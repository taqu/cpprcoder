#ifndef INC_BLKSORT_H_
#define INC_BLKSORT_H_
/**
@file blksort.h

USAGE:
Put '#define BLKSORT_IMPLEMENTATION' before including this file to create the implementation.
*/

#include <cstdint>
#if defined(_MSC_VER)
#    define BLKSORT_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
#    define BLKSORT_RESTRICT __restrict
#else
#    define BLKSORT_RESTRICT
#endif

#ifndef BLKSORT_MALLOC
#    if defined(_MSC_VER)
#        include <malloc.h>
#        define BLKSORT_MALLOC(a, x) _aligned_malloc(x, a)
#    elif defined(__GNUC__) || defined(__clang__)
#        include <cstdlib>
#        define BLKSORT_MALLOC(a, x) std::aligned_alloc(a, x)
#    else
#        error
#    endif
#endif

#ifndef BLKSORT_FREE
#    if defined(_MSC_VER)
#        include <malloc.h>
#        define BLKSORT_FREE(x) _aligned_free(x)
#    elif defined(__GNUC__) || defined(__clang__)
#        include <cstdlib>
#        define BLKSORT_FREE(x) std::free(x)
#    else
#        error
#    endif
#endif

#ifndef BLKSORT_ALIGN
#    if defined(_MSC_VER)
#        define BLKSORT_ALIGN(x) __declspec(align(x))
#    elif defined(__GNUC__) || defined(__clang__)
#        define BLKSORT_ALIGN(x) __attribute__((aligned(x)))
#    else
#        error
#    endif
#endif

#define BLOCKSORT_DEBUG (0)
#define BLOCKSORT_PERF (1)
#define BLOCKSORT_MTF (0)

namespace blksort
{
struct Item
{
    uint8_t* str_;
    uint16_t id_;
    uint16_t head2() const
    {
        return (static_cast<uint16_t>(str_[1]) << 8) | static_cast<uint16_t>(str_[0]);
    }
};

int32_t strcmp(const uint8_t* BLKSORT_RESTRICT x0, const uint8_t* BLKSORT_RESTRICT x1, uint32_t depth);
void insertionsort(uint32_t size, Item* BLKSORT_RESTRICT v, uint32_t depth);
void heapsort(uint32_t n, Item* BLKSORT_RESTRICT v, uint32_t depth);
void sort(uint32_t size, Item* BLKSORT_RESTRICT data, uint32_t depth);

void counting_sort(uint32_t size, uint16_t* BLKSORT_RESTRICT dst, const uint8_t* BLKSORT_RESTRICT key, const uint16_t* BLKSORT_RESTRICT value);

class BlkSort
{
public:
    inline static constexpr uint32_t Align = 16;
    inline static constexpr uint32_t BlockSize = 1024 * 32;
    inline static constexpr uint32_t BlockShift = 15;
    inline static constexpr uint32_t BlockMask = BlockSize - 1;
    inline static constexpr uint32_t EncodedSize = BlockSize + 2;

    BlkSort();
    ~BlkSort();
    static uint32_t encodeBound(uint32_t size);
    static uint32_t decodeBound(uint32_t size);

    void encode(uint32_t size, uint8_t* BLKSORT_RESTRICT dst, const uint8_t* BLKSORT_RESTRICT data);
    void decode(uint32_t size, uint8_t* BLKSORT_RESTRICT dst, uint8_t* BLKSORT_RESTRICT data);

private:
    BlkSort(const BlkSort&) = delete;
    BlkSort& operator=(const BlkSort&) = delete;

    void encode_internal(uint8_t* BLKSORT_RESTRICT dst, const uint8_t* BLKSORT_RESTRICT data);
    void decode_internal(uint8_t* BLKSORT_RESTRICT dst, uint8_t* BLKSORT_RESTRICT data);

    void mtf_init(uint8_t* BLKSORT_RESTRICT id);
    uint8_t mtf_find(const uint8_t* BLKSORT_RESTRICT table, uint8_t x);
    void mtf_move_to_front_one(uint8_t* BLKSORT_RESTRICT table, uint8_t x);
    void mtf_encode(uint32_t size, uint8_t* BLKSORT_RESTRICT data);
    void mtf_decode(uint32_t size, uint8_t* BLKSORT_RESTRICT data);

    const uint32_t size_;
    uint8_t* buffer_;
};

} // namespace blksort
#endif // INC_BLKSORT_H_

#ifdef BLKSORT_IMPLEMENTATION
#include <cassert>
#include <cstring>

#include <algorithm>
#ifdef __AVX__
#    define BLKSORT_AVX (1)
#    include <immintrin.h>
#endif

#ifdef __ARM_NEON
#    define BLKSORT_NEON (1)
#    include <arm_neon.h>
#endif

#if BLOCKSORT_PERF
#    include <chrono>
#endif

#if BLOCKSORT_PERF
double encode_setup = 0.0;
double encode_sort = 0.0;
double encode_finalize = 0.0;
double decode_setup = 0.0;
double decode_sort = 0.0;
double decode_finalize = 0.0;
double mtf_encode_setup = 0.0;
double mtf_encode_main = 0.0;
double mtf_decode_setup = 0.0;
double mtf_decode_main = 0.0;
#endif

namespace blksort
{
namespace
{
#if BLOCKSORT_DEBUG
    void print(uint32_t size, const uint8_t* data)
    {
        for(uint32_t i = 0; i < size; ++i) {
            printf("%02X ", data[i]);
        }
        printf("\n");
    }
#endif

#if 0
    uint32_t bitscan(uint32_t mask)
    {
        unsigned long index = 0;
        return 0 == _BitScanForward(&index, mask) ? 0x0U : index;
    }
#endif

    Item median(uint32_t size, const Item* data)
    {
        uint32_t m = size >> 2;
        uint8_t x0 = data[m].str_[0];
        uint32_t m2 = m + m;
        uint8_t x1 = data[m2].str_[0];
        uint32_t m3 = m + m2;
        uint8_t x2 = data[m3].str_[0];
        if(x0 < x1) {
            return x1 < x2 ? data[m2] : (x0 < x2 ? data[m3] : data[m]);
        } else {
            return x0 < x2 ? data[m] : (x1 < x2 ? data[m3] : data[m2]);
        }
    }

    bool less(const uint8_t* x0, const uint8_t* x1, uint32_t depth)
    {
        assert((depth & 15) == 0);
#if 0
        for(uint32_t d=0; d<depth; d+=16){
            __m128i m0 = _mm_loadu_si128((const __m128i*)&x0[d]);
            __m128i m1 = _mm_loadu_si128((const __m128i*)&x1[d]);
            uint32_t mask0 = (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(m0, m1));
            if(0xFFFFUL != mask0) {
                mask0 = (~mask0) & 0xFFFFUL;
                uint32_t mask1 = (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(m0, _mm_min_epu8(m0, m1)));
                uint32_t msb0 = bitscan(mask0);
                return (mask1 >> msb0) & 0x01U;
            }
            d += 16;
        }
#else
        for(uint32_t d = 0; d < depth; ++d) {
            if(x0[d] != x1[d]) {
                return x0[d] < x1[d];
            }
        }
#endif
        return false;
    }

} // namespace

int32_t strcmp(const uint8_t* x0, const uint8_t* x1, uint32_t depth)
{
    assert((depth & 15) == 0);
    for(uint32_t d = 0; d < depth; ++d) {
        if(x0[d] != x1[d]) {
            return x0[d] < x1[d] ? -1 : 1;
        }
    }
    return 0;
}

void insertionsort(uint32_t size, Item* v, uint32_t depth)
{
    for(uint32_t i = 1; i < size; ++i) {
        Item x = v[i];
        int64_t j;
        for(j = i - 1; 0 <= j && less(x.str_, v[j].str_, depth); --j) {
            v[j + 1] = v[j];
        }
        v[j + 1] = x;
    }
}

void heapsort(uint32_t n, Item* v, uint32_t depth)
{
    --v; // set index = 1
    int32_t i, j;
    Item x;
    for(int32_t k = n >> 1; k >= 1; --k) {
        i = k;
        x = v[k];
        while((j = i << 1) <= n) {
            if(j < n && less(v[j].str_, v[j + 1].str_, depth)) {
                ++j;
            }

            if(!less(x.str_, v[j].str_, depth)) {
                break;
            }
            v[i] = v[j];
            i = j;
        }
        v[i] = x;
    }

    while(n > 1) {
        x = v[n];
        v[n] = v[1];
        --n;
        i = 1;
        while((j = i << 1) <= n) {
            if(j < n && less(v[j].str_, v[j + 1].str_, depth)) {
                ++j;
            }

            if(!less(x.str_, v[j].str_, depth)) {
                break;
            }
            v[i] = v[j];
            i = j;
        }
        v[i] = x;
    }
}

void mqsort(uint32_t size, Item* data, uint32_t d, uint32_t depth, int32_t level)
{
    static constexpr uint32_t SwitchN = 37;
    if(level <= 0) {
        heapsort(size, data, depth);
        return;
    }
    while(d < depth) {
        if(size < SwitchN) {
            insertionsort(size, data, depth);
            return;
        }
        Item pivot = median(size, data);
        const uint8_t p = pivot.str_[d];

        int32_t h = static_cast<int32_t>(size) - 1;
        int32_t i0 = 0;
        int32_t i1 = h;
        int32_t m0 = i0;
        int32_t m1 = i1;
        for(;;) {
            while(i0 <= i1) {
                uint8_t c = data[i0].str_[d];
                if(p < c) {
                    break;
                }
                if(p == c) {
                    std::swap(data[i0], data[m0]);
                    ++m0;
                }
                ++i0;
            }

            while(i0 <= i1) {
                uint8_t c = data[i1].str_[d];
                if(c < p) {
                    break;
                }
                if(p == c) {
                    std::swap(data[i1], data[m1]);
                    --m1;
                }
                --i1;
            }
            if(i1 < i0) {
                break;
            }
            std::swap(data[i0], data[i1]);
            ++i0;
            --i1;
        }
        int32_t r0 = std::min(m0, i0 - m0);
        for(int32_t i = 0; i < r0; ++i) {
            std::swap(data[i], data[i1 - i]);
        }
        m0 = i0 - m0;
        int32_t r1 = std::min(h - m1, m1 - i1);
        for(int32_t i = 0; i < r1; ++i) {
            std::swap(data[i0 + i], data[h - i]);
        }
        m1 = h - (m1 - i1) + 1;
        if(0 < (m0 - 1)) {
            mqsort(static_cast<uint32_t>(m0), data, d, depth, level - 1);
        }
        if(m1 < h) {
            mqsort(static_cast<uint32_t>(static_cast<int32_t>(size) - m1), data + m1, d, depth, level - 1);
        }
        if(m1 <= m0) {
            break;
        }
        data += m0;
        size = static_cast<uint32_t>(m1 - m0);
        ++d;
    }
}

void sort(uint32_t size, Item* data, uint32_t depth)
{
#if 0
    int32_t level = 0;
    uint32_t t = size;
    while(1 < t) {
        ++level;
        t >>= 1;
    }
#else
    const int32_t level = 11;
#endif
    mqsort(size, data, 0, depth, level);
}

void counting_sort(uint32_t size, uint16_t* dst, const uint8_t* key, const uint16_t* value)
{
    assert(0 == (size & 15));
    BLKSORT_ALIGN(16)
    uint16_t count[259];
    ::memset(count, 0, 256 * sizeof(uint16_t));
    for(uint32_t i = 0; i < size; i += 4) {
        count[key[i + 0]] += 1;
        count[key[i + 1]] += 1;
        count[key[i + 2]] += 1;
        count[key[i + 3]] += 1;
    }

    for(uint32_t i = 1; i < 256; i += 4) {
        uint16_t c0 = count[i - 1];
        uint16_t c1 = count[i];
        uint16_t c2 = count[i + 1];
        uint16_t c3 = count[i + 2];
        count[i] += c0;
        c1 += c0;
        count[i + 1] += c1;
        c2 += c1;
        count[i + 2] += c2;
        c3 += c2;
        count[i + 3] += c3;
    }

    for(int32_t i = static_cast<int32_t>(size - 1); 0 <= i; --i) {
        uint8_t j0 = key[i];
        --count[j0];
        dst[count[j0]] = value[i];
    }
}

BlkSort::BlkSort()
    : size_(BlockSize)
    , buffer_(nullptr)
{
    size_t encodeSize = (sizeof(uint8_t) * 2 + sizeof(Item)) * size_;
    size_t decodeSize = sizeof(uint16_t) * 2 * size_;
    size_t mtfSize = sizeof(uint8_t) * 256;
    size_t workSize = std::max(encodeSize, std::max(decodeSize, mtfSize));
    buffer_ = (uint8_t*)BLKSORT_MALLOC(Align, workSize);
}

BlkSort::~BlkSort()
{
    BLKSORT_FREE(buffer_);
    buffer_ = nullptr;
}

uint32_t BlkSort::encodeBound(uint32_t size)
{
    uint32_t blocks = size >> BlockShift;
    uint32_t rawSize = size - (blocks << BlockShift);
    return blocks * EncodedSize + rawSize;
}

uint32_t BlkSort::decodeBound(uint32_t size)
{
    uint32_t blocks = size >> BlockShift;
    uint32_t rawSize = size - (blocks << BlockShift);
    return blocks * BlockSize + rawSize;
}

void BlkSort::encode(uint32_t size, uint8_t* BLKSORT_RESTRICT dst, const uint8_t* BLKSORT_RESTRICT data)
{
    uint32_t blocks = size >> BlockShift;
    uint32_t rawSize = size - blocks * BlockSize;
    for(uint32_t i = 0; i < blocks; ++i) {
        encode_internal(dst, data);
        dst += EncodedSize;
        data += BlockSize;
    }
    ::memcpy(dst, data, rawSize);
}

void BlkSort::decode(uint32_t size, uint8_t* BLKSORT_RESTRICT dst, uint8_t* BLKSORT_RESTRICT data)
{
    uint32_t blocks = size / EncodedSize;
    uint32_t rawSize = size - blocks * EncodedSize;

    for(uint32_t i = 0; i < blocks; ++i) {
        decode_internal(dst, data);
        dst += BlockSize;
        data += EncodedSize;
    }
    ::memcpy(dst, data, rawSize);
}

void BlkSort::encode_internal(uint8_t* BLKSORT_RESTRICT dst, const uint8_t* BLKSORT_RESTRICT data)
{
#if BLOCKSORT_PERF
    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();
#endif
    uint8_t* buffer = buffer_;
    ::memcpy(buffer, data, size_);
    ::memcpy(buffer + size_, data, size_);

    Item* strings = (Item*)(buffer_ + size_ + size_);
    {
        for(int32_t i = 0; i < size_; i += 8) {
            strings[i + 0].id_ = i + 0;
            strings[i + 1].id_ = i + 1;
            strings[i + 2].id_ = i + 2;
            strings[i + 3].id_ = i + 3;
            strings[i + 0].str_ = &buffer[i + 0];
            strings[i + 1].str_ = &buffer[i + 1];
            strings[i + 2].str_ = &buffer[i + 2];
            strings[i + 3].str_ = &buffer[i + 3];

            strings[i + 4].id_ = i + 4;
            strings[i + 5].id_ = i + 5;
            strings[i + 6].id_ = i + 6;
            strings[i + 7].id_ = i + 7;
            strings[i + 4].str_ = &buffer[i + 4];
            strings[i + 5].str_ = &buffer[i + 5];
            strings[i + 6].str_ = &buffer[i + 6];
            strings[i + 7].str_ = &buffer[i + 7];
        }
    }
#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    encode_setup += std::chrono::duration<double>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
#endif

    sort(size_, strings, size_);

#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    encode_sort += std::chrono::duration<double>(end - start).count();
#endif
#if BLOCKSORT_DEBUG
    for(uint32_t i = 1; i < size_; ++i) {
        assert(szsort::strcmp(strings[i - 1].str_, strings[i].str_, size_) < 0);
    }
    for(uint32_t i = 0; i < size_; ++i) {
        for(uint32_t j = i + 1; j < size_; ++j) {
            assert(strings[i].id_ != strings[j].id_);
        }
    }
#endif
#if BLOCKSORT_PERF
    start = std::chrono::high_resolution_clock::now();
#endif
    uint16_t pos = 0;
    for(uint32_t i = 0; i < size_; ++i) {
        dst[i] = strings[i].str_[size_ - 1];
        if(0 == strings[i].id_) {
            assert(0 == ::memcmp(strings[i].str_, data, size_));
            pos = i;
        }
    }
    ::memcpy(dst + size_, &pos, sizeof(uint16_t));
#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    encode_finalize += std::chrono::duration<double>(end - start).count();
#endif

#if BLOCKSORT_DEBUG
    print(16, data);
    for(uint32_t i = 0; i < 16; ++i) {
        const szsort::Item& item = strings[i];
        printf("[%d] ", i);
        print(16, item.str_ + Size - 16);
    }
    print(16, dst);
    printf("%d\n", pos);
#endif
#if BLOCKSORT_MTF
    mtf_encode(BlockSize, dst);
#endif
}

void BlkSort::decode_internal(uint8_t* BLKSORT_RESTRICT dst, uint8_t* BLKSORT_RESTRICT data)
{
#if BLOCKSORT_MTF
    mtf_decode(BlockSize, data);
#endif
#if BLOCKSORT_PERF
    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();
#endif
    uint16_t* id = (uint16_t*)buffer_;
    // clang-format off
    BLKSORT_ALIGN(Align) static const uint16_t ID0[8] = {0,1,2,3,4,5,6,7};
    BLKSORT_ALIGN(Align) static const uint16_t ID1[8] = {8,9,10,11,12,13,14,15};
    BLKSORT_ALIGN(Align) static const uint16_t ID2[8] = {16,17,18,19,20,21,22,23};
    BLKSORT_ALIGN(Align) static const uint16_t ID3[8] = {24,25,26,27,28,29,30,31};
    // clang-format on
#if defined(BLKSORT_AVX)
#    if 0
        if (size_ <= 32) {
            __m128i c0 = _mm_load_si128((const __m128i*)ID0);
            __m128i c1 = _mm_load_si128((const __m128i*)ID1);
            __m128i add = _mm_set1_epi16(16);
            for(uint32_t i = 0; i < size_; i += 16) {
                _mm_store_si128((__m128i*)&id[i], c0);
                c0 = _mm_adds_epi16(c0, add);
                _mm_store_si128((__m128i*)&id[i + 8], c1);
                c1 = _mm_adds_epi16(c1, add);
            }

        } else
#    else
    __m128i c0 = _mm_load_si128((const __m128i*)ID0);
    __m128i c1 = _mm_load_si128((const __m128i*)ID1);
    __m128i c2 = _mm_load_si128((const __m128i*)ID2);
    __m128i c3 = _mm_load_si128((const __m128i*)ID3);
    __m128i add = _mm_set1_epi16(32);
    for(uint32_t i = 0; i < size_; i += 32) {
        _mm_store_si128((__m128i*)&id[i], c0);
        c0 = _mm_adds_epi16(c0, add);
        _mm_store_si128((__m128i*)&id[i + 8], c1);
        c1 = _mm_adds_epi16(c1, add);
        _mm_store_si128((__m128i*)&id[i + 16], c2);
        c2 = _mm_adds_epi16(c2, add);
        _mm_store_si128((__m128i*)&id[i + 24], c3);
        c3 = _mm_adds_epi16(c3, add);
    }
#    endif
#elif defined(BLKSORT_NEON)
    uint16x8_t c0 = vld1q_u16(ID0);
    uint16x8_t c1 = vld1q_u16(ID1);
    uint16x8_t c2 = vld1q_u16(ID2);
    uint16x8_t c3 = vld1q_u16(ID3);
    uint16x8_t add = vmovq_n_u16(32);
    for(uint32_t i = 0; i < size_; i += 32) {
        vst1q_u16(&id[i], c0);
        c0 = vaddq_u16(c0, add);
        vst1q_u16(&id[i + 8], c1);
        c1 = vaddq_u16(c1, add);
        vst1q_u16(&id[i + 16], c2);
        c2 = vaddq_u16(c2, add);
        vst1q_u16(&id[i + 24], c3);
        c3 = vaddq_u16(c3, add);
    }

#else
    for(uint32_t i = 0; i < size_; i += 4) {
        id[i + 0] = i + 0;
        id[i + 1] = i + 1;
        id[i + 2] = i + 2;
        id[i + 3] = i + 3;
    }
#endif

#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    decode_setup += std::chrono::duration<double>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
#endif

    uint16_t* out_id = (uint16_t*)(buffer_ + size_ + size_);
    counting_sort(size_, out_id, data, id);

#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    decode_sort += std::chrono::duration<double>(end - start).count();
#endif
#if BLOCKSORT_DEBUG
    for(uint32_t i = 0; i < Size; ++i) {
        for(uint32_t j = i + 1; j < Size; ++j) {
            assert(out_id[i] != out_id[j]);
        }
    }
    for(uint32_t i = 1; i < Size; ++i) {
        assert(data[out_id[i - 1]] <= data[out_id[i]]);
    }
#endif

#if BLOCKSORT_PERF
    start = std::chrono::high_resolution_clock::now();
#endif
    uint16_t top;
    ::memcpy(&top, data + size_, sizeof(uint16_t));

    uint16_t p = out_id[top];
    for(uint32_t i = 0; i < size_; ++i) {
        dst[i] = data[p];
        p = out_id[p];
    }
#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    decode_finalize += std::chrono::duration<double>(end - start).count();
#endif
#if BLOCKSORT_DEBUG
    print(16, data);
    print(16, dst);
#endif
}

void BlkSort::mtf_init(uint8_t* BLKSORT_RESTRICT id)
{
    static BLKSORT_ALIGN(Align) const uint8_t ID0[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    static BLKSORT_ALIGN(Align) const uint8_t ID1[16] = {16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
    static BLKSORT_ALIGN(Align) const uint8_t ID2[16] = {32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47};
    static BLKSORT_ALIGN(Align) const uint8_t ID3[16] = {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
#ifdef BLKSORT_AVX
    __m128i c0 = _mm_load_si128((const __m128i*)ID0);
    __m128i c1 = _mm_load_si128((const __m128i*)ID1);
    __m128i c2 = _mm_load_si128((const __m128i*)ID2);
    __m128i c3 = _mm_load_si128((const __m128i*)ID3);
    __m128i add = _mm_set1_epi8(64);
    for(uint32_t i = 0; i < 256; i += 64) {
        _mm_store_si128((__m128i*)&id[i], c0);
        c0 = _mm_add_epi8(c0, add);
        _mm_store_si128((__m128i*)&id[i + 16], c1);
        c1 = _mm_add_epi8(c1, add);
        _mm_store_si128((__m128i*)&id[i + 32], c2);
        c2 = _mm_add_epi8(c2, add);
        _mm_store_si128((__m128i*)&id[i + 48], c3);
        c3 = _mm_add_epi8(c3, add);
    }
#else
    for(uint32_t i = 0; i < 256; i += 4) {
        id[i + 0] = i + 0;
        id[i + 1] = i + 1;
        id[i + 2] = i + 2;
        id[i + 3] = i + 3;
    }
#endif
}

uint8_t BlkSort::mtf_find(const uint8_t* BLKSORT_RESTRICT table, uint8_t x)
{
#ifdef BLKSORT_AVX
    __m128i c = _mm_set1_epi8(*(char*)&x);
    for(uint32_t i = 0; i < 256; i += 16) {
        __m128i str0 = _mm_load_si128((const __m128i*)&table[i]);
        int32_t r0 = _mm_cmpestri(c, 1, str0, 16, 0);
        if(r0 <= 0xF) {
            assert(x == table[i + r0]);
            assert((i + r0) < 256);
            return static_cast<uint8_t>(i + r0);
        }
    }
#else
    for(uint32_t i = 0; i < 256; ++i) {
        if(table[i] == x) {
            return static_cast<uint8_t>(i);
        }
    }
#endif
    return -1;
}

void BlkSort::mtf_move_to_front_one(uint8_t* BLKSORT_RESTRICT table, uint8_t x)
{
    assert(1 < x);
    uint8_t c = table[x];
    ::memmove(table + 2, table + 1, x - 1);
    table[1] = c;
}

void BlkSort::mtf_encode(uint32_t size, uint8_t* BLKSORT_RESTRICT data)
{
#if BLOCKSORT_PERF
    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();
#endif
    BLKSORT_ALIGN(Align) uint8_t table[256];
    mtf_init(table);
#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    mtf_encode_setup += std::chrono::duration<double>(end - start).count();
    start = std::chrono::high_resolution_clock::now();
#endif
    uint8_t prev = 1;
    for(uint32_t i = 0; i < size; ++i) {
        uint8_t c = data[i];
        uint8_t x = mtf_find(table, c);
        if(1 == x) {
            if(0 != prev) {
                table[1] = table[0];
                table[0] = c;
            }
        } else if(1 < x) {
            mtf_move_to_front_one(table, x);
        }
        data[i] = x;
        prev = x;
    }
#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    mtf_encode_main += std::chrono::duration<double>(end - start).count();
#endif
}

void BlkSort::mtf_decode(uint32_t size, uint8_t* BLKSORT_RESTRICT data)
{
#if BLOCKSORT_PERF
    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();
#endif
    BLKSORT_ALIGN(Align) uint8_t table[256];
    mtf_init(table);
#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    mtf_decode_setup += std::chrono::duration<double>(end - start).count();
    start = std::chrono::high_resolution_clock::now();
#endif

    uint8_t prev = 1;
    for(uint32_t i = 0; i < size; ++i) {
        uint8_t x0 = data[i + 0];
        uint8_t c0 = table[x0];
        data[i + 0] = c0;
        if(1 == x0) {
            if(0 != prev) {
                table[1] = table[0];
                table[0] = c0;
            }
        } else if(1 < x0) {
            mtf_move_to_front_one(table, x0);
        }
        prev = x0;
    }
#if BLOCKSORT_PERF
    end = std::chrono::high_resolution_clock::now();
    mtf_decode_main += std::chrono::duration<double>(end - start).count();
#endif
}

} // namespace blksort
#endif // BLKSORT_IMPLEMENTATION
