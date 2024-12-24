#ifndef INC_CPPASE_H_
#define INC_CPPASE_H_
/*
USAGE:
Put '#define CPPASE_IMPLEMENTATION' before including this file to create the implementation.
*/
//#define CPPASE_IMPLEMENTATION(1)

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <immintrin.h>
#include <random>

using int64 = int64_t;

using uint8 = uint8_t;
using uint16 = uint16_t;
using uint32 = uint32_t;
using uint64 = uint64_t;

#if defined(__cplusplus) && 202002L <= __cplusplus && !defined(LIKELY_)
#    define LIKELY_ [[likely]]
#else
#    define LIKELY_
#endif

#if defined(__cplusplus) && 202002L <= __cplusplus && !defined(UNLIKELY_)
#    define UNLIKELY_ [[unlikely]]
#else
#    define UNLIKELY_
#endif

#if defined(__cplusplus) && 202002L <= __cplusplus
#    include <bit>
template<typename To, typename From>
inline To bitcast(const From& x)
{
    return std::bit_cast<To, From>(x);
}
#else
#    include <intrin.h>
template<typename To, typename From>
inline To bitcast(const From& x)
{
    To r;
    ::memcpy(&r, &x, sizeof(To));
    return r;
}
#endif

#ifdef _MSC_VER
uint32 bitscan_lsb(uint32 mask)
{
    unsigned long index = 0;
    if(0 != _BitScanForward(&index, mask)) {
        return index + 1;
    }
    return 0;
}
#else
uint32 bitscan_lsb(uint32 mask)
{
    return 0 != mask ? __builtin_ctz(mask) + 1 : 0;
}
#endif
namespace cppase
{

class ASE
{
public:
    inline static constexpr uint32 EncodeBits = 6;
    inline static constexpr uint32 TableSize = 1 << EncodeBits;

    class Stream
    {
    public:
        Stream(uint8* data)
            : pos_(0)
            , r_(8)
            , data_(data)
        {
            data_[0] = 0;
        }
        void put(uint32 bits, uint32 size);
        uint32 pos_;
        uint32 r_;
        uint8* data_;
    };

    class IStream
    {
    public:
        IStream(const uint8* data)
            : pos_(0)
            , r_(8)
            , data_(data)
        {
        }
        uint32 get(uint32 size);
        uint32 pos_;
        uint32 r_;
        const uint8* data_;
    };

    inline static constexpr uint32 Invalid = 0xFFFF'FFFFUL;

    ASE();
    ~ASE();
    void clear();
    bool empty() const;
    bool full() const;
    /**
     * @brief
     * @return 1 base index
     */
    uint32 find(uint8 x) const;
    void push_back(uint8 x);
    void remove_at(uint32 index);
    uint32 entropy() const;
    static int64 encodeBound(int64 size);
    int64 encode(int64 size, uint8* dst, const uint8* src);
    int64 decode(int64 size, int64 original_size, uint8* dst, const uint8* src);

private:
    ASE(const ASE&) = delete;
    ASE& operator=(const ASE&) = delete;
    uint32 size_;
    uint8 symbols_[TableSize + 16];
};
}
#endif ////INC_CPPASE_H_

#ifdef CPPASE_IMPLEMENTATION
namespace cppase
{
void ASE::Stream::put(uint32 bits, uint32 size)
{
    assert(size <= 32);
    uint32 t = 0;
    while(t < size) {
        if(r_ <= 0) {
            ++pos_;
            data_[pos_] = 0;
            r_ = 8;
        }
        uint32 s = (size - t) < r_ ? (size - t) : r_;
        data_[pos_] |= (bits & ((1UL << s) - 1)) << (8 - r_);
        assert(s <= r_);
        r_ -= s;
        t += s;
        bits >>= s;
    }
}

uint32 ASE::IStream::get(uint32 size)
{
    assert(size <= 32);
    uint32 x = 0;
    uint32 t = 0;
    while(t < size) {
        if(r_ <= 0) {
            ++pos_;
            r_ = 8;
        }
        uint32 s = (size - t) < r_ ? (size - t) : r_;
        x |= ((data_[pos_] >> (8 - r_)) & ((1UL << s) - 1)) << t;
        assert(s <= r_);
        t += s;
        r_ -= s;
    }
    return x;
}

ASE::ASE()
    : size_(0)
{
}

ASE::~ASE()
{
}

void ASE::clear()
{
    size_ = 0;
}

bool ASE::empty() const
{
    return size_ <= 0;
}

bool ASE::full() const
{
    return TableSize <= size_;
}

uint32 ASE::find(uint8 x) const
{
    if(size_ <= 0)
        UNLIKELY_
        {
            return Invalid;
        }
    __m128i m = _mm_set1_epi8(bitcast<char, uint8>(x));
    for(uint32 i = 0; i < size_; i += 16) {
        __m128i m0 = _mm_loadu_si128((const __m128i*)&symbols_[i + 0]);
        m0 = _mm_cmpeq_epi8(m0, m);
        uint32 mask = (uint32)_mm_movemask_epi8(m0);
        uint32 index = bitscan_lsb(mask);
        if(index != 0) {
            index += i;
            return index < size_ ? index : Invalid;
        }
    }
    return Invalid;
}

void ASE::push_back(uint8 x)
{
    assert(!full());
    symbols_[size_] = x;
    ++size_;
}

void ASE::remove_at(uint32 index)
{
    assert(index < size_);
    for(uint32 i = index + 1; i < size_; i += 16) {
        __m128i m0 = _mm_loadu_si128((const __m128i*)&symbols_[i + 0]);
        _mm_storeu_si128((__m128i*)&symbols_[i - 1], m0);
    }
    --size_;
}

uint32 ASE::entropy() const
{
    assert(size_ <= TableSize);
    uint32 e = 0;
    uint32 x = 0;
    uint32 size = size_;
    while(x < size) {
        ++e;
        x = 1UL << e;
    }
    return e;
}

int64 ASE::encodeBound(int64 size)
{
    return (static_cast<int64>(std::ceil(size * (9.0 / 8.0))) + 0x04LL) & ~0x03LL;
}

int64 ASE::encode(int64 size, uint8* dst, const uint8* src)
{
    clear();
    Stream stream(dst);
    uint32 bits = 0;
    for(int64 i = 0; i < size; ++i) {
        uint32 index = find(src[i]);
        if(index != Invalid) {
            assert(1 <= index && index <= size_);
            uint32 out = size_ - index;
            --index;
            assert(src[i] == symbols_[index]);
            remove_at(index);
            push_back(src[i]);
            uint32 c = (out << 1) | 0x01U;
            stream.put(c, bits + 1);
        } else {
            if(TableSize <= size_) {
                remove_at(0);
                push_back(src[i]);
            } else {
                push_back(src[i]);
                bits = entropy();
            }
            uint32 c = static_cast<uint32>(src[i]) << 1;
            stream.put(c, 9);
        }
    }
    return static_cast<int64>(stream.pos_);
}

int64 ASE::decode(int64 size, int64 original_size, uint8* dst, const uint8* src)
{
    clear();
    IStream stream(src);
    uint32 s = static_cast<uint32>(size);
    uint32 bits = 0;
    int64 pos = 0;
    while(stream.pos_ <= s && pos < original_size) {
        uint32 cmark = stream.get(1);
        if(cmark) {
            uint32 index = stream.get(bits);
            assert(index < size_);
            index = size_ - index - 1;
            assert(index < size_);
            dst[pos] = symbols_[index];
            remove_at(index);
            push_back(dst[pos]);
            ++pos;
        } else {
            uint8 c = static_cast<uint8>(stream.get(8));
#if _DEBUG
            uint32 index = find(c);
            assert(Invalid == index);
#endif
            if(TableSize <= size_) {
                remove_at(0);
                push_back(c);
            } else {
                push_back(c);
                bits = entropy();
            }
            dst[pos] = c;
            ++pos;
        }
    }
    return pos;
}
}
#endif
