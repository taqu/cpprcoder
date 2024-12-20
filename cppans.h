#ifndef INC_CPPANS_H_
#define INC_CPPANS_H_
/*
USAGE:
Put '#define CPPANS_IMPLEMENTATION' before including this file to create the implementation.
*/
#define CPPANS_IMPLEMENTATION (1)
#include <cstdint>
#include <immintrin.h>

namespace cppans
{
using s8 = int8_t;
using s16 = int16_t;
using s32 = int32_t;
using s64 = int64_t;

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

class rANS
{
public:
    inline static constexpr u32 MaxSize = 0x7FFF'FFFFUL;
    inline static constexpr u32 ProbBits = 14;
    inline static constexpr u32 ProbScale = 1 << ProbBits;
    inline static constexpr u32 rANSByteLowBounds = 1UL << 23;
    inline static constexpr u32 WordLowBounds = 1UL<<16;
    inline static constexpr u32 WordScaleBits = 12;
    inline static constexpr u32 WordM = 1 << WordScaleBits;

    using State = u32;

    struct EncSymbol
    {
        u32 x_max_; //!< upper bound of pre-normalization interval
        u32 rcp_freq_;
        u32 bias_;
        u16 cmpl_freq_; //!< Complemt of frequency: (!<<scale_bits) - freq
        u16 rcp_shift_;
    };

    struct DecSymbol
    {
        u16 start_;
        u16 freq_;
    };

    union WordSlot
    {
        u32 x_;
        struct{
            u16 freq_;
            u16 bias_;
        };
    };
    struct WordTables
    {
        WordSlot slots_[WordM];
        u8 slot2symbol_[WordM];
    };

    union SimdState
    {
        __m128i simd_;
        u32 lane_[4];
    };

    static u64 calc_encoded_size(u32 size);
    static u32 encode(u32 dst_size, u8* dst, u32 src_size, const u8* src);
    static u32 decode(u32 dst_size, u8* dst, u32 src_size, const u8* src);

    static u32 encode_simd(u32 dst_size, u8* dst, u32 src_size, const u8* src);
    static u32 decode_simd(u32 dst_size, u8* dst, u32 src_size, const u8* src);
private:
    rANS(const rANS&) = delete;
    rANS& operator=(const rANS&) = delete;
};
} // namespace cppans

#ifdef CPPANS_IMPLEMENTATION
#    include <cassert>
#    include <cstring>
#ifdef _DEBUG
#include <stdio.h>
#endif
#    if defined(_MSC_VER)
#        define CPPANS_RESTRICT __restrict
#    elif defined(__gnuc__)
#        define CPPANS_RESTRICT __restrict
#    elif defined(__clang__)
#        define CPPANS_RESTRICT __restrict
#    else
#    endif

namespace cppans
{
namespace
{
    void count(u32* CPPANS_RESTRICT dst, u32 size, const u8* CPPANS_RESTRICT src)
    {
        ::memset(dst, 0, 256 * sizeof(u32));
        u32 s = size & ~0x0FUL;
        for(u32 i = 0; i < s; i += 16) {
            ++dst[src[i + 0]];
            ++dst[src[i + 1]];
            ++dst[src[i + 2]];
            ++dst[src[i + 3]];
            ++dst[src[i + 4]];
            ++dst[src[i + 5]];
            ++dst[src[i + 6]];
            ++dst[src[i + 7]];

            ++dst[src[i + 8]];
            ++dst[src[i + 9]];
            ++dst[src[i + 10]];
            ++dst[src[i + 11]];
            ++dst[src[i + 12]];
            ++dst[src[i + 13]];
            ++dst[src[i + 14]];
            ++dst[src[i + 15]];
        }
        for(u32 i = s; i < size; ++i) {
            ++dst[src[i]];
        }
    }

    void cumulative(u32* CPPANS_RESTRICT dst, const u32* CPPANS_RESTRICT src)
    {
        dst[0] = 0;
        for(u32 i = 0; i < 256; ++i) {
            dst[i + 1] = dst[i] + src[i];
        }
    }

    void normalize(u32* CPPANS_RESTRICT freqs, u32* CPPANS_RESTRICT cum_freqs, u64 target_total)
    {
        u32 current_total = cum_freqs[256];
        for(u32 i = 1; i < 257; ++i) {
            cum_freqs[i] = (target_total * cum_freqs[i]) / current_total;
        }
        for(u32 i = 0; i < 256; ++i) {
            if(freqs[i] && cum_freqs[i + 1] == cum_freqs[i]) {
                u32 best_freq = ~0UL;
                s32 best_steal = -1;
                for(s32 j = 0; j < 256; ++j) {
                    u32 freq = cum_freqs[j + 1] - cum_freqs[j];
                    if(1 < freq && freq < best_freq) {
                        best_freq = freq;
                        best_steal = j;
                    }
                }
                assert(-1 != best_steal);
                if(static_cast<u32>(best_steal) < i) {
                    for(s32 j = best_steal + 1; j <= i; ++j) {
                        --cum_freqs[j];
                    }
                } else {
                    assert(i < best_steal);
                    for(s32 j = i + 1; j <= best_steal; ++j) {
                        ++cum_freqs[j];
                    }
                }
            }
        }
        for(u32 i = 0; i < 256; ++i) {
#    if _DEBUG
            if(0 == freqs[i]) {
                assert(cum_freqs[i] == cum_freqs[i + 1]);
            } else {
                assert(cum_freqs[i] < cum_freqs[i + 1]);
            }
#    endif
            freqs[i] = cum_freqs[i + 1] - cum_freqs[i];
        }
    }

    void init(rANS::EncSymbol& symbol, u32 start, u32 freq, u32 scale_bits)
    {
        assert(scale_bits <= 16);
        assert(start <= (1u << scale_bits));
        assert(freq <= (1u << scale_bits) - start);

        // Say M := 1 << scale_bits.
        //
        // The original encoder does:
        //   x_new = (x/freq)*M + start + (x%freq)
        //
        // The fast encoder does (schematically):
        //   q     = mul_hi(x, rcp_freq) >> rcp_shift   (division)
        //   r     = x - q*freq                         (remainder)
        //   x_new = q*M + bias + r                     (new x)
        // plugging in r into x_new yields:
        //   x_new = bias + x + q*(M - freq)
        //        =: bias + x + q*cmpl_freq             (*)
        //
        // and we can just precompute cmpl_freq. Now we just need to
        // set up our parameters such that the original encoder and
        // the fast encoder agree.

        symbol.x_max_ = ((rANS::rANSByteLowBounds >> scale_bits) << 8) * freq;
        symbol.cmpl_freq_ = static_cast<u16>((1 << scale_bits) - freq);
        if(freq < 2) {
            // freq=0 symbols are never valid to encode, so it doesn't matter what
            // we set our values to.
            //
            // freq=1 is tricky, since the reciprocal of 1 is 1; unfortunately,
            // our fixed-point reciprocal approximation can only multiply by values
            // smaller than 1.
            //
            // So we use the "next best thing": rcp_freq=0xffffffff, rcp_shift=0.
            // This gives:
            //   q = mul_hi(x, rcp_freq) >> rcp_shift
            //     = mul_hi(x, (1<<32) - 1)) >> 0
            //     = floor(x - x/(2^32))
            //     = x - 1 if 1 <= x < 2^32
            // and we know that x>0 (x=0 is never in a valid normalization interval).
            //
            // So we now need to choose the other parameters such that
            //   x_new = x*M + start
            // plug it in:
            //     x*M + start                   (desired result)
            //   = bias + x + q*cmpl_freq        (*)
            //   = bias + x + (x - 1)*(M - 1)    (plug in q=x-1, cmpl_freq)
            //   = bias + 1 + (x - 1)*M
            //   = x*M + (bias + 1 - M)
            //
            // so we have start = bias + 1 - M, or equivalently
            //   bias = start + M - 1.
            symbol.rcp_freq_ = ~0u;
            symbol.rcp_shift_ = 0;
            symbol.bias_ = start + (1 << scale_bits) - 1;
        } else {
            // Alverson, "Integer Division using reciprocals"
            // shift=ceil(log2(freq))
            u32 shift = 0;
            while(freq > (1UL << shift)) {
                shift++;
            }

            symbol.rcp_freq_ = static_cast<u32>(((1ULL << (shift + 31)) + freq - 1) / freq);
            symbol.rcp_shift_ = shift - 1;

            // With these values, 'q' is the correct quotient, so we
            // have bias=start.
            symbol.bias_ = start;
        }
    }

    void init(rANS::DecSymbol& symbol, u32 start, u32 freq)
    {
        assert(start <= (1 << 16));
        assert(freq <= (1 << 16) - start);
        symbol.start_ = static_cast<u16>(start);
        symbol.freq_ = static_cast<u16>(freq);
    }

    inline void init(rANS::State& state)
    {
        state = rANS::rANSByteLowBounds;
    }

    void put(rANS::State& r, u8*& dst, const rANS::EncSymbol& symbol)
    {
        assert(0 < symbol.x_max_);

        // renormalize
        u32 x = r;
        u32 x_max = symbol.x_max_;
        if(x_max <= x) {
            u8* ptr = dst;
            do {
                *--ptr = static_cast<u8>(x & 0xFFUL);
                x >>= 8;
            } while(x_max <= x);
            dst = ptr;
        }

        // x = C(s,x)
        // NOTE: written this way so we get a 32-bit "multiply high" when
        // available. If you're on a 64-bit platform with cheap multiplies
        // (e.g. x64), just bake the +32 into rcp_shift.
        u32 q = static_cast<u32>((static_cast<uint64_t>(x) * symbol.rcp_freq_) >> 32) >> symbol.rcp_shift_;
        r = x + symbol.bias_ + q * symbol.cmpl_freq_;
    }

    void flush(u8*& dst, const rANS::State& r)
    {
        u32 x = r;
        u8* ptr = dst;
        ptr -= 4;
        ptr[0] = static_cast<u8>(x >> 0);
        ptr[1] = static_cast<u8>(x >> 8);
        ptr[2] = static_cast<u8>(x >> 16);
        ptr[3] = static_cast<u8>(x >> 24);
        dst = ptr;
    }

    // Initializes a rANS decoder.
    // Unlike the encoder, the decoder works forwards as you'd expect.
    void init_decode(rANS::State& r, const u8*& ptr)
    {
        r = ptr[0] << 0;
        r |= ptr[1] << 8;
        r |= ptr[2] << 16;
        r |= ptr[3] << 24;
        ptr += 4;
    }

    // Returns the current cumulative frequency (map it to a symbol yourself!)
    inline u32 get(rANS::State& r, u32 scale_bits)
    {
        return r & ((1UL << scale_bits) - 1);
    }

    // Advances in the bit stream by "popping" a single symbol with range start
    // "start" and frequency "freq". All frequencies are assumed to sum to "1 << scale_bits",
    // and the resulting bytes get written to ptr (which is updated).
    void advance(rANS::State& r, const u8*& ptr, u32 start, u32 freq, u32 scale_bits)
    {
        u32 mask = (1UL << scale_bits) - 1;
        // s, x = D(x)
        u32 x = r;
        x = freq * (x >> scale_bits) + (x & mask) - start;
        // renormalize
        if(x < rANS::rANSByteLowBounds) {
            do {
                x = (x << 8) | *ptr++;
            } while(x < rANS::rANSByteLowBounds);
        }
        r = x;
    }

    rANS::State wordEncInit()
    {
        return rANS::WordLowBounds;
    }

    // Initialize slots for a symbol in the table
void initSymbols(rANS::WordTables& tables, u8 sym, u32 start, u32 freq)
{
    for (u32 i=0; i < freq; ++i) {
        u32 slot = start + i;
        assert(slot<rANS::WordM);
        tables.slot2symbol_[slot] = sym;
        tables.slots_[slot].freq_ = static_cast<u16>(freq);
        tables.slots_[slot].bias_ = static_cast<u16>(i);
    }
}

void wordEncPut(rANS::State& r, u16*& ptr, u32 start, u32 freq)
{
    // renormalize
    u32 x = r;
    if (((rANS::WordLowBounds >> rANS::WordScaleBits) << 16) * freq <= x) {
        ptr -= 1;
        *ptr = static_cast<u16>(x & 0xFFFFUL);
        x >>= 16;
    }
    // x = C(s,x)
    r = ((x / freq) << rANS::WordScaleBits) + (x % freq) + start;
}

// Flushes the rANS encoder
static inline void wordEncFlush(rANS::State& r, u16*& ptr)
{
    u32 x = r;
    ptr -= 2;
    ptr[0] = static_cast<u16>(x >> 0);
    ptr[1] = static_cast<u16>(x >> 16);
}

// Initializes a rANS decoder.
void wordDecInit(rANS::State& r, const u16*& ptr)
{
    r  = ptr[0] << 0;
    r |= ptr[1] << 16;
    ptr += 2;
}

// Decodes a symbol using the given tables.
u8 wordDecSym(rANS::State& r, const rANS::WordTables& table)
{
    u32 x = r;
    u32 slot = x & (rANS::WordM - 1);

    // s, x = D(x)
    r = table.slots_[slot].freq_ * (x >> rANS::WordScaleBits) + table.slots_[slot].bias_;
    return table.slot2symbol_[slot];
}

// Renormalize after decoding a symbol.
void wordDecRenorm(rANS::State& r, const u16*& ptr)
{
    u32 x = r;
    if (x < rANS::WordLowBounds) {
        r = (x << 16) | *ptr;
        ptr += 1;
    }
}

// Initializes a SIMD rANS decoder.
void simdDecInit(rANS::SimdState& r, const u16*& ptr)
{
    r.simd_ = _mm_loadu_si128((const __m128i*)ptr);
    ptr += 2*4;
}

// Decodes a four symbols in parallel using the given tables.
u32 simdDecSym(rANS::SimdState& r, const rANS::WordTables& tables)
{
    __m128i freq_bias_lo, freq_bias_hi, freq_bias;
    __m128i freq, bias;
    __m128i xscaled;
    __m128i x = r.simd_;
    __m128i slots = _mm_and_si128(x, _mm_set1_epi32(rANS::WordM - 1));
    u32 i0 = (u32) _mm_cvtsi128_si32(slots);
    u32 i1 = (u32) _mm_extract_epi32(slots, 1);
    u32 i2 = (u32) _mm_extract_epi32(slots, 2);
    u32 i3 = (u32) _mm_extract_epi32(slots, 3);

    // symbol
    u32 s = tables.slot2symbol_[i0] | (tables.slot2symbol_[i1] << 8) | (tables.slot2symbol_[i2] << 16) | (tables.slot2symbol_[i3] << 24);

    // gather freq_bias
    freq_bias_lo = _mm_cvtsi32_si128(tables.slots_[i0].x_);
    freq_bias_lo = _mm_insert_epi32(freq_bias_lo, tables.slots_[i1].x_, 1);
    freq_bias_hi = _mm_cvtsi32_si128(tables.slots_[i2].x_);
    freq_bias_hi = _mm_insert_epi32(freq_bias_hi, tables.slots_[i3].x_, 1);
    freq_bias = _mm_unpacklo_epi64(freq_bias_lo, freq_bias_hi);

    // s, x = D(x)
    xscaled = _mm_srli_epi32(x, rANS::WordScaleBits);
    freq = _mm_and_si128(freq_bias, _mm_set1_epi32(0xffff));
    bias = _mm_srli_epi32(freq_bias, 16);
    r.simd_ = _mm_add_epi32(_mm_mullo_epi32(xscaled, freq), bias);
    return s;
}

// Renormalize after decoding a symbol.
static inline void simdDecRenorm(rANS::SimdState& r, const u16*& ptr)
{
    static alignas(16) const int8_t shuffles[16][16] = {
#define _ -1 // for readability
        { _,_,_,_, _,_,_,_, _,_,_,_, _,_,_,_ }, // 0000
        { 0,1,_,_, _,_,_,_, _,_,_,_, _,_,_,_ }, // 0001
        { _,_,_,_, 0,1,_,_, _,_,_,_, _,_,_,_ }, // 0010
        { 0,1,_,_, 2,3,_,_, _,_,_,_, _,_,_,_ }, // 0011
        { _,_,_,_, _,_,_,_, 0,1,_,_, _,_,_,_ }, // 0100
        { 0,1,_,_, _,_,_,_, 2,3,_,_, _,_,_,_ }, // 0101
        { _,_,_,_, 0,1,_,_, 2,3,_,_, _,_,_,_ }, // 0110
        { 0,1,_,_, 2,3,_,_, 4,5,_,_, _,_,_,_ }, // 0111
        { _,_,_,_, _,_,_,_, _,_,_,_, 0,1,_,_ }, // 1000
        { 0,1,_,_, _,_,_,_, _,_,_,_, 2,3,_,_ }, // 1001
        { _,_,_,_, 0,1,_,_, _,_,_,_, 2,3,_,_ }, // 1010
        { 0,1,_,_, 2,3,_,_, _,_,_,_, 4,5,_,_ }, // 1011
        { _,_,_,_, _,_,_,_, 0,1,_,_, 2,3,_,_ }, // 1100
        { 0,1,_,_, _,_,_,_, 2,3,_,_, 4,5,_,_ }, // 1101
        { _,_,_,_, 0,1,_,_, 2,3,_,_, 4,5,_,_ }, // 1110
        { 0,1,_,_, 2,3,_,_, 4,5,_,_, 6,7,_,_ }, // 1111
#undef _
    };
    static uint8_t const numbits[16] = {
        0,1,1,2, 1,2,2,3, 1,2,2,3, 2,3,3,4
    };

    __m128i x = r.simd_;

    // NOTE: SSE2+ only offer a signed 32-bit integer compare, while we
    // need unsigned. So we subtract 0x80000000 before the compare,
    // which converts unsigned integers to signed integers in an
    // order-preserving manner.
    __m128i x_biased = _mm_xor_si128(x, _mm_set1_epi32((s32) 0x80000000));
    __m128i greater = _mm_cmpgt_epi32(_mm_set1_epi32(rANS::WordLowBounds - 0x80000000), x_biased);
    unsigned int mask = _mm_movemask_ps(_mm_castsi128_ps(greater));

    // NOTE: this will read slightly past the end of the input buffer.
    // In practice, either pad the input buffer by 8 bytes at the end,
    // or switch to the non-SIMD version once you get close to the end.
    __m128i memvals = _mm_loadl_epi64((const __m128i*)ptr);
    __m128i xshifted = _mm_slli_epi32(x, 16);
    __m128i shufmask = _mm_load_si128((const __m128i*)shuffles[mask]);
    __m128i newx = _mm_or_si128(xshifted, _mm_shuffle_epi8(memvals, shufmask));
    r.simd_ = _mm_blendv_epi8(x, newx, greater);
    ptr += numbits[mask];
}

} // namespace

u64 rANS::calc_encoded_size(u32 size)
{
    return static_cast<u64>(size) * 2 + sizeof(u32) * 258;
}

u32 rANS::encode(u32 dst_size, u8* dst, u32 src_size, const u8* src)
{
    assert(0 < dst_size);
    assert(nullptr != dst);
    assert(0 < src_size);
    assert(nullptr != src);

    u32 freqs[256];
    count(freqs, src_size, src);
    u32 cum_freqs[257];
    cumulative(cum_freqs, freqs);
    normalize(freqs, cum_freqs, ProbScale);
    EncSymbol symbols[256];
    for(u32 i = 0; i < 256; ++i) {
        init(symbols[i], cum_freqs[i], freqs[i], ProbBits);
    }
    State rans;
    init(rans);
    u8* ptr = dst + dst_size;
    for(u32 i = src_size; 0 < i; --i) {
        u8 s = src[i - 1];
        put(rans, ptr, symbols[s]);
    }
    flush(ptr, rans);
    ptr -= sizeof(u32) * 258;
    if(ptr < dst) {
        return 0;
    }
    u32* u32ptr = reinterpret_cast<u32*>(ptr);
    ::memcpy(u32ptr+1, cum_freqs, sizeof(u32)*257);
    ::memcpy(u32ptr, &src_size, sizeof(u32));
    u32 encoded_size = static_cast<u32>(dst + dst_size - ptr);
    return encoded_size;
}

u32 rANS::decode(u32 dst_size, u8* dst, u32 src_size, const u8* src)
{
    assert(0 < dst_size);
    assert(nullptr != dst);
    assert(0 < src_size);
    assert(nullptr != src);
    assert(258 * sizeof(u32) <= src_size);
    u32 original_size;
    ::memcpy(&original_size, src, sizeof(u32));
    if(dst_size<original_size){
        return 0;
    }
    const u32* cum_freqs = reinterpret_cast<const u32*>(src) + 1;
    u8 cum2sym[ProbScale] = {};
    for(u32 s = 0; s < 256; ++s) {
        for(u32 i = cum_freqs[s]; i < cum_freqs[s + 1]; ++i) {
            cum2sym[i] = s;
        }
    }

    const u8* ptr = src + sizeof(u32) * 258;
    State rans;
    init_decode(rans, ptr);

    for(u32 i = 0; i < original_size; ++i) {
        u8 s = cum2sym[get(rans, ProbBits)];
        dst[i] = s;
        u32 freq = cum_freqs[s+1] - cum_freqs[s];
        advance(rans, ptr, cum_freqs[s], freq, ProbBits);
    }
    u32 decoded_size = static_cast<u32>(ptr - (src+sizeof(u32) * 258));
    return decoded_size;
}


u32 rANS::encode_simd(u32 dst_size, u8* dst, u32 src_size, const u8* src)
{
    assert(0 < dst_size);
    assert(nullptr != dst);
    assert(0 < src_size);
    assert(nullptr != src);

    u32 freqs[256];
    count(freqs, src_size, src);
    u32 cum_freqs[257];
    cumulative(cum_freqs, freqs);
    normalize(freqs, cum_freqs, rANS::WordM);

    WordTables tables;
    for(u32 s=0; s<256; ++s){
        initSymbols(tables, static_cast<u8>(s), cum_freqs[s], freqs[s]);
    }

    State rans[8];
    for(u32 i=0; i<8; ++i){
        rans[i] = wordEncInit();
    }

    u16* ptr = reinterpret_cast<u16*>(dst + dst_size);
    for(u32 i = src_size; 0 < i; --i) {
        u8 s = src[i - 1];
        wordEncPut(rans[(i - 1) & 7], ptr, cum_freqs[s], freqs[s]);
    }
    for(u32 i=8; 0<i; --i){
        wordEncFlush(rans[i-1], ptr);
    }
    u8* header = reinterpret_cast<u8*>(ptr) - sizeof(u32) * 258;
    if(header < dst) {
        return 0;
    }
    u32* u32ptr = reinterpret_cast<u32*>(header);
    ::memcpy(u32ptr+1, cum_freqs, sizeof(u32)*257);
    ::memcpy(u32ptr, &src_size, sizeof(u32));
    u32 encoded_size = static_cast<u32>(dst + dst_size - header);
    return encoded_size;
}

u32 rANS::decode_simd(u32 dst_size, u8* dst, u32 src_size, const u8* src)
{
    assert(0 < dst_size);
    assert(nullptr != dst);
    assert(0 < src_size);
    assert(nullptr != src);
    assert(258 * sizeof(u32) <= src_size);
    u32 original_size;
    ::memcpy(&original_size, src, sizeof(u32));
    if(dst_size<original_size){
        return 0;
    }
    const u32* cum_freqs = reinterpret_cast<const u32*>(src) + 1;

    WordTables tables;
    for(u32 s=0; s<256; ++s){
        initSymbols(tables, static_cast<u8>(s), cum_freqs[s], cum_freqs[s+1]-cum_freqs[s]);
    }

    const u16* ptr = reinterpret_cast<const u16*>(src + sizeof(u32) * 258);
    SimdState rans0,rans1;
    simdDecInit(rans0, ptr);
    simdDecInit(rans1, ptr);

    u32 simdSize = original_size & ~0x07UL;

    for(u32 i = 0; i < simdSize; i += 8) {
        u32 s03 = simdDecSym(rans0, tables);
        u32 s47 = simdDecSym(rans1, tables);
        *(u32*)(dst + i) = s03;
        *(u32*)(dst + i + 4) = s47;
        simdDecRenorm(rans0, ptr);
        simdDecRenorm(rans1, ptr);
    }
    for(u32 i = simdSize; i < original_size; ++i) {
        SimdState* which = (i & 4) != 0 ? &rans1 : &rans0;
        uint8_t s = wordDecSym(which->lane_[i & 3], tables);
        dst[i] = s;
    }
    return original_size;
}

} // namespace cppans
#endif
#endif INC_CPPANS_H_
