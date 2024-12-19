#ifndef INC_CPPANS_H_
#define INC_CPPANS_H_
/*
USAGE:
Put '#define CPPANS_IMPLEMENTATION' before including this file to create the implementation.
*/
#define CPPANS_IMPLEMENTATION (1)
#include <cstdint>

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
    using State = u32;

    struct EncSymbol
    {
        u32 x_max_; //!< upper bound of pre-normalization interval
        u32 rcp_freq_;
        u32 bias_;
        u16 cmpl_freq_; //!< Complemt of frequency: (!<<scale_bits) - freq
        u16 rcp_shift_;
    };

    static u64 calc_encoded_size(u32 size);
    static u32 encode(u32 dst_size, u8* dst, u32 src_size, const u8* src);
    static u32 decode(u32 dst_size, u8* dst, u32 src_size, const u8* src);

private:
    rANS(const rANS&) = delete;
    rANS& operator=(const rANS&) = delete;
};
} // namespace cppans

#ifdef CPPANS_IMPLEMENTATION
#    include <cassert>
#    include <cstring>
#    include <immintrin.h>

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
    u32* u32ptr = reinterpret_cast<u32*>(ptr) + 1;
    for(u32 i = 0; i < 257; ++i) {
        u32ptr[i] = cum_freqs[i];
    }
    u32ptr[-1] = src_size;
    u32 encoded_size = static_cast<u32>(dst + dst_size - ptr);
    return encoded_size;
}

u32 rANS::decode(u32 dst_size, u8* dst, u32 src_size, const u8* src)
{
    assert(0 < dst_size);
    assert(nullptr != dst);
    assert(0 < src_size);
    assert(nullptr != src);
    assert(257 * sizeof(u32) <= src_size);
    u32 original_size;
    ::memcpy(&original_size, src, sizeof(u32));
    const u32* cum_freqs = reinterpret_cast<const u32*>(src) + 1;
    u8 cum2sym[ProbScale] = {};
    for(u32 s = 0; s < 256; ++s) {
        for(u32 i = cum_freqs[s]; i < cum_freqs[s + 1]; ++i) {
            cum2sym[i] = s;
        }
    }

    const u8* ptr = src + sizeof(u32) * 257;
    State rans;
    init_decode(rans, ptr);

    for(u32 i = 0; i < original_size; ++i) {
        u8 s = cum2sym[get(rans, ProbBits)];
        dst[i] = s;
        u32 freq = cum_freqs[s+1] - cum_freqs[s];
        advance(rans, ptr, cum_freqs[s], freq, ProbBits);
    }
    u32 decoded_size = static_cast<u32>(ptr - (src+sizeof(u32) * 257));
    return original_size;
}
} // namespace cppans
#endif
#endif INC_CPPANS_H_
