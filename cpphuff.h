#ifndef INC_CPPHUFF_H_
#define INC_CPPHUFF_H_
/*
USAGE:
Put '#define CPPHUFF_IMPLEMENTATION' before including this file to create the implementation.
*/
#define CPPHUFF_IMPLEMENTATION (1)
#include <cstdint>
#include <immintrin.h>

namespace cpphuff
{
using s8 = int8_t;
using s16 = int16_t;
using s32 = int32_t;
using s64 = int64_t;

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

class Huff
{
public:
    inline static constexpr u32 SymbolMax = 255;

    struct CompressTables
    {
        u32 counts_[SymbolMax+1];
    };

    u32 compress(u32 dst_size, u8* dst, u32 src_size, const u8* src);
private:
    Huff(const Huff&) = delete;
    Huff& operator=(const Huff&) = delete;
};
} // namespace cpphuff

#ifdef CPPHUFF_IMPLEMENTATION
namespace cpphuff
{
    u32 Huff::compress(u32 dst_size, u8* dst, u32 src_size, const u8* src)
{
}
} // namespace cpphuff
#    endif
