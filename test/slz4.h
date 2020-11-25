#ifndef INC_SLZ4_H_
#define INC_SLZ4_H_
/**
@file slz4.h
@author t-sakai
@date 2018/11/21 create

# License
This software is distributed under two licenses, choose whichever you like.

## MIT License
Copyright (c) 2020 Takuro Sakai

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Public Domain
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org>

USAGE:
Put '#define SLZ4_IMPLEMENTATION' before including this file to create the implementation.
*/
#ifdef __cplusplus
#    include <cassert>
#    include <cstddef>
#    include <cstdint>
#else
#    include <assert.h>
#    include <stdbool.h>
#    include <stddef.h>
#    include <stdint.h>
#endif

namespace slz4
{
#ifndef SLZ4_NULL
#    ifdef __cplusplus
#        if 201103L <= __cplusplus || 1700 <= _MSC_VER
#            define SLZ4_NULL nullptr
#        else
#            define SLZ4_NULL 0
#        endif
#    else
#        define SLZ4_NULL (void*)0
#    endif
#endif

typedef int8_t s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float f32;
typedef double f64;

typedef ::size_t size_t;
typedef ::ptrdiff_t ptrdiff_t;

#ifdef _NDEBUG
#    define SLZ4_ASSERT(exp)
#else
#    define SLZ4_ASSERT(exp) assert(exp)
#endif

inline s32 minimum(s32 x0, s32 x1)
{
    return x0 < x1 ? x0 : x1;
}

inline u16 minimum(u16 x0, u16 x1)
{
    return x0 < x1 ? x0 : x1;
}

inline const u8* minimum(const u8* x0, const u8* x1)
{
    return x0 < x1 ? x0 : x1;
}

inline const u8* maximum(const u8* x0, const u8* x1)
{
    return x0 < x1 ? x1 : x0;
}

static const u32 MAX_BLOCK_SIZE = 0x7E000000;
static const u16 MAX_DISTANCE = 65535;
static const u32 MIN_MATCH_LENGTH = 4;
static const u32 END_LITERALS = 5;
static const u32 LOOKAHEAD_LENGTH = 12;
static const u32 MAX_IN_TOKEN = 15;

static const u32 DICTIONARY_SIZE = 0x10000U / 4;

struct SLZ4Context
{
    s32 entries_[DICTIONARY_SIZE];
};

struct LZSSMatch
{
    s32 distance_;
    s32 length_;
};

/**
    @return The number of bytes written into 'dst' or a negative value if this function fails.
    @param capacity ... Byte size of a destination buffer 'dst'. This is at least more than the size which is reterned by the compressBound.
    @param dst ... A destination buffer.
    @param size ... Byte size of a source data. The maximum supported value for this is MAX_BLOCK_SIZE.
    @param src ... A source data.
    */
s32 compress(SLZ4Context& context, u32 capacity, u8* dst, u32 size, const u8* src);

/**
    @return The maximum number of bytes that compression may write into a destination buffer.
    @param size ... Byte size of a source data. The maximum supported value for this is MAX_BLOCK_SIZE.
    */
s32 compressBound(u32 size);

/**
    @return The number of bytes written into 'dst' or a negative value if this function fails.
    @param capacity ... Byte size of a destination buffer 'dst'. The maximum supported value for this is MAX_BLOCK_SIZE.
    @param dst ... A destination buffer.
    @param size ... Byte size of a compressed data 'src'.
    @param src ... A source compressed data.
    */
s32 decompress(u32 capacity, u8* dst, u32 size, const u8* src);

} //namespace slz4
#endif //INC_SLZ4_H_

#ifdef SLZ4_IMPLEMENTATION
#include <cstring>
#include <immintrin.h>

#define XXH_INLINE_ALL
#define XXH_STATIC_LINKING_ONLY   /* access advanced declarations */
#define XXH_IMPLEMENTATION   /* access definitions */
#include <xxhash.h>

namespace slz4
{
namespace
{
    u32 pack(const u8* ptr)
    {
        u32 result;
        u8* u = reinterpret_cast<u8*>(&result);
        u[0] = ptr[0];
        u[1] = ptr[1];
        u[2] = ptr[2];
        u[3] = ptr[3];
        return result;
        //return ptr[0] | (ptr[1] << 8) | (ptr[2] << 16) | (ptr[3] << 24);
    }

    u32 push(u32 code, u32 u)
    {
        return (code >> 8) | (u<<24);
    }

    u32 hash(u32 code)
    {
#if 1
        return XXH32_u32(code, 0x811C9DC5U);
#else
        u32 hash = 0x811C9DC5U;

        hash ^= ((code >> 0) & 0xFFU);
        hash *= 0x01000193U;

        hash ^= ((code >> 8) & 0xFFU);
        hash *= 0x01000193U;

        hash ^= ((code >> 16) & 0xFFU);
        hash *= 0x01000193U;

        hash ^= ((code >> 24) & 0xFFU);
        hash *= 0x01000193U;

        return hash;
#endif
    }

    //-------------------------------------------------------------------
    LZSSMatch findLongestMatch(SLZ4Context& context, u32 code, const u8* start, s32 position, s32 end)
    {
        SLZ4_ASSERT(MIN_MATCH_LENGTH <= (end - position));

        u32 index = hash(code) & (DICTIONARY_SIZE - 1);
        s32 result = context.entries_[index];
        context.entries_[index] = position;
        if(result < 0) {
            return {0, 0};
        }
        if(code != pack(start + result)
           || MAX_DISTANCE < (position - result)) {
            return {0, 0};
        }

        SLZ4_ASSERT(result < position);
        s32 distance = position - result;
        s32 match = result + 4;
        position += 4;
        SLZ4_ASSERT(start[match - 4] == start[position - 4]);
        SLZ4_ASSERT(start[match - 3] == start[position - 3]);
        SLZ4_ASSERT(start[match - 2] == start[position - 2]);
        SLZ4_ASSERT(start[match - 1] == start[position - 1]);

        while(position < end) {
            if(*reinterpret_cast<const u32*>(&start[match]) != *reinterpret_cast<const u32*>(&start[position])) {
                while(position < end) {
                    if(start[match] != start[position]) {
                        break;
                    }
                    ++match;
                    ++position;
                }
                break;
            }
            match += 4;
            position += 4;
        }
        return {distance, match - result};
    }

    //-------------------------------------------------------------------
    bool writeLength(u32& pos, u32 capacity, u8* dst, u32 length)
    {
        SLZ4_ASSERT(0 <= length);
        if(length < MAX_IN_TOKEN) {
            return true;
        }
        length -= MAX_IN_TOKEN;

        while(255 <= length) {
            if(capacity <= pos) {
                return false;
            }
            dst[pos++] = 255;
            length -= 255;
        }
        if(capacity <= pos) {
            return false;
        }
        dst[pos++] = static_cast<u8>(length);
        return true;
    }

    //-------------------------------------------------------------------
    bool writeLiterals(u32& pos, u32 capacity, u8* dst, const u8* src, u32 literalLength)
    {
        SLZ4_ASSERT(0 <= literalLength);
        if(capacity <= pos) {
            return false;
        }

        //Write a sequence header
        u32 firstLiteralBits = literalLength < MAX_IN_TOKEN ? literalLength : MAX_IN_TOKEN;

        dst[pos++] = static_cast<u8>(firstLiteralBits << 4);

        //Write literal length
        if(!writeLength(pos, capacity, dst, literalLength)) {
            return false;
        }

        //Write literals
        if(capacity <= (pos + literalLength)) {
            return false;
        }
        for(u32 i = 0; i < literalLength; ++i, ++src) {
            dst[pos++] = *src;
        }
        return true;
    }

    //-------------------------------------------------------------------
    bool writeLiterals(u32& pos, u32 capacity, u8* dst, const u8* src, u32 literalLength, const LZSSMatch& match)
    {
        SLZ4_ASSERT(0 <= literalLength);
        SLZ4_ASSERT(0 <= match.distance_ && match.distance_ <= MAX_DISTANCE);
        SLZ4_ASSERT(MIN_MATCH_LENGTH <= match.length_);
        if(capacity <= pos) {
            return false;
        }

        u32 matchLength = match.length_ - MIN_MATCH_LENGTH;

        //Write a sequence header
        u32 firstLiteralBits = literalLength < MAX_IN_TOKEN ? literalLength : MAX_IN_TOKEN;
        u32 firstMatchBits = matchLength < MAX_IN_TOKEN ? matchLength : MAX_IN_TOKEN;

        dst[pos++] = static_cast<u8>((firstLiteralBits << 4) | firstMatchBits);

        //Write literal length
        if(!writeLength(pos, capacity, dst, literalLength)) {
            return false;
        }

        //Write literals
        if(capacity <= (pos + literalLength)) {
            return false;
        }
        for(u32 i = 0; i < literalLength; ++i, ++src) {
            dst[pos++] = *src;
        }

        //Write match offset
        if(capacity <= (pos + 1)) {
            return false;
        }
        dst[pos++] = static_cast<u8>(match.distance_ & 0xFFU);
        dst[pos++] = static_cast<u8>((match.distance_ >> 8) & 0xFFU);
        if(!writeLength(pos, capacity, dst, matchLength)) {
            return false;
        }
        return true;
    }

    //-------------------------------------------------------------------
    u32 decodeLength(const u8*& current, const u8* end)
    {
        u32 length = 0;
        while(current < end) {
            s32 len = current[0];
            length += len;
            ++current;
            if(len < 255) {
                break;
            }
        }
        return length;
    }

    //-------------------------------------------------------------------
    void copy(u8* dst, const u8* src, u32 size)
    {
        u32 s = size>>6;
        for(u32 i=0; i<s; ++i){
            _mm_storeu_ps(reinterpret_cast<f32*>(dst), _mm_loadu_ps(reinterpret_cast<const f32*>(src)));
            _mm_storeu_ps(reinterpret_cast<f32*>(dst+16), _mm_loadu_ps(reinterpret_cast<const f32*>(src+16)));
            _mm_storeu_ps(reinterpret_cast<f32*>(dst+32), _mm_loadu_ps(reinterpret_cast<const f32*>(src+32)));
            _mm_storeu_ps(reinterpret_cast<f32*>(dst+48), _mm_loadu_ps(reinterpret_cast<const f32*>(src+48)));
            dst += 64;
            src += 64;
        }
        size -= (s<<6);
        s = size>>4;
        for(u32 i=0; i<s; ++i){
            _mm_storeu_ps(reinterpret_cast<f32*>(dst), _mm_loadu_ps(reinterpret_cast<const f32*>(src)));
            dst += 16;
            src += 16;
        }

        s = size - (s<<4);
        for(u32 i=0; i<s; ++i){
            dst[i] = src[i];
        }
    }
} // namespace

//-------------------------------------------------------------------
s32 compress(SLZ4Context& context, u32 capacity, u8* dst, u32 size, const u8* src)
{
    SLZ4_ASSERT(0 <= size && size <= MAX_BLOCK_SIZE);
    if(MAX_BLOCK_SIZE < size) {
        return -1;
    }

    u32 dstPos = 0;
    if(size < END_LITERALS) {
        return writeLiterals(dstPos, capacity, dst, src, size)
                   ? dstPos
                   : -1;
    }

    memset(context.entries_, -1, sizeof(s32) * DICTIONARY_SIZE);
    { //Add the first code to our dictionary
        u32 code = pack(src);
        u32 index = hash(code) & (DICTIONARY_SIZE - 1);
        context.entries_[index] = 0;
    }

    u32 end = size;
    u32 endMatch = (LOOKAHEAD_LENGTH < size) ? size - LOOKAHEAD_LENGTH : 0;
    u32 pending = 0;
    u32 position = 4;
    while(position < endMatch) {
        u32 code = pack(src+position);
        LZSSMatch match = findLongestMatch(context, code, src, position, end);

        if(MIN_MATCH_LENGTH <= match.length_) {
            u32 literalLength = position - pending;
            if(!writeLiterals(dstPos, capacity, dst, src + pending, literalLength, match)) {
                return -1;
            }
            position += match.length_;
            pending = position;
            continue;
        }
        ++position;
    }

    if(pending < end) {
        u32 literalLength = end - pending;
        if(!writeLiterals(dstPos, capacity, dst, src + pending, literalLength)) {
            return -1;
        }
    }
    return dstPos;
}

//-------------------------------------------------------------------
s32 compressBound(u32 size)
{
    SLZ4_ASSERT(0 <= size && size <= MAX_BLOCK_SIZE);
    return (size <= MAX_BLOCK_SIZE) ? size + (size / 255) + 16 : -1;
}

//-------------------------------------------------------------------
s32 decompress(u32 capacity, u8* dst, u32 size, const u8* src)
{
    SLZ4_ASSERT(0 <= size);
    SLZ4_ASSERT(0 <= capacity && capacity <= MAX_BLOCK_SIZE);
    if(MAX_BLOCK_SIZE<capacity){
        return -1;
    }
    const u8* current = src;
    const u8* end0 = src + size;
    const u8* end1 = end0 - END_LITERALS;
    const u8* dend = dst + capacity;
    u8* d = dst;

    for(;;) {
        //Decode token
        if(end0 <= current) {
            return -1;
        }
        u32 literalLength = (current[0] >> 4);
        u32 matchLength = (current[0]) & 0xFU;
        ++current;

        //Decode literal length
        if(MAX_IN_TOKEN <= literalLength) {
            literalLength += decodeLength(current, end0);
        }

        //Read literals
        if(end0<(current + literalLength)){
            return -1;
        }
        if(dend<(d + literalLength)){
            return -1;
        }
        //memcpy(d, current, literalLength);
        copy(d, current, literalLength);
        d += literalLength;
        current += literalLength;

        if(end1 <= current) {
            break;
        }

        //Decode match offset
        s32 offset = current[0] | (current[1] << 8);
        current += 2;

        //Decode match length
        if(MAX_IN_TOKEN <= matchLength) {
            matchLength += decodeLength(current, end0);
        }

        //Copy match
        if(static_cast<s32>(d-dst) < offset) {
            return -1;
        }
        matchLength += MIN_MATCH_LENGTH;
        if(dend<(d + matchLength)){
            return -1;
        }
        if(16 <= offset) {
            copy(d, d-offset, matchLength);
            d += matchLength;
        } else {
            while(0 < matchLength) {
                d[0] = d[-offset];
                ++d;
                --matchLength;
            }
        }
    }
    return static_cast<s32>(d-dst);
}
} //namespace slz4
#endif //SLZ4_IMPLEMENTATION
