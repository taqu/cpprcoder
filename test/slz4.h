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
#include <cstdint>
#include <cassert>
#include <cstddef>
#else
#include <stdint.h>
#include <assert.h>
#include <stddef.h>
#include <stdbool.h>
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
#define SLZ4_ASSERT(exp)
#else
#define SLZ4_ASSERT(exp) assert(exp)
#endif

    inline s32 minimum(s32 x0, s32 x1)
    {
        return x0<x1? x0 : x1;
    }

    inline u16 minimum(u16 x0, u16 x1)
    {
        return x0<x1? x0 : x1;
    }

    inline const u8* minimum(const u8* x0, const u8* x1)
    {
        return x0<x1? x0 : x1;
    }

    inline const u8* maximum(const u8* x0, const u8* x1)
    {
        return x0<x1? x1 : x0;
    }

    static const s32 MAX_BLOCK_SIZE = 0x7E000000;
    static const u16 MAX_DISTANCE = 65535;
    static const s32 MIN_MATCH_LENGTH = 4;
    static const s32 LAZY_MATCH_LENGTH = 8;
    static const s32 END_LITERALS = 5;
    static const s32 LOOKAHEAD_LENGTH = 12;
    static const s32 MAX_IN_TOKEN = 15;

    struct LZSSMatch
    {
        u16 distance_;
        s32 length_;
    };

    /**
    @return The number of bytes written into 'dst' or a negative value if this function fails.
    @param capacity ... Byte size of a destination buffer 'dst'. This is at least more than the size which is reterned by the compressBound.
    @param dst ... A destination buffer.
    @param size ... Byte size of a source data. The maximum supported value for this is MAX_BLOCK_SIZE.
    @param src ... A source data.
    */
    s32 compress(s32 capacity, u8* dst, s32 size, const u8* src);

    /**
    @return The maximum number of bytes that compression may write into a destination buffer.
    @param size ... Byte size of a source data. The maximum supported value for this is MAX_BLOCK_SIZE.
    */
    s32 compressBound(s32 size);

    /**
    @return The number of bytes written into 'dst' or a negative value if this function fails.
    @param capacity ... Byte size of a destination buffer 'dst'. The maximum supported value for this is MAX_BLOCK_SIZE.
    @param dst ... A destination buffer.
    @param size ... Byte size of a compressed data 'src'.
    @param src ... A source compressed data.
    */
    s32 decompress(s32 capacity, u8* dst, s32 size, const u8* src);

} //namespace slz4
#endif //INC_SLZ4_H_

#ifdef SLZ4_IMPLEMENTATION
namespace slz4
{
namespace
{
    u32 pack(const u8* ptr)
    {
        return ptr[0] | (ptr[1]<<8) | (ptr[2]<<16) | (ptr[3]<<24);
    }

    u32 push(u32 code, const u8* ptr)
    {
        return (code <<8) | ptr[0];
    }

    //-------------------------------------------------------------------
    LZSSMatch findLongestMatch(const u8* current, const u8* start, const u8* end)
    {
#if !defined(_NDEBUG)
        SLZ4_ASSERT(static_cast<s32>(end-start)<=MAX_BLOCK_SIZE);
        SLZ4_ASSERT(MIN_MATCH_LENGTH<=static_cast<s32>(end-current));
#endif
        u32 curCode = pack(current);
        u32 preCode = curCode;
        s32 length = 0;
        const u8* match = current;

        for(const u8* pos = current-1; start <= pos; --pos){
            preCode = push(preCode, pos);
            if(preCode == curCode){
                const u8* s = pos + 4;
                const u8* c = current + 4;
                for(; c<end; ++s,++c){
                    if(*c != *s){
                        break;
                    }
                }
                s32 l = static_cast<s32>(s-pos);
                if(length<l){
                    length = l;
                    match = pos;
                }
            }
        }
        return {static_cast<u16>(current-match), length};
    }

    //-------------------------------------------------------------------
    bool writeLength(s32& pos, s32 capacity, u8* dst, s32 length)
    {
        SLZ4_ASSERT(0<=length);
        if(length<MAX_IN_TOKEN){
            return true;
        }
        length -= MAX_IN_TOKEN;

        while(255<=length){
            if(capacity<=pos){
                return false;
            }
            dst[pos++] = 255;
            length -= 255;
        }
        if(capacity<=pos){
            return false;
        }
        dst[pos++] = static_cast<u8>(length);
        return true;
    }

    //-------------------------------------------------------------------
    bool writeLiterals(s32& pos, s32 capacity, u8* dst, const u8* src, s32 literalLength, const LZSSMatch& match)
    {
        SLZ4_ASSERT(0<=literalLength);
        if(capacity<=pos){
            return false;
        }

        s32 matchLength = (MIN_MATCH_LENGTH<=match.length_)? (match.length_-MIN_MATCH_LENGTH) : 0;

        //Write a sequence header
        s32 firstLiteralBits = literalLength<MAX_IN_TOKEN? literalLength : MAX_IN_TOKEN;
        s32 firstMatchBits = matchLength<MAX_IN_TOKEN? matchLength : MAX_IN_TOKEN;

        dst[pos++] = static_cast<u8>((firstLiteralBits<<4) | firstMatchBits);

        //Write literal length
        if(!writeLength(pos, capacity, dst, literalLength)){
            return false;
        }

        //Write literals
        if(capacity<=(pos+literalLength)){
            return false;
        }
        for(s32 i = 0; i<literalLength; ++i, ++src){
            dst[pos++] = *src;
        }

        //Write match offset
        if(match.length_<=0){
            return true;
        }
        if(capacity<=(pos+1)){
            return false;
        }
        dst[pos++] = static_cast<u8>(match.distance_ & 0xFFU);
        dst[pos++] = static_cast<u8>((match.distance_>>8) & 0xFFU);
        if(!writeLength(pos, capacity, dst, matchLength)){
            return false;
        }
        return true;
    }

    //-------------------------------------------------------------------
    void decodeLength(s32& length, const u8*& current, const u8* end)
    {
        if(length<MAX_IN_TOKEN){
            return;
        }
        while(current<end){
            s32 len = current[0];
            length += len;
            ++current;
            if(len<255){
                break;
            }
        }
    }
}

    //-------------------------------------------------------------------
    s32 compress(s32 capacity, u8* dst, s32 size, const u8* src)
    {
        SLZ4_ASSERT(0<=size && size<=MAX_BLOCK_SIZE);
        if(MAX_BLOCK_SIZE<size){
            return -1;
        }

        s32 dstPos = 0;

        const u8* current = src;
        const u8* start = current;
        const u8* end = current + size;
        const u8* endMatch = (LOOKAHEAD_LENGTH<size)? current + size - LOOKAHEAD_LENGTH : current;
        const u8* endLiterals = (END_LITERALS<=size)? current + size - END_LITERALS : current;
        const u8* pending = current;
        while(current<endMatch){
            start = maximum(current-MAX_DISTANCE, src);
            const u8* ahead = start + MAX_BLOCK_SIZE;
            const u8* end0 = ahead<endLiterals? ahead : endLiterals;
            LZSSMatch match = findLongestMatch(current, start, end0);

            if(MIN_MATCH_LENGTH<=match.length_){
                s32 literalLength = static_cast<s32>(current-pending);
                if(!writeLiterals(dstPos, capacity, dst, pending, literalLength, match)){
                    return -1;
                }
                current += match.length_;
                pending = current;
                continue;
            }
            ++current;
        }

        if(pending<end){
            s32 literalLength = static_cast<s32>(end-pending);
            if(!writeLiterals(dstPos, capacity, dst, pending, literalLength, {0,})){
                return -1;
            }
            pending = SLZ4_NULL;
        }
        return dstPos;
    }

    //-------------------------------------------------------------------
    s32 compressBound(s32 size)
    {
        SLZ4_ASSERT(0<=size && size<=MAX_BLOCK_SIZE);
        return (size<=MAX_BLOCK_SIZE)? size + (size/255) + 16 : -1;
    }

    //-------------------------------------------------------------------
    s32 decompress(s32 capacity, u8* dst, s32 size, const u8* src)
    {
        SLZ4_ASSERT(0<=size);
        SLZ4_ASSERT(0<=capacity && capacity<=MAX_BLOCK_SIZE);
        s32 dstPos = 0;

        const u8* current = src;
        const u8* end0 = src + size;
        const u8* end1 = end0 - END_LITERALS;

        while(current<end0 && dstPos<capacity){
            //Decode token
            s32 literalLength = (current[0]>>4);
            s32 matchLength = (current[0]) & 0xFU;
            ++current;

            //Decode literal length
            decodeLength(literalLength, current, end0);

            //Read literals
            const u8* litEnd = minimum(end0, current+literalLength);
            while(current<litEnd && dstPos<capacity){
                dst[dstPos++] = current[0];
                ++current;
            }

            if(end1<=current){
                break;
            }

            //Decode match offset
            s32 offset = current[0] | (current[1]<<8);
            current += 2;

            //Decode match length
            decodeLength(matchLength, current, end0);

            //Copy match
            if(dstPos<offset){
                return -1;
            }
            matchLength += MIN_MATCH_LENGTH;
            while(0<matchLength && dstPos<capacity){
                dst[dstPos++] = dst[dstPos-offset];
                --matchLength;
            }
        }
        return dstPos;
    }
} //namespace slz4
#endif //SLZ4_IMPLEMENTATION
