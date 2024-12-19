#ifndef INC_CPPRCODER_H_
#define INC_CPPRCODER_H_
/**
@author t-sakai
@date 2018/09/04 create

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
Put '#define CPPRCODER_IMPLEMENTATION' before including this file to create the implementation.

References:
[1] http://www.nct9.ne.jp/m_hiroi/light/pyalgo37.html
*/
#include <cassert>
#include <cstdint>
#include <cstring>

#define CPPRCODER_USE_SIMD

#ifdef CPPRCODER_USE_SIMD
#    include <immintrin.h>
#endif

#define CPPRCODER_CREATE_MEMORYSTREAM

#ifndef CPPRCODER_ASSERT
#    define CPPRCODER_ASSERT(exp) assert(exp)
#endif

namespace cpprcoder
{
#ifndef CPPRCODER_NULL
#    ifdef __cplusplus
#        if 201103L <= __cplusplus || 1700 <= _MSC_VER
#            define CPPRCODER_NULL nullptr
#        else
#            define CPPRCODER_NULL 0
#        endif
#    else
#        define CPPRCODER_NULL (void*)0
#    endif
#endif

#ifndef CPPRCODER_TYPES
#    define CPPRCODER_TYPES
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

typedef char Char;

using ::size_t;
using ::uintptr_t;
#endif

#ifndef CPPRCODER_OFF_T
#    define CPPRCODER_OFF_T
typedef s64 off_t;
#endif

inline u32 maximum(u32 x0, u32 x1)
{
    return x0 < x1 ? x1 : x0;
}

enum Status
{
    Status_Success = 0,
    Status_Pending = 1,
    Status_Error = -1,
};

struct Result
{
    Status status_;
    u32 requestSize_;
};

//----------------------------------------------
//---
//--- IStream
//---
//----------------------------------------------
template<class T>
class IStream
{
public:
    s32 read(s32 size, u8* bytes)
    {
        return T::read(size, bytes);
    }
    bool readByte(u8& byte)
    {
        return T::readByte(byte);
    }

    s32 write(s32 size, const u8* bytes)
    {
        return T::write(size, bytes);
    }

    bool writeByte(u8 byte)
    {
        return T::writeByte(byte);
    }

protected:
    IStream()
    {
    }
    ~IStream()
    {
    }

private:
    IStream(const IStream&) = delete;
    IStream(IStream&&) = delete;
    IStream& operator=(const IStream&) = delete;
    IStream& operator=(IStream&&) = delete;
};

#ifdef CPPRCODER_CREATE_MEMORYSTREAM

#    ifndef CPPRCODER_MALLOC
#        define CPPRCODER_MALLOC(size) malloc(size)
#    endif

#    ifndef CPPRCODER_FREE
#        define CPPRCODER_FREE(ptr) \
            free((ptr)); \
            (ptr) = CPPRCODER_NULL
#    endif

//----------------------------------------------
//---
//--- MemoryStream
//---
//----------------------------------------------
class MemoryStream: public IStream<MemoryStream>
{
public:
    MemoryStream();
    explicit MemoryStream(s32 capacity);
    ~MemoryStream();

    inline s32 capacity() const;
    inline s32 size() const;

    inline const u8* get() const;
    inline const u8& operator[](s32 index) const;
    inline u8& operator[](s32 index);

    void reserve(s32 capacity);
    void resize(s32 size);

    s32 read(s32 size, u8* bytes);
    bool readByte(u8& byte);
    s32 write(s32 size, const u8* bytes);
    bool writeByte(u8 byte);

private:
    MemoryStream(const MemoryStream&) = delete;
    MemoryStream(MemoryStream&&) = delete;
    MemoryStream& operator=(const MemoryStream&) = delete;
    MemoryStream& operator=(MemoryStream&&) = delete;

    static const s32 EXPAND_LIMIT_SIZE = 4096 * 4;

    bool expand(s32 size);

    s32 capacity_;
    s32 size_;
    u8* buffer_;
};

inline s32 MemoryStream::capacity() const
{
    return capacity_;
}

inline s32 MemoryStream::size() const
{
    return size_;
}

inline const u8* MemoryStream::get() const
{
    return buffer_;
}

inline const u8& MemoryStream::operator[](s32 index) const
{
    CPPRCODER_ASSERT(0 <= index && index < size_);
    return buffer_[index];
}

inline u8& MemoryStream::operator[](s32 index)
{
    CPPRCODER_ASSERT(0 <= index && index < size_);
    return buffer_[index];
}

#endif

//----------------------------------------------
//---
//--- AdaptiveFrequencyTable
//---
//----------------------------------------------
class AdaptiveFrequencyTable
{
public:
    static const u32 MAXRANGE = 0xFFFFFF00U;
    static const u32 MINRANGE = 0x1000000U;

    static const u32 CHUNK_SIZE = 16;
    static const u32 SIZE = 256;
    static const u32 CHUNKS = (SIZE / CHUNK_SIZE);
    static const uintptr_t ALIGN = 16;
    static const uintptr_t ALIGN_MASK = ALIGN - 1;

    AdaptiveFrequencyTable();
    ~AdaptiveFrequencyTable();

    void initialize();

    void update(u8 b);
    inline u32 size() const;
    inline u32 total() const;
    inline u32 operator[](u32 index) const;
    u32 cumulative(u8 byte) const;
    void find(u32& count, u8& code, u32 target) const;

private:
    AdaptiveFrequencyTable(const AdaptiveFrequencyTable&) = delete;
    AdaptiveFrequencyTable(AdaptiveFrequencyTable&&) = delete;
    AdaptiveFrequencyTable& operator=(const AdaptiveFrequencyTable&) = delete;
    AdaptiveFrequencyTable& operator=(AdaptiveFrequencyTable&&) = delete;

    void countChunks();

    u32 total_;
#if defined(CPPRCODER_USE_SIMD)
    u32* prefix_;
    u32* frequencies_;
    u32 prefix_buffer_[CHUNKS + 4];
    u32 frequencies_buffer_[SIZE + 4];
#else
    u32 prefix_[CHUNKS];
    u32 frequencies_[SIZE];
#endif
};

inline u32 AdaptiveFrequencyTable::size() const
{
    return SIZE;
}

inline u32 AdaptiveFrequencyTable::total() const
{
    return total_;
}

inline u32 AdaptiveFrequencyTable::operator[](u32 index) const
{
    CPPRCODER_ASSERT(0 <= index && index < SIZE);
    return frequencies_[index];
}

//----------------------------------------------
//---
//--- RangeEncoder
//---
//----------------------------------------------
template<class T = MemoryStream>
class RangeEncoder
{
public:
    static const u32 SHIFT = 24;
    static const u32 MAX_RANGE = 0xFFFFFFFFU;
    static const u32 MIN_RANGE = 0x01U << SHIFT;
    static const u32 MASK = 0x00FFFFFFU;
    static const u32 MAX_SIZE = 0x7FFFFFFFU;
    static const u32 FREQUENCY_SIZE = 256;
    static const u32 HeaderSize16 = sizeof(u32) + sizeof(u16) * FREQUENCY_SIZE;

    RangeEncoder();
    ~RangeEncoder();

    bool encode(T& stream, u32 size, const u8* bytes);
    bool decode(T& stream, u32 size, const u8* bytes);

private:
    RangeEncoder(const RangeEncoder&) = delete;
    RangeEncoder(RangeEncoder&&) = delete;
    RangeEncoder& operator=(const RangeEncoder&) = delete;
    RangeEncoder& operator=(RangeEncoder&&) = delete;

    u8 find(u32 target) const;
    u32 count(u32 index) const;
    /**
        @brief
        @return mode
        */
    void count(u32 size, const u8* bytes);
    void calcCumulatives();
    bool write16(T& stream);
    bool read16(const u8*& bytes, const u8* end);

    u32 cumulatives_[FREQUENCY_SIZE + 1];
    u32 range_;
    u32 low_;
    u32 count_;
    u32 buffer_;
};

template<class T>
RangeEncoder<T>::RangeEncoder()
    : range_(MAX_RANGE)
    , low_(0)
{
}

template<class T>
RangeEncoder<T>::~RangeEncoder()
{
}

template<class T>
bool RangeEncoder<T>::encode(T& stream, u32 size, const u8* bytes)
{
    CPPRCODER_ASSERT(size <= MAX_SIZE);
    CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);

    count(size, bytes);
    range_ = MAX_RANGE;
    low_ = 0;
    count_ = 0;
    buffer_ = 0;
    u8 output[4];
    output[0] = (size >> 0) & 0xFFU;
    output[1] = (size >> 8) & 0xFFU;
    output[2] = (size >> 16) & 0xFFU;
    output[3] = (size >> 24) & 0xFFU;
    if(stream.write(4, output) <= 0) {
        return false;
    }
    //write frequencies
    if(!write16(stream)) {
        return false;
    }
    calcCumulatives();

    for(u32 i = 0; i < size; ++i) {
        u32 t = range_ / cumulatives_[FREQUENCY_SIZE];
        u32 newLow = low_ + cumulatives_[bytes[i]] * t;
        range_ = count(bytes[i]) * t;

        if(newLow < low_) {
            ++buffer_;
            if(0 < count_) {
                for(; 0 != count_; --count_) {
                    if(!stream.writeByte(static_cast<u8>(buffer_))) {
                        return false;
                    }
                    buffer_ = 0;
                }
            }
        }
        low_ = newLow;

        while(range_ < MIN_RANGE) {
            if(low_ < (0xFFU << SHIFT)) {
                if(!stream.writeByte(static_cast<u8>(buffer_))) {
                    return false;
                }
                for(; 0 != count_; --count_) {
                    if(!stream.writeByte(0xFFU)) {
                        return false;
                    }
                }
                buffer_ = low_ >> SHIFT;

            } else {
                ++count_;
            }
            low_ <<= 8;
            range_ <<= 8;
        }
    }

    //
    u8 c = 0xFFU;
    if(MAX_RANGE <= low_) {
        ++buffer_;
        c = 0;
    }
    if(!stream.writeByte(static_cast<u8>(buffer_))) {
        return false;
    }
    for(; 0 != count_; --count_) {
        if(!stream.writeByte(c)) {
            return false;
        }
    }

    output[0] = (low_ >> 24) & 0xFFU;
    output[1] = (low_ >> 16) & 0xFFU;
    output[2] = (low_ >> 8) & 0xFFU;
    output[3] = (low_ >> 0) & 0xFFU;
    return 0 < stream.write(4, output);
}

template<class T>
bool RangeEncoder<T>::decode(T& stream, u32 size, const u8* bytes)
{
    CPPRCODER_ASSERT(size <= MAX_SIZE);
    CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);

    range_ = MAX_RANGE;

    if(size < 1) {
        return false;
    }

    const u8* end = bytes + size;
    u32 uncompressedSize;
    if(size < HeaderSize16) {
        return false;
    }
    uncompressedSize = bytes[3];
    uncompressedSize = (uncompressedSize << 8) | bytes[2];
    uncompressedSize = (uncompressedSize << 8) | bytes[1];
    uncompressedSize = (uncompressedSize << 8) | bytes[0];
    if(uncompressedSize <= 0) {
        return true;
    }
    bytes += 4;
    //read frequencies
    if(!read16(bytes, end)) {
        return false;
    }
    calcCumulatives();

    if(static_cast<u32>(end - bytes) < 5) {
        return false;
    }
    low_ = bytes[1];
    low_ = (low_ << 8) | bytes[2];
    low_ = (low_ << 8) | bytes[3];
    low_ = (low_ << 8) | bytes[4];
    bytes += 5;

    for(u32 i = 0; i < uncompressedSize; ++i) {
        u32 t = range_ / cumulatives_[FREQUENCY_SIZE];
        u8 c = find(low_ / t);
        low_ -= cumulatives_[c] * t;
        range_ = count(c) * t;

        while(range_ < MIN_RANGE) {
            if(end <= bytes) {
                return false;
            }
            range_ <<= 8;
            low_ = (low_ << 8) | bytes[0];
            ++bytes;
        }
        if(!stream.writeByte(c)) {
            return false;
        }
    }
    return true;
}

template<class T>
u8 RangeEncoder<T>::find(u32 target) const
{
    u32 left = 0;
    u32 right = FREQUENCY_SIZE - 1;
    while(left < right) {
        u32 m = (left + right) >> 1;
        if(cumulatives_[m + 1] <= target) {
            left = m + 1;
        } else {
            right = m;
        }
    }
    return static_cast<u8>(left);
}

template<class T>
u32 RangeEncoder<T>::count(u32 index) const
{
    return cumulatives_[index + 1] - cumulatives_[index];
}

template<class T>
void RangeEncoder<T>::count(u32 size, const u8* bytes)
{
    CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);
    memset(cumulatives_, 0, sizeof(u32) * (FREQUENCY_SIZE + 1));
    for(u32 i = 0; i < size; ++i) {
        if(0xFFFFU <= cumulatives_[bytes[i]]) {
            for(u32 j = 0; j < FREQUENCY_SIZE; ++j) {
                if(0 < cumulatives_[j]) {
                    cumulatives_[j] = (cumulatives_[j] >> 1) | 0x01U;
                }
            }
        }

        cumulatives_[bytes[i]] += 1;
    }

    //Sum of frecuencies equals to the size, and should be less than or equals to MIN_RANGE.
    if(MIN_RANGE < size) {
        u32 n = 0;
        while(MIN_RANGE < size) {
            size >>= 1;
            ++n;
        }
        for(u32 i = 1; i < FREQUENCY_SIZE; ++i) {
            cumulatives_[i] = (0 != cumulatives_[i]) ? ((cumulatives_[i] >> n) | 0x01U) : 0;
        }
    }
}

template<class T>
void RangeEncoder<T>::calcCumulatives()
{
    u32 sum = cumulatives_[0];
    cumulatives_[0] = 0;
    for(u32 i = 0; i < FREQUENCY_SIZE; ++i) {
        u32 next = sum + cumulatives_[i + 1];
        cumulatives_[i + 1] = sum;
        sum = next;
    }
}

template<class T>
bool RangeEncoder<T>::read16(const u8*& bytes, const u8* end)
{
    CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);

    static const u32 STEP = 64;
    u16 buffer[STEP];
    for(u32 i = 0; i < FREQUENCY_SIZE; i += STEP) {
        memcpy(buffer, bytes, sizeof(u16) * STEP);
        bytes += sizeof(u16) * STEP;
        for(u32 j = 0; j < STEP; ++j) {
            cumulatives_[i + j] = buffer[j];
        }
    }

    cumulatives_[FREQUENCY_SIZE] = 0;
    return bytes < end;
}

template<class T>
bool RangeEncoder<T>::write16(T& stream)
{
    static const u32 STEP = 64;
    u16 buffer[STEP];
    for(u32 i = 0; i < FREQUENCY_SIZE; i += STEP) {
        for(u32 j = 0; j < STEP; ++j) {
            u32 x = cumulatives_[i + j];
            buffer[j] = static_cast<u16>(x & 0xFFFFU);
        }
        if(stream.write(sizeof(u16) * STEP, reinterpret_cast<const u8*>(buffer)) <= 0) {
            return false;
        }
    }
    return true;
}

//----------------------------------------------
//---
//--- AdaptiveRangeEncoder
//---
//----------------------------------------------
template<class T = MemoryStream>
class AdaptiveRangeEncoder
{
public:
    static const u32 MAXRANGE = 0xFFFFFF00U;
    static const u32 MINRANGE = 0x1000000U;

    AdaptiveRangeEncoder();
    ~AdaptiveRangeEncoder();

    bool initialize(T& stream, u32 umcompressedSize);
    Result encode(s32 size, const u8* bytes);
    Result encode(u8 byte);

private:
    AdaptiveRangeEncoder(const AdaptiveRangeEncoder&) = delete;
    AdaptiveRangeEncoder(AdaptiveRangeEncoder&&) = delete;
    AdaptiveRangeEncoder& operator=(const AdaptiveRangeEncoder&) = delete;
    AdaptiveRangeEncoder& operator=(AdaptiveRangeEncoder&&) = delete;

    static const s32 SHIFT = 24;

    bool finish();
    bool normalize(u32 prevLow);

    T* stream_;
    AdaptiveFrequencyTable frequencies_;
    u32 umcompressedSize_;
    u32 inSize_;
    u32 range_;
    u32 low_;
    u8 buffer_;
    u32 carry_;
};

template<class T>
AdaptiveRangeEncoder<T>::AdaptiveRangeEncoder()
    : stream_(CPPRCODER_NULL)
    , umcompressedSize_(0)
    , inSize_(0)
    , range_(MAXRANGE)
    , low_(0)
    , buffer_(0)
    , carry_(0)
{
}

template<class T>
AdaptiveRangeEncoder<T>::~AdaptiveRangeEncoder()
{
}

template<class T>
bool AdaptiveRangeEncoder<T>::initialize(T& stream, u32 umcompressedSize)
{
    frequencies_.initialize();
    stream_ = &stream;
    range_ = MAXRANGE;
    umcompressedSize_ = umcompressedSize;
    inSize_ = 0;
    low_ = 0;
    buffer_ = 0;
    carry_ = 0;
    u8 bytes[4];
    bytes[0] = (umcompressedSize_ >> 0) & 0xFFU;
    bytes[1] = (umcompressedSize_ >> 8) & 0xFFU;
    bytes[2] = (umcompressedSize_ >> 16) & 0xFFU;
    bytes[3] = (umcompressedSize_ >> 24) & 0xFFU;
    return 0 < stream_->write(4, bytes);
}

template<class T>
Result AdaptiveRangeEncoder<T>::encode(s32 size, const u8* bytes)
{
    CPPRCODER_ASSERT((inSize_ + size) <= umcompressedSize_);

    for(s32 i = 0; i < size; ++i) {
        u32 t = range_ / frequencies_.total();
        u8 c = bytes[i];
        u32 prevLow = low_;
        low_ += frequencies_.cumulative(c) * t;
        range_ = frequencies_[c] * t;
        if(!normalize(prevLow)) {
            inSize_ += i;
            return {Status_Pending, umcompressedSize_ - inSize_};
        }
        frequencies_.update(c);
    }
    inSize_ += size;
    if(umcompressedSize_ <= inSize_) {
        finish();
        return {Status_Success, 0};
    }
    return {Status_Pending, umcompressedSize_ - inSize_};
}

template<class T>
Result AdaptiveRangeEncoder<T>::encode(u8 byte)
{
    CPPRCODER_ASSERT((inSize_ + 1) <= umcompressedSize_);

    u32 t = range_ / frequencies_.total();
    u32 prevLow = low_;
    low_ += frequencies_.cumulative(byte) * t;
    range_ = frequencies_[byte] * t;
    if(!normalize(prevLow)) {
        inSize_ += 1;
        return {Status_Pending, umcompressedSize_ - inSize_};
    }
    frequencies_.update(byte);
    ++inSize_;
    if(umcompressedSize_ <= inSize_) {
        finish();
        return {Status_Success, 0};
    }
    return {Status_Pending, umcompressedSize_ - inSize_};
}

template<class T>
bool AdaptiveRangeEncoder<T>::finish()
{
    if(!stream_->writeByte(buffer_)) {
        return false;
    }
    for(u32 i = 0; i < carry_; ++i) {
        if(!stream_->writeByte(0xFFU)) {
            return false;
        }
    }

    u8 bytes[4];
    bytes[0] = (low_ >> 24) & 0xFFU;
    bytes[1] = (low_ >> 16) & 0xFFU;
    bytes[2] = (low_ >> 8) & 0xFFU;
    bytes[3] = (low_ >> 0) & 0xFFU;
    return 0 < stream_->write(4, bytes);
}

template<class T>
bool AdaptiveRangeEncoder<T>::normalize(u32 prevLow)
{
    if(low_ < prevLow) {
        buffer_ += 1;
        if(0 < carry_) {
            if(!stream_->writeByte(buffer_)) {
                return false;
            }
            for(u32 i = 1; i < carry_; ++i) {
                if(!stream_->writeByte(0)) {
                    return false;
                }
            }
            buffer_ = 0;
            carry_ = 0;
        }
    }

    while(range_ < MINRANGE) {
        if(low_ < (0xFFU << SHIFT)) {
            if(!stream_->writeByte(buffer_)) {
                return false;
            }
            for(u32 i = 0; i < carry_; ++i) {
                if(!stream_->writeByte(0xFFU)) {
                    return false;
                }
            }
            buffer_ = (low_ >> SHIFT) & 0xFFU;
            carry_ = 0;
        } else {
            carry_ += 1;
        }
        low_ <<= 8;
        range_ <<= 8;
    }
    return true;
}

//----------------------------------------------
//---
//--- AdaptiveRangeDecoder
//---
//----------------------------------------------
template<class T = MemoryStream>
class AdaptiveRangeDecoder
{
public:
    static const u32 MAXRANGE = 0x00FFFFFFU;
    static const u32 MINRANGE = 0x1000000U;

    AdaptiveRangeDecoder();
    ~AdaptiveRangeDecoder();

    bool initialize(T& stream);
    Result decode(s32 size, const u8* bytes);

private:
    AdaptiveRangeDecoder(const AdaptiveRangeDecoder&) = delete;
    AdaptiveRangeDecoder(AdaptiveRangeDecoder&&) = delete;
    AdaptiveRangeDecoder& operator=(const AdaptiveRangeDecoder&) = delete;
    AdaptiveRangeDecoder& operator=(AdaptiveRangeDecoder&&) = delete;

    static const s32 State_Init = 0;
    static const s32 State_Decode = 1;

    bool normalize(s32& size, const u8*& bytes);

    T* stream_;
    AdaptiveFrequencyTable frequencies_;
    u32 umcompressedSize_;
    u32 outSize_;
    s32 state_;
    u32 range_;
    u32 low_;
    u8 lastCode_;
};

template<class T>
AdaptiveRangeDecoder<T>::AdaptiveRangeDecoder()
    : stream_(CPPRCODER_NULL)
    , umcompressedSize_(0)
    , outSize_(0)
    , state_(State_Init)
    , range_(MAXRANGE)
    , low_(0)
{
}

template<class T>
AdaptiveRangeDecoder<T>::~AdaptiveRangeDecoder()
{
}

template<class T>
bool AdaptiveRangeDecoder<T>::initialize(T& stream)
{
    frequencies_.initialize();
    stream_ = &stream;
    umcompressedSize_ = 0;
    outSize_ = 0;
    state_ = State_Init;
    range_ = MAXRANGE;
    low_ = 0;
    return true;
}

template<class T>
Result AdaptiveRangeDecoder<T>::decode(s32 size, const u8* bytes)
{
    u32 t, count;
    switch(state_) {
    case State_Init:
        if(size < 8) {
            return {Status_Pending, 8};
        }

        umcompressedSize_ = bytes[3];
        umcompressedSize_ = (umcompressedSize_ << 8) | bytes[2];
        umcompressedSize_ = (umcompressedSize_ << 8) | bytes[1];
        umcompressedSize_ = (umcompressedSize_ << 8) | bytes[0];
        bytes += 4;
        size -= 4;

        low_ = bytes[0];
        low_ = (low_ << 8) | bytes[1];
        low_ = (low_ << 8) | bytes[2];
        low_ = (low_ << 8) | bytes[3];
        bytes += 4;
        size -= 4;
        state_ = State_Decode;
        goto STATE_DECODE;

    case State_Decode:
    STATE_DECODE:
        for(;;) {
            if(!normalize(size, bytes)) {
                return {Status_Pending, umcompressedSize_ - outSize_};
            }
            t = range_ / frequencies_.total();
            frequencies_.find(count, lastCode_, low_ / t);
            low_ -= t * count;
            range_ = t * frequencies_[lastCode_];

            if(!stream_->writeByte(lastCode_)) {
                return {Status_Pending, umcompressedSize_ - outSize_};
            }
            if(umcompressedSize_ <= ++outSize_) {
                return {Status_Success, 0};
            }

            frequencies_.update(lastCode_);
        }
        break;
    default:
        CPPRCODER_ASSERT(false);
        break;
    }
    return {Status_Error, 0};
}

template<class T>
bool AdaptiveRangeDecoder<T>::normalize(s32& size, const u8*& bytes)
{
    while(range_ < MINRANGE) {
        if(size <= 0) {
            return false;
        }
        range_ <<= 8;
        u8 byte = bytes[0];
        --size;
        ++bytes;
        low_ = (low_ << 8) + byte;
    }
    return true;
}

} // namespace cpprcoder

#ifdef CPPRCODER_IMPLEMENTATION

#    include <cassert>

#    if defined(CPPRCODER_CREATE_MEMORYSTREAM) || defined(_MSCVER)
#        include <stdlib.h>
#    endif

#    ifdef CPPRCODER_CREATE_MEMORYSTREAM
#        include <string.h>
#    endif

namespace cpprcoder
{
#    ifdef CPPRCODER_CREATE_MEMORYSTREAM
//----------------------------------------------
//---
//--- MemoryStream
//---
//----------------------------------------------
MemoryStream::MemoryStream()
    : capacity_(0)
    , size_(0)
    , buffer_(CPPRCODER_NULL)
{
}

MemoryStream::MemoryStream(s32 capacity)
    : capacity_(capacity)
    , size_(0)
{
    capacity_ = (capacity_ <= 0) ? 16 : (capacity_ + 15U) & (~15U);
    CPPRCODER_ASSERT(0 <= capacity_);
    buffer_ = reinterpret_cast<u8*>(CPPRCODER_MALLOC(capacity_));
}

MemoryStream::~MemoryStream()
{
    CPPRCODER_FREE(buffer_);
}

void MemoryStream::reserve(s32 capacity)
{
    capacity = (capacity + 15U) & (~15U);
    if(capacity < capacity_) {
        return;
    }
    CPPRCODER_FREE(buffer_);
    capacity_ = capacity;
    buffer_ = reinterpret_cast<u8*>(CPPRCODER_MALLOC(capacity_));
}

void MemoryStream::resize(s32 size)
{
    CPPRCODER_ASSERT(0 <= size);
    if(capacity_ < size) {
        reserve(size);
    }
    size_ = size;
}

s32 MemoryStream::read(s32 size, u8* bytes)
{
    CPPRCODER_ASSERT(0 <= size);
    CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);

    s32 end = size_ + size;
    if(capacity_ < end) {
        return -1;
    }
    for(s32 i = 0; size_ < end; ++i, ++size_) {
        bytes[i] = buffer_[size_];
    }
    return size;
}

bool MemoryStream::readByte(u8& byte)
{
    CPPRCODER_ASSERT(CPPRCODER_NULL != buffer_);
    s32 end = size_ + 1;
    if(capacity_ < end) {
        return false;
    }
    byte = buffer_[size_++];
    return true;
}

s32 MemoryStream::write(s32 size, const u8* bytes)
{
    CPPRCODER_ASSERT(0 <= size);
    CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);
    s32 end = size_ + size;
    if(capacity_ < end) {
        if(!expand(end)) {
            return -1;
        }
    }
    for(s32 i = 0; size_ < end; ++i, ++size_) {
        buffer_[size_] = bytes[i];
    }
    return size;
}

bool MemoryStream::writeByte(u8 byte)
{
    if(capacity_ <= size_) {
        return false;
    }
    buffer_[size_++] = byte;
    return true;
}

bool MemoryStream::expand(s32 size)
{
    s32 prev = capacity_;
    s32 capacity = 0;
    do {
        if(prev <= 0) {
            prev = capacity = 1024;
        } else if(prev < EXPAND_LIMIT_SIZE) {
            prev = capacity = prev << 1;
        } else {
            prev = capacity = prev + EXPAND_LIMIT_SIZE;
        }
    } while(capacity < size);

    capacity = (capacity + 15U) & (~15U);
    u8* buffer = reinterpret_cast<u8*>(CPPRCODER_MALLOC(capacity));
    memcpy(buffer, buffer_, capacity_);
    CPPRCODER_FREE(buffer_);
    capacity_ = capacity;
    buffer_ = buffer;
    return CPPRCODER_NULL != buffer_;
}
#    endif //CPPRCODER_CREATE_MEMORYSTREAM

//----------------------------------------------
//---
//--- AdaptiveFrequencyTable
//---
//----------------------------------------------
AdaptiveFrequencyTable::AdaptiveFrequencyTable()
    : total_(0)
{
}

AdaptiveFrequencyTable::~AdaptiveFrequencyTable()
{
}

void AdaptiveFrequencyTable::initialize()
{
    total_ = SIZE;
#    if defined(CPPRCODER_USE_SIMD)
    prefix_ = reinterpret_cast<u32*>((reinterpret_cast<uintptr_t>(prefix_buffer_) + ALIGN_MASK) & ~ALIGN_MASK);
    frequencies_ = reinterpret_cast<u32*>((reinterpret_cast<uintptr_t>(frequencies_buffer_) + ALIGN_MASK) & ~ALIGN_MASK);

    __m128i one = _mm_set1_epi32(1);
    __m128i* p = reinterpret_cast<__m128i*>(frequencies_);
    __m128i* end = p + (SIZE / 4);
    while(p < end) {
        _mm_store_si128(p, one);
        ++p;
    }

    u32 prefix[] =
        {
            CHUNK_SIZE * 1,
            CHUNK_SIZE * 2,
            CHUNK_SIZE * 3,
            CHUNK_SIZE * 4,
        };
    __m128i base = _mm_set1_epi32(prefix[3]);
    __m128i total = _mm_loadu_si128(reinterpret_cast<const __m128i*>(prefix));
    for(s32 i = 0; i < CHUNKS; i += 4) {
        _mm_store_si128(reinterpret_cast<__m128i*>(&prefix_[i]), total);
        total = _mm_add_epi32(total, base);
    }
#    else
    for(u32 i = 0; i < SIZE; ++i) {
        frequencies_[i] = 1;
    }

    prefix_[0] = CHUNK_SIZE;
    for(s32 i = 1; i < CHUNKS; ++i) {
        prefix_[i] = prefix_[i - 1] + CHUNK_SIZE;
    }
#    endif
}

void AdaptiveFrequencyTable::update(u8 b)
{
#    if defined(CPPRCODER_USE_SIMD)
    ++frequencies_[b];
    if(MINRANGE <= (++total_)) {
        __m128i one = _mm_set1_epi32(1);
        __m128i total = _mm_setzero_si128();
        __m128i* p = reinterpret_cast<__m128i*>(frequencies_);
        __m128i* end = p + (SIZE / 4);
        while(p < end) {
            __m128i v = _mm_load_si128(p);
            v = _mm_srli_epi32(v, 1);
            v = _mm_or_si128(v, one);
            total = _mm_add_epi32(total, v);
            _mm_store_si128(p, v);
            ++p;
        }
        u32 count[4];
        _mm_storeu_si128(reinterpret_cast<__m128i*>(count), total);
        total_ = count[0] + count[1] + count[2] + count[3];
        countChunks();
    } else {
        b >>= 4;
        for(u8 i = b; i < CHUNKS; ++i) {
            ++prefix_[i];
        }
    }
#    else
    ++frequencies_[b];
    if(MINRANGE <= (++total_)) {
        total_ = 0;
        for(u32 i = 0; i < SIZE; ++i) {
            frequencies_[i] = (frequencies_[i] >> 1) | 0x01U;
            total_ += frequencies_[i];
        }
        countChunks();
    } else {
        b >>= 4;
        for(u8 i = b; i < CHUNKS; ++i) {
            ++prefix_[i];
        }
    }
#    endif
}

u32 AdaptiveFrequencyTable::cumulative(u8 byte) const
{
    u32 chunk = byte >> 4;
    u32 count = (0 < chunk) ? prefix_[chunk - 1] : 0;
    for(u32 i = (chunk << 4); i < byte; ++i) {
        count += frequencies_[i];
    }
    return count;
}

void AdaptiveFrequencyTable::find(u32& count, u8& code, u32 target) const
{
#    if 0 //#if defined(CPPRCODER_USE_SIMD)
        __m128i target128 = _mm_set1_epi32(target);
        __m128i one = _mm_set1_epi32(1);

        code = 0;

        __m128i* pprefix = reinterpret_cast<__m128i*>(prefix_);
        __m128i tmp = _mm_load_si128(&pprefix[0]);

        __m128i chunk_total128 = _mm_and_si128(_mm_cmplt_epi32(target128, _mm_load_si128(&pprefix[0])), one);
        chunk_total128 = _mm_add_epi32(chunk_total128, _mm_and_si128(_mm_cmplt_epi32(target128, _mm_load_si128(&pprefix[1])), one));
        chunk_total128 = _mm_add_epi32(chunk_total128, _mm_and_si128(_mm_cmplt_epi32(target128, _mm_load_si128(&pprefix[2])), one));
        chunk_total128 = _mm_add_epi32(chunk_total128, _mm_and_si128(_mm_cmplt_epi32(target128, _mm_load_si128(&pprefix[3])), one));

        u32 tmp_chunk_total[4];
        _mm_storeu_si128(reinterpret_cast<__m128i*>(tmp_chunk_total), chunk_total128);
        u32 chunk = tmp_chunk_total[0] + tmp_chunk_total[1] + tmp_chunk_total[2] + tmp_chunk_total[3];
        chunk = (0<chunk)? CHUNKS - chunk : chunk;
        count = (0<chunk)? prefix_[chunk-1] : 0;

        code = chunk<<4;

        for(u32 i = code; i<SIZE; ++i){
            if(target<(count+frequencies_[i])){
                code = static_cast<u8>(i);
                return;
            }
            count += frequencies_[i];
        }
#    else
    count = 0;
    code = 0;

    u8 chunk = 0;
    for(u8 i = 0; i < CHUNKS; ++i) {
        if(target < prefix_[i]) {
            chunk = i;
            count = (0 < chunk) ? prefix_[chunk - 1] : 0;
            break;
        }
    }
    code = chunk << 4;

    for(u32 i = code; i < SIZE; ++i) {
        u32 next = count + frequencies_[i];
        if(target < next) {
            code = static_cast<u8>(i);
            return;
        }
        count = next;
    }
#    endif
}

void AdaptiveFrequencyTable::countChunks()
{
    u32 sum = frequencies_[0];
    for(s32 j = 1; j < CHUNK_SIZE; ++j) {
        sum += frequencies_[j];
    }
    prefix_[0] = sum;

    for(s32 i = 1; i < CHUNKS; ++i) {
        s32 freq = i << 4;
        sum = frequencies_[freq];
        for(s32 j = 1; j < CHUNK_SIZE; ++j) {
            sum += frequencies_[freq + j];
        }
        prefix_[i] = prefix_[i - 1] + sum;
    }
}

} // namespace cpprcoder
#endif

#endif //INC_CPPRCODER_H_
