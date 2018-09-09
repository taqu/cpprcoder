#ifndef INC_CPPRCODER_H_
#define INC_CPPRCODER_H_
/*
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
*/
/**
@author t-sakai
@date 2018/09/04 create

USAGE:
Put '#define CPPRCODER_IMPLEMENTATION' before including this file to create the implementation.
*/
#include <cassert>
#include <cstdint>

#define CPPRCODER_CREATE_FILESTREAM
#define CPPRCODER_CREATE_MEMORYSTREAM


#ifdef CPPRCODER_CREATE_FILESTREAM
#include <stdio.h>
#endif

#if defined(CPPRCODER_CREATE_MEMORYSTREAM) || defined(_MSCVER)
#include <stdlib.h>
#endif

//#define CPPRCODER_IMPLEMENTATION

namespace cpprcoder
{
#ifdef __cplusplus
#   if 201103L<=__cplusplus || 1900<=_MSC_VER
#       define CPPRCODER_CPP11 1
#   endif
#endif

#ifdef __cplusplus
#   ifdef CPPRCODER_CPP11
#       define CPPRCODER_NULL nullptr
#   else
#       define CPPRCODER_NULL 0
#   endif
#else
#   define CPPRCODER_NULL (void*)0
#endif

#ifndef CPPRCODER_TYPES
#define CPPRCODER_TYPES
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

#ifndef CPPRCODER_ASSERT
#define CPPRCODER_ASSERT(exp) assert(exp)
#endif

#ifndef CPPRCODER_OFF_T
#define CPPRCODER_OFF_T
    typedef s64 off_t;
#endif


#ifndef CPPRCODER_MINRANGE
#define CPPRCODER_MINRANGE (0x1000000U)
#endif

#ifndef CPPRCODER_MAXRANGE_ENCODE
#define CPPRCODER_MAXRANGE_ENCODE (0xFFFFFF00U)
#endif

#ifndef CPPRCODER_MAXRANGE_DECODE
#define CPPRCODER_MAXRANGE_DECODE (0x00FFFFFFU)
#endif

#if defined(__AVX__)
#define CPPRCODER_USE_AVX (1)
#elif defined(__SSE2__)
#define CPPRCODER_USE_SSE2 (1)
#endif


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
    class IStream
    {
    public:
        virtual s32 read(s32 size, u8* bytes) =0;
        virtual bool readByte(u8& byte) =0;
        virtual s32 write(s32 size, const u8* bytes) =0;
        virtual bool writeByte(u8 byte) =0;

    protected:
        IStream()
        {}
        virtual ~IStream()
        {}

    private:
        IStream(const IStream&) = delete;
        IStream(IStream&&) = delete;
        IStream& operator=(const IStream&) = delete;
        IStream& operator=(IStream&&) = delete;
    };

#ifdef CPPRCODER_CREATE_FILESTREAM
    enum IOS
    {
        IOS_Read,
        IOS_Write,
    };

    //----------------------------------------------
    //---
    //--- FileStream
    //---
    //----------------------------------------------
    class FileStream : public IStream
    {
    public:
        FileStream();
        virtual ~FileStream();

        bool open(const Char* path, IOS mode);
        void close();

        bool is_open() const;

        virtual s32 read(s32 size, u8* bytes);
        virtual bool readByte(u8& byte);
        virtual s32 write(s32 size, const u8* bytes);
        virtual bool writeByte(u8 byte);

    private:
        FileStream(const FileStream&) = delete;
        FileStream(FileStream&&) = delete;
        FileStream& operator=(const FileStream&) = delete;
        FileStream& operator=(FileStream&&) = delete;

        FILE* file_;
    };
#endif

#ifdef CPPRCODER_CREATE_MEMORYSTREAM

#ifndef CPPRCODER_MALLOC
#define CPPRCODER_MALLOC(size) malloc(size)
#endif

#ifndef CPPRCODER_FREE
#define CPPRCODER_FREE(ptr) free((ptr));(ptr) = CPPRCODER_NULL
#endif
    //----------------------------------------------
    //---
    //--- MemoryStream
    //---
    //----------------------------------------------
    class MemoryStream : public IStream
    {
    public:
        MemoryStream();
        explicit MemoryStream(s32 capacity);
        virtual ~MemoryStream();

        s32 capacity() const;
        s32 size() const;

        const u8& operator[](s32 index) const;
        u8& operator[](s32 index);

        void reserve(s32 capacity);
        void resize(s32 size);

        virtual s32 read(s32 size, u8* bytes);
        virtual bool readByte(u8& byte);
        virtual s32 write(s32 size, const u8* bytes);
        virtual bool writeByte(u8 byte);

    private:
        MemoryStream(const MemoryStream&) = delete;
        MemoryStream(MemoryStream&&) = delete;
        MemoryStream& operator=(const MemoryStream&) = delete;
        MemoryStream& operator=(MemoryStream&&) = delete;

        static const s32 DOUBLE_LIMIT_SIZE = 4096*4;

        bool expand(s32 size);

        s32 capacity_;
        s32 size_;
        u8* buffer_;
    };
#endif

    //----------------------------------------------
    //---
    //--- FrequencyTable
    //---
    //----------------------------------------------
    class FrequencyTable
    {
    public:
        static const u8 CHUNK_SIZE = 16;
        static const u32 SIZE = 256;
        static const u8 CHUNKS = (SIZE/CHUNK_SIZE);
        static const uintptr_t ALIGN = 16;
        static const uintptr_t ALIGN_MASK = ALIGN-1;

        FrequencyTable();
        ~FrequencyTable();

        void update(u8 b);
        u32 size() const;
        u32 total() const;
        u32 operator[](u32 index) const;
        u32 cumulative(u32 size) const;
        void find(u32& count, u8& code, u32 target) const;
    private:
        FrequencyTable(const FrequencyTable&) = delete;
        FrequencyTable(FrequencyTable&&) = delete;
        FrequencyTable& operator=(const FrequencyTable&) = delete;
        FrequencyTable& operator=(FrequencyTable&&) = delete;

        void countChunks();

        u32 total_;
#if defined(CPPRCODER_USE_SSE2) || defined(CPPRCODER_USE_AVX)
        u32* prefix_;
        u32* frequencies_;
        u32 prefix_buffer_[CHUNKS + 4];
        u32 frequencies_buffer_[SIZE + 4];
#else
        u32 prefix_[CHUNKS];
        u32 frequencies_[SIZE];
#endif
    };

#if 0
    //----------------------------------------------
    //---
    //--- BinaryIndexedTree
    //---
    //----------------------------------------------
    class BinaryIndexedTree
    {
    public:
        static const u32 SIZE = 256;

        BinaryIndexedTree();
        ~BinaryIndexedTree();

        /**
        @brief set all frequencies to value
        @param value ... frequency for all
        */
        void initialize(u32 value);

        void update(u32 index);
        void add(u32 index, s32 value);
        u32 size() const;
        u32 total() const;
        u32 operator[](u32 index) const;
        u32 cumulative(u32 index) const;
        void find(u32& count, u8& code, u32 target) const;
    private:
        BinaryIndexedTree(const BinaryIndexedTree&) = delete;
        BinaryIndexedTree(BinaryIndexedTree&&) = delete;
        BinaryIndexedTree& operator=(const BinaryIndexedTree&) = delete;
        BinaryIndexedTree& operator=(BinaryIndexedTree&&) = delete;

        static inline u32 LSB(u32 x)
        {
            return x & (x-1);
        }

        u32 total_;
        u32 mid_;
        u32 frequencies_[SIZE];
    };
#endif

    //----------------------------------------------
    //---
    //--- AdaptiveRangeEncoder
    //---
    //----------------------------------------------
    class AdaptiveRangeEncoder
    {
    public:
        AdaptiveRangeEncoder();
        ~AdaptiveRangeEncoder();

        bool initialize(IStream& stream, u32 umcompressedSize);
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

        IStream* stream_;
        FrequencyTable frequencies_;
        u32 umcompressedSize_;
        u32 inSize_;
        u32 range_;
        u32 low_;
        u8 buffer_;
        u32 carry_;
    };

    //----------------------------------------------
    //---
    //--- AdaptiveRangeDecoder
    //---
    //----------------------------------------------
    class AdaptiveRangeDecoder
    {
    public:
        AdaptiveRangeDecoder();
        ~AdaptiveRangeDecoder();

        bool initialize(IStream& stream);
        Result decode(s32 size, const u8* bytes);
    private:
        AdaptiveRangeDecoder(const AdaptiveRangeDecoder&) = delete;
        AdaptiveRangeDecoder(AdaptiveRangeDecoder&&) = delete;
        AdaptiveRangeDecoder& operator=(const AdaptiveRangeDecoder&) = delete;
        AdaptiveRangeDecoder& operator=(AdaptiveRangeDecoder&&) = delete;

        static const s32 State_Init = 0;
        static const s32 State_Decode = 1;

        bool normalize(s32& size, const u8*& bytes);

        IStream* stream_;
        FrequencyTable frequencies_;
        u32 umcompressedSize_;
        u32 outSize_;
        s32 state_;
        u32 range_;
        u32 low_;
        u8 lastCode_;
    };
}

#ifdef CPPRCODER_IMPLEMENTATION

#ifdef CPPRCODER_USE_SSE2
#include <immintrin.h>
#include <emmintrin.h>
#endif

#ifdef CPPRCODER_USE_AVX
#include <immintrin.h>
#endif

#ifdef CPPRCODER_CREATE_MEMORYSTREAM
#include <string.h>
#endif

namespace cpprcoder
{
#ifdef CPPRCODER_CREATE_FILESTREAM
    //----------------------------------------------
    //---
    //--- FileStream
    //---
    //----------------------------------------------
    FileStream::FileStream()
        :file_(CPPRCODER_NULL)
    {
    }

    FileStream::~FileStream()
    {
        close();
    }

    bool FileStream::open(const Char* path, IOS mode)
    {
        CPPRCODER_ASSERT(CPPRCODER_NULL != path);
        CPPRCODER_ASSERT(IOS_Read<=mode && mode<=IOS_Write);

        if(CPPRCODER_NULL != file_){
            fclose(file_);
            file_ = CPPRCODER_NULL;
        }
        
        static const Char* smode[] = {"rb","wb"};
#if _MSC_VER
        fopen_s(&file_, path, smode[mode]);
#else
        file_ = fopen(path, smode[mode]);
#endif
        return CPPRCODER_NULL != file_;
    }

    void FileStream::close()
    {
        if(CPPRCODER_NULL != file_){
            fclose(file_);
            file_ = CPPRCODER_NULL;
        }
    }

    bool FileStream::is_open() const
    {
        return CPPRCODER_NULL != file_;
    }

    s32 FileStream::read(s32 size, u8* bytes)
    {
        CPPRCODER_ASSERT(CPPRCODER_NULL != file_);
        CPPRCODER_ASSERT(0<=size);
        CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);
        return (0<fread(bytes, size, 1, file_))? size : -1;
    }

    bool FileStream::readByte(u8& byte)
    {
        CPPRCODER_ASSERT(CPPRCODER_NULL != file_);
        s32 r = fgetc(file_);
        if(r==EOF){
            return false;
        }
        byte = static_cast<u8>(r);
        return true;
    }

    s32 FileStream::write(s32 size, const u8* bytes)
    {
        CPPRCODER_ASSERT(CPPRCODER_NULL != file_);
        CPPRCODER_ASSERT(0<=size);
        CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);
        return (0<fwrite(bytes, size, 1, file_))? size : -1;
    }

    bool FileStream::writeByte(u8 byte)
    {
        CPPRCODER_ASSERT(CPPRCODER_NULL != file_);
        return EOF != fputc(byte, file_);
    }
#endif

#ifdef CPPRCODER_CREATE_MEMORYSTREAM
    //----------------------------------------------
    //---
    //--- MemoryStream
    //---
    //----------------------------------------------
    MemoryStream::MemoryStream()
        :capacity_(0)
        ,size_(0)
        ,buffer_(CPPRCODER_NULL)
    {
    }

    MemoryStream::MemoryStream(s32 capacity)
        :capacity_(capacity)
        ,size_(0)
    {
        capacity_ = (capacity_<=0)? 16 : (capacity_+15U) & (~15U);
        CPPRCODER_ASSERT(0<=capacity_);
        buffer_ = reinterpret_cast<u8*>(CPPRCODER_MALLOC(capacity_));
    }

    MemoryStream::~MemoryStream()
    {
        CPPRCODER_FREE(buffer_);
    }

    s32 MemoryStream::capacity() const
    {
        return capacity_;
    }

    s32 MemoryStream::size() const
    {
        return size_;
    }

    const u8& MemoryStream::operator[](s32 index) const
    {
        CPPRCODER_ASSERT(0<=index && index<size_);
        return buffer_[index];
    }

    u8& MemoryStream::operator[](s32 index)
    {
        CPPRCODER_ASSERT(0<=index && index<size_);
        return buffer_[index];
    }

    void MemoryStream::reserve(s32 capacity)
    {
        capacity = (capacity + 15U) & (~15U);
        if(capacity<capacity_){
            return;
        }
        CPPRCODER_FREE(buffer_);
        capacity_ = capacity;
        buffer_ = reinterpret_cast<u8*>(CPPRCODER_MALLOC(capacity_));
    }

    void MemoryStream::resize(s32 size)
    {
        CPPRCODER_ASSERT(0<=size);
        if(capacity_<size){
            reserve(size);
        }
        size_ = size;
    }

    s32 MemoryStream::read(s32 size, u8* bytes)
    {
        CPPRCODER_ASSERT(0<=size);
        CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);

        s32 end = size_ + size;
        if(capacity_<end){
            return -1;
        }
        for(s32 i=0; size_<end; ++i, ++size_){
            bytes[i] = buffer_[size_];
        }
        return size;
    }

    bool MemoryStream::readByte(u8& byte)
    {
        CPPRCODER_ASSERT(CPPRCODER_NULL != buffer_);
        s32 end = size_ + 1;
        if(capacity_<end){
            return false;
        }
        byte = buffer_[size_++];
        return true;
    }

    s32 MemoryStream::write(s32 size, const u8* bytes)
    {
        CPPRCODER_ASSERT(0<=size);
        CPPRCODER_ASSERT(CPPRCODER_NULL != bytes);
        s32 end = size_ + size;
        if(capacity_<end){
            if(!expand(end)){
                return -1;
            }
        }
        for(s32 i=0; size_<end; ++i, ++size_){
            buffer_[size_] = bytes[i];
        }
        return size;
    }

    bool MemoryStream::writeByte(u8 byte)
    {
        s32 end = size_ + 1;
        if(capacity_<end){
            if(!expand(end)){
                return false;
            }
        }
        buffer_[size_++] = byte;
        return true;
    }

    bool MemoryStream::expand(s32 size)
    {
        s32 prev = capacity_;
        s32 capacity = 0;
        do{
            if(prev<=0){
                prev = capacity = 1024;
            } else if(prev<DOUBLE_LIMIT_SIZE){
                prev = capacity = prev<<1;
            } else {
                prev = capacity = prev + DOUBLE_LIMIT_SIZE;
            }
        } while(capacity<size);

        capacity = (capacity+15U) & (~15U);
        u8* buffer = reinterpret_cast<u8*>(CPPRCODER_MALLOC(capacity));
        memcpy(buffer, buffer_, capacity_);
        CPPRCODER_FREE(buffer_);
        capacity_ = capacity;
        buffer_ = buffer;
        return CPPRCODER_NULL != buffer_;
    }
#endif

    //----------------------------------------------
    //---
    //--- FrequencyTable
    //---
    //----------------------------------------------
    FrequencyTable::FrequencyTable()
        :total_(SIZE)
    {
#if defined(CPPRCODER_USE_SSE2) || defined(CPPRCODER_USE_AVX)
        prefix_ = reinterpret_cast<u32*>((reinterpret_cast<uintptr_t>(prefix_buffer_) + ALIGN_MASK) & ~ALIGN_MASK);
        frequencies_ = reinterpret_cast<u32*>((reinterpret_cast<uintptr_t>(frequencies_buffer_) + ALIGN_MASK) & ~ALIGN_MASK);

        __m128i one = _mm_set1_epi32(1);
        __m128i* p = reinterpret_cast<__m128i*>(frequencies_);
        __m128i* end = p + (SIZE/4);
        while(p<end){
            _mm_store_si128(p, one);
            ++p;
        }

        u32 prefix[] =
        {
            CHUNK_SIZE*1, CHUNK_SIZE*2, CHUNK_SIZE*3, CHUNK_SIZE*4,
        };
        __m128i base = _mm_set1_epi32(prefix[3]);
        __m128i total = _mm_loadu_si128(reinterpret_cast<const __m128i*>(prefix));
        for(s32 i=0; i<CHUNKS; i+=4){
            _mm_store_si128(reinterpret_cast<__m128i*>(&prefix_[i]), total);
            total = _mm_add_epi32(total, base);
        }
#else
        for(u32 i=0; i<SIZE; ++i){
            frequencies_[i] = 1;
        }

        for(s32 i=1; i<CHUNKS; ++i){
            prefix_[i] = prefix_[i-1] + CHUNK_SIZE;
        }
#endif
    }

    FrequencyTable::~FrequencyTable()
    {}

    void FrequencyTable::update(u8 b)
    {
#if defined(CPPRCODER_USE_SSE2) || defined(CPPRCODER_USE_AVX)
        ++frequencies_[b];
        if(CPPRCODER_MINRANGE<=(++total_)){
            __m128i one = _mm_set1_epi32(1);
            __m128i total = _mm_setzero_si128();
            __m128i* p = reinterpret_cast<__m128i*>(frequencies_);
            __m128i* end = p + (SIZE/4);
            while(p<end){
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
        }else{
            b >>= 4;
            for(u8 i=b; i<CHUNKS; ++i){
                ++prefix_[i];
            }
        }
#else
        ++frequencies_[b];
        if(CPPRCODER_MINRANGE<=(++total_)){
            total_ = 0;
            for(u32 i=0; i<SIZE; ++i){
                frequencies_[i] = (frequencies_[i]>>1) | 0x01U;
                total_ += frequencies_[i];
            }
            countChunks();
        }else{
            b >>= 4;
            for(u8 i=b; i<CHUNKS; ++i){
                ++prefix_[i];
            }
        }
#endif
    }

    u32 FrequencyTable::size() const
    {
        return SIZE;
    }

    u32 FrequencyTable::total() const
    {
        return total_;
    }

    u32 FrequencyTable::operator[](u32 index) const
    {
        CPPRCODER_ASSERT(0<=index && index<SIZE);
        return frequencies_[index];
    }

    u32 FrequencyTable::cumulative(u32 size) const
    {
        u32 chunk = size>>4;
        u32 count = (0<chunk)? prefix_[chunk-1] : 0;
        for(u32 i=(chunk<<4); i<size; ++i){
            count += frequencies_[i];
        }
        return count;
    }

    void FrequencyTable::find(u32& count, u8& code, u32 target) const
    {
#if 0 //#if defined(CPPRCODER_USE_SSE2) || defined(CPPRCODER_USE_AVX)
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

        for(u32 i=code; i<SIZE; ++i){
            if(target<(count+frequencies_[i])){
                code = static_cast<u8>(i);
                return;
            }
            count += frequencies_[i];
        }
#else
        count = 0;
        code = 0;

        u8 chunk = 0;
        for(u8 i=0; i<CHUNKS; ++i){
            if(target<prefix_[i]){
                chunk = i;
                count = (0<chunk)? prefix_[chunk-1] : 0;
                break;
            }
        }
        code = chunk<<4;

        for(u32 i=code; i<SIZE; ++i){
            u32 next = count + frequencies_[i];
            if(target<next){
                code = static_cast<u8>(i);
                return;
            }
            count = next;
        }
#endif
    }

    void FrequencyTable::countChunks()
    {
        u32 sum = frequencies_[0];
        for(s32 j=1; j<CHUNK_SIZE; ++j){
            sum += frequencies_[j];
        }
        prefix_[0] = sum;

        for(s32 i=1; i<CHUNKS; ++i){
            s32 freq = i<<4;
            sum = frequencies_[freq];
            for(s32 j=1; j<CHUNK_SIZE; ++j){
                sum += frequencies_[freq+j];
            }
            prefix_[i] = prefix_[i-1] + sum;
        }
    }

#if 0
    //----------------------------------------------
    //---
    //--- BinaryIndexedTree
    //---
    //----------------------------------------------
    BinaryIndexedTree::BinaryIndexedTree()
        :total_(0)
    {
        initialize(1);
    }

    BinaryIndexedTree::~BinaryIndexedTree()
    {}

    void BinaryIndexedTree::initialize(u32 value)
    {
#if defined(CPPRCODER_USE_SSE2) || defined(CPPRCODER_USE_AVX)
        __m128i one = _mm_set1_epi32(0);
        __m128i* p = reinterpret_cast<__m128i*>(frequencies_);
        __m128i* end = p + (SIZE/4);
        while(p<end){
            _mm_storeu_si128(p, one);
            ++p;
        }
#else
        for(u32 i = 0; i<SIZE; ++i){
            frequencies_[i] = 0;
        }
#endif
        total_ = value*SIZE;
        frequencies_[0] = value;
        for(u32 i=1; i<SIZE; ++i){
            u32 index = i;
            do{
                frequencies_[index] += value;
                index += index & -index;
            }while(index<SIZE);
        }
        for(mid_=1; (mid_<<1)<(SIZE-1); mid_<<=1);
    }

    void BinaryIndexedTree::update(u32 index)
    {
        add(index, 1);
        if(CPPRCODER_MINRANGE<=total_){
            for(u32 i=0; i<SIZE; ++i){
                s32 frequency = static_cast<s32>((*this)[i] >> 1);
                if(0<frequency){
                    add(i, -frequency);
                }
            }
        }
    }

    void BinaryIndexedTree::add(u32 index, s32 value)
    {
        if(0<index){
            do{
                frequencies_[index] += value;
                index += index & -index;
            }while(index<SIZE);
        }else{
            frequencies_[index] += value;
        }
        total_ += value;
    }

    u32 BinaryIndexedTree::size() const
    {
        return SIZE;
    }

    u32 BinaryIndexedTree::total() const
    {
        return total_;
    }

    u32 BinaryIndexedTree::operator[](u32 index) const
    {
        CPPRCODER_ASSERT(0<=index && index<SIZE);
        if((index&0x01U) || (index<=0)){
            return frequencies_[index];
        }
        u32 frequency = frequencies_[index];
        u32 end=LSB(index);
        index -=1;
        for(; end!=index; index = LSB(index)){
            frequency -= frequencies_[index];
        }
        return frequency;
    }

    u32 BinaryIndexedTree::cumulative(u32 index) const
    {
        if(index<=0){
            return 0;
        }
        u32 count = frequencies_[0];
        index -= 1;
        while(0<index){
            count += frequencies_[index];
            index = LSB(index);
        }
        return count;
    }

    void BinaryIndexedTree::find(u32& count, u8& code, u32 target) const
    {
        code = 0;
        if(target<frequencies_[0]){
            count = 0;
            return;
        }
        u32 mid = mid_;
        count = frequencies_[0];
        while(0<mid){
            u32 t;
            if(((code + mid)<SIZE) && ((t = frequencies_[code + mid])<=target)){
                code += mid;
                target -= t;
            }
            mid >>= 1;
        }
        code += 1;
    }
#endif

    //----------------------------------------------
    //---
    //--- AdaptiveRangeEncoder
    //---
    //----------------------------------------------
    AdaptiveRangeEncoder::AdaptiveRangeEncoder()
        :stream_(CPPRCODER_NULL)
        ,umcompressedSize_(0)
        ,inSize_(0)
        ,range_(CPPRCODER_MAXRANGE_ENCODE)
        ,low_(0)
        ,buffer_(0)
        ,carry_(0)
    {
    }

    AdaptiveRangeEncoder::~AdaptiveRangeEncoder()
    {
    }

    bool AdaptiveRangeEncoder::initialize(IStream& stream, u32 umcompressedSize)
    {
        stream_ = &stream;
        range_ = CPPRCODER_MAXRANGE_ENCODE;
        umcompressedSize_ = umcompressedSize;
        inSize_ = 0;
        low_ = 0;
        buffer_ = 0;
        carry_ = 0;
        u8 bytes[4];
        bytes[0] = (umcompressedSize_>> 0) & 0xFFU;
        bytes[1] = (umcompressedSize_>> 8) & 0xFFU;
        bytes[2] = (umcompressedSize_>>16) & 0xFFU;
        bytes[3] = (umcompressedSize_>>24) & 0xFFU;
        return stream_->write(4, bytes);
    }

    Result AdaptiveRangeEncoder::encode(s32 size, const u8* bytes)
    {
        CPPRCODER_ASSERT((inSize_+size)<=umcompressedSize_);

        for(s32 i=0; i<size; ++i){
            u32 t = range_ / frequencies_.total();
            u8 c = bytes[i];
            u32 prevLow = low_;
            low_ += frequencies_.cumulative(c) * t;
            range_ = frequencies_[c] * t;
            if(!normalize(prevLow)){
                inSize_ += i;
                return {Status_Pending, umcompressedSize_-inSize_};
            }
            frequencies_.update(c);
        }
        inSize_ += size;
        if(umcompressedSize_<=inSize_){
            finish();
            return {Status_Success, 0};
        }
        return {Status_Pending, umcompressedSize_-inSize_};
    }

    Result AdaptiveRangeEncoder::encode(u8 byte)
    {
        CPPRCODER_ASSERT((inSize_+1)<=umcompressedSize_);

        u32 t = range_ / frequencies_.total();
        u32 prevLow = low_;
        low_ += frequencies_.cumulative(byte) * t;
        range_ = frequencies_[byte] * t;
        if(!normalize(prevLow)){
            inSize_ += 1;
            return {Status_Pending, umcompressedSize_-inSize_};
        }
        frequencies_.update(byte);
        ++inSize_;
        if(umcompressedSize_<=inSize_){
            finish();
            return {Status_Success, 0};
        }
        return {Status_Pending, umcompressedSize_-inSize_};
    }

    bool AdaptiveRangeEncoder::finish()
    {
        if(!stream_->writeByte(buffer_)){
            return false;
        }
        for(u32 i=0; i<carry_; ++i){
            if(!stream_->writeByte(0xFFU)){
                return false;
            }
        }

        u8 bytes[4];
        bytes[0] = (low_>>24) & 0xFFU;
        bytes[1] = (low_>>16) & 0xFFU;
        bytes[2] = (low_>> 8) & 0xFFU;
        bytes[3] = (low_>> 0) & 0xFFU;
        return 0<stream_->write(4, bytes);
    }

    bool AdaptiveRangeEncoder::normalize(u32 prevLow)
    {
        if(low_<prevLow){
            buffer_ += 1;
            if(0<carry_){
                if(!stream_->writeByte(buffer_)){
                    return false;
                }
                for(u32 i=1; i<carry_; ++i){
                    if(!stream_->writeByte(0)){
                        return false;
                    }
                }
                buffer_ = 0;
                carry_ = 0;
            }
        }

        while(range_<CPPRCODER_MINRANGE){
            if(low_<(0xFFU<<SHIFT)){
                if(!stream_->writeByte(buffer_)){
                    return false;
                }
                for(u32 i=0; i<carry_; ++i){
                    if(!stream_->writeByte(0xFFU)){
                        return false;
                    }
                }
                buffer_ = (low_>>SHIFT) & 0xFFU;
                carry_ = 0;
            }else{
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
    AdaptiveRangeDecoder::AdaptiveRangeDecoder()
        :stream_(CPPRCODER_NULL)
        ,umcompressedSize_(0)
        ,outSize_(0)
        ,state_(State_Init)
        ,range_(CPPRCODER_MAXRANGE_DECODE)
        ,low_(0)
    {
    }

    AdaptiveRangeDecoder::~AdaptiveRangeDecoder()
    {
    }

    bool AdaptiveRangeDecoder::initialize(IStream& stream)
    {
        stream_ = &stream;
        umcompressedSize_ = 0;
        outSize_ = 0;
        state_ = State_Init;
        range_ = CPPRCODER_MAXRANGE_DECODE;
        low_ = 0;
        return true;
    }

    Result AdaptiveRangeDecoder::decode(s32 size, const u8* bytes)
    {
        u32 t,count;
        switch(state_){
        case State_Init:
            if(size<8){
                return {Status_Pending, 8};
            }

            umcompressedSize_ = bytes[3];
            umcompressedSize_ = (umcompressedSize_<<8) | bytes[2];
            umcompressedSize_ = (umcompressedSize_<<8) | bytes[1];
            umcompressedSize_ = (umcompressedSize_<<8) | bytes[0];
            bytes += 4;
            size -= 4;

            low_ = bytes[0];
            low_ = (low_<<8) | bytes[1];
            low_ = (low_<<8) | bytes[2];
            low_ = (low_<<8) | bytes[3];
            bytes += 4;
            size -= 4;
            state_ = State_Decode;
            goto STATE_DECODE;

        case State_Decode:
STATE_DECODE:
            for(;;){
                if(!normalize(size, bytes)){
                    return {Status_Pending, umcompressedSize_-outSize_};
                }
                t = range_ / frequencies_.total();
                frequencies_.find(count, lastCode_, low_/t);
                low_ -= t * count;
                range_ = t * frequencies_[lastCode_];

                if(!stream_->writeByte(lastCode_)){
                    return {Status_Pending, umcompressedSize_-outSize_};
                }
                if(umcompressedSize_<=++outSize_){
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

    bool AdaptiveRangeDecoder::normalize(s32& size, const u8*& bytes)
    {
        while(range_<CPPRCODER_MINRANGE){
            if(size<=0){
                return false;
            }
            range_ <<= 8;
            u8 byte = bytes[0];
            --size;
            ++bytes;
            low_ = (low_<<8) + byte;
        }
        return true;
    }
}
#endif

#endif //INC_CPPRCODER_H_
