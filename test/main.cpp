#define CPPRCODER_IMPLEMENTATION
#include "../cpprcoder.h"

#include <chrono>
#include <fstream>
#include <memory>
#include <random>

#include "slz4.h"

#ifdef _MSC_VER
#    include <zlib/zlib.h>
#else
#    include <zlib.h>
#endif

#ifdef _MSC_VER
#    define USE_LZ4
#    include <lz4/lz4.h>
#endif

#define USE_RC
#define USE_ADAPTIVE
#define USE_SLZ4
#define USE_ZLIB

//#undef USE_LZ4

class Timer
{
public:
    void start();
    void stop();
    long long microseconds() const;

private:
    std::chrono::high_resolution_clock::time_point start_;
    std::chrono::high_resolution_clock::time_point end_;
};

void Timer::start()
{
    start_ = std::chrono::high_resolution_clock::now();
}

void Timer::stop()
{
    end_ = std::chrono::high_resolution_clock::now();
}

long long Timer::microseconds() const
{
    return std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count();
}

#ifdef USE_SLZ4
int def_slz4(cpprcoder::MemoryStream& outStream, cpprcoder::u32 srcSize, cpprcoder::u8* src)
{
    slz4::SLZ4Context context;
    cpprcoder::s32 result = slz4::compress(context, outStream.capacity(), reinterpret_cast<slz4::u8*>(&outStream[0]), srcSize, reinterpret_cast<slz4::u8*>(src));
    if(0 <= result) {
        outStream.resize(result);
    }
    return result;
}

int inf_slz4(cpprcoder::MemoryStream& outStream, cpprcoder::u32 srcSize, cpprcoder::u8* src)
{
    cpprcoder::s32 result = slz4::decompress(outStream.capacity(), reinterpret_cast<slz4::u8*>(&outStream[0]), srcSize, reinterpret_cast<slz4::u8*>(src));
    if(0 <= result) {
        outStream.resize(result);
    }
    return result;
}
#endif

#ifdef USE_ZLIB
int def_zlib(cpprcoder::MemoryStream& outStream, cpprcoder::u32 srcSize, cpprcoder::u8* src)
{
    int ret, flush;
    z_stream stream;
    stream.zalloc = NULL;
    stream.zfree = NULL;
    stream.opaque = NULL;
    ret = deflateInit(&stream, Z_DEFAULT_COMPRESSION);

    if(Z_OK != ret) {
        return ret;
    }

    outStream.reserve(srcSize);

    const int Chunk = 16384;
    cpprcoder::u8 out[Chunk];
    cpprcoder::u32 count = 0;
    int outCount = 0;
    do {
        if(srcSize <= count) {
            deflateEnd(&stream);
            return outCount;
        }
        cpprcoder::u32 size = srcSize - count;
        stream.avail_in = size;
        stream.next_in = src + count;
        count += size;
        flush = (srcSize <= count) ? Z_FINISH : Z_NO_FLUSH;
        do {
            stream.avail_out = Chunk;
            stream.next_out = out;
            ret = deflate(&stream, flush);
            int s = Chunk - stream.avail_out;
            outStream.write(s, out);
            outCount += s;
        } while(stream.avail_out == 0);
        CPPRCODER_ASSERT(stream.avail_in <= 0);
    } while(flush != Z_FINISH);
    deflateEnd(&stream);
    return outCount;
}

int inf_zlib(cpprcoder::MemoryStream& outStream, cpprcoder::u32 srcSize, cpprcoder::u8* src)
{
    int ret;
    z_stream stream;
    stream.zalloc = NULL;
    stream.zfree = NULL;
    stream.opaque = NULL;
    ret = inflateInit(&stream);
    if(Z_OK != ret) {
        return ret;
    }

    const int Chunk = 16384;
    cpprcoder::u8 out[Chunk];
    cpprcoder::u32 count = 0;
    int outCount = 0;
    do {
        if(srcSize <= count) {
            inflateEnd(&stream);
            return Z_OK;
        }
        cpprcoder::u32 size = srcSize - count;
        stream.avail_in = size;
        stream.next_in = src + count;
        count += size;

        do {
            stream.avail_out = Chunk;
            stream.next_out = out;
            ret = inflate(&stream, Z_NO_FLUSH);
            switch(ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;
                break;
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                inflateEnd(&stream);
                return ret;
            }
            int s = Chunk - stream.avail_out;
            outStream.write(s, out);
            outCount += s;
        } while(stream.avail_out == 0);
        CPPRCODER_ASSERT(stream.avail_in <= 0);
    } while(ret != Z_STREAM_END);
    inflateEnd(&stream);
    return ret == Z_STREAM_END ? outCount : Z_DATA_ERROR;
}
#endif

#ifdef USE_LZ4
int def_lz4(cpprcoder::MemoryStream& outStream, cpprcoder::u32 srcSize, cpprcoder::u8* src)
{
    cpprcoder::s32 result = LZ4_compress_default(reinterpret_cast<char*>(src), reinterpret_cast<char*>(&outStream[0]), srcSize, outStream.size());
    if(0 <= result) {
        outStream.resize(result);
    }
    return result;
}

int inf_lz4(cpprcoder::MemoryStream& outStream, cpprcoder::u32 srcSize, cpprcoder::u8* src)
{
    cpprcoder::s32 result = LZ4_decompress_safe(reinterpret_cast<const char*>(src), reinterpret_cast<char*>(&outStream[0]), srcSize, outStream.size());
    return result;
}
#endif

#ifdef USE_RC
void run_rangecoder(const char* filepath)
{
    Timer timer;
    long long deflateTime, inflateTime;
    std::ifstream file(filepath, std::ios::binary);
    if(!file.is_open()) {
        return;
    }
    file.seekg(0, std::ios::end);
    cpprcoder::u32 size = static_cast<cpprcoder::u32>(file.tellg());
    file.seekg(0, std::ios::beg);

    cpprcoder::u8* src = new cpprcoder::u8[size];
    file.read(reinterpret_cast<char*>(src), size);
    file.close();

    cpprcoder::MemoryStream encstream(size);
    cpprcoder::MemoryStream decstream(size);

    timer.start();
    cpprcoder::RangeEncoder<> encoder;
    if(!encoder.encode(encstream, size, src)) {
        delete[] src;
        return;
    }
    timer.stop();
    deflateTime = timer.microseconds();

    timer.start();
    if(!encoder.decode(decstream, encstream.size(), encstream.get())) {
        delete[] src;
        return;
    }
    timer.stop();
    inflateTime = timer.microseconds();

    double ratio = (double)encstream.size() / size;
    //printf("|%s|%d|%d|%f|%lld|%lld|\n", filepath, size, encstream.size(), ratio, deflateTime, inflateTime);
    printf("|%s|%f|%lld|%lld|\n", filepath, ratio, deflateTime, inflateTime);
    for(cpprcoder::u32 i = 0; i < size; ++i) {
        if(decstream[i] != src[i]) {
            printf("[%d] %d != %d\n", i, decstream[i], src[i]);
        }
    }
    delete[] src;
}
#endif

#ifdef USE_ADAPTIVE
void run_adaptive(const char* filepath)
{
    Timer timer;
    long long deflateTime, inflateTime;
    std::ifstream file(filepath, std::ios::binary);
    if(!file.is_open()) {
        return;
    }
    file.seekg(0, std::ios::end);
    cpprcoder::u32 size = static_cast<cpprcoder::u32>(file.tellg());
    file.seekg(0, std::ios::beg);

    cpprcoder::u8* src = new cpprcoder::u8[size];
    file.read(reinterpret_cast<char*>(src), size);
    file.close();

    cpprcoder::MemoryStream encstream(size);
    cpprcoder::MemoryStream decstream(size);

    timer.start();
    cpprcoder::AdaptiveRangeEncoder<> encoder;
    if(!encoder.initialize(encstream, size)) {
        delete[] src;
        return;
    }
    cpprcoder::Result result0 = encoder.encode(size, src);
    if(cpprcoder::Status_Success != result0.status_) {
        delete[] src;
        return;
    }
    timer.stop();
    deflateTime = timer.microseconds();

    timer.start();
    cpprcoder::AdaptiveRangeDecoder<> decoder;
    if(!decoder.initialize(decstream)) {
        delete[] src;
        return;
    }
    cpprcoder::Result result1 = decoder.decode(encstream.size(), &encstream[0]);
    if(cpprcoder::Status_Success != result1.status_) {
        delete[] src;
        return;
    }
    timer.stop();
    inflateTime = timer.microseconds();

    double ratio = (double)encstream.size() / size;
    //printf("|%s|%d|%d|%f|%lld|%lld|\n", filepath, size, encstream.size(), ratio, deflateTime, inflateTime);
    printf("|%s|%f|%lld|%lld|\n", filepath, ratio, deflateTime, inflateTime);
    for(cpprcoder::u32 i = 0; i < size; ++i) {
        if(decstream[i] != src[i]) {
            printf("[%d] %d != %d\n", i, decstream[i], src[i]);
        }
    }
    delete[] src;
}
#endif

#ifdef USE_SLZ4
void run_slz4(const char* filepath)
{
    Timer timer;
    long long deflateTime, inflateTime;

    std::ifstream file(filepath, std::ios::binary);
    if(!file.is_open()) {
        return;
    }
    file.seekg(0, std::ios::end);
    cpprcoder::u32 size = static_cast<cpprcoder::u32>(file.tellg());
    file.seekg(0, std::ios::beg);

    cpprcoder::u8* src = new cpprcoder::u8[size];
    file.read(reinterpret_cast<char*>(src), size);
    file.close();

    cpprcoder::MemoryStream encstream(size * 2);
    encstream.resize(size * 2);
    cpprcoder::MemoryStream decstream(size);
    decstream.resize(size);

    timer.start();
    if(def_slz4(encstream, size, src) < 0) {
        printf("Error in def_slz4\n");
        delete[] src;
        return;
    }
    timer.stop();
    deflateTime = timer.microseconds();

    timer.start();
    if(inf_slz4(decstream, encstream.size(), &encstream[0]) < 0) {
        printf("Error in inf_slz4\n");
        delete[] src;
        return;
    }
    timer.stop();
    inflateTime = timer.microseconds();

    for(int i=0; i<decstream.size(); ++i){
        if(src[i] != decstream[i]){
            printf("Error: %d != %d\n", src[i], decstream[i]);
            return;
        }
    }
    double ratio = (double)encstream.size() / size;
    //printf("|%s|%d|%d|%f|%lld|%lld|\n", filepath, size, encstream.size(), ratio, deflateTime, inflateTime);
    printf("|%s|%f|%lld|%lld|\n", filepath, ratio, deflateTime, inflateTime);
    delete[] src;
}
#endif

#ifdef USE_ZLIB
void run_zlib(const char* filepath)
{
    Timer timer;
    long long deflateTime, inflateTime;

    std::ifstream file(filepath, std::ios::binary);
    if(!file.is_open()) {
        return;
    }
    file.seekg(0, std::ios::end);
    cpprcoder::u32 size = static_cast<cpprcoder::u32>(file.tellg());
    file.seekg(0, std::ios::beg);

    cpprcoder::u8* src = new cpprcoder::u8[size];
    file.read(reinterpret_cast<char*>(src), size);
    file.close();

    cpprcoder::MemoryStream encstream(size);
    cpprcoder::MemoryStream decstream(size);

    timer.start();
    if(def_zlib(encstream, size, src) < 0) {
        delete[] src;
        return;
    }
    timer.stop();
    deflateTime = timer.microseconds();

    timer.start();
    if(inf_zlib(decstream, encstream.size(), &encstream[0]) < 0) {
        delete[] src;
        return;
    }
    timer.stop();
    inflateTime = timer.microseconds();

    double ratio = (double)encstream.size() / size;
    //printf("|%s|%d|%d|%f|%lld|%lld|\n", filepath, size, encstream.size(), ratio, deflateTime, inflateTime);
    printf("|%s|%f|%lld|%lld|\n", filepath, ratio, deflateTime, inflateTime);
    delete[] src;
}
#endif

#ifdef USE_LZ4
void run_lz4(const char* filepath)
{
    Timer timer;
    long long deflateTime, inflateTime;

    std::ifstream file(filepath, std::ios::binary);
    if(!file.is_open()) {
        return;
    }
    file.seekg(0, std::ios::end);
    cpprcoder::u32 size = static_cast<cpprcoder::u32>(file.tellg());
    file.seekg(0, std::ios::beg);

    cpprcoder::u8* src = new cpprcoder::u8[size];
    file.read(reinterpret_cast<char*>(src), size);
    file.close();

    cpprcoder::MemoryStream encstream(size * 2);
    encstream.resize(size * 2);
    cpprcoder::MemoryStream decstream(size);
    decstream.resize(size);

    timer.start();
    if(def_lz4(encstream, size, src) < 0) {
        delete[] src;
        return;
    }
    timer.stop();
    deflateTime = timer.microseconds();
    
    timer.start();
    if(inf_lz4(decstream, encstream.size(), &encstream[0]) < 0) {
        delete[] src;
        return;
    }
    timer.stop();
    inflateTime = timer.microseconds();

    double ratio = (double)encstream.size() / size;
    //printf("|%s|%d|%d|%f|%lld|%lld|\n", filepath, size, encstream.size(), ratio, deflateTime, inflateTime);
    printf("|%s|%f|%lld|%lld|\n", filepath, ratio, deflateTime, inflateTime);
    delete[] src;
}
#endif

#ifdef USE_RC
bool test_rangecoder()
{
    std::mt19937 mt;
    std::random_device rand;
    mt.seed(rand());
    static const int Size = 127;
    std::unique_ptr<cpprcoder::u8[]> src(new cpprcoder::u8[Size]);
    for(int i = 0; i < Size; ++i) {
        src[i] = static_cast<cpprcoder::u8>(mt() & 0x0FU);
    }

    cpprcoder::RangeEncoder<> encoder;
    cpprcoder::MemoryStream encstream(Size);
    if(!encoder.encode(encstream, Size, src.get())) {
        return false;
    }
    cpprcoder::MemoryStream decstream(Size);
    if(!encoder.decode(decstream, encstream.size(), encstream.get())) {
        return false;
    }
    for(int i = 0; i < Size; ++i) {
        if(decstream[i] != src[i]) {
            printf("[%d] %d != %d\n", i, decstream[i], src[i]);
            return false;
        }
    }
    return true;
}
#endif

#ifdef USE_ADAPTIVE
bool test_adaptive()
{
    std::mt19937 mt;
    std::random_device rand;
    mt.seed(rand());
    static const int Size = 128 * 1024 * 1024;
    std::unique_ptr<cpprcoder::u8[]> src(new cpprcoder::u8[Size]);
    for(int i = 0; i < Size; ++i) {
        src[i] = static_cast<cpprcoder::u8>(mt() & 0xFFU);
    }

    cpprcoder::AdaptiveRangeEncoder<> encoder;
    cpprcoder::MemoryStream encstream(Size);
    if(!encoder.initialize(encstream, Size)) {
        return false;
    }
    cpprcoder::Result result0 = encoder.encode(Size, src.get());
    if(cpprcoder::Status_Success != result0.status_) {
        return false;
    }
    cpprcoder::AdaptiveRangeDecoder<> decoder;
    cpprcoder::MemoryStream decstream(Size);
    if(!decoder.initialize(decstream)) {
        return false;
    }
    cpprcoder::Result result1 = decoder.decode(encstream.size(), &encstream[0]);
    if(cpprcoder::Status_Success != result1.status_) {
        return false;
    }
    for(int i = 0; i < Size; ++i) {
        if(decstream[i] != src[i]) {
            printf("[%d] %d != %d\n", i, decstream[i], src[i]);
            return false;
        }
    }
    return true;
}
#endif

int main(int /*argc*/, char** /*argv*/)
{
    //test_rangecoder();
    //test_adaptive();

    const char* files[] =
        {
            "../cantrbry/alice29.txt",
            "../cantrbry/asyoulik.txt",
            "../cantrbry/cp.html",
            "../cantrbry/fields.c",
            "../cantrbry/grammar.lsp",
            "../cantrbry/kennedy.xls",
            "../cantrbry/lcet10.txt",
            "../cantrbry/plrabn12.txt",
            "../cantrbry/ptt5",
            "../cantrbry/sum",
            "../cantrbry/xargs.1",
        };
    static const int numFiles = sizeof(files) / (sizeof(files[0]));

    #ifdef USE_RC
    printf("Range Coder\n");
    printf("-------------------------------------------\n");
    for(int i = 0; i < numFiles; ++i) {
        run_rangecoder(files[i]);
    }
    #endif

    #ifdef USE_ADAPTIVE
    printf("Adaptive Range Coder\n");
    printf("-------------------------------------------\n");
    for(int i = 0; i < numFiles; ++i) {
        run_adaptive(files[i]);
    }
    #endif

    #ifdef USE_SLZ4
    printf("SLZ4\n");
    printf("-------------------------------------------\n");
    for(int i = 0; i < numFiles; ++i) {
        run_slz4(files[i]);
    }
    #endif

#ifdef USE_ZLIB
    printf("ZLib\n");
    printf("-------------------------------------------\n");
    for(int i = 0; i < numFiles; ++i) {
        run_zlib(files[i]);
    }
#endif

#ifdef USE_LZ4
    printf("LZ4\n");
    printf("-------------------------------------------\n");
    for(int i = 0; i < numFiles; ++i) {
        run_lz4(files[i]);
    }
#endif
    return 0;
}
