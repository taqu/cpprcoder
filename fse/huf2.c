#include "huf2.h"
#include "huf.h"
#include <assert.h>
#include <string.h>
#ifdef _DEBUG
#    include <stdio.h>
#endif

#if defined(__cplusplus)
extern "C"
{
#endif
    uint32_t HUF_compressBlocksBound(uint32_t size)
    {
        assert(0 < size);
        uint32_t num_blocks = (size + ((1UL << 17) - 1)) >> 17;
        uint32_t block_bound = HUF_compressBound(HUF_BLOCKSIZE_MAX);
        uint32_t upper_bits = (2 * num_blocks + 0x07UL) >> 3;
        return num_blocks * (sizeof(HUF_BlockHeader) + block_bound) + sizeof(HUF_Header) + upper_bits;
    }

    inline static void HUF_set_header(uint32_t index, HUF_BlockHeader* headers, uint8_t* upper_bits, uint32_t original, uint32_t compressed)
    {
        assert(0<original && original <= (1UL << 17));
        assert(0<compressed && compressed <= (1UL << 17));
        --original;
        --compressed;
        headers[index].size_ = original | (compressed << 17);

        index <<= 1;
        uint32_t i = index >> 3;
        uint32_t offset = index - (i << 3);
        upper_bits[i] |= (compressed >> 15) << offset;
    }

    inline static void HUF_get_header(uint32_t index, uint32_t* original, uint32_t* compressed, const HUF_BlockHeader* headers, const uint8_t* upper_bits)
    {
        *original = headers[index].size_ & ((1UL<<17)-1);
        *original += 1;
        *compressed = headers[index].size_ >> 17;

        index <<= 1;
        uint32_t i = index >> 3;
        uint32_t offset = index - (i << 3);
        uint32_t bits = (upper_bits[i]>>offset) & 0x3UL;
        *compressed |= bits<<15;
        *compressed += 1;
    }

    uint32_t HUF_compressBlocks(void* dst, uint32_t dstCapacity, const void* src, uint32_t srcSize)
    {
        assert(NULL != dst);
        assert(0 < dstCapacity);
        assert(NULL != src);
        assert(0 < srcSize);
        uint32_t num_blocks = (uint32_t)((srcSize + ((1UL << 17) - 1)) >> 17);
        uint8_t* d = (uint8_t*)dst;
        HUF_Header* h = (HUF_Header*)d;
        const uint8_t* s = (const uint8_t*)src;
        {
            HUF_Header header;
            header.size_ = (uint32_t)srcSize;
            memcpy(d, &header, sizeof(HUF_Header));
            d += sizeof(HUF_Header);
        }
        HUF_BlockHeader* headers = (HUF_BlockHeader*)d;
        d += sizeof(HUF_BlockHeader) * num_blocks;
        uint8_t* upper_bits = d;
        {
            uint32_t bytes = (2 * num_blocks + 0x07UL) >> 3;
            memset(upper_bits, 0, bytes);
            d += bytes;
        }
        uint32_t block_bound = HUF_compressBound(HUF_BLOCKSIZE_MAX);
        for(uint32_t i = 0; i < num_blocks; ++i) {
            uint32_t size = (HUF_BLOCKSIZE_MAX < srcSize) ? HUF_BLOCKSIZE_MAX : srcSize;
            assert(0 < size);
            uint32_t r = HUF_compress(d, block_bound, s, size);
            if(r <= 0) {
                r = size;
                memcpy(d, s, size);
            } else if(HUF_isError(r)) {
#ifdef _DEBUG
                printf("%s\n", HUF_getErrorName(r));
#endif
                return 0;
            }
            d += r;
            s += size;
            srcSize -= size;

            assert(size <= (1UL << 17));
            assert(r <= (1UL << 17));
            HUF_set_header(i, headers, upper_bits, size, r);
        }
        return (uint32_t)(d - (uint8_t*)dst);
    }

    uint32_t HUF_decompressBlocks(void* dst, uint32_t originalSize, const void* src, uint32_t srcSize)
    {
        assert(NULL != dst);
        assert(0 < originalSize);
        assert(NULL != src);
        assert(sizeof(HUF_Header) <= srcSize);
        HUF_Header header = {0};
        memcpy(&header, src, sizeof(HUF_Header));
        if(header.size_ != originalSize) {
            return 0;
        }
        uint32_t num_blocks = (uint32_t)((originalSize + ((1UL << 17) - 1)) >> 17);

        uint8_t* d = (uint8_t*)dst;
        const uint8_t* s = (const uint8_t*)src + sizeof(HUF_Header);
        const uint8_t* send = (const uint8_t*)src + srcSize;

        const HUF_BlockHeader* headers = (const HUF_BlockHeader*)s;
        s += sizeof(HUF_BlockHeader) * num_blocks;
        const uint8_t* upper_bits = s;
        {
            uint32_t bytes = (2 * num_blocks + 0x07UL) >> 3;
            s += bytes;
        }
        for(uint32_t i=0; i<num_blocks; ++i){
            uint32_t osize;
            uint32_t csize;
            HUF_get_header(i, &osize, &csize, headers, upper_bits);
            uint32_t r = HUF_decompress(d, osize, s, csize);
            if(HUF_isError(r)) {
#ifdef _DEBUG
                printf("%s\n", HUF_getErrorName(r));
#endif
                return 0;
            }
            s += csize;
            d += osize;
        }
        return (uint32_t)(d - (uint8_t*)dst) == originalSize ? originalSize : 0;
    }

#if defined(__cplusplus)
}
#endif
