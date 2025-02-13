#include "huf2.h"
#include "huf.h"
#include <assert.h>
#include <string.h>

#if defined(__cplusplus)
extern "C"
{
#endif
    size_t HUF_compressBlocksBound(size_t size)
    {
        assert(0 < size);
        size_t num_blocks = (size + ((1UL << 17) - 1)) >> 17;
        size_t block_bound = HUF_compressBound(HUF_BLOCKSIZE_MAX);
        return num_blocks * (sizeof(HUF_BlockHeader) + block_bound) + sizeof(HUF_Header);
    }

    size_t HUF_compressBlocks(void* dst, size_t dstCapacity, const void* src, size_t srcSize)
    {
        assert(NULL != dst);
        assert(0 < dstCapacity);
        assert(NULL != src);
        assert(0 < srcSize);
        uint32_t num_blocks = (uint32_t)((srcSize + ((1UL << 17) - 1)) >> 17);
        uint8_t* d = (uint8_t*)dst;
        const uint8_t* s = (const uint8_t*)src;
        {
            HUF_Header header;
            header.size_ = (uint32_t)srcSize;
            memcpy(d, &header, sizeof(HUF_Header));
            d += sizeof(HUF_Header);
        }
        size_t block_bound = HUF_compressBound(HUF_BLOCKSIZE_MAX);
        for(uint32_t i = 0; i < num_blocks; ++i) {
            size_t size = (HUF_BLOCKSIZE_MAX < srcSize) ? HUF_BLOCKSIZE_MAX : srcSize;
            assert(0 < size);
            uint8_t* bheader = d;
            d += sizeof(HUF_BlockHeader);
            size_t r = HUF_compress(d, block_bound, s, size);
            if(r <= 0) {
                return 0;
            }
            d += r;
            s += size;

            size -= 1;
            r -= 1;
            assert(size < (1UL << 17));
            assert(r < (1UL << 17));
            uint32_t lower = (uint32_t)(((r & ((1UL << 15) - 1)) << 17) | size);
            uint8_t upper = (uint8_t)(r >> 15UL);
            memcpy(bheader, &lower, sizeof(uint32_t));
            bheader[4] = upper;
        }
        return (size_t)(d - (uint8_t*)dst);
    }

    size_t HUF_decompressBlocks(void* dst, size_t originalSize, const void* src, size_t srcSize)
    {
        assert(NULL != dst);
        assert(0 < originalSize);
        assert(NULL != src);
        assert(sizeof(HUF_Header) <= srcSize);
        HUF_Header header = {0};
        memcpy(&header, src, sizeof(HUF_Header));
        if(header.size_<=0 || header.size_<originalSize){
            return 0;
        }
        uint8_t* d = (uint8_t*)dst;
        const uint8_t* s = (const uint8_t*)src + sizeof(HUF_Header);
        const uint8_t* send = (const uint8_t*)src + srcSize;
        while(s<send){
            uint32_t osize;
            memcpy(&osize, s, 3);
            osize &= ((1UL<<17)-1);
            uint32_t csize;
            memcpy(&csize, &s[2], 3);
            csize >>= 2;
            s += sizeof(HUF_BlockHeader);
            size_t r = HUF_decompress(d, osize, s, csize);
            if(r != osize){
                return 0;
            }
            s += csize;
            d += osize;
        }
        return (size_t)(d - (uint8_t*)dst) == originalSize? originalSize : 0;
    }
#if defined(__cplusplus)
}
#endif
