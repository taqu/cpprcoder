#ifndef INC_HUF2_H_
#define INC_HUF2_H_
#include <stdint.h>
#include <stddef.h>

#if defined (__cplusplus)
extern "C" {
#endif
	typedef struct HUF_Header_t
    {
		uint32_t size_;
    } HUF_Header;

	typedef struct HUF_BlockHeader_t
    {
		uint32_t size_;
    } HUF_BlockHeader;


	uint32_t HUF_compressBlocksBound(uint32_t size);
	uint32_t HUF_compressBlocks(void* dst, uint32_t dstCapacity, const void* src, uint32_t srcSize);
	uint32_t HUF_decompressBlocks(void* dst, uint32_t originalSize, const void* src, uint32_t srcSize);

#if defined (__cplusplus)
}
#endif
#endif //INC_HUF2_H_
