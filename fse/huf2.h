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
		uint8_t value_[5];
    } HUF_BlockHeader;

	size_t HUF_compressBlocksBound(size_t size);
	size_t HUF_compressBlocks(void* dst, size_t dstCapacity, const void* src, size_t srcSize);
	size_t HUF_decompressBlocks(void* dst,  size_t originalSize, const void* src, size_t srcSize);
#if defined (__cplusplus)
}
#endif
#endif //INC_HUF2_H_
