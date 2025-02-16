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

	void HUF_set_header(uint32_t index, HUF_BlockHeader* headers, uint8_t* upper_bits, uint32_t original, uint32_t compressed);
	void HUF_get_header(uint32_t index, uint32_t* original, uint32_t* compressed, const HUF_BlockHeader* headers, const uint8_t* upper_bits);

	uint32_t HUF_compressBlocksBound(uint32_t size);
	uint32_t HUF_compressBlocks(void* dst, uint32_t dstCapacity, const void* src, uint32_t srcSize);
	uint32_t HUF_decompressBlocks(void* dst, uint32_t originalSize, const void* src, uint32_t srcSize);

	uint32_t HUF_decompressBlocksWorkSize(uint32_t size);
	uint32_t HUF_decompressBlocksInplace(void* dst, uint32_t originalSize, const void* src, uint32_t srcSize, void* work);
#if defined (__cplusplus)
}
#endif
#endif //INC_HUF2_H_
