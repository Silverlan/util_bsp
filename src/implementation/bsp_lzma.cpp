// SPDX-FileCopyrightText: (c) 2024 Silverlan <opensource@pragma-engine.com>
// SPDX-License-Identifier: MIT

module;

#include "LzmaLib.h"

module source_engine.bsp;

namespace source_engine::bsp {
	bool lzma_uncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t *srcLen, const unsigned char *props, size_t propsSize) { return LzmaUncompress(dest, destLen, src, srcLen, props, propsSize) == SZ_OK; }
};
