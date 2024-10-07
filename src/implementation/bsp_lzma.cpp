/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

module;

#include "LzmaLib.h"

module source_engine.bsp;

namespace source_engine::bsp {
	bool lzma_uncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t *srcLen, const unsigned char *props, size_t propsSize) { return LzmaUncompress(dest, destLen, src, srcLen, props, propsSize) == SZ_OK; }
};
