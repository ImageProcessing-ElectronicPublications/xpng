/* gcc xpng.c -o xpng -lpng -lz -lm*/
/* Requires libpng and zlib */

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <zlib.h>
#include <png.h>

#ifndef __XPNG_H
#define __XPNG_H

#define XPNG_VERSION "1.5"
#define XPNG_BPP 3

void xpng(png_bytep buffer, int32_t w, int32_t h, size_t pngsize, uint8_t clevel, int32_t radius);

#endif /* __XPNG_H */
