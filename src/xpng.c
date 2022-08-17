/* gcc xpng.c -o xpng -lpng -lm*/
/* Requires libpng and zlib */

#include "xpng.h"

int16_t byteclamp16_t(int16_t c)
{
    int16_t buff[3] = {0, c, 255};
    return buff[ (int)(c > 0) + (int)(c > 255) ];
}

/* Returns pointer to a certain address in a 24 bit datafield */
png_bytep pix(int32_t i, int32_t j, uint8_t c, int32_t w, png_bytep data)
{
    return &(data[i*w*3+j*3+c]);
}

/* Returns the difference between two bytes for error calculation */
uint8_t pix_err(png_byte x1, png_byte x2)
{
    uint8_t a1 = x1, a2 = x2;
    return (a1 > a2) ? (a1 - a2) : (a2 - a1);
}

/* Average predictor for PNG */
png_byte pix_avg(int32_t i, int32_t j, uint8_t c, int32_t w, png_bytep data)
{
    int16_t p;
    png_byte left = (j > 0) ? *pix(i,j-1,c,w,data) : 0;
    png_byte up = (i > 0) ? *pix(i-1,j,c,w,data) : 0;

    p = up;
    p += left;
    p /= 2;
    return p;
}

/* Top predictor for PNG - unused */
png_byte pix_top(int32_t i, int32_t j, uint8_t c, int32_t w, png_bytep data)
{
    return ((i > 0) ? *pix(i-1,j,c,w,data) : 0);
}

/* Left predictor for PNG - used in first row */
png_byte pix_sub(int32_t i, int32_t j, uint8_t c, int32_t w, png_bytep data)
{
    return ((j > 0) ? *pix(i,j-1,c,w,data) : 0);
}

/* Paeth predictor for PNG - unused */
png_byte pix_paeth(int32_t i, int32_t j, uint8_t c, int32_t w, png_bytep data)
{
    uint8_t dist = 0;
    int16_t p;

    png_byte left = (j > 0) ? *pix(i,j-1,c,w,data) : 0;
    png_byte up = (i > 0) ? *pix(i-1,j,c,w,data) : 0;
    png_byte topleft = (i > 0 && j > 0) ? *pix(i-1,j-1,c,w,data) : 0;

    p = left;
    p += up;
    p -= topleft;
    p = byteclamp16_t(p);

    png_byte base;

    base = left;
    dist = pix_err((png_byte) p, left);

    if (pix_err((png_byte) p, up) < dist)
    {
        dist = pix_err((png_byte) p, up);
        base = up;
    }

    if (pix_err((png_byte) p, topleft) < dist)
    {
        dist = pix_err((png_byte) p, topleft);
        base = topleft;
    }

    return base;
}

/* Blured images: Ib = M(I, w, w); w - window (+-) radius */
png_byte pix_blur(int32_t i, int32_t j, uint8_t c, int32_t h, int32_t w, int32_t radius, png_bytep data)
{
    int16_t M = 0;
    int32_t ii, jj, N = 0;
    int64_t p0, p, X = 0;

    p0 = *pix(i, j, c, w, data);
    for (ii = i - radius; ii <= i + radius; ii++)
        if ((ii >= 0) && (ii < h))
            for (jj = j - radius; jj <= j + radius; jj++)
                if ((jj >= 0) && (jj < w))
                {
                    p = *pix(ii, jj, c, w, data);
                    X += p;
                    N++;
                }
    X -= p0;
    N--;
    M = (N > 0) ? (X / N) : p0;
    M = byteclamp16_t(M);
    return M;
}

/* Calculates noisiness - used for adaptive quantization
 * i = I - Ib; Ib - blur area
 * e = Sum(i*i) / n - M(i)*M(i);
 *  */
png_byte pix_calc_noise(int32_t i, int32_t j, uint8_t c, int32_t h, int32_t w, int32_t radius, png_bytep data, png_bytep blur)
{
    int16_t err = 0;
    int32_t ii, jj, N = 0;
    int64_t p, b, X = 0, X2 = 0;

    for (ii = i - radius; ii <= i + radius; ii++)
        for (jj = j - radius; jj <= j + radius; jj++)
            if (ii >= 0 && ii < h && jj >= 0 && jj < w)
            {
                p = *pix(ii, jj, c, w, data);
                b = *pix(ii, jj, c, w, blur);
                p -= b;
                X += p;
                X2 += p * p;
                N++;
            }
    if (N > 0)
    {
        X /= N;
        X2 /= N;
        err = sqrt(X2 - X * X);
        err = byteclamp16_t(err);
    }
    return err;
}

png_bytep xpng(png_bytep buffer, int32_t w, int32_t h, size_t pngsize, uint8_t clevel, int32_t radius)
{
    uint8_t c, d, delta;
    uint16_t jumpsize, qlevel;
    int32_t i, j, refdist;
    uint64_t f, freq[256] = {};
    size_t rd;
    png_byte base, target, ltarget, rtarget, approx, a;

    png_bytep buffer2; /* Output image buffer */
    png_bytep diff; /* Residuals from predictor */
    png_bytep noise; /* Local noisiness */

    qlevel = clevel;
    qlevel = sqrt(qlevel * 256 * 256);

    buffer2 = malloc(pngsize);
    diff = malloc(pngsize);
    noise = malloc(pngsize);

    if (buffer == NULL || buffer2 == NULL || diff == NULL || noise == NULL)
    {
        fprintf(stderr, "xpng: error: insufficient memory\n");
        return NULL;
    }

    memset(buffer2, 0, pngsize);

    /* Calculate local noisiness */
    rd = 0;
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++)
            for (c = 0; c < 3; c++)
            {
                diff[rd] = pix_blur(i,j,c,h,w,radius,buffer);
                rd++;
            }
    rd = 0;
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++)
            for (c = 0; c < 3; c++)
            {
                noise[rd] = pix_calc_noise(i,j,c,h,w,radius,buffer,diff);
                rd++;
            }

    /* Specify the preference levels of each quantization */
    for (jumpsize = 1; jumpsize <= 128; jumpsize *= 2)
        for (i = -128; i <= 128; i += jumpsize)
            freq[(uint8_t)i] = jumpsize;

    /* Go through image and adaptively quantize for noise masking */
    for (i = 0; i < h; i++)
    {
        refdist = 1;
        for (j = 0; j < w; j++)
        {
            for (c = 0; c < 3; c++)
            {
                /* Left predictor for first row, avg. for others */
                if (i == 0 && j > 0)
                {
                    target = *pix(i,j,c,w,buffer);
                    base = pix_sub(i,j,c,w,buffer2);
                }
                else
                {
                    target = *pix(i,j,c,w,buffer);
                    base = pix_avg(i,j,c,w,buffer2);
                }

                /* Calculate tolerance using comp. level and local noise */
                delta = (qlevel + (uint32_t) *pix(i,j,c,w,noise) * 256 * qlevel / (2047 + qlevel)) / 256;

                /* Calculate interval based on tolerance */
                ltarget = (target > delta) ? (target - delta) : 0;
                rtarget = ((255 - target) > delta) ? (target + delta) : 255;

                /* Find best quantization within allowable tolerance */
                *pix(i,j,c,w,buffer2) = target;
                approx = target;
                f = 0;
                a = ltarget;
                do
                {
                    d = a-base;
                    if (freq[d] > f)
                    {
                        approx = a;
                        f = freq[d];
                    }
                }
                while (a++ != rtarget);

                *pix(i,j,c,w,diff) = approx - base;
                /*freq[*pix(i,j,c,diff)]++;*/
                *pix(i,j,c,w,buffer2) = approx;

                /* Use recent value for runlength encoding if possible */
                if (pix_err(target,base+*(pix(i,j,c,w,diff)-refdist)) < delta)
                {
                    *pix(i,j,c,w,diff) = *(pix(i,j,c,w,diff)-refdist);
                    *pix(i,j,c,w,buffer2) = base + *pix(i,j,c,w,diff);
                    continue;
                }
                for (rd = 1; rd <= j*3+c && rd <= 4; rd++)
                {
                    if (pix_err(target,base+*(pix(i,j,c,w,diff)-rd)) < delta)
                    {
                        refdist = rd;
                        *pix(i,j,c,w,diff) = *(pix(i,j,c,w,diff)-rd);
                        *pix(i,j,c,w,buffer2) = base + *pix(i,j,c,w,diff);
                        break;
                    }
                }
            }
        }
    }
    free(diff);
    free(noise);

    return buffer2;
}

int main(int argc, const char **argv)
{
    uint8_t clevel = 32;
    int32_t h, w, radius = 2;
    size_t pngsize;

    png_image image; /* The control structure used by libpng */
    png_bytep buffer; /* Buffer for original image */
    png_bytep buffer2; /* Output image buffer */

    if (argc < 3)
    {
        fprintf(stderr, "XPNG version %s (lossy PNG encoder)\n", XPNG_VERSION);
        fprintf(stderr, "Usage: %s input.png output.png [level=%d, {0..255}] [radius=%d, {1..16}]\n", argv[0], clevel, radius);
        exit(1);
    }
    if (argc > 3) clevel = atoi(argv[3]);
    if (clevel > 255)
    {
        fprintf(stderr, "xpng: error: level=%d ! {0..255}]\n", clevel);
        exit(1);
    }
    if (argc > 4) radius = atoi(argv[4]);
    if (radius < 1 || radius > 16)
    {
        fprintf(stderr, "xpng: error: radius=%d ! {1..16}]\n", radius);
        exit(1);
    }

    /* Initialize png_image structure */
    memset(&image, 0, (sizeof image));
    image.version = PNG_IMAGE_VERSION;

    /* The first argument is the input file */
    if (png_image_begin_read_from_file(&image, argv[1]) == 0)
    {
        fprintf(stderr, "xpng: error: %s\n", image.message);
        exit (1);
    }

    image.format = PNG_FORMAT_RGB;
    pngsize = PNG_IMAGE_SIZE(image);
    buffer = malloc(pngsize);

    if (buffer == NULL)
    {
        fprintf(stderr, "xpng: error: insufficient memory\n");
        exit(1);
    }
    memset(buffer, 0, PNG_IMAGE_SIZE(image));

    if (buffer == NULL ||
            png_image_finish_read(&image, NULL/*bg*/, buffer,
                                  0/*row_stride*/, NULL/*colormap*/) == 0)
    {
        fprintf(stderr, "xpng: error: %s\n", image.message);
        exit (1);
    }

    h = image.height;
    w = image.width;

    buffer2 = xpng(buffer, w, h, pngsize, clevel, radius);

    free(buffer);

    if (buffer2 != NULL)
    {
        if (png_image_write_to_file(&image, argv[2], 0/*already_8bit*/,
                                    buffer2, 0/*row_stride*/, NULL/*colormap*/) == 0)
        {
            fprintf(stderr, "xpng: error: %s\n", image.message);
            exit (1);
        }
        free(buffer2);
    }

    exit(0);
}
