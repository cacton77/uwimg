#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    for (int c = 0; c < im.c; c++)
    {
        float sum = 0.0;
        for (int y = 0; y < im.h; y++)
        {
            for (int x = 0; x < im.w; x++)
            {
                float v = get_pixel(im, x, y, c);
                sum += v;
            }
        }
        sum = sum == 0.0 ? 1.0 : sum;
        for (int y = 0; y < im.h; y++)
        {
            for (int x = 0; x < im.w; x++)
            {
                float v = get_pixel(im, x, y, c);
                set_pixel(im, x, y, c, v / sum);
            }
        }
    }
}

image make_box_filter(int w)
{
    image bf = make_image(w, w, 1);
    for (int y = 0; y < bf.h; y++)
    {
        for (int x = 0; x < bf.w; x++)
        {
            set_pixel(bf, x, y, 0, 1.0);
        }
    }
    l1_normalize(bf);
    return bf;
}

image convolve_image(image im, image filter, int preserve)
{
    assert(filter.c == 1 || filter.c == im.c);
    if (filter.c == im.c)
    {
        image im_conv = make_image(im.w, im.h, im.c);
        int vx = (filter.w - 1) / 2;
        int vy = (filter.h - 1) / 2;
        float qxy, fuv, imuv;
        for (int c = 0; c < im.c; c++)
        {
            for (int y = 0; y < im.h; y++)
            {
                for (int x = 0; x < im.w; x++)
                {
                    qxy = 0.0;
                    for (int v = -vy; v <= vy; v++)
                    {
                        for (int u = -vx; u <= vx; u++)
                        {
                            fuv = get_pixel(filter, vx + u, vy + v, c);
                            imuv = get_pixel(im, x + u, y + v, c);
                            qxy = qxy + fuv * imuv;
                        }
                    }
                    set_pixel(im_conv, x, y, c, qxy);
                }
            }
        }
        return im_conv;
    }
    else if (preserve == 1)
    {
        image im_conv = make_image(im.w, im.h, im.c);
        int vx = (filter.w - 1) / 2;
        int vy = (filter.h - 1) / 2;
        float qxy, fuv, imuv;
        for (int c = 0; c < im.c; c++)
        {
            for (int y = 0; y < im.h; y++)
            {
                for (int x = 0; x < im.w; x++)
                {
                    qxy = 0.0;
                    for (int v = -vy; v <= vy; v++)
                    {
                        for (int u = -vx; u <= vx; u++)
                        {
                            fuv = get_pixel(filter, vx + u, vy + v, 0);
                            imuv = get_pixel(im, x + u, y + v, c);
                            qxy = qxy + fuv * imuv;
                        }
                    }
                    set_pixel(im_conv, x, y, c, qxy);
                }
            }
        }
        return im_conv;
    }
    else
    {
        image im_conv = make_image(im.w, im.h, 1);
        int vx = (filter.w - 1) / 2;
        int vy = (filter.h - 1) / 2;
        float qxy_sum, qxy, fuv, imuv;
        for (int y = 0; y < im.h; y++)
        {
            for (int x = 0; x < im.w; x++)
            {
                qxy_sum = 0.0;
                for (int c = 0; c < im.c; c++)
                {
                    qxy = 0.0;
                    for (int v = -vy; v <= vy; v++)
                    {
                        for (int u = -vx; u <= vx; u++)
                        {
                            fuv = get_pixel(filter, vx + u, vy + v, 0);
                            imuv = get_pixel(im, x + u, y + v, c);
                            qxy = qxy + fuv * imuv;
                        }
                    }
                    qxy_sum += qxy;
                }
                set_pixel(im_conv, x, y, 0, qxy_sum);
            }
        }
        return im_conv;
    }
}

image make_highpass_filter()
{
    image hp = make_image(3, 3, 1);
    set_pixel(hp, 0, 0, 0, 0.0);
    set_pixel(hp, 1, 0, 0, -1.0);
    set_pixel(hp, 2, 0, 0, 0.0);
    set_pixel(hp, 0, 1, 0, -1.0);
    set_pixel(hp, 1, 1, 0, 4.0);
    set_pixel(hp, 2, 1, 0, -1.0);
    set_pixel(hp, 0, 2, 0, 0.0);
    set_pixel(hp, 1, 2, 0, -1.0);
    set_pixel(hp, 2, 2, 0, 0.0);
    return hp;
}

image make_sharpen_filter()
{
    image shrp = make_image(3, 3, 1);
    set_pixel(shrp, 0, 0, 0, 0.0);
    set_pixel(shrp, 1, 0, 0, -1.0);
    set_pixel(shrp, 2, 0, 0, 0.0);
    set_pixel(shrp, 0, 1, 0, -1.0);
    set_pixel(shrp, 1, 1, 0, 5.0);
    set_pixel(shrp, 2, 1, 0, -1.0);
    set_pixel(shrp, 0, 2, 0, 0.0);
    set_pixel(shrp, 1, 2, 0, -1.0);
    set_pixel(shrp, 2, 2, 0, 0.0);
    return shrp;
}

image make_emboss_filter()
{
    image emb = make_image(3, 3, 1);
    set_pixel(emb, 0, 0, 0, -2.0);
    set_pixel(emb, 1, 0, 0, -1.0);
    set_pixel(emb, 2, 0, 0, 0.0);
    set_pixel(emb, 0, 1, 0, -1.0);
    set_pixel(emb, 1, 1, 0, 1.0);
    set_pixel(emb, 2, 1, 0, 1.0);
    set_pixel(emb, 0, 2, 0, 0.0);
    set_pixel(emb, 1, 2, 0, 1.0);
    set_pixel(emb, 2, 2, 0, 2.0);
    return emb;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: We should use preserve with the sharpen and emboss filters and not use preserve with the highpass filter. If we don't preserve the rgb values for the sharpen filter, we end up with the same output as the highpass filter (leaving only grayscale information about high frequency content of the image), thus rgb data is important for this image process. Without preserving color data for emboss, the output loses a large amount of information about the image, and the "embossing" effect is lost.

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: After applying all of the above filters, clamping must be performed to keep the pixel values within the range of 0-1.

image make_gaussian_filter(float sigma)
{
    int s = (int)ceil(6 * sigma);
    int w = s % 2 == 1 ? s : s + 1;
    image gf = make_image(w, w, 1);
    int cx = (w - 1) / 2;
    int cy = cx;
    float Gxy;
    for (int y = -cy; y <= cy; y++)
    {
        for (int x = -cx; x <= cx; x++)
        {
            Gxy = 1.0 / (2.0 * M_PI * sigma * sigma) * exp(-((float)x * x + (float)y * y) / (2 * sigma * sigma));
            set_pixel(gf, cx + x, cy + y, 0, Gxy);
        }
    }
    l1_normalize(gf);
    return gf;
}

image add_image(image a, image b)
{
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image apb = make_image(a.w, a.h, a.c);
    for (int c = 0; c < apb.c; c++)
    {
        for (int y = 0; y < apb.h; y++)
        {
            for (int x = 0; x < apb.w; x++)
            {
                float axy = get_pixel(a, x, y, c);
                float bxy = get_pixel(b, x, y, c);
                set_pixel(apb, x, y, c, axy + bxy);
            }
        }
    }
    return apb;
}

image sub_image(image a, image b)
{
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image amb = make_image(a.w, a.h, a.c);
    for (int c = 0; c < amb.c; c++)
    {
        for (int y = 0; y < amb.h; y++)
        {
            for (int x = 0; x < amb.w; x++)
            {
                float axy = get_pixel(a, x, y, c);
                float bxy = get_pixel(b, x, y, c);
                set_pixel(amb, x, y, c, axy - bxy);
            }
        }
    }
    return amb;
}

image make_gx_filter()
{
    image gx = make_image(3, 3, 1);
    set_pixel(gx, 0, 0, 0, -1.0);
    set_pixel(gx, 1, 0, 0, 0.0);
    set_pixel(gx, 2, 0, 0, 1.0);
    set_pixel(gx, 0, 1, 0, -2.0);
    set_pixel(gx, 1, 1, 0, 0.0);
    set_pixel(gx, 2, 1, 0, 2.0);
    set_pixel(gx, 0, 2, 0, -1.0);
    set_pixel(gx, 1, 2, 0, 0.0);
    set_pixel(gx, 2, 2, 0, 1.0);
    return gx;
}

image make_gy_filter()
{
    image gy = make_image(3, 3, 1);
    set_pixel(gy, 0, 0, 0, -1.0);
    set_pixel(gy, 1, 0, 0, -2.0);
    set_pixel(gy, 2, 0, 0, -1.0);
    set_pixel(gy, 0, 1, 0, 0.0);
    set_pixel(gy, 1, 1, 0, 0.0);
    set_pixel(gy, 2, 1, 0, 0.0);
    set_pixel(gy, 0, 2, 0, 1.0);
    set_pixel(gy, 1, 2, 0, 2.0);
    set_pixel(gy, 2, 2, 0, 1.0);
    return gy;
}

void feature_normalize(image im)
{
    float max, min, v;
    for (int c = 0; c < im.c; c++)
    {
        for (int y = 0; y < im.h; y++)
        {
            for (int x = 0; x < im.w; x++)
            {
                v = get_pixel(im, x, y, c);
                max = v > max ? v : max;
                min = v < min ? v : min;
            }
        }
    }
    float s = max - min == 0 ? 1 : 1.0 / (max - min);
    for (int c = 0; c < im.c; c++)
    {
        shift_image(im, c, -min);
        scale_image(im, c, s);
    }
}

image *sobel_image(image im)
{
    image im_g = rgb_to_grayscale(im);
    image gy = make_gy_filter();
    image gx = make_gx_filter();
    image im_gy = convolve_image(im, gy, 0);
    image im_gx = convolve_image(im, gx, 0);
    image G = make_image(im.w, im.h, 1);
    image T = make_image(im.w, im.h, 1);
    float Gy, Gx, Gxy, Txy;
    for (int y = 0; y < G.w; y++)
    {
        for (int x = 0; x < G.w; x++)
        {
            Gy = get_pixel(im_gy, x, y, 0);
            Gx = get_pixel(im_gx, x, y, 0);
            Gxy = sqrt(pow(Gx, 2) + pow(Gy, 2));
            Txy = atan2(Gy, Gx);
            set_pixel(G, x, y, 0, Gxy);
            set_pixel(T, x, y, 0, Txy);
        }
    }

    image *sobel_images = calloc(2, sizeof(im_g));
    sobel_images[0] = G;
    sobel_images[1] = T;
    return sobel_images;
}

image colorize_sobel(image im)
{
    image *res = sobel_image(im);
    image mag = res[0];
    feature_normalize(mag);
    smooth_image(mag, 6);
    image theta = res[1];
    image im_color = make_image(im.w, im.h, im.c);
    float h, s, v;
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            h = get_pixel(theta, x, y, 0);
            s = get_pixel(mag, x, y, 0);
            v = s;
            set_pixel(im_color, x, y, 0, h);
            set_pixel(im_color, x, y, 1, s);
            set_pixel(im_color, x, y, 2, v);
        }
    }
    hsv_to_rgb(im_color);
    return im_color;
}
