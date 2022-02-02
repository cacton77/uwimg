#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // Perform clamping if necessary
    x = x < 0 ? 0 : (x >= im.w ? im.w - 1 : x);
    y = y < 0 ? 0 : (y >= im.h ? im.h - 1 : y);
    c = c < 0 ? 0 : (c >= im.c ? im.c - 1 : c);

    float v = im.data[c * (im.w * im.h) + y * im.w + x];
    return v;
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // Perform clamping if necessary
    x = x < 0 ? 0 : (x >= im.w ? im.w - 1 : x);
    y = y < 0 ? 0 : (y >= im.h ? im.h - 1 : y);
    c = c < 0 ? 0 : (c >= im.c ? im.c - 1 : c);

    im.data[c * (im.w * im.h) + y * im.w + x] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, (im.w * im.h * im.c) * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            float r = im.data[y * im.w + x];
            float g = im.data[im.w * im.h + y * im.w + x];
            float b = im.data[2 * im.w * im.h + y * im.w + x];
            float v = 0.299 * r + 0.587 * g + 0.114 * b;
            gray.data[y * im.w + x] = v;
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            im.data[c * im.w * im.h + y * im.w + x] += v;
        }
    }
}

void scale_image(image im, int c, float v)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            im.data[c * im.w * im.h + y * im.w + x] *= v;
        }
    }
}

void clamp_image(image im)
{
    for (int c = 0; c < im.c; c++)
    {
        for (int y = 0; y < im.h; y++)
        {
            for (int x = 0; x < im.w; x++)
            {
                float v = im.data[c * im.w * im.h + y * im.w + x];
                v = v > 1.0 ? 1.0 : (v < 0.0 ? 0.0 : v);
                im.data[c * im.w * im.h + y * im.w + x] = v;
            }
        }
    }
}

// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

void rgb_to_hsv(image im)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            float r = im.data[y * im.w + x];
            float g = im.data[im.w * im.h + y * im.w + x];
            float b = im.data[2 * im.w * im.h + y * im.w + x];
            float c, h, s, v;
            // Value
            v = three_way_max(r, g, b);
            // Saturation
            c = v - three_way_min(r, g, b);
            s = v == 0.0 ? 0.0 : c / v;
            // Hue
            if (c == 0)
            {
                h = 0;
            }
            else if (v == r)
            {
                h = (g - b) / c;
            }
            else if (v == g)
            {
                h = (b - r) / c + 2;
            }
            else // if v == b
            {
                h = (r - g) / c + 4;
            }
            if (h < 0)
            {
                h = h / 6 + 1;
            }
            else
            {
                h = h / 6;
            }
            im.data[y * im.w + x] = h;
            im.data[im.w * im.h + y * im.w + x] = s;
            im.data[2 * im.w * im.h + y * im.w + x] = v;
        }
    }
}

void hsv_to_rgb(image im)
{
    for (int y = 0; y < im.h; y++)
    {
        for (int x = 0; x < im.w; x++)
        {
            float h = 6 * im.data[y * im.w + x];
            float s = im.data[im.w * im.h + y * im.w + x];
            float v = im.data[2 * im.w * im.h + y * im.w + x];
            float c = s * v;
            float p = c * (1 - fabs(fmod(h, 2) - 1));
            float m = v - c;
            float r, g, b;
            if (h >= 0 && h < 1)
            {
                r = c, g = p, b = 0;
            }
            else if (h >= 1 && h < 2)
            {
                r = p, g = c, b = 0;
            }
            else if (h >= 2 && h < 3)
            {
                r = 0, g = c, b = p;
            }
            else if (h >= 3 && h < 4)
            {
                r = 0, g = p, b = c;
            }
            else if (h >= 4 && h < 5)
            {
                r = p, g = 0, b = c;
            }
            else
            {
                r = c, g = 0, b = p;
            }
            r += m, g += m, b += m;

            im.data[y * im.w + x] = r;
            im.data[im.w * im.h + y * im.w + x] = g;
            im.data[2 * im.w * im.h + y * im.w + x] = b;
        }
    }
    // TODO Fill this in
}
