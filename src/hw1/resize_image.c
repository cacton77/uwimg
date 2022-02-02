#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    int x_nn = (int)round(x);
    int y_nn = (int)round(y);
    float v = get_pixel(im, x_nn, y_nn, c);
    return v;
}

image nn_resize(image im, int w, int h)
{
    image im_rszd = make_image(w, h, im.c);
    float ax = (float)im.w / w;
    float ay = (float)im.h / h;
    float bx = -0.5 + 0.5 * ax;
    float by = -0.5 + 0.5 * ay;
    for (int c = 0; c < w; c++)
    {
        for (int y = 0; y < w; y++)
        {
            float y_mapd = ay * y + by;
            for (int x = 0; x < w; x++)
            {
                float x_mapd = ax * x + bx;
                float v = nn_interpolate(im, x_mapd, y_mapd, c);
                set_pixel(im_rszd, x, y, c, v);
            }
        }
    }
    return im_rszd;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    // V1_ _ _ _ _ _ _ _ _V2
    // |          |       |
    // |    A4    |   A3  |
    // |          |       |
    // |_ _ _ _ _ Q _ _ _ |
    // |    A2    |   A1  |
    // V3_ _ _ _ _|_ _ _ _V4

    float Q, V1, V2, V3, V4;
    int x_lwr = (int)floor(x);
    int x_upr = (int)ceil(x);
    int y_lwr = (int)floor(y);
    int y_upr = (int)ceil(y);

    V1 = get_pixel(im, x_lwr, y_lwr, c);
    V2 = get_pixel(im, x_upr, y_lwr, c);
    V3 = get_pixel(im, x_lwr, y_upr, c);
    V4 = get_pixel(im, x_upr, y_upr, c);

    float A1 = (ceil(x) - x) * (ceil(y) - y);
    float A2 = (x - floor(x)) * (ceil(y) - y);
    float A3 = (ceil(x) - x) * (y - floor(y));
    float A4 = (x - floor(x)) * (y - floor(y));

    Q = V1 * A1 + V2 * A2 + V3 * A3 + V4 * A4;

    return Q;
}

image bilinear_resize(image im, int w, int h)
{
    image im_rszd = make_image(w, h, im.c);
    float ax = (float)im.w / w;
    float ay = (float)im.h / h;
    float bx = -0.5 + 0.5 * ax;
    float by = -0.5 + 0.5 * ay;
    for (int c = 0; c < im.c; c++)
    {
        for (int y = 0; y < h; y++)
        {
            float y_mapd = ay * y + by;
            for (int x = 0; x < w; x++)
            {
                float x_mapd = ax * x + bx;
                float v = bilinear_interpolate(im, x_mapd, y_mapd, c);
                set_pixel(im_rszd, x, y, c, v);
            }
        }
    }
    return im_rszd;
}
