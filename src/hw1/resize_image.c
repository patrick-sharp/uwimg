#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    return get_pixel(im, (int) round(x), (int) round(y), c);
}

image nn_resize(image im, int w, int h)
{
    image resized = make_image(w, h, im.c);
    for (int c = 0; c < im.c; c++) {
        for (int x = 0; x < w; x++) {
            for (int y = 0; y < h; y++) {
                float new_x = (float) im.w / w * (x + 0.5) - 0.5;
                float new_y = (float) im.h / h * (y + 0.5) - 0.5;
                set_pixel(resized, x, y, c, nn_interpolate(im, new_x, new_y, c));
            }
        }
    }
    return resized;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    float flx = floor(x);
    float fly = floor(y);

    float d1 = x - flx;
    float d2 = ceil(x) - x;
    float d3 = y - fly;
    float d4 = ceil(y) - y;

    if (d2 == 0.0) d2 = 1.0;
    if (d4 == 0.0) d4 = 1.0;

    int i_flx = (int) flx;
    int i_fly = (int) fly;

    float v1 = get_pixel(im, i_flx, i_fly, c);
    float v2 = get_pixel(im, i_flx + 1, i_fly, c);
    float v3 = get_pixel(im, i_flx, i_fly + 1, c);
    float v4 = get_pixel(im, i_flx + 1, i_fly + 1, c);

    float q1 = v1 * d2 + v2 * d1;
    float q2 = v3 * d2 + v4 * d1;
    return q1 * d4 + q2 * d3;
}

image bilinear_resize(image im, int w, int h)
{
    image resized = make_image(w, h, im.c);
    for (int c = 0; c < im.c; c++) {
        for (int x = 0; x < w; x++) {
            for (int y = 0; y < h; y++) {
                float new_x = (float) im.w / w * (x + 0.5) - 0.5;
                float new_y = (float) im.h / h * (y + 0.5) - 0.5;
                set_pixel(resized, x, y, c, bilinear_interpolate(im, new_x, new_y, c));
            }
        }
    }
    return resized;
}

