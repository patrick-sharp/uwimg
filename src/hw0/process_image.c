#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

void clamp_pad(image im, int *x, int *y, int *c) {
    *c = MAX(*c, 0);
    *c = MIN(*c, im.c - 1);
    *x = MAX(*x, 0);
    *x = MIN(*x, im.w - 1);
    *y = MAX(*y, 0);
    *y = MIN(*y, im.h - 1);
}

int get_index(image im, int x, int y, int c) {
    return im.w * im.h * c + y * im.w + x;
}

float get_pixel(image im, int x, int y, int c)
{
    clamp_pad(im, &x, &y, &c);
    return im.data[get_index(im, x, y, c)];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if (x >= im.w || x < 0 || y >= im.h || y < 0 || c >= im.c || c < 0) {
        return;
    }
    im.data[get_index(im, x, y, c)] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, copy.w * copy.h * copy.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    for (int i = 0; i < im.w * im.h; i++) {
        gray.data[i] = 0.299 * im.data[i] + 0.587 * im.data[im.w * im.h + i] + .114 * im.data[im.w * im.h * 2 + i];
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    for (int i = 0; i < im.h * im.w; i++) {
        im.data[c * im.h * im.w + i] += v;
    }
}

void clamp_image(image im)
{
    for (int i = 0; i < im.c * im.h * im.w; i++) {
        im.data[i] = MAX(0.0, im.data[i]);
        im.data[i] = MIN(1.0, im.data[i]);
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    for (int i = 0; i < im.h * im.w; i++) {
        float R = im.data[i];
        float G = im.data[im.h * im.w + i];
        float B = im.data[2 * im.h * im.w + i];
        float V = three_way_max(R, G, B);
        float m = three_way_min(R, G, B);
        float C = V - m;
        float S = V != 0 ? C / V : 0;
        float Hprime;
        float H;
        if (C == 0) {
            H = 0;
        } else {
            if (V == R) {
                Hprime = (G - B) / C;
                Hprime = Hprime < 0 ? Hprime + 6 : Hprime;
            } else if (V == G) {
                Hprime = (B - R) / C + 2;
            } else {
                Hprime = (R - G) / C + 4;
            }
            H = Hprime / 6;
        }
        im.data[i] = H;
        im.data[im.h * im.w + i] = S;
        im.data[2 * im.h * im.w + i] = V;
    }
}

void hsv_to_rgb(image im)
{
    for (int i = 0; i < im.h * im.w; i++) {
        float H = im.data[i];
        float S = im.data[im.h * im.w + i];
        float V = im.data[2 * im.h * im.w + i];
        float C;
        if (V == 0.0) {
            C = 0.0;
        } else {
            C = S * V;
        }
        float m = V - C;
        float R;
        float G;
        float B;
        float Hprime = H * 6;
        if (0 <= Hprime && Hprime < 1) {
            R = V;
            B = m;
            G = B + Hprime * C;
        } else if (1 <= Hprime && Hprime < 2) {
            G = V;
            B = m;
            R = B - (Hprime - 2) * C;
        } else if (2 <= Hprime && Hprime < 3) {
            R = m;
            G = V;
            B = R + (Hprime - 2) * C;
        } else if (3 <= Hprime && Hprime < 4) {
            R = m;
            B = V;
            G = R - (Hprime - 4) * C;
        } else if (4 <= Hprime && Hprime < 5) {
            G = m;
            B = V;
            R = G + (Hprime - 4) * C;
        } else {
            R = V;
            G = m;
            B = G - (Hprime - 6) * C;
        }
        im.data[i] = R;
        im.data[im.h * im.w + i] = G;
        im.data[2 * im.h * im.w + i] = B;
    }
}
