#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    float sum = 0.0;
    for (int i = 0; i < im.c; i++) {
        for (int j = 0; j < im.w; j++) {
            for (int k = 0; k < im.w; k++) {
                sum += get_pixel(im, j, k, i);
            }
        }
    }
    for (int i = 0; i < im.c; i++) {
        for (int j = 0; j < im.w; j++) {
            for (int k = 0; k < im.w; k++) {
                set_pixel(im, j, k, i, get_pixel(im, j, k, i) / sum);
            }
        }
    }
}

image make_box_filter(int w)
{
    image filter = make_image(w, w, 1);
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            set_pixel(filter, i, j, 0, 1.0);
        }
    }
    l1_normalize(filter);
    return filter;
}

float convolve_pixel(image im, image filter, int x, int y, int c) {
    float weighted_sum = 0;
    for (int i = 0; i < filter.w; i++) {
        for (int j = 0; j < filter.w; j++) {
            float filter_pixel = get_pixel(filter, i, j, c);
            float im_pixel = get_pixel(im, x - filter.w / 2 + i, y - filter.w / 2 + j, c);
            weighted_sum += filter_pixel * im_pixel;
        }
    }
    return weighted_sum;
}

float convolve_pixel_to_monochannel(image im, image filter, int x, int y) {
    float sum = 0;
    for (int i = 0; i < im.c; i++) {
        sum += convolve_pixel(im, filter, x, y, i);
    }
    return sum;
}

image convolve_image(image im, image filter, int preserve)
{
    image convolved;
    if (im.c == filter.c) {
        if (preserve == 1) {
                        convolved = make_image(im.w, im.h, im.c);       
            for (int im_x = 0; im_x < im.w; im_x++) {
                for (int im_y = 0; im_y < im.h; im_y++) {
                    for (int im_c = 0; im_c < im.c; im_c++) {
                        float weighted_sum = 0.0;
                        for (int f_x = 0; f_x < filter.w; f_x++) {
                            for (int f_y = 0; f_y < filter.w; f_y++) {
                                float filter_pixel = get_pixel(filter, f_x, f_y, im_c);
                                float im_pixel = get_pixel(im, im_x - filter.w / 2 + f_x, im_y - filter.w / 2 + f_y, im_c);
                                weighted_sum += filter_pixel * im_pixel;
                            }
                        }
                        set_pixel(convolved, im_x, im_y, im_c, weighted_sum);
                    }
                }
            }
        } else {
                        convolved = make_image(im.w, im.h, 1);
            for (int im_x = 0; im_x < im.w; im_x++) {
                for (int im_y = 0; im_y < im.h; im_y++) {
                    float weighted_sum = 0.0;
                    for (int im_c = 0; im_c < im.c; im_c++) {
                        for (int f_x = 0; f_x < filter.w; f_x++) {
                            for (int f_y = 0; f_y < filter.w; f_y++) {
                                float filter_pixel = get_pixel(filter, f_x, f_y, im_c);
                                float im_pixel = get_pixel(im, im_x - filter.w / 2 + f_x, im_y - filter.w / 2 + f_y, im_c);
                                weighted_sum += filter_pixel * im_pixel;
                            }
                        }
                    }
                    set_pixel(convolved, im_x, im_y, 0, weighted_sum);
                }
            }
        }
    } else {
        if (preserve == 1) {
            convolved = make_image(im.w, im.h, im.c);
            for (int im_x = 0; im_x < im.w; im_x++) {
                for (int im_y = 0; im_y < im.h; im_y++) {
                    for (int im_c = 0; im_c < im.c; im_c++) {
                        float weighted_sum = 0.0;
                        for (int f_x = 0; f_x < filter.w; f_x++) {
                            for (int f_y = 0; f_y < filter.w; f_y++) {
                                float filter_pixel = get_pixel(filter, f_x, f_y, 0);
                                float im_pixel = get_pixel(im, im_x - filter.w / 2 + f_x, im_y - filter.w / 2 + f_y, im_c);
                                weighted_sum += filter_pixel * im_pixel;
                            }
                        }
                        set_pixel(convolved, im_x, im_y, im_c, weighted_sum);
                    }
                }
            }
        } else {
            convolved = make_image(im.w, im.h, 1);
            for (int im_x = 0; im_x < im.w; im_x++) {
                for (int im_y = 0; im_y < im.h; im_y++) {
                    float weighted_sum = 0.0;
                    for (int im_c = 0; im_c < im.c; im_c++) {
                        for (int f_x = 0; f_x < filter.w; f_x++) {
                            for (int f_y = 0; f_y < filter.w; f_y++) {
                                float filter_pixel = get_pixel(filter, f_x, f_y, 0);
                                float im_pixel = get_pixel(im, im_x - filter.w / 2 + f_x, im_y - filter.w / 2 + f_y, im_c);
                                weighted_sum += filter_pixel * im_pixel;
                            }
                        }
                    }
                    set_pixel(convolved, im_x, im_y, 0, weighted_sum);
                }
            }
        }
    }
    return convolved;
}

image make_highpass_filter()
{
    image f = make_image(3,3,1);
    set_pixel(f, 1, 0, 0, -1.0);
    set_pixel(f, 0, 1, 0, -1.0);
    set_pixel(f, 2, 1, 0, -1.0);
    set_pixel(f, 1, 2, 0, -1.0);
    set_pixel(f, 1, 1, 0, 4.0);
    return f;
}

image make_sharpen_filter()
{
    image f = make_image(3,3,1);
    set_pixel(f, 1, 0, 0, -1.0);
    set_pixel(f, 0, 1, 0, -1.0);
    set_pixel(f, 2, 1, 0, -1.0);
    set_pixel(f, 1, 2, 0, -1.0);
    set_pixel(f, 1, 1, 0, 5.0);
    return f;
}

image make_emboss_filter()
{
    image f = make_image(3,3,1);
    set_pixel(f, 0, 0, 0, -2.0);
    set_pixel(f, 1, 0, 0, -1.0);
    set_pixel(f, 0, 1, 0, -1.0);
    set_pixel(f, 1, 1, 0, 1.0);
    set_pixel(f, 2, 1, 0, 1.0);
    set_pixel(f, 1, 2, 0, 1.0);
    set_pixel(f, 2, 2, 0, 2.0);
    return f;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    int filter_size = (int) ceil(6 * sigma);
    filter_size += filter_size % 2 == 0 ? 1 : 0;
    image filter = make_image(filter_size, filter_size, 1);
    float shift = filter_size / 2.0 - 0.5;
    for (int i = 0; i < filter_size; i++) {
        for (int j = 0; j < filter_size; j++) {
            float x = i - shift;
            float y = j - shift;
            float v = 1.0 / (TWOPI * sigma * sigma) * expf(-(x * x + y * y) / (2 * sigma * sigma));
            set_pixel(filter, i, j, 0, v);
        }
    }
    l1_normalize(filter);
    return filter;
}

image add_image(image a, image b)
{
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image sum = make_image(a.w, b.h, a.c);
    for (int x = 0; x < a.w; x++) {
        for (int y = 0; y < a.h; y++) {
            for (int c = 0; c < a.c; c++) {
                float a_v = get_pixel(a, x, y, c);
                float b_v = get_pixel(b, x, y, c);
                set_pixel(sum, x, y, c, a_v + b_v);
            }
        }
    }
    return sum;
}

image sub_image(image a, image b)
{
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image difference = make_image(a.w, b.h, a.c);
    for (int x = 0; x < a.w; x++) {
        for (int y = 0; y < a.h; y++) {
            for (int c = 0; c < a.c; c++) {
                float a_v = get_pixel(a, x, y, c);
                float b_v = get_pixel(b, x, y, c);
                set_pixel(difference, x, y, c, a_v - b_v);
            }
        }
    }
    return difference;
}

image make_gx_filter()
{
    image f = make_image(3,3,1);
    set_pixel(f, 0, 0, 0, -1.0);
    set_pixel(f, 2, 0, 0, 1.0);
    set_pixel(f, 0, 1, 0, -2.0);
    set_pixel(f, 2, 1, 0, 2.0);
    set_pixel(f, 0, 2, 0, -1.0);
    set_pixel(f, 2, 2, 0, 1.0);
    return f;
}

image make_gy_filter()
{
    image f = make_image(3,3,1);
    set_pixel(f, 0, 0, 0, -1.0);
    set_pixel(f, 0, 2, 0, 1.0);
    set_pixel(f, 1, 0, 0, -2.0);
    set_pixel(f, 1, 2, 0, 2.0);
    set_pixel(f, 2, 0, 0, -1.0);
    set_pixel(f, 2, 2, 0, 1.0);
    return f;
}

void feature_normalize(image im)
{
    float min = get_pixel(im, 0,0,0);
    float max = min;
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            for (int c = 0; c < im.c; c++) {
                float v = get_pixel(im, x, y, c);
                min = v < min ? v : min;
                max = v > max ? v : max;
            }
        }
    }
    float range = max - min;
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            for (int c = 0; c < im.c; c++) {
                float v = range <= 0.0 ? 0.0 : (get_pixel(im, x, y, c) - min) / range;
                set_pixel(im, x, y, c, v);
            }
        }
    }
}

image *sobel_image(image im)
{
    image gx = make_gx_filter();
    image gy = make_gy_filter();
    image Gx = convolve_image(im, gx, 0);
    image Gy = convolve_image(im, gy, 0);

    image *mag_dir = calloc(2, sizeof(image));
    mag_dir[0] = make_image(im.w, im.h, 1);
    mag_dir[1] = make_image(im.w, im.h, 1);
    image mag = mag_dir[0];
    image dir = mag_dir[1];
    for (int x = 0; x < mag.w; x++) {
        for (int y = 0; y < mag.h; y++) {
            float Gx_v = get_pixel(Gx, x, y, 0);
            float Gy_v = get_pixel(Gy, x, y, 0);

            float mag_v = sqrtf(Gx_v * Gx_v + Gy_v * Gy_v);
            set_pixel(mag, x, y, 0, mag_v);

            // float dir_v = atanf(Gx_v == 0 ? 0.0 : (Gy_v / Gx_v));
            float dir_v = atan2f(Gy_v, Gx_v);
            set_pixel(dir, x, y, 0, dir_v);
        }
    }
    return mag_dir;
}

image colorize_sobel(image im)
{
    image *mag_dir = sobel_image(im);
    image mag = mag_dir[0];
    image dir = mag_dir[1];
    feature_normalize(mag);
    feature_normalize(dir);
    image result = make_image(im.w, im.h, 3);
    for (int x = 0; x < im.w; x++) {
        for (int y = 0; y < im.h; y++) {
            float dir_v = get_pixel(dir, x, y, 0);
            set_pixel(result, x, y, 0, dir_v);

            float mag_v = get_pixel(mag, x, y, 0);
            set_pixel(result, x, y, 1, mag_v);
            set_pixel(result, x, y, 2, mag_v);
        }
    }
    free_image(mag);
    free_image(dir);
    hsv_to_rgb(result);
    return result;
}
