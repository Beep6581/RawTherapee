#ifndef CALC_DISTORTION__H
#define CALC_DISTORTION__H
extern "C" {
    double calcDistortion (unsigned char* img1, unsigned char* img2, int ncols, int nrows);
}
#endif
