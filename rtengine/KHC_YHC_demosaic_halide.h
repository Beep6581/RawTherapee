//
// Created by Lucas Cooper on 3/11/2015.
//

#ifndef RAWTHERAPEE_KHC_YHC_DEMOSAIC_HALIDE_H_H
#define RAWTHERAPEE_KHC_YHC_DEMOSAIC_HALIDE_H_H

#endif //RAWTHERAPEE_KHC_YHC_DEMOSAIC_HALIDE_H_H

#include <Halide.h>
#include <HalideRuntime.h>

Halide::Func make_demosaic_func(Halide::ImageParam param, Halide::Type out_type);