# Halide

The files in this directory support the use of Halide based image processing pipelines.

The files in this root directory are copied from the [Halide project](https://github.com/halide/Halide) repository and
are under the license in HALIDE_LICENSE.

The files in the generators directory are part of RawTherapee and are licensed under the LICENSE.txt in the root of this
repository.

## Using

To build RawTherapee with Halide support, the following steps are required:

- Acquire and build [Halide](https://github.com/halide/Halide). Make note of the build output directory.
- Configure RawTherapee with CMake, setting `HALIDE_PATH` to the Halide build directory (the directory that contains
  `lib/libHalide.<platform specific extension>` and `include/Halide.h`) and, if GPU acceleration is desired, set
  `OPTION_HALIDE_OPENCL` to `ON` (the default, as it supports CPU as well).
- Build.

## Generators

The current Halide accelerated processes are:

- [IGD demosaic](http://www.eie.polyu.edu.hk/~enyhchan/J-JEI-Low_complexity_color_demosaicing_algorithm_based_on_IG.pdf) (used with permission of the authors)