//
// A wrapper around the Halide implementation of "A Low Complexity Color Demosaicing Algorithm Based on Integrated
// Gradient" by King-Hong Chung and Yuk-Hee Chan:
//
// http://www.eie.polyu.edu.hk/~enyhchan/J-JEI-Low_complexity_color_demosaicing_algorithm_based_on_IG.pdf
//

#include <cmath>
#include "rawimagesource.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"
#include "opthelper.h"
#include "igd_demosaic.h"
#include "StopWatch.h"
#include <ctime>
#include <stdio.h>
#include <iostream>


#ifdef HALIDE_OPENCL
#include "halide/HalideRuntime.h"
#include "halide/HalideRuntimeOpenCL.h"
#endif

using namespace std;
using namespace rtengine;

void RawImageSource::igd_halide_demosaic(int winx, int winy, int winw, int winh){
    StopWatch measure("IGD");
    double progress = 0;

    const bool plistenerActive = (const bool) plistener;

    if (plistener) {
        plistener->setProgressStr(
                Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"),
                                       RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::igd]));

        plistener->setProgress(progress);
    }

    buffer_t input_buf = {0};
    input_buf.host = (uint8_t *)(const float *)rawData;

    input_buf.stride[0] = 1;
    input_buf.stride[1] = W;

    input_buf.extent[0] = W;
    input_buf.extent[1] = H;

    input_buf.elem_size = sizeof(float);

    buffer_t output_buf = {0};

    output_buf.stride[0] = 1;
    output_buf.stride[1] = W;
    output_buf.stride[2] = W * H;
    output_buf.stride[3] = W * H * 3;

    output_buf.extent[0] = W;
    output_buf.extent[1] = H;
    output_buf.extent[2] = 3;
    output_buf.elem_size = sizeof(float);

    uint8_t *output_data = new uint8_t[sizeof(float) * W * H * 3];
    output_buf.host = output_data;
    output_buf.host_dirty = false;
    output_buf.dev_dirty = false;
    output_buf.dev = 0;

    uint8_t layout = 0;
    if (FC(0, 0) == 0) {
        layout = 0;
    } else if (FC(0, 1) == 0) {
        layout = 1;
    } else if (FC(1, 0) == 0) {
        layout = 2;
    } else {
        layout = 3;
    }

#ifdef HALIDE_OPENCL
    input_buf.host_dirty = true;
    halide_copy_to_device(NULL, &input_buf, halide_opencl_device_interface());

    buffer_t input_no_host = *(&input_buf);
    input_no_host.host = NULL;

    buffer_t output_no_host = *(&output_buf);
    output_no_host.host = (uint8_t *) 1;

    igd_demosaic(&input_no_host, layout, &output_no_host);

    output_no_host.host = output_data;
    halide_copy_to_host(NULL, &output_no_host);
#else
    igd_demosaic(&input_buf, layout, &output_buf);
#endif
    float *output_start = (float *) output_buf.host;

    float *redstart = output_start + (0 * output_buf.stride[2]);
    float *greenstart = output_start + (1 * output_buf.stride[2]);
    float *bluestart = output_start + (2 * output_buf.stride[2]);

    int rowstride = output_buf.stride[1];

    for(int i=0;i<H;i++) {
        ((float **)red)[i] = redstart + rowstride * i;
        ((float **)green)[i] = greenstart + rowstride * i;
        ((float **)blue)[i] = bluestart + rowstride * i;
    }

    if(plistenerActive) {
        plistener->setProgress(1.00);
    }
}