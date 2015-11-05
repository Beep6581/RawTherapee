//
// A Halide implementation of "A Low Complexity Color Demosaicing Algorithm Based on Integrated Gradient" by
// King-Hong Chung and Yuk-Hee Chan:
//
// http://www.eie.polyu.edu.hk/~enyhchan/J-JEI-Low_complexity_color_demosaicing_algorithm_based_on_IG.pdf
//

#include <cmath>
#include "rawimagesource.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"
#include "opthelper.h"
#include "halide_debayer.h"

using namespace std;
using namespace rtengine;

void RawImageSource::khc_yhc_demosaic(int winx, int winy, int winw, int winh){
    double progress = 0;

    const bool plistenerActive = (const bool) plistener;

    if (plistener) {
        plistener->setProgressStr(
                Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"),
                                       RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::khc_yhc]));

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

    halide_debayer(&input_buf, layout, &output_buf);

    uint64_t plane_size = output_buf.extent[0] * output_buf.extent[1];
    float *output_start = (float *) output_buf.host;

    float *redaddr[H];
    float *greenaddr[H];
    float *blueaddr[H];

    float *redstart = output_start + (0 * output_buf.stride[2]);
    float *greenstart = output_start + (1 * output_buf.stride[2]);
    float *bluestart = output_start + (2 * output_buf.stride[2]);

    int rowstride = output_buf.stride[1];

    for(int i=0;i<H;i++) {

        redaddr[i] = redstart + rowstride * i;
        greenaddr[i] = greenstart + rowstride * i;
        blueaddr[i] = bluestart + rowstride * i;
    }

    red = *(new array2D<float>(W, H, redaddr, 0));
    green = *(new array2D<float>(W, H, greenaddr, 0));
    blue = *(new array2D<float>(W, H, blueaddr, 0));

    if(plistenerActive) {
        plistener->setProgress(1.00);
    }
}