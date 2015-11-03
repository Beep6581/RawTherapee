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
#include "KHC_YHC_demosaic_halide.h"

using namespace std;
using namespace rtengine;

void RawImageSource::khc_yhc_demosaic(int winx, int winy, int winw, int winh){
    double progress = 0;

    const bool plistenerActive = (const bool) plistener;

    if (plistener) {
        plistener->setProgressStr (
            Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"),
            RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::khc_yhc]));

        plistener->setProgress (progress);
    }

    buffer_t input_buf = {0};
    input_buf.host = (uint8_t *) &rawData[0][0];
    input_buf.stride[0] = 1;
    input_buf.stride[1] = W;
    input_buf.extent[0] = W;
    input_buf.extent[1] = H;
    input_buf.elem_size = sizeof(rawData[0][0]);

    Halide::Image<float> input_image = Halide::Image<float>(&input_buf, "image");
    Halide::ImageParam param = Halide::ImageParam(Halide::Float(32), 2);
    Halide::Func output = make_demosaic_func(param);
    param.set(input_image);
    Halide::Image<float> output_image = output.realize(input_image.width(), input_image.height(), 3);
    output_image.copy_to_host();

    float *redaddr = (float *)output_image.address_of(0, 0, 0);
    float *greenaddr = (float *)output_image.address_of(0, 0, 1);
    float *blueaddr = (float *)output_image.address_of(0, 0, 2);

    red = *(new array2D<float>(W, H, (float **) &redaddr, 0));
    green = *(new array2D<float>(W, H, (float **) &greenaddr, 0));
    blue = *(new array2D<float>(W, H, (float **) &blueaddr, 0));

    if(plistenerActive) {
        plistener->setProgress(1.00);
    }
}