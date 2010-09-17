/*
 * improcfuns.cc
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#include "improcfuns.h"

namespace rtengine {

void ImProcFunctions::getAutoExp  (unsigned int* histogram, int histcompr, double expcomp, double clip, double& br, int& bl) {

    double sum = 0;
    for (int i=0; i<65536>>histcompr; i++)
        sum += histogram[i];

    // compute clipping points based on the original histograms (linear, without exp comp.)
    int clippable = (int)(sum * clip);
    int clipped = 0;
    int aw = (65536>>histcompr) - 1;
    while (aw>1 && histogram[aw]+clipped <= clippable) {
        clipped += histogram[aw];
        aw--;
    }

    clipped = 0;
    int shc = 0;
    while (shc<aw-1 && histogram[shc]+clipped <= clippable) {
        clipped += histogram[shc];
        shc++;
    }

    aw <<= histcompr;
    shc <<= histcompr;

    double corr = pow(2.0, expcomp);

    // black point selection is based on the linear result (yielding better visual results)
    bl = (int)(shc * corr);
    // compute the white point of the exp. compensated gamma corrected image
    double awg = (int)(CurveFactory::gamma2 (aw * corr / 65536.0) * 65536.0);

    // compute average intensity of the exp compensated, gamma corrected image
    double gavg = 0;
    for (int i=0; i<65536>>histcompr; i++)
        gavg += histogram[i] * CurveFactory::gamma2((int)(corr*(i<<histcompr)<65535 ? corr*(i<<histcompr) : 65535)) / sum;


    if (bl < gavg) {
        int maxaw = (gavg - bl) * 4 / 3 + bl; // dont let aw be such large that the histogram average goes above 3/4
        double mavg = 65536.0 / (awg-bl) * (gavg - bl);
        if (awg < maxaw)
            awg = maxaw;
    }

    br = log(65535.0 / (awg-bl)) / log(2.0);
    if (br<0)
        br = 0;
}
}
