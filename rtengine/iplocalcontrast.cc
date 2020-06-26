/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Ported from G'MIC by Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  The original implementation in G'MIC was authored by Arto Huotari, and was
 *  released under the CeCILL free software license (see
 *  http://www.cecill.info/licences/Licence_CeCILL_V2-en.html)
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "array2D.h"
#include "gauss.h"
#include "labimage.h"
#include "improcfun.h"
#include "procparams.h"
#include "settings.h"

namespace rtengine
{

void ImProcFunctions::localContrast(LabImage *lab, float **destination, const rtengine::procparams::LocalContrastParams &localContrastParams, bool fftwlc, double scale)
{
    if (!localContrastParams.enabled) {
        return;
    }

    const int width = lab->W;
    const int height = lab->H;
    const float a = localContrastParams.amount;
    const float dark = localContrastParams.darkness;
    const float light = localContrastParams.lightness;
    array2D<float> buf(width, height);
    float sigma = localContrastParams.radius / scale;
    //printf("wi%i he=%i am=%f da=%f li=%f si=%f\n", width, height, a, dark, light, sigma);
    if(!fftwlc) {
#ifdef _OPENMP
        #pragma omp parallel if(multiThread)
#endif
        gaussianBlur(lab->L, buf, width, height, sigma);
    } else {
        float kr = 1.f;
        //emprical adjustment between FFTW radius and Gaussainblur
        //under 50 ==> 10.f
        //above 400 ==> 1.f
        if(settings->fftwsigma == false) {//empirical formula
            float ak = -9.f / 350.f;
            float bk = 10.f - 50.f * ak;
            kr = ak * sigma + bk;
            if(sigma < 50.f) kr = 10.f;
            if(sigma > 400.f) kr = 1.f;
        } else {//sigma *= sigma
            kr = sigma;
        }
        //OPENMP disabled
        ImProcFunctions::fftw_convol_blur2(lab->L, buf, width, height, kr * sigma, 0, 0);
    }
#ifdef _OPENMP
    #pragma omp parallel for if(multiThread)
#endif

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float bufval = (lab->L[y][x] - buf[y][x]) * a;

            if (dark != 1 || light != 1) {
                bufval *= (bufval > 0.f) ? light : dark;
            }

            destination[y][x] = LIM(lab->L[y][x] + bufval, 0.0001f, 32767.f);
        }
    }
}

} // namespace rtengine
