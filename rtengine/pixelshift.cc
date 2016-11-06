////////////////////////////////////////////////////////////////
//
//  simple pentax pixelshift algorithm
//  copyright (c) Ingo Weyrich 2016
//
//
//  pixelshift.cc is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include <cmath>
#include "rawimagesource.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"
#include "opthelper.h"
#define BENCHMARK
#include "StopWatch.h"
using namespace std;
using namespace rtengine;

void RawImageSource::pixelshift_simple(int winx, int winy, int winw, int winh)
{

BENCHFUN
    double progress = 0.0;
    const bool plistenerActive = plistener;

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::pixelshift_simple]));
        plistener->setProgress (progress);
    }

    const int bord = 4;

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int i = bord; i < winh - bord; ++i) {
        float *greenDest = green[i];
        float *nonGreenDest0 = red[i];
        float *nonGreenDest1 = blue[i];
        int j = bord;
        int c = FC(i,j);
        if (c == 2 || ((c&1) && FC(i,j+1) == 2)) {
            std::swap(nonGreenDest0, nonGreenDest1);
        }
        if(c&1) {
            greenDest[j] = (riFrames[0]->data[i][j] + riFrames[2]->data[i+1][j+1]) / 2.f;
            nonGreenDest0[j] = riFrames[3]->data[i][j+1];
            nonGreenDest1[j] = riFrames[1]->data[i+1][j];
            j++;
        }
        for(; j< winw - bord; j+=2) {
            nonGreenDest0[j] = riFrames[0]->data[i][j];
            greenDest[j] = (riFrames[3]->data[i][j+1] + riFrames[1]->data[i+1][j] ) / 2.f;
            nonGreenDest1[j] = riFrames[2]->data[i+1][j+1];
            greenDest[j+1] = (riFrames[0]->data[i][j+1] + riFrames[2]->data[i+1][j+2]) / 2.f;
            nonGreenDest0[j+1] = riFrames[3]->data[i][j+2];
            nonGreenDest1[j+1] = riFrames[1]->data[i+1][j+1];
        }
    }

    if(plistenerActive) {
        plistener->setProgress(1.00);
    }

}
#undef TS
#undef CLF
