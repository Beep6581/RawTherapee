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
//#define BENCHMARK
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

    const int bord = 2;

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int i = bord; i < winh - bord; ++i) {
        for(int j = bord; j< winw - bord; ++j) {
            int c = FC(i,j);
            if(c == 0) {
                red[i][j] = riFrames[0]->data[i][j];
                green[i][j] = riFrames[3]->data[i][j+1];
                blue[i][j] = riFrames[2]->data[i+1][j+1];
            } else if(c & 1) {
                green[i][j] = riFrames[0]->data[i][j];
                if(FC(i,j+1) == 0) {
                    red[i][j] = riFrames[3]->data[i][j+1];
                    blue[i][j] = riFrames[1]->data[i+1][j];
                } else {
                    blue[i][j] = riFrames[3]->data[i][j+1];
                    red[i][j] = riFrames[1]->data[i+1][j];
                }
            } else {
                blue[i][j] = riFrames[0]->data[i][j];
                red[i][j] = riFrames[2]->data[i+1][j+1];
                green[i][j] = riFrames[3]->data[i][j+1];
            }
        }
    }

    if(plistenerActive) {
        plistener->setProgress(1.00);
    }

}
#undef TS
#undef CLF
