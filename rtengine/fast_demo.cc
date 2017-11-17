////////////////////////////////////////////////////////////////
//
//      Fast demosaicing algorythm
//
//      copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: August 26, 2010
//
//  fast_demo.cc is free software: you can redistribute it and/or modify
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

using namespace std;
using namespace rtengine;

#define TS 224
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*
LUTf RawImageSource::initInvGrad()
{
    LUTf invGrad (0x10000);

    //set up directional weight function
    for (int i=0; i<0x10000; i++)
        invGrad[i] = 1.0/SQR(1.0+i);

    return invGrad;
}
*/
#define INVGRAD(i) (16.0f/SQR(4.0f+i))
#ifdef __SSE2__
#define INVGRADV(i) (c16v*_mm_rcp_ps(SQRV(fourv+i)))
#endif
//LUTf RawImageSource::invGrad = RawImageSource::initInvGrad();

SSEFUNCTION void RawImageSource::fast_demosaic(int winx, int winy, int winw, int winh)
{

    double progress = 0.0;
    const bool plistenerActive = plistener;

    //int winx=0, winy=0;
    //int winw=W, winh=H;

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST)));
        plistener->setProgress (progress);
    }


    const int bord = 5;

    float clip_pt = 4 * 65535 * initialGain;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

        char (*buffer);
        float (*greentile);
        float (*redtile);
        float (*bluetile);
#define CLF 1
        // assign working space
        buffer = (char *) calloc(3 * sizeof(float) * TS * TS + 3 * CLF * 64 + 63, 1);
        char    *data;
        data = (char*)( ( uintptr_t(buffer) + uintptr_t(63)) / 64 * 64);

        greentile = (float (*))         data; //pointers to array
        redtile   = (float (*))         ((char*)greentile + sizeof(float) * TS * TS + CLF * 64);
        bluetile  = (float (*))         ((char*)redtile + sizeof(float) * TS * TS + CLF * 64);

#ifdef _OPENMP
        #pragma omp sections
#endif
        {
#ifdef _OPENMP
            #pragma omp section
#endif
            {

                //first, interpolate borders using bilinear
                for (int i = 0; i < H; i++)
                {

                    float sum[6];
                    int imin = max(0, i - 1);
                    int imax = min(i + 2, H);

                    for (int j = 0; j < bord; j++) { //first few columns
                        for (int c = 0; c < 6; c++) {
                            sum[c] = 0;
                        }

                        int jmin = max(0, j - 1);

                        for (int i1 = imin; i1 < imax; i1++)
                            for (int j1 = jmin; j1 < j + 2; j1++) {
                                int c = FC(i1, j1);
                                sum[c] += rawData[i1][j1];
                                sum[c + 3]++;
                            }

                        int c = FC(i, j);

                        if (c == 1) {
                            red[i][j] = sum[0] / sum[3];
                            green[i][j] = rawData[i][j];
                            blue[i][j] = sum[2] / sum[5];
                        } else {
                            green[i][j] = sum[1] / sum[4];

                            if (c == 0) {
                                red[i][j] = rawData[i][j];
                                blue[i][j] = sum[2] / sum[5];
                            } else {
                                red[i][j] = sum[0] / sum[3];
                                blue[i][j] = rawData[i][j];
                            }
                        }
                    }//j

                    for (int j = W - bord; j < W; j++) { //last few columns
                        for (int c = 0; c < 6; c++) {
                            sum[c] = 0;
                        }

                        int jmax = min(j + 2, W);

                        for (int i1 = imin; i1 < imax; i1++)
                            for (int j1 = j - 1; j1 < jmax; j1++) {
                                int c = FC(i1, j1);
                                sum[c] += rawData[i1][j1];
                                sum[c + 3]++;
                            }

                        int c = FC(i, j);

                        if (c == 1) {
                            red[i][j] = sum[0] / sum[3];
                            green[i][j] = rawData[i][j];
                            blue[i][j] = sum[2] / sum[5];
                        } else {
                            green[i][j] = sum[1] / sum[4];

                            if (c == 0) {
                                red[i][j] = rawData[i][j];
                                blue[i][j] = sum[2] / sum[5];
                            } else {
                                red[i][j] = sum[0] / sum[3];
                                blue[i][j] = rawData[i][j];
                            }
                        }
                    }//j
                }//i

            }
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef _OPENMP
            #pragma omp section
#endif
            {

                for (int j = bord; j < W - bord; j++)
                {
                    float sum[6];

                    for (int i = 0; i < bord; i++) { //first few rows
                        for (int c = 0; c < 6; c++) {
                            sum[c] = 0;
                        }

                        for (int i1 = max(0, i - 1); i1 < i + 2; i1++)
                            for (int j1 = j - 1; j1 < j + 2; j1++) {
                                int c = FC(i1, j1);
                                sum[c] += rawData[i1][j1];
                                sum[c + 3]++;
                            }

                        int c = FC(i, j);

                        if (c == 1) {
                            red[i][j] = sum[0] / sum[3];
                            green[i][j] = rawData[i][j];
                            blue[i][j] = sum[2] / sum[5];
                        } else {
                            green[i][j] = sum[1] / sum[4];

                            if (c == 0) {
                                red[i][j] = rawData[i][j];
                                blue[i][j] = sum[2] / sum[5];
                            } else {
                                red[i][j] = sum[0] / sum[3];
                                blue[i][j] = rawData[i][j];
                            }
                        }
                    }//i

                    for (int i = H - bord; i < H; i++) { //last few rows
                        for (int c = 0; c < 6; c++) {
                            sum[c] = 0;
                        }

                        for (int i1 = i - 1; i1 < min(i + 2, H); i1++)
                            for (int j1 = j - 1; j1 < j + 2; j1++) {
                                int c = FC(i1, j1);
                                sum[c] += rawData[i1][j1];
                                sum[c + 3]++;
                            }

                        int c = FC(i, j);

                        if (c == 1) {
                            red[i][j] = sum[0] / sum[3];
                            green[i][j] = rawData[i][j];
                            blue[i][j] = sum[2] / sum[5];
                        } else {
                            green[i][j] = sum[1] / sum[4];

                            if (c == 0) {
                                red[i][j] = rawData[i][j];
                                blue[i][j] = sum[2] / sum[5];
                            } else {
                                red[i][j] = sum[0] / sum[3];
                                blue[i][j] = rawData[i][j];
                            }
                        }
                    }//i
                }//j

            }
        }
#ifdef _OPENMP
        #pragma omp single
#endif
        {
            if(plistenerActive) {
                progress += 0.1;
                plistener->setProgress(progress);
            }
        }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        int progressCounter = 0;
        const double progressInc = 16.0 * (1.0 - progress) / ((H * W) / ((TS - 4) * (TS - 4)));

#ifdef _OPENMP
        #pragma omp for nowait
#endif

        for (int top = bord - 2; top < H - bord + 2; top += TS - (4))
            for (int left = bord - 2; left < W - bord + 2; left += TS - (4)) {
                int bottom = min(top + TS, H - bord + 2);
                int right  = min(left + TS, W - bord + 2);

#ifdef __SSE2__
                int j, cc;
                __m128 wtuv, wtdv, wtlv, wtrv;
                __m128 greenv, tempv, absv, abs2v;
                __m128 c16v = _mm_set1_ps( 16.0f );
                __m128 fourv = _mm_set1_ps( 4.0f );
                vmask selmask;
                vmask andmask = _mm_set_epi32( 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff );

                if(FC(top, left) == 1) {
                    selmask = _mm_set_epi32( 0, 0xffffffff, 0, 0xffffffff );
                } else {
                    selmask = _mm_set_epi32( 0xffffffff, 0, 0xffffffff, 0 );
                }

#endif

                // interpolate G using gradient weights
                for (int i = top, rr = 0; i < bottom; i++, rr++) {
                    float   wtu, wtd, wtl, wtr;
#ifdef __SSE2__
                    selmask = (vmask)_mm_andnot_ps( (__m128)selmask, (__m128)andmask);

                    for (j = left, cc = 0; j < right - 3; j += 4, cc += 4) {
                        tempv = LVFU(rawData[i][j]);
                        absv = vabsf(LVFU(rawData[i - 1][j]) - LVFU(rawData[i + 1][j]));
                        wtuv = INVGRADV(absv + vabsf(tempv - LVFU(rawData[i - 2][j])) + vabsf(LVFU(rawData[i - 1][j]) - LVFU(rawData[i - 3][j])));
                        wtdv = INVGRADV(absv + vabsf(tempv - LVFU(rawData[i + 2][j])) + vabsf(LVFU(rawData[i + 1][j]) - LVFU(rawData[i + 3][j])));
                        abs2v = vabsf(LVFU(rawData[i][j - 1]) - LVFU(rawData[i][j + 1]));
                        wtlv = INVGRADV(abs2v + vabsf(tempv - LVFU(rawData[i][j - 2])) + vabsf(LVFU(rawData[i][j - 1]) - LVFU(rawData[i][j - 3])));
                        wtrv = INVGRADV(abs2v + vabsf(tempv - LVFU(rawData[i][j + 2])) + vabsf(LVFU(rawData[i][j + 1]) - LVFU(rawData[i][j + 3])));
                        greenv = (wtuv * LVFU(rawData[i - 1][j]) + wtdv * LVFU(rawData[i + 1][j]) + wtlv * LVFU(rawData[i][j - 1]) + wtrv * LVFU(rawData[i][j + 1])) / (wtuv + wtdv + wtlv + wtrv);
                        _mm_store_ps(&greentile[rr * TS + cc], vself(selmask, greenv, tempv));
                        _mm_store_ps(&redtile[rr * TS + cc], tempv);
                        _mm_store_ps(&bluetile[rr * TS + cc], tempv);
                    }

                    for (; j < right; j++, cc++) {

                        if (FC(i, j) == 1) {
                            greentile[rr * TS + cc] = rawData[i][j];

                        } else {
                            //compute directional weights using image gradients
                            wtu = INVGRAD((abs(rawData[i + 1][j] - rawData[i - 1][j]) + abs(rawData[i][j] - rawData[i - 2][j]) + abs(rawData[i - 1][j] - rawData[i - 3][j])));
                            wtd = INVGRAD((abs(rawData[i - 1][j] - rawData[i + 1][j]) + abs(rawData[i][j] - rawData[i + 2][j]) + abs(rawData[i + 1][j] - rawData[i + 3][j])));
                            wtl = INVGRAD((abs(rawData[i][j + 1] - rawData[i][j - 1]) + abs(rawData[i][j] - rawData[i][j - 2]) + abs(rawData[i][j - 1] - rawData[i][j - 3])));
                            wtr = INVGRAD((abs(rawData[i][j - 1] - rawData[i][j + 1]) + abs(rawData[i][j] - rawData[i][j + 2]) + abs(rawData[i][j + 1] - rawData[i][j + 3])));

                            //store in rgb array the interpolated G value at R/B grid points using directional weighted average
                            greentile[rr * TS + cc] = (wtu * rawData[i - 1][j] + wtd * rawData[i + 1][j] + wtl * rawData[i][j - 1] + wtr * rawData[i][j + 1]) / (wtu + wtd + wtl + wtr);
                        }

                        redtile[rr * TS + cc] = rawData[i][j];
                        bluetile[rr * TS + cc] = rawData[i][j];
                    }

#else

                    for (int j = left, cc = 0; j < right; j++, cc++) {
                        if (FC(i, j) == 1) {
                            greentile[rr * TS + cc] = rawData[i][j];
                        } else {
                            //compute directional weights using image gradients
                            wtu = INVGRAD((abs(rawData[i + 1][j] - rawData[i - 1][j]) + abs(rawData[i][j] - rawData[i - 2][j]) + abs(rawData[i - 1][j] - rawData[i - 3][j])));
                            wtd = INVGRAD((abs(rawData[i - 1][j] - rawData[i + 1][j]) + abs(rawData[i][j] - rawData[i + 2][j]) + abs(rawData[i + 1][j] - rawData[i + 3][j])));
                            wtl = INVGRAD((abs(rawData[i][j + 1] - rawData[i][j - 1]) + abs(rawData[i][j] - rawData[i][j - 2]) + abs(rawData[i][j - 1] - rawData[i][j - 3])));
                            wtr = INVGRAD((abs(rawData[i][j - 1] - rawData[i][j + 1]) + abs(rawData[i][j] - rawData[i][j + 2]) + abs(rawData[i][j + 1] - rawData[i][j + 3])));

                            //store in rgb array the interpolated G value at R/B grid points using directional weighted average
                            greentile[rr * TS + cc] = (wtu * rawData[i - 1][j] + wtd * rawData[i + 1][j] + wtl * rawData[i][j - 1] + wtr * rawData[i][j + 1]) / (wtu + wtd + wtl + wtr);
                        }

                        redtile[rr * TS + cc] = rawData[i][j];
                        bluetile[rr * TS + cc] = rawData[i][j];
                    }

#endif
                }

#ifdef __SSE2__
                __m128 zd25v = _mm_set1_ps(0.25f);
                __m128 clip_ptv = _mm_set1_ps( clip_pt );
#endif

                for (int i = top + 1, rr = 1; i < bottom - 1; i++, rr++) {
                    if (FC(i, left + (FC(i, 2) & 1) + 1) == 0)
#ifdef __SSE2__
                        for (int j = left + 1, cc = 1; j < right - 1; j += 4, cc += 4) {
                            //interpolate B/R colors at R/B sites
                            _mm_storeu_ps(&bluetile[rr * TS + cc], LVFU(greentile[rr * TS + cc]) - zd25v * ((LVFU(greentile[(rr - 1)*TS + (cc - 1)]) + LVFU(greentile[(rr - 1)*TS + (cc + 1)]) + LVFU(greentile[(rr + 1)*TS + cc + 1]) + LVFU(greentile[(rr + 1)*TS + cc - 1])) -
                                          _mm_min_ps(clip_ptv, LVFU(rawData[i - 1][j - 1]) + LVFU(rawData[i - 1][j + 1]) + LVFU(rawData[i + 1][j + 1]) + LVFU(rawData[i + 1][j - 1]))));
                        }

#else

                        for (int cc = (FC(i, 2) & 1) + 1, j = left + cc; j < right - 1; j += 2, cc += 2) {
                            //interpolate B/R colors at R/B sites
                            bluetile[rr * TS + cc] = greentile[rr * TS + cc] - 0.25f * ((greentile[(rr - 1) * TS + (cc - 1)] + greentile[(rr - 1) * TS + (cc + 1)] + greentile[(rr + 1) * TS + cc + 1] + greentile[(rr + 1) * TS + cc - 1]) -
                                                     min(clip_pt, rawData[i - 1][j - 1] + rawData[i - 1][j + 1] + rawData[i + 1][j + 1] + rawData[i + 1][j - 1]));
                        }

#endif
                    else
#ifdef __SSE2__
                        for (int j = left + 1, cc = 1; j < right - 1; j += 4, cc += 4) {
                            //interpolate B/R colors at R/B sites
                            _mm_storeu_ps(&redtile[rr * TS + cc], LVFU(greentile[rr * TS + cc]) - zd25v * ((LVFU(greentile[(rr - 1)*TS + cc - 1]) + LVFU(greentile[(rr - 1)*TS + cc + 1]) + LVFU(greentile[(rr + 1)*TS + cc + 1]) + LVFU(greentile[(rr + 1)*TS + cc - 1])) -
                                          _mm_min_ps(clip_ptv, LVFU(rawData[i - 1][j - 1]) + LVFU(rawData[i - 1][j + 1]) + LVFU(rawData[i + 1][j + 1]) + LVFU(rawData[i + 1][j - 1]))));
                        }

#else

                        for (int cc = (FC(i, 2) & 1) + 1, j = left + cc; j < right - 1; j += 2, cc += 2) {
                            //interpolate B/R colors at R/B sites
                            redtile[rr * TS + cc] = greentile[rr * TS + cc] - 0.25f * ((greentile[(rr - 1) * TS + cc - 1] + greentile[(rr - 1) * TS + cc + 1] + greentile[(rr + 1) * TS + cc + 1] + greentile[(rr + 1) * TS + cc - 1]) -
                                                    min(clip_pt, rawData[i - 1][j - 1] + rawData[i - 1][j + 1] + rawData[i + 1][j + 1] + rawData[i + 1][j - 1]));
                        }

#endif
                }


#ifdef __SSE2__
                __m128 temp1v, temp2v, greensumv;
                selmask = _mm_set_epi32( 0xffffffff, 0, 0xffffffff, 0 );
#endif

                // interpolate R/B using color differences
                for (int i = top + 2, rr = 2; i < bottom - 2; i++, rr++) {
#ifdef __SSE2__

                    for (int cc = 2 + (FC(i, 2) & 1), j = left + cc; j < right - 2; j += 4, cc += 4) {
                        // no need to take care about the borders of the tile. There's enough free space.
                        //interpolate R and B colors at G sites
                        greenv = LVFU(greentile[rr * TS + cc]);
                        greensumv = LVFU(greentile[(rr - 1) * TS + cc]) + LVFU(greentile[(rr + 1) * TS + cc]) + LVFU(greentile[rr * TS + cc - 1]) + LVFU(greentile[rr * TS + cc + 1]);

                        temp1v = LVFU(redtile[rr * TS + cc]);
                        temp2v = greenv - zd25v * (greensumv - LVFU(redtile[(rr - 1) * TS + cc]) - LVFU(redtile[(rr + 1) * TS + cc]) - LVFU(redtile[rr * TS + cc - 1]) - LVFU(redtile[rr * TS + cc + 1]));

//              temp2v = greenv - zd25v*((LVFU(greentile[(rr-1)*TS+cc])-LVFU(redtile[(rr-1)*TS+cc]))+(LVFU(greentile[(rr+1)*TS+cc])-LVFU(redtile[(rr+1)*TS+cc]))+
//                                                         (LVFU(greentile[rr*TS+cc-1])-LVFU(redtile[rr*TS+cc-1]))+(LVFU(greentile[rr*TS+cc+1])-LVFU(redtile[rr*TS+cc+1])));
                        _mm_storeu_ps( &redtile[rr * TS + cc], vself(selmask, temp1v, temp2v));

                        temp1v = LVFU(bluetile[rr * TS + cc]);

                        temp2v = greenv - zd25v * (greensumv - LVFU(bluetile[(rr - 1) * TS + cc]) - LVFU(bluetile[(rr + 1) * TS + cc]) - LVFU(bluetile[rr * TS + cc - 1]) - LVFU(bluetile[rr * TS + cc + 1]));

//              temp2v = greenv - zd25v*((LVFU(greentile[(rr-1)*TS+cc])-LVFU(bluetile[(rr-1)*TS+cc]))+(LVFU(greentile[(rr+1)*TS+cc])-LVFU(bluetile[(rr+1)*TS+cc]))+
//                                                          (LVFU(greentile[rr*TS+cc-1])-LVFU(bluetile[rr*TS+cc-1]))+(LVFU(greentile[rr*TS+cc+1])-LVFU(bluetile[rr*TS+cc+1])));
                        _mm_storeu_ps( &bluetile[rr * TS + cc], vself(selmask, temp1v, temp2v));
                    }

#else

                    for (int cc = 2 + (FC(i, 2) & 1), j = left + cc; j < right - 2; j += 2, cc += 2) {
                        //interpolate R and B colors at G sites
                        redtile[rr * TS + cc] = greentile[rr * TS + cc] - 0.25f * ((greentile[(rr - 1) * TS + cc] - redtile[(rr - 1) * TS + cc]) + (greentile[(rr + 1) * TS + cc] - redtile[(rr + 1) * TS + cc]) +
                                                (greentile[rr * TS + cc - 1] - redtile[rr * TS + cc - 1]) + (greentile[rr * TS + cc + 1] - redtile[rr * TS + cc + 1]));
                        bluetile[rr * TS + cc] = greentile[rr * TS + cc] - 0.25f * ((greentile[(rr - 1) * TS + cc] - bluetile[(rr - 1) * TS + cc]) + (greentile[(rr + 1) * TS + cc] - bluetile[(rr + 1) * TS + cc]) +
                                                 (greentile[rr * TS + cc - 1] - bluetile[rr * TS + cc - 1]) + (greentile[rr * TS + cc + 1] - bluetile[rr * TS + cc + 1]));
                    }

#endif
                }


                for (int i = top + 2, rr = 2; i < bottom - 2; i++, rr++) {
#ifdef __SSE2__

                    for (j = left + 2, cc = 2; j < right - 5; j += 4, cc += 4) {
                        _mm_storeu_ps(&red[i][j], LVFU(redtile[rr * TS + cc]));
                        _mm_storeu_ps(&green[i][j], LVFU(greentile[rr * TS + cc]));
                        _mm_storeu_ps(&blue[i][j], LVFU(bluetile[rr * TS + cc]));
                    }

                    for (; j < right - 2; j++, cc++) {
                        red[i][j] = redtile[rr * TS + cc];
                        green[i][j] = greentile[rr * TS + cc];
                        blue[i][j] = bluetile[rr * TS + cc];
                    }

#else

                    for (int j = left + 2, cc = 2; j < right - 2; j++, cc++) {
                        red[i][j] = redtile[rr * TS + cc];
                        green[i][j] = greentile[rr * TS + cc];
                        blue[i][j] = bluetile[rr * TS + cc];
                    }

#endif


                }

                if(plistenerActive && ((++progressCounter) % 16 == 0)) {
#ifdef _OPENMP
                    #pragma omp critical (updateprogress)
#endif
                    {
                        progress += progressInc;
                        progress = min(1.0, progress);
                        plistener->setProgress (progress);
                    }
                }

            }

        free(buffer);
    } // End of parallelization

    if(plistenerActive) {
        plistener->setProgress(1.00);
    }



}
#undef TS
#undef CLF
