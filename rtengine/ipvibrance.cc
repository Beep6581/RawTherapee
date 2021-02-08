/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2011 Jacques Desmis  <jdesmis@gmail.com>
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

#include "rt_math.h"

#include "rtengine.h"
#include "improcfun.h"
#include "iccstore.h"
#include "labimage.h"
#include "curves.h"
#include "color.h"
#include "procparams.h"
#include "StopWatch.h"

using namespace std;

namespace rtengine
{

using namespace procparams;

void fillCurveArrayVib (DiagonalCurve* diagCurve, LUTf &outCurve)
{

    if (diagCurve) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i <= 0xffff; i++ ) {
            // change to [0,1] range
            // apply custom/parametric/NURBS curve, if any
            // and store result in a temporary array
            outCurve[i] = 65535.0 * diagCurve->getVal(i / 65535.0);
        }
    } else {
        outCurve.makeIdentity();
    }
}


/*
 * Vibrance correction
 * copyright (c)2011  Jacques Desmis <jdesmis@gmail.com> and Jean-Christophe Frisch <natureh@free.fr>
 *
 */
void ImProcFunctions::vibrance (LabImage* lab, const procparams::VibranceParams &vibranceParams, bool highlight, const Glib::ustring &workingProfile)
{
    if (!vibranceParams.enabled) {
        return;
    }
    BENCHFUN

//  int skip=1; //scale==1 ? 1 : 16;
    bool skinCurveIsSet = false;
    DiagonalCurve* dcurve = new DiagonalCurve (vibranceParams.skintonescurve, CURVES_MIN_POLY_POINTS);

    if (dcurve) {
        if (!dcurve->isIdentity()) {
            skinCurveIsSet = true;
        } else {
            delete dcurve;
            dcurve = nullptr;
        }
    }

    if (!skinCurveIsSet && !vibranceParams.pastels && !vibranceParams.saturated) {
        if (dcurve) {
            delete dcurve;
            dcurve = nullptr;
        }

        return;
    }

    const int width = lab->W;
    const int height = lab->H;

    // skin hue curve
    // I use diagonal because I think it's better

    const float chromaPastel = vibranceParams.pastels / 100.f;
    const float chromaSatur = vibranceParams.saturated / 100.f;
    constexpr float p00 = 0.07f;
    const float limitpastelsatur = (static_cast<float>(vibranceParams.psthreshold.getTopLeft()) / 100.f) * (1.f - p00) + p00;
    const float maxdp = (limitpastelsatur - p00) / 4.f;
    const float maxds = (1.f - limitpastelsatur) / 4.f;
    const float p0 = p00 + maxdp;
    const float p1 = p00 + 2.f * maxdp;
    const float p2 = p00 + 3.f * maxdp;
    const float s0 = limitpastelsatur + maxds;
    const float s1 = limitpastelsatur + 2.f * maxds;
    const float s2 = limitpastelsatur + 3.f * maxds;
    const float transitionweighting = static_cast<float>(vibranceParams.psthreshold.getBottomLeft()) / 100.f;
    float chromamean = 0.f;

    if (chromaPastel != chromaSatur) {
        //if sliders pastels and saturated are different: transition with a double linear interpolation: between p2 and limitpastelsatur, and between limitpastelsatur and s0
        //modify the "mean" point in function of double threshold  => differential transition
        chromamean = maxdp * (chromaSatur - chromaPastel) / (s0 - p2) + chromaPastel;

        // move chromaMean up or down depending on transitionCtrl
        if (transitionweighting > 0.f) {
            chromamean = (chromaSatur - chromamean) * transitionweighting + chromamean;
        } else if (transitionweighting < 0.f) {
            chromamean = (chromamean - chromaPastel)  * transitionweighting + chromamean;
        }
    }

    const float chromaPastel_a = (chromaPastel - chromamean) / (p2 - limitpastelsatur);
    const float chromaPastel_b = chromaPastel - chromaPastel_a * p2;

    const float chromaSatur_a = (chromaSatur - chromamean) / (s0 - limitpastelsatur);
    const float chromaSatur_b = chromaSatur - chromaSatur_a * s0;

    constexpr float dhue = 0.15f; //hue transition
    constexpr float dchr = 20.f; //chroma transition
    constexpr float skbeg = -0.05f; //begin hue skin
    constexpr float skend = 1.60f; //end hue skin
    constexpr float xx = 0.5f; //soft : between 0.3 and 1.0
    constexpr float ask = 65535.f / (skend - skbeg);
    constexpr float bsk0 = -skbeg;
    constexpr float bsk = -skbeg * ask;

    LUTf skin_curve (65536, 0);

    if (skinCurveIsSet) {
        fillCurveArrayVib (dcurve, skin_curve);
        skin_curve /= ask;
//        skin_curve *= 2.f;
    }

    if (dcurve) {
        delete dcurve;
        dcurve = nullptr;
    }

    const bool protectskins = vibranceParams.protectskins;
    const bool avoidcolorshift = vibranceParams.avoidcolorshift;

    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (workingProfile);
    //inverse matrix user select
    const float wip[3][3] = {
        {static_cast<float>(wiprof[0][0]), static_cast<float>(wiprof[0][1]), static_cast<float>(wiprof[0][2])},
        {static_cast<float>(wiprof[1][0]), static_cast<float>(wiprof[1][1]), static_cast<float>(wiprof[1][2])},
        {static_cast<float>(wiprof[2][0]), static_cast<float>(wiprof[2][1]), static_cast<float>(wiprof[2][2])}
    };

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {

#ifdef __SSE2__
        float HHbuffer[width] ALIGNED16;
        float CCbuffer[width] ALIGNED16;
#endif
        float sathue[5], sathue2[4]; // adjust sat in function of hue

#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 16)
#endif

        for (int i = 0; i < height; i++) {
#ifdef __SSE2__
            // vectorized per row calculation of HH and CC
            vfloat c327d68v = F2V(327.68f);
            int k = 0;
            for (; k < width - 3; k += 4) {
                vfloat av = LVFU(lab->a[i][k]);
                vfloat bv = LVFU(lab->b[i][k]);
                STVF(HHbuffer[k], xatan2f(bv, av));
                STVF(CCbuffer[k], vsqrtf(SQRV(av) + SQRV(bv)) / c327d68v);
            }
            for (; k < width; k++) {
                HHbuffer[k] = xatan2f (lab->b[i][k], lab->a[i][k]);
                CCbuffer[k] = sqrt (SQR (lab->a[i][k]) + SQR (lab->b[i][k])) / 327.68f;
            }
#endif
            for (int j = 0; j < width; j++) {
                float LL = lab->L[i][j] / 327.68f;
#ifdef __SSE2__
                float HH = HHbuffer[j];
                float CC = CCbuffer[j];
#else
                float HH = xatan2f (lab->b[i][j], lab->a[i][j]);
                float CC = sqrt (SQR (lab->a[i][j]) + SQR (lab->b[i][j])) / 327.68f;
#endif

                // here we work on Chromaticity and Hue
                // variation of Chromaticity  ==> saturation via RGB
                // Munsell correction, then conversion to Lab
                float Lprov = LL;
                float Chprov = CC;
                float2 sincosval;

                if (CC == 0.f) {
                    sincosval.y = 1.f;
                    sincosval.x = 0.f;
                } else {
                    sincosval.y = lab->a[i][j] / (CC * 327.68f);
                    sincosval.x = lab->b[i][j] / (CC * 327.68f);
                }

                //gamut control : Lab values are in gamut
                float saturation;
                Color::gamutLchonly(HH, sincosval, Lprov, Chprov, saturation, wip, highlight, 0.15f, 0.98f);

                if (Chprov > 6.f) {
                    float satredu = 1.f; //reduct sat in function of skin

                    if (protectskins) {
                        Color::SkinSat (LL, HH, CC, satredu);// for skin colors
                    }

                    if (saturation > 0.f) {
                        if (satredu != 1.f) {
                            // for skin, no differentiation
                            sathue [0] = sathue [1] = sathue [2] = sathue [3] = sathue[4] = 1.f;
                            sathue2[0] = sathue2[1] = sathue2[2] = sathue2[3]          = 1.f;
                        } else {
                            //double pyramid: LL and HH
                            //I try to take into account: Munsell response (human vision) and Gamut..(less response for red): preferably using Prophoto or WideGamut
                            //blue: -1.80 -3.14  green = 2.1 3.14   green-yellow=1.4 2.1  red:0 1.4  blue-purple:-0.7  -1.4   purple: 0 -0.7
                            //these values allow a better and differential response
                            if (LL < 20.0f) { //more for blue-purple, blue and red modulate
                                sathue[4] = 0.4f;
                                sathue2[3] = 1.f;
                                if     (/*HH> -3.1415f &&*/ HH < -1.5f   ) {
                                    sathue[0] = 1.3f;    //blue
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue2[0] = 1.05f;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f   ) {
                                    sathue[0] = 1.6f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.15f;
                                    sathue2[2] = 1.1f ;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f   ) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=1.1f;sathue[1]=1.1f;sathue[2]=1.1f;sathue[3]=1.0f ;sathue[4]=0.4f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f   ) {
                                    sathue[0] = 1.3f;    //red   0.8 0.7
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.0f ;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f   ) {
                                    sathue[0] = 1.0f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.4f;    //green
                                    sathue[1] = 1.3f;
                                    sathue[2] = 1.2f;
                                    sathue[3] = 1.15f;
                                    sathue2[0] = 1.15f;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                }
                            } else if (LL < 50.0f) { //more for blue and green, less for red and green-yellow
                                sathue[4] = 0.4f;
                                if     (/*HH> -3.1415f &&*/ HH < -1.5f   ) {
                                    sathue[0] = 1.5f;    //blue
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f   ) {
                                    sathue[0] = 1.3f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue2[0] = 1.05f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f   ) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 1.0f ;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f;sathue[3]=0.8f ;sathue[4]=0.4f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f   ) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f ;
                                    sathue2[0] = 0.8f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f   ) {
                                    sathue[0] = 1.1f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.1f;
                                    sathue[2] = 1.1f;
                                    sathue[3] = 1.05f;
                                    sathue2[0] = 0.9f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.7f ;
                                    sathue2[3] = 0.6f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.5f;    //green
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f;
                                    sathue[3] = 1.2f ;
                                    sathue2[0] = 1.2f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                }

                            } else if (LL < 80.0f) { //more for green, less for red and green-yellow
                                sathue[4] = 0.3f;
                                if     (/*HH> -3.1415f &&*/ HH < -1.5f   ) {
                                    sathue[0] = 1.3f;    //blue
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.15f;
                                    sathue[3] = 1.1f ;
                                    sathue2[0] = 1.1f ;
                                    sathue2[1] = 1.1f ;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f   ) {
                                    sathue[0] = 1.3f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.15f;
                                    sathue[3] = 1.1f ;
                                    sathue2[0] = 1.1f ;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f   ) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f ;
                                    sathue[3] = 1.0f ;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 1.0f ;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f ;sathue[3]=0.8f ;sathue[4]=0.3f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f   ) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f ;
                                    sathue[3] = 0.8f ;
                                    sathue2[0] = 0.8f ;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f   ) {
                                    sathue[0] = 1.3f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.2f;
                                    sathue[2] = 1.1f ;
                                    sathue[3] = 1.05f;
                                    sathue2[0] = 1.0f ;
                                    sathue2[1] = 0.9f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.7f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.6f;    //green - even with Prophoto green are too "little"  1.5 1.3
                                    sathue[1] = 1.4f;
                                    sathue[2] = 1.3f ;
                                    sathue[3] = 1.25f;
                                    sathue2[0] = 1.25f;
                                    sathue2[1] = 1.2f ;
                                    sathue2[2] = 1.15f;
                                    sathue2[3] = 1.05f;
                                }
                            } else { /*if (LL>=80.0f)*/ //more for green-yellow, less for red and purple
                                sathue[4] = 0.2f;
                                if     (/*HH> -3.1415f &&*/ HH < -1.5f   ) {
                                    sathue[0] = 1.0f;    //blue
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>=-1.5f    &&*/ HH < -0.7f   ) {
                                    sathue[0] = 1.0f;    //blue purple  1.2 1.1
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>=-0.7f    &&*/ HH <  0.0f   ) {
                                    sathue[0] = 1.2f;    //purple
                                    sathue[1] = 1.0f;
                                    sathue[2] = 1.0f;
                                    sathue[3] = 0.9f;
                                    sathue2[0] = 0.9f;
                                    sathue2[1] = 0.9f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                }
                                //          else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f;sathue[3]=0.8f;sathue[4]=0.2f;sathue2[0]=0.8f;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
                                else if (/*HH>= 0.0f    &&*/ HH <= 1.4f   ) {
                                    sathue[0] = 1.1f;    //red   0.8 0.7
                                    sathue[1] = 1.0f;
                                    sathue[2] = 0.9f;
                                    sathue[3] = 0.8f;
                                    sathue2[0] = 0.8f;
                                    sathue2[1] = 0.8f ;
                                    sathue2[2] = 0.8f ;
                                    sathue2[3] = 0.8f;
                                } else if (/*HH>  1.4f    &&*/ HH <= 2.1f   ) {
                                    sathue[0] = 1.6f;    //green yellow 1.2 1.1
                                    sathue[1] = 1.5f;
                                    sathue[2] = 1.4f;
                                    sathue[3] = 1.2f;
                                    sathue2[0] = 1.1f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.0f ;
                                    sathue2[3] = 1.0f;
                                } else { /*if(HH>  2.1f    && HH<= 3.1415f)*/
                                    sathue[0] = 1.4f;    //green
                                    sathue[1] = 1.3f;
                                    sathue[2] = 1.2f;
                                    sathue[3] = 1.1f;
                                    sathue2[0] = 1.1f;
                                    sathue2[1] = 1.05f;
                                    sathue2[2] = 1.05f;
                                    sathue2[3] = 1.0f;
                                }
                            }
                        }

                        float chmod = 0.f;
                        // variables to improve transitions
                        float pa, pb;// transition = pa*saturation + pb

                        // We handle only positive values here ;  improve transitions
                        if      (saturation < p00) {
                            pa = 0.f;
                            pb = chromaPastel * sathue[4];
                        } else if (saturation < p0 )               {
                            float chl00 = chromaPastel * sathue[4];
                            float chl0  = chromaPastel * sathue[0];
                            pa = (chl00 - chl0) / (p00 - p0);
                            pb = chl00 - pa * p00;
                        } else if (saturation < p1)                {
                            float chl0  = chromaPastel * sathue[0];
                            float chl1  = chromaPastel * sathue[1];
                            pa = (chl0 - chl1) / (p0 - p1);
                            pb = chl0 - pa * p0;
                        } else if (saturation < p2)                {
                            float chl1  = chromaPastel * sathue[1];
                            float chl2  = chromaPastel * sathue[2];
                            pa = (chl1 - chl2) / (p1 - p2);
                            pb = chl1 - pa * p1;
                        } else if (saturation < limitpastelsatur)  {
                            float chl2  = chromaPastel *  sathue[2];
                            float chl3  = chromaPastel * sathue[3];
                            pa = (chl2 - chl3) / (p2 - limitpastelsatur);
                            pb = chl2 - pa * p2;
                        } else if (saturation < s0)                {
                            float chl3  = chromaPastel * sathue[3];
                            float chs0  = chromaSatur * sathue2[0];
                            pa = (chl3 - chs0) / (limitpastelsatur - s0) ;
                            pb = chl3 - pa * limitpastelsatur;
                        } else if (saturation < s1)                {
                            float chs0  = chromaSatur * sathue2[0];
                            float chs1  = chromaSatur * sathue2[1];
                            pa = (chs0 - chs1) / (s0 - s1);
                            pb = chs0 - pa * s0;
                        } else if (saturation < s2)                {
                            float chs1  = chromaSatur * sathue2[1];
                            float chs2  = chromaSatur * sathue2[2];
                            pa = (chs1 - chs2) / (s1 - s2);
                            pb = chs1 - pa * s1;
                        } else                                     {
                            float chs2  = chromaSatur * sathue2[2];
                            float chs3  = chromaSatur * sathue2[3];
                            pa = (chs2 - chs3) / (s2 - 1.f);
                            pb = chs2 - pa * s2;
                        }
                        chmod = pa * saturation + pb;
                        chmod *= satredu;

                        if (chromaPastel != chromaSatur) {

                            // Pastels
                            if (saturation > p2 && saturation < limitpastelsatur) {
                                float newchromaPastel = chromaPastel_a * saturation + chromaPastel_b;
                                chmod = newchromaPastel * satredu * sathue[3];
                            }

                            // Saturated
                            if (saturation < s0 && saturation >= limitpastelsatur) {
                                float newchromaSatur = chromaSatur_a * saturation + chromaSatur_b;
                                chmod = newchromaSatur * satredu * sathue2[0];
                            }
                        }// end transition

                        if (saturation <= limitpastelsatur) {
                            chmod = rtengine::LIM(chmod, -0.93f, 2.f);
                            Chprov *= 1.0f + chmod;

                        } else { //if (saturation > limitpastelsatur)
                            chmod = rtengine::LIM(chmod, -0.93f, 1.8f);
                            Chprov *= 1.0f + chmod;

                        }
                        Chprov = rtengine::max(Chprov, 6.f);
                    }
                }

                bool hhModified = false;

                // Vibrance's Skin curve
                if (skinCurveIsSet) {
                    if (HH > skbeg && HH < skend) {
                        if (Chprov < 60.f) { //skin hue  : todo ==> transition
                            float HHsk = ask * HH + bsk;
                            float Hn = skin_curve[HHsk] - bsk0;
                            float Hc = Hn * xx + HH * (1.f - xx);
                            HH = Hc;
                            hhModified = true;
                        } else if (Chprov < (60.f + dchr)) { //transition chroma
                            float HHsk = ask * HH + bsk;
                            float Hn = skin_curve[HHsk] - bsk0;
                            float Hc = Hn * xx + HH * (1.f - xx);
                            float aa = (HH - Hc) / dchr ;
                            float bb = HH - (60.f + dchr) * aa;
                            HH = aa * Chprov + bb;
                            hhModified = true;
                        }
                    }
                    //transition hue
                    else if (HH > (skbeg - dhue) && HH <= skbeg && Chprov < (60.f + dchr * 0.5f)) {
                        float HHsk = ask * skbeg + bsk;
                        float Hn = skin_curve[HHsk] - bsk0;
                        float Hcc = Hn * xx + skbeg * (1.f - xx);
                        float adh = (Hcc - (skbeg - dhue)) / dhue;
                        float bdh = Hcc - adh * skbeg;
                        HH = adh * HH + bdh;
                        hhModified = true;
                    } else if (HH >= skend && HH < (skend + dhue) && Chprov < (60.f + dchr * 0.5f)) {
                        float HHsk = ask * skend + bsk;
                        float Hn = skin_curve[HHsk] - bsk0;
                        float Hcc = Hn * xx + skend * (1.f - xx);
                        float adh = (skend + dhue - Hcc) / dhue;
                        float bdh = Hcc - adh * skend;
                        HH = adh * HH + bdh;
                        hhModified = true;
                    }
                } // end skin hue

                //Munsell correction
                if (!avoidcolorshift && hhModified) {
                    sincosval = xsincosf (HH);
                }

                float aprovn, bprovn;
                bool inGamut;

                const float fyy = Color::c1By116 * Lprov + Color::c16By116;
                const float yy_ = (Lprov > static_cast<float>(Color::epskap)) ? fyy * fyy*fyy : Lprov / Color::kappaf;
                float ChprovOld = std::numeric_limits<float>::min();
                do {
                    inGamut = true;

                    if (avoidcolorshift) {
                        float correctionHue = 0.0f;

                        Color::AllMunsellLch(Lprov, HH, Chprov, CC, correctionHue);

                        if (correctionHue != 0.f || hhModified) {
                            sincosval = xsincosf(HH + correctionHue);
                            hhModified = false;
                        }
                    }
                    aprovn = Chprov * sincosval.y;
                    bprovn = Chprov * sincosval.x;

                    if (Chprov == ChprovOld) { // avoid endless loop
                        break;
                    } else {
                        ChprovOld = Chprov;
                    }

                    float fxx = 0.002f * aprovn + fyy;
                    float fzz = fyy - 0.005f * bprovn;
                    float xx_ = Color::f2xyz(fxx) * Color::D50x;
                    float zz_ = Color::f2xyz(fzz) * Color::D50z;

                    float R, G, B;
                    Color::xyz2rgb (xx_, yy_, zz_, R, G, B, wip);

                    if (rtengine::min(R, G, B) < 0.0f) {
                        Chprov *= 0.98f;
                        inGamut = false;
                    }

                    // if "highlight reconstruction" enabled don't control Gamut for highlights
                    if (!highlight && max(R, G, B) > 1.f && min(R, G, B) <= 1.f) {
                        Chprov *= 0.98f;
                        inGamut = false;
                    }
                } while (!inGamut);

                //put new values in Lab
                lab->L[i][j] = Lprov * 327.68f;
                lab->a[i][j] = aprovn * 327.68f;
                lab->b[i][j] = bprovn * 327.68f;
            }
        }
    } // end of parallelization
}


}
