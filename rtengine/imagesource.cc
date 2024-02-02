/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2024 RawTherapee team
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

#include "imagesource.h"
#include "procparams.h"


namespace rtengine
{

void ImageSource::getAutoWBMultipliersItcGreen(
        procparams::ProcParams &params,
        bool forcewbgrey,
        int kcam,
        double greenitc,
        bool extra,
        float &temp0,
        float &delta,
        int &bia,
        int &dread,
        int nocam,
        float &studgood,
        float &minchrom,
        int &kmin,
        float &minhist,
        float &maxhist,
        int fh,
        int fw,
        ColorTemp &currWB,
        int tempnotisraw,
        double greennotisraw,
        bool skipRecalculate,
        ColorTemp &autoWB,
        double &rm,
        double &gm,
        double &bm
        )
{
    float tem = 5000.f;
    float gre  = 1.f;
    double tempref0bias = 5000.;
    double tempitc = 5000.f;
    bool autowb1 = true;
    double green_thres = 0.8;

    if (isRAW()) {// only with Raw files

        auto currWBitc = getWB();

        double greenref = currWBitc.getGreen();
        double tempref0bias0 = currWBitc.getTemp();

        if (greenref > green_thres && params.wb.itcwb_prim == "srgb") {
            forcewbgrey = true;
        }

        if (!forcewbgrey && (tempref0bias0 < 3300.f)  && (greenref < 1.13f && greenref > 0.88f)) { //seems good with temp and green...To fixe...limits 1.13 and 0.88
            if (settings->verbose) {
                printf("Keep camera settings temp=%f green=%f\n", tempref0bias0, greenref);
            }

            autowb1 = true;
            kcam = 1;
        }

        if (autowb1) {
            //alternative to camera if camera settings out, using autowb grey to find new ref, then mixed with camera
            // kcam = 0;
            params.wb.method = "autold";
            tempitc = 5000.f;
            greenitc = 1.;
            currWBitc = getWB();
            tempref0bias = currWBitc.getTemp();
            double greenref = currWBitc.getGreen();
            bool pargref = true;
            bool pargre = true;

            if ((greenref > 1.5f || tempref0bias < 3300.f || tempref0bias > 7700.f || forcewbgrey) && kcam != 1 && !params.wb.itcwb_sampling) { //probably camera out to adjust...
                getAutoWBMultipliersitc(extra, tempref0bias, greenref, tempitc, greenitc, temp0, delta, bia, dread, kcam, nocam, studgood, minchrom, kmin, minhist, maxhist, 0, 0, fh, fw, 0, 0, fh, fw, rm, gm, bm,  params.wb, params.icm, params.raw, params.toneCurve);
                wbMul2Camera(rm, gm, bm);
                wbCamera2Mul(rm, gm, bm);
                ColorTemp ct(rm, gm, bm, 1.0, currWB.getObserver());
                tem = ct.getTemp();
                gre  = ct.getGreen();

                if (gre > 1.3f) {
                    pargre = false;
                }

                if (greenref > 1.3f) {
                    pargref = false;
                }

                double deltemp = tem - tempref0bias;

                if (gre > 1.5f && !forcewbgrey) { //probable wrong value
                    tem = 0.3 * tem + 0.7 * tempref0bias;//find a mixed value
                    gre = 0.5f + 0.5f * LIM(gre, 0.9f, 1.1f);//empirical formula in case  system out
                } else {
                    if (!forcewbgrey) {
                        gre = 0.2f + 0.8f * LIM(gre, 0.85f, 1.15f);
                        tem = 0.3 * tem + 0.7 * tempref0bias;//find a mixed value
                        nocam = 0;
                    } else {//set temp and green to init itcwb algorithm
                        double grepro = LIM(greenref, green_thres, 1.15);
                        gre = 0.5f * grepro + 0.5f * LIM(gre, 0.9f, 1.1f);//empirical green between green camera and autowb grey

                        if (abs(deltemp) < 400.) { //arbitraries thresholds to refine
                            tem = 0.3 * tem + 0.7 * tempref0bias;//find a mixed value between camera and auto grey

                            if (deltemp > 0.) {
                                nocam = 1;
                            } else {
                                nocam = 2;
                            }
                        } else if (abs(deltemp) < 900.) { //other arbitrary threshold
                            tem = 0.4 * tem + 0.6 * tempref0bias;//find a mixed value between camera and auto grey

                            if (deltemp > 0.) {
                                nocam = 3;
                            } else {
                                nocam = 4;
                            }
                        } else if (abs(deltemp) < 1500. && tempref0bias < 4500.f) {
                            if ((pargre && pargref) || (!pargre && !pargref)) {
                                tem = 0.45 * tem + 0.55 * tempref0bias;//find a mixed value between camera and auto grey
                            }

                            if (pargre && !pargref) {
                                tem = 0.7 * tem + 0.3 * tempref0bias;//find a mixed value between camera and auto grey
                            }

                            if (!pargre && pargref) {
                                tem = 0.3 * tem + 0.7 * tempref0bias;//find a mixed value between camera and auto grey
                            }

                            nocam = 5;
                        } else if (abs(deltemp) < 1500. && tempref0bias >= 4500.f) {
                            if ((pargre && pargref) || (!pargre && !pargref)) {
                                tem = 0.45 * tem + 0.55 * tempref0bias;//find a mixed value between camera and auto grey
                            }

                            if (pargre && !pargref) {
                                tem = 0.7 * tem + 0.3 * tempref0bias;//find a mixed value between camera and auto grey
                            }

                            if (!pargre && pargref) {
                                tem = 0.3 * tem + 0.7 * tempref0bias;//find a mixed value between camera and auto grey
                            }

                            nocam = 6;
                        } else if (abs(deltemp) >= 1500. && tempref0bias < 5500.f) {
                            if (tem >= 4500.f) {
                                if ((pargre && pargref) || (!pargre && !pargref)) {
                                    tem = 0.7 * tem + 0.3 * tempref0bias;//find a mixed value between camera and auto grey
                                }

                                if (pargre && !pargref) {
                                    tem = 0.8 * tem + 0.2 * tempref0bias;//find a mixed value between camera and auto grey
                                }

                                if (!pargre && pargref) {
                                    tem = 0.3 * tem + 0.7 * tempref0bias;//find a mixed value between camera and auto grey
                                }

                                nocam = 7;
                            } else {
                                tem = 0.3 * tem + 0.7 * tempref0bias;//find a mixed value between camera and auto grey
                                nocam = 8;
                            }
                        } else if (abs(deltemp) >= 1500. && tempref0bias >= 5500.f) {
                            if (tem >= 10000.f) {
                                tem = 0.99 * tem + 0.01 * tempref0bias;//find a mixed value between camera and auto grey
                                nocam = 9;
                            } else {
                                if ((pargre && pargref) || (!pargre && !pargref)) {
                                    tem = 0.45 * tem + 0.55 * tempref0bias;//find a mixed value between camera and auto grey
                                }

                                if (pargre && !pargref) {
                                    tem = 0.7 * tem + 0.3 * tempref0bias;//find a mixed value between camera and auto grey
                                }

                                if (!pargre && pargref) {
                                    tem = 0.3 * tem + 0.7 * tempref0bias;//find a mixed value between camera and auto grey
                                }

                                nocam = 10;
                            }
                        } else {
                            tem = 0.4 * tem + 0.6 * tempref0bias;
                            nocam = 11;
                        }
                    }
                }

                tempitc = tem ;

                extra = true;

                if (settings->verbose) {
                    printf("Using new references AWB grey or mixed  Enable Extra - temgrey=%f gregrey=%f tempitc=%f nocam=%i\n", (double) tem, (double) gre, (double) tempitc, nocam);
                }
            }
        }

        params.wb.method = "autitcgreen";

    } else if (!isRAW()) { // Itcwb and no raw
        params.wb.temperature = tempnotisraw;
        params.wb.green = greennotisraw;
        params.wb.equal = 1.;
    }
    float greenitc_low = 1.f;
    float tempitc_low = 5000.f;
    //raw files and autitcgreen
    if (isRAW() || !skipRecalculate) {
        greenitc = 1.;
        auto currWBitc = getWB();
        currWBitc = currWBitc.convertObserver(params.wb.observer);//change the temp/green couple with the same multipliers

        double tempref = currWBitc.getTemp() * (1. + params.wb.tempBias);
        double greenref = currWBitc.getGreen();
        greenitc = greenref;

        if ((greenref > 1.5f || tempref0bias < 3300.f || tempref0bias > 7700.f || forcewbgrey) && autowb1 && kcam != 1 && !params.wb.itcwb_sampling) { //probably camera out to adjust = greenref ? tempref0bias ?
            tempref = tem * (1. + params.wb.tempBias);
            greenref = gre;
        } else {

        }

        if(params.wb.itcwb_sampling) {
            greenitc_low = greenref;
            tempitc_low = tempref;
        }

        if (settings->verbose) {
            printf("tempref=%f greref=%f tempitc=%f greenitc=%f\n", tempref, greenref, tempitc, greenitc);
        }

        getAutoWBMultipliersitc(extra, tempref, greenref, tempitc, greenitc, temp0, delta,  bia, dread, kcam, nocam, studgood, minchrom, kmin, minhist, maxhist, 0, 0, fh, fw, 0, 0, fh, fw, rm, gm, bm,  params.wb, params.icm, params.raw, params.toneCurve);

        params.wb.temperature = tempitc;
        params.wb.green = greenitc;
        if(params.wb.itcwb_sampling) {
            params.wb.temperature = tempitc_low;
            params.wb.green = greenitc_low;
        }

        currWB = ColorTemp(params.wb.temperature, params.wb.green, 1., params.wb.method, params.wb.observer);
        currWB.getMultipliers(rm, gm, bm);
        autoWB.update(rm, gm, bm, params.wb.equal, params.wb.observer, 0.); //params.wb.tempBias already used before

    }
}

} // namespace rtengine

