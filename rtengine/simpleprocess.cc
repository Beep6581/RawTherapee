/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "rtengine.h"
#include "colortemp.h"
#include "imagesource.h"
#include "improcfun.h"
#include "curves.h"
#include "iccstore.h"
#include "clutstore.h"
#include "processingjob.h"
#include <glibmm.h>
#include "../rtgui/options.h"
#include "rawimagesource.h"
#include "../rtgui/multilangmgr.h"
#include "mytime.h"
#undef THREAD_PRIORITY_NORMAL

namespace rtengine
{
extern const Settings* settings;

IImage16* processImage (ProcessingJob* pjob, int& errorCode, ProgressListener* pl, bool tunnelMetaData, bool flush)
{

    errorCode = 0;

    ProcessingJobImpl* job = static_cast<ProcessingJobImpl*>(pjob);

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_PROCESSING");
        pl->setProgress (0.0);
    }

    InitialImage* ii = job->initialImage;

    if (!ii) {
        ii = InitialImage::load (job->fname, job->isRaw, &errorCode);

        if (errorCode) {
            delete job;
            return NULL;
        }
    }

    procparams::ProcParams& params = job->pparams;

    // acquire image from imagesource
    ImageSource* imgsrc = ii->getImageSource ();

    int tr = getCoarseBitMask(params.coarse);
    int fw, fh;
    imgsrc->getFullSize (fw, fh, tr);

    // check the crop params
    if (params.crop.x > fw || params.crop.y > fh) {
        // the crop is completely out of the image, so we disable the crop
        params.crop.enabled = false;
        // and we set the values to the defaults
        params.crop.x = 0;
        params.crop.y = 0;
        params.crop.w = fw;
        params.crop.h = fh;
    } else {
        if (params.crop.x < 0) {
            params.crop.x = 0;
        }

        if (params.crop.y < 0) {
            params.crop.y = 0;
        }

        if ((params.crop.x + params.crop.w) > fw) {
            // crop overflow in the width dimension ; we trim it
            params.crop.w = fw - params.crop.x;
        }

        if ((params.crop.y + params.crop.h) > fh) {
            // crop overflow in the height dimension ; we trim it
            params.crop.h = fh - params.crop.y;
        }
    }

//    MyTime t1,t2;
//    t1.set();

    ImProcFunctions ipf (&params, true);

    PreviewProps pp (0, 0, fw, fh, 1);
    imgsrc->preprocess( params.raw, params.lensProf, params.coarse, params.dirpyrDenoise.enabled);

    if (params.toneCurve.autoexp) {// this enabled HLRecovery
        LUTu histRedRaw(256), histGreenRaw(256), histBlueRaw(256);
        imgsrc->getRAWHistogram(histRedRaw, histGreenRaw, histBlueRaw);

        if (ToneCurveParams::HLReconstructionNecessary(histRedRaw, histGreenRaw, histBlueRaw) && !params.toneCurve.hrenabled) {
            params.toneCurve.hrenabled = true;
            // WARNING: Highlight Reconstruction is being forced 'on', should we force a method here too?
        }
    }

    if (pl) {
        pl->setProgress (0.20);
    }

    imgsrc->demosaic( params.raw);

    if (pl) {
        pl->setProgress (0.30);
    }

    if(params.retinex.enabled) { //enabled Retinex
        LUTf cdcurve (65536, 0);
        LUTf mapcurve (65536, 0);
        LUTu dummy;
        RetinextransmissionCurve dehatransmissionCurve;
        RetinexgaintransmissionCurve dehagaintransmissionCurve;
        bool dehacontlutili = false;
        bool mapcontlutili = false;
        bool useHsl = false;
//        multi_array2D<float, 3> conversionBuffer(1, 1);
        multi_array2D<float, 4> conversionBuffer(1, 1);
        imgsrc->retinexPrepareBuffers(params.icm, params.retinex, conversionBuffer, dummy);
        imgsrc->retinexPrepareCurves(params.retinex, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, dehacontlutili, mapcontlutili, useHsl, dummy, dummy );
        float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
        imgsrc->retinex( params.icm, params.retinex, params.toneCurve, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, conversionBuffer, dehacontlutili, mapcontlutili, useHsl, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax, dummy);
    }

    if (pl) {
        pl->setProgress (0.40);
    }

    imgsrc->HLRecovery_Global( params.toneCurve );


    if (pl) {
        pl->setProgress (0.45);
    }

    // set the color temperature
    ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

    if (params.wb.method == "Camera") {
        currWB = imgsrc->getWB ();
    } else if (params.wb.method == "Auto") {
        double rm, gm, bm;
        imgsrc->getAutoWBMultipliers(rm, gm, bm);
        currWB.update(rm, gm, bm, params.wb.equal);
    }

    NoiseCurve noiseLCurve;
    NoiseCurve noiseCCurve;
    Imagefloat *calclum = NULL ;
    params.dirpyrDenoise.getCurves(noiseLCurve, noiseCCurve);
    float autoNR = (float) settings->nrauto;//
    float autoNRmax = (float) settings->nrautomax;//
    int tilesize;
    int overlap;

    if(settings->leveldnti == 0) {
        tilesize = 1024;
        overlap = 128;
    }

    if(settings->leveldnti == 1) {
        tilesize = 768;
        overlap = 96;
    }

    //  const int tilesize = 768;
    //  const int overlap = 96;
    int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;
    ipf.Tile_calc (tilesize, overlap, 2, fw, fh, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);
    int nbtl = numtiles_W * numtiles_H;

    if((settings->leveldnautsimpl == 1 && params.dirpyrDenoise.Cmethod == "AUT") || (settings->leveldnautsimpl == 0 && params.dirpyrDenoise.C2method == "AUTO")) {
        nbtl = 9;
    }

    float *ch_M = new float [nbtl];//allocate memory
    float *max_r = new float [nbtl];
    float *max_b = new float [nbtl];
    float *min_b = new float [9];
    float *min_r = new float [9];
    float *lumL = new float [nbtl];
    float *chromC = new float [nbtl];
    float *ry = new float [nbtl];
    float *sk = new float [nbtl];
    float *pcsk = new float [nbtl];

    //  printf("expert=%d\n",settings->leveldnautsimpl);
    if(settings->leveldnautsimpl == 1 && params.dirpyrDenoise.Cmethod == "PON") {
        MyTime t1pone, t2pone;
        t1pone.set();
        int crW, crH;

        if(settings->leveldnv == 0) {
            crW = 100;
            crH = 100;
        }

        if(settings->leveldnv == 1) {
            crW = 250;
            crH = 250;
        }

        if(settings->leveldnv == 2) {
            crW = int(tileWskip / 2);
            crH = int(tileHskip / 2);
        }

        //  if(settings->leveldnv ==2) {crW=int(tileWskip/2);crH=int(1.15f*(tileWskip/2));}//adapted to scale of preview
        if(settings->leveldnv == 3) {
            crW = tileWskip - 10;
            crH = tileHskip - 10;
        }

        float lowdenoise = 1.f;
        int levaut = settings->leveldnaut;

        if(levaut == 1) { //Standard
            lowdenoise = 0.7f;
        }

        //  int crW=tileWskip-10;//crop noise width
        //  int crH=tileHskip-10;//crop noise height
//      Imagefloat *origCropPart;//init auto noise
//          origCropPart = new Imagefloat (crW, crH);//allocate memory
        if (params.dirpyrDenoise.enabled) {//evaluate Noise
            LUTf gamcurve(65536, 0);
            float gam, gamthresh, gamslope;
            ipf.RGB_denoise_infoGamCurve(params.dirpyrDenoise, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope);
            #pragma omp parallel
            {
                Imagefloat *origCropPart;//init auto noise
                origCropPart = new Imagefloat (crW, crH);//allocate memory
                Imagefloat *provicalc = new Imagefloat ((crW + 1) / 2, (crH + 1) / 2); //for denoise curves
                int skipP = 1;
                #pragma omp for schedule(dynamic) collapse(2) nowait

                for(int wcr = 0; wcr < numtiles_W; wcr++) {
                    for(int hcr = 0; hcr < numtiles_H; hcr++) {
                        int beg_tileW = wcr * tileWskip + tileWskip / 2.f - crW / 2.f;
                        int beg_tileH = hcr * tileHskip + tileHskip / 2.f - crH / 2.f;
                        PreviewProps ppP (beg_tileW , beg_tileH, crW, crH, skipP);
                        imgsrc->getImage (currWB, tr, origCropPart, ppP, params.toneCurve, params.icm, params.raw );

                        // we only need image reduced to 1/4 here
                        for(int ii = 0; ii < crH; ii += 2) {
                            for(int jj = 0; jj < crW; jj += 2) {
                                provicalc->r(ii >> 1, jj >> 1) = origCropPart->r(ii, jj);
                                provicalc->g(ii >> 1, jj >> 1) = origCropPart->g(ii, jj);
                                provicalc->b(ii >> 1, jj >> 1) = origCropPart->b(ii, jj);
                            }
                        }

                        imgsrc->convertColorSpace(provicalc, params.icm, currWB);//for denoise luminance curve
                        float maxr = 0.f;
                        float maxb = 0.f;
                        float pondcorrec = 1.0f;
                        float chaut, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, nresi, highresi, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc;
                        int Nb;
                        chaut = 0.f;
                        redaut = 0.f;
                        blueaut = 0.f;
                        maxredaut = 0.f;
                        maxblueaut = 0.f;
                        chromina = 0.f;
                        sigma = 0.f;
                        ipf.RGB_denoise_info(origCropPart, provicalc, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope, params.dirpyrDenoise, imgsrc->getDirPyrDenoiseExpComp(), chaut, Nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, nresi, highresi, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc);
                        float multip = 1.f;
                        float adjustr = 1.f;

                        if      (params.icm.working == "ProPhoto")   {
                            adjustr = 1.f;   //
                        } else if (params.icm.working == "Adobe RGB")  {
                            adjustr = 1.f / 1.3f;
                        } else if (params.icm.working == "sRGB")       {
                            adjustr = 1.f / 1.3f;
                        } else if (params.icm.working == "WideGamut")  {
                            adjustr = 1.f / 1.1f;
                        } else if (params.icm.working == "Rec2020")  {
                            adjustr = 1.f / 1.1f;
                        } else if (params.icm.working == "Beta RGB")   {
                            adjustr = 1.f / 1.2f;
                        } else if (params.icm.working == "BestRGB")    {
                            adjustr = 1.f / 1.2f;
                        } else if (params.icm.working == "BruceRGB")   {
                            adjustr = 1.f / 1.2f;
                        }

                        if(!imgsrc->isRAW()) {
                            multip = 2.f;    //take into account gamma for TIF / JPG approximate value...not good fot gamma=1
                        }

                        float maxmax = max(maxredaut, maxblueaut);
                        float delta;
                        int mode = 2;
                        int lissage = settings->leveldnliss;
                        ipf.calcautodn_info (chaut, delta, Nb, levaut, maxmax, lumema, chromina, mode, lissage, redyel, skinc, nsknc);

                        //    printf("PROCESS cha=%f red=%f bl=%f redM=%f bluM=%f chrom=%f sigm=%f lum=%f sigL=%f\n",chaut,redaut,blueaut, maxredaut, maxblueaut, chromina, sigma, lumema, sigma_L);
                        if(maxredaut > maxblueaut) {
                            maxr = (delta) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);

                            if(minblueaut <= minredaut  && minblueaut < chaut) {
                                maxb = (-chaut + minblueaut) / (autoNRmax * multip * adjustr * lowdenoise);
                            }
                        } else {
                            maxb = (delta) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);

                            if(minredaut <= minblueaut  && minredaut < chaut) {
                                maxr = (-chaut + minredaut) / (autoNRmax * multip * adjustr * lowdenoise);
                            }
                        }//maxb mxr - empirical evaluation red / blue

                        ch_M[hcr * numtiles_W + wcr] = pondcorrec * chaut / (autoNR * multip * adjustr * lowdenoise);
                        max_r[hcr * numtiles_W + wcr] = pondcorrec * maxr;
                        max_b[hcr * numtiles_W + wcr] = pondcorrec * maxb;
                        lumL[hcr * numtiles_W + wcr] = lumema;
                        chromC[hcr * numtiles_W + wcr] = chromina;
                        ry[hcr * numtiles_W + wcr] = redyel;
                        sk[hcr * numtiles_W + wcr] = skinc;
                        pcsk[hcr * numtiles_W + wcr] = nsknc;

                    }
                }

                delete provicalc;
                delete origCropPart;
            }

            int liss = settings->leveldnliss; //smooth result around mean

            if(liss == 2 || liss == 3) {
                // I smooth only mean and not delta (max)
                float nchm = 0.f;
                float koef = 0.4f; //between 0.1 to 0.9

                if(liss == 3) {
                    koef = 0.0f;    //quasi auto for mean Ch
                }

                for(int wcr = 0; wcr < numtiles_W; wcr++) {
                    for(int hcr = 0; hcr < numtiles_H; hcr++) {
                        nchm += ch_M[hcr * numtiles_W + wcr];
                    }
                }

                nchm /= (numtiles_H * numtiles_W);

                for(int wcr = 0; wcr < numtiles_W; wcr++) {
                    for(int hcr = 0; hcr < numtiles_H; hcr++) {
                        ch_M[hcr * numtiles_W + wcr] = nchm + (ch_M[hcr * numtiles_W + wcr] - nchm) * koef;
                    }
                }
            }

            if(liss == 3) { //same as auto but with much cells
                float MaxR = 0.f;
                float MaxB = 0.f;
                float MaxRMoy = 0.f;
                float MaxBMoy = 0.f;

                for(int k = 0; k < nbtl; k++) {
                    MaxBMoy += max_b[k];
                    MaxRMoy += max_r[k];

                    if(max_r[k] > MaxR) {
                        MaxR = max_r[k];
                    }

                    if(max_b[k] > MaxB) {
                        MaxB = max_b[k];
                    }

                }

                MaxBMoy /= nbtl;
                MaxRMoy /= nbtl;

                for(int k = 0; k < nbtl; k++) {
                    if(MaxR > MaxB) {
                        max_r[k] = MaxRMoy + (MaxR - MaxRMoy) * 0.66f; //#std Dev
                        //max_b[k]=MinB;
                        max_b[k] = MaxBMoy + (MaxB - MaxBMoy) * 0.66f;

                    } else {
                        max_b[k] = MaxBMoy + (MaxB - MaxBMoy) * 0.66f;
                        //max_r[k]=MinR;
                        max_r[k] = MaxRMoy + (MaxR - MaxRMoy) * 0.66f;

                    }
                }
            }

            if (settings->verbose) {
                t2pone.set();
                printf("Info denoise ponderated performed in %d usec:\n", t2pone.etime(t1pone));
            }

        }
    }


    if((settings->leveldnautsimpl == 1 && params.dirpyrDenoise.Cmethod == "AUT")  || (settings->leveldnautsimpl == 0 && params.dirpyrDenoise.C2method == "AUTO")) {
        MyTime t1aue, t2aue;
        t1aue.set();
        int crW, crH;

        if(settings->leveldnv == 0) {
            crW = 100;
            crH = 100;
        }

        if(settings->leveldnv == 1) {
            crW = 250;
            crH = 250;
        }

        if(settings->leveldnv == 2) {
            crW = int(tileWskip / 2);
            crH = int(tileHskip / 2);
        }

        //  if(settings->leveldnv ==2) {crW=int(tileWskip/2);crH=int(1.15f*(tileWskip/2));}//adapted to scale of preview
        if(settings->leveldnv == 3) {
            crW = tileWskip - 10;
            crH = tileHskip - 10;
        }

        float lowdenoise = 1.f;
        int levaut = settings->leveldnaut;

        if(levaut == 1) { //Standard
            lowdenoise = 0.7f;
        }

        if (params.dirpyrDenoise.enabled) {//evaluate Noise
            LUTf gamcurve(65536, 0);
            float gam, gamthresh, gamslope;
            ipf.RGB_denoise_infoGamCurve(params.dirpyrDenoise, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope);
            int Nb[9];
            int  coordW[3];//coordonate of part of image to mesure noise
            int  coordH[3];
            int begW = 50;
            int begH = 50;
            coordW[0] = begW;
            coordW[1] = fw / 2 - crW / 2;
            coordW[2] = fw - crW - begW;
            coordH[0] = begH;
            coordH[1] = fh / 2 - crH / 2;
            coordH[2] = fh - crH - begH;
            #pragma omp parallel
            {
                Imagefloat *origCropPart;//init auto noise
                origCropPart = new Imagefloat (crW, crH);//allocate memory
                Imagefloat *provicalc = new Imagefloat ((crW + 1) / 2, (crH + 1) / 2); //for denoise curves

                #pragma omp for schedule(dynamic) collapse(2) nowait

                for(int wcr = 0; wcr <= 2; wcr++) {
                    for(int hcr = 0; hcr <= 2; hcr++) {
                        PreviewProps ppP (coordW[wcr] , coordH[hcr], crW, crH, 1);
                        imgsrc->getImage (currWB, tr, origCropPart, ppP, params.toneCurve, params.icm, params.raw);

                        // we only need image reduced to 1/4 here
                        for(int ii = 0; ii < crH; ii += 2) {
                            for(int jj = 0; jj < crW; jj += 2) {
                                provicalc->r(ii >> 1, jj >> 1) = origCropPart->r(ii, jj);
                                provicalc->g(ii >> 1, jj >> 1) = origCropPart->g(ii, jj);
                                provicalc->b(ii >> 1, jj >> 1) = origCropPart->b(ii, jj);
                            }
                        }

                        imgsrc->convertColorSpace(provicalc, params.icm, currWB);//for denoise luminance curve
                        int nb = 0;
                        float chaut = 0.f, redaut = 0.f, blueaut = 0.f, maxredaut = 0.f, maxblueaut = 0.f, minredaut = 0.f, minblueaut = 0.f, nresi = 0.f, highresi = 0.f, chromina = 0.f, sigma = 0.f, lumema = 0.f, sigma_L = 0.f, redyel = 0.f, skinc = 0.f, nsknc = 0.f;
                        ipf.RGB_denoise_info(origCropPart, provicalc, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope,  params.dirpyrDenoise, imgsrc->getDirPyrDenoiseExpComp(), chaut, nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, nresi, highresi, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc);
                        Nb[hcr * 3 + wcr] = nb;
                        ch_M[hcr * 3 + wcr] = chaut;
                        max_r[hcr * 3 + wcr] = maxredaut;
                        max_b[hcr * 3 + wcr] = maxblueaut;
                        min_r[hcr * 3 + wcr] = minredaut;
                        min_b[hcr * 3 + wcr] = minblueaut;
                        lumL[hcr * 3 + wcr] = lumema;
                        chromC[hcr * 3 + wcr] = chromina;
                        ry[hcr * 3 + wcr] = redyel;
                        sk[hcr * 3 + wcr] = skinc;
                        pcsk[hcr * 3 + wcr] = nsknc;
                    }
                }

                delete provicalc;
                delete origCropPart;
            }
            float chM = 0.f;
            float MaxR = 0.f;
            float MaxB = 0.f;
            float MinR = 100000000.f;
            float MinB = 100000000.f;
            float maxr = 0.f;
            float maxb = 0.f;
            float multip = 1.f;
            float adjustr = 1.f;
            float Max_R[9] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
            float Max_B[9] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
            float Min_R[9];
            float Min_B[9];
            float MaxRMoy = 0.f;
            float MaxBMoy = 0.f;
            float MinRMoy = 0.f;
            float MinBMoy = 0.f;

            if      (params.icm.working == "ProPhoto")   {
                adjustr = 1.f;
            } else if (params.icm.working == "Adobe RGB")  {
                adjustr = 1.f / 1.3f;
            } else if (params.icm.working == "sRGB")       {
                adjustr = 1.f / 1.3f;
            } else if (params.icm.working == "WideGamut")  {
                adjustr = 1.f / 1.1f;
            } else if (params.icm.working == "Rec2020")  {
                adjustr = 1.f / 1.1f;
            } else if (params.icm.working == "Beta RGB")   {
                adjustr = 1.f / 1.2f;
            } else if (params.icm.working == "BestRGB")    {
                adjustr = 1.f / 1.2f;
            } else if (params.icm.working == "BruceRGB")   {
                adjustr = 1.f / 1.2f;
            }

            if(!imgsrc->isRAW()) {
                multip = 2.f;    //take into account gamma for TIF / JPG approximate value...not good fot gamma=1
            }

            float delta[9];
            int mode = 1;
            int lissage = settings->leveldnliss;

            for(int k = 0; k < 9; k++) {
                float maxmax = max(max_r[k], max_b[k]);
                ipf.calcautodn_info (ch_M[k], delta[k], Nb[k], levaut, maxmax, lumL[k], chromC[k], mode, lissage, ry[k], sk[k], pcsk[k] );
                //  printf("ch_M=%f delta=%f\n",ch_M[k], delta[k]);
            }

            for(int k = 0; k < 9; k++) {
                if(max_r[k] > max_b[k]) {
                    //printf("R delta=%f  koef=%f\n",delta[k],autoNRmax*multip*adjustr*lowdenoise);
                    Max_R[k] = (delta[k]) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);
                    Min_B[k] = -(ch_M[k] - min_b[k]) / (autoNRmax * multip * adjustr * lowdenoise);
                    Max_B[k] = 0.f;
                    Min_R[k] = 0.f;
                } else {
                    //printf("B delta=%f  koef=%f\n",delta[k],autoNRmax*multip*adjustr*lowdenoise);
                    Max_B[k] = (delta[k]) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);
                    Min_R[k] = - (ch_M[k] - min_r[k])   / (autoNRmax * multip * adjustr * lowdenoise);
                    Min_B[k] = 0.f;
                    Max_R[k] = 0.f;
                }
            }

            for(int k = 0; k < 9; k++) {
                //  printf("ch_M= %f Max_R=%f Max_B=%f min_r=%f min_b=%f\n",ch_M[k],Max_R[k], Max_B[k],Min_R[k], Min_B[k]);
                chM += ch_M[k];
                MaxBMoy += Max_B[k];
                MaxRMoy += Max_R[k];
                MinRMoy += Min_R[k];
                MinBMoy += Min_B[k];

                if(Max_R[k] > MaxR) {
                    MaxR = Max_R[k];
                }

                if(Max_B[k] > MaxB) {
                    MaxB = Max_B[k];
                }

                if(Min_R[k] < MinR) {
                    MinR = Min_R[k];
                }

                if(Min_B[k] < MinB) {
                    MinB = Min_B[k];
                }

            }

            chM /= 9;
            MaxBMoy /= 9;
            MaxRMoy /= 9;
            MinBMoy /= 9;
            MinRMoy /= 9;

            if(MaxR > MaxB) {
                maxr = MaxRMoy + (MaxR - MaxRMoy) * 0.66f; //#std Dev
                //  maxb=MinB;
                maxb = MinBMoy + (MinB - MinBMoy) * 0.66f;

            } else {
                maxb = MaxBMoy + (MaxB - MaxBMoy) * 0.66f;
                //  maxr=MinR;
                maxr = MinRMoy + (MinR - MinRMoy) * 0.66f;

            }

//              printf("SIMPL cha=%f red=%f bl=%f \n",chM,maxr,maxb);

            params.dirpyrDenoise.chroma = chM / (autoNR * multip * adjustr);
            params.dirpyrDenoise.redchro = maxr;
            params.dirpyrDenoise.bluechro = maxb;
        }

        if (settings->verbose) {
            t2aue.set();
            printf("Info denoise auto performed in %d usec:\n", t2aue.etime(t1aue));
        }

        //end evaluate noise
    }





    Imagefloat* baseImg = new Imagefloat (fw, fh);
    imgsrc->getImage (currWB, tr, baseImg, pp, params.toneCurve, params.icm, params.raw);

    if (pl) {
        pl->setProgress (0.50);
    }

//  LUTf Noisecurve (65536,0);
//!!!// auto exposure!!!
    double expcomp = params.toneCurve.expcomp;
    int    bright = params.toneCurve.brightness;
    int    contr = params.toneCurve.contrast;
    int    black = params.toneCurve.black;
    int    hlcompr = params.toneCurve.hlcompr;
    int    hlcomprthresh = params.toneCurve.hlcomprthresh;

    if (params.toneCurve.autoexp) {
        LUTu aehist;
        int aehistcompr;
        imgsrc->getAutoExpHistogram (aehist, aehistcompr);
        ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, expcomp, bright, contr, black, hlcompr, hlcomprthresh);
    }

    // at this stage, we can flush the raw data to free up quite an important amount of memory
    // commented out because it makes the application crash when batch processing...
    // TODO: find a better place to flush rawData and rawRGB
    if(flush) {
        imgsrc->flushRawData();
        imgsrc->flushRGB();
    }

    // perform luma/chroma denoise
//  CieImage *cieView;
//  NoisCurve noiseLCurve;
//    bool lldenoiseutili=false;
//  Imagefloat *calclum ;
//    params.dirpyrDenoise.getCurves(noiseLCurve, lldenoiseutili);
//  if (params.dirpyrDenoise.enabled  && lldenoiseutili) {

    DirPyrDenoiseParams denoiseParams = params.dirpyrDenoise;   // make a copy because we cheat here

    if(denoiseParams.Lmethod == "CUR") {
        if(noiseLCurve) {
            denoiseParams.luma = 0.5f;
        } else {
            denoiseParams.luma = 0.0f;
        }
    } else if(denoiseParams.Lmethod == "SLI") {
        noiseLCurve.Reset();
    }

    if (denoiseParams.enabled  && (noiseLCurve || noiseCCurve )) {
        // we only need image reduced to 1/4 here
        calclum = new Imagefloat ((fw + 1) / 2, (fh + 1) / 2); //for luminance denoise curve
        #pragma omp parallel for

        for(int ii = 0; ii < fh; ii += 2) {
            for(int jj = 0; jj < fw; jj += 2) {
                calclum->r(ii >> 1, jj >> 1) = baseImg->r(ii, jj);
                calclum->g(ii >> 1, jj >> 1) = baseImg->g(ii, jj);
                calclum->b(ii >> 1, jj >> 1) = baseImg->b(ii, jj);
            }
        }

        imgsrc->convertColorSpace(calclum, params.icm, currWB);
    }

    if (denoiseParams.enabled) {
        // CurveFactory::denoiseLL(lldenoiseutili, denoiseParams.lcurve, Noisecurve,1);
        //denoiseParams.getCurves(noiseLCurve);
//      ipf.RGB_denoise(baseImg, baseImg, calclum, imgsrc->isRAW(), denoiseParams, params.defringe, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, lldenoiseutili);
        float chaut, redaut, blueaut, maxredaut, maxblueaut, nresi, highresi;
        int kall = 2;
        ipf.RGB_denoise(kall, baseImg, baseImg, calclum, ch_M, max_r, max_b, imgsrc->isRAW(), denoiseParams, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, noiseCCurve, chaut, redaut, blueaut, maxredaut, maxblueaut, nresi, highresi);

    }

//  delete calclum;
    delete [] ch_M;
    delete [] max_r;
    delete [] max_b;
    delete [] min_r;
    delete [] min_b;
    delete [] lumL;
    delete [] chromC;
    delete [] ry;
    delete [] sk;
    delete [] pcsk;

    imgsrc->convertColorSpace(baseImg, params.icm, currWB);

    // perform first analysis
    LUTu hist16 (65536);

    ipf.firstAnalysis (baseImg, params, hist16);

    // perform transform (excepted resizing)
    if (ipf.needsTransform()) {
        Imagefloat* trImg = new Imagefloat (fw, fh);
        ipf.transform (baseImg, trImg, 0, 0, 0, 0, fw, fh, fw, fh, imgsrc->getMetaData()->getFocalLen(), imgsrc->getMetaData()->getFocalLen35mm(),
                       imgsrc->getMetaData()->getFocusDist(), imgsrc->getRotateDegree(), true);
        delete baseImg;
        baseImg = trImg;
    }


    if (params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled) {
        const int W = baseImg->getWidth();
        const int H = baseImg->getHeight();
        LabImage labcbdl(W, H);
        ipf.rgb2lab(*baseImg, labcbdl, params.icm.working);
        ipf.dirpyrequalizer (&labcbdl, 1);
        ipf.lab2rgb(labcbdl, *baseImg, params.icm.working);
    }

    // update blurmap
    SHMap* shmap = NULL;

    if (params.sh.enabled) {
        shmap = new SHMap (fw, fh, true);
        double radius = sqrt (double(fw * fw + fh * fh)) / 2.0;
        double shradius = params.sh.radius;

        if (!params.sh.hq) {
            shradius *= radius / 1800.0;
        }

        shmap->update (baseImg, shradius, ipf.lumimul, params.sh.hq, 1);
    }

    if (params.spot.enabled && !params.spot.entries.empty()) {
        ipf.removeSpots(baseImg, params.spot.entries, pp);
    }

    // RGB processing

    LUTf curve1 (65536);
    LUTf curve2 (65536);
    LUTf curve (65536, 0);
    LUTf satcurve (65536, 0);
    LUTf lhskcurve (65536, 0);
    LUTf lumacurve(32770, 0); // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
    LUTf clcurve (65536, 0);
    LUTf clToningcurve;
    LUTf cl2Toningcurve;
    LUTf wavclCurve (65536, 0);

    LUTf rCurve;
    LUTf gCurve;
    LUTf bCurve;
    LUTu dummy;

    ToneCurve customToneCurve1, customToneCurve2;
    ColorGradientCurve ctColorCurve;
    OpacityCurve ctOpacityCurve;
    ColorAppearance customColCurve1, customColCurve2, customColCurve3 ;
    ToneCurve customToneCurvebw1;
    ToneCurve customToneCurvebw2;
    //if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;

    CurveFactory::complexCurve (expcomp, black / 65535.0, hlcompr, hlcomprthresh, params.toneCurve.shcompr, bright, contr,
                                params.toneCurve.curveMode, params.toneCurve.curve, params.toneCurve.curveMode2, params.toneCurve.curve2,
                                hist16, curve1, curve2, curve, dummy, customToneCurve1, customToneCurve2 );

    CurveFactory::RGBCurve (params.rgbCurves.rcurve, rCurve, 1);
    CurveFactory::RGBCurve (params.rgbCurves.gcurve, gCurve, 1);
    CurveFactory::RGBCurve (params.rgbCurves.bcurve, bCurve, 1);

    bool opautili = false;

    if(params.colorToning.enabled) {
        TMatrix wprof = iccStore->workingSpaceMatrix (params.icm.working);
        double wp[3][3] = {
            {wprof[0][0], wprof[0][1], wprof[0][2]},
            {wprof[1][0], wprof[1][1], wprof[1][2]},
            {wprof[2][0], wprof[2][1], wprof[2][2]}
        };
        TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params.icm.working);
        double wip[3][3] = {
            {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
            {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
            {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
        };
        params.colorToning.getCurves(ctColorCurve, ctOpacityCurve, wp, wip, opautili);
        clToningcurve (65536, 0);
        CurveFactory::curveToning(params.colorToning.clcurve, clToningcurve, 1);
        cl2Toningcurve (65536, 0);
        CurveFactory::curveToning(params.colorToning.cl2curve, cl2Toningcurve, 1);
    }

    LabImage* labView = new LabImage (fw, fh);

    if(params.blackwhite.enabled) {
        CurveFactory::curveBW (params.blackwhite.beforeCurve, params.blackwhite.afterCurve, hist16, dummy, customToneCurvebw1, customToneCurvebw2, 1);
    }

    double rrm, ggm, bbm;
    float autor, autog, autob;
    float satLimit = float(params.colorToning.satProtectionThreshold) / 100.f * 0.7f + 0.3f;
    float satLimitOpacity = 1.f - (float(params.colorToning.saturatedOpacity) / 100.f);

    if(params.colorToning.enabled  && params.colorToning.autosat) { //for colortoning evaluation of saturation settings
        float moyS = 0.f;
        float eqty = 0.f;
        ipf.moyeqt (baseImg, moyS, eqty);//return image : mean saturation and standard dev of saturation
        float satp = ((moyS + 1.5f * eqty) - 0.3f) / 0.7f; //1.5 sigma ==> 93% pixels with high saturation -0.3 / 0.7 convert to Hombre scale

        if(satp >= 0.92f) {
            satp = 0.92f;    //avoid values too high (out of gamut)
        }

        if(satp <= 0.15f) {
            satp = 0.15f;    //avoid too low values
        }

        satLimit = 100.f * satp;

        satLimitOpacity = 100.f * (moyS - 0.85f * eqty); //-0.85 sigma==>20% pixels with low saturation
    }

    autor = -9000.f; // This will ask to compute the "auto" values for the B&W tool (have to be inferior to -5000)
    DCPProfile::ApplyState as;
    DCPProfile *dcpProf = imgsrc->getDCP(params.icm, currWB, as);

    ipf.rgbProc (baseImg, labView, NULL, curve1, curve2, curve, shmap, params.toneCurve.saturation, rCurve, gCurve, bCurve, satLimit , satLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2, customToneCurvebw1, customToneCurvebw2, rrm, ggm, bbm, autor, autog, autob, expcomp, hlcompr, hlcomprthresh, dcpProf, as);

    if (settings->verbose) {
        printf("Output image / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", autor, autog, autob);
    }

    // if clut was used and size of clut cache == 1 we free the memory used by the clutstore (default clut cache size = 1 for 32 bit OS)
    if ( params.filmSimulation.enabled && !params.filmSimulation.clutFilename.empty() && options.clutCacheSize == 1) {
        CLUTStore::getInstance().clearCache();
    }

    // freeing up some memory
    customToneCurve1.Reset();
    customToneCurve2.Reset();
    ctColorCurve.Reset();
    ctOpacityCurve.Reset();
    noiseLCurve.Reset();
    noiseCCurve.Reset();
    customToneCurvebw1.Reset();
    customToneCurvebw2.Reset();

    // Freeing baseImg because not used anymore
    delete baseImg;
    baseImg = NULL;

    if (shmap) {
        delete shmap;
    }

    shmap = NULL;

    if (pl) {
        pl->setProgress (0.55);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // start tile processing...???


    if(params.labCurve.contrast != 0) { //only use hist16 for contrast
        hist16.clear();

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            LUTu hist16thr (hist16.getSize());  // one temporary lookup table per thread
            hist16thr.clear();
#ifdef _OPENMP
            #pragma omp for schedule(static) nowait
#endif

            for (int i = 0; i < fh; i++)
                for (int j = 0; j < fw; j++) {
                    hist16thr[(int)((labView->L[i][j]))]++;
                }

            #pragma omp critical
            {
                hist16 += hist16thr;
            }
        }
    }

    bool utili;
    CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, hist16, lumacurve, dummy, 1, utili);

    bool clcutili;
    CurveFactory::curveCL(clcutili, params.labCurve.clcurve, clcurve, 1);

    bool autili, butili, ccutili, cclutili;
    CurveFactory::complexsgnCurve (autili, butili, ccutili, cclutili, params.labCurve.acurve, params.labCurve.bcurve, params.labCurve.cccurve,
                                   params.labCurve.lccurve, curve1, curve2, satcurve, lhskcurve, 1);

    ipf.chromiLuminanceCurve (NULL, 1, labView, labView, curve1, curve2, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, dummy, dummy);

    if((params.colorappearance.enabled && !params.colorappearance.tonecie) || (!params.colorappearance.enabled)) {
        ipf.EPDToneMap(labView, 5, 1);
    }


    ipf.vibrance(labView);

    if((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
        ipf.impulsedenoise (labView);
    }

    // for all treatments Defringe, Sharpening, Contrast detail ,Microcontrast they are activated if "CIECAM" function are disabled

    if((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
        ipf.defringe (labView);
    }

    if (params.sharpenEdge.enabled) {
        ipf.MLsharpen(labView);
    }

    if (params.sharpenMicro.enabled) {
        if((params.colorappearance.enabled && !settings->autocielab) ||  (!params.colorappearance.enabled)) {
            ipf.MLmicrocontrast (labView);    //!params.colorappearance.sharpcie
        }
    }

    if(((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) && params.sharpening.enabled) {

        float **buffer = new float*[fh];

        for (int i = 0; i < fh; i++) {
            buffer[i] = new float[fw];
        }

        ipf.sharpening (labView, (float**)buffer, params.sharpening);

        for (int i = 0; i < fh; i++) {
            delete [] buffer[i];
        }

        delete [] buffer;
    }

    WaveletParams WaveParams = params.wavelet;
    WavCurve wavCLVCurve;
    WavOpacityCurveRG waOpacityCurveRG;
    WavOpacityCurveBY waOpacityCurveBY;
    WavOpacityCurveW waOpacityCurveW;
    WavOpacityCurveWL waOpacityCurveWL;

    params.wavelet.getCurves(wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL );


    // directional pyramid wavelet
    if(params.dirpyrequalizer.cbdlMethod == "aft") {
        if((params.colorappearance.enabled && !settings->autocielab)  || !params.colorappearance.enabled) {
            ipf.dirpyrequalizer (labView, 1);    //TODO: this is the luminance tonecurve, not the RGB one
        }
    }

    int kall = 2;
    bool wavcontlutili = false;

    CurveFactory::curveWavContL(wavcontlutili, params.wavelet.wavclCurve, wavclCurve,/* hist16C, dummy,*/ 1);

    if((params.wavelet.enabled)) {
        ipf.ip_wavelet(labView, labView, kall, WaveParams, wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW,  waOpacityCurveWL, wavclCurve, wavcontlutili, 1);
    }

    wavCLVCurve.Reset();

    //Colorappearance and tone-mapping associated

    int f_w = 1, f_h = 1;
    int begh = 0, endh = fh;

    if(params.colorappearance.tonecie || params.colorappearance.enabled) {
        f_w = fw;
        f_h = fh;
    }

    CieImage *cieView = new CieImage (f_w, (f_h));
    begh = 0;
    endh = fh;
    CurveFactory::curveLightBrightColor (
        params.colorappearance.curve,
        params.colorappearance.curve2,
        params.colorappearance.curve3,
        hist16, dummy,
        dummy, dummy,
        customColCurve1,
        customColCurve2,
        customColCurve3,
        1);

    if(params.colorappearance.enabled) {
        double adap;
        float fnum = imgsrc->getMetaData()->getFNumber  ();// F number
        float fiso = imgsrc->getMetaData()->getISOSpeed () ;// ISO
        float fspeed = imgsrc->getMetaData()->getShutterSpeed () ;//speed
        float fcomp = imgsrc->getMetaData()->getExpComp  ();//compensation + -

        if(fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) {
            adap = 2000.;
        }//if no exif data or wrong
        else {
            float E_V = fcomp + log2 ((fnum * fnum) / fspeed / (fiso / 100.f));
            E_V += params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
            E_V += log2(params.raw.expos);// exposure raw white point ; log2 ==> linear to EV
            adap = powf(2.f, E_V - 3.f); //cd / m2
        }

        LUTf CAMBrightCurveJ;
        LUTf CAMBrightCurveQ;
        float CAMMean = NAN;

        if (params.sharpening.enabled) {
            float d;
            double dd;

            int sk = 1;

            if(settings->ciecamfloat) {
                ipf.ciecam_02float (cieView, float(adap), begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, d, sk, 1);
            } else {
                ipf.ciecam_02 (cieView, adap, begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, dd, 1, 1);
            }
        } else {
            float d;

            double dd;
            int sk = 1;

            if(settings->ciecamfloat) {
                ipf.ciecam_02float (cieView, float(adap), begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, d, sk, 1);
            } else {
                ipf.ciecam_02 (cieView, adap, begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, dd, 1, 1);
            }
        }
    }

    delete cieView;
    cieView = NULL;




    // end tile processing...???
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (pl) {
        pl->setProgress (0.60);
    }

    int imw, imh;
    double tmpScale = ipf.resizeScale(&params, fw, fh, imw, imh);
    bool labResize = params.resize.enabled && params.resize.method != "Nearest" && tmpScale != 1.0;
    LabImage *tmplab;

    // crop and convert to rgb16
    int cx = 0, cy = 0, cw = labView->W, ch = labView->H;

    if (params.crop.enabled) {
        cx = params.crop.x;
        cy = params.crop.y;
        cw = params.crop.w;
        ch = params.crop.h;

        if(labResize) { // crop lab data
            tmplab = new LabImage(cw, ch);

            for(int row = 0; row < ch; row++) {
                for(int col = 0; col < cw; col++) {
                    tmplab->L[row][col] = labView->L[row + cy][col + cx];
                    tmplab->a[row][col] = labView->a[row + cy][col + cx];
                    tmplab->b[row][col] = labView->b[row + cy][col + cx];
                }
            }

            delete labView;
            labView = tmplab;
            cx = 0;
            cy = 0;
        }
    }

    if (labResize) { // resize lab data
        // resize image
        tmplab = new LabImage(imw, imh);
        ipf.Lanczos (labView, tmplab, tmpScale);
        delete labView;
        labView = tmplab;
        cw = labView->W;
        ch = labView->H;

        if(params.prsharpening.enabled) {
            for(int i = 0; i < ch; i++)
                for(int j = 0; j < cw; j++) {
                    labView->L[i][j] = labView->L[i][j] < 0.f ? 0.f : labView->L[i][j];
                }

            float **buffer = new float*[ch];

            for (int i = 0; i < ch; i++) {
                buffer[i] = new float[cw];
            }

            ipf.sharpening (labView, (float**)buffer, params.prsharpening);

            for (int i = 0; i < ch; i++) {
                delete [] buffer[i];
            }

            delete [] buffer;
        }
    }

    Image16* readyImg = NULL;
    cmsHPROFILE jprof = NULL;
    bool customGamma = false;
    bool useLCMS = false;

    if(params.icm.gamma != "default" || params.icm.freegamma) { // if select gamma output between BT709, sRGB, linear, low, high, 2.2 , 1.8
        cmsMLU *DescriptionMLU, *CopyrightMLU, *DmndMLU, *DmddMLU;// for modification TAG

        cmsToneCurve* GammaTRC[3] = { NULL, NULL, NULL };
        cmsFloat64Number Parameters[7];
        double ga0, ga1, ga2, ga3, ga4, ga5, ga6;
        //  if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;
        readyImg = ipf.lab2rgb16b (labView, cx, cy, cw, ch, params.icm.output, params.icm.outputIntent, params.icm.working, params.icm.gamma, params.icm.freegamma, params.icm.gampos, params.icm.slpos, ga0, ga1, ga2, ga3, ga4, ga5, ga6, params.blackwhite.enabled );
        customGamma = true;

        //or selected Free gamma
        useLCMS = false;
        bool pro = false;
        Glib::ustring chpro, outProfile;
        bool present_space[10] = {false, false, false, false, false, false, false, false, false, false};
        std::vector<Glib::ustring> opnames = iccStore->getProfiles ();

        //test if files are in system
        for (int j = 0; j < 10; j++) {
            // one can modify "option" [Color Management] to adapt the profile's name if they are different for windows, MacOS, Linux ??
            // some of them are actually provided by RT, thanks to Jacques Desmis
            if     (j == 0) {
                chpro = options.rtSettings.prophoto;
            } else if(j == 1) {
                chpro = options.rtSettings.adobe;
            } else if(j == 2) {
                chpro = options.rtSettings.widegamut;
            } else if(j == 3) {
                chpro = options.rtSettings.beta;
            } else if(j == 4) {
                chpro = options.rtSettings.best;
            } else if(j == 5) {
                chpro = options.rtSettings.bruce;
            } else if(j == 6) {
                chpro = options.rtSettings.srgb;
            } else if(j == 7) {
                chpro = options.rtSettings.srgb10;    //gamma 1.0
            } else if(j == 8) {
                chpro = options.rtSettings.prophoto10;    //gamma 1.0
            } else if(j == 9) {
                chpro = options.rtSettings.rec2020;
            }

            for (unsigned int i = 0; i < opnames.size(); i++) {
                if(chpro.compare(opnames[i]) == 0) {
                    present_space[j] = true;
                }
            }

            if (!present_space[j] && settings->verbose) {
                printf("Missing file: %s\n", chpro.c_str());
            }
        }

        if (params.icm.freegamma && params.icm.gampos < 1.35) {
            pro = true;    //select profil with gammaTRC modified :
        } else if (params.icm.gamma == "linear_g1.0" || (params.icm.gamma == "High_g1.3_s3.35")) {
            pro = true;    //pro=0  RT_sRGB || Prophoto
        }

        // Check that output profiles exist, otherwise use LCMS2
        // Use the icc/icm profiles associated to possible working profiles, set in "options"
        if      (params.icm.working == "ProPhoto"  && present_space[0] && !pro) {
            outProfile = options.rtSettings.prophoto;
        } else if (params.icm.working == "Adobe RGB" && present_space[1]        ) {
            outProfile = options.rtSettings.adobe;
        } else if (params.icm.working == "WideGamut" && present_space[2]        ) {
            outProfile = options.rtSettings.widegamut;
        } else if (params.icm.working == "Beta RGB"  && present_space[3]        ) {
            outProfile = options.rtSettings.beta;
        } else if (params.icm.working == "BestRGB"   && present_space[4]        ) {
            outProfile = options.rtSettings.best;
        } else if (params.icm.working == "BruceRGB"  && present_space[5]        ) {
            outProfile = options.rtSettings.bruce;
        } else if (params.icm.working == "sRGB"      && present_space[6] && !pro) {
            outProfile = options.rtSettings.srgb;
        } else if (params.icm.working == "sRGB"      && present_space[7] &&  pro) {
            outProfile = options.rtSettings.srgb10;
        } else if (params.icm.working == "ProPhoto"  && present_space[8] &&  pro) {
            outProfile = options.rtSettings.prophoto10;
        } else if (params.icm.working == "Rec2020"  && present_space[9]) {
            outProfile = options.rtSettings.rec2020;
        } else {
            // Should not occurs
            if (settings->verbose) {
                printf("\"%s\": unknown working profile! - use LCMS2 substitution\n", params.icm.working.c_str() );
            }

            useLCMS = true;
        }

        //begin adaptation rTRC gTRC bTRC
        //"jprof" profile has the same characteristics than RGB values, but TRC are adapted... for applying profile
        if (!useLCMS) {
            if (settings->verbose) {
                printf("Output Gamma - profile: \"%s\"\n", outProfile.c_str()  );    //c_str()
            }

            jprof = iccStore->getProfile(outProfile); //get output profile

            if (jprof == NULL) {
                useLCMS = true;

                if (settings->verbose) {
                    printf("\"%s\" ICC output profile not found!\n", outProfile.c_str());
                }
            } else {
                Parameters[0] = ga0;
                Parameters[1] = ga1;
                Parameters[2] = ga2;
                Parameters[3] = ga3;
                Parameters[4] = ga4;
                Parameters[5] = ga5;
                Parameters[6] = ga6;
                // 7 parameters for smoother curves
                //change desc Tag , to "free gamma", or "BT709", etc.
                cmsContext ContextID = cmsGetProfileContextID(jprof);//modification TAG
                DescriptionMLU  = cmsMLUalloc(ContextID, 1);
                CopyrightMLU    = cmsMLUalloc(ContextID, 1);//for ICC
                DmndMLU = cmsMLUalloc(ContextID, 1); //for ICC
                DmddMLU = cmsMLUalloc(ContextID, 1); // for ICC


                // instruction with //ICC are used for generate icc profile
                if (DescriptionMLU == NULL) {
                    printf("Description error\n");
                }

                cmsMLUsetWide(CopyrightMLU, "en", "US", L"General Public License - AdobeRGB compatible")    ;//adapt to profil
                cmsMLUsetWide(DmndMLU,      "en", "US", L"RawTherapee") ;
                cmsMLUsetWide(DmddMLU,      "en", "US", L"RTMedium")    ;   //adapt to profil

                //display Tag desc with : selection of gamma and Primaries
                if (!params.icm.freegamma) {
                    std::wstring gammaStr;

                    if(params.icm.gamma == "High_g1.3_s3.35") {
                        gammaStr = std::wstring(L"GammaTRC: High g=1.3 s=3.35");
                    } else if (params.icm.gamma == "Low_g2.6_s6.9") {
                        gammaStr = std::wstring(L"GammaTRC: Low g=2.6 s=6.9");
                    } else if (params.icm.gamma == "sRGB_g2.4_s12.92") {
                        gammaStr = std::wstring(L"GammaTRC: sRGB g=2.4 s=12.92");
                    } else if (params.icm.gamma == "BT709_g2.2_s4.5") {
                        gammaStr = std::wstring(L"GammaTRC: BT709 g=2.2 s=4.5");
                    } else if (params.icm.gamma == "linear_g1.0") {
                        gammaStr = std::wstring(L"GammaTRC: Linear g=1.0");
                    } else if (params.icm.gamma == "standard_g2.2") {
                        gammaStr = std::wstring(L"GammaTRC: g=2.2");
                    } else if (params.icm.gamma == "standard_g1.8") {
                        gammaStr = std::wstring(L"GammaTRC: g=1.8");
                    }

                    cmsMLUsetWide(DescriptionMLU,  "en", "US", gammaStr.c_str());

                    //for elaboration ICC profiles
                    //  else if (params.icm.gamma== "sRGB_g2.4_s12.92" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_Medium gamma sRGB(AdobeRGB compatible)");
                    //  else if (params.icm.gamma== "BT709_g2.2_s4.5" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_sRGB gamma BT709(IEC61966 equivalent)");
                    //  else if (params.icm.gamma== "sRGB_g2.4_s12.92" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_sRGB gamma sRGB(IEC61966 equivalent)");
                    //  else if (params.icm.gamma== "linear_g1.0" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_sRGB gamma Linear1.0(IEC61966 equivalent)");
                    //else if (params.icm.gamma==   "BT709_g2.2_s4.5" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_Large gamma BT709(Prophoto compatible)");
                    //  else if (params.icm.gamma== "sRGB_g2.4_s12.92" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_Large gamma sRGB(Prophoto compatible)");
                    //  else if (params.icm.gamma== "linear_g1.0" && !params.icm.freegamma) cmsMLUsetWide(DescriptionMLU,  "en", "US", L"RT_Large gamma Linear1.0(Prophoto compatible)");
                } else {
                    // create description with gamma + slope + primaries
                    std::wostringstream gammaWs;
                    gammaWs.precision(2);
                    gammaWs << "Manual GammaTRC: g=" << (float)params.icm.gampos << " s=" << (float)params.icm.slpos;
                    cmsMLUsetWide(DescriptionMLU,  "en", "US", gammaWs.str().c_str());
                }

                cmsWriteTag(jprof, cmsSigProfileDescriptionTag,  DescriptionMLU);//desc changed
                //  cmsWriteTag(jprof, cmsSigCopyrightTag,           CopyrightMLU);
                //  cmsWriteTag(jprof, cmsSigDeviceMfgDescTag, DmndMLU);
                //  cmsWriteTag(jprof, cmsSigDeviceModelDescTag, DmddMLU);

                // Calculate output profile's rTRC bTRC gTRC
                GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(NULL, 5, Parameters);
                cmsWriteTag(jprof, cmsSigGreenTRCTag, (void*)GammaTRC[1] );
                cmsWriteTag(jprof, cmsSigRedTRCTag,   (void*)GammaTRC[0] );
                cmsWriteTag(jprof, cmsSigBlueTRCTag,  (void*)GammaTRC[2] );
                //for generation ICC profiles : here Prophoto ==> Large
                //  if(params.icm.gamma==   "BT709_g2.2_s4.5") cmsSaveProfileToFile(jprof, "RT_sRGB_gBT709.icm");
                //  else if (params.icm.gamma== "sRGB_g2.4_s12.92") cmsSaveProfileToFile(jprof, "RT_Medium_gsRGB.icc");
                //  else if (params.icm.gamma== "linear_g1.0") cmsSaveProfileToFile(jprof, "RT_Large_g10.icc");


            }
        }

        if (GammaTRC[0]) {
            cmsFreeToneCurve(GammaTRC[0]);
        }
    } else {
        // if Default gamma mode: we use the profile selected in the "Output profile" combobox;
        // gamma come from the selected profile, otherwise it comes from "Free gamma" tool

        //  readyImg = ipf.lab2rgb16 (labView, cx, cy, cw, ch, params.icm.output, params.blackwhite.enabled);
        bool bwonly = params.blackwhite.enabled &&  !params.colorToning.enabled ;

        if(autili || butili ) {
            bwonly = false;
        }

        readyImg = ipf.lab2rgb16 (labView, cx, cy, cw, ch, params.icm.output, params.icm.outputIntent, bwonly);

        if (settings->verbose) {
            printf("Output profile_: \"%s\"\n", params.icm.output.c_str());
        }
    }

    delete labView;
    labView = NULL;



    if(!autili && !butili ) {
        if(params.blackwhite.enabled &&  !params.colorToning.enabled ) {//force BW r=g=b
            if (settings->verbose) {
                printf("Force BW\n");
            }

            for (int ccw = 0; ccw < cw; ccw++) {
                for (int cch = 0; cch < ch; cch++) {
                    readyImg->r(cch, ccw) = readyImg->g(cch, ccw);
                    readyImg->b(cch, ccw) = readyImg->g(cch, ccw);
                }
            }
        }
    }

    if (pl) {
        pl->setProgress (0.70);
    }

    if (tmpScale != 1.0 && params.resize.method == "Nearest") { // resize rgb data (gamma applied)
        Image16* tempImage = new Image16 (imw, imh);
        ipf.resize (readyImg, tempImage, tmpScale);
        delete readyImg;
        readyImg = tempImage;
    }

    if (tunnelMetaData) {
        readyImg->setMetadata (ii->getMetaData()->getExifData ());
    } else {
        readyImg->setMetadata (ii->getMetaData()->getExifData (), params.exif, params.iptc);
    }


    // Setting the output curve to readyImg
    if (customGamma) {
        if (!useLCMS) {
            // use corrected sRGB profile in order to apply a good TRC if present, otherwise use LCMS2 profile generated by lab2rgb16b
            ProfileContent pc(jprof);
            readyImg->setOutputProfile (pc.data, pc.length);
        }
    } else {
        // use RT_sRGB.icm profile if present, otherwise use LCMS2 profile generate by lab2rgb16b
        Glib::ustring outputProfile;

        if (params.icm.output != "" && params.icm.output != ColorManagementParams::NoICMString) {
            outputProfile = params.icm.output;

            /*  if we'd wanted the RT_sRGB profile we would have selected it
            else {
            // use RT_sRGB.icm profile if present, otherwise use LCMS2 profile generate by lab2rgb16b
            if (settings->verbose) printf("No output profiles set ; looking for the default sRGB profile (\"%s\")...\n", options.rtSettings.srgb.c_str());
            outputProfile = options.rtSettings.srgb;
            }*/

            // if iccStore->getProfile send back an object, then iccStore->getContent will do too
            cmsHPROFILE jprof = iccStore->getProfile(outputProfile); //get outProfile

            if (jprof == NULL) {
                if (settings->verbose) {
                    printf("\"%s\" ICC output profile not found!\n - use LCMS2 substitution\n", outputProfile.c_str());
                }
            } else {
                if (settings->verbose) {
                    printf("Using \"%s\" output profile\n", outputProfile.c_str());
                }

                ProfileContent pc = iccStore->getContent (outputProfile);
                readyImg->setOutputProfile (pc.data, pc.length);
            }
        } else {
            // No ICM
            readyImg->setOutputProfile (NULL, 0);
        }
    }

//    t2.set();
//    if( settings->verbose )
//           printf("Total:- %d usec\n", t2.etime(t1));

    if (!job->initialImage) {
        ii->decreaseRef ();
    }

    delete job;

    if (pl) {
        pl->setProgress (0.75);
    }

    /*  curve1.reset();curve2.reset();
        curve.reset();
        satcurve.reset();
        lhskcurve.reset();

        rCurve.reset();
        gCurve.reset();
        bCurve.reset();
        hist16.reset();
        hist16C.reset();
    */
    return readyImg;
}

void batchProcessingThread (ProcessingJob* job, BatchProcessingListener* bpl, bool tunnelMetaData)
{

    ProcessingJob* currentJob = job;

    while (currentJob) {
        int errorCode;
        IImage16* img = processImage (currentJob, errorCode, bpl, tunnelMetaData, true);

        if (errorCode) {
            bpl->error (M("MAIN_MSG_CANNOTLOAD"));
            currentJob = NULL;
        } else {
            try {
                currentJob = bpl->imageReady (img);
            } catch (Glib::Exception& ex) {
                bpl->error (ex.what());
                currentJob = NULL;
            }
        }
    }
}

void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl, bool tunnelMetaData)
{

    if (bpl)
#if __GNUC__ == 4 && __GNUC_MINOR__ == 8 && defined( WIN32 ) && defined(__x86_64__)
        // See Issue 2384 "Very bad response time on win7/64 using gcc 4.8 when queue is running"
        Glib::Thread::create(sigc::bind(sigc::ptr_fun(batchProcessingThread), job, bpl, tunnelMetaData), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);

#else
        Glib::Thread::create(sigc::bind(sigc::ptr_fun(batchProcessingThread), job, bpl, tunnelMetaData), 0, true, true, Glib::THREAD_PRIORITY_LOW);
#endif

}

}
