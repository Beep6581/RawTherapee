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
#include <iostream>
#include <fstream>
#include <string>
#include "../rtgui/md5helper.h"
#include "../rtgui/thresholdselector.h"


#undef THREAD_PRIORITY_NORMAL

namespace rtengine
{
extern const Settings* settings;

namespace
{

template <typename T>
void adjust_radius(const T &default_param, double scale_factor, T &param)
{
    const double delta = (param - default_param) * scale_factor;
    param = default_param + delta;
}


class ImageProcessor
{
public:
    ImageProcessor(
        ProcessingJob* pjob,
        int& errorCode,
        ProgressListener* pl,
        bool flush
    ) :
        job(static_cast<ProcessingJobImpl*>(pjob)),
        errorCode(errorCode),
        pl(pl),
        flush(flush),
        // internal state
        ii(nullptr),
        imgsrc(nullptr),
        fw(0),
        fh(0),
        tr(0),
        pp(0, 0, 0, 0, 0),
        calclum(nullptr),
        autoNR(0.f),
        autoNRmax(0.f),
        tilesize(0),
        overlap(0),
        ch_M(nullptr),
        max_r(nullptr),
        max_b(nullptr),
        min_b(nullptr),
        min_r(nullptr),
        lumL(nullptr),
        chromC(nullptr),
        ry(nullptr),
        sk(nullptr),
        pcsk(nullptr),
        expcomp(0.0),
        bright(0),
        contr(0),
        black(0),
        hlcompr(0),
        hlcomprthresh(0),
        baseImg(nullptr),
        labView(nullptr),
        autili(false),
        butili(false)
    {
    }

    Imagefloat *operator()()
    {
        if (!job->fast) {
            return normal_pipeline();
        } else {
            return fast_pipeline();
        }
    }

private:
    Imagefloat *normal_pipeline()
    {
        if (!stage_init()) {
            return nullptr;
        }

        stage_denoise();
        stage_transform();
        return stage_finish();
    }

    Imagefloat *fast_pipeline()
    {
        if (!job->pparams.resize.enabled) {
            return normal_pipeline();
        }

        pl = nullptr;

        if (!stage_init()) {
            return nullptr;
        }

        stage_transform();
        stage_early_resize();
        stage_denoise();
        return stage_finish();
    }

    bool stage_init()
    {
        errorCode = 0;

        if (pl) {
            pl->setProgressStr("PROGRESSBAR_PROCESSING");
            pl->setProgress(0.0);
        }

        ii = job->initialImage;

        if (!ii) {
            ii = InitialImage::load(job->fname, job->isRaw, &errorCode);

            if (errorCode) {
                delete job;
                return false; //return nullptr;
            }
        }

        procparams::ProcParams& params = job->pparams;

        // acquire image from imagesource
        imgsrc = ii->getImageSource();

        tr = getCoarseBitMask(params.coarse);
        imgsrc->getFullSize(fw, fh, tr);

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

        ipf_p.reset(new ImProcFunctions(&params, true));
        ImProcFunctions &ipf = * (ipf_p.get());

        pp = PreviewProps(0, 0, fw, fh, 1);
        imgsrc->setCurrentFrame(params.raw.bayersensor.imageNum);
        imgsrc->preprocess(params.raw, params.lensProf, params.coarse, params.dirpyrDenoise.enabled);

        if (params.toneCurve.autoexp) {// this enabled HLRecovery
            LUTu histRedRaw(256), histGreenRaw(256), histBlueRaw(256);
            imgsrc->getRAWHistogram(histRedRaw, histGreenRaw, histBlueRaw);

            if (ToneCurveParams::HLReconstructionNecessary(histRedRaw, histGreenRaw, histBlueRaw) && !params.toneCurve.hrenabled) {
                params.toneCurve.hrenabled = true;
                // WARNING: Highlight Reconstruction is being forced 'on', should we force a method here too?
            }
        }

        if (pl) {
            pl->setProgress(0.20);
        }

        imgsrc->demosaic(params.raw);

        if (pl) {
            pl->setProgress(0.30);
        }

        if (params.retinex.enabled) { //enabled Retinex
            LUTf cdcurve(65536, 0);
            LUTf mapcurve(65536, 0);
            LUTu dummy;
            RetinextransmissionCurve dehatransmissionCurve;
            RetinexgaintransmissionCurve dehagaintransmissionCurve;
            bool dehacontlutili = false;
            bool mapcontlutili = false;
            bool useHsl = false;
//        multi_array2D<float, 3> conversionBuffer(1, 1);
            multi_array2D<float, 4> conversionBuffer(1, 1);
            imgsrc->retinexPrepareBuffers(params.icm, params.retinex, conversionBuffer, dummy);
            imgsrc->retinexPrepareCurves(params.retinex, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, dehacontlutili, mapcontlutili, useHsl, dummy, dummy);
            float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            imgsrc->retinex(params.icm, params.retinex, params.toneCurve, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, conversionBuffer, dehacontlutili, mapcontlutili, useHsl, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax, dummy);
        }

        if (pl) {
            pl->setProgress(0.40);
        }

        imgsrc->HLRecovery_Global(params.toneCurve);


        if (pl) {
            pl->setProgress(0.45);
        }

        // set the color temperature
        currWB = ColorTemp(params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

        if (!params.wb.enabled) {
            currWB = ColorTemp();
        } else if (params.wb.method == "Camera") {
            currWB = imgsrc->getWB();
        } else if (params.wb.method == "Auto") {
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers(rm, gm, bm);
            currWB.update(rm, gm, bm, params.wb.equal, params.wb.tempBias);
        }

        calclum = nullptr ;
        params.dirpyrDenoise.getCurves(noiseLCurve, noiseCCurve);
        autoNR = (float) settings->nrauto;//
        autoNRmax = (float) settings->nrautomax;//

        if (settings->leveldnti == 0) {
            tilesize = 1024;
            overlap = 128;
        }

        if (settings->leveldnti == 1) {
            tilesize = 768;
            overlap = 96;
        }

        //  const int tilesize = 768;
        //  const int overlap = 96;
        int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;
        ipf.Tile_calc(tilesize, overlap, 2, fw, fh, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);
        int nbtl = numtiles_W * numtiles_H;

        if ((settings->leveldnautsimpl == 1 && params.dirpyrDenoise.Cmethod == "AUT") || (settings->leveldnautsimpl == 0 && params.dirpyrDenoise.C2method == "AUTO")) {
            nbtl = 9;
        }

        ch_M = new float [nbtl];//allocate memory
        max_r = new float [nbtl];
        max_b = new float [nbtl];
        min_b = new float [9];
        min_r = new float [9];
        lumL = new float [nbtl];
        chromC = new float [nbtl];
        ry = new float [nbtl];
        sk = new float [nbtl];
        pcsk = new float [nbtl];

        //  printf("expert=%d\n",settings->leveldnautsimpl);
        if (settings->leveldnautsimpl == 1 && params.dirpyrDenoise.Cmethod == "PON") {
            MyTime t1pone, t2pone;
            t1pone.set();
            int crW = 100; // settings->leveldnv == 0
            int crH = 100; // settings->leveldnv == 0

            if (settings->leveldnv == 1) {
                crW = 250;
                crH = 250;
            }

            if (settings->leveldnv == 2) {
                crW = int (tileWskip / 2);
                crH = int (tileHskip / 2);
            }

            //  if(settings->leveldnv ==2) {crW=int(tileWskip/2);crH=int(1.15f*(tileWskip/2));}//adapted to scale of preview
            if (settings->leveldnv == 3) {
                crW = tileWskip - 10;
                crH = tileHskip - 10;
            }

            float lowdenoise = 1.f;
            int levaut = settings->leveldnaut;

            if (levaut == 1) { //Standard
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
                    origCropPart = new Imagefloat(crW, crH); //allocate memory
                    Imagefloat *provicalc = new Imagefloat((crW + 1) / 2, (crH + 1) / 2);  //for denoise curves
                    int skipP = 1;
                    #pragma omp for schedule(dynamic) collapse(2) nowait

                    for (int wcr = 0; wcr < numtiles_W; wcr++) {
                        for (int hcr = 0; hcr < numtiles_H; hcr++) {
                            int beg_tileW = wcr * tileWskip + tileWskip / 2.f - crW / 2.f;
                            int beg_tileH = hcr * tileHskip + tileHskip / 2.f - crH / 2.f;
                            PreviewProps ppP(beg_tileW, beg_tileH, crW, crH, skipP);
                            imgsrc->getImage(currWB, tr, origCropPart, ppP, params.toneCurve, params.raw);
                            //baseImg->getStdImage(currWB, tr, origCropPart, ppP, true, params.toneCurve);

                            // we only need image reduced to 1/4 here
                            for (int ii = 0; ii < crH; ii += 2) {
                                for (int jj = 0; jj < crW; jj += 2) {
                                    provicalc->r(ii >> 1, jj >> 1) = origCropPart->r(ii, jj);
                                    provicalc->g(ii >> 1, jj >> 1) = origCropPart->g(ii, jj);
                                    provicalc->b(ii >> 1, jj >> 1) = origCropPart->b(ii, jj);
                                }
                            }

                            imgsrc->convertColorSpace(provicalc, params.icm, currWB);  //for denoise luminance curve
                            float maxr = 0.f;
                            float maxb = 0.f;
                            float pondcorrec = 1.0f;
                            float chaut, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc;
                            int Nb;
                            chaut = 0.f;
                            redaut = 0.f;
                            blueaut = 0.f;
                            maxredaut = 0.f;
                            maxblueaut = 0.f;
                            chromina = 0.f;
                            sigma = 0.f;
                            ipf.RGB_denoise_info(origCropPart, provicalc, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope, params.dirpyrDenoise, imgsrc->getDirPyrDenoiseExpComp(), chaut, Nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc);
                            float multip = 1.f;
                            float adjustr = 1.f;

                            if (params.icm.working == "ProPhoto")   {
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

                            if (!imgsrc->isRAW()) {
                                multip = 2.f;    //take into account gamma for TIF / JPG approximate value...not good fot gamma=1
                            }

                            float maxmax = max(maxredaut, maxblueaut);
                            float delta;
                            int mode = 2;
                            int lissage = settings->leveldnliss;
                            ipf.calcautodn_info(chaut, delta, Nb, levaut, maxmax, lumema, chromina, mode, lissage, redyel, skinc, nsknc);

                            //    printf("PROCESS cha=%f red=%f bl=%f redM=%f bluM=%f chrom=%f sigm=%f lum=%f sigL=%f\n",chaut,redaut,blueaut, maxredaut, maxblueaut, chromina, sigma, lumema, sigma_L);
                            if (maxredaut > maxblueaut) {
                                maxr = (delta) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);

                                if (minblueaut <= minredaut  && minblueaut < chaut) {
                                    maxb = (-chaut + minblueaut) / (autoNRmax * multip * adjustr * lowdenoise);
                                }
                            } else {
                                maxb = (delta) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);

                                if (minredaut <= minblueaut  && minredaut < chaut) {
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

                if (liss == 2 || liss == 3) {
                    // I smooth only mean and not delta (max)
                    float nchm = 0.f;
                    float koef = 0.4f; //between 0.1 to 0.9

                    if (liss == 3) {
                        koef = 0.0f;    //quasi auto for mean Ch
                    }

                    for (int wcr = 0; wcr < numtiles_W; wcr++) {
                        for (int hcr = 0; hcr < numtiles_H; hcr++) {
                            nchm += ch_M[hcr * numtiles_W + wcr];
                        }
                    }

                    nchm /= (numtiles_H * numtiles_W);

                    for (int wcr = 0; wcr < numtiles_W; wcr++) {
                        for (int hcr = 0; hcr < numtiles_H; hcr++) {
                            ch_M[hcr * numtiles_W + wcr] = nchm + (ch_M[hcr * numtiles_W + wcr] - nchm) * koef;
                        }
                    }
                }

                if (liss == 3) { //same as auto but with much cells
                    float MaxR = 0.f;
                    float MaxB = 0.f;
                    float MaxRMoy = 0.f;
                    float MaxBMoy = 0.f;

                    for (int k = 0; k < nbtl; k++) {
                        MaxBMoy += max_b[k];
                        MaxRMoy += max_r[k];

                        if (max_r[k] > MaxR) {
                            MaxR = max_r[k];
                        }

                        if (max_b[k] > MaxB) {
                            MaxB = max_b[k];
                        }

                    }

                    MaxBMoy /= nbtl;
                    MaxRMoy /= nbtl;

                    for (int k = 0; k < nbtl; k++) {
                        if (MaxR > MaxB) {
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


        if ((settings->leveldnautsimpl == 1 && params.dirpyrDenoise.Cmethod == "AUT")  || (settings->leveldnautsimpl == 0 && params.dirpyrDenoise.C2method == "AUTO")) {
            MyTime t1aue, t2aue;
            t1aue.set();
            int crW, crH;

            if (settings->leveldnv == 0) {
                crW = 100;
                crH = 100;
            }

            if (settings->leveldnv == 1) {
                crW = 250;
                crH = 250;
            }

            if (settings->leveldnv == 2) {
                crW = int (tileWskip / 2);
                crH = int (tileHskip / 2);
            }

            //  if(settings->leveldnv ==2) {crW=int(tileWskip/2);crH=int(1.15f*(tileWskip/2));}//adapted to scale of preview
            if (settings->leveldnv == 3) {
                crW = tileWskip - 10;
                crH = tileHskip - 10;
            }

            float lowdenoise = 1.f;
            int levaut = settings->leveldnaut;

            if (levaut == 1) { //Standard
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
                    origCropPart = new Imagefloat(crW, crH); //allocate memory
                    Imagefloat *provicalc = new Imagefloat((crW + 1) / 2, (crH + 1) / 2);  //for denoise curves

                    #pragma omp for schedule(dynamic) collapse(2) nowait

                    for (int wcr = 0; wcr <= 2; wcr++) {
                        for (int hcr = 0; hcr <= 2; hcr++) {
                            PreviewProps ppP(coordW[wcr], coordH[hcr], crW, crH, 1);
                            imgsrc->getImage(currWB, tr, origCropPart, ppP, params.toneCurve, params.raw);
                            //baseImg->getStdImage(currWB, tr, origCropPart, ppP, true, params.toneCurve);


                            // we only need image reduced to 1/4 here
                            for (int ii = 0; ii < crH; ii += 2) {
                                for (int jj = 0; jj < crW; jj += 2) {
                                    provicalc->r(ii >> 1, jj >> 1) = origCropPart->r(ii, jj);
                                    provicalc->g(ii >> 1, jj >> 1) = origCropPart->g(ii, jj);
                                    provicalc->b(ii >> 1, jj >> 1) = origCropPart->b(ii, jj);
                                }
                            }

                            imgsrc->convertColorSpace(provicalc, params.icm, currWB);  //for denoise luminance curve
                            int nb = 0;
                            float chaut = 0.f, redaut = 0.f, blueaut = 0.f, maxredaut = 0.f, maxblueaut = 0.f, minredaut = 0.f, minblueaut = 0.f, chromina = 0.f, sigma = 0.f, lumema = 0.f, sigma_L = 0.f, redyel = 0.f, skinc = 0.f, nsknc = 0.f;
                            ipf.RGB_denoise_info(origCropPart, provicalc, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope,  params.dirpyrDenoise, imgsrc->getDirPyrDenoiseExpComp(), chaut, nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc);
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

                if (params.icm.working == "ProPhoto")   {
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

                if (!imgsrc->isRAW()) {
                    multip = 2.f;    //take into account gamma for TIF / JPG approximate value...not good fot gamma=1
                }

                float delta[9];
                int mode = 1;
                int lissage = settings->leveldnliss;

                for (int k = 0; k < 9; k++) {
                    float maxmax = max(max_r[k], max_b[k]);
                    ipf.calcautodn_info(ch_M[k], delta[k], Nb[k], levaut, maxmax, lumL[k], chromC[k], mode, lissage, ry[k], sk[k], pcsk[k]);
                    //  printf("ch_M=%f delta=%f\n",ch_M[k], delta[k]);
                }

                for (int k = 0; k < 9; k++) {
                    if (max_r[k] > max_b[k]) {
                        //printf("R delta=%f  koef=%f\n",delta[k],autoNRmax*multip*adjustr*lowdenoise);
                        Max_R[k] = (delta[k]) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);
                        Min_B[k] = - (ch_M[k] - min_b[k]) / (autoNRmax * multip * adjustr * lowdenoise);
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

                for (int k = 0; k < 9; k++) {
                    //  printf("ch_M= %f Max_R=%f Max_B=%f min_r=%f min_b=%f\n",ch_M[k],Max_R[k], Max_B[k],Min_R[k], Min_B[k]);
                    chM += ch_M[k];
                    MaxBMoy += Max_B[k];
                    MaxRMoy += Max_R[k];
                    MinRMoy += Min_R[k];
                    MinBMoy += Min_B[k];

                    if (Max_R[k] > MaxR) {
                        MaxR = Max_R[k];
                    }

                    if (Max_B[k] > MaxB) {
                        MaxB = Max_B[k];
                    }

                    if (Min_R[k] < MinR) {
                        MinR = Min_R[k];
                    }

                    if (Min_B[k] < MinB) {
                        MinB = Min_B[k];
                    }

                }

                chM /= 9;
                MaxBMoy /= 9;
                MaxRMoy /= 9;
                MinBMoy /= 9;
                MinRMoy /= 9;

                if (MaxR > MaxB) {
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

        baseImg = new Imagefloat(fw, fh);
        imgsrc->getImage(currWB, tr, baseImg, pp, params.toneCurve, params.raw);

        if (pl) {
            pl->setProgress(0.50);
        }

//  LUTf Noisecurve (65536,0);
//!!!// auto exposure!!!
        expcomp = params.toneCurve.expcomp;
        bright = params.toneCurve.brightness;
        contr = params.toneCurve.contrast;
        black = params.toneCurve.black;
        hlcompr = params.toneCurve.hlcompr;
        hlcomprthresh = params.toneCurve.hlcomprthresh;


        if (params.toneCurve.autoexp) {
            LUTu aehist;
            int aehistcompr;
            imgsrc->getAutoExpHistogram(aehist, aehistcompr);
            ipf.getAutoExp(aehist, aehistcompr, params.toneCurve.clip, expcomp, bright, contr, black, hlcompr, hlcomprthresh);
        }

        // at this stage, we can flush the raw data to free up quite an important amount of memory
        // commented out because it makes the application crash when batch processing...
        // TODO: find a better place to flush rawData and rawRGB
        if (flush) {
            imgsrc->flushRawData();
            imgsrc->flushRGB();
        }

        return true;
    }

    void stage_denoise()
    {
        procparams::ProcParams& params = job->pparams;
        //ImProcFunctions ipf (&params, true);
        ImProcFunctions &ipf = * (ipf_p.get());

        // perform luma/chroma denoise
//  CieImage *cieView;
//  NoisCurve noiseLCurve;
//    bool lldenoiseutili=false;
//  Imagefloat *calclum ;
//    params.dirpyrDenoise.getCurves(noiseLCurve, lldenoiseutili);
//  if (params.dirpyrDenoise.enabled  && lldenoiseutili) {

        DirPyrDenoiseParams denoiseParams = params.dirpyrDenoise;   // make a copy because we cheat here

        if (denoiseParams.Lmethod == "CUR") {
            if (noiseLCurve) {
                denoiseParams.luma = 0.5f;
            } else {
                denoiseParams.luma = 0.0f;
            }
        } else if (denoiseParams.Lmethod == "SLI") {
            noiseLCurve.Reset();
        }

        if (denoiseParams.enabled  && (noiseLCurve || noiseCCurve)) {
            // we only need image reduced to 1/4 here
            calclum = new Imagefloat((fw + 1) / 2, (fh + 1) / 2);  //for luminance denoise curve
            #pragma omp parallel for

            for (int ii = 0; ii < fh; ii += 2) {
                for (int jj = 0; jj < fw; jj += 2) {
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
            float nresi, highresi;
            int kall = 2;
            ipf.RGB_denoise(kall, baseImg, baseImg, calclum, ch_M, max_r, max_b, imgsrc->isRAW(), denoiseParams, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, noiseCCurve, nresi, highresi);

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
    }

    void stage_transform()
    {
        procparams::ProcParams& params = job->pparams;
        //ImProcFunctions ipf (&params, true);
        ImProcFunctions &ipf = * (ipf_p.get());

        imgsrc->convertColorSpace(baseImg, params.icm, currWB);

        // perform first analysis
        hist16(65536);

        ipf.firstAnalysis(baseImg, params, hist16);

        if (params.fattal.enabled) {
            ipf.ToneMapFattal02(baseImg);
        }

        // perform transform (excepted resizing)
        if (ipf.needsTransform()) {
            Imagefloat* trImg = nullptr;

            if (ipf.needsLuminanceOnly()) {
                trImg = baseImg;
            } else {
                trImg = new Imagefloat(fw, fh);
            }

            ipf.transform(baseImg, trImg, 0, 0, 0, 0, fw, fh, fw, fh,
                          imgsrc->getMetaData(), imgsrc->getRotateDegree(), true);

            if (trImg != baseImg) {
                delete baseImg;
                baseImg = trImg;
            }
        }
    }

    Imagefloat *stage_finish()
    {
        procparams::ProcParams& params = job->pparams;
        //ImProcFunctions ipf (&params, true);
        ImProcFunctions &ipf = * (ipf_p.get());

        if (params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled) {
            const int W = baseImg->getWidth();
            const int H = baseImg->getHeight();
            LabImage labcbdl(W, H);
            ipf.rgb2lab(*baseImg, labcbdl, params.icm.working);
            ipf.dirpyrequalizer(&labcbdl, 1);
            ipf.lab2rgb(labcbdl, *baseImg, params.icm.working);
        }

        // update blurmap
        SHMap* shmap = nullptr;

        if (params.sh.enabled) {
            shmap = new SHMap(fw, fh, true);
            double radius = sqrt(double (fw * fw + fh * fh)) / 2.0;
            double shradius = params.sh.radius;

            if (!params.sh.hq) {
                shradius *= radius / 1800.0;
            }

            shmap->update(baseImg, shradius, ipf.lumimul, params.sh.hq, 1);
        }

        // RGB processing

        curve1(65536);
        curve2(65536);
        curve(65536, 0);
        satcurve(65536, 0);
        lhskcurve(65536, 0);
        lumacurve(32770, 0);  // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
        clcurve(65536, 0);
        wavclCurve(65536, 0);

        //if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;

        CurveFactory::complexCurve(expcomp, black / 65535.0, hlcompr, hlcomprthresh, params.toneCurve.shcompr, bright, contr,
                                   params.toneCurve.curve, params.toneCurve.curve2,
                                   hist16, curve1, curve2, curve, dummy, customToneCurve1, customToneCurve2);

        CurveFactory::RGBCurve(params.rgbCurves.rcurve, rCurve, 1);
        CurveFactory::RGBCurve(params.rgbCurves.gcurve, gCurve, 1);
        CurveFactory::RGBCurve(params.rgbCurves.bcurve, bCurve, 1);

        bool opautili = false;

        if (params.colorToning.enabled) {
            TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix(params.icm.working);
            double wp[3][3] = {
                {wprof[0][0], wprof[0][1], wprof[0][2]},
                {wprof[1][0], wprof[1][1], wprof[1][2]},
                {wprof[2][0], wprof[2][1], wprof[2][2]}
            };
            params.colorToning.getCurves(ctColorCurve, ctOpacityCurve, wp, opautili);
            clToningcurve(65536, 0);
            CurveFactory::curveToning(params.colorToning.clcurve, clToningcurve, 1);
            cl2Toningcurve(65536, 0);
            CurveFactory::curveToning(params.colorToning.cl2curve, cl2Toningcurve, 1);
        }

        labView = new LabImage(fw, fh);
        reservView = new LabImage(fw, fh);

        if (params.blackwhite.enabled) {
            CurveFactory::curveBW(params.blackwhite.beforeCurve, params.blackwhite.afterCurve, hist16, dummy, customToneCurvebw1, customToneCurvebw2, 1);
        }

        double rrm, ggm, bbm;
        float autor, autog, autob;
        float satLimit = float (params.colorToning.satProtectionThreshold) / 100.f * 0.7f + 0.3f;
        float satLimitOpacity = 1.f - (float (params.colorToning.saturatedOpacity) / 100.f);

        if (params.colorToning.enabled  && params.colorToning.autosat && params.colorToning.method != "LabGrid") { //for colortoning evaluation of saturation settings
            float moyS = 0.f;
            float eqty = 0.f;
            ipf.moyeqt(baseImg, moyS, eqty); //return image : mean saturation and standard dev of saturation
            float satp = ((moyS + 1.5f * eqty) - 0.3f) / 0.7f; //1.5 sigma ==> 93% pixels with high saturation -0.3 / 0.7 convert to Hombre scale

            if (satp >= 0.92f) {
                satp = 0.92f;    //avoid values too high (out of gamut)
            }

            if (satp <= 0.15f) {
                satp = 0.15f;    //avoid too low values
            }

            satLimit = 100.f * satp;

            satLimitOpacity = 100.f * (moyS - 0.85f * eqty); //-0.85 sigma==>20% pixels with low saturation
        }

        autor = -9000.f; // This will ask to compute the "auto" values for the B&W tool (have to be inferior to -5000)
        DCPProfile::ApplyState as;
        DCPProfile *dcpProf = imgsrc->getDCP(params.icm, as);

        LUTu histToneCurve;

        ipf.rgbProc(baseImg, labView, nullptr, curve1, curve2, curve, shmap, params.toneCurve.saturation, rCurve, gCurve, bCurve, satLimit, satLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2, customToneCurvebw1, customToneCurvebw2, rrm, ggm, bbm, autor, autog, autob, expcomp, hlcompr, hlcomprthresh, dcpProf, as, histToneCurve);

        if (settings->verbose) {
            printf("Output image / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", autor, autog, autob);
        }

        // if clut was used and size of clut cache == 1 we free the memory used by the clutstore (default clut cache size = 1 for 32 bit OS)
        if (params.filmSimulation.enabled && !params.filmSimulation.clutFilename.empty() && options.clutCacheSize == 1) {
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
        baseImg = nullptr;

        if (shmap) {
            delete shmap;
        }

        shmap = nullptr;

        if (pl) {
            pl->setProgress(0.55);
        }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // start tile processing...???


        if (params.labCurve.contrast != 0) { //only use hist16 for contrast
            hist16.clear();

#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                LUTu hist16thr(hist16.getSize());   // one temporary lookup table per thread
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
        CurveFactory::complexLCurve(params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, hist16, lumacurve, dummy, 1, utili);

        bool clcutili;
        CurveFactory::curveCL(clcutili, params.labCurve.clcurve, clcurve, 1);

        bool ccutili, cclutili;
        CurveFactory::complexsgnCurve(autili, butili, ccutili, cclutili, params.labCurve.acurve, params.labCurve.bcurve, params.labCurve.cccurve,
                                      params.labCurve.lccurve, curve1, curve2, satcurve, lhskcurve, 1);


        //     bool locallutili = false;
        //     bool localcutili = false;
        reservView->CopyFrom(labView);

        if (params.locallab.enabled) {
            MyTime t1, t2;
            t1.set();

            std::string mdfive = getMD5(imgsrc->getFileName());

            Glib::ustring pop = options.cacheBaseDir + "/mip/";

            Glib::ustring datalab;

            if (options.mip == MI_opt) {
                datalab = pop + Glib::path_get_basename(imgsrc->getFileName() + "." + mdfive + ".mip");
            }

            if (options.mip == MI_prev) {
                datalab = imgsrc->getFileName() + ".mip";
            }


            LocretigainCurve locRETgainCurve;
            LocLHCurve loclhCurve;
            LocHHCurve lochhCurve;

            LocretigainCurverab locRETgainCurverab;
            LUTf lllocalcurve(65536, 0);
            LUTf cclocalcurve(65536, 0);
            LUTf sklocalcurve(65536, 0);
            LUTf hltonecurveloc(32768, 0);
            LUTf shtonecurveloc(32768, 0);
            LUTf tonecurveloc(32768, 0);
            LUTf exlocalcurve(32768, 0);
            //    int realspot = params.locallab.nbspot;
            int maxspot = settings->nspot + 1;
            ifstream fic0(datalab, ios::in);
            float** shbuffer = nullptr;
            int versionmip = 0;
            std::string delim[69] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
                                     "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
                                     "&", "#", "{", "[", "]", "}", "$", "*", "?", ">", "!", ";", "<"",(", ")", "+", "-"
                                    };

            if (params.locallab.inverssha) {
                shbuffer   = new float*[fh];

                for (int i = 0; i < fh; i++) {
                    shbuffer[i] = new float[fw];
                }
            }

            if (fic0) {//normally we don't use here but ??
                //find the version mip
                string line;
                string spotline;
                //   int cont = 0;

                while (getline(fic0, line)) {
                    spotline = line;
                    std::size_t pos = spotline.find("=");
                    std::size_t posend = spotline.find("@");  //in case of for futur use

                    if (spotline.substr(0, pos) == "Mipversion") {
                        string strversion = spotline.substr(pos + 1, (posend - pos));
                        versionmip = std::stoi(strversion.c_str());
                    }

                    if (spotline.substr(0, pos) == "Spot") {
                        //        cont = 0;
                    }


                }

                fic0.close();
            }

            ifstream fich(datalab, ios::in);
            int maxdata = 102; //101 10023 //100 10022 //99 10021 // 91 10021 //88 10019 //87 10018//86 10017//85 10016 //82;//78;//73 10011

            if (fich && versionmip != 0) {
                std::string inser;

                int **dataspots;
                dataspots = new int*[maxdata];

                for (int i = 0; i < maxdata; i++) {
                    dataspots[i] = new int[maxspot];
                }


                std::string *retistrs;
                retistrs = new std::string[maxspot];
                std::string *llstrs;

                llstrs = new std::string[maxspot];

                std::string *lhstrs;
                lhstrs = new std::string[maxspot];

                std::string *ccstrs;
                ccstrs = new std::string[maxspot];

                std::string *hhstrs;
                hhstrs = new std::string[maxspot];

                std::string *skinstrs;
                skinstrs = new std::string[maxspot];

                std::string *pthstrs;
                pthstrs = new std::string[maxspot];

                std::string *exstrs;
                exstrs = new std::string[maxspot];

                {
                    dataspots[2][0] =  params.locallab.circrad;
                    dataspots[3][0] =  params.locallab.locX;
                    dataspots[4][0] =  params.locallab.locY;
                    dataspots[5][0] =  params.locallab.locYT;
                    dataspots[6][0] =  params.locallab.locXL;
                    dataspots[7][0] =  params.locallab.centerX;
                    dataspots[8][0] =  params.locallab.centerY;
                    dataspots[9][0] =  params.locallab.lightness;
                    dataspots[10][0] =  params.locallab.contrast;
                    dataspots[11][0] =  params.locallab.chroma;
                    dataspots[12][0] =  params.locallab.sensi;
                    dataspots[13][0] =  params.locallab.transit;

                    if (!params.locallab.invers) {
                        dataspots[14][0] =  0;
                    } else {
                        dataspots[14][0] =  1;
                    }

                    if (params.locallab.Smethod == "IND") {
                        dataspots[15][0] =  0;
                    } else if (params.locallab.Smethod == "SYM") {
                        dataspots[15][0] =  1;
                    } else if (params.locallab.Smethod == "INDSL") {
                        dataspots[15][0] =  2;
                    } else if (params.locallab.Smethod == "SYMSL") {
                        dataspots[15][0] =  3;
                    }

                    dataspots[17][0] =  params.locallab.radius;
                    dataspots[18][0] =  params.locallab.strength;
                    dataspots[19][0] =  params.locallab.sensibn;

                    if (!params.locallab.inversrad) {
                        dataspots[20][0] =  0;
                    } else {
                        dataspots[20][0] =  1;
                    }

                    dataspots[21][0] = params.locallab.str;
                    dataspots[22][0] = params.locallab.chrrt;
                    dataspots[23][0] = params.locallab.neigh;
                    dataspots[24][0] = params.locallab.vart;
                    dataspots[25][0] = params.locallab.sensih;

                    if (!params.locallab.inversret) {
                        dataspots[26][0] =  0;
                    } else {
                        dataspots[26][0] =  1;
                    }

                    if (params.locallab.retinexMethod == "low") {
                        dataspots[27][0] =  0;
                    } else if (params.locallab.retinexMethod == "uni") {
                        dataspots[27][0] =  1;
                    } else if (params.locallab.retinexMethod == "high") {
                        dataspots[27][0] =  2;
                    }


                    dataspots[28][0] = params.locallab.sharradius;
                    dataspots[29][0] = params.locallab.sharamount;
                    dataspots[30][0] = params.locallab.shardamping;
                    dataspots[31][0] = params.locallab.shariter;
                    dataspots[32][0] = params.locallab.sensisha;

                    if (!params.locallab.inverssha) {
                        dataspots[33][0] =  0;
                    } else {
                        dataspots[33][0] =  1;
                    }

                    if (params.locallab.qualityMethod == "std") {
                        dataspots[34][0] =  0;
                    } else if (params.locallab.qualityMethod == "enh") {
                        dataspots[34][0] =  1;
                    } else if (params.locallab.qualityMethod == "enhden") {
                        dataspots[34][0] =  2;
                    }

                    dataspots[35][0] = params.locallab.thres;
                    dataspots[36][0] = params.locallab.proxi;

                    dataspots[37][0] = params.locallab.noiselumf;
                    dataspots[38][0] = params.locallab.noiselumc;
                    dataspots[39][0] = params.locallab.noisechrof;
                    dataspots[40][0] = params.locallab.noisechroc;

                    dataspots[41][0] = params.locallab.mult[0];
                    dataspots[42][0] = params.locallab.mult[1];
                    dataspots[43][0] = params.locallab.mult[2];
                    dataspots[44][0] = params.locallab.mult[3];
                    dataspots[45][0] = params.locallab.mult[4];
                    dataspots[46][0] = params.locallab.threshold;
                    dataspots[47][0] = params.locallab.sensicb;

                    if (!params.locallab.activlum) {
                        dataspots[48][0] =  0;
                    } else {
                        dataspots[48][0] =  1;
                    }

                    dataspots[49][0] = params.locallab.stren;
                    dataspots[50][0] = params.locallab.gamma;
                    dataspots[51][0] = params.locallab.estop;
                    dataspots[52][0] = params.locallab.scaltm;
                    dataspots[53][0] = params.locallab.rewei;
                    dataspots[54][0] = params.locallab.sensitm;
                    dataspots[55][0] = params.locallab.retrab;

                    if (!params.locallab.curvactiv) {
                        dataspots[56][0] =  0;
                    } else {
                        dataspots[56][0] =  1;
                    }

                    if (params.locallab.qualitycurveMethod == "none") {
                        dataspots[57][0] =  0;
                    } else if (params.locallab.qualitycurveMethod == "std") {
                        dataspots[57][0] =  1;
                    } else if (params.locallab.qualitycurveMethod == "enh") {
                        dataspots[57][0] =  2;
                    }

                    dataspots[58][0] = params.locallab.sensiv;
                    dataspots[59][0] = params.locallab.pastels;
                    dataspots[60][0] = params.locallab.saturated;

                    if (!params.locallab.protectskins) {
                        dataspots[61][0] = 0;
                    } else {
                        dataspots[61][0] = 1;
                    }

                    if (!params.locallab.avoidcolorshift) {
                        dataspots[62][0] = 0;
                    } else {
                        dataspots[62][0] = 1;
                    }

                    if (!params.locallab.pastsattog) {
                        dataspots[63][0] = 0;
                    } else {
                        dataspots[63][0] = 1;
                    }

                    dataspots[64][0] = params.locallab.expcomp;
                    dataspots[65][0] = params.locallab.black;
                    dataspots[66][0] = params.locallab.hlcompr;
                    dataspots[67][0] = params.locallab.hlcomprthresh;
                    dataspots[68][0] = params.locallab.shcompr;
                    dataspots[69][0] = params.locallab.sensiex;

                    dataspots[70][0] = params.locallab.centerXbuf;
                    dataspots[71][0] = params.locallab.centerYbuf;
                    dataspots[72][0] = params.locallab.adjblur;

                    if (!params.locallab.cutpast) {
                        dataspots[73][0] = 0;
                    } else {
                        dataspots[73][0] = 1;
                    }

                    dataspots[74][0] = params.locallab.chromacbdl;

                    if (!params.locallab.lastdust) {
                        dataspots[75][0] = 0;
                    } else {
                        dataspots[75][0] = 1;
                    }

                    if (params.locallab.blurMethod == "norm") {
                        dataspots[76][0] =  0;
                    } else if (params.locallab.blurMethod == "inv") {
                        dataspots[76][0] =  1;
                    } else if (params.locallab.blurMethod == "sym") {
                        dataspots[76][0] =  2;
                    }

                    if (params.locallab.dustMethod == "cop") {
                        dataspots[77][0] =  0;
                    } else if (params.locallab.dustMethod == "mov") {
                        dataspots[77][0] =  1;
                    } else if (params.locallab.dustMethod == "pas") {
                        dataspots[77][0] =  2;
                    }


                    if (params.locallab.Exclumethod == "norm") {
                        dataspots[78][0] =  0;
                    } else if (params.locallab.Exclumethod == "exc") {
                        dataspots[78][0] =  1;
                    }

                    dataspots[79][0] = params.locallab.sensiexclu;
                    dataspots[80][0] = params.locallab.struc;
                    dataspots[81][0] = params.locallab.warm;
                    dataspots[82][0] = params.locallab.noiselumdetail;
                    dataspots[83][0] = params.locallab.noisechrodetail;
                    dataspots[84][0] = params.locallab.sensiden;

                    if (!params.locallab.expdenoi) {
                        dataspots[85][0] =  0;
                    } else {
                        dataspots[85][0] =  1;
                    }

                    if (!params.locallab.expcolor) {
                        dataspots[86][0] =  0;
                    } else {
                        dataspots[86][0] =  1;
                    }

                    if (!params.locallab.expvibrance) {
                        dataspots[87][0] =  0;
                    } else {
                        dataspots[87][0] =  1;
                    }

                    if (!params.locallab.expblur) {
                        dataspots[88][0] =  0;
                    } else {
                        dataspots[88][0] =  1;
                    }

                    if (!params.locallab.exptonemap) {
                        dataspots[89][0] =  0;
                    } else {
                        dataspots[89][0] =  1;
                    }

                    if (!params.locallab.expreti) {
                        dataspots[90][0] =  0;
                    } else {
                        dataspots[90][0] =  1;
                    }

                    if (!params.locallab.expsharp) {
                        dataspots[91][0] =  0;
                    } else {
                        dataspots[91][0] =  1;
                    }

                    if (!params.locallab.expcbdl) {
                        dataspots[92][0] =  0;
                    } else {
                        dataspots[92][0] =  1;
                    }

                    if (!params.locallab.expexpose) {
                        dataspots[93][0] =  0;
                    } else {
                        dataspots[93][0] =  1;
                    }

                    dataspots[94][0] = params.locallab.bilateral;
                    dataspots[95][0] = params.locallab.noiselequal;

                    if (params.locallab.shapemethod == "ELI") {
                        dataspots[96][0] =  0;
                    } else if (params.locallab.shapemethod == "RECT") {
                        dataspots[96][0] =  1;
                    }

                    dataspots[maxdata - 5][0] = 100.f * params.locallab.huerefblur;
                    dataspots[maxdata - 4][0] = 100.f * params.locallab.hueref;
                    dataspots[maxdata - 3][0] = params.locallab.chromaref;
                    dataspots[maxdata - 2][0] = params.locallab.lumaref;
                    dataspots[maxdata - 1][0] = params.locallab.sobelref;

                    //curve Reti local
                    int siz = params.locallab.localTgaincurve.size();

                    if (siz > 69) {
                        siz = 69;//avoid crash
                    }

                    //   int s_cur[siz + 1];
                    int s_datcur[siz + 1];

                    for (int j = 0; j < siz; j++) {
                        s_datcur[j] = (int)(1000. * params.locallab.localTgaincurve[j]);
                    }

                    std::string cur_str = "";

                    for (int j = 0; j < siz; j++) {
                        cur_str = cur_str + std::to_string(s_datcur[j]) + delim[j];
                    }

                    inser = retistrs[0] = cur_str + "@";

                    int sizl = params.locallab.llcurve.size();

                    if (sizl > 69) {
                        sizl = 69;
                    }

                    //    int s_curl[sizl + 1];
                    int s_datcurl[sizl + 1];

                    for (int j = 0; j < sizl; j++) {
                        s_datcurl[j] = (int)(1000. * params.locallab.llcurve[j]);
                    }

                    std::string ll_str = "";

                    for (int j = 0; j < sizl; j++) {
                        ll_str = ll_str + std::to_string(s_datcurl[j]) + delim[j];
                    }

                    llstrs[0] = ll_str + "@";


                    int sizc = params.locallab.cccurve.size();

                    if (sizc > 69) {
                        sizc = 69;
                    }

                    //        int s_curc[sizc + 1];
                    int s_datcurc[sizc + 1];

                    for (int j = 0; j < sizc; j++) {
                        s_datcurc[j] = (int)(1000. * params.locallab.cccurve[j]);
                    }

                    std::string cc_str = "";

                    for (int j = 0; j < sizc; j++) {
                        cc_str = cc_str + std::to_string(s_datcurc[j]) + delim[j];
                    }

                    ccstrs[0] = cc_str + "@";

                    //

                    int sizh = params.locallab.LHcurve.size();

                    if (sizh > 69) {
                        sizh = 69;
                    }

                    //     int s_curh[sizh + 1];
                    int s_datcurh[sizh + 1];

                    for (int j = 0; j < sizh; j++) {
                        s_datcurh[j] = (int)(1000. * params.locallab.LHcurve[j]);
                    }

                    std::string lh_str = "";

                    for (int j = 0; j < sizh; j++) {
                        lh_str = lh_str + std::to_string(s_datcurh[j]) + delim[j];
                    }

                    lhstrs[0] = lh_str + "@";

                    int sizhh = params.locallab.HHcurve.size();

                    if (sizhh > 69) {
                        sizhh = 69;
                    }

                    //     int s_curh[sizh + 1];
                    int s_datcurhh[sizhh + 1];

                    for (int j = 0; j < sizhh; j++) {
                        s_datcurhh[j] = (int)(1000. * params.locallab.HHcurve[j]);
                    }

                    std::string hh_str = "";

                    for (int j = 0; j < sizhh; j++) {
                        hh_str = hh_str + std::to_string(s_datcurhh[j]) + delim[j];
                    }

                    hhstrs[0] = hh_str + "@";

                    //Skin curve
                    int sizsk = params.locallab.skintonescurve.size();

                    if (sizsk > 69) {
                        sizsk = 69;//to avoid crash
                    }


                    int s_datcursk[sizsk + 1];

                    for (int j = 0; j < sizsk; j++) {
                        s_datcursk[j] = (int)(1000. * params.locallab.skintonescurve[j]);
                    }

                    std::string sk_str = "";

                    for (int j = 0; j < sizsk; j++) {
                        sk_str = sk_str + std::to_string(s_datcursk[j]) + delim[j];
                    }

                    skinstrs[0] = sk_str + "@";

                    //end local skin
                    //PSThreshold
                    int sizps = 2;
                    int s_datps[sizps + 1];
                    s_datps[1] =  static_cast<int>(params.locallab.psthreshold.getTopLeft());

                    s_datps[0] =  static_cast<int>(params.locallab.psthreshold.getBottomLeft());

                    std::string ps_str = "";

                    for (int j = 0; j < sizps; j++) {
                        ps_str = ps_str + std::to_string(s_datps[j]) + delim[j];
                    }

                    pthstrs[0] = ps_str + "@";

                    //end local ps

                    //expos
                    //Skin curve
                    int sizex = params.locallab.excurve.size();

                    if (sizex > 69) {
                        sizex = 69;//to avoid crash
                    }


                    int s_datcurex[sizsk + 1];

                    for (int j = 0; j < sizex; j++) {
                        s_datcurex[j] = (int)(1000. * params.locallab.excurve[j]);
                    }

                    std::string ex_str = "";

                    for (int j = 0; j < sizex; j++) {
                        ex_str = ex_str + std::to_string(s_datcurex[j]) + delim[j];
                    }

                    exstrs[0] = ex_str + "@";




                }


                int ns = 0;

                if (fich) {

                    std::string line;
                    std::string spotline;
                    int cont = 0;

                    while (getline(fich, line)) {
                        spotline = line;
                        std::size_t pos = spotline.find("=");
                        std::size_t posend = spotline.find("@");  //in case of for futur use

                        if (spotline.substr(0, pos) == "Mipversion") {
                            std::string strversion = spotline.substr(pos + 1, (posend - pos));
                            versionmip = std::stoi(strversion.c_str());
                        }

                        if (spotline.substr(0, pos) == "Spot") {
                            cont = 0;
                        }

                        cont++;
                        std::string str3 = spotline.substr(pos + 1, (posend - pos));

                        if (cont == 1) {
                            ns =  std::stoi(str3.c_str());
                        }

                        if (cont >= 2  && cont < 16) {
                            dataspots[cont][ns] = std::stoi(str3.c_str());

                        }

                        if (spotline.substr(0, pos) == "Currentspot") {
                            dataspots[16][0] = std::stoi(str3.c_str());
                        }

                        if (cont > 16  && cont < maxdata) {
                            dataspots[cont][ns] = std::stoi(str3.c_str());

                        }

                        if (spotline.substr(0, pos) == "curveReti") {
                            retistrs[ns] = str3;
                        }

                        if (spotline.substr(0, pos) == "curveLL") {
                            llstrs[ns] = str3;
                        }

                        if (spotline.substr(0, pos) == "curveLH") {
                            lhstrs[ns] = str3;
                        }

                        if (spotline.substr(0, pos) == "curveHH") {
                            hhstrs[ns] = str3;
                        }

                        if (spotline.substr(0, pos) == "curveCC") {
                            ccstrs[ns] = str3;
                        }

                        if (spotline.substr(0, pos) == "curveskin") {
                            skinstrs[ns] = str3;
                        }

                        if (spotline.substr(0, pos) == "pthres") {
                            pthstrs[ns] = str3;
                        }

                        if (spotline.substr(0, pos) == "curveex") {
                            exstrs[ns] = str3;
                        }


                    }

                    fich.close();
                }


                for (int sp = 1; sp < maxspot; sp++) { //spots default
                    params.locallab.huerefblur = INFINITY;
                    params.locallab.hueref = INFINITY;
                    params.locallab.chromaref = INFINITY;
                    params.locallab.lumaref = INFINITY;
                    params.locallab.sobelref = INFINITY;

                    params.locallab.circrad = dataspots[2][sp];
                    params.locallab.locX = dataspots[3][sp];
                    params.locallab.locY = dataspots[4][sp];
                    params.locallab.locYT = dataspots[5][sp];
                    params.locallab.locXL = dataspots[6][sp];
                    params.locallab.centerX = dataspots[7][sp];
                    params.locallab.centerY = dataspots[8][sp];
                    params.locallab.lightness = dataspots[9][sp];
                    params.locallab.contrast = dataspots[10][sp];
                    params.locallab.chroma = dataspots[11][sp];
                    params.locallab.sensi = dataspots[12][sp];
                    params.locallab.transit = dataspots[13][sp];

                    if (dataspots[14][sp] ==  0) {
                        params.locallab.invers = false;
                    } else {
                        params.locallab.invers = true;
                    }

                    if (dataspots[15][sp] ==  0) {
                        params.locallab.Smethod = "IND" ;
                    } else if (dataspots[15][sp] ==  1) {
                        params.locallab.Smethod = "SYM" ;
                    } else if (dataspots[15][sp] ==  2) {
                        params.locallab.Smethod = "INDSL";
                    } else if (dataspots[15][sp] ==  3) {
                        params.locallab.Smethod = "SYMSL";
                    }

                    params.locallab.radius = dataspots[17][sp];
                    params.locallab.strength = dataspots[18][sp];
                    params.locallab.sensibn = dataspots[19][sp];

                    if (dataspots[20][sp] ==  0) {
                        params.locallab.inversrad = false;
                    } else {
                        params.locallab.inversrad = true;
                    }

                    params.locallab.str = dataspots[21][sp];
                    params.locallab.chrrt = dataspots[22][sp];
                    params.locallab.neigh = dataspots[23][sp];
                    params.locallab.vart = dataspots[24][sp];
                    params.locallab.sensih = dataspots[25][sp];

                    if (dataspots[26][sp] ==  0) {
                        params.locallab.inversret = false;
                    } else {
                        params.locallab.inversret = true;
                    }

                    if (dataspots[27][sp] ==  0) {
                        params.locallab.retinexMethod = "low" ;
                    } else if (dataspots[27][sp] ==  1) {
                        params.locallab.retinexMethod = "uni" ;
                    } else if (dataspots[27][sp] ==  2) {
                        params.locallab.retinexMethod = "high";
                    }

                    params.locallab.sharradius = dataspots[28][sp];
                    params.locallab.sharamount = dataspots[29][sp];
                    params.locallab.shardamping = dataspots[30][sp];
                    params.locallab.shariter = dataspots[31][sp];
                    params.locallab.sensisha = dataspots[32][sp];

                    if (dataspots[33][sp] ==  0) {
                        params.locallab.inverssha = false;
                    } else {
                        params.locallab.inverssha = true;
                    }

                    if (dataspots[34][sp] ==  0) {
                        params.locallab.qualityMethod = "std" ;
                    } else if (dataspots[34][sp] ==  1) {
                        params.locallab.qualityMethod = "enh" ;
                    } else if (dataspots[34][sp] ==  2) {
                        params.locallab.qualityMethod = "enhden" ;
                    }

                    params.locallab.thres = dataspots[35][sp];
                    params.locallab.proxi = dataspots[36][sp];

                    params.locallab.noiselumf = dataspots[37][sp];
                    params.locallab.noiselumc = dataspots[38][sp];
                    params.locallab.noisechrof = dataspots[39][sp];
                    params.locallab.noisechroc = dataspots[40][sp];

                    params.locallab.mult[0] = dataspots[41][sp];
                    params.locallab.mult[1] = dataspots[42][sp];
                    params.locallab.mult[2] = dataspots[43][sp];
                    params.locallab.mult[3] = dataspots[44][sp];
                    params.locallab.mult[4] = dataspots[45][sp];
                    params.locallab.threshold = dataspots[46][sp];
                    params.locallab.sensicb = dataspots[47][sp];

                    if (dataspots[48][sp] ==  0) {
                        params.locallab.activlum = false;
                    } else {
                        params.locallab.activlum = true;
                    }

                    params.locallab.stren = dataspots[49][sp];
                    params.locallab.gamma = dataspots[50][sp];
                    params.locallab.estop = dataspots[51][sp];
                    params.locallab.scaltm = dataspots[52][sp];
                    params.locallab.rewei = dataspots[53][sp];
                    params.locallab.sensitm = dataspots[54][sp];
                    params.locallab.retrab = dataspots[55][sp];

                    if (dataspots[56][sp] ==  0) {
                        params.locallab.curvactiv = false;
                    } else {
                        params.locallab.curvactiv = true;
                    }

                    if (dataspots[57][sp] ==  0) {
                        params.locallab.qualitycurveMethod = "none" ;
                    } else if (dataspots[57][sp] ==  1) {
                        params.locallab.qualitycurveMethod = "std" ;
                    } else if (dataspots[57][sp] ==  2) {
                        params.locallab.qualitycurveMethod = "enh" ;
                    }


                    params.locallab.sensiv = dataspots[58][sp];
                    params.locallab.pastels = dataspots[59][sp];
                    params.locallab.saturated = dataspots[60][sp];

                    if (dataspots[61][sp] ==  0) {
                        params.locallab.protectskins = false;
                    } else {
                        params.locallab.protectskins  = true;
                    }

                    if (dataspots[62][sp] ==  0) {
                        params.locallab.avoidcolorshift = false;
                    } else {
                        params.locallab.avoidcolorshift  = true;
                    }

                    if (dataspots[63][sp] ==  0) {
                        params.locallab.pastsattog = false;
                    } else {
                        params.locallab.pastsattog  = true;
                    }

                    params.locallab.expcomp = dataspots[64][sp];
                    params.locallab.black = dataspots[65][sp];
                    params.locallab.hlcompr = dataspots[66][sp];
                    params.locallab.hlcomprthresh = dataspots[67][sp];
                    params.locallab.shcompr =  dataspots[68][sp];
                    params.locallab.sensiex =  dataspots[69][sp];

                    params.locallab.centerXbuf  = dataspots[70][sp];
                    params.locallab.centerYbuf =  dataspots[71][sp];
                    params.locallab.adjblur = dataspots[72][sp];

                    if (dataspots[73][sp] ==  0) {
                        params.locallab.cutpast = false;
                    } else {
                        params.locallab.cutpast  = true;
                    }

                    params.locallab.chromacbdl = dataspots[74][sp];

                    if (dataspots[75][sp] ==  0) {
                        params.locallab.lastdust = false;
                    } else {
                        params.locallab.lastdust  = true;
                    }

                    if (dataspots[76][sp] ==  0) {
                        params.locallab.blurMethod = "norm" ;
                    } else if (dataspots[76][sp] ==  1) {
                        params.locallab.blurMethod = "inv" ;
                    } else if (dataspots[76][sp] ==  2) {
                        params.locallab.blurMethod = "sym" ;
                    }

                    if (dataspots[77][sp] ==  0) {
                        params.locallab.dustMethod = "cop" ;
                    } else if (dataspots[77][sp] ==  1) {
                        params.locallab.dustMethod = "mov" ;
                    } else if (dataspots[77][sp] ==  2) {
                        params.locallab.dustMethod = "pas" ;
                    }


                    if (dataspots[78][sp] ==  0) {
                        params.locallab.Exclumethod = "norm" ;
                    } else if (dataspots[78][sp] ==  1) {
                        params.locallab.Exclumethod = "exc" ;
                    }

                    params.locallab.sensiexclu = dataspots[79][sp];
                    params.locallab.struc = dataspots[80][sp];
                    params.locallab.warm = dataspots[81][sp];
                    params.locallab.noiselumdetail = dataspots[82][sp];
                    params.locallab.noisechrodetail = dataspots[83][sp];
                    params.locallab.sensiden = dataspots[84][sp];

                    if (dataspots[85][sp] ==  0) {
                        params.locallab.expdenoi = false;
                    } else {
                        params.locallab.expdenoi = true;
                    }

                    if (dataspots[86][sp] ==  0) {
                        params.locallab.expcolor = false;
                    } else {
                        params.locallab.expcolor = true;
                    }

                    if (dataspots[87][sp] ==  0) {
                        params.locallab.expvibrance = false;
                    } else {
                        params.locallab.expvibrance = true;
                    }

                    if (dataspots[88][sp] ==  0) {
                        params.locallab.expblur = false;
                    } else {
                        params.locallab.expblur = true;
                    }

                    if (dataspots[89][sp] ==  0) {
                        params.locallab.exptonemap = false;
                    } else {
                        params.locallab.exptonemap = true;
                    }

                    if (dataspots[90][sp] ==  0) {
                        params.locallab.expreti = false;
                    } else {
                        params.locallab.expreti = true;
                    }

                    if (dataspots[91][sp] ==  0) {
                        params.locallab.expsharp = false;
                    } else {
                        params.locallab.expsharp = true;
                    }

                    if (dataspots[92][sp] ==  0) {
                        params.locallab.expcbdl = false;
                    } else {
                        params.locallab.expcbdl = true;
                    }

                    if (dataspots[93][sp] ==  0) {
                        params.locallab.expexpose = false;
                    } else {
                        params.locallab.expexpose = true;
                    }

                    params.locallab.bilateral = dataspots[94][sp];
                    params.locallab.noiselequal = dataspots[95][sp];

                    if (dataspots[96][sp] ==  0) {
                        params.locallab.shapemethod = "ELI" ;
                    } else if (dataspots[96][sp] ==  1) {
                        params.locallab.shapemethod = "RECT" ;
                    }

                    params.locallab.huerefblur = ((float) dataspots[maxdata - 5][sp]) / 100.f;
                    params.locallab.hueref = ((float) dataspots[maxdata - 4][sp]) / 100.f;
                    params.locallab.chromaref = dataspots[maxdata - 3][sp];
                    params.locallab.lumaref = dataspots[maxdata - 2][sp];
                    params.locallab.sobelref = dataspots[maxdata - 1][sp];


                    int *s_datc;
                    s_datc = new int[70];
                    int siz;

                    ipf.strcurv_data(retistrs[sp], s_datc, siz);
                    std::vector<double>   cretiend;

                    for (int j = 0; j < siz; j++) {
                        cretiend.push_back((double)(s_datc[j]) / 1000.);
                    }

                    delete [] s_datc;

                    int *s_datcl;
                    s_datcl = new int[70];
                    int sizl;

                    ipf.strcurv_data(llstrs[sp], s_datcl, sizl);


                    std::vector<double>   cllend;

                    for (int j = 0; j < sizl; j++) {
                        cllend.push_back((double)(s_datcl[j]) / 1000.);
                    }

                    delete [] s_datcl;

                    int *s_datcc;
                    s_datcc = new int[70];
                    int sizc;

                    ipf.strcurv_data(ccstrs[sp], s_datcc, sizc);


                    std::vector<double>   cccend;

                    for (int j = 0; j < sizc; j++) {
                        cccend.push_back((double)(s_datcc[j]) / 1000.);
                    }

                    delete [] s_datcc;

                    int *s_datch;
                    s_datch = new int[70];
                    int sizh;

                    ipf.strcurv_data(lhstrs[sp], s_datch, sizh);


                    std::vector<double>   clhend;

                    for (int j = 0; j < sizh; j++) {
                        clhend.push_back((double)(s_datch[j]) / 1000.);
                    }

                    int *s_datchh;
                    s_datchh = new int[70];
                    int sizhh;

                    ipf.strcurv_data(hhstrs[sp], s_datchh, sizhh);


                    std::vector<double>   chhend;

                    for (int j = 0; j < sizhh; j++) {
                        chhend.push_back((double)(s_datchh[j]) / 1000.);
                    }

                    delete [] s_datchh;

                    int *s_datcsk;
                    s_datcsk = new int[70];
                    int sizsk;

                    ipf.strcurv_data(skinstrs[sp], s_datcsk, sizsk);


                    std::vector<double>   cskend;

                    for (int j = 0; j < sizsk; j++) {
                        cskend.push_back((double)(s_datcsk[j]) / 1000.);
                    }

                    delete [] s_datcsk;

                    //PSThreshold + 1
                    int sizps = 2;
                    int s_datcps[sizps + 1];
                    ipf.strcurv_data(pthstrs[sp], s_datcps, sizps);

                    params.locallab.psthreshold.setValues(s_datcps[0], s_datcps[1]);

                    //expos
                    int *s_datcex;
                    s_datcex = new int[70];
                    int sizex;

                    ipf.strcurv_data(exstrs[sp], s_datcex, sizex);


                    std::vector<double>   cexend;

                    for (int j = 0; j < sizex; j++) {
                        cexend.push_back((double)(s_datcex[j]) / 1000.);
                    }

                    delete [] s_datcex;

                    params.locallab.localTgaincurve.clear();
                    params.locallab.llcurve.clear();
                    params.locallab.LHcurve.clear();
                    params.locallab.cccurve.clear();
                    params.locallab.HHcurve.clear();
                    params.locallab.skintonescurve.clear();
                    params.locallab.excurve.clear();

                    params.locallab.localTgaincurve = cretiend;
                    params.locallab.llcurve = cllend;
                    params.locallab.LHcurve = clhend;
                    params.locallab.cccurve = cccend;
                    params.locallab.HHcurve = chhend;
                    params.locallab.skintonescurve = cskend;
                    params.locallab.excurve = cexend;

                    bool LHutili = false;
                    bool HHutili = false;
                    std::string t_curvhhref = "1000A0B500C350D350E166F500G350H350I333J500K350L350M500N500O350P350Q666R500S350T350U833V500W350X350Y@";

                    if (lhstrs[sp].c_str() != t_curvhhref) {
                        LHutili = true;
                    }

                    if (hhstrs[sp].c_str() != t_curvhhref) {
                        HHutili = true;
                    }


                    params.locallab.getCurves(locRETgainCurve, locRETgainCurverab, loclhCurve, lochhCurve, LHutili, HHutili);
                    bool locallutili = false;
                    bool localcutili = false;
                    bool localskutili = false;
                    bool localexutili = false;
                    std::string t_curvskinref = "3000A0B0C1000D1000E@";
                    std::string t_none = "0A@";

                    if (skinstrs[sp].c_str() != t_curvskinref  && skinstrs[sp].c_str() != t_none) {
                        localskutili = true;
                    }

                    std::string t_curvexref = "3000A0B0C1000D1000E@";

                    if (exstrs[sp].c_str() != t_curvexref  && exstrs[sp].c_str() != t_none) {
                        localexutili = true;
                    }

                    CurveFactory::curveLocal(locallutili, params.locallab.llcurve, lllocalcurve, 1);
                    CurveFactory::curveCCLocal(localcutili, params.locallab.cccurve, cclocalcurve, 1);
                    CurveFactory::curveskLocal(localskutili, params.locallab.skintonescurve, sklocalcurve, 1);

                    CurveFactory::curveexLocal(localexutili, params.locallab.excurve, exlocalcurve, 1);
                    //provisory
                    double br = 0.;
                    double contr = 0.;
                    double ecomp = params.locallab.expcomp;
                    double black = params.locallab.black;
                    double hlcompr = params.locallab.hlcompr;
                    double hlcomprthresh = params.locallab.hlcomprthresh;
                    double shcompr = params.locallab.shcompr;

                    CurveFactory::complexCurvelocal(ecomp, black / 65535., hlcompr, hlcomprthresh, shcompr, br, contr,
                                                    hist16, hltonecurveloc, shtonecurveloc, tonecurveloc,
                                                    1);

                    double huere, chromare, lumare, huerefblu;
                    double sobelre;

                    ipf.calc_ref(labView, labView, 0, 0, fw, fh, 1, huerefblu, huere, chromare, lumare, sobelre);

                    params.locallab.huerefblur = huerefblu;
                    params.locallab.hueref = huere;
                    params.locallab.chromaref = chromare;
                    params.locallab.lumaref = lumare;
                    params.locallab.sobelref = sobelre;

                    ipf.Lab_Local(2, (float**)shbuffer, labView, labView, reservView, 0, 0, fw, fh,  1, locRETgainCurve, lllocalcurve, loclhCurve, lochhCurve,
                                  LHutili, HHutili, cclocalcurve, localskutili, sklocalcurve, localexutili, exlocalcurve, hltonecurveloc, shtonecurveloc, tonecurveloc, params.locallab.huerefblur, params.locallab.hueref, params.locallab.chromaref, params.locallab.lumaref, params.locallab.sobelref);
                    lllocalcurve.clear();
                    cclocalcurve.clear();
                    sklocalcurve.clear();
                    exlocalcurve.clear();

                }



                for (int i = 0; i < maxdata; i++) {
                    delete [] dataspots[i];
                }

                delete [] dataspots;




                delete [] retistrs;
                delete [] llstrs;
                delete [] lhstrs;
                delete [] ccstrs;
                delete [] hhstrs;
                delete [] skinstrs;
                delete [] exstrs;

                if (params.locallab.inverssha) {

                    for (int i = 0; i < fh; i++) {
                        delete [] shbuffer[i];
                    }

                    delete [] shbuffer;
                }


            }

            t2.set();

            if (settings->verbose) {
                printf("Total local:- %d usec\n", t2.etime(t1));
            }

        }

        delete reservView;
        reservView = nullptr;

        ipf.chromiLuminanceCurve(nullptr, 1, labView, labView, curve1, curve2, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, dummy, dummy);

        if ((params.colorappearance.enabled && !params.colorappearance.tonecie) || (!params.colorappearance.enabled)) {
            ipf.EPDToneMap(labView, 5, 1);
        }


        ipf.vibrance(labView);

        if ((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
            ipf.impulsedenoise(labView);
        }

        // for all treatments Defringe, Sharpening, Contrast detail ,Microcontrast they are activated if "CIECAM" function are disabled

        if ((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
            ipf.defringe(labView);
        }

        if (params.sharpenEdge.enabled) {
            ipf.MLsharpen(labView);
        }

        if (params.sharpenMicro.enabled) {
            if ((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
                ipf.MLmicrocontrast(labView);     //!params.colorappearance.sharpcie
            }
        }

        if (((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) && params.sharpening.enabled) {

            float **buffer = new float*[fh];

            for (int i = 0; i < fh; i++) {
                buffer[i] = new float[fw];
            }

            ipf.sharpening(labView, (float**)buffer, params.sharpening);

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

        params.wavelet.getCurves(wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL);


        // directional pyramid wavelet
        if (params.dirpyrequalizer.cbdlMethod == "aft") {
            if ((params.colorappearance.enabled && !settings->autocielab)  || !params.colorappearance.enabled) {
                ipf.dirpyrequalizer(labView, 1);     //TODO: this is the luminance tonecurve, not the RGB one
            }
        }

        bool wavcontlutili = false;

        CurveFactory::curveWavContL(wavcontlutili, params.wavelet.wavclCurve, wavclCurve,/* hist16C, dummy,*/ 1);

        if (params.wavelet.enabled) {
            ipf.ip_wavelet(labView, labView, 2, WaveParams, wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW,  waOpacityCurveWL, wavclCurve, 1);
        }

        wavCLVCurve.Reset();

        //Colorappearance and tone-mapping associated

        int f_w = 1, f_h = 1;

        if (params.colorappearance.tonecie || params.colorappearance.enabled) {
            f_w = fw;
            f_h = fh;
        }

        CieImage *cieView = new CieImage(f_w, (f_h));

        CurveFactory::curveLightBrightColor(
            params.colorappearance.curve,
            params.colorappearance.curve2,
            params.colorappearance.curve3,
            hist16, dummy,
            dummy, dummy,
            customColCurve1,
            customColCurve2,
            customColCurve3,
            1);

        if (params.colorappearance.enabled) {
            double adap;
            int imgNum = 0;

            if (imgsrc->getSensorType() == ST_BAYER) {
                imgNum = params.raw.bayersensor.imageNum;
            } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                //imgNum = params.raw.xtranssensor.imageNum;
            }

            float fnum = imgsrc->getMetaData()->getFNumber(imgNum);          // F number
            float fiso = imgsrc->getMetaData()->getISOSpeed(imgNum) ;        // ISO
            float fspeed = imgsrc->getMetaData()->getShutterSpeed(imgNum) ;  //speed
            float fcomp = imgsrc->getMetaData()->getExpComp(imgNum);         //compensation + -

            if (fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) {
                adap = 2000.;
            }//if no exif data or wrong
            else {
                float E_V = fcomp + log2((fnum * fnum) / fspeed / (fiso / 100.f));
                E_V += params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
                E_V += log2(params.raw.expos);  // exposure raw white point ; log2 ==> linear to EV
                adap = powf(2.f, E_V - 3.f);  //cd / m2
            }

            LUTf CAMBrightCurveJ;
            LUTf CAMBrightCurveQ;
            float CAMMean = NAN;

            if (params.sharpening.enabled) {
                if (settings->ciecamfloat) {
                    float d, dj, yb;
                    ipf.ciecam_02float(cieView, float (adap), 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, d, dj, yb, 1);
                } else {
                    double dd, dj;
                    ipf.ciecam_02(cieView, adap, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, dd, dj, 1);
                }
            } else {
                if (settings->ciecamfloat) {
                    float d, dj, yb;
                    ipf.ciecam_02float(cieView, float (adap), 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, d, dj, yb, 1);
                } else {
                    double dd, dj;
                    ipf.ciecam_02(cieView, adap, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, dd, dj, 1);
                }
            }
        }

        delete cieView;
        cieView = nullptr;




        // end tile processing...???
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (pl) {
            pl->setProgress(0.60);
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

            if (labResize) { // crop lab data
                tmplab = new LabImage(cw, ch);

                for (int row = 0; row < ch; row++) {
                    for (int col = 0; col < cw; col++) {
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
            ipf.Lanczos(labView, tmplab, tmpScale);
            delete labView;
            labView = tmplab;
            cw = labView->W;
            ch = labView->H;

            if (params.prsharpening.enabled) {
                for (int i = 0; i < ch; i++)
                    for (int j = 0; j < cw; j++) {
                        labView->L[i][j] = labView->L[i][j] < 0.f ? 0.f : labView->L[i][j];
                    }

                float **buffer = new float*[ch];

                for (int i = 0; i < ch; i++) {
                    buffer[i] = new float[cw];
                }

                ipf.sharpening(labView, (float**)buffer, params.prsharpening);

                for (int i = 0; i < ch; i++) {
                    delete [] buffer[i];
                }

                delete [] buffer;
            }
        }

        Imagefloat* readyImg = nullptr;
        cmsHPROFILE jprof = nullptr;
        bool customGamma = false;
        bool useLCMS = false;
        bool bwonly = params.blackwhite.enabled && !params.colorToning.enabled && !autili && !butili ;

        if (params.icm.gamma != "default" || params.icm.freegamma) { // if select gamma output between BT709, sRGB, linear, low, high, 2.2 , 1.8

            GammaValues ga;
            //  if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;
            readyImg = ipf.lab2rgbOut (labView, cx, cy, cw, ch, params.icm, &ga);
            customGamma = true;

            //or selected Free gamma
            useLCMS = false;

            if ((jprof = ICCStore::getInstance()->createCustomGammaOutputProfile(params.icm, ga)) == nullptr) {
                useLCMS = true;
            }

        } else {
            // if Default gamma mode: we use the profile selected in the "Output profile" combobox;
            // gamma come from the selected profile, otherwise it comes from "Free gamma" tool

            readyImg = ipf.lab2rgbOut (labView, cx, cy, cw, ch, params.icm);

            if (settings->verbose) {
                printf("Output profile_: \"%s\"\n", params.icm.output.c_str());
            }
        }

        delete labView;
        labView = nullptr;

//       delete reservView;
//       reservView = nullptr;


        if (bwonly) { //force BW r=g=b
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

        if (pl) {
            pl->setProgress(0.70);
        }

        if (tmpScale != 1.0 && params.resize.method == "Nearest") { // resize rgb data (gamma applied)
            Imagefloat* tempImage = new Imagefloat (imw, imh);
            ipf.resize(readyImg, tempImage, tmpScale);
            delete readyImg;
            readyImg = tempImage;
        }

        switch (params.metadata.mode) {
            case MetaDataParams::TUNNEL:
                // Sending back the whole first root, which won't necessarily be the selected frame number
                // and may contain subframe depending on initial raw's hierarchy
                readyImg->setMetadata(ii->getMetaData()->getRootExifData());
                break;

            case MetaDataParams::EDIT:
                // ask for the correct frame number, but may contain subframe depending on initial raw's hierarchy
                readyImg->setMetadata(ii->getMetaData()->getBestExifData(imgsrc, &params.raw), params.exif, params.iptc);
                break;

            default: // case MetaDataParams::STRIP
                // nothing to do
                break;
        }


        // Setting the output curve to readyImg
        if (customGamma) {
            if (!useLCMS) {
                // use corrected sRGB profile in order to apply a good TRC if present, otherwise use LCMS2 profile generated by lab2rgb16 w/ gamma
                ProfileContent pc(jprof);
                readyImg->setOutputProfile(pc.getData().c_str(), pc.getData().size());
            }
        } else {
            // use the selected output profile if present, otherwise use LCMS2 profile generate by lab2rgb16 w/ gamma

            if (params.icm.output != "" && params.icm.output != ColorManagementParams::NoICMString) {

                // if ICCStore::getInstance()->getProfile send back an object, then ICCStore::getInstance()->getContent will do too
                cmsHPROFILE jprof = ICCStore::getInstance()->getProfile(params.icm.output);  //get outProfile

                if (jprof == nullptr) {
                    if (settings->verbose) {
                        printf("\"%s\" ICC output profile not found!\n - use LCMS2 substitution\n", params.icm.output.c_str());
                    }
                } else {
                    if (settings->verbose) {
                        printf("Using \"%s\" output profile\n", params.icm.output.c_str());
                    }

                    ProfileContent pc = ICCStore::getInstance()->getContent(params.icm.output);
                    readyImg->setOutputProfile(pc.getData().c_str(), pc.getData().size());
                }
            } else {
                // No ICM
                readyImg->setOutputProfile(nullptr, 0);
            }
        }

//    t2.set();
//    if( settings->verbose )
//           printf("Total:- %d usec\n", t2.etime(t1));

        if (!job->initialImage) {
            ii->decreaseRef();
        }

        delete job;

        if (pl) {
            pl->setProgress(0.75);
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

    void stage_early_resize()
    {
        procparams::ProcParams& params = job->pparams;
        //ImProcFunctions ipf (&params, true);
        ImProcFunctions &ipf = * (ipf_p.get());

        int imw, imh;
        double scale_factor = ipf.resizeScale(&params, fw, fh, imw, imh);

        std::unique_ptr<LabImage> tmplab(new LabImage(fw, fh));
        ipf.rgb2lab(*baseImg, *tmplab, params.icm.working);

        if (params.crop.enabled) {
            int cx = params.crop.x;
            int cy = params.crop.y;
            int cw = params.crop.w;
            int ch = params.crop.h;

            std::unique_ptr<LabImage> cropped(new LabImage(cw, ch));

            for (int row = 0; row < ch; row++) {
                for (int col = 0; col < cw; col++) {
                    cropped->L[row][col] = tmplab->L[row + cy][col + cx];
                    cropped->a[row][col] = tmplab->a[row + cy][col + cx];
                    cropped->b[row][col] = tmplab->b[row + cy][col + cx];
                }
            }

            tmplab = std::move(cropped);
        }

        assert(params.resize.enabled);

        // resize image
        {
            std::unique_ptr<LabImage> resized(new LabImage(imw, imh));
            ipf.Lanczos(tmplab.get(), resized.get(), scale_factor);
            tmplab = std::move(resized);
        }

        adjust_procparams(scale_factor);

        fw = imw;
        fh = imh;

        delete baseImg;
        baseImg = new Imagefloat(fw, fh);
        ipf.lab2rgb(*tmplab, *baseImg, params.icm.working);
    }

    void adjust_procparams(double scale_factor)
    {
        procparams::ProcParams &params = job->pparams;
        procparams::ProcParams defaultparams;

        params.resize.enabled = false;
        params.crop.enabled = false;

        if (params.prsharpening.enabled) {
            params.sharpening = params.prsharpening;
        } else {
            params.sharpening.radius *= scale_factor;
        }

        params.impulseDenoise.thresh *= scale_factor;

        if (scale_factor < 0.5) {
            params.impulseDenoise.enabled = false;
        }

        params.wavelet.strength *= scale_factor;
        params.dirpyrDenoise.luma *= scale_factor * scale_factor;
        //params.dirpyrDenoise.Ldetail += (100 - params.dirpyrDenoise.Ldetail) * scale_factor;
        auto &lcurve = params.dirpyrDenoise.lcurve;

        for (size_t i = 2; i < lcurve.size(); i += 4) {
            lcurve[i] *= min (scale_factor * scale_factor, 1.0);
        }

        noiseLCurve.Set(lcurve);
        const char *medmethods[] = { "soft", "33", "55soft", "55", "77", "99" };

        if (params.dirpyrDenoise.median) {
            auto &key = params.dirpyrDenoise.methodmed == "RGB" ? params.dirpyrDenoise.rgbmethod : params.dirpyrDenoise.medmethod;

            for (int i = 1; i < int (sizeof(medmethods) / sizeof(const char *)); ++i) {
                if (key == medmethods[i]) {
                    int j = i - int (1.0 / scale_factor);

                    if (j < 0) {
                        params.dirpyrDenoise.median = false;
                    } else {
                        key = medmethods[j];
                    }

                    break;
                }
            }
        }

        params.epd.scale *= scale_factor;
        //params.epd.edgeStopping *= scale_factor;

        const double dirpyreq_scale = min(scale_factor * 1.5, 1.0);

        for (int i = 0; i < 6; ++i) {
            adjust_radius(defaultparams.dirpyrequalizer.mult[i], dirpyreq_scale,
                          params.dirpyrequalizer.mult[i]);
        }

        params.dirpyrequalizer.threshold *= scale_factor;

        adjust_radius(defaultparams.defringe.radius, scale_factor,
                      params.defringe.radius);
        params.sh.radius *= scale_factor;
        params.localContrast.radius *= scale_factor;

        if (params.raw.xtranssensor.method == procparams::RAWParams::XTransSensor::getMethodString(procparams::RAWParams::XTransSensor::Method::THREE_PASS)) {
            params.raw.xtranssensor.method = procparams::RAWParams::XTransSensor::getMethodString(procparams::RAWParams::XTransSensor::Method::ONE_PASS);
        }

        if (params.raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::PIXELSHIFT)) {
            params.raw.bayersensor.method = procparams::RAWParams::BayerSensor::getMethodString(params.raw.bayersensor.pixelShiftLmmse ? procparams::RAWParams::BayerSensor::Method::LMMSE : procparams::RAWParams::BayerSensor::Method::AMAZE);
        }
    }

private:
    ProcessingJobImpl* job;
    int& errorCode;
    ProgressListener* pl;
    bool flush;

    // internal state
    std::unique_ptr<ImProcFunctions> ipf_p;
    InitialImage *ii;
    ImageSource *imgsrc;
    int fw;
    int fh;

    int tr;
    PreviewProps pp;

    NoiseCurve noiseLCurve;
    NoiseCurve noiseCCurve;
    Imagefloat *calclum;
    float autoNR;
    float autoNRmax;
    int tilesize;
    int overlap;

    float *ch_M;
    float *max_r;
    float *max_b;
    float *min_b;
    float *min_r;
    float *lumL;
    float *chromC;
    float *ry;
    float *sk;
    float *pcsk;

    double expcomp;
    int bright;
    int contr;
    int black;
    int hlcompr;
    int hlcomprthresh;

    ColorTemp currWB;
    Imagefloat *baseImg;
    LabImage* labView;
    LabImage* reservView;

    LUTu hist16;

    LUTf curve1;
    LUTf curve2;
    LUTf curve;
    LUTf satcurve;
    LUTf lhskcurve;
    LUTf lumacurve;
    LUTf clcurve;
    LUTf clToningcurve;
    LUTf cl2Toningcurve;
    LUTf wavclCurve;

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

    bool autili, butili;
};

} // namespace


IImagefloat* processImage (ProcessingJob* pjob, int& errorCode, ProgressListener* pl, bool flush)
{
    ImageProcessor proc(pjob, errorCode, pl, flush);
    return proc();
}

void batchProcessingThread(ProcessingJob* job, BatchProcessingListener* bpl)
{

    ProcessingJob* currentJob = job;

    while (currentJob) {
        int errorCode;
        IImagefloat* img = processImage (currentJob, errorCode, bpl, true);

        if (errorCode) {
            bpl->error(M("MAIN_MSG_CANNOTLOAD"));
            currentJob = nullptr;
        } else {
            try {
                currentJob = bpl->imageReady(img);
            } catch (Glib::Exception& ex) {
                bpl->error(ex.what());
                currentJob = nullptr;
            }
        }
    }
}

void startBatchProcessing(ProcessingJob* job, BatchProcessingListener* bpl)
{

    if (bpl) {
        Glib::Thread::create(sigc::bind(sigc::ptr_fun(batchProcessingThread), job, bpl), 0, true, true, Glib::THREAD_PRIORITY_LOW);
    }

}

}
