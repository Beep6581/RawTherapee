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

namespace
{

template <typename T>
void adjust_radius (const T &default_param, double scale_factor, T &param)
{
    const double delta = (param - default_param) * scale_factor;
    param = default_param + delta;
}


class ImageProcessor
{
public:
    ImageProcessor (ProcessingJob* pjob, int& errorCode,
                    ProgressListener* pl, bool tunnelMetaData, bool flush):
        job (static_cast<ProcessingJobImpl*> (pjob)),
        errorCode (errorCode),
        pl (pl),
        tunnelMetaData (tunnelMetaData),
        flush (flush),
        // internal state
        ipf_p (nullptr),
        ii (nullptr),
        imgsrc (nullptr),
        fw (-1),
        fh (-1),
        pp (0, 0, 0, 0, 0)
    {
    }

    Image16 *operator()()
    {
        if (!job->fast) {
            return normal_pipeline();
        } else {
            return fast_pipeline();
        }
    }

private:
    Image16 *normal_pipeline()
    {
        if (!stage_init()) {
            return nullptr;
        }

        stage_denoise();
        stage_transform();
        return stage_finish();
    }

    Image16 *fast_pipeline()
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
            pl->setProgressStr ("PROGRESSBAR_PROCESSING");
            pl->setProgress (0.0);
        }

        ii = job->initialImage;

        if (!ii) {
            ii = InitialImage::load (job->fname, job->isRaw, &errorCode);

            if (errorCode) {
                delete job;
                return false; //return nullptr;
            }
        }

        procparams::ProcParams& params = job->pparams;

        // acquire image from imagesource
        imgsrc = ii->getImageSource ();

        tr = getCoarseBitMask (params.coarse);
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

        ipf_p.reset (new ImProcFunctions (&params, true));
        ImProcFunctions &ipf = * (ipf_p.get());

        pp = PreviewProps (0, 0, fw, fh, 1);
        imgsrc->setCurrentFrame (params.raw.bayersensor.imageNum);
        imgsrc->preprocess ( params.raw, params.lensProf, params.coarse, params.dirpyrDenoise.enabled);

        if (params.toneCurve.autoexp) {// this enabled HLRecovery
            LUTu histRedRaw (256), histGreenRaw (256), histBlueRaw (256);
            imgsrc->getRAWHistogram (histRedRaw, histGreenRaw, histBlueRaw);

            if (ToneCurveParams::HLReconstructionNecessary (histRedRaw, histGreenRaw, histBlueRaw) && !params.toneCurve.hrenabled) {
                params.toneCurve.hrenabled = true;
                // WARNING: Highlight Reconstruction is being forced 'on', should we force a method here too?
            }
        }

        if (pl) {
            pl->setProgress (0.20);
        }

        imgsrc->demosaic ( params.raw);

        if (pl) {
            pl->setProgress (0.30);
        }

        if (params.retinex.enabled) { //enabled Retinex
            LUTf cdcurve (65536, 0);
            LUTf mapcurve (65536, 0);
            LUTu dummy;
            RetinextransmissionCurve dehatransmissionCurve;
            RetinexgaintransmissionCurve dehagaintransmissionCurve;
            bool dehacontlutili = false;
            bool mapcontlutili = false;
            bool useHsl = false;
//        multi_array2D<float, 3> conversionBuffer(1, 1);
            multi_array2D<float, 4> conversionBuffer (1, 1);
            imgsrc->retinexPrepareBuffers (params.icm, params.retinex, conversionBuffer, dummy);
            imgsrc->retinexPrepareCurves (params.retinex, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, dehacontlutili, mapcontlutili, useHsl, dummy, dummy );
            float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            imgsrc->retinex ( params.icm, params.retinex, params.toneCurve, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, conversionBuffer, dehacontlutili, mapcontlutili, useHsl, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax, dummy);
        }

        if (pl) {
            pl->setProgress (0.40);
        }

        imgsrc->HLRecovery_Global ( params.toneCurve );


        if (pl) {
            pl->setProgress (0.45);
        }

        // set the color temperature
        currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

        if (params.wb.method == "Camera") {
            currWB = imgsrc->getWB ();
        } else if (params.wb.method == "Auto") {
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers (rm, gm, bm);
            currWB.update (rm, gm, bm, params.wb.equal, params.wb.tempBias);
        }

        calclum = nullptr ;
        params.dirpyrDenoise.getCurves (noiseLCurve, noiseCCurve);
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
        ipf.Tile_calc (tilesize, overlap, 2, fw, fh, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);
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
                LUTf gamcurve (65536, 0);
                float gam, gamthresh, gamslope;
                ipf.RGB_denoise_infoGamCurve (params.dirpyrDenoise, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope);
                #pragma omp parallel
                {
                    Imagefloat *origCropPart;//init auto noise
                    origCropPart = new Imagefloat (crW, crH);//allocate memory
                    Imagefloat *provicalc = new Imagefloat ((crW + 1) / 2, (crH + 1) / 2); //for denoise curves
                    int skipP = 1;
                    #pragma omp for schedule(dynamic) collapse(2) nowait

                    for (int wcr = 0; wcr < numtiles_W; wcr++) {
                        for (int hcr = 0; hcr < numtiles_H; hcr++) {
                            int beg_tileW = wcr * tileWskip + tileWskip / 2.f - crW / 2.f;
                            int beg_tileH = hcr * tileHskip + tileHskip / 2.f - crH / 2.f;
                            PreviewProps ppP (beg_tileW, beg_tileH, crW, crH, skipP);
                            imgsrc->getImage (currWB, tr, origCropPart, ppP, params.toneCurve, params.icm, params.raw );
                            //baseImg->getStdImage(currWB, tr, origCropPart, ppP, true, params.toneCurve);

                            // we only need image reduced to 1/4 here
                            for (int ii = 0; ii < crH; ii += 2) {
                                for (int jj = 0; jj < crW; jj += 2) {
                                    provicalc->r (ii >> 1, jj >> 1) = origCropPart->r (ii, jj);
                                    provicalc->g (ii >> 1, jj >> 1) = origCropPart->g (ii, jj);
                                    provicalc->b (ii >> 1, jj >> 1) = origCropPart->b (ii, jj);
                                }
                            }

                            imgsrc->convertColorSpace (provicalc, params.icm, currWB); //for denoise luminance curve
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
                            ipf.RGB_denoise_info (origCropPart, provicalc, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope, params.dirpyrDenoise, imgsrc->getDirPyrDenoiseExpComp(), chaut, Nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc);
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

                            if (!imgsrc->isRAW()) {
                                multip = 2.f;    //take into account gamma for TIF / JPG approximate value...not good fot gamma=1
                            }

                            float maxmax = max (maxredaut, maxblueaut);
                            float delta;
                            int mode = 2;
                            int lissage = settings->leveldnliss;
                            ipf.calcautodn_info (chaut, delta, Nb, levaut, maxmax, lumema, chromina, mode, lissage, redyel, skinc, nsknc);

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
                    printf ("Info denoise ponderated performed in %d usec:\n", t2pone.etime (t1pone));
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
                LUTf gamcurve (65536, 0);
                float gam, gamthresh, gamslope;
                ipf.RGB_denoise_infoGamCurve (params.dirpyrDenoise, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope);
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

                    for (int wcr = 0; wcr <= 2; wcr++) {
                        for (int hcr = 0; hcr <= 2; hcr++) {
                            PreviewProps ppP (coordW[wcr], coordH[hcr], crW, crH, 1);
                            imgsrc->getImage (currWB, tr, origCropPart, ppP, params.toneCurve, params.icm, params.raw);
                            //baseImg->getStdImage(currWB, tr, origCropPart, ppP, true, params.toneCurve);


                            // we only need image reduced to 1/4 here
                            for (int ii = 0; ii < crH; ii += 2) {
                                for (int jj = 0; jj < crW; jj += 2) {
                                    provicalc->r (ii >> 1, jj >> 1) = origCropPart->r (ii, jj);
                                    provicalc->g (ii >> 1, jj >> 1) = origCropPart->g (ii, jj);
                                    provicalc->b (ii >> 1, jj >> 1) = origCropPart->b (ii, jj);
                                }
                            }

                            imgsrc->convertColorSpace (provicalc, params.icm, currWB); //for denoise luminance curve
                            int nb = 0;
                            float chaut = 0.f, redaut = 0.f, blueaut = 0.f, maxredaut = 0.f, maxblueaut = 0.f, minredaut = 0.f, minblueaut = 0.f, chromina = 0.f, sigma = 0.f, lumema = 0.f, sigma_L = 0.f, redyel = 0.f, skinc = 0.f, nsknc = 0.f;
                            ipf.RGB_denoise_info (origCropPart, provicalc, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope,  params.dirpyrDenoise, imgsrc->getDirPyrDenoiseExpComp(), chaut, nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc);
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

                if (!imgsrc->isRAW()) {
                    multip = 2.f;    //take into account gamma for TIF / JPG approximate value...not good fot gamma=1
                }

                float delta[9];
                int mode = 1;
                int lissage = settings->leveldnliss;

                for (int k = 0; k < 9; k++) {
                    float maxmax = max (max_r[k], max_b[k]);
                    ipf.calcautodn_info (ch_M[k], delta[k], Nb[k], levaut, maxmax, lumL[k], chromC[k], mode, lissage, ry[k], sk[k], pcsk[k] );
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
                printf ("Info denoise auto performed in %d usec:\n", t2aue.etime (t1aue));
            }

            //end evaluate noise
        }

        baseImg = new Imagefloat (fw, fh);
        imgsrc->getImage (currWB, tr, baseImg, pp, params.toneCurve, params.icm, params.raw);

        if (pl) {
            pl->setProgress (0.50);
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
            imgsrc->getAutoExpHistogram (aehist, aehistcompr);
            ipf.getAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), params.toneCurve.clip, expcomp, bright, contr, black, hlcompr, hlcomprthresh);
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

        if (denoiseParams.enabled  && (noiseLCurve || noiseCCurve )) {
            // we only need image reduced to 1/4 here
            calclum = new Imagefloat ((fw + 1) / 2, (fh + 1) / 2); //for luminance denoise curve
            #pragma omp parallel for

            for (int ii = 0; ii < fh; ii += 2) {
                for (int jj = 0; jj < fw; jj += 2) {
                    calclum->r (ii >> 1, jj >> 1) = baseImg->r (ii, jj);
                    calclum->g (ii >> 1, jj >> 1) = baseImg->g (ii, jj);
                    calclum->b (ii >> 1, jj >> 1) = baseImg->b (ii, jj);
                }
            }

            imgsrc->convertColorSpace (calclum, params.icm, currWB);
        }

        if (denoiseParams.enabled) {
            // CurveFactory::denoiseLL(lldenoiseutili, denoiseParams.lcurve, Noisecurve,1);
            //denoiseParams.getCurves(noiseLCurve);
//      ipf.RGB_denoise(baseImg, baseImg, calclum, imgsrc->isRAW(), denoiseParams, params.defringe, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, lldenoiseutili);
            float chaut, redaut, blueaut, maxredaut, maxblueaut, nresi, highresi;
            int kall = 2;
            ipf.RGB_denoise (kall, baseImg, baseImg, calclum, ch_M, max_r, max_b, imgsrc->isRAW(), denoiseParams, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, noiseCCurve, chaut, redaut, blueaut, maxredaut, maxblueaut, nresi, highresi);

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

        imgsrc->convertColorSpace (baseImg, params.icm, currWB);

        // perform first analysis
        hist16 (65536);

        ipf.firstAnalysis (baseImg, params, hist16);

        // perform transform (excepted resizing)
        if (ipf.needsTransform()) {
            Imagefloat* trImg = nullptr;
            if (ipf.needsLuminanceOnly()) {
                trImg = baseImg;
            } else {
                trImg = new Imagefloat (fw, fh);
            }
            ipf.transform (baseImg, trImg, 0, 0, 0, 0, fw, fh, fw, fh,
                           imgsrc->getMetaData(), imgsrc->getRotateDegree(), true);
            if(trImg != baseImg) {
                delete baseImg;
                baseImg = trImg;
            }
        }
    }

    Image16 *stage_finish()
    {
        procparams::ProcParams& params = job->pparams;
        //ImProcFunctions ipf (&params, true);
        ImProcFunctions &ipf = * (ipf_p.get());

        if (params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled) {
            const int W = baseImg->getWidth();
            const int H = baseImg->getHeight();
            LabImage labcbdl (W, H);
            ipf.rgb2lab (*baseImg, labcbdl, params.icm.working);
            ipf.dirpyrequalizer (&labcbdl, 1);
            ipf.lab2rgb (labcbdl, *baseImg, params.icm.working);
        }

        // update blurmap
        SHMap* shmap = nullptr;

        if (params.sh.enabled) {
            shmap = new SHMap (fw, fh, true);
            double radius = sqrt (double (fw * fw + fh * fh)) / 2.0;
            double shradius = params.sh.radius;

            if (!params.sh.hq) {
                shradius *= radius / 1800.0;
            }

            shmap->update (baseImg, shradius, ipf.lumimul, params.sh.hq, 1);
        }

        // RGB processing

        curve1 (65536);
        curve2 (65536);
        curve (65536, 0);
        satcurve (65536, 0);
        lhskcurve (65536, 0);
        lumacurve (32770, 0); // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
        clcurve (65536, 0);
        wavclCurve (65536, 0);

        //if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;

        CurveFactory::complexCurve (expcomp, black / 65535.0, hlcompr, hlcomprthresh, params.toneCurve.shcompr, bright, contr,
                                    params.toneCurve.curveMode, params.toneCurve.curve, params.toneCurve.curveMode2, params.toneCurve.curve2,
                                    hist16, curve1, curve2, curve, dummy, customToneCurve1, customToneCurve2 );

        CurveFactory::RGBCurve (params.rgbCurves.rcurve, rCurve, 1);
        CurveFactory::RGBCurve (params.rgbCurves.gcurve, gCurve, 1);
        CurveFactory::RGBCurve (params.rgbCurves.bcurve, bCurve, 1);

        bool opautili = false;

        if (params.colorToning.enabled) {
            TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params.icm.working);
            double wp[3][3] = {
                {wprof[0][0], wprof[0][1], wprof[0][2]},
                {wprof[1][0], wprof[1][1], wprof[1][2]},
                {wprof[2][0], wprof[2][1], wprof[2][2]}
            };
            TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params.icm.working);
            double wip[3][3] = {
                {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
                {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
                {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
            };
            params.colorToning.getCurves (ctColorCurve, ctOpacityCurve, wp, wip, opautili);
            clToningcurve (65536, 0);
            CurveFactory::curveToning (params.colorToning.clcurve, clToningcurve, 1);
            cl2Toningcurve (65536, 0);
            CurveFactory::curveToning (params.colorToning.cl2curve, cl2Toningcurve, 1);
        }

        labView = new LabImage (fw, fh);

        if (params.blackwhite.enabled) {
            CurveFactory::curveBW (params.blackwhite.beforeCurve, params.blackwhite.afterCurve, hist16, dummy, customToneCurvebw1, customToneCurvebw2, 1);
        }

        double rrm, ggm, bbm;
        float autor, autog, autob;
        float satLimit = float (params.colorToning.satProtectionThreshold) / 100.f * 0.7f + 0.3f;
        float satLimitOpacity = 1.f - (float (params.colorToning.saturatedOpacity) / 100.f);

        if (params.colorToning.enabled  && params.colorToning.autosat) { //for colortoning evaluation of saturation settings
            float moyS = 0.f;
            float eqty = 0.f;
            ipf.moyeqt (baseImg, moyS, eqty);//return image : mean saturation and standard dev of saturation
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
        DCPProfile *dcpProf = imgsrc->getDCP (params.icm, currWB, as);

        LUTu histToneCurve;

        ipf.rgbProc (baseImg, labView, nullptr, curve1, curve2, curve, shmap, params.toneCurve.saturation, rCurve, gCurve, bCurve, satLimit, satLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2, customToneCurvebw1, customToneCurvebw2, rrm, ggm, bbm, autor, autog, autob, expcomp, hlcompr, hlcomprthresh, dcpProf, as, histToneCurve);

        if (settings->verbose) {
            printf ("Output image / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", autor, autog, autob);
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
        baseImg = nullptr;

        if (shmap) {
            delete shmap;
        }

        shmap = nullptr;

        if (pl) {
            pl->setProgress (0.55);
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
                LUTu hist16thr (hist16.getSize());  // one temporary lookup table per thread
                hist16thr.clear();
#ifdef _OPENMP
                #pragma omp for schedule(static) nowait
#endif

                for (int i = 0; i < fh; i++)
                    for (int j = 0; j < fw; j++) {
                        hist16thr[ (int) ((labView->L[i][j]))]++;
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
        CurveFactory::curveCL (clcutili, params.labCurve.clcurve, clcurve, 1);

        bool ccutili, cclutili;
        CurveFactory::complexsgnCurve (autili, butili, ccutili, cclutili, params.labCurve.acurve, params.labCurve.bcurve, params.labCurve.cccurve,
                                       params.labCurve.lccurve, curve1, curve2, satcurve, lhskcurve, 1);

        ipf.chromiLuminanceCurve (nullptr, 1, labView, labView, curve1, curve2, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, dummy, dummy);

        if (params.fattal.enabled) {
            ipf.ToneMapFattal02(labView, 3);
        }

        if ((params.colorappearance.enabled && !params.colorappearance.tonecie) || (!params.colorappearance.enabled)) {
            ipf.EPDToneMap (labView, 5, 1);
        }


        ipf.vibrance (labView);

        if ((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
            ipf.impulsedenoise (labView);
        }

        // for all treatments Defringe, Sharpening, Contrast detail ,Microcontrast they are activated if "CIECAM" function are disabled

        if ((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
            ipf.defringe (labView);
        }

        if (params.sharpenEdge.enabled) {
            ipf.MLsharpen (labView);
        }

        if (params.sharpenMicro.enabled) {
            if ((params.colorappearance.enabled && !settings->autocielab) ||  (!params.colorappearance.enabled)) {
                ipf.MLmicrocontrast (labView);    //!params.colorappearance.sharpcie
            }
        }

        if (((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) && params.sharpening.enabled) {

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

        params.wavelet.getCurves (wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL );


        // directional pyramid wavelet
        if (params.dirpyrequalizer.cbdlMethod == "aft") {
            if ((params.colorappearance.enabled && !settings->autocielab)  || !params.colorappearance.enabled) {
                ipf.dirpyrequalizer (labView, 1);    //TODO: this is the luminance tonecurve, not the RGB one
            }
        }

        bool wavcontlutili = false;

        CurveFactory::curveWavContL (wavcontlutili, params.wavelet.wavclCurve, wavclCurve,/* hist16C, dummy,*/ 1);

        if (params.wavelet.enabled) {
            ipf.ip_wavelet (labView, labView, 2, WaveParams, wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW,  waOpacityCurveWL, wavclCurve, wavcontlutili, 1);
        }

        wavCLVCurve.Reset();

        //Colorappearance and tone-mapping associated

        int f_w = 1, f_h = 1;
        int begh = 0, endh = fh;

        if (params.colorappearance.tonecie || params.colorappearance.enabled) {
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

        if (params.colorappearance.enabled) {
            double adap;
            int imgNum = 0;
            if (imgsrc->getSensorType() == ST_BAYER) {
                imgNum = params.raw.bayersensor.imageNum;
            } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                //imgNum = params.raw.xtranssensor.imageNum;
            }
            float fnum = imgsrc->getMetaData()->getFNumber (imgNum);         // F number
            float fiso = imgsrc->getMetaData()->getISOSpeed (imgNum) ;       // ISO
            float fspeed = imgsrc->getMetaData()->getShutterSpeed (imgNum) ; //speed
            float fcomp = imgsrc->getMetaData()->getExpComp (imgNum);        //compensation + -

            if (fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) {
                adap = 2000.;
            }//if no exif data or wrong
            else {
                float E_V = fcomp + log2 ((fnum * fnum) / fspeed / (fiso / 100.f));
                E_V += params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
                E_V += log2 (params.raw.expos); // exposure raw white point ; log2 ==> linear to EV
                adap = powf (2.f, E_V - 3.f); //cd / m2
            }

            LUTf CAMBrightCurveJ;
            LUTf CAMBrightCurveQ;
            float CAMMean = NAN;

            if (params.sharpening.enabled) {
                if (settings->ciecamfloat) {
                    float d, dj, yb;
                    ipf.ciecam_02float (cieView, float (adap), begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, d, dj, yb, 1);
                } else {
                    double dd, dj, yb;
                    ipf.ciecam_02 (cieView, adap, begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, dd, dj, yb, 1);
                }
            } else {
                if (settings->ciecamfloat) {
                    float d, dj, yb;
                    ipf.ciecam_02float (cieView, float (adap), begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, d, dj, yb, 1);
                } else {
                    double dd, dj, yb;
                    ipf.ciecam_02 (cieView, adap, begh, endh, 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, dd, dj, yb, 1);
                }
            }
        }

        delete cieView;
        cieView = nullptr;




        // end tile processing...???
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (pl) {
            pl->setProgress (0.60);
        }

        int imw, imh;
        double tmpScale = ipf.resizeScale (&params, fw, fh, imw, imh);
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
                tmplab = new LabImage (cw, ch);

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
            tmplab = new LabImage (imw, imh);
            ipf.Lanczos (labView, tmplab, tmpScale);
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

                ipf.sharpening (labView, (float**)buffer, params.prsharpening);

                for (int i = 0; i < ch; i++) {
                    delete [] buffer[i];
                }

                delete [] buffer;
            }
        }

        Image16* readyImg = nullptr;
        cmsHPROFILE jprof = nullptr;
        bool customGamma = false;
        bool useLCMS = false;
        bool bwonly = params.blackwhite.enabled && !params.colorToning.enabled && !autili && !butili ;

        if (params.icm.gamma != "default" || params.icm.freegamma) { // if select gamma output between BT709, sRGB, linear, low, high, 2.2 , 1.8

            GammaValues ga;
            //  if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;
            readyImg = ipf.lab2rgb16 (labView, cx, cy, cw, ch, params.icm, bwonly, &ga);
            customGamma = true;

            //or selected Free gamma
            useLCMS = false;

            if ((jprof = ICCStore::getInstance()->createCustomGammaOutputProfile (params.icm, ga)) == nullptr) {
                useLCMS = true;
            }

        } else {
            // if Default gamma mode: we use the profile selected in the "Output profile" combobox;
            // gamma come from the selected profile, otherwise it comes from "Free gamma" tool

            readyImg = ipf.lab2rgb16 (labView, cx, cy, cw, ch, params.icm, bwonly);

            if (settings->verbose) {
                printf ("Output profile_: \"%s\"\n", params.icm.output.c_str());
            }
        }

        delete labView;
        labView = nullptr;



        if (bwonly) { //force BW r=g=b
            if (settings->verbose) {
                printf ("Force BW\n");
            }

            for (int ccw = 0; ccw < cw; ccw++) {
                for (int cch = 0; cch < ch; cch++) {
                    readyImg->r (cch, ccw) = readyImg->g (cch, ccw);
                    readyImg->b (cch, ccw) = readyImg->g (cch, ccw);
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
            // Sending back the whole first root, which won't necessarily be the selected frame number
            // and may contain subframe depending on initial raw's hierarchy
            readyImg->setMetadata (ii->getMetaData()->getRootExifData ());
        } else {
            // ask for the correct frame number, but may contain subframe depending on initial raw's hierarchy
            readyImg->setMetadata (ii->getMetaData()->getBestExifData(imgsrc, &params.raw), params.exif, params.iptc);
        }


        // Setting the output curve to readyImg
        if (customGamma) {
            if (!useLCMS) {
                // use corrected sRGB profile in order to apply a good TRC if present, otherwise use LCMS2 profile generated by lab2rgb16 w/ gamma
                ProfileContent pc (jprof);
                readyImg->setOutputProfile (pc.getData().c_str(), pc.getData().size());
            }
        } else {
            // use the selected output profile if present, otherwise use LCMS2 profile generate by lab2rgb16 w/ gamma

            if (params.icm.output != "" && params.icm.output != ColorManagementParams::NoICMString) {

                // if ICCStore::getInstance()->getProfile send back an object, then ICCStore::getInstance()->getContent will do too
                cmsHPROFILE jprof = ICCStore::getInstance()->getProfile (params.icm.output); //get outProfile

                if (jprof == nullptr) {
                    if (settings->verbose) {
                        printf ("\"%s\" ICC output profile not found!\n - use LCMS2 substitution\n", params.icm.output.c_str());
                    }
                } else {
                    if (settings->verbose) {
                        printf ("Using \"%s\" output profile\n", params.icm.output.c_str());
                    }

                    ProfileContent pc = ICCStore::getInstance()->getContent (params.icm.output);
                    readyImg->setOutputProfile (pc.getData().c_str(), pc.getData().size());
                }
            } else {
                // No ICM
                readyImg->setOutputProfile (nullptr, 0);
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

    void stage_early_resize()
    {
        procparams::ProcParams& params = job->pparams;
        //ImProcFunctions ipf (&params, true);
        ImProcFunctions &ipf = * (ipf_p.get());

        int imw, imh;
        double scale_factor = ipf.resizeScale (&params, fw, fh, imw, imh);

        std::unique_ptr<LabImage> tmplab (new LabImage (fw, fh));
        ipf.rgb2lab (*baseImg, *tmplab, params.icm.working);

        if (params.crop.enabled) {
            int cx = params.crop.x;
            int cy = params.crop.y;
            int cw = params.crop.w;
            int ch = params.crop.h;

            std::unique_ptr<LabImage> cropped (new LabImage (cw, ch));

            for (int row = 0; row < ch; row++) {
                for (int col = 0; col < cw; col++) {
                    cropped->L[row][col] = tmplab->L[row + cy][col + cx];
                    cropped->a[row][col] = tmplab->a[row + cy][col + cx];
                    cropped->b[row][col] = tmplab->b[row + cy][col + cx];
                }
            }

            tmplab = std::move (cropped);
        }

        assert (params.resize.enabled);

        // resize image
        {
            std::unique_ptr<LabImage> resized (new LabImage (imw, imh));
            ipf.Lanczos (tmplab.get(), resized.get(), scale_factor);
            tmplab = std::move (resized);
        }

        adjust_procparams (scale_factor);

        fw = imw;
        fh = imh;

        delete baseImg;
        baseImg = new Imagefloat (fw, fh);
        ipf.lab2rgb (*tmplab, *baseImg, params.icm.working);
    }

    void adjust_procparams (double scale_factor)
    {
        procparams::ProcParams &params = job->pparams;
        procparams::ProcParams defaultparams;

        params.resize.enabled = false;
        params.crop.enabled = false;

        if (params.prsharpening.enabled) {
            params.sharpening = params.prsharpening;
        } else {
            adjust_radius (defaultparams.sharpening.radius, scale_factor,
                           params.sharpening.radius);
        }

        params.impulseDenoise.thresh *= scale_factor;

        if (scale_factor < 0.5) {
            params.impulseDenoise.enabled = false;
        }

        params.wavelet.strength *= scale_factor;
        params.dirpyrDenoise.luma *= scale_factor;
        //params.dirpyrDenoise.Ldetail += (100 - params.dirpyrDenoise.Ldetail) * scale_factor;
        auto &lcurve = params.dirpyrDenoise.lcurve;

        for (size_t i = 2; i < lcurve.size(); i += 4) {
            lcurve[i] *= min (scale_factor * 2, 1.0);
        }

        noiseLCurve.Set (lcurve);
        const char *medmethods[] = { "soft", "33", "55soft", "55", "77", "99" };

        if (params.dirpyrDenoise.median) {
            auto &key = params.dirpyrDenoise.methodmed == "RGB" ? params.dirpyrDenoise.rgbmethod : params.dirpyrDenoise.medmethod;

            for (int i = 1; i < int (sizeof (medmethods) / sizeof (const char *)); ++i) {
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

        const double dirpyreq_scale = min (scale_factor * 1.5, 1.0);

        for (int i = 0; i < 6; ++i) {
            adjust_radius (defaultparams.dirpyrequalizer.mult[i], dirpyreq_scale,
                           params.dirpyrequalizer.mult[i]);
        }

        params.dirpyrequalizer.threshold *= scale_factor;

        adjust_radius (defaultparams.defringe.radius, scale_factor,
                       params.defringe.radius);
        adjust_radius (defaultparams.sh.radius, scale_factor, params.sh.radius);

        if (params.raw.xtranssensor.method ==
                procparams::RAWParams::XTransSensor::methodstring[
             procparams::RAWParams::XTransSensor::threePass]) {
            params.raw.xtranssensor.method =
                procparams::RAWParams::XTransSensor::methodstring[
            procparams::RAWParams::XTransSensor::onePass];
        }

        if (params.raw.bayersensor.method == procparams::RAWParams::BayerSensor::methodstring[procparams::RAWParams::BayerSensor::pixelshift]) {
            params.raw.bayersensor.method = procparams::RAWParams::BayerSensor::methodstring[params.raw.bayersensor.pixelShiftLmmse ? procparams::RAWParams::BayerSensor::lmmse : procparams::RAWParams::BayerSensor::amaze];
        }
    }

private:
    ProcessingJobImpl* job;
    int& errorCode;
    ProgressListener* pl;
    bool tunnelMetaData;
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


IImage16* processImage (ProcessingJob* pjob, int& errorCode, ProgressListener* pl, bool tunnelMetaData, bool flush)
{
    ImageProcessor proc (pjob, errorCode, pl, tunnelMetaData, flush);
    return proc();
}

void batchProcessingThread (ProcessingJob* job, BatchProcessingListener* bpl, bool tunnelMetaData)
{

    ProcessingJob* currentJob = job;

    while (currentJob) {
        int errorCode;
        IImage16* img = processImage (currentJob, errorCode, bpl, tunnelMetaData, true);

        if (errorCode) {
            bpl->error (M ("MAIN_MSG_CANNOTLOAD"));
            currentJob = nullptr;
        } else {
            try {
                currentJob = bpl->imageReady (img);
            } catch (Glib::Exception& ex) {
                bpl->error (ex.what());
                currentJob = nullptr;
            }
        }
    }
}

void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl, bool tunnelMetaData)
{

    if (bpl) {
        Glib::Thread::create (sigc::bind (sigc::ptr_fun (batchProcessingThread), job, bpl, tunnelMetaData), 0, true, true, Glib::THREAD_PRIORITY_LOW);
    }

}

}
