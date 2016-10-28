////////////////////////////////////////////////////////////////
//
//
//
//
//  code dated: December , 2014
//
//  Ipwaveletcc is free software: you can redistribute it and/or modify
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
// *  2014 Jacques Desmis <jdesmis@gmail.com>
// *  2014 Ingo Weyrich <heckflosse@i-weyrich.de>

//
////////////////////////////////////////////////////////////////

#include <cassert>
#include <cmath>

#include "../rtgui/threadutils.h"

#include "rtengine.h"
#include "improcfun.h"
#include "LUT.h"
#include "array2D.h"
#include "boxblur.h"
#include "rt_math.h"
#include "mytime.h"
#include "sleef.c"
#include "opthelper.h"
#include "median.h"
#include "EdgePreservingDecomposition.h"
#include "iccstore.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cplx_wavelet_dec.h"

#define TS 64       // Tile size
#define offset 25   // shift between tiles
#define fTS ((TS/2+1))  // second dimension of Fourier tiles
#define blkrad 1    // radius of block averaging

#define epsilon 0.001f/(TS*TS) //tolerance


namespace rtengine
{

extern const Settings* settings;

struct cont_params {
    float mul[10];
    int chrom;
    int chro;
    int contrast;
    float th;
    float thH;
    float conres;
    float conresH;
    float chrores;
    float hueres;
    float sky;
    float b_l, t_l, b_r, t_r;
    float b_ly, t_ly, b_ry, t_ry;
    float b_lsl, t_lsl, b_rsl, t_rsl;
    float b_lhl, t_lhl, b_rhl, t_rhl;
    float edg_low, edg_mean, edg_sd, edg_max;
    float lev0s, lev0n, lev1s, lev1n, lev2s, lev2n, lev3s, lev3n;
    float b_lpast, t_lpast, b_rpast, t_rpast;
    float b_lsat, t_lsat, b_rsat, t_rsat;
    int rad;
    int val;
    int til;
    int numlevH, numlevS;
    float mulC[9];
    float mulopaRG[9];
    float mulopaBY[9];
    bool curv;
    bool opaBY;
    bool opaRG;
    bool edgcurv;
    bool diagcurv;
    int CHmet;
    int CHSLmet;
    int EDmet;
    bool HSmet;
    bool avoi;
    float strength;
    int reinforce;
    bool detectedge;
    int backm;
    float eddet;
    float eddetthr;
    float eddetthrHi;
    bool link;
    bool lip3;
    bool tonemap;
    bool diag;
    int TMmeth;
    float tmstrength;
    float balan;
    int ite;
    int contmet;
    bool opaW;
    int BAmet;
    bool bam;
    float blhigh;
    float grhigh;
    float blmed;
    float grmed;
    float bllow;
    float grlow;
    bool cbena;
    bool contena;
    bool chromena;
    bool edgeena;
    bool resena;
    bool finena;
    bool toningena;
    bool noiseena;
    int maxilev;
    float edgsens;
    float edgampl;
    int neigh;
    bool lipp;
};

int wavNestedLevels = 1;


SSEFUNCTION void ImProcFunctions::ip_wavelet(LabImage * lab, LabImage * dst, int kall, const procparams::WaveletParams & waparams, const WavCurve & wavCLVCcurve, const WavOpacityCurveRG & waOpacityCurveRG, const WavOpacityCurveBY & waOpacityCurveBY,  const WavOpacityCurveW & waOpacityCurveW, const WavOpacityCurveWL & waOpacityCurveWL, LUTf &wavclCurve, bool wavcontlutili, int skip)


{
#ifdef _DEBUG
    // init variables to display Munsell corrections
    MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif
    TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
    double wip[3][3] = {
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };
    const short int imheight = lab->H, imwidth = lab->W;
    struct cont_params cp;
    cp.avoi = params->wavelet.avoid;

    if(params->wavelet.Medgreinf == "more") {
        cp.reinforce = 1;
    }

    if(params->wavelet.Medgreinf == "none") {
        cp.reinforce = 2;
    }

    if(params->wavelet.Medgreinf == "less") {
        cp.reinforce = 3;
    }

    if(params->wavelet.NPmethod == "none") {
        cp.lip3 = false;
    }

    if(params->wavelet.NPmethod == "low") {
        cp.lip3 = true;
        cp.neigh = 0;
    }

    if(params->wavelet.NPmethod == "high") {
        cp.lip3 = true;
        cp.neigh = 1;
    }

    cp.lipp = params->wavelet.lipst;
    cp.diag = params->wavelet.tmr;
    cp.balan = (float)params->wavelet.balance;
    cp.ite = params->wavelet.iter;
    cp.tonemap = false;
    cp.bam = false;

    if(params->wavelet.tmrs == 0) {
        cp.tonemap = false;
    } else {
        cp.tonemap = true;
    }

    if(params->wavelet.TMmethod == "cont") {
        cp.contmet = 1;
    } else if(params->wavelet.TMmethod == "tm") {
        cp.contmet = 2;
    }

    if(params->wavelet.BAmethod != "none") {
        cp.bam = true;
    }

    if(params->wavelet.BAmethod == "sli") {
        cp.BAmet = 1;
    }

    if(params->wavelet.BAmethod == "cur") {
        cp.BAmet = 2;
    }

    cp.tmstrength = params->wavelet.tmrs;
    //cp.tonemap = params->wavelet.tmr;
    cp.contena = params->wavelet.expcontrast;
    cp.chromena = params->wavelet.expchroma;
    cp.edgeena = params->wavelet.expedge;
    cp.resena = params->wavelet.expresid;
    cp.finena = params->wavelet.expfinal;
    cp.toningena = params->wavelet.exptoning;
    cp.noiseena = params->wavelet.expnoise;

    if(params->wavelet.Backmethod == "black") {
        cp.backm = 0;
    }

    if(params->wavelet.Backmethod == "grey") {
        cp.backm = 1;
    }

    if(params->wavelet.Backmethod == "resid") {
        cp.backm = 2;
    }

    cp.link = params->wavelet.linkedg;
    cp.eddet = (float) params->wavelet.edgedetect;
    cp.eddetthr = (float) params->wavelet.edgedetectthr;
    cp.eddetthrHi = (float) params->wavelet.edgedetectthr2;

    cp.edgsens = 60.f;
    cp.edgampl = 10.f;

    if(cp.lipp) {
        cp.edgsens = (float) params->wavelet.edgesensi;
        cp.edgampl = (float) params->wavelet.edgeampli;
    }

    int N = imheight * imwidth;
    int maxmul = params->wavelet.thres;
    cp.maxilev = maxmul;
    static const float scales[10] = {1.f, 2.f, 4.f, 8.f, 16.f, 32.f, 64.f, 128.f, 256.f, 512.f};
    float scaleskip[10];

    for(int sc = 0; sc < 10; sc++) {
        scaleskip[sc] = scales[sc] / skip;
    }

    float atten0 = 0.40f;
    float atten123 = 0.90f;

    //int DaubLen = settings->daubech ? 8 : 6;
    int DaubLen;

    if(params->wavelet.daubcoeffmethod == "2_") {
        DaubLen = 4;
    }

    if(params->wavelet.daubcoeffmethod == "4_") {
        DaubLen = 6;
    }

    if(params->wavelet.daubcoeffmethod == "6_") {
        DaubLen = 8;
    }

    if(params->wavelet.daubcoeffmethod == "10_") {
        DaubLen = 12;
    }

    if(params->wavelet.daubcoeffmethod == "14_") {
        DaubLen = 16;
    }

    cp.CHSLmet = 1;
//  if(params->wavelet.CHSLmethod=="SL")    cp.CHSLmet=1;
//  if(params->wavelet.CHSLmethod=="CU")    cp.CHSLmet=2;
    cp.EDmet = 1;

    if(params->wavelet.EDmethod == "SL") {
        cp.EDmet = 1;
    }

    if(params->wavelet.EDmethod == "CU") {
        cp.EDmet = 2;
    }

    cp.cbena = params->wavelet.cbenab;
    cp.blhigh = (float)params->wavelet.bluehigh;
    cp.grhigh = (float)params->wavelet.greenhigh;
    cp.blmed = (float)params->wavelet.bluemed;
    cp.grmed = (float)params->wavelet.greenmed;
    cp.bllow = (float)params->wavelet.bluelow;
    cp.grlow = (float)params->wavelet.greenlow;
    cp.curv = false;
    cp.edgcurv = false;
    cp.diagcurv = false;
    cp.opaRG = false;
    cp.opaBY = false;
    cp.opaW = false;
    cp.CHmet = 0;
    cp.HSmet = false;

    if(params->wavelet.CHmethod == "with") {
        cp.CHmet = 1;
    }

    if(params->wavelet.CHmethod == "link") {
        cp.CHmet = 2;
    }

    if(params->wavelet.HSmethod == "with") {
        cp.HSmet = true;
    }

    cp.strength = min(1.f, max(0.f, ((float)params->wavelet.strength / 100.f)));

    for(int m = 0; m < maxmul; m++) {
        cp.mulC[m] = waparams.ch[m];
    }

    if(waOpacityCurveRG) {
        cp.opaRG = true;
    }

    if(cp.opaRG) {
        cp.mulopaRG[0] = 200.f * (waOpacityCurveRG[0] - 0.5f);
        cp.mulopaRG[1] = 200.f * (waOpacityCurveRG[62] - 0.5f);
        cp.mulopaRG[2] = 200.f * (waOpacityCurveRG[125] - 0.5f);
        cp.mulopaRG[3] = 200.f * (waOpacityCurveRG[187] - 0.5f);
        cp.mulopaRG[4] = 200.f * (waOpacityCurveRG[250] - 0.5f);
        cp.mulopaRG[5] = 200.f * (waOpacityCurveRG[312] - 0.5f);
        cp.mulopaRG[6] = 200.f * (waOpacityCurveRG[375] - 0.5f);
        cp.mulopaRG[7] = 200.f * (waOpacityCurveRG[438] - 0.5f);
        cp.mulopaRG[8] = 200.f * (waOpacityCurveRG[500] - 0.5f);
    } else {
        for(int level = 0; level < 9; level++) {
            cp.mulopaRG[level] = 0.f;
        }
    }

    if(waOpacityCurveBY) {
        cp.opaBY = true;
    }

    if(cp.opaBY) {
        cp.mulopaBY[0] = 200.f * (waOpacityCurveBY[0] - 0.5f);
        cp.mulopaBY[1] = 200.f * (waOpacityCurveBY[62] - 0.5f);
        cp.mulopaBY[2] = 200.f * (waOpacityCurveBY[125] - 0.5f);
        cp.mulopaBY[3] = 200.f * (waOpacityCurveBY[187] - 0.5f);
        cp.mulopaBY[4] = 200.f * (waOpacityCurveBY[250] - 0.5f);
        cp.mulopaBY[5] = 200.f * (waOpacityCurveBY[312] - 0.5f);
        cp.mulopaBY[6] = 200.f * (waOpacityCurveBY[375] - 0.5f);
        cp.mulopaBY[7] = 200.f * (waOpacityCurveBY[438] - 0.5f);
        cp.mulopaBY[8] = 200.f * (waOpacityCurveBY[500] - 0.5f);
    } else {
        for(int level = 0; level < 9; level++) {
            cp.mulopaBY[level] = 0.f;
        }
    }

    if(wavCLVCcurve) {
        cp.edgcurv = true;
    }

    if(waOpacityCurveWL) {
        cp.diagcurv = true;
    }

    for(int m = 0; m < maxmul; m++) {
        cp.mul[m] = waparams.c[m];
    }

    cp.mul[9] = (float) waparams.sup;

    for(int sc = 0; sc < 10; sc++) { //reduce strength if zoom < 100%  for contrast
        if(sc == 0) {
            if(scaleskip[sc] < 1.f) {
                cp.mul[sc] *= (atten0 * scaleskip[sc]);
            }
        } else {
            if(scaleskip[sc] < 1.f) {
                cp.mul[sc] *= (atten123 * scaleskip[sc]);
            }
        }
    }

//  if(settings->verbose) printf("Wav mul 0=%f 1=%f 2=%f 3=%f 4=%f 5=%f 6=%f 7=%f 8=%f 9=%f\n",cp.mul[0],cp.mul[1],cp.mul[2],cp.mul[3],cp.mul[4],cp.mul[5],cp.mul[6],cp.mul[7],cp.mul[8],cp.mul[9]);
    for(int sc = 0; sc < 9; sc++) { //reduce strength if zoom < 100%  for chroma and tuning
        if(sc == 0) {
            if(scaleskip[sc] < 1.f) {
                cp.mulC[sc] *= (atten0 * scaleskip[sc]);
                cp.mulopaRG[sc] *= (atten0 * scaleskip[sc]);
                cp.mulopaBY[sc] *= (atten0 * scaleskip[sc]);
            }
        } else {
            if(scaleskip[sc] < 1.f) {
                cp.mulC[sc] *= (atten123 * scaleskip[sc]);
                cp.mulopaRG[sc] *= (atten123 * scaleskip[sc]);
                cp.mulopaBY[sc] *= (atten123 * scaleskip[sc]);
            }
        }
    }

    cp.chro = waparams.chro;
    cp.chrom = waparams.chroma;
    cp.contrast = waparams.contrast;
    cp.rad = waparams.edgrad;
    cp.val = waparams.edgval;
    cp.til = waparams.edgthresh;

    cp.conres = waparams.rescon;
    cp.conresH = waparams.resconH;
    cp.chrores = waparams.reschro;
    //cp.hueres=waparams.reshue;
    cp.hueres = 2.f;
    cp.th = float(waparams.thr);
    cp.thH = float(waparams.thrH);
    cp.sky = waparams.sky;
    //skin
    cp.b_l = static_cast<float>(params->wavelet.hueskin.value[0]) / 100.0f;
    cp.t_l = static_cast<float>(params->wavelet.hueskin.value[1]) / 100.0f;
    cp.b_r = static_cast<float>(params->wavelet.hueskin.value[2]) / 100.0f;
    cp.t_r = static_cast<float>(params->wavelet.hueskin.value[3]) / 100.0f;

    cp.b_ly = static_cast<float>(params->wavelet.hueskin2.value[0]) / 100.0f;
    cp.t_ly = static_cast<float>(params->wavelet.hueskin2.value[1]) / 100.0f;
    cp.b_ry = static_cast<float>(params->wavelet.hueskin2.value[2]) / 100.0f;
    cp.t_ry = static_cast<float>(params->wavelet.hueskin2.value[3]) / 100.0f;
    cp.numlevH = params->wavelet.threshold;

    //shadows
    cp.b_lsl = static_cast<float>(params->wavelet.bllev.value[0]);
    cp.t_lsl = static_cast<float>(params->wavelet.bllev.value[1]);
    cp.b_rsl = static_cast<float>(params->wavelet.bllev.value[2]);
    cp.t_rsl = static_cast<float>(params->wavelet.bllev.value[3]);
    cp.numlevS = params->wavelet.threshold2;
    int maxlevS = 9 - cp.numlevH;
    cp.numlevS = MIN(cp.numlevS, maxlevS);
    //printf("levHigh=%d levShad=%d\n",cp.numlevH,cp.numlevS);
    //highlight
    cp.b_lhl = static_cast<float>(params->wavelet.hllev.value[0]);
    cp.t_lhl = static_cast<float>(params->wavelet.hllev.value[1]);
    cp.b_rhl = static_cast<float>(params->wavelet.hllev.value[2]);
    cp.t_rhl = static_cast<float>(params->wavelet.hllev.value[3]);
    //printf("BL=%f TL=%f BR=%f TR=%f\n",cp.b_lhl,cp.t_lhl,cp.b_rhl,cp.t_rhl);
    //pastel
    cp.b_lpast = static_cast<float>(params->wavelet.pastlev.value[0]);
    cp.t_lpast = static_cast<float>(params->wavelet.pastlev.value[1]);
    cp.b_rpast = static_cast<float>(params->wavelet.pastlev.value[2]);
    cp.t_rpast = static_cast<float>(params->wavelet.pastlev.value[3]);
    //saturated
    cp.b_lsat = static_cast<float>(params->wavelet.satlev.value[0]);
    cp.t_lsat = static_cast<float>(params->wavelet.satlev.value[1]);
    cp.b_rsat = static_cast<float>(params->wavelet.satlev.value[2]);
    cp.t_rsat = static_cast<float>(params->wavelet.satlev.value[3]);
    //edge local contrast
    cp.edg_low = static_cast<float>(params->wavelet.edgcont.value[0]);
    cp.edg_mean = static_cast<float>(params->wavelet.edgcont.value[1]);
    cp.edg_max = static_cast<float>(params->wavelet.edgcont.value[2]);
    cp.edg_sd = static_cast<float>(params->wavelet.edgcont.value[3]);
    //level noise
    cp.lev0s = static_cast<float>(params->wavelet.level0noise.value[0]);
    cp.lev0n = static_cast<float>(params->wavelet.level0noise.value[1]);
    cp.lev1s = static_cast<float>(params->wavelet.level1noise.value[0]);
    cp.lev1n = static_cast<float>(params->wavelet.level1noise.value[1]);
    cp.lev2s = static_cast<float>(params->wavelet.level2noise.value[0]);
    cp.lev2n = static_cast<float>(params->wavelet.level2noise.value[1]);
    cp.lev3s = static_cast<float>(params->wavelet.level3noise.value[0]);
    cp.lev3n = static_cast<float>(params->wavelet.level3noise.value[1]);

    cp.detectedge = params->wavelet.medianlev;
    //printf("low=%f mean=%f sd=%f max=%f\n",cp.edg_low,cp.edg_mean,cp.edg_sd,cp.edg_max);
    int minwin = min(imwidth, imheight);
    int maxlevelcrop = 9;

    if(cp.mul[9] != 0) {
        maxlevelcrop = 10;
    }

    // adap maximum level wavelet to size of crop
    if(minwin * skip < 1024) {
        maxlevelcrop = 9;    //sampling wavelet 512
    }

    if(minwin * skip < 512) {
        maxlevelcrop = 8;    //sampling wavelet 256
    }

    if(minwin * skip < 256) {
        maxlevelcrop = 7;    //sampling 128
    }

    if(minwin * skip < 128) {
        maxlevelcrop = 6;
    }

    if(minwin < 64) {
        maxlevelcrop = 5;
    }

    //  printf("minwin=%d maxcrop=%d\n",minwin, maxlevelcrop);

    int levwav = params->wavelet.thres;

    if(levwav == 9 && cp.mul[9] != 0) {
        levwav = 10;
    }

    levwav = min(maxlevelcrop, levwav);

    // determine number of levels to process.
    //  for(levwav=min(maxlevelcrop,levwav);levwav>0;levwav--)
    //      if(cp.mul[levwav-1]!=0.f  || cp.curv)
    //  if(cp.mul[levwav-1]!=0.f)
    //          break;
    // I suppress this fonctionality ==> crash for level < 3
    if(levwav < 1) {
        return;    // nothing to do
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // begin tile processing of image

    //output buffer
    int realtile;

    if(params->wavelet.Tilesmethod == "big") {
        realtile = 22;
    }

    if(params->wavelet.Tilesmethod == "lit") {
        realtile = 12;
    }

    int tilesize = 128 * realtile;
    int overlap = (int) tilesize * 0.125f;
    int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;

    if(params->wavelet.Tilesmethod == "full") {
        kall = 0;
    }

    Tile_calc (tilesize, overlap, kall, imwidth, imheight, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);

    const int numtiles = numtiles_W * numtiles_H;
    LabImage * dsttmp;

    if(numtiles == 1) {
        dsttmp = dst;
    } else {
        dsttmp = new LabImage(imwidth, imheight);

        for (int n = 0; n < 3 * imwidth * imheight; n++) {
            dsttmp->data[n] = 0;
        }
    }

    //now we have tile dimensions, overlaps
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int minsizetile = min(tilewidth, tileheight);
    int maxlev2 = 10;

    if(minsizetile < 1024 && levwav == 10) {
        maxlev2 = 9;
    }

    if(minsizetile < 512) {
        maxlev2 = 8;
    }

    if(minsizetile < 256) {
        maxlev2 = 7;
    }

    if(minsizetile < 128) {
        maxlev2 = 6;
    }

    levwav = min(maxlev2, levwav);

    //printf("levwav = %d\n",levwav);

#ifdef _OPENMP
    int numthreads = 1;
    int maxnumberofthreadsforwavelet = 0;

    //reduce memory for big tile size
    if(kall != 0) {
        if(realtile <= 22) {
            maxnumberofthreadsforwavelet = 2;
        }

        if(realtile <= 20) {
            maxnumberofthreadsforwavelet = 3;
        }

        if(realtile <= 18) {
            maxnumberofthreadsforwavelet = 4;
        }

        if(realtile <= 16) {
            maxnumberofthreadsforwavelet = 6;
        }

        if(realtile <= 14) {
            maxnumberofthreadsforwavelet = 8;
        }

        //printf("maxNRT=%d\n",maxnumberofthreadsforwavelet);
        if((maxnumberofthreadsforwavelet == 6 || maxnumberofthreadsforwavelet == 8)  && levwav == 10) {
            maxnumberofthreadsforwavelet -= 2;
        }

        if(levwav <= 7 && maxnumberofthreadsforwavelet == 8) {
            maxnumberofthreadsforwavelet = 0;
        }
    }

    //printf("maxthre=%d\n",maxnumberofthreadsforwavelet);


    // Calculate number of tiles. If less than omp_get_max_threads(), then limit num_threads to number of tiles
    if( options.rgbDenoiseThreadLimit > 0) {
        maxnumberofthreadsforwavelet = min(max(options.rgbDenoiseThreadLimit / 2, 1), maxnumberofthreadsforwavelet);
    }

    numthreads = MIN(numtiles, omp_get_max_threads());

    if(maxnumberofthreadsforwavelet > 0) {
        numthreads = MIN(numthreads, maxnumberofthreadsforwavelet);
    }

#ifdef _RT_NESTED_OPENMP
    wavNestedLevels = omp_get_max_threads() / numthreads;
    bool oldNested = omp_get_nested();

    if(wavNestedLevels < 2) {
        wavNestedLevels = 1;
    } else {
        omp_set_nested(true);
    }

    if(maxnumberofthreadsforwavelet > 0)
        while(wavNestedLevels * numthreads > maxnumberofthreadsforwavelet) {
            wavNestedLevels--;
        }

#endif

    if(settings->verbose) {
        printf("Ip Wavelet uses %d main thread(s) and up to %d nested thread(s) for each main thread\n", numthreads, wavNestedLevels);
    }

    #pragma omp parallel num_threads(numthreads)
#endif
    {
        float *mean = new float [9];
        float *meanN = new float [9];
        float *sigma = new float [9];
        float *sigmaN = new float [9];
        float *MaxP = new float [9];
        float *MaxN = new float [9];

        float** varhue = new float*[tileheight];

        for (int i = 0; i < tileheight; i++) {
            varhue[i] = new float[tilewidth];
        }

        float** varchro = new float*[tileheight];

        for (int i = 0; i < tileheight; i++) {
            varchro[i] = new float[tilewidth];
        }

#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif

        for (int tiletop = 0; tiletop < imheight; tiletop += tileHskip) {
            for (int tileleft = 0; tileleft < imwidth ; tileleft += tileWskip) {
                int tileright = MIN(imwidth, tileleft + tilewidth);
                int tilebottom = MIN(imheight, tiletop + tileheight);
                int width  = tileright - tileleft;
                int height = tilebottom - tiletop;
                LabImage * labco;
                float **Lold;
                float *LoldBuffer = nullptr;

                if(numtiles == 1) { // untiled processing => we can use output buffer for labco
                    labco = dst;

                    if(cp.avoi) { // we need a buffer to hold a copy of the L channel
                        Lold = new float*[tileheight];
                        LoldBuffer = new float[tilewidth * tileheight];
                        memcpy(LoldBuffer, lab->L[0], tilewidth * tileheight * sizeof(float));

                        for (int i = 0; i < tileheight; i++) {
                            Lold[i] = LoldBuffer + i * tilewidth;
                        }
                    }

                } else {
                    labco = new LabImage(width, height);
                    Lold = lab->L;
                }

#ifdef _RT_NESTED_OPENMP
                #pragma omp parallel for num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif

                for (int i = tiletop; i < tilebottom; i++) {
                    int i1 = i - tiletop;
                    int j;
#ifdef __SSE2__
                    __m128 c327d68v = _mm_set1_ps(327.68f);
                    __m128 av, bv, huev, chrov;

                    for (j = tileleft; j < tileright - 3; j += 4) {
                        int j1 = j - tileleft;
                        av = LVFU(lab->a[i][j]);
                        bv = LVFU(lab->b[i][j]);
                        huev = xatan2f(bv, av);
                        chrov = _mm_sqrt_ps(SQRV(av) + SQRV(bv)) / c327d68v;
                        _mm_storeu_ps(&varhue[i1][j1], huev);
                        _mm_storeu_ps(&varchro[i1][j1], chrov);

                        if(labco != lab) {
                            _mm_storeu_ps(&(labco->L[i1][j1]), LVFU(lab->L[i][j]));
                            _mm_storeu_ps(&(labco->a[i1][j1]), av);
                            _mm_storeu_ps(&(labco->b[i1][j1]), bv);
                        }
                    }

#else
                    j = tileleft;
#endif

                    for (; j < tileright; j++) {
                        int j1 = j - tileleft;
                        float a = lab->a[i][j];
                        float b = lab->b[i][j];
                        varhue[i1][j1] = xatan2f(b, a);
                        varchro[i1][j1] = (sqrtf(a * a + b * b)) / 327.68f;

                        if(labco != lab) {
                            labco->L[i1][j1] = lab->L[i][j];
                            labco->a[i1][j1] = a;
                            labco->b[i1][j1] = b;
                        }
                    }
                }

                //to avoid artifacts in blue sky
                if(params->wavelet.median) {
                    float** tmL;
                    int wid = labco->W;
                    int hei = labco->H;
                    int borderL = 1;
                    tmL = new float*[hei];

                    for (int i = 0; i < hei; i++) {
                        tmL[i] = new float[wid];
                    }

                    for(int i = borderL; i < hei - borderL; i++ ) {
                        for(int j = borderL; j < wid - borderL; j++) {
                            tmL[i][j] = labco->L[i][j];
                        }
                    }

#ifdef _RT_NESTED_OPENMP
                    #pragma omp parallel for num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif

                    for (int i = 1; i < hei - 1; i++) {
                        for (int j = 1; j < wid - 1; j++) {
                            if((varhue[i][j] < -1.3f && varhue[i][j] > - 2.5f)  && (varchro[i][j] > 15.f && varchro[i][j] < 55.f) && labco->L[i][j] > 6000.f) { //blue sky + med3x3  ==> after for more effect use denoise
                                tmL[i][j] = median(labco->L[i][j] , labco->L[i - 1][j], labco->L[i + 1][j] , labco->L[i][j + 1], labco->L[i][j - 1], labco->L[i - 1][j - 1], labco->L[i - 1][j + 1], labco->L[i + 1][j - 1], labco->L[i + 1][j + 1]);    //3x3
                            }
                        }
                    }

                    for(int i = borderL; i < hei - borderL; i++ ) {
                        for(int j = borderL; j < wid - borderL; j++) {
                            labco->L[i][j] = tmL[i][j];
                        }
                    }

                    for (int i = 0; i < hei; i++) {
                        delete [] tmL[i];
                    }

                    delete [] tmL;
                    // end blue sky
                }

                if(numtiles == 1) {
                    // reduce the varhue array to get faster access in following processing and reduce peak memory usage
                    float temphue[(tilewidth + 1) / 2] ALIGNED64;

                    for (int i = 0; i < (tileheight + 1) / 2; i++) {
                        for (int j = 0; j < (tilewidth + 1) / 2; j++) {
                            temphue[j] = varhue[i * 2][j * 2];
                        }

                        delete [] varhue[i];
                        varhue[i] = new float[(tilewidth + 1) / 2];
                        memcpy(varhue[i], temphue, ((tilewidth + 1) / 2) * sizeof(float));
                    }

                    for(int i = (tileheight + 1) / 2; i < tileheight; i++) {
                        delete [] varhue[i];
                        varhue[i] = nullptr;
                    }
                } else { // reduce the varhue array to get faster access in following processing
                    for (int i = 0; i < (tileheight + 1) / 2; i++) {
                        for (int j = 0; j < (tilewidth + 1) / 2; j++) {
                            varhue[i][j] = varhue[i * 2][j * 2];
                        }
                    }
                }

                int datalen = labco->W * labco->H;

                int levwavL = levwav;
                bool ref0 = false;

                if((cp.lev0s > 0.f || cp.lev1s > 0.f || cp.lev2s > 0.f || cp.lev3s > 0.f) && cp.noiseena) {
                    ref0 = true;
                }

                //  printf("LevwavL before: %d\n",levwavL);
                if(cp.contrast == 0.f && !cp.tonemap && cp.conres == 0.f && cp.conresH == 0.f && cp.val == 0  && !ref0 && params->wavelet.CLmethod == "all") { // no processing of residual L  or edge=> we probably can reduce the number of levels
                    while(levwavL > 0 && cp.mul[levwavL - 1] == 0.f) { // cp.mul[level] == 0.f means no changes to level
                        levwavL--;
                    }
                }

                //  printf("LevwavL after: %d\n",levwavL);
                //  if(cp.noiseena){
                if(levwavL < 4 ) {
                    levwavL = 4;    //to allow edge  => I always allocate 3 (4) levels..because if user select wavelet it is to do something !!
                }

                //  }
                //  else {
                //      if(levwavL < 3) levwavL=3;//to allow edge  => I always allocate 3 (4) levels..because if user select wavelet it is to do something !!
                //  }
                if(levwavL > 0) {
                    wavelet_decomposition* Ldecomp = new wavelet_decomposition (labco->data, labco->W, labco->H, levwavL, 1, skip, max(1, wavNestedLevels), DaubLen );

                    if(!Ldecomp->memoryAllocationFailed) {

                        float madL[8][3];
#ifdef _RT_NESTED_OPENMP
                        #pragma omp parallel for schedule(dynamic) collapse(2) num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif

                        for (int lvl = 0; lvl < 4; lvl++) {
                            for (int dir = 1; dir < 4; dir++) {
                                int Wlvl_L = Ldecomp->level_W(lvl);
                                int Hlvl_L = Ldecomp->level_H(lvl);

                                float ** WavCoeffs_L = Ldecomp->level_coeffs(lvl);

                                madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                            }
                        }

                        int ind = 0;
                        bool ref = false;

                        if((cp.lev0s > 0.f || cp.lev1s > 0.f || cp.lev2s > 0.f || cp.lev3s > 0.f) && cp.noiseena) {
                            ref = true;
                        }

                        bool contr = false;

                        for(int f = 0; f < levwavL; f++) {
                            if(cp.mul[f] != 0.f) {
                                contr = true;
                            }
                        }

                        if(cp.val > 0 || ref || contr) {//edge
                            Evaluate2(*Ldecomp, cp, ind, mean, meanN, sigma, sigmaN, MaxP, MaxN, madL);
                        }

                        //init for edge and denoise
                        float vari[4];

                        vari[0] = 8.f * SQR((cp.lev0n / 125.0) * (1.0 + cp.lev0n / 25.0));
                        vari[1] = 8.f * SQR((cp.lev1n / 125.0) * (1.0 + cp.lev1n / 25.0));
                        vari[2] = 8.f * SQR((cp.lev2n / 125.0) * (1.0 + cp.lev2n / 25.0));
                        vari[3] = 8.f * SQR((cp.lev3n / 125.0) * (1.0 + cp.lev3n / 25.0));

                        if((cp.lev0n > 0.1f || cp.lev1n > 0.1f || cp.lev2n > 0.1f || cp.lev3n > 0.1f) && cp.noiseena) {
                            int edge = 1;
                            vari[0] = max(0.0001f, vari[0]);
                            vari[1] = max(0.0001f, vari[1]);
                            vari[2] = max(0.0001f, vari[2]);
                            vari[3] = max(0.0001f, vari[3]);
                            float* noisevarlum = nullptr;  // we need a dummy to pass it to WaveletDenoiseAllL

                            WaveletDenoiseAllL(*Ldecomp, noisevarlum, madL, vari, edge);
                        }

                        ind = 1;
                        //Flat curve for Contrast=f(H) in levels
                        FlatCurve* ChCurve = nullptr;//curve C=f(H)
                        bool Chutili = false;
                        ChCurve = new FlatCurve(params->wavelet.Chcurve);

                        if (!ChCurve || ChCurve->isIdentity()) {
                            if (ChCurve) {
                                delete ChCurve;
                                ChCurve = nullptr;
                            }
                        } else {
                            Chutili = true;
                        }


                        WaveletcontAllL(labco, varhue, varchro, *Ldecomp, cp, skip, mean, meanN, sigma, sigmaN, MaxP, MaxN, wavCLVCcurve, waOpacityCurveW, waOpacityCurveWL, ChCurve, Chutili);

                        if(cp.val > 0 || ref || contr  || cp.diagcurv) {//edge
                            Evaluate2(*Ldecomp, cp, ind, mean, meanN, sigma, sigmaN, MaxP, MaxN, madL);
                        }

                        WaveletcontAllLfinal(*Ldecomp, cp, mean, sigma, MaxP, waOpacityCurveWL);
                        //Evaluate2(*Ldecomp, cp, ind, mean, meanN, sigma, sigmaN, MaxP, MaxN, madL);

                        Ldecomp->reconstruct(labco->data, cp.strength);
                    }

                    delete Ldecomp;
                }

                //Flat curve for H=f(H) in residual image
                FlatCurve* hhCurve = nullptr;//curve H=f(H)
                bool hhutili = false;
                hhCurve = new FlatCurve(params->wavelet.hhcurve);

                if (!hhCurve || hhCurve->isIdentity()) {
                    if (hhCurve) {
                        delete hhCurve;
                        hhCurve = nullptr;
                    }
                } else {
                    hhutili = true;
                }


                if(!hhutili) {//always a or b
                    int levwava = levwav;

                    //  printf("Levwava before: %d\n",levwava);
                    if(cp.chrores == 0.f && params->wavelet.CLmethod == "all" && !cp.cbena) { // no processing of residual ab => we probably can reduce the number of levels
                        while(levwava > 0 && !cp.diag && (((cp.CHmet == 2 && (cp.chro == 0.f || cp.mul[levwava - 1] == 0.f )) || (cp.CHmet != 2 && (levwava == 10 || (!cp.curv  || cp.mulC[levwava - 1] == 0.f))))) && (!cp.opaRG || levwava == 10 || (cp.opaRG && cp.mulopaRG[levwava - 1] == 0.f)) && ((levwava == 10 || (cp.CHSLmet == 1 && cp.mulC[levwava - 1] == 0.f)))) {
                            levwava--;
                        }
                    }

                    //printf("Levwava after: %d\n",levwava);
                    if(levwava > 0) {
                        wavelet_decomposition* adecomp = new wavelet_decomposition (labco->data + datalen, labco->W, labco->H, levwava, 1, skip, max(1, wavNestedLevels), DaubLen );

                        if(!adecomp->memoryAllocationFailed) {
                            WaveletcontAllAB(labco, varhue, varchro, *adecomp, waOpacityCurveW, cp, true);
                            adecomp->reconstruct(labco->data + datalen, cp.strength);
                        }

                        delete adecomp;
                    }

                    int levwavb = levwav;

                    //printf("Levwavb before: %d\n",levwavb);
                    if(cp.chrores == 0.f && params->wavelet.CLmethod == "all" && !cp.cbena) { // no processing of residual ab => we probably can reduce the number of levels
                        while(levwavb > 0 &&  !cp.diag && (((cp.CHmet == 2 && (cp.chro == 0.f || cp.mul[levwavb - 1] == 0.f )) || (cp.CHmet != 2 && (levwavb == 10 || (!cp.curv || cp.mulC[levwavb - 1] == 0.f))))) && (!cp.opaBY || levwavb == 10 || (cp.opaBY && cp.mulopaBY[levwavb - 1] == 0.f)) && ((levwavb == 10 || (cp.CHSLmet == 1 && cp.mulC[levwavb - 1] == 0.f)))) {
                            levwavb--;
                        }
                    }

                    //  printf("Levwavb after: %d\n",levwavb);
                    if(levwavb > 0) {
                        wavelet_decomposition* bdecomp = new wavelet_decomposition (labco->data + 2 * datalen, labco->W, labco->H, levwavb, 1, skip, max(1, wavNestedLevels), DaubLen );

                        if(!bdecomp->memoryAllocationFailed) {
                            WaveletcontAllAB(labco, varhue, varchro, *bdecomp, waOpacityCurveW, cp, false);
                            bdecomp->reconstruct(labco->data + 2 * datalen, cp.strength);
                        }

                        delete bdecomp;
                    }
                } else {// a and b
                    int levwavab = levwav;

                    //  printf("Levwavab before: %d\n",levwavab);
                    if(cp.chrores == 0.f && !hhutili && params->wavelet.CLmethod == "all") { // no processing of residual ab => we probably can reduce the number of levels
                        while(levwavab > 0 && (((cp.CHmet == 2 && (cp.chro == 0.f || cp.mul[levwavab - 1] == 0.f )) || (cp.CHmet != 2 && (levwavab == 10 || (!cp.curv  || cp.mulC[levwavab - 1] == 0.f))))) && (!cp.opaRG || levwavab == 10 || (cp.opaRG && cp.mulopaRG[levwavab - 1] == 0.f)) && ((levwavab == 10 || (cp.CHSLmet == 1 && cp.mulC[levwavab - 1] == 0.f)))) {
                            levwavab--;
                        }
                    }

                    //  printf("Levwavab after: %d\n",levwavab);
                    if(levwavab > 0) {
                        wavelet_decomposition* adecomp = new wavelet_decomposition (labco->data + datalen, labco->W, labco->H, levwavab, 1, skip, max(1, wavNestedLevels), DaubLen );
                        wavelet_decomposition* bdecomp = new wavelet_decomposition (labco->data + 2 * datalen, labco->W, labco->H, levwavab, 1, skip, max(1, wavNestedLevels), DaubLen );

                        if(!adecomp->memoryAllocationFailed && !bdecomp->memoryAllocationFailed) {
                            WaveletcontAllAB(labco, varhue, varchro, *adecomp, waOpacityCurveW, cp, true);
                            WaveletcontAllAB(labco, varhue, varchro, *bdecomp, waOpacityCurveW, cp, false);
                            WaveletAandBAllAB(labco, varhue, varchro, *adecomp, *bdecomp, cp, waOpacityCurveW, hhCurve, hhutili );

                            adecomp->reconstruct(labco->data + datalen, cp.strength);
                            bdecomp->reconstruct(labco->data + 2 * datalen, cp.strength);

                        }

                        delete adecomp;
                        delete bdecomp;
                    }
                }

                if (hhCurve) {
                    delete hhCurve;
                }

                if(numtiles > 1 || (numtiles == 1 /*&& cp.avoi*/)) {//in all case since I add contrast curve
                    //calculate mask for feathering output tile overlaps
                    float Vmask[height + overlap] ALIGNED16;
                    float Hmask[width + overlap] ALIGNED16;

                    if(numtiles > 1) {
                        for (int i = 0; i < height; i++) {
                            Vmask[i] = 1;
                        }

                        for (int j = 0; j < width; j++) {
                            Hmask[j] = 1;
                        }

                        for (int i = 0; i < overlap; i++) {
                            float mask = SQR(sin((M_PI * i) / (2 * overlap)));

                            if (tiletop > 0) {
                                Vmask[i] = mask;
                            }

                            if (tilebottom < imheight) {
                                Vmask[height - i] = mask;
                            }

                            if (tileleft > 0) {
                                Hmask[i] = mask;
                            }

                            if (tileright < imwidth) {
                                Hmask[width - i] = mask;
                            }
                        }
                    }

                    bool highlight = params->toneCurve.hrenabled;

#ifdef _RT_NESTED_OPENMP
                    #pragma omp parallel for schedule(dynamic,16) num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif

                    for (int i = tiletop; i < tilebottom; i++) {
                        int i1 = i - tiletop;
                        float L, a, b;
#ifdef __SSE2__
                        int rowWidth = tileright - tileleft;
                        float atan2Buffer[rowWidth] ALIGNED64;
                        float chprovBuffer[rowWidth] ALIGNED64;
                        float xBuffer[rowWidth] ALIGNED64;
                        float yBuffer[rowWidth] ALIGNED64;

                        if(cp.avoi) {
                            int col;
                            __m128 av, bv;
                            __m128 cv, yv, xv;
                            __m128 zerov = _mm_setzero_ps();
                            __m128 onev = _mm_set1_ps(1.f);
                            __m128 c327d68v = _mm_set1_ps(327.68f);
                            vmask xyMask;

                            for(col = 0; col < rowWidth - 3; col += 4) {
                                av = LVFU(labco->a[i1][col]);
                                bv = LVFU(labco->b[i1][col]);
                                STVF(atan2Buffer[col], xatan2f(bv, av));

                                cv = _mm_sqrt_ps(SQRV(av) + SQRV(bv));
                                yv = av / cv;
                                xv = bv / cv;
                                xyMask = vmaskf_eq(zerov, cv);
                                yv = vself(xyMask, onev, yv);
                                xv = vself(xyMask, zerov, xv);
                                STVF(yBuffer[col], yv);
                                STVF(xBuffer[col], xv);
                                STVF(chprovBuffer[col], cv / c327d68v);

                            }

                            for(; col < rowWidth; col++) {
                                float a = labco->a[i1][col];
                                float b = labco->b[i1][col];
                                atan2Buffer[col] = xatan2f(b, a);
                                float Chprov1 = sqrtf(SQR(a) + SQR(b));
                                yBuffer[col] = (Chprov1 == 0.f) ? 1.f : a / Chprov1;
                                xBuffer[col] = (Chprov1 == 0.f) ? 0.f : b / Chprov1;
                                chprovBuffer[col] = Chprov1 / 327.68;
                            }
                        }

#endif

                        for (int j = tileleft; j < tileright; j++) {
                            int j1 = j - tileleft;

                            if(cp.avoi) { //Gamut and Munsell
#ifdef __SSE2__
                                float HH = atan2Buffer[j1];
                                float Chprov1 = chprovBuffer[j1];
                                float2 sincosv;
                                sincosv.y = yBuffer[j1];
                                sincosv.x = xBuffer[j1];
#else
                                a = labco->a[i1][j1];
                                b = labco->b[i1][j1];
                                float HH = xatan2f(b, a);
                                float Chprov1 = sqrtf(SQR(a) + SQR(b));
                                float2 sincosv;
                                sincosv.y = (Chprov1 == 0.0f) ? 1.f : a / (Chprov1);
                                sincosv.x = (Chprov1 == 0.0f) ? 0.f : b / (Chprov1);
                                Chprov1 /= 327.68f;
#endif
                                L = labco->L[i1][j1];
                                const float Lin = labco->L[i1][j1];

                                if(wavclCurve  && cp.finena) {
                                    labco->L[i1][j1] = (0.5f * Lin  + 1.5f * wavclCurve[Lin]) / 2.f;   //apply contrast curve
                                }

                                L = labco->L[i1][j1];

                                float Lprov1 = L / 327.68f;
                                float Lprov2 = Lold[i][j] / 327.68f;
                                float memChprov = varchro[i1][j1];
                                float R, G, B;
#ifdef _DEBUG
                                bool neg = false;
                                bool more_rgb = false;
                                Color::gamutLchonly(HH, sincosv, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                                Color::gamutLchonly(HH, sincosv, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif
                                L = Lprov1 * 327.68f;

                                a = 327.68f * Chprov1 * sincosv.y; //gamut
                                b = 327.68f * Chprov1 * sincosv.x; //gamut
                                float correctionHue = 0.0f; // Munsell's correction
                                float correctlum = 0.0f;
                                Lprov1 = L / 327.68f;
                                float Chprov = sqrtf(SQR(a) + SQR(b)) / 327.68f;
#ifdef _DEBUG
                                Color::AllMunsellLch(true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum, MunsDebugInfo);
#else
                                Color::AllMunsellLch(true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum);
#endif

                                if(correctionHue != 0.f || correctlum != 0.f) { // only calculate sin and cos if HH changed
                                    if(fabs(correctionHue) < 0.015f) {
                                        HH += correctlum;    // correct only if correct Munsell chroma very little.
                                    }

                                    sincosv = xsincosf(HH + correctionHue);
                                }

                                a = 327.68f * Chprov * sincosv.y; // apply Munsell
                                b = 327.68f * Chprov * sincosv.x; //aply Munsell
                            } else {//general case
                                L = labco->L[i1][j1];
                                const float Lin = labco->L[i1][j1];

                                if(wavclCurve  && cp.finena) {
                                    labco->L[i1][j1] = (0.5f * Lin + 1.5f * wavclCurve[Lin]) / 2.f;   //apply contrast curve
                                }

                                L = labco->L[i1][j1];
                                a = labco->a[i1][j1];
                                b = labco->b[i1][j1];
                            }

                            if(numtiles > 1) {
                                float factor = Vmask[i1] * Hmask[j1];
                                dsttmp->L[i][j] += factor * L;
                                dsttmp->a[i][j] += factor * a;
                                dsttmp->b[i][j] += factor * b;
                            } else {
                                dsttmp->L[i][j] = L;
                                dsttmp->a[i][j] = a;
                                dsttmp->b[i][j] = b;

                            }
                        }
                    }
                }

                if(LoldBuffer != nullptr) {
                    delete [] LoldBuffer;
                    delete [] Lold;
                }

                if(numtiles > 1) {
                    delete labco;
                }
            }
        }

        for (int i = 0; i < tileheight; i++)
            if(varhue[i] != nullptr) {
                delete [] varhue[i];
            }

        delete [] varhue;

        for (int i = 0; i < tileheight; i++) {
            delete [] varchro[i];
        }

        delete [] varchro;

        delete [] mean;
        delete [] meanN;
        delete [] sigma;
        delete [] sigmaN;

    }
#ifdef _RT_NESTED_OPENMP
    omp_set_nested(oldNested);
#endif

    if(numtiles > 1) {
        dst->CopyFrom(dsttmp);
        delete dsttmp;
    }

}

#undef TS
#undef fTS
#undef offset
#undef epsilon

void ImProcFunctions::Aver( float *  RESTRICT DataList, int datalen, float &averagePlus, float &averageNeg, float &max, float &min)
{

    //find absolute mean
    int countP = 0, countN = 0;
    float averaP = 0.f, averaN = 0.f;

    float thres = 5.f;//different fom zero to take into account only data large enough
    max = 0.f;
    min = 0.f;
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
    {
        float lmax = 0.f, lmin = 0.f;
#ifdef _RT_NESTED_OPENMP
        #pragma omp for reduction(+:averaP,averaN,countP,countN) nowait
#endif

        for(int i = 0; i < datalen; i++) {
            if(DataList[i] >= thres) {
                averaP += DataList[i];

                if(DataList[i] > lmax) {
                    lmax = DataList[i];
                }

                countP++;
            } else if(DataList[i] < -thres) {
                averaN += DataList[i];

                if(DataList[i] < lmin) {
                    lmin = DataList[i];
                }

                countN++;
            }
        }

#ifdef _RT_NESTED_OPENMP
        #pragma omp critical
#endif
        {
            max = max > lmax ? max : lmax;
            min = min < lmin ? min : lmin;
        }
    }

    if(countP > 0) {
        averagePlus = averaP / countP;
    } else {
        averagePlus = 0;
    }

    if(countN > 0) {
        averageNeg = averaN / countN;
    } else {
        averageNeg = 0;
    }

}


void ImProcFunctions::Sigma( float *  RESTRICT DataList, int datalen, float averagePlus, float averageNeg, float &sigmaPlus, float &sigmaNeg)
{
    int countP = 0, countN = 0;
    float variP = 0.f, variN = 0.f;
    float thres = 5.f;//different fom zero to take into account only data large enough

#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for reduction(+:variP,variN,countP,countN) num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif

    for(int i = 0; i < datalen; i++) {
        if(DataList[i] >= thres) {
            variP += SQR(DataList[i] - averagePlus);
            countP++;
        } else if(DataList[i] <= -thres) {
            variN += SQR(DataList[i] - averageNeg);
            countN++;
        }
    }

    if(countP > 0) {
        sigmaPlus = sqrt(variP / countP);
    } else {
        sigmaPlus = 0;
    }

    if(countN > 0) {
        sigmaNeg = sqrt(variN / countN);
    } else {
        sigmaNeg = 0;
    }

}

void ImProcFunctions::Evaluate2(wavelet_decomposition &WaveletCoeffs_L,
                                const struct cont_params& cp, int ind, float *mean, float *meanN, float *sigma, float *sigmaN, float *MaxP, float *MaxN, float madL[8][3])
{
//StopWatch Stop1("Evaluate2");
    int maxlvl = WaveletCoeffs_L.maxlevel();

    for (int lvl = 0; lvl < maxlvl; lvl++) {

        int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
        int Hlvl_L = WaveletCoeffs_L.level_H(lvl);

        int skip_L = WaveletCoeffs_L.level_stride(lvl);

        float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);

        Eval2 (WavCoeffs_L, lvl, cp, Wlvl_L, Hlvl_L, skip_L,  ind, mean, meanN, sigma, sigmaN, MaxP, MaxN, madL[lvl]);
    }

}
void ImProcFunctions::Eval2 (float ** WavCoeffs_L,  int level, const struct cont_params& cp,
                             int W_L, int H_L, int skip_L, int ind, float *mean, float *meanN, float *sigma, float *sigmaN, float *MaxP, float *MaxN, float *madL)
{

    float avLP[4], avLN[4];
    float maxL[4], minL[4];
    float sigP[4], sigN[4];
    float AvL, AvN, SL, SN, maxLP, maxLN;

    for (int dir = 1; dir < 4; dir++) {
        Aver(WavCoeffs_L[dir], W_L * H_L,  avLP[dir], avLN[dir], maxL[dir], minL[dir]);
        Sigma(WavCoeffs_L[dir], W_L * H_L, avLP[dir], avLN[dir], sigP[dir], sigN[dir]);
    }

    AvL = 0.f;
    AvN = 0.f;
    SL = 0.f;
    SN = 0.f;
    maxLP = 0.f;
    maxLN = 0.f;

    for (int dir = 1; dir < 4; dir++) {
        AvL += avLP[dir];
        AvN += avLN[dir];
        SL += sigP[dir];
        SN += sigN[dir];
        maxLP += maxL[dir];
        maxLN += minL[dir];
    }

    AvL /= 3;
    AvN /= 3;
    SL /= 3;
    SN /= 3;
    maxLP /= 3;
    maxLN /= 3;

    mean[level] = AvL;
    meanN[level] = AvN;
    sigma[level] = SL;
    sigmaN[level] = SN;
    MaxP[level] = maxLP;
    MaxN[level] = maxLN;
}

float *ImProcFunctions::ContrastDR(float *Source, int skip, struct cont_params &cp, int W_L, int H_L, float Compression, float DetailBoost, float max0, float min0, float ave, float ah, float bh, float al, float bl, float factorx, float *Contrast)
{
    int n = W_L * H_L;

    if(Contrast == nullptr) {
        Contrast = new float[n];
    }

    memcpy(Contrast, Source, n * sizeof(float));
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for
#endif

    for (int i = 0; i < W_L * H_L; i++) { //contrast
        Contrast[i] = Source[i] ;
    }

    return Contrast;
}

SSEFUNCTION float *ImProcFunctions::CompressDR(float *Source, int skip, struct cont_params &cp, int W_L, int H_L, float Compression, float DetailBoost, float max0, float min0, float ave, float ah, float bh, float al, float bl, float factorx, float *Compressed)
{

    const float eps = 0.000001f;
    int n = W_L * H_L;

#ifdef __SSE2__
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel
#endif
    {
        __m128 epsv = _mm_set1_ps( eps );
#ifdef _RT_NESTED_OPENMP
        #pragma omp for
#endif

        for(int ii = 0; ii < n - 3; ii += 4) {
            _mm_storeu_ps( &Source[ii], xlogf(LVFU(Source[ii]) + epsv));
        }
    }

    for(int ii = n - (n % 4); ii < n; ii++) {
        Source[ii] = xlogf(Source[ii] + eps);
    }

#else
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for
#endif

    for(int ii = 0; ii < n; ii++) {
        Source[ii] = xlogf(Source[ii] + eps);
    }

#endif

    float *ucr = ContrastDR(Source, skip, cp, W_L, H_L, Compression, DetailBoost, max0, min0, ave, ah, bh, al, bl, factorx);

    if(Compressed == nullptr) {
        Compressed = ucr;
    }

    float temp;

    if(DetailBoost > 0.f && DetailBoost < 0.05f ) {
        float betemp = expf(-(2.f - DetailBoost + 0.694f)) - 1.f; //0.694 = log(2)
        temp = 1.2f * xlogf( -betemp);
        temp /= 20.f;
    } else if(DetailBoost >= 0.05f && DetailBoost < 0.25f ) {
        float betemp = expf(-(2.f - DetailBoost + 0.694f)) - 1.f; //0.694 = log(2)
        temp = 1.2f * xlogf( -betemp);
        temp /= (-75.f * DetailBoost + 23.75f);
    } else if(DetailBoost >= 0.25f) {
        float betemp = expf(-(2.f - DetailBoost + 0.694f)) - 1.f; //0.694 = log(2)
        temp = 1.2f * xlogf( -betemp);
        temp /= (-2.f * DetailBoost + 5.5f);
    }

    else {
        temp = (Compression - 1.0f) / 20.f;
    }

#ifdef __SSE2__
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel
#endif
    {
        __m128 cev, uev, sourcev;
        __m128 epsv = _mm_set1_ps( eps );
        __m128 DetailBoostv = _mm_set1_ps( DetailBoost );
        __m128 tempv = _mm_set1_ps( temp );
#ifdef _RT_NESTED_OPENMP
        #pragma omp for
#endif

        for(int i = 0; i < n - 3; i += 4) {
            cev = xexpf(LVFU(Source[i]) + LVFU(ucr[i]) * (tempv)) - epsv;
            uev = xexpf(LVFU(ucr[i])) - epsv;
            sourcev = xexpf(LVFU(Source[i])) - epsv;
            _mm_storeu_ps( &Source[i], sourcev);
            _mm_storeu_ps( &Compressed[i], cev + DetailBoostv * (sourcev - uev) );
        }
    }

    for(int i = n - (n % 4); i < n; i++) {
        float ce = xexpf(Source[i] + ucr[i] * (temp)) - eps;
        float ue = xexpf(ucr[i]) - eps;
        Source[i] = xexpf(Source[i]) - eps;
        Compressed[i] = ce + DetailBoost * (Source[i] - ue);
    }

#else
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for
#endif

    for(int i = 0; i < n; i++) {
        float ce = xexpf(Source[i] + ucr[i] * (temp)) - eps;
        float ue = xexpf(ucr[i]) - eps;
        Source[i] = xexpf(Source[i]) - eps;
        Compressed[i] = ce + DetailBoost * (Source[i] - ue);
    }

#endif

    if(Compressed != ucr) {
        delete[] ucr;
    }

    return Compressed;


}

void ImProcFunctions::ContrastResid(float * WavCoeffs_L0,  unsigned int Iterates, int skip, struct cont_params &cp, int W_L, int H_L, float max0, float min0, float ave, float ah, float bh, float al, float bl, float factorx)
{
    float stren = cp.tmstrength;
    float gamm = params->wavelet.gamma;
    cp.TMmeth = 2; //default after testing

    if(cp.TMmeth == 1) {
        min0 = 0.0f;
        max0 = 32768.f;
    } else if (cp.TMmeth == 2) {
        min0 = 0.0f;
    }

#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for
#endif

    for(int i = 0; i < W_L * H_L; i++) {
        WavCoeffs_L0[i] = (WavCoeffs_L0[i] - min0) / max0;
        WavCoeffs_L0[i] *= gamm;
    }

    float Compression = expf(-stren);       //This modification turns numbers symmetric around 0 into exponents.
    float DetailBoost = stren;

    if(stren < 0.0f) {
        DetailBoost = 0.0f;    //Go with effect of exponent only if uncompressing.
    }


    CompressDR(WavCoeffs_L0, skip, cp, W_L, H_L, Compression, DetailBoost, max0, min0, ave, ah, bh, al, bl, factorx, WavCoeffs_L0);


#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for            // removed schedule(dynamic,10)
#endif

    for(int ii = 0; ii < W_L * H_L; ii++) {
        WavCoeffs_L0[ii] = WavCoeffs_L0[ii] * max0 * (1.f / gamm) + min0;
    }
}




void ImProcFunctions::EPDToneMapResid(float * WavCoeffs_L0,  unsigned int Iterates, int skip, struct cont_params& cp, int W_L, int H_L, float max0, float min0)
{


    float stren = cp.tmstrength;
    float edgest = params->epd.edgeStopping;
    float sca = params->epd.scale;
    float gamm = params->wavelet.gamma;
    float rew = params->epd.reweightingIterates;
    EdgePreservingDecomposition epd2(W_L, H_L);
    cp.TMmeth = 2; //default after testing

    if(cp.TMmeth == 1) {
        min0 = 0.0f;
        max0 = 32768.f;
    } else if (cp.TMmeth == 2) {
        min0 = 0.0f;
    }

    //  max0=32768.f;
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for
#endif

    for(int i = 0; i < W_L * H_L; i++) {
        WavCoeffs_L0[i] = (WavCoeffs_L0[i] - min0) / max0;
        WavCoeffs_L0[i] *= gamm;
    }

    float Compression = expf(-stren);       //This modification turns numbers symmetric around 0 into exponents.
    float DetailBoost = stren;

    if(stren < 0.0f) {
        DetailBoost = 0.0f;    //Go with effect of exponent only if uncompressing.
    }

    //Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
    if(Iterates == 0) {
        Iterates = (unsigned int)(edgest * 15.0f);
    }


    epd2.CompressDynamicRange(WavCoeffs_L0, (float)sca / skip, edgest, Compression, DetailBoost, Iterates, rew, WavCoeffs_L0);

    //Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel for            // removed schedule(dynamic,10)
#endif

    for(int ii = 0; ii < W_L * H_L; ii++) {
        WavCoeffs_L0[ii] = WavCoeffs_L0[ii] * max0 * (1.f / gamm) + min0;
    }
}

void ImProcFunctions::WaveletcontAllLfinal(wavelet_decomposition &WaveletCoeffs_L, struct cont_params &cp, float *mean, float *sigma, float *MaxP, const WavOpacityCurveWL & waOpacityCurveWL)
{
    int maxlvl = WaveletCoeffs_L.maxlevel();
    float * WavCoeffs_L0 = WaveletCoeffs_L.coeff0;

    for (int dir = 1; dir < 4; dir++) {
        for (int lvl = 0; lvl < maxlvl; lvl++) {
            int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
            int Hlvl_L = WaveletCoeffs_L.level_H(lvl);
            float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);
            finalContAllL (WavCoeffs_L, WavCoeffs_L0, lvl, dir, cp, Wlvl_L, Hlvl_L, mean, sigma, MaxP, waOpacityCurveWL);
        }
    }
}


void ImProcFunctions::WaveletcontAllL(LabImage * labco, float ** varhue, float **varchrom, wavelet_decomposition &WaveletCoeffs_L,
                                      struct cont_params &cp, int skip, float *mean, float *meanN, float *sigma, float *sigmaN, float *MaxP, float *MaxN, const WavCurve & wavCLVCcurve, const WavOpacityCurveW & waOpacityCurveW, const WavOpacityCurveWL & waOpacityCurveWL, FlatCurve* ChCurve, bool Chutili)
{
    int maxlvl = WaveletCoeffs_L.maxlevel();
    int W_L = WaveletCoeffs_L.level_W(0);
    int H_L = WaveletCoeffs_L.level_H(0);
    float * WavCoeffs_L0 = WaveletCoeffs_L.coeff0;

    float maxh = 2.5f; //amplification contrast above mean
    float maxl = 2.5f; //reduction contrast under mean
    float contrast = cp.contrast;
    float multL = (float)contrast * (maxl - 1.f) / 100.f + 1.f;
    float multH = (float) contrast * (maxh - 1.f) / 100.f + 1.f;
    double avedbl = 0.f; // use double precision for big summations
    float max0 = 0.f;
    float min0 = FLT_MAX;

    if(contrast != 0.f || (cp.tonemap  && cp.resena)) { // contrast = 0.f means that all will be multiplied by 1.f, so we can skip this step
#ifdef _RT_NESTED_OPENMP
        #pragma omp parallel for reduction(+:avedbl) num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif

        for (int i = 0; i < W_L * H_L; i++) {
            avedbl += WavCoeffs_L0[i];
        }

#ifdef _RT_NESTED_OPENMP
        #pragma omp parallel num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
        {
            float lminL = FLT_MAX;
            float lmaxL = 0.f;

#ifdef _RT_NESTED_OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < W_L * H_L; i++) {
                if(WavCoeffs_L0[i] < lminL) {
                    lminL = WavCoeffs_L0[i];
                }

                if(WavCoeffs_L0[i] > lmaxL) {
                    lmaxL = WavCoeffs_L0[i];
                }

            }

#ifdef _RT_NESTED_OPENMP
            #pragma omp critical
#endif
            {
                if(lminL < min0) {
                    min0 = lminL;
                }

                if(lmaxL > max0) {
                    max0 = lmaxL;
                }
            }

        }

    }

    //      printf("MAXmax0=%f MINmin0=%f\n",max0,min0);

//tone mapping
    if(cp.tonemap && cp.contmet == 2  && cp.resena) {
        //iterate = 5
        EPDToneMapResid(WavCoeffs_L0, 5, skip, cp, W_L, H_L, max0, min0);

    }

//end tonemapping


    max0 /= 327.68f;
    min0 /= 327.68f;
    float ave = avedbl / (double)(W_L * H_L);
    float av = ave / 327.68f;
    float ah = (multH - 1.f) / (av - max0); //
    float bh = 1.f - max0 * ah;
    float al = (multL - 1.f) / (av - min0);
    float bl = 1.f - min0 * al;
    float factorx = 1.f;
//      float *koeLi[9];
//      float maxkoeLi[9];
    float *koeLi[12];
    float maxkoeLi[12];

    float *koeLibuffer = nullptr;

    for(int y = 0; y < 12; y++) {
        maxkoeLi[y] = 0.f;    //9
    }

    koeLibuffer = new float[12 * H_L * W_L]; //12

    for (int i = 0; i < 12; i++) { //9
        koeLi[i] = &koeLibuffer[i * W_L * H_L];
    }

    for(int j = 0; j < 12; j++) //9
        for (int i = 0; i < W_L * H_L; i++) {
            koeLi[j][i] = 0.f;
        }

#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
    {
        if(contrast != 0.f  && cp.resena) { // contrast = 0.f means that all will be multiplied by 1.f, so we can skip this step
            {
#ifdef _RT_NESTED_OPENMP
                #pragma omp for
#endif

                for (int i = 0; i < W_L * H_L; i++) { //contrast
                    if(WavCoeffs_L0[i] < 32768.f) {
                        float prov;

                        if( WavCoeffs_L0[i] > ave) {
                            float kh = ah * (WavCoeffs_L0[i] / 327.68f) + bh;
                            prov = WavCoeffs_L0[i];
                            WavCoeffs_L0[i] = ave + kh * (WavCoeffs_L0[i] - ave);
                        } else {
                            float kl = al * (WavCoeffs_L0[i] / 327.68f) + bl;
                            prov = WavCoeffs_L0[i];
                            WavCoeffs_L0[i] = ave - kl * (ave - WavCoeffs_L0[i]);
                        }

                        float diflc = WavCoeffs_L0[i] - prov;
                        diflc *= factorx;
                        WavCoeffs_L0[i] =  prov + diflc;
                    }
                }
            }
        }

        if(cp.tonemap && cp.contmet == 1  && cp.resena) {
            float maxp = max0 * 256.f;
            float minp = min0 * 256.f;
#ifdef _RT_NESTED_OPENMP
            #pragma omp single
#endif
            ContrastResid(WavCoeffs_L0, 5, skip, cp, W_L, H_L, maxp, minp, ave, ah, bh, al, bl, factorx );
        }

#ifdef _RT_NESTED_OPENMP
        #pragma omp barrier
#endif

        if((cp.conres != 0.f || cp.conresH != 0.f) && cp.resena) { // cp.conres = 0.f and cp.comresH = 0.f means that all will be multiplied by 1.f, so we can skip this step
#ifdef _RT_NESTED_OPENMP
            #pragma omp for nowait
#endif

            for (int i = 0; i < W_L * H_L; i++) {
                float LL = WavCoeffs_L0[i];
                float LL100 = LL / 327.68f;
                float tran = 5.f;//transition
                //shadow
                float alp = 3.f; //increase contrast sahdow in lowlights  between 1 and ??

                if(cp.th > (100.f - tran)) {
                    tran = 100.f - cp.th;
                }

                if(LL100 < cp.th) {
                    float aalp = (1.f - alp) / cp.th; //no changes for LL100 = cp.th
                    float kk = aalp * LL100 + alp;
                    WavCoeffs_L0[i] *= (1.f + kk * cp.conres / 200.f);
                } else if(LL100 < cp.th + tran) {
                    float ath = -cp.conres / tran;
                    float bth = cp.conres - ath * cp.th;
                    WavCoeffs_L0[i] *= (1.f + (LL100 * ath + bth) / 200.f);
                }

                //highlight
                tran = 5.f;

                if(cp.thH < (tran)) {
                    tran = cp.thH;
                }

                if(LL100 > cp.thH) {
                    WavCoeffs_L0[i] *= (1.f + cp.conresH / 200.f);
                } else if(LL100 > (cp.thH - tran)) {
                    float athH = cp.conresH / tran;
                    float bthH = cp.conresH - athH * cp.thH;
                    WavCoeffs_L0[i] *= (1.f + (LL100 * athH + bthH) / 200.f);
                }
            }
        }

        //enabled Lipschitz..replace simple by complex edge detection
        // I found this concept on the web (doctoral thesis on medical Imaging)
        // I was inspired by the principle of Canny and Lipschitz (continuity and derivability)
        // I adapted the principle but have profoundly changed the algorithm
        // One can 1) change all parameters and found good parameters;
        //one can also chnage in calckoe
        float edd = settings->ed_detec;
        float eddlow = settings->ed_low; //5 to 40
        //  float eddlipinfl=settings->ed_lipinfl;
        //  float eddlipampl=settings->ed_lipampl;
        float eddlipinfl = 0.005f * cp.edgsens + 0.4f;
        float eddlipampl = 1.f + cp.edgampl / 50.f;
        //  float eddlow=5.f + cp.edgampl/2.f;//settings->ed_low;//5 to 40


        if(cp.detectedge) { //enabled Lipschitz control...more memory..more time...
            float *tmCBuffer = new float[H_L * W_L];
            float *tmC[H_L];

            for (int i = 0; i < H_L; i++) {
                tmC[i] = &tmCBuffer[i * W_L];
            }

#ifdef _RT_NESTED_OPENMP
            #pragma omp for schedule(dynamic) collapse(2)
#endif

            for (int lvl = 0; lvl < 4; lvl++) {
                for (int dir = 1; dir < 4; dir++) {
                    int W_L = WaveletCoeffs_L.level_W(lvl);
                    int H_L = WaveletCoeffs_L.level_H(lvl);

                    float ** WavCoeffs_LL = WaveletCoeffs_L.level_coeffs(lvl);
                    calckoe(WavCoeffs_LL, cp, koeLi, lvl , dir, W_L, H_L, edd, maxkoeLi, tmC);
                    // return convolution KoeLi and maxkoeLi of level 0 1 2 3 and Dir Horiz, Vert, Diag
                }
            }

            delete [] tmCBuffer;

            float aamp = 1.f + cp.eddetthrHi / 100.f;

            for (int lvl = 0; lvl < 4; lvl++) {
#ifdef _RT_NESTED_OPENMP
                #pragma omp for schedule(dynamic,16)
#endif

                for (int i = 1; i < H_L - 1; i++) {
                    for (int j = 1; j < W_L - 1; j++) {
                        //treatment of koeLi and maxkoeLi
                        float interm = 0.f;

                        if(cp.lip3 && cp.lipp) {
                            // comparaison between pixel and neighbours
                            const auto neigh = cp.neigh == 1;
                            const auto kneigh = neigh ? 28.f : 38.f;
                            const auto somm = neigh ? 40.f : 50.f;

                            for (int dir = 1; dir < 4; dir++) { //neighbours proxi
                                koeLi[lvl * 3 + dir - 1][i * W_L + j] = (kneigh * koeLi[lvl * 3 + dir - 1][i * W_L + j] + 2.f * koeLi[lvl * 3 + dir - 1][(i - 1) * W_L + j] + 2.f * koeLi[lvl * 3 + dir - 1][(i + 1) * W_L + j]
                                                                        + 2.f * koeLi[lvl * 3 + dir - 1][i * W_L + j + 1] + 2.f * koeLi[lvl * 3 + dir - 1][i * W_L + j - 1] + koeLi[lvl * 3 + dir - 1][(i - 1) * W_L + j - 1]
                                                                        + koeLi[lvl * 3 + dir - 1][(i - 1) * W_L + j + 1] + koeLi[lvl * 3 + dir - 1][(i + 1) * W_L + j - 1] + koeLi[lvl * 3 + dir - 1][(i + 1) * W_L + j + 1]) / somm;
                            }
                        }

                        for (int dir = 1; dir < 4; dir++) {
                            //here I evaluate combinaison of vert / diag / horiz...we are with multiplicators of the signal
                            interm += SQR(koeLi[lvl * 3 + dir - 1][i * W_L + j]);
                        }

                        interm = sqrt(interm);

//                  interm /= 1.732f;//interm = pseudo variance koeLi
                        interm *= 0.57736721f;
                        float kampli = 1.f;
                        float eps = 0.0001f;
                        // I think this double ratio (alph, beta) is better than arctg

                        float alph = koeLi[lvl * 3][i * W_L + j] / (koeLi[lvl * 3 + 1][i * W_L + j] + eps); //ratio between horizontal and vertical
                        float beta = koeLi[lvl * 3 + 2][i * W_L + j] / (koeLi[lvl * 3 + 1][i * W_L + j] + eps); //ratio between diagonal and horizontal

                        float alipinfl = (eddlipampl - 1.f) / (1.f - eddlipinfl);
                        float blipinfl = eddlipampl - alipinfl;

                        //alph evaluate the direction of the gradient regularity Lipschitz
                        // if = 1 we are on an edge
                        // if 0 we are not
                        // we can change and use log..or Arctg but why ?? we can change if need ...
                        //Liamp=1 for eddlipinfl
                        //liamp > 1 for alp >eddlipinfl and alph < 1
                        //Liamp < 1 for alp < eddlipinfl and alph > 0
                        if(alph > 1.f) {
                            alph = 1.f / alph;
                        }

                        if(beta > 1.f) {
                            beta = 1.f / beta;
                        }

                        //take into account diagonal
                        //if in same value OK
                        //if not no edge or reduction
                        float bet = 1.f;

                        //if(cp.lip3) {//enhance algorithm
                        if(alph > eddlipinfl && beta < 0.85f * eddlipinfl) { //0.85 arbitrary value ==> eliminate from edge if H V D too different
                            bet = beta;
                        }

                        //}
                        float AmpLip = 1.f;

                        if(alph > eddlipinfl) {
                            AmpLip = alipinfl * alph + blipinfl;    //If beta low reduce kampli
                            kampli = SQR(bet) * AmpLip * aamp;
                        } else {
                            AmpLip = (1.f / eddlipinfl) * SQR(SQR(alph * bet));    //Strong Reduce if beta low
                            kampli = AmpLip / aamp;
                        }

                        // comparaison betwwen pixel and neighbours to do ==> I think 3 dir above is better
                        /*      if(cp.lip3){
                                koeLi[lvl*3][i*W_L + j] = (koeLi[lvl*3][i*W_L + j] + koeLi[lvl*3][(i-1)*W_L + j] + koeLi[lvl*3][(i+1)*W_L + j]
                                        + koeLi[lvl*3][i*W_L + j+1] + koeLi[lvl*3][i*W_L + j-1] + koeLi[lvl*3][(i-1)*W_L + j-1]
                                        + koeLi[lvl*3][(i-1)*W_L + j+1] +koeLi[lvl*3][(i+1)*W_L + j-1] +koeLi[lvl*3][(i+1)*W_L + j+1])/9.f;
                                }
                        */
                        // apply to each direction Wavelet level : horizontal / vertiacle / diagonal
                        //interm += SQR(koeLi[lvl*3 + dir-1][i*W_L + j]);

                        interm *= kampli;

                        if(interm < cp.eddetthr / eddlow) {
                            interm = 0.01f;    //eliminate too low values
                        }

                        //we can change this part of algo==> not equal but ponderate
                        koeLi[lvl * 3][i * W_L + j] = koeLi[lvl * 3 + 1][i * W_L + j] = koeLi[lvl * 3 + 2][i * W_L + j] = interm; //new value
                        //here KoeLi contains values where gradient is high and coef high, and eliminate low values...
                    }
                }
            }

            // end
        }

#ifdef _RT_NESTED_OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif

        for (int dir = 1; dir < 4; dir++) {
            for (int lvl = 0; lvl < maxlvl; lvl++) {

                int Wlvl_L = WaveletCoeffs_L.level_W(lvl);
                int Hlvl_L = WaveletCoeffs_L.level_H(lvl);

                float ** WavCoeffs_L = WaveletCoeffs_L.level_coeffs(lvl);

                ContAllL (koeLi, maxkoeLi, true, maxlvl, labco,  varhue, varchrom, WavCoeffs_L, WavCoeffs_L0, lvl, dir, cp, Wlvl_L, Hlvl_L, skip, mean, meanN, sigma, sigmaN, MaxP, MaxN, wavCLVCcurve, waOpacityCurveW, ChCurve, Chutili);


            }
        }
    }

    //delete edge detection
    if(koeLibuffer) {
        delete [] koeLibuffer;
    }
}

void ImProcFunctions::WaveletAandBAllAB(LabImage * labco, float ** varhue, float **varchrom, wavelet_decomposition &WaveletCoeffs_a, wavelet_decomposition &WaveletCoeffs_b,
                                        struct cont_params &cp, const WavOpacityCurveW & waOpacityCurveW, FlatCurve* hhCurve, bool hhutili)
{
    //   StopWatch Stop1("WaveletAandBAllAB");
    if (hhutili  && cp.resena) {  // H=f(H)
        int W_L = WaveletCoeffs_a.level_W(0);
        int H_L = WaveletCoeffs_a.level_H(0);

        float * WavCoeffs_a0 = WaveletCoeffs_a.coeff0;
        float * WavCoeffs_b0 = WaveletCoeffs_b.coeff0;
#ifdef _RT_NESTED_OPENMP
        #pragma omp parallel num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
        {
#ifdef __SSE2__
            float huebuffer[W_L] ALIGNED64;
            float chrbuffer[W_L] ALIGNED64;
#endif // __SSE2__
#ifdef _RT_NESTED_OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int i = 0; i < H_L; i++) {
#ifdef __SSE2__
                // precalculate hue and chr
                int k;

                for (k = 0; k < W_L - 3; k += 4) {
                    __m128 av = LVFU(WavCoeffs_a0[i * W_L + k]);
                    __m128 bv = LVFU(WavCoeffs_b0[i * W_L + k]);
                    __m128 huev = xatan2f(bv, av);
                    __m128 chrv = _mm_sqrt_ps(SQRV(av) + SQRV(bv));
                    STVF(huebuffer[k], huev);
                    STVF(chrbuffer[k], chrv);
                }

                for(; k < W_L; k++) {
                    huebuffer[k] = xatan2f(WavCoeffs_b0[i * W_L + k], WavCoeffs_a0[i * W_L + k]);
                    chrbuffer[k] = sqrtf(SQR(WavCoeffs_b0[i * W_L + k]) + SQR(WavCoeffs_a0[i * W_L + k])) / 327.68f;
                }

#endif // __SSE2__

                for (int j = 0; j < W_L; j++) {

#ifdef __SSE2__
                    float hueR = huebuffer[j];
                    float chR = chrbuffer[j];
#else
                    float hueR = xatan2f(WavCoeffs_b0[i * W_L + j], WavCoeffs_a0[i * W_L + j]);
                    float chR = sqrtf(SQR(WavCoeffs_b0[i * W_L + j]) + SQR(WavCoeffs_a0[i * W_L + j]));
#endif
                    /*      if (editID == EUID_WW_HHCurve) {//H pipette
                                            float valpar =Color::huelab_to_huehsv2(hueR);
                                            editWhatever->v(i,j) = valpar;
                                    }
                    */
                    float valparam = float((hhCurve->getVal(Color::huelab_to_huehsv2(hueR)) - 0.5f) * 1.7f) + hueR; //get H=f(H)  1.7 optimisation !
                    float2 sincosval = xsincosf(valparam);
                    WavCoeffs_a0[i * W_L + j] = chR * sincosval.y;
                    WavCoeffs_b0[i * W_L + j] = chR * sincosval.x;
                }
            }
        }
    }

}

void ImProcFunctions::WaveletcontAllAB(LabImage * labco, float ** varhue, float **varchrom, wavelet_decomposition &WaveletCoeffs_ab, const WavOpacityCurveW & waOpacityCurveW,
                                       struct cont_params &cp, const bool useChannelA)
{

    int maxlvl = WaveletCoeffs_ab.maxlevel();
    int W_L = WaveletCoeffs_ab.level_W(0);
    int H_L = WaveletCoeffs_ab.level_H(0);

    float * WavCoeffs_ab0 = WaveletCoeffs_ab.coeff0;

#ifdef _RT_NESTED_OPENMP
    #pragma omp parallel num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif
    {
        if(cp.chrores != 0.f  && cp.resena) { // cp.chrores == 0.f means all will be multiplied by 1.f, so we can skip the processing of residual

#ifdef _RT_NESTED_OPENMP
            #pragma omp for nowait
#endif

            for (int i = 0; i < W_L * H_L; i++) {
                const float skyprot = cp.sky;
                //chroma
                int ii = i / W_L;
                int jj = i - ii * W_L;
                float modhue = varhue[ii][jj];
                float scale = 1.f;

                if(skyprot > 0.f) {
                    if((modhue < cp.t_ry && modhue > cp.t_ly)) {
                        scale = (100.f - cp.sky) / 100.1f;
                    } else if((modhue >= cp.t_ry && modhue < cp.b_ry)) {
                        scale = (100.f - cp.sky) / 100.1f;
                        float ar = (scale - 1.f) / (cp.t_ry - cp.b_ry);
                        float br = scale - cp.t_ry * ar;
                        scale = ar * modhue + br;
                    } else if((modhue > cp.b_ly && modhue < cp.t_ly)) {
                        scale = (100.f - cp.sky) / 100.1f;
                        float al = (scale - 1.f) / (-cp.b_ly + cp.t_ly);
                        float bl = scale - cp.t_ly * al;
                        scale = al * modhue + bl;
                    }
                } else if(skyprot < 0.f) {
                    if((modhue > cp.t_ry || modhue < cp.t_ly)) {
                        scale = (100.f + cp.sky) / 100.1f;
                    }

                    /*  else if((modhue >= cp.t_ry && modhue < cp.b_ry)) {
                            scale=(100.f+cp.sky)/100.1f;
                            float ar=(scale-1.f)/(cp.t_ry- cp.b_ry);
                            float br=scale-cp.t_ry*ar;
                            scale=ar*modhue+br;
                        }
                        else if((modhue > cp.b_ly && modhue < cp.t_ly)) {
                            scale=(100.f+cp.sky)/100.1f;
                            float al=(scale-1.f)/(-cp.b_ly + cp.t_ly);
                            float bl=scale-cp.t_ly*al;
                            scale=al*modhue+bl;
                        }
                    */
                }

                WavCoeffs_ab0[i] *= (1.f + cp.chrores * (scale) / 100.f);

            }
        }

        if(cp.cbena  && cp.resena) {//if user select Toning and color balance

#ifdef _RT_NESTED_OPENMP
            #pragma omp for nowait
#endif

            for (int i = 0; i < W_L * H_L; i++) {
                int ii = i / W_L;
                int jj = i - ii * W_L;
                float LL = (labco->L[ii * 2][jj * 2]) / 327.68f; //I use labco but I can use also WavCoeffs_L0 (more exact but more memory)

                float sca = 1.f; //amplifer - reducter...about 1, but perhaps 0.6 or 1.3

                if(useChannelA) {//green red (little magenta)
                    //transition to avoid artifacts with 6 between 30 to 36 and  63 to 69
                    float aa = (cp.grmed - cp.grlow) / 6.f;
                    float bb = cp.grlow - 30.f * aa;
                    float aaa = (cp.grhigh - cp.grmed) / 6.f;
                    float bbb = cp.grmed - 63.f * aaa;

                    if(LL < 30.f) { //shadows
                        WavCoeffs_ab0[i] += cp.grlow * (sca) * 300.f;
                    } else if(LL >= 30.f && LL < 36.f) { //transition
                        float tr = aa * LL + bb;
                        WavCoeffs_ab0[i] += tr * (sca) * 300.f;
                    } else if(LL >= 36.f && LL < 63.f) { //midtones
                        WavCoeffs_ab0[i] += cp.grmed * (sca) * 300.f;
                    } else if(LL >= 63.f && LL < 69.f) { //transition
                        float trh = aaa * LL + bbb;
                        WavCoeffs_ab0[i] += trh * (sca) * 300.f;
                    } else if(LL >= 69.f) { //highlights
                        WavCoeffs_ab0[i] += cp.grhigh * (sca) * 300.f;
                    }
                } else { //blue yellow
                    //transition with 6 between 30 to 36 and 63 to 69
                    float aa1 = (cp.blmed - cp.bllow) / 6.f;
                    float bb1 = cp.bllow - 30.f * aa1;
                    float aaa1 = (cp.blhigh - cp.blmed) / 6.f;
                    float bbb1 = cp.blmed - 63.f * aaa1;

                    if(LL < 30.f) {
                        WavCoeffs_ab0[i] += cp.bllow * (sca) * 300.f;
                    } else if(LL >= 30.f && LL < 36.f) {
                        float tr1 = aa1 * LL + bb1;
                        WavCoeffs_ab0[i] += tr1 * (sca) * 300.f;
                    } else if(LL >= 36.f && LL < 63.f) {
                        WavCoeffs_ab0[i] += cp.blmed * (sca) * 300.f;
                    } else if(LL >= 63.f && LL < 69.f) {
                        float trh1 = aaa1 * LL + bbb1;
                        WavCoeffs_ab0[i] += trh1 * (sca) * 300.f;
                    } else if(LL >= 69.f) {
                        WavCoeffs_ab0[i] += cp.blhigh * (sca) * 300.f;
                    }
                }
            }
        }

#ifdef _RT_NESTED_OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif

        for (int dir = 1; dir < 4; dir++) {
            for (int lvl = 0; lvl < maxlvl; lvl++) {

                int Wlvl_ab = WaveletCoeffs_ab.level_W(lvl);
                int Hlvl_ab = WaveletCoeffs_ab.level_H(lvl);

                float ** WavCoeffs_ab = WaveletCoeffs_ab.level_coeffs(lvl);
                ContAllAB (labco,  maxlvl, varhue, varchrom, WavCoeffs_ab, WavCoeffs_ab0, lvl, dir, waOpacityCurveW, cp, Wlvl_ab, Hlvl_ab, useChannelA);
            }
        }


    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void ImProcFunctions::calckoe(float ** WavCoeffs_LL, const struct cont_params& cp, float *koeLi[12], int level, int dir, int W_L, int H_L, float edd, float *maxkoeLi, float **tmC)
{
    int borderL = 2;

//  printf("cpedth=%f\n",cp.eddetthr);
    if(cp.eddetthr < 30.f) {
        borderL = 1;

        // I calculate coefficients with r size matrix 3x3 r=1 ; 5x5 r=2; 7x7 r=3
        /*
        float k[2*r][2*r];
        for(int i=1;i<=(2*r+1);i++) {
                    for(int j=1;j<=(2*r+1);j++) {
                        k[i][j]=(1.f/6.283*sigma*sigma)*exp(-SQR(i-r-1)+SQR(j-r-1)/2.f*SQR(sigma));
                    }
        }
        //I could also use Gauss.h for 3x3
        // If necessary I can put a 7x7 matrix
        */
        for (int i = 1; i < H_L - 1; i++) { //sigma=0.55
            for (int j = 1; j < W_L - 1; j++) {
                tmC[i][j] = (8.94f * WavCoeffs_LL[dir][i * W_L + j] + 1.71f * (WavCoeffs_LL[dir][(i - 1) * W_L + j] + 1.71f * WavCoeffs_LL[dir][(i + 1) * W_L + j]
                             + 1.71f * WavCoeffs_LL[dir][i * W_L + j + 1] + 1.71f * WavCoeffs_LL[dir][i * W_L + j - 1]) + 0.33f * WavCoeffs_LL[dir][(i - 1) * W_L + j - 1]
                             + 0.33f * WavCoeffs_LL[dir][(i - 1) * W_L + j + 1] + 0.33f * WavCoeffs_LL[dir][(i + 1) * W_L + j - 1] + 0.33f * WavCoeffs_LL[dir][(i + 1) * W_L + j + 1]) * 0.0584795f;
                // apply to each direction Wavelet level : horizontal / vertiacle / diagonal


            }
        }
    } else if(cp.eddetthr >= 30.f && cp.eddetthr < 50.f) {
        borderL = 1;

        for (int i = 1; i < H_L - 1; i++) { //sigma=0.85
            for (int j = 1; j < W_L - 1; j++) {
                tmC[i][j] = (4.0091f * WavCoeffs_LL[dir][i * W_L + j] + 2.0068f * (WavCoeffs_LL[dir][(i - 1) * W_L + j] + 2.0068f * WavCoeffs_LL[dir][(i + 1) * W_L + j]
                             + 2.0068f * WavCoeffs_LL[dir][i * W_L + j + 1] + 2.0068f * WavCoeffs_LL[dir][i * W_L + j - 1]) + 1.0045f * WavCoeffs_LL[dir][(i - 1) * W_L + j - 1]
                             + 1.0045f * WavCoeffs_LL[dir][(i - 1) * W_L + j + 1] + 1.0045f * WavCoeffs_LL[dir][(i + 1) * W_L + j - 1] + 1.0045f * WavCoeffs_LL[dir][(i + 1) * W_L + j + 1]) * 0.062288f;
                // apply to each direction Wavelet level : horizontal / vertiacle / diagonal


            }
        }
    }


    else if(cp.eddetthr >= 50.f && cp.eddetthr < 75.f) {
        borderL = 1;

        for (int i = 1; i < H_L - 1; i++) {
            for (int j = 1; j < W_L - 1; j++) { //sigma=1.1
                tmC[i][j] = (3.025f * WavCoeffs_LL[dir][i * W_L + j] + 2.001f * (WavCoeffs_LL[dir][(i - 1) * W_L + j] + 2.001f * WavCoeffs_LL[dir][(i + 1) * W_L + j]
                             + 2.001f * WavCoeffs_LL[dir][i * W_L + j + 1] + 2.001f * WavCoeffs_LL[dir][i * W_L + j - 1]) + 1.323f * WavCoeffs_LL[dir][(i - 1) * W_L + j - 1]
                             + 1.323f * WavCoeffs_LL[dir][(i - 1) * W_L + j + 1] + 1.323f * WavCoeffs_LL[dir][(i + 1) * W_L + j - 1] + 1.323f * WavCoeffs_LL[dir][(i + 1) * W_L + j + 1]) * 0.06127f;
            }
        }
    }

    else if(cp.eddetthr >= 75.f) {
        borderL = 2;

        //if(cp.lip3 && level > 1) {
        if(level > 1) {// do not activate 5x5 if level 0 or 1

            for (int i = 2; i < H_L - 2; i++) {
                for (int j = 2; j < W_L - 2; j++) {
                    // Gaussian 1.1
                    // 0.5 2 3 2 0.5
                    // 2 7 10 7 2
                    // 3 10 15 10 3
                    // 2 7 10 7 2
                    // 0.5 2 3 2 0.5
                    // divi 113
                    //Gaussian 1.4
                    // 2 4 5 4 2
                    // 4 9 12 9 4
                    // 5 12 15 12 5
                    // 4 9 12 9 4
                    // 2 4 5 4 2
                    // divi 159
                    if(cp.eddetthr < 85.f) { //sigma=1.1
                        tmC[i][j] = (15.f * WavCoeffs_LL[dir][i * W_L + j]  + 10.f * WavCoeffs_LL[dir][(i - 1) * W_L + j] + 10.f * WavCoeffs_LL[dir][(i + 1) * W_L + j]
                                     + 10.f * WavCoeffs_LL[dir][i * W_L + j + 1] + 10.f * WavCoeffs_LL[dir][i * W_L + j - 1] + 7.f * WavCoeffs_LL[dir][(i - 1) * W_L + j - 1]
                                     + 7.f * WavCoeffs_LL[dir][(i - 1) * W_L + j + 1] + 7.f * WavCoeffs_LL[dir][(i + 1) * W_L + j - 1] + 7.f * WavCoeffs_LL[dir][(i + 1) * W_L + j + 1]
                                     + 3.f * WavCoeffs_LL[dir][(i - 2) * W_L + j] + 3.f * WavCoeffs_LL[dir][(i + 2) * W_L + j] + 3.f * WavCoeffs_LL[dir][i * W_L + j - 2] + 3.f * WavCoeffs_LL[dir][i * W_L + j + 2]
                                     + 2.f * WavCoeffs_LL[dir][(i - 2) * W_L + j - 1] + 2.f * WavCoeffs_LL[dir][(i - 2) * W_L + j + 1] + 2.f * WavCoeffs_LL[dir][(i + 2) * W_L + j + 1] + 2.f * WavCoeffs_LL[dir][(i + 2) * W_L + j - 1]
                                     + 2.f * WavCoeffs_LL[dir][(i - 1) * W_L + j - 2] + 2.f * WavCoeffs_LL[dir][(i - 1) * W_L + j + 2] + 2.f * WavCoeffs_LL[dir][(i + 1) * W_L + j + 2] + 2.f * WavCoeffs_LL[dir][(i + 1) * W_L + j - 2]
                                     + 0.5f * WavCoeffs_LL[dir][(i - 2) * W_L + j - 2] + 0.5f * WavCoeffs_LL[dir][(i - 2) * W_L + j + 2] + 0.5f * WavCoeffs_LL[dir][(i + 2) * W_L + j - 2] + 0.5f * WavCoeffs_LL[dir][(i + 2) * W_L + j + 2]
                                    ) * 0.0088495f;

                    }

                    else {//sigma=1.4
                        tmC[i][j] = (15.f * WavCoeffs_LL[dir][i * W_L + j] + 12.f * WavCoeffs_LL[dir][(i - 1) * W_L + j] + 12.f * WavCoeffs_LL[dir][(i + 1) * W_L + j]
                                     + 12.f * WavCoeffs_LL[dir][i * W_L + j + 1] + 12.f * WavCoeffs_LL[dir][i * W_L + j - 1] + 9.f * WavCoeffs_LL[dir][(i - 1) * W_L + j - 1]
                                     + 9.f * WavCoeffs_LL[dir][(i - 1) * W_L + j + 1] + 9.f * WavCoeffs_LL[dir][(i + 1) * W_L + j - 1] + 9.f * WavCoeffs_LL[dir][(i + 1) * W_L + j + 1]
                                     + 5.f * WavCoeffs_LL[dir][(i - 2) * W_L + j] + 5.f * WavCoeffs_LL[dir][(i + 2) * W_L + j] + 5.f * WavCoeffs_LL[dir][i * W_L + j - 2] + 5.f * WavCoeffs_LL[dir][i * W_L + j + 2]
                                     + 4.f * WavCoeffs_LL[dir][(i - 2) * W_L + j - 1] + 4.f * WavCoeffs_LL[dir][(i - 2) * W_L + j + 1] + 4.f * WavCoeffs_LL[dir][(i + 2) * W_L + j + 1] + 4.f * WavCoeffs_LL[dir][(i + 2) * W_L + j - 1]
                                     + 4.f * WavCoeffs_LL[dir][(i - 1) * W_L + j - 2] + 4.f * WavCoeffs_LL[dir][(i - 1) * W_L + j + 2] + 4.f * WavCoeffs_LL[dir][(i + 1) * W_L + j + 2] + 4.f * WavCoeffs_LL[dir][(i + 1) * W_L + j - 2]
                                     + 2.f * WavCoeffs_LL[dir][(i - 2) * W_L + j - 2] + 2.f * WavCoeffs_LL[dir][(i - 2) * W_L + j + 2] + 2.f * WavCoeffs_LL[dir][(i + 2) * W_L + j - 2] + 2.f * WavCoeffs_LL[dir][(i + 2) * W_L + j + 2]
                                    ) * 0.0062893f;
                    }


                    // apply to each direction Wavelet level : horizontal / vertiacle / diagonal
                }
            }
        }

    }


    /*
    // I suppress these 2 convolutions ==> lees good results==> probably because structure data different and also I compare to original value which have + and -
        for(int i = borderL; i < H_L-borderL; i++ ) {//[-1 0 1] x==>j
            for(int j = borderL; j < W_L-borderL; j++) {
            tmC[i][j]=- WavCoeffs_LL[dir][(i)*W_L + j-1] +  WavCoeffs_LL[dir][(i)*W_L + j+1];
            }
        }
        for(int i = borderL; i < H_L-borderL; i++ ) {//[1 0 -1] y==>i
            for(int j = borderL; j < W_L-borderL; j++) {
            tmC[i][j]= - WavCoeffs_LL[dir][(i-1)*W_L + j] + WavCoeffs_LL[dir][(i+1)*W_L + j];
            }
        }
    */

    float thr = 40.f; //avoid artifact eg. noise...to test
    float thr2 = 1.5f * edd; //edd can be modified in option ed_detect
    thr2 += cp.eddet / 30.f; //to test
    float diffFactor = (cp.eddet / 100.f);

    for(int i = 0; i < H_L; i++ ) {
        for(int j = 0; j < W_L; j++) {
            koeLi[level * 3 + dir - 1][i * W_L + j] = 1.f;
        }
    }

    for(int i = borderL; i < H_L - borderL; i++ ) {
        for(int j = borderL; j < W_L - borderL; j++) {
            // my own algo : probably a little false, but simpler as Lipschitz !
            // Thr2 = maximum of the function ==> Lipsitch says = probably edge
//                              float temp = WavCoeffs_LL[dir][i*W_L + j];
//                              if(temp>=0.f &&  temp < thr) temp = thr;
//                              if(temp < 0.f &&  temp > -thr) temp = -thr;
            float temp = max(fabsf(WavCoeffs_LL[dir][i * W_L + j]), thr );
            koeLi[level * 3 + dir - 1][i * W_L + j] = min(thr2, fabs(tmC[i][j] / temp)); // limit maxi

            //it will be more complicated to calculate both Wh and Wv, but we have also Wd==> pseudo Lipschitz
            if(koeLi[level * 3 + dir - 1][i * W_L + j] > maxkoeLi[level * 3 + dir - 1]) {
                maxkoeLi[level * 3 + dir - 1] = koeLi[level * 3 + dir - 1][i * W_L + j];
            }

            float diff = maxkoeLi[level * 3 + dir - 1] - koeLi[level * 3 + dir - 1][i * W_L + j];
            diff *= diffFactor;
            koeLi[level * 3 + dir - 1][i * W_L + j] = maxkoeLi[level * 3 + dir - 1] - diff;
        }
    }

}

void ImProcFunctions::finalContAllL (float ** WavCoeffs_L, float * WavCoeffs_L0, int level, int dir, struct cont_params &cp,
                                     int W_L, int H_L, float *mean, float *sigma, float *MaxP, const WavOpacityCurveWL & waOpacityCurveWL)
{
    if(cp.diagcurv  && cp.finena && MaxP[level] > 0.f && mean[level] != 0.f && sigma[level] != 0.f ) { //curve
        float insigma = 0.666f; //SD
        float logmax = log(MaxP[level]); //log Max
        float rapX = (mean[level] + sigma[level]) / MaxP[level]; //rapport between sD / max
        float inx = log(insigma);
        float iny = log(rapX);
        float rap = inx / iny; //koef
        float asig = 0.166f / sigma[level];
        float bsig = 0.5f - asig * mean[level];
        float amean = 0.5f / mean[level];

#ifdef _RT_NESTED_OPENMP
        #pragma omp parallel for schedule(dynamic, W_L * 16) num_threads(wavNestedLevels) if(wavNestedLevels>1)
#endif

        for (int i = 0; i < W_L * H_L; i++) {
            float absciss;

            if(fabsf(WavCoeffs_L[dir][i]) >= (mean[level] + sigma[level])) { //for max
                float valcour = xlogf(fabsf(WavCoeffs_L[dir][i]));
                float valc = valcour - logmax;
                float vald = valc * rap;
                absciss = xexpf(vald);
            } else if(fabsf(WavCoeffs_L[dir][i]) >= mean[level]) {
                absciss = asig * fabsf(WavCoeffs_L[dir][i]) + bsig;
            } else {
                absciss = amean * fabsf(WavCoeffs_L[dir][i]);
            }

            float kc = waOpacityCurveWL[absciss * 500.f] - 0.5f;
            float reduceeffect = kc <= 0.f ? 1.f : 1.5f;

            float kinterm = 1.f + reduceeffect * kc;
            kinterm = kinterm <= 0.f ? 0.01f : kinterm;

            WavCoeffs_L[dir][i] *=  kinterm;
        }
    }

    int choicelevel = atoi(params->wavelet.Lmethod.data()) - 1;
    choicelevel = choicelevel == -1 ? 4 : choicelevel;

    int choiceClevel = 0;

    if(params->wavelet.CLmethod == "one") {
        choiceClevel = 0;
    } else if(params->wavelet.CLmethod == "inf") {
        choiceClevel = 1;
    } else if(params->wavelet.CLmethod == "sup") {
        choiceClevel = 2;
    } else if(params->wavelet.CLmethod == "all") {
        choiceClevel = 3;
    }

    int choiceDir = 0;

    if(params->wavelet.Dirmethod == "one") {
        choiceDir = 1;
    } else if(params->wavelet.Dirmethod == "two") {
        choiceDir = 2;
    } else if(params->wavelet.Dirmethod == "thr") {
        choiceDir = 3;
    } else if(params->wavelet.Dirmethod == "all") {
        choiceDir = 0;
    }

    int dir1 = (choiceDir == 2) ? 1 : 2;
    int dir2 = (choiceDir == 3) ? 1 : 3;

    if(choiceClevel < 3) { // not all levels visible, paint residual
        if(level == 0) {
            if(cp.backm != 2) { // nothing to change when residual is used as background
                float backGroundColor = (cp.backm == 1) ? 12000.f : 0.f;

                for (int i = 0; i < W_L * H_L; i++) {
                    WavCoeffs_L0[i] = backGroundColor;
                }
            }
        }
    }

    if(choiceClevel == 0) { // Only one level

        if(choiceDir == 0) { // All directions
            if(level != choicelevel) { // zero all for the levels != choicelevel
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_L * H_L; i++) {
                        WavCoeffs_L[dir][i] = 0.f;
                    }
                }
            }
        } else { // zero the unwanted directions for level == choicelevel

            if(choicelevel >= cp.maxilev) {
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_L * H_L; i++) {
                        WavCoeffs_L[dir][i] = 0.f;
                    }
                }
            } else if(level != choicelevel) { // zero all for the levels != choicelevel
                for (int i = 0; i < W_L * H_L; i++) {
                    WavCoeffs_L[dir1][i] = WavCoeffs_L[dir2][i] = 0.f;
                }
            }
        }
    } else if(choiceClevel == 1) { // Only below level
        if(choiceDir == 0) { // All directions
            if(level > choicelevel) {
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_L * H_L; i++) {
                        WavCoeffs_L[dir][i] = 0.f;
                    }
                }
            }
        } else { // zero the unwanted directions for level >= choicelevel
            if(level > choicelevel) {
                for (int i = 0; i < W_L * H_L; i++) {
                    WavCoeffs_L[dir1][i] = WavCoeffs_L[dir2][i] = 0.f;
                }
            }
        }
    } else if(choiceClevel == 2) { // Only above level
        if(choiceDir == 0) { // All directions
            if(level <= choicelevel) {
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_L * H_L; i++) {
                        WavCoeffs_L[dir][i] = 0.f;
                    }
                }
            }
        } else { // zero the unwanted directions for level >= choicelevel
            if(choicelevel >= cp.maxilev) {
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_L * H_L; i++) {
                        WavCoeffs_L[dir][i] = 0.f;
                    }
                }
            }


            else if(level <= choicelevel) {
                for (int i = 0; i < W_L * H_L; i++) {
                    WavCoeffs_L[dir1][i] = WavCoeffs_L[dir2][i] = 0.f;
                }
            }
        }
    }


}

void ImProcFunctions::ContAllL (float *koeLi[12], float *maxkoeLi, bool lipschitz, int maxlvl, LabImage * labco, float ** varhue, float **varchrom, float ** WavCoeffs_L, float * WavCoeffs_L0, int level, int dir, struct cont_params &cp,
                                int W_L, int H_L, int skip, float *mean, float *meanN, float *sigma, float *sigmaN, float *MaxP, float *MaxN, const WavCurve & wavCLVCcurve, const WavOpacityCurveW & waOpacityCurveW, FlatCurve* ChCurve, bool Chutili)
{
    assert (level >= 0);
    assert (maxlvl > level);

    static const float scales[10] = {1.f, 2.f, 4.f, 8.f, 16.f, 32.f, 64.f, 128.f, 256.f, 512.f};
    float scaleskip[10];

    for(int sc = 0; sc < 10; sc++) {
        scaleskip[sc] = scales[sc] / skip;
    }

    float t_r = settings->top_right;
    float t_l = settings->top_left;
    float b_r = settings->bot_right;
    float edd = settings->ed_detec;
    float eddstrength = settings->ed_detecStr;
    float aedstr = (eddstrength - 1.f) / 90.f;
    float bedstr = 1.f - 10.f * aedstr;

    if(cp.val > 0  && cp.edgeena) {
        float * koe;
        float maxkoe = 0.f;

        if(!lipschitz) {
            koe = new float [H_L * W_L];

            for (int i = 0; i < W_L * H_L; i++) {
                koe[i] = 0.f;
            }

            maxkoe = 0.f;

            if(cp.detectedge) {
                float** tmC;
                int borderL = 1;
                tmC = new float*[H_L];

                for (int i = 0; i < H_L; i++) {
                    tmC[i] = new float[W_L];
                }

                {
                    for (int i = 1; i < H_L - 1; i++) {
                        for (int j = 1; j < W_L - 1; j++) {
                            //edge detection wavelet TMC Canny
                            // also possible to detect noise with 5x5 instead of 3x3
                            tmC[i][j] = (4.f * WavCoeffs_L[dir][i * W_L + j] + 2.f * WavCoeffs_L[dir][(i - 1) * W_L + j] + 2.f * WavCoeffs_L[dir][(i + 1) * W_L + j]
                                         + 2.f * WavCoeffs_L[dir][i * W_L + j + 1] + 2.f * WavCoeffs_L[dir][i * W_L + j - 1] + WavCoeffs_L[dir][(i - 1) * W_L + j - 1]
                                         + WavCoeffs_L[dir][(i - 1) * W_L + j + 1] + WavCoeffs_L[dir][(i + 1) * W_L + j - 1] + WavCoeffs_L[dir][(i + 1) * W_L + j + 1]) / 16.f;

                            // apply to each direction Wavelet level : horizontal / vertiacle / diagonal
                        }
                    }
                }




                for(int i = borderL; i < H_L - borderL; i++ ) {
                    for(int j = borderL; j < W_L - borderL; j++) {
                        // my own algo : probably a little false, but simpler as Lipschitz !
                        float thr = 40.f; //avoid artifact eg. noise...to test
                        float thr2 = edd; //edd can be modified in option ed_detect
                        thr2 += cp.eddet / 30.f; //to test
                        float temp = WavCoeffs_L[dir][i * W_L + j];

                        if(temp >= 0.f &&  temp < thr) {
                            temp = thr;
                        }

                        if(temp < 0.f &&  temp > -thr) {
                            temp = -thr;
                        }

                        koe[i * W_L + j] = min(thr2, fabs(tmC[i][j] / temp));

                        if(koe[i * W_L + j] > maxkoe) {
                            maxkoe = koe[i * W_L + j];
                        }

                        float diff = maxkoe - koe[i * W_L + j];
                        diff *= (cp.eddet / 100.f);
                        float interm = maxkoe - diff;

                        if(interm < cp.eddetthr / 30.f) {
                            interm = 0.01f;
                        }

                        koe[i * W_L + j] = interm;

                    }
                }

                for (int i = 0; i < H_L; i++) {
                    delete [] tmC[i];
                }

                delete [] tmC;

            }
        }

        //end detect edge
        float rad = ((float)cp.rad) / 60.f; //radius ==> not too high value to avoid artifacts
        float value = ((float)cp.val) / 8.f; //strength

        if (scaleskip[1] < 1.f) {
            float atten01234 = 0.80f;
            value *= (atten01234 * scaleskip[1]);    //for zoom < 100% reduce strength...I choose level 1...but!!
        }

        float edge = 1.f;
        float lim0 = 20.f; //arbitrary limit for low radius and level between 2 or 3 to 30 maxi
        float lev = float (level);
        float repart = (float)cp.til;
        float brepart;

        if(cp.reinforce == 1) {
            brepart = 3.f;
        }

        if(cp.reinforce == 3) {
            brepart = 0.5f;    //arbitrary value to increase / decrease repart, between 1 and 0
        }

        float arepart = -(brepart - 1.f) / (lim0 / 60.f);

        if (cp.reinforce != 2) {
            if(rad < lim0 / 60.f) {
                repart *= (arepart * rad + brepart);    //linear repartition of repart
            }
        }

        float al0 = 1.f + (repart) / 50.f;
        float al10 = 1.0f; //arbitrary value ==> less = take into account high levels
        //  float ak =-(al0-al10)/10.f;//10 = maximum levels
        float ak = -(al0 - al10) / 10.f; //10 = maximum levels
        float bk = al0;
        float koef = ak * level + bk; //modulate for levels : more levels high, more koef low ==> concentrated action on low levels, without or near for high levels
        float expkoef = -pow(fabs(rad - lev), koef); //reduce effect for high levels

        if (cp.reinforce == 3) {
            if(rad < lim0 / 60.f && level == 0) {
                expkoef *= abs(repart);    //reduce effect for low values of rad and level=0==> quasi only level 1 is effective
            }
        }

        if (cp.reinforce == 1) {
            if(rad < lim0 / 60.f && level == 1) {
                expkoef /= repart;    //increase effect for low values of rad and level=1==> quasi only level 0 is effective
            }
        }

        //take into account local contrast
        float refin = value * exp (expkoef);

        if(cp.link  && cp.noiseena) { //combi
            {
                if(level == 0) {
                    refin *= (1.f + cp.lev0s / 50.f);    // we can change this sensibility!
                }

                if(level == 1) {
                    refin *= (1.f + cp.lev1s / 50.f);
                }

                if(level == 2) {
                    refin *= (1.f + cp.lev2s / 50.f);
                }

                if(level == 3) {
                    refin *= (1.f + cp.lev3s / 50.f);
                }
            }
        }

        float edgePrecalc = 1.f + refin; //estimate edge "pseudo variance"

        if(cp.EDmet == 2) { //curve
            //  if(exa) {//curve
            float insigma = 0.666f; //SD
            float logmax = log(MaxP[level]); //log Max
            float rapX = (mean[level] + sigma[level]) / MaxP[level]; //rapport between sD / max
            float inx = log(insigma);
            float iny = log(rapX);
            float rap = inx / iny; //koef
            float asig = 0.166f / sigma[level];
            float bsig = 0.5f - asig * mean[level];
            float amean = 0.5f / mean[level];
            float absciss;
            float kinterm;
            float kmul;
            int borderL = 1;

            for(int i = borderL; i < H_L - borderL; i++ ) {
                for(int j = borderL; j < W_L - borderL; j++) {
                    int k = i * W_L + j;

                    if(cp.detectedge) {
                        if(!lipschitz) {
                            if(cp.eddet > 10.f) {
                                edge = (aedstr * cp.eddet + bedstr) * (edgePrecalc * (1.f + koe[k])) / (1.f + 0.9f * maxkoe);
                            } else {
                                edge = (edgePrecalc * (1.f + koe[k])) / (1.f + 0.9f * maxkoe);
                            }
                        }

                        if(lipschitz) {
                            if(level < 4) {
                                edge = 1.f + (edgePrecalc - 1.f) * (koeLi[level * 3][k]) / (1.f + 0.9f * maxkoeLi[level * 3 + dir - 1]);
                            } else {
                                edge = edgePrecalc;
                            }
                        }
                    } else {
                        edge = edgePrecalc;
                    }

                    if(cp.edgcurv) {
                        if(fabs(WavCoeffs_L[dir][k]) >= (mean[level] + sigma[level])) { //for max
                            float valcour = log(fabs(WavCoeffs_L[dir][k]));
                            float valc = valcour - logmax;
                            float vald = valc * rap;
                            absciss = exp(vald);

                        } else if(fabs(WavCoeffs_L[dir][k]) >= mean[level] &&  fabs(WavCoeffs_L[dir][k]) < (mean[level] + sigma[level])) {
                            absciss = asig * fabs(WavCoeffs_L[dir][k]) + bsig;
                        } else if(fabs(WavCoeffs_L[dir][k]) < mean[level]) {
                            absciss = amean * fabs(WavCoeffs_L[dir][k]);
                        }

                        // Threshold adjuster settings==> approximative for curve
                        //kmul about average cbrt(3--40 / 10)==>1.5 to 2.5
                        //kmul about SD   10--60  / 35 ==> 2
                        // kmul about low  cbrt((5.f+cp.edg_low)/5.f);==> 1.5
                        // kmul about max ==> 9
                        // we can change these values
                        // result is different not best or bad than threshold slider...but similar
                        float abssd = 4.f; //amplification reference
                        float bbssd = 2.f; //mini ampli
                        float maxamp = 2.5f; //maxi ampli at end
                        float maxampd = 10.f; //maxi ampli at end
                        float a_abssd = (maxamp - abssd) / 0.333f;
                        float b_abssd = maxamp - a_abssd;
                        float da_abssd = (maxampd - abssd) / 0.333f;
                        float db_abssd = maxampd - da_abssd;
                        float am = (abssd - bbssd) / 0.666f;
                        float kmuld = 0.f;

                        if(absciss > 0.666f && absciss < 1.f) {
                            kmul = a_abssd * absciss + b_abssd;    //about max  ==> kinterm
                            kmuld = da_abssd * absciss + db_abssd;
                        } else {
                            kmul = kmuld = absciss * am + bbssd;
                        }

                        kinterm = 1.f;
                        float kc = kmul * (wavCLVCcurve[absciss * 500.f] - 0.5f);
                        float kcd = kmuld * (wavCLVCcurve[absciss * 500.f] - 0.5f);

                        if(kc >= 0.f) {
                            float reduceeffect = 0.6f;
                            kinterm = 1.f + reduceeffect * kmul * (wavCLVCcurve[absciss * 500.f] - 0.5f);    //about 1 to 3 general and big amplification for max (under 0)
                        } else {
                            kinterm = 1.f - (SQR(kcd)) / 10.f;
                        }

                        if(kinterm < 0.f) {
                            kinterm = 0.01f;
                        }

                        edge *= kinterm;

                        if(edge < 1.f) {
                            edge = 1.f;
                        }
                    }

                    WavCoeffs_L[dir][k] *=  edge;
                }
            }
        } else if(cp.EDmet == 1) { //threshold adjuster
            float MaxPCompare = MaxP[level] * SQR(cp.edg_max / 100.f); //100 instead of b_r...case if b_r < 100
            float MaxNCompare = MaxN[level] * SQR(cp.edg_max / 100.f); //always rduce a little edge for near max values
            float edgeSdCompare = (mean[level] + 1.5f * sigma[level]) * SQR(cp.edg_sd / t_r); // 1.5 standard deviation #80%  range between mean 50% and 80%
            float edgeMeanCompare = mean[level] * SQR(cp.edg_mean / t_l);
            float edgeLowCompare = (5.f + SQR(cp.edg_low));
            float edgeMeanFactor = cbrt(cp.edg_mean / t_l);
            float interm;

            if(cp.edg_low < 10.f) {
                interm = cbrt((5.f + cp.edg_low) / 5.f);
            } else {
                interm = 1.437f;    //cbrt(3);
            }

            float edgeLowFactor = interm;
            float edgeSdFactor = cp.edg_sd / t_r;
            float edgeMaxFactor = SQR(cp.edg_max / b_r);
            float edgMaxFsup = (cp.edg_max / b_r); //reduce increase of effect for high values contrast..if slider > b_r

            //for (int i=0; i<W_L*H_L; i++) {
            int borderL = 1;

            for(int i = borderL; i < H_L - borderL; i++ ) {
                for(int j = borderL; j < W_L - borderL; j++) {
                    int k = i * W_L + j;

                    if(cp.detectedge) {
                        if(!lipschitz) {
                            if(cp.eddet > 10.f) {
                                edge = (aedstr * cp.eddet + bedstr) * (edgePrecalc * (1.f + koe[k])) / (1.f + 0.9f * maxkoe);
                            } else {
                                edge = (edgePrecalc * (1.f + koe[k])) / (1.f + 0.9f * maxkoe);
                            }
                        }

                        if(lipschitz) {
                            if(level < 4) {
                                edge = 1.f + (edgePrecalc - 1.f) * (koeLi[level * 3][k]) / (1.f + 0.9f * maxkoeLi[level * 3 + dir - 1]);
                            } else {
                                edge = edgePrecalc;
                            }
                        }
                    } else {
                        edge = edgePrecalc;
                    }

                    //algorithm that take into account local contrast
                    // I use a thresholdadjuster with
                    // Bottom left ==> minimal low value for local contrast (not 0, but 5...we can change)
                    // 0 10*10 35*35 100*100 substantially correspond to the true distribution of low value, mean, standard-deviation and max (ed 5, 50, 400, 4000
                    // Top left ==> mean reference value (for each level), we can change cbrt(cp.edg_mean/10.f)
                    // Top Right==> standard deviation (for each level) we can change (cp.edg_sd/35.f)
                    // bottom right ==> Max for positif and negatif contrast we can change cp.edg_max/100.f
                    // If we move sliders to the left, local contrast is reduced
                    // if we move sliders to the right local contrast is increased
                    // MaxP, MaxN, mean, sigma are calculated if necessary (val > 0) by evaluate2(), eval2(), aver() , sigma()
                    if(b_r < 100.f  && cp.edg_max / b_r > 1.f) { //in case of b_r < 100 and slider move to right
                        if (WavCoeffs_L[dir][k] > MaxPCompare * cp.edg_max / b_r) {
                            edge *= edgMaxFsup;

                            if(edge < 1.f) {
                                edge = 1.f;
                            }
                        } else if (WavCoeffs_L[dir][k] < MaxNCompare * cp.edg_max / b_r) {
                            edge *= edgMaxFsup;

                            if(edge < 1.f) {
                                edge = 1.f;
                            }
                        }
                    }

                    if (WavCoeffs_L[dir][k] > MaxPCompare) {
                        edge *= edgeMaxFactor;

                        if(edge < 1.f) {
                            edge = 1.f;
                        }
                    }//reduce edge if > new max
                    else if (WavCoeffs_L[dir][k] < MaxNCompare) {
                        edge *= edgeMaxFactor;

                        if(edge < 1.f) {
                            edge = 1.f;
                        }
                    }

                    if (fabs(WavCoeffs_L[dir][k]) >= edgeMeanCompare && fabs(WavCoeffs_L[dir][k]) < edgeSdCompare) {
                        //if (fabs(WavCoeffs_L[dir][i]) > edgeSdCompare) {
                        edge *= edgeSdFactor;

                        if(edge < 1.f) {
                            edge = 1.f;
                        }
                    }//mofify effect if sd change

                    if (fabs(WavCoeffs_L[dir][k]) < edgeMeanCompare) {
                        edge *= edgeMeanFactor;

                        if(edge < 1.f) {
                            edge = 1.f;
                        }
                    } // modify effect if mean change

                    if (fabs(WavCoeffs_L[dir][k]) < edgeLowCompare) {
                        edge *= edgeLowFactor;

                        if(edge < 1.f) {
                            edge = 1.f;
                        }
                    }

                    WavCoeffs_L[dir][k] *= edge;
                }
            }
        }

        if(!lipschitz) {
            delete [] koe;
        }
    }


    if(!cp.link && cp.noiseena)   { //used both with denoise 1 2 3
        float refine = 0.f;

        for (int i = 0; i < W_L * H_L; i++) {
            if(level == 0) {
                refine = cp.lev0s / 40.f;
            }

            if(level == 1) {
                refine = cp.lev1s / 40.f;
            }

            if(level == 2) {
                refine = cp.lev2s / 40.f;
            }

            if(level == 3) {
                refine = cp.lev3s / 40.f;
            }

            WavCoeffs_L[dir][i] *= (1.f + refine);
        }
    }


    float cpMul = cp.mul[level];

    if(cpMul != 0.f && cp.contena) { // cpMul == 0.f means all will be multiplied by 1.f, so we can skip this

        const float skinprot = params->wavelet.skinprotect;
        const float skinprotneg = -skinprot;
        const float factorHard = (1.f - skinprotneg / 100.f);

        //to adjust increase contrast with local contrast

        //for each pixel and each level
        float beta;
        float mea[9];
        mea[0] = mean[level] / 6.f;
        mea[1] = mean[level] / 2.f;
        mea[2] = mean[level]; // 50% data
        mea[3] = mean[level] + sigma[level] / 2.f;
        mea[4] = mean[level] + sigma[level]; //66%
        mea[5] = mean[level] + 1.2f * sigma[level];
        mea[6] = mean[level] + 1.5f * sigma[level]; //
        mea[7] = mean[level] + 2.f * sigma[level]; //95%
        mea[8] = mean[level] + 2.5f * sigma[level]; //99%

        bool useChromAndHue = (skinprot != 0.f || cp.HSmet);
        float modchro;

        for (int i = 0; i < W_L * H_L; i++) {
            float kLlev = 1.f;

            if(cpMul < 0.f) {
                beta = 1.f; // disabled for negatives values "less contrast"
            } else {
                float WavCL = fabsf(WavCoeffs_L[dir][i]);

                //reduction amplification: max action between mean / 2 and mean + sigma
                // arbitrary coefficient, we can add a slider !!
                if(WavCL < mea[0]) {
                    beta = 0.6f;    //preserve very low contrast (sky...)
                } else if(WavCL < mea[1]) {
                    beta = 0.8f;
                } else if(WavCL < mea[2]) {
                    beta = 1.f;    //standard
                } else if(WavCL < mea[3]) {
                    beta = 1.f;
                } else if(WavCL < mea[4]) {
                    beta = 0.8f;    //+sigma
                } else if(WavCL < mea[5]) {
                    beta = 0.6f;
                } else if(WavCL < mea[6]) {
                    beta = 0.4f;
                } else if(WavCL < mea[7]) {
                    beta = 0.2f;    // + 2 sigma
                } else if(WavCL < mea[8]) {
                    beta = 0.1f;
                } else {
                    beta = 0.0f;
                }
            }

            float scale = 1.f;
            float scale2 = 1.f;

            float LL100, LL100res, LL100init, kH[maxlvl];

            int ii = i / W_L;
            int jj = i - ii * W_L;
            float LL = labco->L[ii * 2][jj * 2];
            LL100 = LL100init = LL / 327.68f;
            LL100res = WavCoeffs_L0[i] / 327.68f;
            float delta = fabs(LL100init - LL100res) / (maxlvl / 2);

            for(int ml = 0; ml < maxlvl; ml++) {
                if(ml < maxlvl / 2) {
                    kH[ml] = (LL100res + ml * delta) / LL100res;    // fixed a priori max to level middle
                } else {
                    kH[ml] = (LL100init - ml * delta) / LL100res;
                }
            }


            if(useChromAndHue) {
                float modhue = varhue[ii][jj];
                modchro = varchrom[ii * 2][jj * 2];
                // hue chroma skin with initial lab datas
                scale = 1.f;

                if(skinprot > 0.f) {
                    Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprot, scale, true, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 0); //0 for skin and extand
                } else if(skinprot < 0.f) {
                    Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprotneg, scale, false, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 0);

                    if (scale == 1.f) {
                        scale = factorHard;
                    } else {
                        scale = 1.f;
                    }
                }

            }

            if(Chutili) {
                int i_i = i / W_L;
                int j_j = i - i_i * W_L;
                double lr;
                float modhue2 = varhue[i_i][j_j];
                float valparam = float((ChCurve->getVal(lr = Color::huelab_to_huehsv2(modhue2)) - 0.5f)); //get valparam=f(H)

                if(valparam > 0.f) {
                    scale2 = 1.f + 3.f * valparam;    //arbitrary value
                } else {
                    scale2 = 1.f + 1.9f * valparam;    //near 0 but not zero if curve # 0
                }
            }

            //linear transition HL
            float diagacc = 1.f;
            float alpha = (1024.f + 15.f * (float) cpMul * scale * scale2 * beta * diagacc) / 1024.f ;

            if(cp.HSmet  && cp.contena) {
                float aaal = (1.f - alpha) / ((cp.b_lhl - cp.t_lhl) * kH[level]);
                float bbal = 1.f - aaal * cp.b_lhl * kH[level];
                float aaar = (alpha - 1.f) / (cp.t_rhl - cp.b_rhl) * kH[level];
                float bbbr = 1.f - cp.b_rhl * aaar * kH[level];
                //linear transition Shadows
                float aaalS = (1.f - alpha) / (cp.b_lsl - cp.t_lsl);
                float bbalS = 1.f - aaalS * cp.b_lsl;
                float aaarS = (alpha - 1.f) / (cp.t_rsl - cp.b_rsl);
                float bbbrS = 1.f - cp.b_rsl * aaarS;

                if(level <= cp.numlevH) { //in function of levels
                    if((LL100 > cp.t_lhl * kH[level] && LL100 < cp.t_rhl * kH[level])) {
                        kLlev = alpha;
                    } else if((LL100 > cp.b_lhl * kH[level] && LL100 <= cp.t_lhl * kH[level])) {
                        kLlev = aaal * LL100 + bbal;
                    } else if((LL100 > cp.t_rhl * kH[level] && LL100 <= cp.b_rhl * kH[level])) {
                        kLlev = aaar * LL100 + bbbr;
                    } else {
                        kLlev = 1.f;
                    }
                }

                if(level >= (9 - cp.numlevS)) {
                    if((LL100 > cp.t_lsl && LL100 < cp.t_rsl)) {
                        kLlev = alpha;
                    } else if((LL100 > cp.b_lsl && LL100 <= cp.t_lsl)) {
                        kLlev = aaalS * LL100 + bbalS;
                    } else if((LL100 > cp.t_rsl && LL100 <= cp.b_rsl)) {
                        kLlev = aaarS * LL100 + bbbrS;
                    } else {
                        kLlev = 1.f;
                    }
                }

            } else {
                kLlev = alpha;
            }

            WavCoeffs_L[dir][i] *= (kLlev);
        }
    }

    if(waOpacityCurveW) {
        cp.opaW = true;
    }

    if(cp.bam && cp.finena) {
        if(cp.opaW && cp.BAmet == 2) {
            int iteration = cp.ite;
            int itplus = 7 + iteration;
            int itmoins = 7 - iteration;
            int med = maxlvl / 2;
            int it;

            if(level < med) {
                it = itmoins;
            } else if(level == med) {
                it = 7;
            } else /*if(level > med)*/ {
                it = itplus;
            }

            for(int j = 0; j < it; j++) {
                //float bal = cp.balan;//-100 +100
                float kba = 1.f;

                //  if(dir <3) kba= 1.f + bal/600.f;
                //  if(dir==3) kba = 1.f - bal/300.f;
                for (int i = 0; i < W_L * H_L; i++) {
                    int ii = i / W_L;
                    int jj = i - ii * W_L;
                    float LL100 = labco->L[ii * 2][jj * 2] / 327.68f;
                    float k1 = 0.3f * (waOpacityCurveW[6.f * LL100] - 0.5f); //k1 between 0 and 0.5    0.5==> 1/6=0.16
                    float k2 = k1 * 2.f;

                    if(dir < 3) {
                        kba = 1.f + k1;
                    }

                    if(dir == 3) {
                        kba = 1.f - k2;
                    }

                    WavCoeffs_L[dir][i] *= (kba);
                }
            }
        }

        if(cp.BAmet == 1) {
            int iteration = cp.ite;
            int itplus = 7 + iteration;
            int itmoins = 7 - iteration;
            int med = maxlvl / 2;
            int it;

            if(level < med) {
                it = itmoins;
            } else if(level == med) {
                it = 7;
            } else /*if(level > med)*/ {
                it = itplus;
            }

            for(int j = 0; j < it; j++) {
                float bal = cp.balan;//-100 +100
                float kba = 1.f;

                //  if(dir <3) kba= 1.f + bal/600.f;
                //  if(dir==3) kba = 1.f - bal/300.f;
                for (int i = 0; i < W_L * H_L; i++) {
                    int ii = i / W_L;
                    int jj = i - ii * W_L;
                    float k1 = 600.f;
                    float k2 = 300.f;
                    float LL100 = labco->L[ii * 2][jj * 2] / 327.68f;
                    float aa = 4970.f;
                    float bb = -397000.f;
                    float b0 = 100000.f;
                    float a0 = -4970.f;

                    if(LL100 > 80.f) {
                        k1 = aa * LL100 + bb;
                        k2 = 0.5f * k1;
                    }

                    if(LL100 < 20.f) {
                        k1 = a0 * LL100 + b0;
                        k2 = 0.5f * k1;
                    }

                    //k1=600.f;
                    //k2=300.f;
                    //k1=0.3f*(waOpacityCurveW[6.f*LL100]-0.5f);//k1 between 0 and 0.5    0.5==> 1/6=0.16
                    //k2=k1*2.f;
                    if(dir < 3) {
                        kba = 1.f + bal / k1;
                    }

                    if(dir == 3) {
                        kba = 1.f - bal / k2;
                    }

                    WavCoeffs_L[dir][i] *= (kba);
                }
            }
        }

    }

    // to see each level of wavelet ...level from 0 to 8
    int choicelevel = atoi(params->wavelet.Lmethod.data()) - 1;
    choicelevel = choicelevel == -1 ? 4 : choicelevel;
}

void ImProcFunctions::ContAllAB (LabImage * labco, int maxlvl, float ** varhue, float **varchrom, float ** WavCoeffs_ab, float * WavCoeffs_ab0, int level, int dir, const WavOpacityCurveW & waOpacityCurveW, struct cont_params &cp,
                                 int W_ab, int H_ab, const bool useChannelA)
{
    float cpMul = cp.mul[level];

    if(cpMul != 0.f && cp.CHmet == 2 && cp.chro != 0.f  && cp.chromena) { // cpMul == 0.f or cp.chro = 0.f means all will be multiplied by 1.f, so we can skip this
        const float skinprot = params->wavelet.skinprotect;
        const float skinprotneg = -skinprot;
        const float factorHard = (1.f - skinprotneg / 100.f);
        const float cpChrom = cp.chro;

        //to adjust increase contrast with local contrast
        bool useSkinControl = (skinprot != 0.f);
        float alphaC = (1024.f + 15.f * cpMul * cpChrom / 50.f) / 1024.f ;

        for (int i = 0; i < W_ab * H_ab; i++) {
            if(useSkinControl) {
                int ii = i / W_ab;
                int jj = i - ii * W_ab;
                float LL100 = labco->L[ii * 2][jj * 2] / 327.68f;
                float modhue = varhue[ii][jj];
                float modchro = varchrom[ii * 2][jj * 2];
                // hue chroma skin with initial lab datas
                float scale = 1.f;

                if(skinprot > 0.f) {
                    Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprot, scale, true, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 0); //0 for skin and extand
                } else if(skinprot < 0.f) {
                    Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprotneg, scale, false, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 0);
                    scale = (scale == 1.f) ? factorHard : 1.f;
                }

                alphaC = (1024.f + 15.f * cpMul * cpChrom * scale / 50.f) / 1024.f ;
            }

            WavCoeffs_ab[dir][i] *= alphaC;
        }
    }

    //Curve chro

    float cpMulC = cp.mulC[level];

    //  if( (cp.curv || cp.CHSLmet==1) && cp.CHmet!=2 && level < 9 && cpMulC != 0.f) { // cpMulC == 0.f means all will be multiplied by 1.f, so we can skip
    if( cp.CHmet != 2 && level < 9 && cpMulC != 0.f  && cp.chromena) { // cpMulC == 0.f means all will be multiplied by 1.f, so we can skip
        const float skinprot = params->wavelet.skinprotect;
        const float skinprotneg = -skinprot;
        const float factorHard = (1.f - skinprotneg / 100.f);
        bool useSkinControl = (skinprot != 0.f);

        for (int i = 0; i < W_ab * H_ab; i++) {
            int ii = i / W_ab;
            int jj = i - ii * W_ab;
            //WL and W_ab are identical
            float scale = 1.f;
            float modchro = varchrom[ii * 2][jj * 2];

            if(useSkinControl) {
                // hue chroma skin with initial lab datas
                float LL100 = labco->L[ii * 2][jj * 2] / 327.68f;
                float modhue = varhue[ii][jj];

                if(skinprot > 0.f) {
                    Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprot, scale, true, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 1); //1 for curve
                } else if(skinprot < 0.f) {
                    Color::SkinSatCbdl2 (LL100, modhue, modchro, skinprotneg, scale, false, cp.b_l, cp.t_l, cp.t_r, cp.b_r, 1);
                    scale = (scale == 1.f) ? factorHard : 1.f;
                }
            }

            float beta = (1024.f + 20.f * cpMulC * scale) / 1024.f ;

            if(beta < 0.02f) {
                beta = 0.02f;
            }

            float kClev = beta;

            if(cp.CHmet == 1) {
                if(level < cp.chrom) {
                    //linear for saturated
                    if((modchro > cp.t_lsat && modchro < cp.t_rsat)) {
                        kClev = beta;
                    } else if((modchro > cp.b_lsat && modchro <= cp.t_lsat)) {
                        float aaal = (1.f - beta) / (cp.b_lsat - cp.t_lsat);
                        float bbal = 1.f - aaal * cp.b_lsat;
                        kClev = aaal * modchro + bbal;
                    } else if((modchro > cp.t_rsat &&  modchro <= cp.b_rsat)) {
                        float aaar = (beta - 1.f) / (cp.t_rsat - cp.b_rsat);
                        float bbbr = 1.f - cp.b_rsat * aaar;
                        kClev = aaar * modchro + bbbr;
                    } else {
                        kClev = 1.f;
                    }
                } else {
                    //linear for pastel
                    if((modchro > cp.t_lpast && modchro < cp.t_rpast)) {
                        kClev = beta;
                    } else if((modchro > cp.b_lpast && modchro <= cp.t_lpast)) {
                        float aaalS = (1.f - beta) / (cp.b_lpast - cp.t_lpast);
                        float bbalS = 1.f - aaalS * cp.b_lpast;
                        kClev = aaalS * modchro + bbalS;
                    } else if((modchro > cp.t_rpast &&  modchro <= cp.b_rpast)) {
                        float aaarS = (beta - 1.f) / (cp.t_rpast - cp.b_rpast);
                        float bbbrS = 1.f - cp.b_rpast * aaarS;
                        kClev = aaarS * modchro + bbbrS;
                    } else {
                        kClev = 1.f;
                    }
                }
            } else if(cp.CHmet == 0) {
                kClev = beta;
            }

            WavCoeffs_ab[dir][i] *= kClev;
        }
    }

    bool useOpacity;
    float mulOpacity;

    if(useChannelA) {
        useOpacity = cp.opaRG;
        mulOpacity = cp.mulopaRG[level];
    } else {
        useOpacity = cp.opaBY;
        mulOpacity = cp.mulopaBY[level];
    }

    if((useOpacity && level < 9 && mulOpacity != 0.f) && cp.toningena) { //toning

        float beta = (1024.f + 20.f * mulOpacity) / 1024.f ;

        //float beta = (1000.f * mulOpacity);
        for (int i = 0; i < W_ab * H_ab; i++) {
            WavCoeffs_ab[dir][i] *= beta;
        }

        //  WavCoeffs_ab[dir][i] += beta;
    }

    if(waOpacityCurveW) {
        cp.opaW = true;
    }

    if(cp.bam  && cp.diag) {
//printf("OK Chroma\n");
        if(cp.opaW && cp.BAmet == 2) {
            int iteration = cp.ite;
            int itplus = 7 + iteration;
            int itmoins = 7 - iteration;
            int med = maxlvl / 2;
            int it;

            if(level < med) {
                it = itmoins;
            } else if(level == med) {
                it = 7;
            } else /*if(level > med)*/ {
                it = itplus;
            }

            for(int j = 0; j < it; j++) {
                //float bal = cp.balan;//-100 +100
                float kba = 1.f;

                //  if(dir <3) kba= 1.f + bal/600.f;
                //  if(dir==3) kba = 1.f - bal/300.f;
                for (int i = 0; i < W_ab * H_ab; i++) {
                    int ii = i / W_ab;
                    int jj = i - ii * W_ab;
                    float LL100 = labco->L[ii * 2][jj * 2] / 327.68f;
                    float k1 = 0.3f * (waOpacityCurveW[6.f * LL100] - 0.5f); //k1 between 0 and 0.5    0.5==> 1/6=0.16
                    float k2 = k1 * 2.f;

                    if(dir < 3) {
                        kba = 1.f + k1;
                    }

                    if(dir == 3) {
                        kba = 1.f - k2;
                    }

                    WavCoeffs_ab[dir][i] *= (kba);
                }
            }
        }

        if(cp.BAmet == 1) {
            int iteration = cp.ite;
            int itplus = 7 + iteration;
            int itmoins = 7 - iteration;
            int med = maxlvl / 2;
            int it;

            if(level < med) {
                it = itmoins;
            } else if(level == med) {
                it = 7;
            } else /*if(level > med)*/ {
                it = itplus;
            }

            for(int j = 0; j < it; j++) {
                float bal = cp.balan;//-100 +100
                float kba = 1.f;

                //  if(dir <3) kba= 1.f + bal/600.f;
                //  if(dir==3) kba = 1.f - bal/300.f;
                for (int i = 0; i < W_ab * H_ab; i++) {
                    int ii = i / W_ab;
                    int jj = i - ii * W_ab;
                    float k1 = 600.f;
                    float k2 = 300.f;
                    float LL100 = labco->L[ii * 2][jj * 2] / 327.68f;
                    float aa = 4970.f;
                    float bb = -397000.f;
                    float b0 = 100000.f;
                    float a0 = -4970.f;

                    if(LL100 > 80.f) {
                        k1 = aa * LL100 + bb;
                        k2 = 0.5f * k1;
                    }

                    if(LL100 < 20.f) {
                        k1 = a0 * LL100 + b0;
                        k2 = 0.5f * k1;
                    }

                    //k1=600.f;
                    //k2=300.f;
                    //k1=0.3f*(waOpacityCurveW[6.f*LL100]-0.5f);//k1 between 0 and 0.5    0.5==> 1/6=0.16
                    //k2=k1*2.f;
                    if(dir < 3) {
                        kba = 1.f + bal / k1;
                    }

                    if(dir == 3) {
                        kba = 1.f - bal / k2;
                    }

                    WavCoeffs_ab[dir][i] *= (kba);
                }
            }
        }

    }

    // to see each level of wavelet ...level from 0 to 8
    int choicelevel = atoi(params->wavelet.Lmethod.data()) - 1;
    choicelevel = choicelevel == -1 ? 4 : choicelevel;
    int choiceClevel = 0;

    if(params->wavelet.CLmethod == "one") {
        choiceClevel = 0;
    } else if(params->wavelet.CLmethod == "inf") {
        choiceClevel = 1;
    } else if(params->wavelet.CLmethod == "sup") {
        choiceClevel = 2;
    } else if(params->wavelet.CLmethod == "all") {
        choiceClevel = 3;
    }

    int choiceDir = 0;

    if(params->wavelet.Dirmethod == "one") {
        choiceDir = 1;
    } else if(params->wavelet.Dirmethod == "two") {
        choiceDir = 2;
    } else if(params->wavelet.Dirmethod == "thr") {
        choiceDir = 3;
    } else if(params->wavelet.Dirmethod == "all") {
        choiceDir = 0;
    }

    int dir1 = (choiceDir == 2) ? 1 : 2;
    int dir2 = (choiceDir == 3) ? 1 : 3;

    if(choiceClevel < 3) { // not all levels visible, paint residual
        if(level == 0) {
            if(cp.backm != 2) { // nothing to change when residual is used as background
                for (int i = 0; i < W_ab * H_ab; i++) {
                    WavCoeffs_ab0[i] = 0.f;
                }
            }
        }
    }

    if(choiceClevel == 0) { // Only one level
        if(choiceDir == 0) { // All directions
            if(level != choicelevel) { // zero all for the levels != choicelevel
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_ab * H_ab; i++) {
                        WavCoeffs_ab[dir][i] = 0.f;
                    }
                }
            }
        } else { // zero the unwanted directions for level == choicelevel
            if(choicelevel >= cp.maxilev) {
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_ab * H_ab; i++) {
                        WavCoeffs_ab[dir][i] = 0.f;
                    }
                }
            } else if(level != choicelevel) { // zero all for the levels != choicelevel
                for (int i = 0; i < W_ab * H_ab; i++) {
                    WavCoeffs_ab[dir1][i] = WavCoeffs_ab[dir2][i] = 0.f;
                }
            }
        }
    } else if(choiceClevel == 1) { // Only below level
        if(choiceDir == 0) { // All directions
            if(level > choicelevel) {
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_ab * H_ab; i++) {
                        WavCoeffs_ab[dir][i] = 0.f;
                    }
                }
            }
        } else { // zero the unwanted directions for level >= choicelevel
            if(level > choicelevel) {
                for (int i = 0; i < W_ab * H_ab; i++) {
                    WavCoeffs_ab[dir1][i] = WavCoeffs_ab[dir2][i] = 0.f;
                }
            }
        }
    } else if(choiceClevel == 2) { // Only above level
        if(choiceDir == 0) { // All directions
            if(level <= choicelevel) {
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_ab * H_ab; i++) {
                        WavCoeffs_ab[dir][i] = 0.f;
                    }
                }
            }
        } else { // zero the unwanted directions for level >= choicelevel
            if(choicelevel >= cp.maxilev) {
                for (int dir = 1; dir < 4; dir++) {
                    for (int i = 0; i < W_ab * H_ab; i++) {
                        WavCoeffs_ab[dir][i] = 0.f;
                    }
                }
            } else if(level <= choicelevel) {
                for (int i = 0; i < W_ab * H_ab; i++) {
                    WavCoeffs_ab[dir1][i] = WavCoeffs_ab[dir2][i] = 0.f;
                }
            }
        }
    }
}
}
