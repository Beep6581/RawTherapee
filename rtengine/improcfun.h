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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <memory>

#include "coord2d.h"
#include "gamutwarning.h"
#include "jaggedarray.h"
#include "pipettebuffer.h"
#include "array2D.h"

#include <cairomm/cairomm.h>

template<typename T>
class LUT;

using LUTu = LUT<uint32_t>;
using LUTf = LUT<float>;

template<typename T, const size_t num>
class multi_array2D;
namespace rtengine
{

class ColorAppearance;
class ColorGradientCurve;
class DCPProfile;
class DCPProfileApplyState;
class FlatCurve;
class FramesMetaData;
class LensCorrection;
class LocCCmaskCurve;
class LocLLmaskCurve;
class LocHHmaskCurve;
class LocwavCurve;
class LocretigainCurve;
class LocretitransCurve;
class LocLHCurve;
class LocHHCurve;
class NoiseCurve;
class OpacityCurve;
class ToneCurve;
class WavCurve;
class WavOpacityCurveBY;
class WavOpacityCurveRG;
class WavOpacityCurveW;
class WavOpacityCurveWL;

class CieImage;
class Image8;
class Imagefloat;
class LabImage;
class wavelet_decomposition;

namespace procparams
{

class ProcParams;

struct DehazeParams;
struct FattalToneMappingParams;
struct ColorManagementParams;
struct DirPyrDenoiseParams;
struct LocalContrastParams;
struct LocallabParams;
struct SharpeningParams;
struct SoftLightParams;
struct VibranceParams;
struct VignettingParams;
struct WaveletParams;

}

enum RenderingIntent : int;

class ImProcFunctions
{
    cmsHTRANSFORM monitorTransform;
    std::unique_ptr<GamutWarning> gamutWarning;
    Cairo::RefPtr<Cairo::ImageSurface> locImage;

    const procparams::ProcParams* params;
    double scale;
    bool multiThread;

    bool lastcutpast;
    int lastcxbuf;
    int lastcybuf;
    int lastcount;
    LabImage *spotbuffer;

    void calcVignettingParams(int oW, int oH, const procparams::VignettingParams& vignetting, double &w2, double &h2, double& maxRadius, double &v, double &b, double &mul);

    void transformLuminanceOnly(Imagefloat* original, Imagefloat* transformed, int cx, int cy, int oW, int oH, int fW, int fH);
    void transformGeneral(bool highQuality, Imagefloat *original, Imagefloat *transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH, const LensCorrection *pLCPMap);
    void transformLCPCAOnly(Imagefloat *original, Imagefloat *transformed, int cx, int cy, const LensCorrection *pLCPMap);

    bool needsCA() const;
    bool needsDistortion() const;
    bool needsRotation() const;
    bool needsPerspective() const;
    bool needsGradient() const;
    bool needsVignetting() const;
    bool needsLCP() const;
    bool needsLensfun() const;
//   static cmsUInt8Number* Mempro = NULL;

public:
    enum class Median {
        TYPE_3X3_SOFT,
        TYPE_3X3_STRONG,
        TYPE_5X5_SOFT,
        TYPE_5X5_STRONG,
        TYPE_7X7,
        TYPE_9X9
    };

    double lumimul[3];

    explicit ImProcFunctions(const procparams::ProcParams* iparams, bool imultiThread = true)
        : monitorTransform(nullptr), params(iparams), scale(1), multiThread(imultiThread), lumimul{} {}
    ~ImProcFunctions();
    bool needsLuminanceOnly()
    {
        return !(needsCA() || needsDistortion() || needsRotation() || needsPerspective() || needsLCP() || needsLensfun()) && (needsVignetting() || needsPCVignetting() || needsGradient());
    }
    void setScale(double iscale);

    bool needsTransform() const;
    bool needsPCVignetting() const;
    float calcGradientFactor (const struct grad_params& gp, int x, int y);
    void firstAnalysis(const Imagefloat* const working, const procparams::ProcParams &params, LUTu & vhist16);
    void updateColorProfiles(const Glib::ustring& monitorProfile, RenderingIntent monitorIntent, bool softProof, bool gamutCheck);
    void rgbProc(Imagefloat* working, LabImage* lab, PipetteBuffer *pipetteBuffer, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                 int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve, float satLimit, float satLimitOpacity, const ColorGradientCurve & ctColorCurve, const OpacityCurve & ctOpacityCurve, bool opautili, LUTf & clcurve, LUTf & cl2curve, const ToneCurve & customToneCurve1, const ToneCurve & customToneCurve2,
                 const ToneCurve & customToneCurvebw1, const ToneCurve & customToneCurvebw2, double &rrm, double &ggm, double &bbm, float &autor, float &autog, float &autob, DCPProfile *dcpProf, const DCPProfileApplyState &asIn, LUTu &histToneCurve, size_t chunkSize = 1, bool measure = false);
    void rgbProc(Imagefloat* working, LabImage* lab, PipetteBuffer *pipetteBuffer, LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve,
                 int sat, LUTf & rCurve, LUTf & gCurve, LUTf & bCurve, float satLimit, float satLimitOpacity, const ColorGradientCurve & ctColorCurve, const OpacityCurve & ctOpacityCurve, bool opautili, LUTf & clcurve, LUTf & cl2curve, const ToneCurve & customToneCurve1, const ToneCurve & customToneCurve2,
                 const ToneCurve & customToneCurvebw1, const ToneCurve & customToneCurvebw2, double &rrm, double &ggm, double &bbm, float &autor, float &autog, float &autob,
                 double expcomp, int hlcompr, int hlcomprthresh, DCPProfile *dcpProf, const DCPProfileApplyState &asIn, LUTu &histToneCurve, size_t chunkSize = 1, bool measure = false);
    void labtoning(float r, float g, float b, float &ro, float &go, float &bo, int algm, int metchrom, int twoc, float satLimit, float satLimitOpacity, const ColorGradientCurve & ctColorCurve, const OpacityCurve & ctOpacityCurve, LUTf & clToningcurve, LUTf & cl2Toningcurve, float iplow, float iphigh, double wp[3][3], double wip[3][3]);
    void toning2col(float r, float g, float b, float &ro, float &go, float &bo, float iplow, float iphigh, float rl, float gl, float bl, float rh, float gh, float bh, float SatLow, float SatHigh, float balanS, float balanH, float reducac, int mode, int preser, float strProtect);
    void toningsmh(float r, float g, float b, float &ro, float &go, float &bo, float RedLow, float GreenLow, float BlueLow, float RedMed, float GreenMed, float BlueMed, float RedHigh, float GreenHigh, float BlueHigh, float reducac, int mode, float strProtect);
    void toningsmh2(float r, float g, float b, float &ro, float &go, float &bo, float low[3], float satLow, float med[3], float satMed, float high[3], float satHigh, float reducac, int mode, int preser);
    void secondeg_begin(float reducac, float vend, float &aam, float &bbm);
    void secondeg_end(float reducac, float vinf, float &aa, float &bb, float &cc);

    void retreavergb(float &r, float &g, float &b);
    void moyeqt(Imagefloat* working, float &moyS, float &eqty);

    void luminanceCurve(LabImage* lold, LabImage* lnew, LUTf &curve);
    void ciecamloc_02float(int sp, LabImage* lab);

    void ciecam_02float(CieImage* ncie, float adap, int pW, int pwb, LabImage* lab, const procparams::ProcParams* params,
                        const ColorAppearance & customColCurve1, const ColorAppearance & customColCurve, const ColorAppearance & customColCurve3,
                        LUTu &histLCAM, LUTu &histCCAM, LUTf & CAMBrightCurveJ, LUTf & CAMBrightCurveQ, float &mean, int Iterates, int scale, bool execsharp, float &d, float &dj, float &yb, int rtt,
                        bool showSharpMask = false);
    void chromiLuminanceCurve(PipetteBuffer *pipetteBuffer, int pW, LabImage* lold, LabImage* lnew, LUTf &acurve, LUTf &bcurve, LUTf & satcurve, LUTf & satclcurve, LUTf &clcurve, LUTf &curve, bool utili, bool autili, bool butili, bool ccutili, bool cclutili, bool clcutili, LUTu &histCCurve, LUTu &histLurve);
    void vibrance(LabImage* lab, const procparams::VibranceParams &vibranceParams, bool highlight, const Glib::ustring &workingProfile);         //Jacques' vibrance
    void softprocess(const LabImage* bufcolorig, array2D<float> &buflight, /* float ** bufchro, float ** buf_a, float ** buf_b, */ float rad, int bfh, int bfw, double epsilmax, double epsilmin,  float thres, int sk, bool multiThread);
    void softproc(const LabImage* bufcolorig, const LabImage* bufcolfin, float rad, int bfh, int bfw, double epsilmax, double epsilmin, float thres, int sk, bool multiThread, int flag);
//    void colorCurve       (LabImage* lold, LabImage* lnew);
    void sharpening(LabImage* lab, const procparams::SharpeningParams &sharpenParam, bool showMask = false);
    void sharpeningcam(CieImage* ncie, float** buffer, bool showMask = false);
    void transform(Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH, const FramesMetaData *metadata, int rawRotationDeg, bool fullImage);
    float resizeScale(const procparams::ProcParams* params, int fw, int fh, int &imw, int &imh);
    void lab2monitorRgb(LabImage* lab, Image8* image);
    void resize(Imagefloat* src, Imagefloat* dst, float dScale);
    void Lanczos(const LabImage* src, LabImage* dst, float scale);
    void Lanczos(const Imagefloat* src, Imagefloat* dst, float scale);

    void deconvsharpening(float** luminance, float** buffer, const float* const * blend, int W, int H, const procparams::SharpeningParams &sharpenParam, double Scale);
    void deconvsharpeningloc(float** luminance, float** buffer, int W, int H, float** loctemp, int damp, double radi, int ite, int amo, int contrast, double blurrad);

    void MLsharpen(LabImage* lab); // Manuel's clarity / sharpening
    void MLmicrocontrast(float** luminance, int W, int H);   //Manuel's microcontrast
    void MLmicrocontrast(LabImage* lab);   //Manuel's microcontrast
    void MLmicrocontrastcam(CieImage* ncie);   //Manuel's microcontrast

    void impulsedenoise(LabImage* lab);   //Emil's impulse denoise
    void impulsedenoisecam(CieImage* ncie, float **buffers[3]);
    void impulse_nr(LabImage* lab, double thresh);
    void impulse_nrcam(CieImage* ncie, double thresh, float **buffers[3]);

    void dirpyrdenoise(LabImage* src);    //Emil's pyramid denoise
    void dirpyrequalizer(LabImage* lab, int scale);  //Emil's wavelet


    void EPDToneMapResid(float * WavCoeffs_L0, unsigned int Iterates,  int skip, struct cont_params& cp, int W_L, int H_L, float max0, float min0);
    void CompressDR(float *Source, int W_L, int H_L, float Compression, float DetailBoost);
    void ContrastResid(float * WavCoeffs_L0, struct cont_params &cp, int W_L, int H_L, float max0, float min0);

    void EPDToneMap(LabImage *lab, unsigned int Iterates = 0, int skip = 1);
    void EPDToneMaplocal(int sp, LabImage *lab, LabImage *tmp1, unsigned int Iterates, int skip);
    void EPDToneMapCIE(CieImage *ncie, float a_w, float c_, int Wid, int Hei, float minQ, float maxQ, unsigned int Iterates = 0, int skip = 1);

    // pyramid denoise
//    procparams::DirPyrDenoiseParams dnparams;
    void dirpyr(LabImage* data_fine, LabImage* data_coarse, int level, LUTf &rangefn_L, LUTf &rangefn_ab,
                int pitch, int scale, const int luma, int chroma); 
    void idirpyr(LabImage* data_coarse, LabImage* data_fine, int level, LUTf &rangefn_L, LUTf & nrwt_l, LUTf & nrwt_ab,
                 int pitch, int scale, const int luma, const int chroma/*, LUTf & Lcurve, LUTf & abcurve*/);
    //locallab
    void maskcalccol(bool invmask, bool pde, int bfw, int bfh, int xstart, int ystart, int sk, int cx, int cy, LabImage* bufcolorig, LabImage* bufmaskblurcol, LabImage* originalmaskcol, LabImage* original, int inv, const struct local_params & lp,
                 float strumask, bool astool,
                 const LocCCmaskCurve & locccmasCurve, bool & lcmasutili, 
                 const LocLLmaskCurve & locllmasCurve, bool & llmasutili, 
                 const LocHHmaskCurve & lochhmasCurve, bool &lhmasutili, const LocHHmaskCurve & lochhhmasCurve, bool &lhhmasutili, 
                 bool multiThread, bool enaMask, bool showmaske, bool deltaE, bool modmask, bool zero, bool modif, float chrom, float rad, float lap, float gamma, float slope, float blendm, int shado, float amountcd, float anchorcd,
                 LUTf & lmasklocalcurve, bool & localmaskutili,
                 const LocwavCurve & loclmasCurvecolwav, bool & lmasutilicolwav, int level_bl, int level_hl, int level_br, int level_hr,
                 int shortcu, bool delt, const float hueref, const float chromaref, const float lumaref,
                 float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope);
                 
    void deltaEforMask(float **rdE, int bfw, int bfh, LabImage* bufcolorig, const float hueref, const float chromaref, const float lumaref,
                          float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope, float balance);
    void discrete_laplacian_threshold(float * data_out, const float * data_in, size_t nx, size_t ny, float t);
    void rex_poisson_dct(float * data, size_t nx, size_t ny, double m);
    void mean_dt(const float * data, size_t size, double * mean_p, double * dt_p);
    double *cos_table(size_t size);

    void normalize_mean_dt(float *data, const float *ref, size_t size, float mod);
    void retinex_pde(float *datain, float * dataout, int bfw, int bfh, float thresh, float multy, float *dE, int show, int dEenable, int normalize);
    void exposure_pde(float *dataor, float *datain, float * dataout, int bfw, int bfh, float thresh, float mod);
    void fftw_convol_blur(float *input, float *output, int bfw, int bfh, float radius, int fftkern, int algo);
    void fftw_convol_blur2(float **input2, float **output2, int bfw, int bfh, float radius, int fftkern, int algo);
    void fftw_tile_blur(int GW, int GH, int tilssize , int max_numblox_W, int min_numblox_W, float **tmp1, int numThreads, double radius);

    void maskforretinex(int sp, int before, float ** luminance, float ** out, int W_L, int H_L, int skip,
         const LocCCmaskCurve & locccmasretiCurve, bool &lcmasretiutili, const  LocLLmaskCurve & locllmasretiCurve, bool &llmasretiutili, const  LocHHmaskCurve & lochhmasretiCurve, bool & lhmasretiutili,
         int llretiMask, bool retiMasktmap, bool retiMask, float rad, float lap, bool pde, float gamm, float slop, float chro, float blend,
         LUTf & lmaskretilocalcurve, bool & localmaskretiutili,
         LabImage * bufreti, LabImage * bufmask, LabImage * buforig, LabImage * buforigmas, bool multiThread,
         bool delt, const float hueref, const float chromaref, const float lumaref,
         float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope, float balance, float lumask);

    void filmGrain(Imagefloat *rgb, int isogr, int strengr, int scalegr, int bfw, int bfh);

    void MSRLocal(int sp, bool fftw, int lum, float** reducDE, LabImage * bufreti, LabImage * bufmask, LabImage * buforig, LabImage * buforigmas, float** luminance, float** templ, const float* const *originalLuminance,
        const int width, const int height, int bfwr, int bfhr, const procparams::LocallabParams &loc, const int skip, const LocretigainCurve &locRETgainCcurve, const LocretitransCurve &locRETtransCcurve,
        const int chrome, const int scall, const float krad, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax,
        const LocCCmaskCurve & locccmasretiCurve, bool &lcmasretiutili, const  LocLLmaskCurve & locllmasretiCurve, bool &llmasretiutili, const  LocHHmaskCurve & lochhmasretiCurve, bool & lhmasretiutili, int llretiMask,
        LUTf & lmaskretilocalcurve, bool & localmaskretiutili,
        LabImage * transformed, bool retiMasktmap, bool retiMask,
        bool delt, const float hueref, const float chromaref, const float lumaref,
        float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope, float balance, float lumask);


    void calc_ref(int sp, LabImage* original, LabImage* transformed, int cx, int cy, int oW, int oH, int sk, double &huerefblur, double &chromarefblur, double &lumarefblur, double &hueref, double &chromaref, double &lumaref, double &sobelref, float &avg);
    void copy_ref(LabImage* spotbuffer, LabImage* original, LabImage* transformed, int cx, int cy, int sk, const struct local_params & lp, double &huerefspot, double &chromarefspot, double &lumarefspot);
    void paste_ref(LabImage* spotbuffer, LabImage* transformed, int cx, int cy, int sk, const struct local_params & lp);
    void Lab_Local(int call, int sp, float** shbuffer, LabImage* original, LabImage* transformed, LabImage* reserved, LabImage* lastorig, int cx, int cy, int oW, int oH, int sk, const LocretigainCurve & locRETgainCcurve, const LocretitransCurve &locRETtransCcurve, LUTf & lllocalcurve,  bool & locallutili, const LocLHCurve & loclhCurve,  const LocHHCurve & lochhCurve,
                LUTf & lmasklocalcurve, bool & localmaskutili,
                LUTf & lmaskexplocalcurve, bool & localmaskexputili,
                LUTf & lmaskSHlocalcurve, bool & localmaskSHutili,
                LUTf & lmaskviblocalcurve, bool & localmaskvibutili,
                LUTf & lmasktmlocalcurve, bool & localmasktmutili,
                LUTf & lmaskretilocalcurve, bool & localmaskretiutili,
                LUTf & lmaskcblocalcurve, bool & localmaskcbutili,
                LUTf & lmaskbllocalcurve, bool & localmaskblutili,
                const LocCCmaskCurve & locccmasCurve, bool & lcmasutili, const  LocLLmaskCurve & locllmasCurve, bool & llmasutili, const  LocHHmaskCurve & lochhmasCurve, bool & lhmasutili, const  LocHHmaskCurve & lochhhmasCurve, bool & lhhmasutili,
                const LocCCmaskCurve & locccmasexpCurve, bool &lcmasexputili, const  LocLLmaskCurve & locllmasexpCurve, bool &llmasexputili, const  LocHHmaskCurve & lochhmasexpCurve, bool & lhmasexputili, 
                const LocCCmaskCurve & locccmasSHCurve, bool &lcmasSHutili, const  LocLLmaskCurve & locllmasSHCurve, bool &llmasSHutili, const  LocHHmaskCurve & lochhmasSHCurve, bool & lhmasSHutili,
                const LocCCmaskCurve & locccmasvibCurve, bool &lcmasvibutili, const  LocLLmaskCurve & locllmasvibCurve, bool &llmasvibutili, const  LocHHmaskCurve & lochhmasvibCurve, bool & lhmasvibutili,
                const LocCCmaskCurve & locccmascbCurve, bool &lcmascbutili, const  LocLLmaskCurve & locllmascbCurve, bool &llmascbutili, const  LocHHmaskCurve & lochhmascbCurve, bool & lhmascbutili,
                const LocCCmaskCurve & locccmasretiCurve, bool &lcmasretiutili, const  LocLLmaskCurve & locllmasretiCurve, bool &llmasretiutili, const  LocHHmaskCurve & lochhmasretiCurve, bool & lhmasretiutili,
                const LocCCmaskCurve & locccmastmCurve, bool &lcmastmutili, const  LocLLmaskCurve & locllmastmCurve, bool &llmastmutili, const  LocHHmaskCurve & lochhmastmCurve, bool & lhmastmutili,
                const LocCCmaskCurve & locccmasblCurve, bool &lcmasblutili, const  LocLLmaskCurve & locllmasblCurve, bool &llmasblutili, const  LocHHmaskCurve & lochhmasblCurve, bool & lhmasblutili,
                const LocwavCurve & loclmasCurveblwav, bool & lmasutiliblwav,
                const LocwavCurve & loclmasCurvecolwav, bool & lmasutilicolwav,
                const LocwavCurve & locwavCurve, bool & locwavutili,
                bool &LHutili, bool &HHutili, LUTf & cclocalcurve, bool & localcutili, LUTf & rgblocalcurve, bool & localrgbutili, bool & localexutili, LUTf & exlocalcurve, LUTf & hltonecurveloc, LUTf & shtonecurveloc, LUTf & tonecurveloc, LUTf & lightCurveloc,
                double & huerefblur, double &chromarefblur, double & lumarefblur, double &hueref, double &chromaref, double &lumaref, double &sobelref, int &lastsav,
                int llColorMask, int llColorMaskinv, int llExpMask, int llExpMaskinv, int llSHMask, int llSHMaskinv, int llvibMask, int llcbMask, int llretiMask, int llsoftMask, int lltmMask, int llblMask,
                float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax);

    void addGaNoise(LabImage *lab, LabImage *dst, const float mean, const float variance, const int sk);
    void BlurNoise_Localold(int call, const struct local_params& lp, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy);
    void InverseBlurNoise_Local(LabImage * originalmask, float **bufchro, const struct local_params& lp, const float hueref, const float chromaref,  const float lumaref, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy, int sk);
    void InverseReti_Local(const struct local_params& lp, const float hueref, const float chromaref,  const float lumaref, LabImage* original, LabImage* transformed, const LabImage* const tmp1, int cx, int cy, int chro, int sk);
    void BlurNoise_Local(LabImage* tmp1, LabImage * originalmask, float **bufchro, const float hueref, const float chromaref, const float lumaref, const local_params& lp, LabImage* original, LabImage* transformed, int cx, int cy, int sk);
    static void strcurv_data(std::string retistr, int *s_datc, int &siz);
    void blendstruc(int bfw, int bfh, LabImage* bufcolorig, float radius, float stru, array2D<float> & blend2, int sk, bool multiThread);

    void wavcontrast4(float ** tmp, float contrast, int bfw, int bfh, int level_bl, int level_hl, int level_br, int level_hr, int sk, bool numThreads, const LocwavCurve & locwavCurve, bool & locwavutili, int & maxlvl);

    void transit_shapedetect_retinex(int call, int senstype, LabImage * bufexporig, LabImage * bufmask, LabImage * buforigmas, float **buflight, float **bufchro, const float hueref, const float chromaref,  const float lumaref, const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk);
    void transit_shapedetect(int senstype, const LabImage *bufexporig, const LabImage *bufcolfin, LabImage * originalmask, float **buflight, float **bufchro, float **buf_a_cat, float ** buf_b_cat, float ** bufhh, bool HHutili, const float hueref, const float chromaref, const float lumaref, float sobelref, float meansobel, float ** blend2, const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk);
    void exlabLocal(local_params& lp, int bfh, int bfw, LabImage* bufexporig, LabImage* lab,  LUTf & hltonecurve, LUTf & shtonecurve, LUTf & tonecurve, float mean);
    void Exclude_Local(float **deltaso, float hueref, float chromaref, float lumaref, float sobelref, float meansobel, const struct local_params & lp, const LabImage * original, LabImage * transformed, const LabImage * rsv, const LabImage * reserv, int cx, int cy, int sk);

    void DeNoise_Local(int call, const struct local_params& lp, LabImage * originalmask, int levred, float hueref, float lumaref, float chromaref, LabImage* original, LabImage* transformed, LabImage &tmp1, int cx, int cy, int sk);
    void DeNoise(int call, int del,  float * slidL, float * slida, float * slidb, int aut, bool noiscfactiv, const struct local_params& lp, LabImage* originalmaskbl, int levred, float huerefblur, float lumarefblur, float chromarefblur, LabImage* original, LabImage* transformed, int cx, int cy, int sk);


    void fftw_denoise(int GW, int GH, int max_numblox_W, int min_numblox_W, float **tmp1, array2D<float> *Lin,  int numThreads, const struct local_params & lp, int chrom);

    void ColorLight_Local(float moddE, float powdE, int call, LabImage * bufcolorig, LabImage * originalmask, float **buflight, float **bufchro, float **bufchroslid, float ** bufhh, float ** buflightslid, bool &LHutili, bool &HHutili, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const float lumaref, float sobelref, float ** blend2, LUTf & lllocalcurve, const LocLHCurve & loclhCurve, const LocHHCurve & lochhCurve, LUTf & lightCurveloc, const local_params& lp, LabImage* original, LabImage* transformed, int cx, int cy, int sk);
    void InverseColorLight_Local(bool tonequ, bool tonecurv, int sp, int senstype, struct local_params& lp, LabImage * originalmask, LUTf & lightCurveloc, LUTf & hltonecurveloc, LUTf & shtonecurveloc, LUTf & tonecurveloc, LUTf & exlocalcurve, LUTf & cclocalcurve, float adjustr, bool localcutili, LUTf & lllocalcurve, bool locallutili, LabImage* original, LabImage* transformed, int cx, int cy, const float hueref, const float chromaref, const float lumaref, int sk);
    void Sharp_Local(int call, float **loctemp,  int senstype, const float hueref,  const float chromaref, const float lumaref, const local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk);

//    void Sharp_Local(int call, float **loctemp,  int senstype, const float hueplus, const float huemoins, const float hueref, const float dhue, const float chromaref, const local_params& lp, LabImage* original, LabImage* transformed, int cx, int cy, int sk);
    void InverseSharp_Local(float **loctemp, const float hueref, const float lumaref, const float chromaref, const local_params& lp, LabImage* original, LabImage* transformed, int cx, int cy, int sk);


    void Tile_calc(int tilesize, int overlap, int kall, int imwidth, int imheight, int &numtiles_W, int &numtiles_H, int &tilewidth, int &tileheight, int &tileWskip, int &tileHskip);
    void ip_wavelet(LabImage * lab, LabImage * dst, int kall, const procparams::WaveletParams & waparams, const WavCurve & wavCLVCcurve, const WavOpacityCurveRG & waOpacityCurveRG, const WavOpacityCurveBY & waOpacityCurveBY,  const WavOpacityCurveW & waOpacityCurveW, const WavOpacityCurveWL & waOpacityCurveWL, LUTf &wavclCurve, int skip);

    void WaveletcontAllL(LabImage * lab, float **varhue, float **varchrom, wavelet_decomposition &WaveletCoeffs_L,
                         struct cont_params &cp, int skip, float *mean, float *sigma, float *MaxP, float *MaxN,  const WavCurve & wavCLVCcurve, const WavOpacityCurveW & waOpacityCurveW, FlatCurve* ChCurve, bool Chutili);
    void WaveletcontAllLfinal(wavelet_decomposition &WaveletCoeffs_L, struct cont_params &cp, float *mean, float *sigma, float *MaxP, const WavOpacityCurveWL & waOpacityCurveWL);
    void WaveletcontAllAB(LabImage * lab, float **varhue, float **varchrom, wavelet_decomposition &WaveletCoeffs_a, const WavOpacityCurveW & waOpacityCurveW,
                          struct cont_params &cp, const bool useChannelA);
    void WaveletAandBAllAB(wavelet_decomposition &WaveletCoeffs_a, wavelet_decomposition &WaveletCoeffs_b,
                           struct cont_params &cp, FlatCurve* hhcurve, bool hhutili);
    void ContAllL(float **koeLi, float *maxkoeLi, bool lipschitz, int maxlvl, LabImage * lab, float **varhue, float **varchrom, float ** WavCoeffs_L, float * WavCoeffs_L0, int level, int dir, struct cont_params &cp,
                  int W_L, int H_L, int skip, float *mean, float *sigma, float *MaxP, float *MaxN,  const WavCurve & wavCLVCcurve, const WavOpacityCurveW & waOpacityCurveW, FlatCurve* ChCurve, bool Chutili);
    void finalContAllL(float ** WavCoeffs_L, float * WavCoeffs_L0, int level, int dir, struct cont_params &cp,
                       int W_L, int H_L, float *mean, float *sigma, float *MaxP, const WavOpacityCurveWL & waOpacityCurveWL);
    void ContAllAB(LabImage * lab, int maxlvl, float **varhue, float **varchrom, float ** WavCoeffs_a, float * WavCoeffs_a0, int level, int dir, const WavOpacityCurveW & waOpacityCurveW, struct cont_params &cp,
                   int W_ab, int H_ab, const bool useChannelA);
    void Evaluate2(wavelet_decomposition &WaveletCoeffs_L,
                   float *mean, float *meanN, float *sigma, float *sigmaN, float *MaxP, float *MaxN);
    void Eval2(float ** WavCoeffs_L, int level,
               int W_L, int H_L, float *mean, float *meanN, float *sigma, float *sigmaN, float *MaxP, float *MaxN);

    void Aver(float * HH_Coeffs, int datalen, float &averagePlus, float &averageNeg, float &max, float &min);
    void Sigma(float * HH_Coeffs, int datalen, float averagePlus, float averageNeg, float &sigmaPlus, float &sigmaNeg);
    void calckoe(float ** WavCoeffs_LL, const struct cont_params& cp, float ** koeLi, int level, int dir, int W_L, int H_L, float edd, float *maxkoeLi, float **tmC = nullptr);



    void Median_Denoise(float **src, float **dst, int width, int height, Median medianType, int iterations, int numThreads, float **buffer = nullptr);
    void Median_Denoise(float **src, float **dst, float upperBound, int width, int height, Median medianType, int iterations, int numThreads, float **buffer = nullptr);
    void RGB_denoise(int kall, Imagefloat * src, Imagefloat * dst, Imagefloat * calclum, float * ch_M, float *max_r, float *max_b, bool isRAW, const procparams::DirPyrDenoiseParams & dnparams, const double expcomp, const NoiseCurve & noiseLCurve, const NoiseCurve & noiseCCurve, float &nresi, float &highresi);
    void RGB_denoise_infoGamCurve(const procparams::DirPyrDenoiseParams & dnparams, const bool isRAW, LUTf &gamcurve, float &gam, float &gamthresh, float &gamslope);
    void RGB_denoise_info(Imagefloat * src, Imagefloat * provicalc, bool isRAW, const LUTf &gamcurve, float gam, float gamthresh, float gamslope, const procparams::DirPyrDenoiseParams & dnparams, const double expcomp, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float & maxblueaut, float &minredaut, float & minblueaut, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc, bool multiThread = false);
    void RGBtile_denoise(float * fLblox, int hblproc, float noisevar_Ldetail);     //for DCT
    void RGBoutput_tile_row(float *bloxrow_L, float ** Ldetail, float ** tilemask_out, int height, int width, int top);
    bool WaveletDenoiseAllL(const wavelet_decomposition &WaveletCoeffs_L, float *noisevarlum, float madL[8][3], float * vari, int edge, int denoiseNestedLevels);
    bool WaveletDenoiseAllAB(const wavelet_decomposition &WaveletCoeffs_L, const wavelet_decomposition &WaveletCoeffs_ab, float *noisevarchrom, float madL[8][3], float *variC, int local, float noisevar_ab, const bool useNoiseCCurve, bool autoch, bool denoiseMethodRgb, int denoiseNestedLevels);
    bool WaveletDenoiseAll_BiShrinkAB(const wavelet_decomposition &WaveletCoeffs_L, const wavelet_decomposition &WaveletCoeffs_ab, float *noisevarchrom, float madL[8][3], float *variC, int local, float noisevar_ab, const bool useNoiseCCurve,  bool autoch, bool denoiseMethodRgb, int denoiseNestedLevels);

    void WaveletDenoiseAll_info(int levwav, const wavelet_decomposition &WaveletCoeffs_a,
                                const wavelet_decomposition &WaveletCoeffs_b, float **noisevarlum, float **noisevarchrom, float **noisevarhue, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float & minblueaut, int schoice, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc,
                                float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool denoiseMethodRgb);

    bool WaveletDenoiseAll_BiShrinkL(const wavelet_decomposition &WaveletCoeffs_L, float *noisevarlum, float madL[8][3], float * vari, int edge, int denoiseNestedLevels);
    void ShrinkAllL(const wavelet_decomposition &WaveletCoeffs_L, float **buffer, int level, int dir, float *noisevarlum, float * madL, float * vari, int edge);
    void ShrinkAllAB(const wavelet_decomposition &WaveletCoeffs_L, const wavelet_decomposition &WaveletCoeffs_ab, float **buffer, int level, int dir,
                     float *noisevarchrom, float noisevar_ab, const bool useNoiseCCurve, bool autoch, bool denoiseMethodRgb, float * madL, float * variC, int local, float * madaab = nullptr, bool madCalculated = false);
    void ShrinkAll_info(float ** WavCoeffs_a, float ** WavCoeffs_b,
                        int W_ab, int H_ab, float **noisevarlum, float **noisevarchrom, float **noisevarhue, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut, int schoice, int lvl, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc,
                        float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool denoiseMethodRgb);
    void Noise_residualAB(const wavelet_decomposition &WaveletCoeffs_ab, float &chresid, float &chmaxresid, bool denoiseMethodRgb);
    void calcautodn_info(float &chaut, float &delta, int Nb, int levaut, float maxmax, float lumema, float chromina, int mode, int lissage, float redyel, float skinc, float nsknc);
    float Mad(const float * DataList, int datalen);
    float MadRgb(const float * DataList, int datalen);

    // pyramid wavelet
    void cbdl_local_temp(float ** src, float ** loctemp, int srcwidth, int srcheight, const float * mult, float kchro, const double dirpyrThreshold, const float mergeL, const float contres, const float blurcb, const double skinprot, const bool gamutlab, float b_l, float t_l, float t_r, float b_r,  int choice, int scale, bool multiThread);
    void dirpyr_equalizer(const float * const * src, float ** dst, int srcwidth, int srcheight, const float * const * l_a, const float * const * l_b, const double * mult, double dirpyrThreshold, double skinprot, float b_l, float t_l, float t_r, int scale);    //Emil's directional pyramid wavelet
    void dirpyr_equalizercam(const CieImage* ncie, float ** src, float ** dst, int srcwidth, int srcheight, const float * const * h_p, const float * const * C_p,  const double * mult, const double dirpyrThreshold, const double skinprot, float b_l, float t_l, float t_r, int scale);    //Emil's directional pyramid wavelet
    void defringe(LabImage* lab);
    void defringecam(CieImage* ncie);
    void badpixcam(CieImage* ncie, double rad, int thr, int mode, float chrom, bool hotbad);
    void badpixlab(LabImage* lab, double rad, int thr, float chrom);

    void PF_correct_RT(LabImage * lab, double radius, int thresh);
    void PF_correct_RTcam(CieImage * ncie, double radius, int thresh);
    void Badpixelscam(CieImage * ncie, double radius, int thresh, int mode, float chrom, bool hotbad);
    void BadpixelsLab(LabImage * lab, double radius, int thresh, float chrom);

    void dehaze(Imagefloat *rgb, const procparams::DehazeParams &dehazeParams);
    void ToneMapFattal02(Imagefloat *rgb, const procparams::FattalToneMappingParams &fatParams, int detail_level);
    void localContrast(LabImage *lab, float **destination, const procparams::LocalContrastParams &localContrastParams, bool fftwlc, double scale);
    void colorToningLabGrid(LabImage *lab, int xstart, int xend, int ystart, int yend, bool MultiThread);
   // void shadowsHighlights(LabImage *lab);
    void shadowsHighlights(LabImage *lab, bool ena, int labmode, int hightli, int shado, int rad, int scal, int hltonal, int shtonal);
    
    void softLight(LabImage *lab, const procparams::SoftLightParams &softLightParams);
    void labColorCorrectionRegions(LabImage *lab);

    Image8*     lab2rgb(LabImage* lab, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm, bool consider_histogram_settings = true);
    Imagefloat*    lab2rgbOut(LabImage* lab, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm);
    // CieImage *ciec;
    void workingtrc(const Imagefloat* src, Imagefloat* dst, int cw, int ch, int mul, const Glib::ustring &profile, double gampos, double slpos, cmsHTRANSFORM &transform, bool normalizeIn = true, bool normalizeOut = true, bool keepTransForm = false) const;

    bool transCoord(int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef = -1, const LensCorrection *pLCPMap = nullptr) const;
    bool transCoord(int W, int H, const std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef = -1, const LensCorrection *pLCPMap = nullptr) const;
    static void getAutoExp(const LUTu & histogram, int histcompr, double clip, double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh);
    static double getAutoDistor(const Glib::ustring& fname, int thumb_size);
    double getTransformAutoFill(int oW, int oH, const LensCorrection *pLCPMap = nullptr) const;
    void rgb2lab(const Imagefloat &src, LabImage &dst, const Glib::ustring &workingSpace);
    void lab2rgb(const LabImage &src, Imagefloat &dst, const Glib::ustring &workingSpace);
};

}
