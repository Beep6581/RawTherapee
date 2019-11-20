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
#include "pipettebuffer.h"
#include "shmap.h"

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

struct ColorManagementParams;
struct DirPyrDenoiseParams;
struct SharpeningParams;
struct VignettingParams;
struct WaveletParams;

}

enum RenderingIntent : int;

class ImProcFunctions
{
    cmsHTRANSFORM monitorTransform;
    std::unique_ptr<GamutWarning> gamutWarning;

    const procparams::ProcParams* params;
    double scale;
    bool multiThread;

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
    void ciecam_02float(CieImage* ncie, float adap, int pW, int pwb, LabImage* lab, const procparams::ProcParams* params,
                        const ColorAppearance & customColCurve1, const ColorAppearance & customColCurve, const ColorAppearance & customColCurve3,
                        LUTu &histLCAM, LUTu &histCCAM, LUTf & CAMBrightCurveJ, LUTf & CAMBrightCurveQ, float &mean, int Iterates, int scale, bool execsharp, float &d, float &dj, float &yb, int rtt,
                        bool showSharpMask = false);
    void chromiLuminanceCurve(PipetteBuffer *pipetteBuffer, int pW, LabImage* lold, LabImage* lnew, LUTf &acurve, LUTf &bcurve, LUTf & satcurve, LUTf & satclcurve, LUTf &clcurve, LUTf &curve, bool utili, bool autili, bool butili, bool ccutili, bool cclutili, bool clcutili, LUTu &histCCurve, LUTu &histLurve);
    void vibrance(LabImage* lab);         //Jacques' vibrance
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
    void EPDToneMapCIE(CieImage *ncie, float a_w, float c_, int Wid, int Hei, float minQ, float maxQ, unsigned int Iterates = 0, int skip = 1);

    // pyramid denoise
//    procparams::DirPyrDenoiseParams dnparams;
    void dirpyr(LabImage* data_fine, LabImage* data_coarse, int level, LUTf &rangefn_L, LUTf &rangefn_ab,
                int pitch, int scale, const int luma, int chroma);
    void idirpyr(LabImage* data_coarse, LabImage* data_fine, int level, LUTf &rangefn_L, LUTf & nrwt_l, LUTf & nrwt_ab,
                 int pitch, int scale, const int luma, const int chroma/*, LUTf & Lcurve, LUTf & abcurve*/);

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
    bool WaveletDenoiseAllL(const wavelet_decomposition &WaveletCoeffs_L, float *noisevarlum, float madL[8][3], float * vari, int edge);
    bool WaveletDenoiseAllAB(const wavelet_decomposition &WaveletCoeffs_L, const wavelet_decomposition &WaveletCoeffs_ab, float *noisevarchrom, float madL[8][3], float noisevar_ab, const bool useNoiseCCurve, bool autoch, bool denoiseMethodRgb);
    void WaveletDenoiseAll_info(int levwav, const wavelet_decomposition &WaveletCoeffs_a,
                                const wavelet_decomposition &WaveletCoeffs_b, float **noisevarlum, float **noisevarchrom, float **noisevarhue, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float & minblueaut, int schoice, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc,
                                float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool denoiseMethodRgb);

    bool WaveletDenoiseAll_BiShrinkL(const wavelet_decomposition &WaveletCoeffs_L, float *noisevarlum, float madL[8][3]);
    bool WaveletDenoiseAll_BiShrinkAB(const wavelet_decomposition &WaveletCoeffs_L, const wavelet_decomposition &WaveletCoeffs_ab, float *noisevarchrom, float madL[8][3], float noisevar_ab,
                                      const bool useNoiseCCurve,  bool autoch, bool denoiseMethodRgb);
    void ShrinkAllL(const wavelet_decomposition &WaveletCoeffs_L, float **buffer, int level, int dir, float *noisevarlum, float * madL, float * vari, int edge);
    void ShrinkAllAB(const wavelet_decomposition &WaveletCoeffs_L, const wavelet_decomposition &WaveletCoeffs_ab, float **buffer, int level, int dir,
                     float *noisevarchrom, float noisevar_ab, const bool useNoiseCCurve, bool autoch, bool denoiseMethodRgb, float * madL, float * madaab = nullptr, bool madCalculated = false);
    void ShrinkAll_info(float ** WavCoeffs_a, float ** WavCoeffs_b,
                        int W_ab, int H_ab, float **noisevarlum, float **noisevarchrom, float **noisevarhue, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut, int schoice, int lvl, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc,
                        float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool denoiseMethodRgb);
    void Noise_residualAB(const wavelet_decomposition &WaveletCoeffs_ab, float &chresid, float &chmaxresid, bool denoiseMethodRgb);
    void calcautodn_info(float &chaut, float &delta, int Nb, int levaut, float maxmax, float lumema, float chromina, int mode, int lissage, float redyel, float skinc, float nsknc);
    float Mad(const float * DataList, int datalen);
    float MadRgb(const float * DataList, int datalen);

    // pyramid wavelet
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

    void dehaze(Imagefloat *rgb);
    void ToneMapFattal02(Imagefloat *rgb);
    void localContrast(LabImage *lab);
    void colorToningLabGrid(LabImage *lab, int xstart, int xend, int ystart, int yend, bool MultiThread);
    void shadowsHighlights(LabImage *lab);
    void softLight(LabImage *lab);
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
