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
#ifndef _RAWIMAGESOURCE_
#define _RAWIMAGESOURCE_

#include "imagesource.h"
#include <lcms2.h>
#include "dcp.h"
#include "array2D.h"
#include "curves.h"
#include <fftw3.h>
#include "color.h"
#include "iimage.h"
#include "../rtgui/cacheimagedata.h"
#include "../rtgui/threadutils.h"

#define HR_SCALE 2

namespace rtengine {

// these two functions "simulate" and jagged array, but just use two allocs
template<class T> T** allocArray (int W, int H, bool initZero=false) {

    T** t = new T*[H];
    t[0] = new T[H*W];

    if (initZero) memset(t[0],0,sizeof(T)*W*H);

    for (int i=1; i<H; i++)
        t[i] = t[i-1]+W;

    return t;
}

template<class T> void freeArray (T** a, int H) {

    delete [] a[0];
    delete [] a;
}


template<class T> void freeArray2 (T** a, int H) {
  //for (int i=0; i<H; i++)
    delete [] a[0];
}

class RawImageSource : public ImageSource {

    private:
        static DiagonalCurve *phaseOneIccCurve;
        static DiagonalCurve *phaseOneIccCurveInv;
        static LUTf invGrad;  // for fast_demosaic
        static LUTf initInvGrad ();
        static bool findInputProfile(Glib::ustring inProfile, cmsHPROFILE embedded, std::string camName, DCPProfile **dcpProf, cmsHPROFILE& in);
        static void colorSpaceConversion_ (Imagefloat* im, ColorManagementParams &cmp, ColorTemp &wb, double pre_mul[3], const RAWParams &raw, cmsHPROFILE embedded, cmsHPROFILE camprofile, double cam[3][3], const std::string &camName);

    protected:
        MyMutex getImageMutex;  // locks getImage

        int W, H;
        ColorTemp camera_wb;
        ProgressListener* plistener;
        float scale_mul[4]; // multiplier for each color
        float c_black[4]; // copy of cblack Dcraw for black level
        float c_white[4];
        float cblacksom[4];
        float ref_pre_mul[4];
        double refwb_red;
        double refwb_green;
        double refwb_blue;
        double rgb_cam[3][3];
        double cam_rgb[3][3];
        double xyz_cam[3][3];
        double cam_xyz[3][3];
        bool fuji;
        bool d1x;
        int border;
        //char** hpmap;
        float** hrmap[3];   // for color propagation
        char** needhr;      // for color propagation
        int max_3[3];
        float chmax[4],hlmax[4];
        double initialGain; // initial gain calculated after scale_colors
        double camInitialGain;
        double defGain;
        bool full;
        cmsHPROFILE camProfile;
        cmsHPROFILE embProfile;
        bool rgbSourceModified;

        RawImage* ri;  // Copy of raw pixels, NOT corrected for initial gain, blackpoint etc.

        // to accelerate CIELAB conversion:
        double lc00, lc01, lc02, lc10, lc11, lc12, lc20, lc21, lc22;
        double* cache;
        int threshold;

        array2D<float> rawData;  // holds preprocessed pixel values, rowData[i][j] corresponds to the ith row and jth column

        // the interpolated green plane:
        array2D<float> green;
        // the interpolated red plane:
        array2D<float> red;
        // the interpolated blue plane:
        array2D<float> blue;


        void hphd_vertical       (float** hpmap, int col_from, int col_to);
        void hphd_horizontal     (float** hpmap, int row_from, int row_to);
        void hphd_green          (float** hpmap);
        void processFalseColorCorrectionThread (Imagefloat* im, int row_from, int row_to);
        void hlRecovery          (std::string method, float* red, float* green, float* blue, int i, int sx1, int width, int skip, const RAWParams &raw, float* hlmax);
        int  defTransform        (int tran);
        void rotateLine          (float* line, PlanarPtr<float> &channel, int tran, int i, int w, int h);
        void transformRect       (PreviewProps pp, int tran, int &sx1, int &sy1, int &width, int &height, int &fw);
        void transformPosition   (int x, int y, int tran, int& tx, int& ty);

        void updateHLRecoveryMap_ColorPropagation ();
        void HLRecovery_ColorPropagation (float* red, float* green, float* blue, int i, int sx1, int width, int skip);
        unsigned FC(int row, int col){ return ri->FC(row,col); }
        inline void getRowStartEnd (int x, int &start, int &end);
        static void getProfilePreprocParams(cmsHPROFILE in, float& gammafac, float& lineFac, float& lineSum);


    public:
        RawImageSource ();
        ~RawImageSource ();

        int         load        (Glib::ustring fname, bool batch = false);
        void        preprocess  (const RAWParams &raw, const LensProfParams &lensProf, const CoarseTransformParams& coarse);
        void        demosaic    (const RAWParams &raw);
        void        flushRawData      ();
        void        flushRGB          ();
        void        HLRecovery_Global (ToneCurveParams hrp);
        void        refinement_lassus (int PassCount);
        void        refinement(int PassCount);

        bool        IsrgbSourceModified() {return rgbSourceModified;} // tracks whether cached rgb output of demosaic has been modified

		void		processFlatField(const RAWParams &raw, RawImage *riFlatFile, unsigned short black[4]);
        void        copyOriginalPixels(const RAWParams &raw, RawImage *ri, RawImage *riDark, RawImage *riFlatFile  );
        void        cfaboxblur  (RawImage *riFlatFile, float* cfablur, int boxH, int boxW );
        void        scaleColors (int winx,int winy,int winw,int winh, const RAWParams &raw);// raw for cblack

        void        getImage    (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, ToneCurveParams hrp, ColorManagementParams cmp, RAWParams raw);
        eSensorType getSensorType () { return ri!=NULL ? ri->getSensorType() : ST_NONE; }
        ColorTemp   getWB       () { return camera_wb; }
        void        getAutoWBMultipliers (double &rm, double &gm, double &bm);
        ColorTemp   getSpotWB   (std::vector<Coord2D> &red, std::vector<Coord2D> &green, std::vector<Coord2D> &blue, int tran, double equal);
        bool        isWBProviderReady () { return rawData; }

        double      getDefGain  () { return defGain; }

        double      getGamma    () { return Color::sRGBGamma; }

        void        getFullSize (int& w, int& h, int tr = TR_NONE);
        void        getSize     (int tran, PreviewProps pp, int& w, int& h);
        int         getRotateDegree() const { return ri->get_rotateDegree(); }

        ImageData*  getImageData () { return idata; }
        ImageMatrices* getImageMatrices () { return &imatrices; }
        bool        isRAW() const { return true; }

        void        setProgressListener (ProgressListener* pl) { plistener = pl; }
        void        getAutoExpHistogram (LUTu & histogram, int& histcompr);
        void        getRAWHistogram (LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw);

        void convertColorSpace(Imagefloat* image, ColorManagementParams cmp, ColorTemp &wb, RAWParams raw);
        static void colorSpaceConversion   (Imagefloat* im, ColorManagementParams cmp, ColorTemp &wb, double pre_mul[3], RAWParams raw, cmsHPROFILE embedded, cmsHPROFILE camprofile, double cam[3][3], std::string camName) {
            colorSpaceConversion_ (im, cmp, wb, pre_mul, raw, embedded, camprofile, cam, camName);
        }
        static void inverse33 (const double (*coeff)[3], double (*icoeff)[3]);

        void boxblur2(float** src, float** dst, int H, int W, int box );
        void boxblur_resamp(float **src, float **dst, int H, int W, int box, int samp );

        //void boxblur_resamp(float **red, float **green, float **blue, int H, int W, float thresh[3], float max[3],
        //                multi_array2D<float,3> & hfsize, multi_array2D<float,3> & hilite, int box );
        void HLRecovery_inpaint (float** red, float** green, float** blue);
        //void HLRecovery_inpaint ();

        static void HLRecovery_Luminance (float* rin, float* gin, float* bin, float* rout, float* gout, float* bout, int width, float maxval);
        static void HLRecovery_CIELab (float* rin, float* gin, float* bin, float* rout, float* gout, float* bout, int width, float maxval, double cam[3][3], double icam[3][3]);
        static void HLRecovery_blend (float* rin, float* gin, float* bin, int width, float maxval, float* pre_mul, const RAWParams &raw, float* hlmax);
	static void init ();
	static void cleanup ();

    protected:
        typedef unsigned short ushort;
                void processFalseColorCorrection (Imagefloat* i, int steps);
        inline  void convert_row_to_YIQ (float* r, float* g, float* b, float* Y, float* I, float* Q, int W);
        inline  void convert_row_to_RGB (float* r, float* g, float* b, float* Y, float* I, float* Q, int W);

        inline  void convert_to_cielab_row  (float* ar, float* ag, float* ab, float* oL, float* oa, float* ob);
        inline  void interpolate_row_g      (float* agh, float* agv, int i);
        inline  void interpolate_row_rb     (float* ar, float* ab, float* pg, float* cg, float* ng, int i);
        inline  void interpolate_row_rb_mul_pp (float* ar, float* ab, float* pg, float* cg, float* ng, int i, float r_mul, float g_mul, float b_mul, int x1, int width, int skip);

        int  LinEqSolve( int nDim, double* pfMatr, double* pfVect, double* pfSolution);//Emil's CA auto correction
        void CA_correct_RT	(double cared, double cablue);
        void ddct8x8s(int isgn, float a[8][8]);
        void processRawWhitepoint (float expos, float preser);  // exposure before interpolation

        int  interpolateBadPixelsBayer( PixelsMap &bitmapBads );
        int  interpolateBadPixelsNColours( PixelsMap &bitmapBads, const int colours );
        int  interpolateBadPixelsXtrans( PixelsMap &bitmapBads );
        int  findHotDeadPixels( PixelsMap &bpMap, float thresh, bool findHotPixels, bool findDeadPixels );

        void cfa_linedn (float linenoiselevel);//Emil's line denoise
        void cfa_tile_denoise (fftwf_complex * fcfablox, int vblk, int hblk, int numblox_H, int numblox_W, float noisevar, float * rolloff );
        void cfa_output_tile_row (float *cfabloxrow, float ** cfahipassdn, float ** tilemask_out, int height, int width, int top, int blkrad);

        void green_equilibrate (float greenthresh);//Emil's green equilibration

        void nodemosaic(bool bw);
        void eahd_demosaic();
        void hphd_demosaic();
        void vng4_demosaic();
        void ppg_demosaic();
        void jdl_interpolate_omp();
        void igv_interpolate(int winw, int winh);
        void lmmse_interpolate_omp(int winw, int winh, int iterations);

        void amaze_demosaic_RT(int winx, int winy, int winw, int winh);//Emil's code for AMaZE
        void fast_demosaic(int winx, int winy, int winw, int winh );//Emil's code for fast demosaicing
        void dcb_demosaic(int iterations, bool dcb_enhance);
        void ahd_demosaic(int winx, int winy, int winw, int winh);
        void border_interpolate(unsigned int border, float (*image)[4], unsigned int start = 0, unsigned int end = 0);
        void border_interpolate2(int winw, int winh, int lborders);
        void dcb_initTileLimits(int &colMin, int &rowMin, int &colMax, int &rowMax, int x0, int y0, int border);
        void fill_raw( float (*cache )[4], int x0, int y0, float** rawData);
        void fill_border( float (*cache )[4], int border, int x0, int y0);
        void copy_to_buffer(float (*image2)[3], float (*image)[4]);
        void dcb_hid(float (*image)[4], float (*bufferH)[3], float (*bufferV)[3], int x0, int y0);
        void dcb_color(float (*image)[4], int x0, int y0);
        void dcb_hid2(float (*image)[4], int x0, int y0);
        void dcb_map(float (*image)[4], int x0, int y0);
        void dcb_correction(float (*image)[4], int x0, int y0);
        void dcb_pp(float (*image)[4], int x0, int y0);
        void dcb_correction2(float (*image)[4], int x0, int y0);
        void restore_from_buffer(float (*image)[4], float (*image2)[3]);
        void dcb_refinement(float (*image)[4], int x0, int y0);
        void dcb_color_full(float (*image)[4], int x0, int y0, float (*chroma)[2]);
        void cielab (const float (*rgb)[3], float* l, float* a, float *b, const int width, const int height, const int labWidth, const float xyz_cam[3][3]);
        void xtransborder_interpolate (int border);
        void xtrans_interpolate (int passes, bool useCieLab);
        void fast_xtrans_interpolate ();
        void    transLine   (float* red, float* green, float* blue, int i, Imagefloat* image, int tran, int imw, int imh, int fw);
        void    hflip       (Imagefloat* im);
        void    vflip       (Imagefloat* im);

};
}
#endif
