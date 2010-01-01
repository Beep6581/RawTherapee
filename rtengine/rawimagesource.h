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

#include <imagesource.h>
#include <lcms.h>
#define HR_SCALE 2

namespace rtengine {

template<class T> void freeArray (T** a, int H) {
  for (int i=0; i<H; i++)
    delete [] a[i];
  delete [] a;
}

template<class T> T** allocArray (int W, int H) {

    T** t = new T*[H];
    for (int i=0; i<H; i++)
        t[i] = new T[W];
    return t;
}


template<class T> void freeArray2 (T** a, int H) {
  for (int i=0; i<H; i++)
    delete [] a[i];
}

class RawImageSource : public ImageSource {

    protected:
        Glib::Mutex isrcMutex;
    
        int W, H;
        ColorTemp wb;
        ProgressListener* plistener;
        double coeff[3][3];
        double icoeff[3][3];
        double cam[3][3];
        double icam[3][3];
        bool fuji;
        bool d1x;
        int border;
        char** hpmap;
        float** hrmap[3];   // for color propagation
        char** needhr;      // for color propagation
        int max[3];
        double defGain;
        bool full;
		Glib::ustring oldmethod;
		cmsHPROFILE camProfile;
		cmsHPROFILE embProfile;

        RawImage* ri;
        
        // to accelerate CIELAB conversion:
        double lc00, lc01, lc02, lc10, lc11, lc12, lc20, lc21, lc22;
        double* cache;
        int threshold;
       
        // the interpolated green plane:
        unsigned short** green; 
    
        void hphd_vertical       (float** hpmap, int col_from, int col_to);
        void hphd_horizontal     (float** hpmap, int row_from, int row_to);
        void hphd_green          (int row_from, int row_to);
        void correction_YIQ_LQ_  (Image16* im, int row_from, int row_to);
        void hlRecovery          (std::string method, unsigned short* red, unsigned short* green, unsigned short* blue, int i, int sx1, int width, int skip);
        int  defTransform        (int tran);
        void rotateLine          (unsigned short* line, unsigned short** channel, int tran, int i, int w, int h);
        void transformRect       (PreviewProps pp, int tran, int &sx1, int &sy1, int &width, int &height, int &fw);
        void transformPosition   (int x, int y, int tran, int& tx, int& ty);

        void updateHLRecoveryMap (std::string method, double rm, double gm, double bm);
        void updateHLRecoveryMap_ColorPropagation ();
        void HLRecovery_ColorPropagation (unsigned short* red, unsigned short* green, unsigned short* blue, int i, int sx1, int width, int skip);
        
    public:
        RawImageSource ();
        ~RawImageSource ();
    
        int         load        (Glib::ustring fname);
        void        getImage    (ColorTemp ctemp, int tran, Image16* image, PreviewProps pp, HRecParams hrp, ColorManagementParams cmp);
        ColorTemp   getWB       () { return wb; }
        ColorTemp   getAutoWB   ();
        ColorTemp   getSpotWB   (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue, int tran);

        double      getDefGain  () { return defGain; }
        double      getGamma    () { return 2.2; }
        
        void        getFullSize (int& w, int& h, int tr = TR_NONE);
        void        getSize     (int tran, PreviewProps pp, int& w, int& h);

        ImageData*  getImageData () { return idata; }
        void        setProgressListener (ProgressListener* pl) { plistener = pl; }
        int         getAEHistogram (int* histogram, int& histcompr);

        static void colorSpaceConversion (Image16* im, ColorManagementParams cmp, cmsHPROFILE embedded, cmsHPROFILE camprofile, double cam[3][3], double& defgain);
        static void inverse33 (double (*coeff)[3], double (*icoeff)[3]);

        static void HLRecovery_Luminance (unsigned short* rin, unsigned short* gin, unsigned short* bin, unsigned short* rout, unsigned short* gout, unsigned short* bout, int width, int maxval);
        static void HLRecovery_CIELab (unsigned short* rin, unsigned short* gin, unsigned short* bin, unsigned short* rout, unsigned short* gout, unsigned short* bout, int width, int maxval, double cam[3][3], double icam[3][3]);

    protected:
                void correction_YIQ_LQ  (Image16* i, int times);
        inline  void convert_row_to_YIQ (unsigned short* r, unsigned short* g, unsigned short* b, int* Y, int* I, int* Q, int W);
        inline  void convert_row_to_RGB (unsigned short* r, unsigned short* g, unsigned short* b, int* Y, int* I, int* Q, int W);

        inline  void convert_to_cielab_row  (unsigned short* ar, unsigned short* ag, unsigned short* ab, short* oL, short* oa, short* ob);
        inline  void interpolate_row_g      (unsigned short* agh, unsigned short* agv, int i);
        inline  void interpolate_row_rb     (unsigned short* ar, unsigned short* ab, unsigned short* pg, unsigned short* cg, unsigned short* ng, int i);
        inline  void interpolate_row_rb_mul_pp (unsigned short* ar, unsigned short* ab, unsigned short* pg, unsigned short* cg, unsigned short* ng, int i, double r_mul, double g_mul, double b_mul, int x1, int width, int skip);

        void    eahd_demosaic   ();
        void    hphd_demosaic   ();
        void    vng4_demosaic   ();

        void    transLine   (unsigned short* red, unsigned short* green, unsigned short* blue, int i, Image16* image, int tran, int imw, int imh, int fw);
        void    hflip       (Image16* im);
        void    vflip       (Image16* im);      
        
};
};
#endif
