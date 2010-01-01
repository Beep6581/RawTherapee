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
#ifndef _IMPROCFUN_H_
#define _IMPROCFUN_H_

#include <image16.h>
#include <image8.h>
#include <procparams.h>
#include <shmap.h>
#include <coord2d.h>

namespace rtengine {

using namespace procparams;

class LabImage {
    private:
        bool fromImage;

    public:
        int W, H;
        unsigned short** L;
        short** a;
        short** b;
        
     LabImage (int w, int h);
     LabImage (Image16* im);
    ~LabImage ();
};

class ImProcFunctions {

    protected:
        struct STemp {
            int cx, cy, sx, sy, oW, oH;
        };
		cmsHTRANSFORM monitorTransform;

        void transform_         (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to);
        void simpltransform_    (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to);
        void vignetting_        (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to);
        void transform_sep_     (Image16* original, Image16* transformed, const ProcParams* params, STemp sizes, int row_from, int row_to);
        void rgbProc_           (Image16* working, LabImage* lab, const ProcParams* params, int* tonecurve, SHMap* shmap, int row_from, int row_to);
        void lab2rgb_           (LabImage* lab, Image8* image, int row_from, int row_to);
        void colorCurve_        (LabImage* lold, LabImage* lnew, const ProcParams* params, int row_from, int row_to, double* cmultiplier);
        void sharpenHaloCtrl    (LabImage* lab, const ProcParams* params, unsigned short** blurmap, unsigned short** base, int W, int row_from, int row_to);
        void firstAnalysis_     (Image16* original, Glib::ustring wprofile, int* histogram, int* chroma_radius, int row_from, int row_to);
        void resize_            (Image16* src, Image16* dst, ResizeParams params, int row_from, int row_to);
        void damping_           (float** aI, unsigned short** aO, float damping, int W, int rowfrom, int rowto);

    public:

        static int* cacheL;
        static int* cachea;
        static int* cacheb;
        static int* xcache;
        static int* ycache;
        static int* zcache;
    
        int chroma_scale;
        int chroma_radius;

        double lumimul[3];
        static unsigned short gamma2curve[65536];
        

        static void initCache ();
        
        ImProcFunctions () : monitorTransform(NULL) {}
		void release    ();
		
		
        void firstAnalysis  (Image16* working, const ProcParams* params, int* vhist16, double gamma);
    
        void rgbProc        (Image16* working, LabImage* lab, const ProcParams* params, int* tonecurve, SHMap* shmap);
        void luminanceCurve (LabImage* lold, LabImage* lnew, int* curve, int row_from, int row_to);
        void colorCurve     (LabImage* lold, LabImage* lnew, const ProcParams* params);
        void sharpening     (LabImage* lab, const ProcParams* params, double scale, unsigned short** buffer);
        void lumadenoise    (LabImage* lab, const ProcParams* params, double scale, int** buffer);
        void colordenoise   (LabImage* lab, const ProcParams* params, double scale, int** buffer);
        void transform      (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int sx, int sy, int oW, int oH);
        void simpltransform (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int sx, int sy, int oW, int oH);
        void vignetting     (Image16* original, Image16* transformed, const ProcParams* params, int cx, int cy, int oW, int oH);
        void lab2rgb        (LabImage* lab, Image8* image);
        void resize         (Image16* src, Image16* dst, ResizeParams params);

        bool transCoord     (const ProcParams* params, int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv);
        bool transCoord     (const ProcParams* params, int W, int H, std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue);
        void deconvsharpening(LabImage* lab, const ProcParams* params, double scale, unsigned short** buffer);

        Image8*     lab2rgb     (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile);
        Image16*    lab2rgb16   (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile);

        static void getAutoExp  (int* histogram, int histcompr, double expcomp, double clip, double& br, int& bl);
};
};
#endif
