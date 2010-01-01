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
#include <stdimagesource.h>
#include <mytime.h>
#include <iccstore.h>
#define MAXVAL  0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)
#include <curves.h>

#undef THREAD_PRIORITY_NORMAL

namespace rtengine {

extern const Settings* settings;

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

#define HR_SCALE 2
StdImageSource::StdImageSource () : ImageSource(), plistener(NULL), img(NULL) {

    hrmap[0] = NULL;
    hrmap[1] = NULL;
    hrmap[2] = NULL;
	needhr = NULL;
	embProfile = NULL;
    idata = NULL;
 }

StdImageSource::~StdImageSource () {

    delete idata;
    
    if (hrmap[0]!=NULL) {
        int dh = img->height/HR_SCALE;
        freeArray<float>(hrmap[0], dh);
        freeArray<float>(hrmap[1], dh);
        freeArray<float>(hrmap[2], dh);
    }       

    delete img;

	if (needhr)
        freeArray<char>(needhr, img->height);
}

int StdImageSource::load (Glib::ustring fname) {

    fileName = fname;

    img = new Image16 ();
    if (plistener) {
        plistener->setProgressStr ("Loading...");
        plistener->setProgress (0.0);
        img->setProgressListener (plistener);
    }

    int error = img->load (fname);
    if (error) {
        delete img;
        img = NULL;
        return error;
    }

    embProfile = img->getEmbeddedProfile ();
    idata = new ImageData (fname); 

    if (plistener) {
        plistener->setProgressStr ("Ready.");
        plistener->setProgress (1.0);
    }

    wb = ColorTemp (1.0, 1.0, 1.0);

    return 0;
}

void StdImageSource::transform (PreviewProps pp, int tran, int &sx1, int &sy1, int &sx2, int &sy2) {

    int W = img->width;
    int H = img->height;
    int sw = W, sh = H;  
    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = H;
        sh = W;
    }
    int ppx = pp.x, ppy = pp.y;
    if (tran & TR_HFLIP) 
        ppx = sw - pp.x - pp.w;
    if (tran & TR_VFLIP) 
        ppy = sh - pp.y - pp.h;
    
    sx1 = ppx;
    sy1 = ppy;
    sx2 = ppx + pp.w;
    sy2 = ppy + pp.h;
    
    if ((tran & TR_ROT) == TR_R180) {
        sx1 = W - ppx - pp.w;
        sy1 = H - ppy - pp.h;
        sx2 = sx1 + pp.w;
        sy2 = sy1 + pp.h;
    }
    else if ((tran & TR_ROT) == TR_R90) {
        sx1 = ppy;
        sy1 = H - ppx - pp.w;
        sx2 = sx1 + pp.h;
        sy2 = sy1 + pp.w;
    }
    else if ((tran & TR_ROT) == TR_R270) {
        sx1 = W - ppy - pp.h;
        sy1 = ppx;
        sx2 = sx1 + pp.h;
        sy2 = sy1 + pp.w;
    }   
    printf ("ppx %d ppy %d ppw %d pph %d s: %d %d %d %d\n",pp.x, pp.y,pp.w,pp.h,sx1,sy1,sx2,sy2);
}

void StdImageSource::getImage_ (ColorTemp ctemp, int tran, Image16* image, PreviewProps pp, bool first, HRecParams hrp) {

    // compute channel multipliers
    double rm, gm, bm;
    ctemp.getMultipliers (rm, gm, bm);
    rm = 1.0 / rm;
    gm = 1.0 / gm;
    bm = 1.0 / bm;
    double mul_lum = 0.299*rm + 0.587*gm + 0.114*bm;
    rm /= mul_lum;
    gm /= mul_lum;
    bm /= mul_lum;    

    int sx1, sy1, sx2, sy2;
    transform (pp, tran, sx1, sy1, sx2, sy2);
    
    int imwidth  = (sx2 - sx1) / pp.skip + ((sx2 - sx1) % pp.skip > 0);
    int imheight = (sy2 - sy1) / pp.skip + ((sy2 - sy1) % pp.skip > 0);

    int istart = sy1, iend = sy2, ix = 0;
    if (first) {
        iend = istart;
        while (iend<(sy1+sy2)/2) {
            iend+=pp.skip;
        }
    }
    else {
        while (istart<(sy1+sy2)/2) {
            ix++;
            istart+=pp.skip;
        }
    }

    int mtran = tran;
    int skip = pp.skip;

    unsigned short* red  = new unsigned short[imwidth];
    unsigned short* grn  = new unsigned short[imwidth];
    unsigned short* blue = new unsigned short[imwidth];
        
    for (int i=istart; i<iend; i+=skip, ix++) {
    for (int j=0,jx=sx1; j<imwidth; j++,jx+=pp.skip) {
            red[j]  = img->r[i][jx];
            grn[j]  = img->g[i][jx];
            blue[j] = img->b[i][jx];
        }        
//        if (hrp.enabled)
//            hlRecovery (red, grn, blue, i, sx1, sx2, pp.skip);
    
        if ((mtran & TR_ROT) == TR_R180) 
            for (int j=0; j<imwidth; j++) {
                image->r[imheight-1-ix][imwidth-1-j] = (int)CLIP(rm*red[j]);
                image->g[imheight-1-ix][imwidth-1-j] = (int)CLIP(gm*grn[j]);
                image->b[imheight-1-ix][imwidth-1-j] = (int)CLIP(bm*blue[j]);
            }
        else if ((mtran & TR_ROT) == TR_R90) 
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {
                image->r[j][imheight-1-ix] = (int)CLIP(rm*red[j]);
                image->g[j][imheight-1-ix] = (int)CLIP(gm*grn[j]);
                image->b[j][imheight-1-ix] = (int)CLIP(bm*blue[j]);
        }
        else if ((mtran & TR_ROT) == TR_R270) 
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {
                image->r[imwidth-1-j][ix] = (int)CLIP(rm*red[j]);
                image->g[imwidth-1-j][ix] = (int)CLIP(gm*grn[j]);
                image->b[imwidth-1-j][ix] = (int)CLIP(bm*blue[j]);
            }
        else {
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {
                image->r[ix][j] = (int)CLIP(rm*red[j]);
                image->g[ix][j] = (int)CLIP(gm*grn[j]);
                image->b[ix][j] = (int)CLIP(bm*blue[j]);
            }
        }
    }
    
    delete [] red;
    delete [] grn;
    delete [] blue;
}

void StdImageSource::getImage (ColorTemp ctemp, int tran, Image16* image, PreviewProps pp, HRecParams hrp, ColorManagementParams cmp) {

    MyTime t1,t2;

    t1.set ();

//    if (hrp.enabled==true && hrmap[0]==NULL) 
//        updateHLRecoveryMap ();
    if (settings->dualThreadEnabled) {
        Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &StdImageSource::getImage_), ctemp, tran, image, pp, true, hrp), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &StdImageSource::getImage_), ctemp, tran, image, pp, false, hrp), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        thread1->join ();
        thread2->join ();
    }
    else {
        getImage_ (ctemp, tran, image, pp, true, hrp);
        getImage_ (ctemp, tran, image, pp, false, hrp);
    }

    colorSpaceConversion (image, cmp, embProfile);
    
    // Flip if needed
    if (tran & TR_HFLIP)
        hflip (image);
    if (tran & TR_VFLIP)
        vflip (image);

    t2.set ();
}

void StdImageSource::colorSpaceConversion (Image16* im, ColorManagementParams cmp, cmsHPROFILE embedded) {

    cmsHPROFILE in;
    cmsHPROFILE out = iccStore.workingSpace (cmp.working);
    if (cmp.input=="(embedded)" || cmp.input=="" || cmp.input=="(camera)") {
        if (embedded)
            in = embedded;
        else
            in = iccStore.getsRGBProfile ();
    }
    else if (cmp.input!="(none)") {
        in = iccStore.getProfile (cmp.input);
        if (in==NULL && embedded)
            in = embedded;
        else if (in==NULL)
            in = iccStore.getsRGBProfile ();
        else if (cmp.gammaOnInput) 
            for (int i=0; i<im->height; i++)
                for (int j=0; j<im->width; j++) {
                    im->r[i][j] = CurveFactory::gamma (im->r[i][j]);
                    im->g[i][j] = CurveFactory::gamma (im->g[i][j]);
                    im->b[i][j] = CurveFactory::gamma (im->b[i][j]);
                }
    }

    if (cmp.input!="(none)") {
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (in, TYPE_RGB_16_PLANAR, out, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, 0);
        lcmsMutex->unlock ();
        cmsDoTransform (hTransform, im->data, im->data, im->planestride/2);
        cmsDeleteTransform(hTransform);
    }
}

void StdImageSource::getFullSize (int& w, int& h, int tr) {

    w = img->width;
    h = img->height;
    if ((tr & TR_ROT) == TR_R90 || (tr & TR_ROT) == TR_R270) {
        w = img->height;
        h = img->width;
    }
}
 
void StdImageSource::getSize (int tran, PreviewProps pp, int& w, int& h) {

    w = pp.w / pp.skip + (pp.w % pp.skip > 0);
    h = pp.h / pp.skip + (pp.h % pp.skip > 0);
}

void StdImageSource::hflip (Image16* image) {
    int width  = image->width;
    int height = image->height;

    unsigned short* rowr = new unsigned short[width];
    unsigned short* rowg = new unsigned short[width];
    unsigned short* rowb = new unsigned short[width];
    for (int i=0; i<height; i++) {
      for (int j=0; j<width; j++) {
        rowr[j] = image->r[i][width-1-j];
        rowg[j] = image->g[i][width-1-j];
        rowb[j] = image->b[i][width-1-j];
      }
      memcpy (image->r[i], rowr, width*sizeof(unsigned short));
      memcpy (image->g[i], rowg, width*sizeof(unsigned short));
      memcpy (image->b[i], rowb, width*sizeof(unsigned short));
    }
    delete [] rowr;
    delete [] rowg;
    delete [] rowb;
}

void StdImageSource::vflip (Image16* image) {
    int width  = image->width;
    int height = image->height;

    register unsigned short tmp;
    for (int i=0; i<height/2; i++) 
      for (int j=0; j<width; j++) {
        tmp = image->r[i][j]; 
        image->r[i][j] = image->r[height-1-i][j];
        image->r[height-1-i][j] = tmp;
        tmp = image->g[i][j]; 
        image->g[i][j] = image->g[height-1-i][j];
        image->g[height-1-i][j] = tmp;
        tmp = image->b[i][j]; 
        image->b[i][j] = image->b[height-1-i][j];
        image->b[height-1-i][j] = tmp;
      }
}
/*    
void hlRecovery (unsigned short* red, unsigned short* green, unsigned short* blue, int H, int W, int i, int sx1, int sx2, int skip, char** needhr, float** hrmap[3]);
void hlmultipliers (int** rec[3], int max[3], int dh, int dw);

void StdImageSource::updateHLRecoveryMap () {

    // detect maximal pixel values
    int maxr = 0, maxg = 0, maxb = 0;
    for (int i=32; i<img->height-32; i++)
        for (int j=32; j<img->width-32; j++) {
            if (img->r[i][j] > maxr) maxr = img->r[i][j];
            if (img->g[i][j] > maxg) maxg = img->g[i][j];
            if (img->b[i][j] > maxb) maxb = img->b[i][j];
        }

    maxr = maxr * 19 / 20;
    maxg = maxg * 19 / 20;
    maxb = maxb * 19 / 20;
    max[0] = maxr;
    max[1] = maxg;
    max[2] = maxb;
    
    // downscale image
    int dw = img->width/HR_SCALE;
    int dh = img->height/HR_SCALE;
    Image16* ds = new Image16 (dw, dh);

    // overburnt areas
    int** rec[3];
    for (int i=0; i<3; i++)
        rec[i] = allocArray<int> (dw, dh);

	if (needhr)
		freeArray<char>(needhr, img->height);
	needhr = allocArray<char> (img->width, img->height);

	for (int i=0; i<img->height; i++)
		for (int j=0; j<img->width; j++)
			if (img->r[i][j]>=max[0] || img->g[i][j]>=max[1] || img->b[i][j]>=max[2])
				needhr[i][j] = 1;
			else
				needhr[i][j] = 0;


    for (int i=0; i<ds->height; i++)
        for (int j=0; j<ds->width; j++) {
            int sumr = 0; int cr = 0;
            int sumg = 0; int cg = 0;
            int sumb = 0; int cb = 0;
            for (int x=0; x<HR_SCALE; x++)
                for (int y=0; y<HR_SCALE; y++) {
                    int ix = HR_SCALE*i+x;
                    int jy = HR_SCALE*j+y;
                    sumr += img->r[ix][jy];
                    if (img->r[ix][jy] < maxr) cr++;
                    sumg += img->g[ix][jy];
                    if (img->g[ix][jy] < maxg) cg++;
                    sumb += img->b[ix][jy];
                    if (img->b[ix][jy] < maxb) cb++;
                }
            if (cr<HR_SCALE*HR_SCALE) rec[0][i][j] = INT_MAX; else rec[0][i][j] = sumr / HR_SCALE/HR_SCALE;
            if (cg<HR_SCALE*HR_SCALE) rec[1][i][j] = INT_MAX; else rec[1][i][j] = sumg / HR_SCALE/HR_SCALE;
            if (cb<HR_SCALE*HR_SCALE) rec[2][i][j] = INT_MAX; else rec[2][i][j] = sumb / HR_SCALE/HR_SCALE;
            ds->r[i][j] = sumr / HR_SCALE/HR_SCALE;
            ds->g[i][j] = sumg / HR_SCALE/HR_SCALE;
            ds->b[i][j] = sumb / HR_SCALE/HR_SCALE;
            
        }

    hlmultipliers (rec, max, dh, dw);

    if (hrmap[0]!=NULL) {
        freeArray<float> (hrmap[0], dh);
        freeArray<float> (hrmap[1], dh);
        freeArray<float> (hrmap[2], dh);
    }

    hrmap[0] = allocArray<float> (dw, dh);
    hrmap[1] = allocArray<float> (dw, dh);
    hrmap[2] = allocArray<float> (dw, dh);

    for (int i=0; i<dh; i++)
        for (int j=0; j<dw; j++) {
            hrmap[0][i][j] = ds->r[i][j]>0 ? (double)rec[0][i][j] / ds->r[i][j] : 1.0;
            hrmap[1][i][j] = ds->g[i][j]>0 ? (double)rec[1][i][j] / ds->g[i][j] : 1.0;
            hrmap[2][i][j] = ds->b[i][j]>0 ? (double)rec[2][i][j] / ds->b[i][j] : 1.0;
        }
    
    delete ds;
    
    freeArray<int> (rec[0], dh);
    freeArray<int> (rec[1], dh);
    freeArray<int> (rec[2], dh);
}

void StdImageSource::hlRecovery (unsigned short* red, unsigned short* green, unsigned short* blue, int i, int sx1, int sx2, int skip) {

    rtengine::hlRecovery (red, green, blue, img->height, img->width, i, sx1, sx2, skip, needhr, hrmap);
}
*/
int StdImageSource::getAEHistogram (int* histogram, int& histcompr) {

    histcompr = 3;

    memset (histogram, 0, (65536>>histcompr)*sizeof(int));

    for (int i=0; i<img->height; i++)
        for (int j=0; j<img->width; j++) {
            histogram[CurveFactory::igamma_srgb (img->r[i][j])>>histcompr]++;
            histogram[CurveFactory::igamma_srgb (img->g[i][j])>>histcompr]++;
            histogram[CurveFactory::igamma_srgb (img->b[i][j])>>histcompr]++;
        }
    return 1;
}

ColorTemp StdImageSource::getAutoWB () {

    double avg_r = 0;
    double avg_g = 0;
    double avg_b = 0;
    int n = 0;
    int p = 6;

    for (int i=1; i<img->height-1; i++)
        for (int j=1; j<img->width-1; j++) {
            if (img->r[i][j]>64000 || img->g[i][j]>64000 || img->b[i][j]>64000)
                continue;
            avg_r += intpow((double)img->r[i][j], p);
            avg_g += intpow((double)img->g[i][j], p);
            avg_b += intpow((double)img->b[i][j], p);
            n++;
        }
    return ColorTemp (pow(avg_r/n, 1.0/p), pow(avg_g/n, 1.0/p), pow(avg_b/n, 1.0/p));
}

void StdImageSource::transformPixel (int x, int y, int tran, int& tx, int& ty) {
    
    int W = img->width;
    int H = img->height;
    int sw = W, sh = H;  
    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = H;
        sh = W;
    }

    int ppx = x, ppy = y;
    if (tran & TR_HFLIP) 
        ppx = sw - 1 - x ;
    if (tran & TR_VFLIP) 
        ppy = sh - 1 - y;
    
    tx = ppx;
    ty = ppy;
    
    if ((tran & TR_ROT) == TR_R180) {
        tx = W - 1 - ppx;
        ty = H - 1 - ppy;
    }
    else if ((tran & TR_ROT) == TR_R90) {
        tx = ppy;
        ty = H - 1 - ppx;
    }
    else if ((tran & TR_ROT) == TR_R270) {
        tx = W - 1 - ppy;
        ty = ppx;
    }   
}

ColorTemp StdImageSource::getSpotWB (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue, int tran) {

    int x; int y;
    double reds = 0, greens = 0, blues = 0;
    int rn = 0, gn = 0, bn = 0;
    for (int i=0; i<red.size(); i++) {
        transformPixel (red[i].x, red[i].y, tran, x, y);
        if (x>=0 && y>=0 && x<img->width && y<img->height) {
            reds += img->r[y][x];
//            img->r[y][x]=0;   // debug!!!
            rn++;
        }
        transformPixel (green[i].x, green[i].y, tran, x, y);
        if (x>=0 && y>=0 && x<img->width && y<img->height) {
            greens += img->g[y][x];
//            img->g[y][x]=0; // debug!!!
            gn++;
        }
        transformPixel (blue[i].x, blue[i].y, tran, x, y);
        if (x>=0 && y>=0 && x<img->width && y<img->height) {
            blues += img->b[y][x];
//            img->b[y][x]=0; // debug!!!
            bn++;
        }
    }
    double img_r, img_g, img_b;
    wb.getMultipliers (img_r, img_g, img_b);
    printf ("AVG: %g %g %g\n", reds/rn, greens/gn, blues/bn);

    return ColorTemp (reds/rn*img_r, greens/gn*img_g, blues/bn*img_b);
}
}

