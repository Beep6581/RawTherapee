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
#include "improcfun.h"
#include <glibmm.h>
#include "iccstore.h"
#include "iccmatrices.h"
#include "mytime.h"
#include "../rtgui/icmpanel.h"
#include "../rtgui/options.h"
#include "settings.h"
#include "curves.h"


#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {

#undef CLIP
#undef CLIPTO
#undef CMAXVAL

#define CMAXVAL 0xffff
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIP01(a) ((a)>0?((a)<1?(a):1):0)
	
extern const Settings* settings;
	
const double (*wprof[])[3]  = {xyz_sRGB, xyz_adobe, xyz_prophoto, xyz_widegamut, xyz_bruce, xyz_beta, xyz_best};
const double (*iwprof[])[3] = {sRGB_xyz, adobe_xyz, prophoto_xyz, widegamut_xyz, bruce_xyz, beta_xyz, best_xyz};
const char* wprofnames[] = {"sRGB", "Adobe RGB", "ProPhoto", "WideGamut", "BruceRGB", "Beta RGB", "BestRGB"};
const int numprof = 7;

void ImProcFunctions::lab2rgb (LabImage* lab, Image8* image) {
	//MyTime tBeg,tEnd;
 //   tBeg.set();
	//gamutmap(lab);

	if (monitorTransform) {
        
        // cmsDoTransform is relatively expensive
        #pragma omp parallel for
		for (int i=0; i<lab->H; i++) {
            float buffer[3*lab->W];
            //float g;  // unused

            const int ix = i * 3 * lab->W;
            int iy = 0;

			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];

            float fy,fx,fz,x_,y_,z_;

			for (int j=0; j<lab->W; j++) {
								
				fy = (0.00862069 * rL[j]) / 327.68 + 0.137932; // (L+16)/116
				fx = (0.002 * ra[j]) / 327.68 + fy;
				fz = fy - (0.005 * rb[j]) / 327.68;
				
				x_ = f2xyz(fx)*D50x;
				y_ = f2xyz(fy);
				z_ = f2xyz(fz)*D50z;

                buffer[iy++] = CLIP01(x_);
                buffer[iy++] = CLIP01(y_);
                buffer[iy++] = CLIP01(z_);
			}

            cmsDoTransform (monitorTransform, buffer, image->data + ix, lab->W);
		}
        
	} else {

		#pragma omp parallel for if (multiThread)
		for (int i=0; i<lab->H; i++) {
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			int ix = i * 3 * lab->W;
            //float g;  // unused

			float R,G,B;
            float fy,fx,fz,x_,y_,z_;

			for (int j=0; j<lab->W; j++) {
			
				//float L1=rL[j],a1=ra[j],b1=rb[j];//for testing
				
				fy = (0.00862069 * rL[j]) / 327.68 + 0.137932; // (L+16)/116
				fx = (0.002 * ra[j]) / 327.68 + fy;
				fz = fy - (0.005 * rb[j]) / 327.68;
				
				x_ = 65535.0 * f2xyz(fx)*D50x;
				y_ = 65535.0 * f2xyz(fy);
				z_ = 65535.0 * f2xyz(fz)*D50z;

				xyz2srgb(x_,y_,z_,R,G,B);
				
				/* copy RGB */
				//int R1=((int)gamma2curve[(R)]) 
				image->data[ix++] = ((int)gamma2curve[CLIP(R)]) >> 8;
				image->data[ix++] = ((int)gamma2curve[CLIP(G)]) >> 8;
				image->data[ix++] = ((int)gamma2curve[CLIP(B)]) >> 8;
			}
		}
	}

    //tEnd.set();
    //printf("lab2rgb %i %d\n", lab->W, tEnd.etime(tBeg));
}

Image8* ImProcFunctions::lab2rgb (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile) {

	//gamutmap(lab);
	
    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw>lab->W) cw = lab->W-cx;
    if (cy+ch>lab->H) ch = lab->H-cy;

    Image8* image = new Image8 (cw, ch);

    cmsHPROFILE oprof = iccStore->getProfile (profile);

    if (oprof) {
        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16, oprof, TYPE_RGB_8, settings->colorimetricIntent,
            cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety
        lcmsMutex->unlock ();

        // cmsDoTransform is relatively expensive
		#pragma omp parallel for
        for (int i=cy; i<cy+ch; i++) {
            short buffer [3*cw];
            //float g;  // unused

            const int ix = i * 3 * cw;
            int iy = 0;

            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];

            for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

                buffer[iy++] = CLIP((int)x_);
                buffer[iy++] = CLIP((int)y_);
                buffer[iy++] = CLIP((int)z_);
            }

            cmsDoTransform (hTransform, buffer, image->data + ix, cw);
        }

        cmsDeleteTransform(hTransform);
    } else {
		
		float rgb_xyz[3][3];
		
		for (int i=0; i<numprof; i++) {
			if (profile==wprofnames[i]) {
				for (int m=0; m<3; m++) 
					for (int n=0; n<3; n++) {
						rgb_xyz[m][n] = iwprof[i][m][n];
					}
				break;
			}
		}
		
		#pragma omp parallel for if (multiThread)
        for (int i=cy; i<cy+ch; i++) {
			//float g;  // unused
			float R,G,B;
            float* rL = lab->L[i];
            float* ra = lab->a[i];
            float* rb = lab->b[i];
			int ix = 3*i*cw;
            for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

				xyz2rgb(x_,y_,z_,R,G,B,rgb_xyz);

                image->data[ix++] = (int)gamma2curve[CLIP(R)] >> 8;
                image->data[ix++] = (int)gamma2curve[CLIP(G)] >> 8;
                image->data[ix++] = (int)gamma2curve[CLIP(B)] >> 8;
            }
        }
    }
    return image;
}
// for default (not gamma)
Image16* ImProcFunctions::lab2rgb16 (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile) {
	
	//gamutmap(lab);

    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw>lab->W) cw = lab->W-cx;
    if (cy+ch>lab->H) ch = lab->H-cy;

    Image16* image = new Image16 (cw, ch);
   cmsHPROFILE oprof = iccStore->getProfile (profile);
 
 

    if (oprof) {
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
			//float g;  // unused
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			short* xa = (short*)image->r[i-cy];
			short* ya = (short*)image->g[i-cy];
			short* za = (short*)image->b[i-cy];
			for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

				xa[j-cx] = CLIP((int)x_);
				ya[j-cx] = CLIP((int)y_);
				za[j-cx] = CLIP((int)z_);
			}
		}

        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
		cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16_PLANAR, oprof, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, cmsFLAGS_NOOPTIMIZE);
        lcmsMutex->unlock ();

		cmsDoTransform (hTransform, image->data, image->data, image->planestride);
		cmsDeleteTransform(hTransform);
	} else {
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
			//float g; //unused
			float R,G,B;
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

				xyz2srgb(x_,y_,z_,R,G,B);

				image->r[i-cy][j-cx] = (int)gamma2curve[CLIP(R)];
				image->g[i-cy][j-cx] = (int)gamma2curve[CLIP(G)];
				image->b[i-cy][j-cx] = (int)gamma2curve[CLIP(B)];
			}
		}
	}
    return image;
}


// for gamma options (BT709...sRGB linear...)
Image16* ImProcFunctions::lab2rgb16b (LabImage* lab, int cx, int cy, int cw, int ch, Glib::ustring profile, Glib::ustring profi, Glib::ustring gam,  bool freegamma, double gampos, double slpos, double &ga0, double &ga1, double &ga2, double &ga3, double &ga4, double &ga5, double &ga6) {
	
	//gamutmap(lab);

    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw>lab->W) cw = lab->W-cx;
    if (cy+ch>lab->H) ch = lab->H-cy;

    Image16* image = new Image16 (cw, ch);
	//cmsBool  rc = TRUE;  // unused
	float p1,p2,p3,p4,p5,p6;//primaries
	//double ga0,ga1,ga2,ga3,ga4,ga5=0.0,ga6=0.0;//gamma parameters
	double g_a0,g_a1,g_a2,g_a3,g_a4,g_a5;//gamma parameters
	double pwr;
	double ts;
	ga6=0.0;
	pwr=1.0/gampos;
	ts=slpos;
	int mode=0, imax=0;
	
	int t50;
	int select_temp =1;//5003K
	double eps=0.000000001;// not divide by zero
	//primaries for 7 working profiles ==> output profiles
	// eventually to adapt primaries  if RT used special profiles !
	if(profi=="ProPhoto") 	  {p1=0.7347; p2=0.2653; p3=0.1596; p4=0.8404; p5=0.0366; p6=0.0001;select_temp=1;}//Prophoto primaries
	else if (profi=="WideGamut") {p1=0.7350; p2=0.2650; p3=0.1150; p4=0.8260; p5=0.1570; p6=0.0180;select_temp=1;}//Widegamut primaries
	else if (profi=="Adobe RGB") {p1=0.6400; p2=0.3300; p3=0.2100; p4=0.7100; p5=0.1500; p6=0.0600;select_temp=2;}//Adobe primaries
	else if (profi=="sRGB") {p1=0.6400; p2=0.3300; p3=0.3000; p4=0.6000; p5=0.1500; p6=0.0600;select_temp=2;} // sRGB primaries
	else if (profi=="BruceRGB") {p1=0.6400; p2=0.3300; p3=0.2800; p4=0.6500; p5=0.1500; p6=0.0600;select_temp=2;} // Bruce primaries
	else if (profi=="Beta RGB") {p1=0.6888; p2=0.3112; p3=0.1986; p4=0.7551; p5=0.1265; p6=0.0352;select_temp=1;} // Beta primaries
	else if (profi=="BestRGB") {p1=0.7347; p2=0.2653; p3=0.2150; p4=0.7750; p5=0.1300; p6=0.0350;select_temp=1;} // Best primaries
	if (!freegamma) {//if Free gamma not selected	
	// gamma : ga0,ga1,ga2,ga3,ga4,ga5 by calcul
    if(gam=="BT709_g2.2_s4.5") 		{ga0=2.222;ga1=1./1.099258;ga2=0.099258/1.099258;ga3=1./4.5; ga4=0.01805;ga5=0.0;}//BT709  2.2  4.5  - my prefered as D.Coffin ga4=0.01805	
	else if (gam=="sRGB_g2.4_s12.92")	{ga0=Color::sRGBGammaCurve-0.0001; ga1=1./1.0550; ga2=0.0550/1.0550;ga3=1./12.92;ga4=0.039289;}//sRGB 2.4 12.92  - RT default as Lightroom
	else if (gam=="High_g1.3_s3.35")	{ga0=1.3 ; ga1=1./1.001724; ga2=0.001724/1.001724;ga3=1./3.35;ga4=0.001715;}//for high dynamic images
	else if (gam== "Low_g2.6_s6.9")   {ga0=2.6 ; ga1=1./1.12213; ga2=0.12213/1.12213;ga3=1./6.90;ga4=0.01;} //gamma 2.6 variable : for low contrast images
	else if (gam=="linear_g1.0")   {ga0=1.0; ga1=1.;ga2=0.;ga3=1./eps;ga4=0.;}//gamma=1 linear : for high dynamic images (cf : D.Coffin...)
	else if (gam=="standard_g2.2")   {ga0=2.2; ga1=1.;ga2=0.;ga3=1./eps;ga4=0.;}//gamma=2.2 (as gamma of Adobe, Widegamut...)
	else if (gam=="standard_g1.8")   {ga0=1.8; ga1=1.;ga2=0.;ga3=1./eps;ga4=0.;}//gamma=1.8  (as gamma of Prophoto)
	}
	else //free gamma selected
	{
	if(slpos==0) slpos=eps;
	calcGamma(pwr, ts, mode, imax,g_a0,g_a1,g_a2,g_a3,g_a4,g_a5);// call to calcGamma with selected gamma and slope : return parameters for LCMS2
	ga0=gampos;ga1=1./(1.0+g_a4);ga2=g_a4/(1.0 + g_a4);ga3=1./slpos; ga4=g_a3;ga5=0;
	}


	if(select_temp==1) t50=5003;// for Widegamut, Prophoto Best, Beta   D50
	else if (select_temp==2) t50=6504;// for sRGB, AdobeRGB, Bruce  D65

	cmsCIExyY       xyD;
	cmsCIExyYTRIPLE Primaries = {{p1, p2, 1.0},//red primaries
								{p3, p4, 1.0}, // green
								{p5, p6, 1.0} //blue
								};							   
    cmsToneCurve* GammaTRC[3];
	cmsFloat64Number Parameters[7];
    Parameters[0] = ga0;
    Parameters[1] = ga1;
    Parameters[2] = ga2;
    Parameters[3] = ga3;
    Parameters[4] = ga4;   
	Parameters[5] = ga5;   
	Parameters[6] = ga6;   
// 7 parameters for smoother curves
    cmsWhitePointFromTemp(&xyD, t50);
    GammaTRC[0] = GammaTRC[1] = GammaTRC[2] =   cmsBuildParametricToneCurve(NULL, 5, Parameters);//5 = more smoother than 4
    cmsHPROFILE oprofdef = cmsCreateRGBProfileTHR(NULL, &xyD, &Primaries, GammaTRC); //oprofdef  become Outputprofile

    cmsFreeToneCurve(GammaTRC[0]);

	
    if (oprofdef) {
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
            //float g;  // unused
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			short* xa = (short*)image->r[i-cy];
			short* ya = (short*)image->g[i-cy];
			short* za = (short*)image->b[i-cy];
			for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

				xa[j-cx] = CLIP((int)x_);
				ya[j-cx] = CLIP((int)y_);
				za[j-cx] = CLIP((int)z_);
			}
		}

        cmsHPROFILE iprof = iccStore->getXYZProfile ();
        lcmsMutex->lock ();
		cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16_PLANAR, oprofdef, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, cmsFLAGS_NOOPTIMIZE);
        lcmsMutex->unlock ();

		cmsDoTransform (hTransform, image->data, image->data, image->planestride);
		cmsDeleteTransform(hTransform);
	} else {
	// 
		#pragma omp parallel for if (multiThread)
		for (int i=cy; i<cy+ch; i++) {
            //float g;  // unused
			float R,G,B;
			float* rL = lab->L[i];
			float* ra = lab->a[i];
			float* rb = lab->b[i];
			for (int j=cx; j<cx+cw; j++) {
				
				float fy = (0.00862069 * rL[j])/327.68 + 0.137932; // (L+16)/116
				float fx = (0.002 * ra[j])/327.68 + fy;
				float fz = fy - (0.005 * rb[j])/327.68;
				
				float x_ = 65535.0 * f2xyz(fx)*D50x;
				float y_ = 65535.0 * f2xyz(fy);
				float z_ = 65535.0 * f2xyz(fz)*D50z;

				xyz2srgb(x_,y_,z_,R,G,B);

				image->r[i-cy][j-cx] = (int)gamma2curve[CLIP(R)];
				image->g[i-cy][j-cx] = (int)gamma2curve[CLIP(G)];
				image->b[i-cy][j-cx] = (int)gamma2curve[CLIP(B)];
			}
		}
	}
    return image;
}
	
//#include "sRGBgamutbdy.cc"

}
