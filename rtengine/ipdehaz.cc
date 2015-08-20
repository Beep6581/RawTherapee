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
 
     * D. J. Jobson, Z. Rahman, and G. A. Woodell. A multi-scale   
     * Retinex for bridging the gap between color images and the   
     * human observation of scenes. IEEE Transactions on Image Processing,  
     * 1997, 6(7): 965-976  
     * inspired from 2003 Fabien Pelisson <Fabien.Pelisson@inrialpes.fr>
 
 */





#include <stdlib.h>
#include <stdio.h>
#include <math.h>  
#include <string.h>  
#include <string.h>  
#include "rtengine.h"
#include "gauss.h"

#include "improcfun.h"
#define MAX_DEHAZE_SCALES   6
#define clipdehaz( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )   
 
 namespace rtengine
{

extern const Settings* settings;

static float DehazeScales[MAX_DEHAZE_SCALES];  

void  dehaze_scales( float* scales, int nscales, int mode, int s)   
    {
      if ( nscales == 1 )
      {
          scales[0] =  (float)s / 2.f;
      }
      else if (nscales == 2)   
      {
          scales[0] = (float) s / 2.f;
          scales[1] = (float) s;
      }
      else
      {
          float size_step = (float) s / (float) nscales;
            
            if(mode==0) {
                for (int i = 0; i < nscales; ++i )
                scales[i] = 2.0f + (float)i * size_step;
            }
            else if (mode==1) {
            size_step = (float)log(s - 2.0f) / (float) nscales;
                for (int i = 0; i < nscales; ++i )
                scales[i] = 2.0f + (float)pow (10.f, (i * size_step) / log (10.f));
            }
            else if(mode==2){
            size_step = (float) log(s - 2.0f) / (float) nscales;
                for ( int i = 0; i < nscales; ++i )
                scales[i] = s - (float)pow (10.f, (i * size_step) / log (10.f));
            }
    }
}

void mean_stddv( float **dst, float &mean, float &stddv, int W_L, int H_L )
    {
    float vsquared;
    int i, j;
    
    vsquared = 0.0f;
    mean = 0.0f;
             for (int i = 0; i <H_L; i++ )
                for (int j=0; j<W_L; j++) {
                    mean += dst[i][j];
                    vsquared += (dst[i][j] * dst[i][j]);
                }
    
    mean /= (float) W_L*H_L;
    vsquared /= (float) W_L*H_L;
    stddv = ( vsquared - (mean * mean) );
    stddv = (float)sqrt(stddv);
    }



void ImProcFunctions::MSR(LabImage* lab, int width, int height, int skip)
    
    {   
      float         pond;
      float         mean, stddv;
      float         mini, delta, maxi;
      float eps = 2.f;
      float gain = (float) params->labCurve.gain;//def =1  not use
      float offset  = 1.f;
      float neg=(float) params->labCurve.offs;//def = 0  not use
      float strength = (float) params->labCurve.str;
      int scal =  params->labCurve.scal;//def=3
      int  nei = (int) 2.5f*params->labCurve.neigh;//def = 200
      int vart= params->labCurve.vart;
      int modedehaz;
      if(params->labCurve.dehazmet=="none") modedehaz=-1;//enabled disabled
      if(params->labCurve.dehazmet=="uni") modedehaz=0;
      if(params->labCurve.dehazmet=="low") modedehaz=1;
      if(params->labCurve.dehazmet=="high") modedehaz=2;
      if (modedehaz !=-1) {//enabled
        int H_L=height;
        int W_L=width;
        float** src;
        src = new float*[H_L];
        for (int i = 0; i < H_L; i++) {
            src[i] = new float[W_L];
            memset( src[i], 0, W_L * sizeof (float) );
        }
        float** dst;
        dst = new float*[H_L];
        for (int i = 0; i < H_L; i++) {
            dst[i] = new float[W_L];
            memset( dst[i], 0, W_L * sizeof (float) );
        }
        float** in;
        in = new float*[H_L];
        for (int i = 0; i < H_L; i++) {
            in[i] = new float[W_L];
            memset( in[i], 0, W_L * sizeof (float) );
        }
        float** out;
        out = new float*[H_L];
        for (int i = 0; i < H_L; i++) {
            out[i] = new float[W_L];
            memset( out[i], 0, W_L * sizeof (float) );
        }
                
        for (int i=0; i< H_L; i++) {
            for (int j=0; j<W_L; j++) {
                src[i][j]=lab->L[i][j] + eps;
            }
        }
        dehaze_scales( DehazeScales, scal, modedehaz, nei );
       
        pond = 1.0f / (float) scal;
        
#ifdef _OPENMP
#pragma omp for
#endif
        for (int i = 0; i < H_L ; i++ )
             for (int j=0; j<W_L; j++)
                {
                    in[i][j] = (float)(src[i][j] + eps);   //avoid log(0)
                }
         
        for ( int scale = 0; scale < scal; scale++ )
            {
#ifdef _OPENMP
#pragma omp parallel
#endif 
                {
                AlignedBufferMP<double> buffer(max(W_L,H_L));
                gaussHorizontal<float> (in, out, buffer, W_L, H_L, DehazeScales[scale]);
                gaussVertical<float>   (out, out, buffer,W_L, H_L, DehazeScales[scale]);
                }
                for ( int i=0; i < H_L; i++)
                    for (int j=0; j < W_L; j++)
                    {
                        dst[i][j] +=  pond * (float)( log(src[i][j] + eps) - log(out[i][j]) );
                    }
            }

     for (int i = 0; i < H_L; i++) {
        delete [] in[i];
    }
    delete [] in;
      for (int i = 0; i < H_L; i++) {
        delete [] out[i];
    }
    delete [] out;





float beta=16384.0f;


#ifdef _OPENMP
#pragma omp for
#endif
            for (int i=0; i< H_L; i++ )
                for (int j=0; j<W_L; j++)
                {
                    float logsrc = (float)log( (float) src[i][j] + eps );
                    dst[i][j] = (gain * ((float)(log(beta * (src[i][j] + eps)) - logsrc) * dst[i][j]) + offset);
                }
     
            mean=0.f;stddv=0.f;
            mean_stddv( dst, mean, stddv, W_L, H_L);
            float nstdv=1.25f;
            
            mini = mean - 1.2f*stddv;
            maxi = mean + 1.2f*stddv;
        //    mini = neg;
            delta = maxi - mini;
         //   delta = vart;
                printf("maxi=%f mini=%f mean=%f std=%f delta=%f\n", maxi, mini, mean, stddv, delta);
     
            if ( !delta ) delta = 1.0f;
#ifdef _OPENMP
#pragma omp for
#endif
             for ( int i=0; i < H_L; i ++ )
                for (int j=0; j< W_L; j++) {
                    float cd = vart*32768.f * ( dst[i][j] - mini ) / delta;
                    src[i][j] = clipdehaz( cd, 0.f, 32768.f );
                    lab->L[i][j]=((100.f - strength)* lab->L[i][j] + strength * src[i][j])/100.f;
                }
                
        for (int i = 0; i < H_L; i++) {
            delete [] dst[i];
        }
        delete [] dst;
        for (int i = 0; i < H_L; i++) {
            delete [] src[i];
        }
        delete [] src;       
      } 
    }   

}

 