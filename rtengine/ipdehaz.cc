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
 
  
  # include <stdlib.h>
    # include <stdio.h>
    # include <math.h>  
    # include <string.h>  
    # include <string.h>  
    #include "rtengine.h"
  
    #include "improcfun.h"
    # define MAX_DEHAZE_SCALES    8 
    # define clipdehaz( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )   
 
 namespace rtengine
{

extern const Settings* settings;
    
    static float DehazeScales[MAX_DEHAZE_SCALES];  
    
    typedef struct   
    {   
      int    N;   
      float  sigma;   
      double B;   
      double b[4];   
    } gauss3;   
    
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

void mean_stddv( float *dst, float &mean, float &stddv, int W_L, int H_L )   
    {   
    float vsquared;   
    int i, j;   
       
    vsquared = 0.0f;   
    mean = 0.0f; 
    for (int i=0; i<H_L*W_L; i++) {
       
           mean += dst[i]; 
           vsquared += (dst[i] * dst[i]);   
        }
    
    mean /= (float) W_L*H_L;  
    vsquared /= (float) W_L*H_L;    
    stddv = ( vsquared - (mean * mean) );   
    stddv = (float)sqrt(stddv);   
    }   

//we can probaly chnage this code with Gausshorizontal and Gaussvertical from Gauss.h

void  compute_coefs3( gauss3 * c, float sigma )   
    {   
      float q, q2, q3;   
       
      q = 0.f;   
       
      if ( sigma >= 2.5f )   
      {   
          q = 0.98711f * sigma - 0.96330f;   
      }   
      else if ( (sigma >= 0.5f) && (sigma < 2.5f) )   
      {   
          q = 3.97156f - 4.14554f * (float) sqrt ((double) 1 - 0.26891 * sigma);   
      }   
      else   
      {   
          q = 0.1147705018520355224609375f;   
      }   
       
      q2 = q * q;   
      q3 = q * q2;   
      c->b[0] = (1.57825f+(2.44413f*q)+(1.4281f *q2)+(0.422205f*q3));   
      c->b[1] = (         (2.44413f*q)+(2.85619f*q2)+(1.26661f *q3));   
      c->b[2] = (                     -((1.4281f*q2)+(1.26661f *q3)));   
      c->b[3] = (                                    (0.422205f*q3));   
      c->B = 1.0f-((c->b[1]+c->b[2]+c->b[3])/c->b[0]);   
      c->sigma = sigma;   
      c->N = 3;   
    }   


void  gausssmooth( float *in, float *out, int size, int rowstride, gauss3 *c )   
    {   
      int i,n, bufsize;   
      float *w1,*w2;   
       
      bufsize = size+3;   
      size -= 1;   
      w1 = new float [bufsize];
      w2 = new float [bufsize];
      w1[0] = in[0];   
      w1[1] = in[0];   
      w1[2] = in[0];   
      for ( i = 0 , n=3; i <= size ; i++, n++)   
      {   
         w1[n] = (float)(c->B*in[i*rowstride] +   
                       ((c->b[1]*w1[n-1] +   
                         c->b[2]*w1[n-2] +   
                         c->b[3]*w1[n-3] ) / c->b[0]));   
      }   
       
      w2[size+1]= w1[size+3];   
      w2[size+2]= w1[size+3];   
      w2[size+3]= w1[size+3];   
      for ( i = size, n = i; i >= 0; i--, n-- )   
      {   
         w2[n]= out[i * rowstride] = (float)(c->B*w1[n] +   
                                           ((c->b[1]*w2[n+1] +   
                                             c->b[2]*w2[n+2] +   
                                             c->b[3]*w2[n+3] ) / c->b[0]));   
      }      
    delete [] w1;
    delete [] w2;
    }   

void ImProcFunctions::MSR(LabImage* lab, int width, int height, int skip)
    
    {   
      int           Is;            
      float         weight;   
      gauss3        coef;   
      float         mean, stddv;   
      float         mini, delta, maxi;   
      float eps = 5.f;
      float gain = (float) params->labCurve.gain;//def =1  not use
      float offset  = (float) params->labCurve.offs;//def = 0  not use
      float strength = (float) params->labCurve.str;
      int scal =  params->labCurve.scal;//def=3
      int  nei = (int) 2.5f*params->labCurve.neigh;//def = 200
      int modedehaz;
      if(params->labCurve.dehazmet=="none") modedehaz=-1;//enabled disabled
      if(params->labCurve.dehazmet=="uni") modedehaz=0;
      if(params->labCurve.dehazmet=="low") modedehaz=1;
      if(params->labCurve.dehazmet=="high") modedehaz=2;
      if (modedehaz !=-1) {//enabled
        int H_L=height;
        int W_L=width;
        float *src = new float[H_L*W_L];
        memset( src, 0, H_L*W_L * sizeof (float) );   
        
        float *dst = new float[H_L*W_L];
        memset( dst, 0, H_L*W_L * sizeof (float) );   
       
        float *out = new float[H_L*W_L];
        memset( out, 0, H_L*W_L * sizeof (float) );   
        
        float *in = new float[H_L*W_L];
        memset( in, 0, H_L*W_L * sizeof (float) );   
        
        for (int i=0; i< H_L; i++) {
            for (int j=0; j<W_L; j++) {
                src[i*W_L + j]=lab->L[i][j] + eps; 
            }
        }    
        Is  = width * height ;   
        dehaze_scales( DehazeScales, scal, modedehaz, nei );   
       
        weight = 1.0f / (float) scal;   
       
        int posw = 0;   
        for (int i = 0; i < Is ; i++ )   
            {   
                in[i] = (float)(src[i] + eps);   //avoid log(0)
            }   
        for ( int scale = 0; scale < scal; scale++ )   
            {   
                compute_coefs3( &coef, DehazeScales[scale] );   
                for (int row = 0; row < height; row++ )   
                {   
                    posw =  row * width;   
                    gausssmooth( in + posw, out + posw, width, 1, &coef);   
                }   
                memcpy( in,  out, Is * sizeof(float) );   
                memset( out, 0  , Is * sizeof(float) );   
       
                for (int col = 0; col < width; col++ )   
                {   
                    posw = col;   
                    gausssmooth( in + posw, out + posw, height, width, &coef );   
                }   
       
                for ( int i = 0; i < Is; i++ )   
                {   
                    dst[i] += weight * (float)( log(src[i] + eps) - log(out[i]) );   
                }   
            }   
            delete [] in;
            delete [] out;
 
       
        // Ci(x,y)=log[a Ii(x,y)]-log[ Ei=1-s Ii(x,y)]    
       
            float beta  = 16384.0f;   
            for ( int i = 0; i < Is; i ++ )   
                { 
                    float logsrc = (float)log( (float) src[i] + eps );   
                    dst[i] = gain * ((float)(log(beta * (src[i] + eps)) - logsrc) * dst[i]) + offset;
                }   
     
            mean=0.f;stddv=0.f;
            mean_stddv( dst, mean, stddv, W_L, H_L); 
            float dd=1.25f;
            
            mini = mean - dd*stddv; 
            maxi = mean + dd*stddv;   
            delta = maxi - mini;   
                printf("maxi=%f mini=%f mean=%f std=%f delta=%f\n", maxi, mini, mean, stddv, delta);
     
            if ( !delta ) delta = 1.0f; 
                for(int i=0;i<W_L*H_L;i++) { 
                    int ii = i/W_L;
                    int jj = i-ii*W_L;

                    float cd = 32768.f * ( dst[i] - mini ) / delta;
                    src[i] = clipdehaz( cd, 0.f, 32768.f );                    
                    lab->L[ii][jj]=((100.f - strength)* lab->L[ii][jj] + strength * src[i])/100.f;
                }
       
        delete [] dst;
        delete [] src;
      
      } 
    }   

}