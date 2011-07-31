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
#include <rtengine.h>
#include <improcfun.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <minmax.h>
#include <gauss.h>
#include <bilateral2.h>

namespace rtengine {

#undef CLIP
#undef CMAXVAL
#undef ABS

#define CMAXVAL 0xffff
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define ABS(a) ((a)<0?-(a):(a))


extern Settings* settings;

void ImProcFunctions::dcdamping (float** aI, float** aO, float damping, int W, int H) {
    const float dampingFac=2.0/(damping*damping);
#ifdef _OPENMP
#pragma omp for
#endif
	for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
            float I = aI[i][j];
            float O = aO[i][j];
            if (O<=0.0 || I<=0.0) {
                aI[i][j] = 0.0;
                continue;
            }
            float U = -(O * log(I/O) - I + O) * dampingFac;
            U = MIN(U,1.0);
            U = U*U*U*U*(5.0-U*4.0);
            aI[i][j] = (O - I) / I * U + 1.0;
        }
}

void ImProcFunctions::deconvsharpening (LabImage* lab, float** b2) {

    if (params->sharpening.enabled==false || params->sharpening.deconvamount<1)
        return;

    int W = lab->W, H = lab->H;

    float** tmpI = new float*[H];
    for (int i=0; i<H; i++) {
        tmpI[i] = new float[W];
        for (int j=0; j<W; j++)
            tmpI[i][j] = (float)lab->L[i][j];
    }

    float** tmp = (float**)b2;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(W,H));
    float damping = params->sharpening.deconvdamping / 5.0;
    bool needdamp = params->sharpening.deconvdamping > 0;
    for (int k=0; k<params->sharpening.deconviter; k++) {

    	// apply blur function (gaussian blur)
        gaussHorizontal<float> (tmpI, tmp, buffer, W, H, params->sharpening.deconvradius / scale, multiThread);
        gaussVertical<float>   (tmp, tmp,  buffer, W, H, params->sharpening.deconvradius / scale, multiThread);

    	if (!needdamp) {
#ifdef _OPENMP
#pragma omp for
#endif
            for (int i=0; i<H; i++)
                for (int j=0; j<W; j++)
                    if (tmp[i][j]>0)
                        tmp[i][j] = (float)lab->L[i][j] / tmp[i][j];
        }
        else
			dcdamping (tmp, lab->L, damping, W, H);

        gaussHorizontal<float> (tmp, tmp, buffer, W, H, params->sharpening.deconvradius / scale, multiThread);
        gaussVertical<float>   (tmp, tmp, buffer, W, H, params->sharpening.deconvradius / scale, multiThread);

#ifdef _OPENMP
#pragma omp for
#endif
        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++)
                tmpI[i][j] = tmpI[i][j] * tmp[i][j];
		} // end for
    delete buffer;

    float p2 = params->sharpening.deconvamount /100.0;
    float p1 = 1.0 - p2;

#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++)
            lab->L[i][j] = lab->L[i][j]*p1 + MAX(tmpI[i][j],0)*p2;

} // end parallel

    for (int i=0; i<H; i++)
        delete [] tmpI[i];
    delete [] tmpI;
}

void ImProcFunctions::sharpening (LabImage* lab, float** b2) {

    if (params->sharpening.method=="rld") {
       deconvsharpening (lab, b2);
        return;
    }

    // Rest is UNSHARP MASK
    if (params->sharpening.enabled==false || params->sharpening.amount<1 || lab->W<8 || lab->H<8)
        return;

    int W = lab->W, H = lab->H;
    float** b3;
    if (params->sharpening.edgesonly)
    {
    	b3 = new float*[H];
    	for (int i=0; i<H; i++)
    	   b3[i] = new float[W];
    }
#ifdef _OPENMP
#pragma omp parallel
#endif
    {


    AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(W,H));
    if (params->sharpening.edgesonly==false) {

        gaussHorizontal<float> (lab->L, b2, buffer, W, H, params->sharpening.radius / scale, multiThread);
        gaussVertical<float>   (b2,     b2, buffer, W, H, params->sharpening.radius / scale, multiThread);
    }
    else {
  		bilateral<float, float> (lab->L, (float**)b3, b2, W, H, params->sharpening.edges_radius / scale, params->sharpening.edges_tolerance, multiThread);
		gaussHorizontal<float> (b3, b2, buffer, W, H, params->sharpening.radius / scale, multiThread);
		gaussVertical<float>   (b2, b2, buffer, W, H, params->sharpening.radius / scale, multiThread);
    }
    delete buffer;

    float** base = lab->L;
    if (params->sharpening.edgesonly)
        base = b3;

    if (params->sharpening.halocontrol==false) {
		#pragma omp for
    	for (int i=0; i<H; i++)
            for (int j=0; j<W; j++) {
                float diff = base[i][j] - b2[i][j];
                if (ABS(diff)>params->sharpening.threshold) {
                    lab->L[i][j] = lab->L[i][j] + params->sharpening.amount * diff / 100.f;
                }
            }
    }
    else
		sharpenHaloCtrl (lab, b2, base, W, H);
    } // end parallel

    if (params->sharpening.edgesonly) {
        for (int i=0; i<H; i++)
            delete [] b3[i];
        delete [] b3;
    }
}

void ImProcFunctions::sharpenHaloCtrl (LabImage* lab, float** blurmap, float** base, int W, int H) {

    float scale = (100.f - params->sharpening.halocontrol_amount) * 0.01f;
    float sharpFac = params->sharpening.amount * 0.01f;
    float** nL = base;
#pragma omp parallel for if (multiThread)
    for (int i=2; i<H-2; i++) {
        float max1 = 0, max2 = 0, min1 = 0, min2 = 0, maxn, minn, np1, np2, np3, min, max, labL;
        for (int j=2; j<W-2; j++) {
            // compute 3 iterations, only forward
            np1 = 2.f * (nL[i-2][j] + nL[i-2][j+1] + nL[i-2][j+2] + nL[i-1][j] + nL[i-1][j+1] + nL[i-1][j+2] + nL[i]  [j] + nL[i]  [j+1] + nL[i]  [j+2]) / 27.f + nL[i-1][j+1] / 3.f;
            np2 = 2.f * (nL[i-1][j] + nL[i-1][j+1] + nL[i-1][j+2] + nL[i]  [j] + nL[i]  [j+1] + nL[i]  [j+2] + nL[i+1][j] + nL[i+1][j+1] + nL[i+1][j+2]) / 27.f + nL[i]  [j+1] / 3.f;
            np3 = 2.f * (nL[i]  [j] + nL[i]  [j+1] + nL[i]  [j+2] + nL[i+1][j] + nL[i+1][j+1] + nL[i+1][j+2] + nL[i+2][j] + nL[i+2][j+1] + nL[i+2][j+2]) / 27.f + nL[i+1][j+1] / 3.f;

            // Max/Min of all these deltas and the last two max/min
            MINMAX3(np1,np2,np3,maxn,minn);
            MAX3(max1,max2,maxn,max);
            MIN3(min1,min2,minn,min);

            // Shift the queue
            max1 = max2; max2 = maxn;
            min1 = min2; min2 = minn;
            labL = lab->L[i][j];
            if (max < labL) max = labL;
            if (min > labL) min = labL;

            // deviation from the environment as measurement
            float diff = nL[i][j] - blurmap[i][j];

            if (ABS(diff) > params->sharpening.threshold) {
                float newL = labL + sharpFac * diff;
                // applying halo control
                if (newL > max)
                    newL = max + (newL-max) * scale;
                else if (newL < min)
                    newL = min - (min-newL) * scale;

                lab->L[i][j] = newL;
            }
        }
    }
}
// To the extent possible under law, Manuel Llorens <manuelllorens@gmail.com>[
	// has waived all copyright and related or neighboring rights to this work.
// This work is published from: Spain.

//thanks to Manuel for this excellent job.. (Jacques Desmis JDC or frej83)
void ImProcFunctions::MLsharpen (LabImage* lab) {
// JD: this algorithm maximize clarity of images; it does not play on accutance. It can remove (partialy) the effects of the AA filter)
// I think we can use this algorithm alone in most cases, or first to clarify image and if you want a very little USM (unsharp mask sharpening) after...
    if (params->clarity.enabled==false)
        return;
	MyTime t1e,t2e;
	t1e.set();

	    int offset,c,i,j,p,width2;
		int width = lab->W, height = lab->H;
	    float *L,lumH,lumV,lumD1,lumD2,v,contrast,med,s;
	    float difL,difR,difT,difB,difLT,difRB,difLB,difRT,wH,wV,wD1,wD2,chmax[3];
	    float f1,f2,f3,f4;
	 	float templab;
		int iii,kkk;
	    width2=2*width;
		float strength;
		strength=params->clarity.clstrength / 100.0f;
		if(strength < 0.00001f) return;

		if (settings->verbose) printf ("Clarity strength %f\n", strength);
		
		L = new float[width*height];

	    chmax[0]=8.0f;
	    chmax[1]=3.0f;
	    chmax[2]=3.0f;
		
		int channels;
		if(params->clarity.clthreechannels) channels=0; else channels=2;
		if (settings->verbose) printf ("Clarity channels %d\n", channels);
		
		int passes=params->clarity.clpasses;
		if (settings->verbose) printf ("Clarity passes %d\n", passes);
		
	        for(p=0;p<passes;p++)
				for(c=0;c<=channels;c++) {// c=0 Luminance only

	        #pragma omp parallel for private(offset) shared(L)
	            for(offset=0;offset<width*height;offset++)
					{int ii=offset/width;
					int kk=offset-ii*width;
						if(c==0) L[offset]=lab->L[ii][kk]/327.68f; // adjust to RT and to 0..100
						else if (c==1) L[offset]=lab->a[ii][kk]/327.68f;
						else if (c==2) L[offset]=lab->b[ii][kk]/327.68f;
						}
	        #pragma omp parallel for private(j,i,iii,kkk, templab,offset,wH,wV,wD1,wD2,s,lumH,lumV,lumD1,lumD2,v,contrast,f1,f2,f3,f4,difT,difB,difL,difR,difLT,difLB,difRT,difRB) shared(lab,L,strength)
	            for(j=2;j<height-2;j++)
	                for(i=2,offset=j*width+i;i<width-2;i++,offset++){
	                    // weight functions
	                    wH=fabs(L[offset+1]-L[offset-1]);
	                    wV=fabs(L[offset+width]-L[offset-width]);  
	 
	                    s=1.0+fabs(wH-wV)/2.0;
	                    wD1=fabs(L[offset+width+1]-L[offset-width-1])/s;
	                    wD2=fabs(L[offset+width-1]-L[offset-width+1])/s;
	                    s=wD1;
	                    wD1/=wD2;
	                    wD2/=wD1;
	 
	                    // initial values
						int ii=offset/width;
						int kk=offset-ii*width;
					    if(c==0)lumH=lumV=lumD1=lumD2=v=lab->L[ii][kk]/327.68f;
						else if (c==1) lumH=lumV=lumD1=lumD2=v=lab->a[ii][kk]/327.68f;
						else if (c==2) lumH=lumV=lumD1=lumD2=v=lab->b[ii][kk]/327.68f;

	 
	                    // contrast detection
	                    contrast=sqrt(fabs(L[offset+1]-L[offset-1])*fabs(L[offset+1]-L[offset-1])+fabs(L[offset+width]-L[offset-width])*fabs(L[offset+width]-L[offset-width]))/chmax[c];
	                    if(contrast>1.0) contrast=1.0;
	 
	                    // new possible values
	                    if((L[offset]<L[offset-1])&&(L[offset]>L[offset+1])||(L[offset]>L[offset-1])&&(L[offset]<L[offset+1])){
	                        f1=fabs(L[offset-2]-L[offset-1]);
	                        f2=fabs(L[offset-1]-L[offset]);
	                        f3=fabs(L[offset-1]-L[offset-width])*fabs(L[offset-1]-L[offset+width]);
	                        f4=sqrt(fabs(L[offset-1]-L[offset-width2])*fabs(L[offset-1]-L[offset+width2]));
	                        difL=f1*f2*f2*f3*f3*f4;
	                        f1=fabs(L[offset+2]-L[offset+1]);
	                        f2=fabs(L[offset+1]-L[offset]);
	                        f3=fabs(L[offset+1]-L[offset-width])*fabs(L[offset+1]-L[offset+width]);
	                        f4=sqrt(fabs(L[offset+1]-L[offset-width2])*fabs(L[offset+1]-L[offset+width2]));
	                        difR=f1*f2*f2*f3*f3*f4;
	                        if((difR!=0)&&(difL!=0)){
	                            lumH=(L[offset-1]*difR+L[offset+1]*difL)/(difL+difR);
	                            lumH=v*(1-contrast)+lumH*contrast;
	                        }
	                    }
	 
	                    if((L[offset]<L[offset-width])&&(L[offset]>L[offset+width])||(L[offset]>L[offset-width])&&(L[offset]<L[offset+width])){
	                        f1=fabs(L[offset-width2]-L[offset-width]);
	                        f2=fabs(L[offset-width]-L[offset]);
	                        f3=fabs(L[offset-width]-L[offset-1])*fabs(L[offset-width]-L[offset+1]);
	                        f4=sqrt(fabs(L[offset-width]-L[offset-2])*fabs(L[offset-width]-L[offset+2]));
	                        difT=f1*f2*f2*f3*f3*f4;
	                        f1=fabs(L[offset+width2]-L[offset+width]);
	                        f2=fabs(L[offset+width]-L[offset]);
	                        f3=fabs(L[offset+width]-L[offset-1])*fabs(L[offset+width]-L[offset+1]);
	                        f4=sqrt(fabs(L[offset+width]-L[offset-2])*fabs(L[offset+width]-L[offset+2]));
	                        difB=f1*f2*f2*f3*f3*f4;
	                        if((difB!=0)&&(difT!=0)){
	                            lumV=(L[offset-width]*difB+L[offset+width]*difT)/(difT+difB);
	                            lumV=v*(1-contrast)+lumV*contrast;
	                        }
	                    }
	 
	                    if((L[offset]<L[offset-1-width])&&(L[offset]>L[offset+1+width])||(L[offset]>L[offset-1-width])&&(L[offset]<L[offset+1+width])){
	                        f1=fabs(L[offset-2-width2]-L[offset-1-width]);
	                        f2=fabs(L[offset-1-width]-L[offset]);
	                        f3=fabs(L[offset-1-width]-L[offset-width+1])*fabs(L[offset-1-width]-L[offset+width-1]);
	                        f4=sqrt(fabs(L[offset-1-width]-L[offset-width2+2])*fabs(L[offset-1-width]-L[offset+width2-2]));
	                        difLT=f1*f2*f2*f3*f3*f4;
	                        f1=fabs(L[offset+2+width2]-L[offset+1+width]);
	                        f2=fabs(L[offset+1+width]-L[offset]);
	                        f3=fabs(L[offset+1+width]-L[offset-width+1])*fabs(L[offset+1+width]-L[offset+width-1]);
	                        f4=sqrt(fabs(L[offset+1+width]-L[offset-width2+2])*fabs(L[offset+1+width]-L[offset+width2-2]));
	                        difRB=f1*f2*f2*f3*f3*f4;
	                        if((difLT!=0)&&(difRB!=0)){
	                            lumD1=(L[offset-1-width]*difRB+L[offset+1+width]*difLT)/(difLT+difRB);
	                            lumD1=v*(1-contrast)+lumD1*contrast;
	                        }
	                    }
	 
	                    if((L[offset]<L[offset+1-width])&&(L[offset]>L[offset-1+width])||(L[offset]>L[offset+1-width])&&(L[offset]<L[offset-1+width])){
	                        f1=fabs(L[offset-2+width2]-L[offset-1+width]);
	                        f2=fabs(L[offset-1+width]-L[offset]);
	                        f3=fabs(L[offset-1+width]-L[offset-width-1])*fabs(L[offset-1+width]-L[offset+width+1]);
	                        f4=sqrt(fabs(L[offset-1+width]-L[offset-width2-2])*fabs(L[offset-1+width]-L[offset+width2+2]));
                        difLB=f1*f2*f2*f3*f3*f4;
	                        f1=fabs(L[offset+2-width2]-L[offset+1-width]);
	                        f2=fabs(L[offset+1-width]-L[offset])*fabs(L[offset+1-width]-L[offset]);
	                        f3=fabs(L[offset+1-width]-L[offset+width+1])*fabs(L[offset+1-width]-L[offset-width-1]);
	                        f4=sqrt(fabs(L[offset+1-width]-L[offset+width2+2])*fabs(L[offset+1-width]-L[offset-width2-2]));
	                        difRT=f1*f2*f2*f3*f3*f4;
	                        if((difLB!=0)&&(difRT!=0)){
	                            lumD2=(L[offset+1-width]*difLB+L[offset-1+width]*difRT)/(difLB+difRT);
	                            lumD2=v*(1-contrast)+lumD2*contrast;
	                        }
	                    }
 
	                    s=strength;                
	 
	                    // avoid sharpening diagonals too much
	                    if(((fabs(wH/wV)<0.45f)&&(fabs(wH/wV)>0.05f))||((fabs(wV/wH)<0.45f)&&(fabs(wV/wH)>0.05f))) s=strength/3.0f;                                                                
 
						// final mix
	                    if((wH!=0.0f)&&(wV!=0.0f)&&(wD1!=0.0f)&&(wD2!=0.0f)) {
							iii=offset/width;
							kkk=offset-iii*width;
							templab=v*(1-s)+(lumH*wH+lumV*wV+lumD1*wD1+lumD2*wD2)/(wH+wV+wD1+wD2)*s;
							if(c==0) lab->L[iii][kkk]=fabs(327.68f*templab);// fabs because lab->L always >0
							else if (c==1){lab->a[iii][kkk]=327.68f*templab;}
							else if (c==2)lab->b[iii][kkk]=327.68f*templab;
							}
				 
	                }
					}

					delete [] L;
					
	t2e.set();
	if( settings->verbose )
		printf("Clarity gradient  %d usec\n", t2e.etime(t1e));
					
	        }

 // To the extent possible under law, Manuel Llorens <manuelllorens@gmail.com>  
 // has waived all copyright and related or neighboring rights to this work.  
  // This code is licensed under CC0 v1.0, see license information at  
  // http://creativecommons.org/publicdomain/zero/1.0/  
  // addition from JD : pyramid  + ponderated contrast with matrix 5x5
   void ImProcFunctions::MLmicrocontrast(LabImage* lab){ 
    if (params->clarity.enabledtwo==false)
        return;
	MyTime t1e,t2e;
	t1e.set();
	int k;
	   if(params->clarity.MLmicromatrix == false) k=2; else k=1;
		// k=2 matrix 5x5  k=1 matrix 3x3
		int offset,offset2,c,i,j,col,row,n; 
		float temp,temp2,temp3,temp4,tempL;	
		float *LM,v,s,contrast,w;  
		int signs[25];  
     	int width = lab->W, height = lab->H;
    	float uniform=params->clarity.uniformity;//between 0 to 100
		int unif;
		unif=(int)(uniform/10.0f); //put unif between 0 to 10
	    float strength=params->clarity.mlstrength/1500.0f; //strength 2000.0 quasi no artefacts ==> 1500 = maximum, after artefacts
		if(strength < 0.000001f) return;
		if(k==1) strength*=2.7f;//25/9 if 3x3
		if (settings->verbose) printf ("Microcontrast strength %f\n", strength);		
		if (settings->verbose) printf ("Microcontrast uniformity %i\n",unif);
		//modulation uniformity in function of luminance
		float L98[11]={0.001f,0.0015f,0.002f,0.004f,0.006f,0.008f,0.01f,0.03f,0.05f,0.1f,0.1f};
		float L95[11]={0.0012f,0.002f,0.005f,0.01f,0.02f,0.05f,0.1f,0.12f,0.15f,0.2f,0.25f};
		float L92[11]={0.01f,0.015f,0.02f,0.06f,0.10f,0.13f,0.17f,0.25f,0.3f,0.32f,0.35f};
		float L90[11]={0.015f,0.02f,0.04f,0.08f,0.12f,0.15f,0.2f,0.3f,0.4f,0.5f,0.6f};
		float L87[11]={0.025f,0.03f,0.05f,0.1f,0.15f,0.25f,0.3f,0.4f,0.5f,0.63f,0.75f};
		float L83[11]={0.055f,0.08f,0.1f,0.15f,0.2f,0.3f,0.4f,0.5f,0.6f,0.75f,0.85f};
		float L80[11]={0.15f,0.2f,0.25f,0.3f,0.35f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f};
		float L75[11]={0.22f,0.25f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.85f,0.9f,0.95f};
		float L70[11]={0.35f,0.4f,0.5f,0.6f,0.7f,0.8f,0.97f,1.0f,1.0f,1.0f,1.0f};
		float L63[11]={0.55f,0.6f,0.7f,0.8f,0.85f,0.9f,1.0f,1.0f,1.0f,1.0f,1.0f};
		float L58[11]={0.75f,0.77f,0.8f,0.9f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};
		//default 5
		//modulation contrast
		float Cont0[11]={0.05f,0.1f,0.2f,0.25f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f};
		float Cont1[11]={0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,0.95f,1.0f};
		float Cont2[11]={0.2f,0.40f,0.6f,0.7f,0.8f,0.85f,0.90f,0.95f,1.0f,1.05f,1.10f};
		float Cont3[11]={0.5f,0.6f,0.7f,0.8f,0.85f,0.9f,1.0f,1.0f,1.05f,1.10f,1.20f};
		float Cont4[11]={0.8f,0.85f,0.9f,0.95f,1.0f,1.05f,1.10f,1.150f,1.2f,1.25f,1.40f};
		float Cont5[11]={1.0f,1.1f,1.2f,1.25f,1.3f,1.4f,1.45f,1.50f,1.6f,1.65f,1.80f};
		
		float chmax=8.0f;
		LM = new float[width*height];//allocation for Luminance
       c=0;  
      #pragma omp parallel for private(offset, i,j) shared(LM)  
				       for(j=0;j<height;j++)  
							for(i=0,offset=j*width+i;i<width;i++,offset++){  
								LM[offset]=lab->L[j][i]/327.68f;// adjust to 0.100 and to RT variables
								}

     #pragma omp parallel for private(j,i,offset,s,signs,v,n,row,col,offset2,contrast,temp,w,temp2,temp3,tempL,temp4) shared(lab,LM,strength,chmax,unif,k,L98,L95,L92,L90,L87,L83,L80,L75,L70,L63,L58,Cont0,Cont1,Cont2,Cont3,Cont4,Cont5)  
       for(j=k;j<height-k;j++)  
           for(i=k,offset=j*width+i;i<width-k;i++,offset++){  
				s=strength;  
				v=LM[offset];              
              n=0;  
              for(row=j-k;row<=j+k;row++)  
                  for(col=i-k,offset2=row*width+col;col<=i+k;col++,offset2++){  
                      signs[n]=0;  
                      if(v<LM[offset2]) signs[n]=-1;  
                      if(v>LM[offset2]) signs[n]=1;  
                      n++;  
                 }     
          if(k==1) contrast=sqrt(fabs(LM[offset+1]-LM[offset-1])*fabs(LM[offset+1]-LM[offset-1])+fabs(LM[offset+width]-LM[offset-width])*fabs(LM[offset+width]-LM[offset-width]))/chmax; //for 3x3 
          else if(k==2) contrast=sqrt(fabs(LM[offset+1]-LM[offset-1])*fabs(LM[offset+1]-LM[offset-1])+fabs(LM[offset+width]-LM[offset-width])*fabs(LM[offset+width]-LM[offset-width])\
		  +fabs(LM[offset+2]-LM[offset-2])*fabs(LM[offset+2]-LM[offset-2])+fabs(LM[offset+2*width]-LM[offset-2*width])*fabs(LM[offset+2*width]-LM[offset-2*width]))/(2*chmax);  //for 5x5

			  if(contrast>1.0f) contrast=1.0f;             
					//matrix 5x5		
				temp=lab->L[j][i]/327.68f;  //begin 3x3
				temp +=(v-LM[offset-width-1])*sqrtf(2.0f)*s;
				temp +=(v-LM[offset-width])*s;
				temp +=(v-LM[offset-width+1])*sqrtf(2.0f)*s;
				temp +=(v-LM[offset-1])*s;
				temp +=(v-LM[offset+1])*s;
				temp +=(v-LM[offset+width-1])*sqrtf(2.0f)*s;
				temp +=(v-LM[offset+width])*s;
				temp +=(v-LM[offset+width+1])*sqrtf(2.0f)*s;//end 3x3
				
				// add JD continue 5x5
				if(k==2) {
				temp +=2.0f*(v-LM[offset+2*width])*s;
				temp +=2.0f*(v-LM[offset-2*width])*s;
				temp +=2.0f*(v-LM[offset-2])*s;
				temp +=2.0f*(v-LM[offset+2])*s;
				
				temp +=2.0f*(v-LM[offset+2*width -1])*s*sqrtf(1.25f);// 1.25  = 1*1 + 0.5*0.5
				temp +=2.0f*(v-LM[offset+2*width -2])*s*sqrtf(2.0f);
				temp +=2.0f*(v-LM[offset+2*width+1])*s*sqrtf(1.25f);;
				temp +=2.0f*(v-LM[offset+2*width+2])*s*sqrtf(2.0f);
				temp +=2.0f*(v-LM[offset+ width+2])*s*sqrtf(1.25f);;
				temp +=2.0f*(v-LM[offset+width-2])*s*sqrtf(1.25f);;
				temp +=2.0f*(v-LM[offset-2*width -1])*s*sqrtf(1.25f);
				temp +=2.0f*(v-LM[offset-2*width -2])*s*sqrtf(2.0f);
				temp +=2.0f*(v-LM[offset-2*width+1])*s*sqrtf(1.25f);;
				temp +=2.0f*(v-LM[offset-2*width+2])*s*sqrtf(2.0f);
				temp +=2.0f*(v-LM[offset- width+2])*s*sqrtf(1.25f);;
				temp +=2.0f*(v-LM[offset-width-2])*s*sqrtf(1.25f);;
				}
				if(temp <0.0f) temp=0.0f;
				v=temp;
               
				n=0;

              for(row=j-k;row<=j+k;row++)  
                 for(col=i-k,offset2=row*width+col;col<=i+k;col++,offset2++){  
                     if(((v<LM[offset2])&&(signs[n]>0))||((v>LM[offset2])&&(signs[n]<0)))
							{
							temp =v*0.75f+LM[offset2]*0.25f;// 0.75 0.25
						n++;  
							} 
						}
						if(LM[offset]>95.0f || LM[offset]<5.0f) contrast*=Cont0[unif]; //+ JD : luminance  pyramid to adjust contrast by evaluation of LM[offset]
						else if(LM[offset]>90.0f || LM[offset]<10.0f) contrast*=Cont1[unif];
						else if(LM[offset]>80.0f || LM[offset]<20.0f) contrast*=Cont2[unif];
						else if(LM[offset]>70.0f || LM[offset]<30.0f) contrast*=(2.0f/k)*Cont3[unif];
						else if(LM[offset]>60.0f || LM[offset]<40.0f) contrast*=(2.0f/k)*Cont4[unif];
						else 
						contrast*=(2.0f/k)*Cont5[unif];
						if(contrast>1.0f) {contrast=1.0f;}   
						tempL=327.68f*(temp*(1.0f-contrast)+LM[offset]*contrast); 
						// JD: modulation of microcontrast in function of original Luminance and modulation of luminance
						temp2=tempL/(327.68f*LM[offset]);//for highlights
						if(temp2>1.0f) {
						if(temp2>2.0f) temp2=2.0f;//limit action
						if(LM[offset]>98.0f) {lab->L[j][i]=LM[offset]*327.68f;}
						else if(LM[offset]>95.0f) {temp3=temp2-1.0f;temp=(L95[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>92.0f) {temp3=temp2-1.0f;temp=(L92[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>90.0f) {temp3=temp2-1.0f;temp=(L90[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>87.0f) {temp3=temp2-1.0f;temp=(L87[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>83.0f) {temp3=temp2-1.0f;temp=(L83[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>80.0f) {temp3=temp2-1.0f;temp=(L80[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}						
						else if(LM[offset]>75.0f) {temp3=temp2-1.0f;temp=(L75[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>70.0f) {temp3=temp2-1.0f;temp=(L70[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>63.0f) {temp3=temp2-1.0f;temp=(L63[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>58.0f) {temp3=temp2-1.0f;temp=(L58[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>42.0f) {temp3=temp2-1.0f;temp=(L58[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>37.0f) {temp3=temp2-1.0f;temp=(L63[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}						
						else if(LM[offset]>30.0f) {temp3=temp2-1.0f;temp=(L70[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>25.0f) {temp3=temp2-1.0f;temp=(L75[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>20.0f) {temp3=temp2-1.0f;temp=(L80[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>17.0f) {temp3=temp2-1.0f;temp=(L83[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>13.0f) {temp3=temp2-1.0f;temp=(L87[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>10.0f) {temp3=temp2-1.0f;temp=(L90[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>5.0f) {temp3=temp2-1.0f;temp=(L95[unif]*temp3)+1.0f;lab->L[j][i]=temp*LM[offset]*327.68f;}
						else if(LM[offset]>0.0f) {lab->L[j][i]=LM[offset]*327.68f;}						
						}
						temp4=(327.68f*LM[offset])/tempL;//
						if(temp4>1.0f) {
						if(temp4>2.f) temp4=2.f;//limit action
						if(LM[offset]<2.0f) {temp3=temp4-1.0f;temp=(L98[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<5.0f) {temp3=temp4-1.0f;temp=(L95[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<8.0f) {temp3=temp4-1.0f;temp=(L92[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}						
						else if(LM[offset]<10.0f) {temp3=temp4-1.0f;temp=(L90[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<13.0f) {temp3=temp4-1.0f;temp=(L87[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<17.0f) {temp3=temp4-1.0f;temp=(L83[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}					
						else if(LM[offset]<20.0f) {temp3=temp4-1.0f;temp=(L80[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}						
						else if(LM[offset]<25.0f) {temp3=temp4-1.0f;temp=(L75[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<30.0f) {temp3=temp4-1.0f;temp=(L70[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<37.0f) {temp3=temp4-1.0f;temp=(L63[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<42.0f) {temp3=temp4-1.0f;temp=(L58[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<58.0f) {temp3=temp4-1.0f;temp=(L58[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<63.0f) {temp3=temp4-1.0f;temp=(L63[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}						
						else if(LM[offset]<70.0f) {temp3=temp4-1.0f;temp=(L70[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<75.0f) {temp3=temp4-1.0f;temp=(L75[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<80.0f) {temp3=temp4-1.0f;temp=(L80[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}					
						else if(LM[offset]<83.0f) {temp3=temp4-1.0f;temp=(L83[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}						
						else if(LM[offset]<87.0f) {temp3=temp4-1.0f;temp=(L87[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<90.0f) {temp3=temp4-1.0f;temp=(L90[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<95.0f) {temp3=temp4-1.0f;temp=(L95[unif]*temp3)+1.0f;lab->L[j][i]=(LM[offset]*327.68f)/temp;}
						else if(LM[offset]<100.0f) {lab->L[j][i]=LM[offset]*327.68f;}										
						}
						
			}  
					delete [] LM;
		t2e.set();
		if( settings->verbose )
		printf("Microcontrast  %d usec\n", t2e.etime(t1e));
				
  }  

}
