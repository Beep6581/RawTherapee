////////////////////////////////////////////////////////////////
//
//		Fast demosaicing algorythm
//
//		copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: August 26, 2010
//
//	fast_demo.cc is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include <cmath>
#include "rawimagesource.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"

/*#ifdef _OPENMP
#include <omp.h>
#endif*/

using namespace std;
using namespace rtengine;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LUTf RawImageSource::initInvGrad()
{
    LUTf invGrad (0x10000);

	//set up directional weight function
	for (int i=0; i<0x10000; i++)
		invGrad[i] = 1.0/SQR(1.0+i);

    return invGrad;
}

LUTf RawImageSource::invGrad = RawImageSource::initInvGrad();

void RawImageSource::fast_demosaic(int winx, int winy, int winw, int winh) {
	//int winx=0, winy=0;
	//int winw=W, winh=H;
	
	if (plistener) {
		plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::methodstring[RAWParams::fast]));
		plistener->setProgress (0.0);
	}
	float progress = 0.0;
	

    const int bord=4;
		
	int clip_pt = 4*65535*initialGain;
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
#ifdef _OPENMP
#pragma omp for
#endif
	//first, interpolate borders using bilinear
	for (int i=0; i<H; i++) {

        float sum[6];

		for (int j=0; j<bord; j++) {//first few columns
			for (int c=0; c<6; c++) sum[c]=0;
			for (int i1=i-1; i1<i+2; i1++)
				for (int j1=j-1; j1<j+2; j1++) {
					if ((i1 > -1) && (i1 < H) && (j1 > -1)) {
						int c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			int c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=rawData[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=rawData[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=rawData[i][j];
				}
			}
		}//j
		
		for (int j=W-bord; j<W; j++) {//last few columns
			for (int c=0; c<6; c++) sum[c]=0;
			for (int i1=i-1; i1<i+2; i1++)
				for (int j1=j-1; j1<j+2; j1++) {
					if ((i1 > -1) && (i1 < H ) && (j1 < W)) {
						int c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			int c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=rawData[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=rawData[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=rawData[i][j];
				}
			}
		}//j
	}//i
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef _OPENMP
#pragma omp for
#endif
	for (int j=bord; j<W-bord; j++) {
        float sum[6];

		for (int i=0; i<bord; i++) {//first few rows
			for (int c=0; c<6; c++) sum[c]=0;
			for (int i1=i-1; i1<i+2; i1++)
				for (int j1=j-1; j1<j+2; j1++) {
					if ((j1 > -1) && (j1 < W) && (i1 > -1)) {
						int c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			int c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=rawData[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=rawData[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=rawData[i][j];
				}
			}
		}//i
		
		for (int i=H-bord; i<H; i++) {//last few rows
			for (int c=0; c<6; c++) sum[c]=0;
			for (int i1=i-1; i1<i+2; i1++)
				for (int j1=j-1; j1<j+2; j1++) {
					if  ((j1 > -1) && (j1 < W) && (i1 < H)) {
						int c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			int c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=rawData[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=rawData[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=rawData[i][j];
				}
			}
		}//i
	}//j
	
	if(plistener) plistener->setProgress(0.05);
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef _OPENMP
#pragma omp for 
#endif
		// interpolate G using gradient weights
		for (int i=bord; i< H-bord; i++) {
			float	wtu, wtd, wtl, wtr;
			for (int j=bord; j < W-bord; j++) {
				
				if (FC(i,j)==1) {
					green[i][j] = rawData[i][j];
					//red[i][j] = green[i][j];
					//blue[i][j] = green[i][j];
					
				} else {
					//compute directional weights using image gradients
					wtu=invGrad[(abs(rawData[i+1][j]-rawData[i-1][j])+abs(rawData[i][j]-rawData[i-2][j])+abs(rawData[i-1][j]-rawData[i-3][j])) /4];
					wtd=invGrad[(abs(rawData[i-1][j]-rawData[i+1][j])+abs(rawData[i][j]-rawData[i+2][j])+abs(rawData[i+1][j]-rawData[i+3][j])) /4];
					wtl=invGrad[(abs(rawData[i][j+1]-rawData[i][j-1])+abs(rawData[i][j]-rawData[i][j-2])+abs(rawData[i][j-1]-rawData[i][j-3])) /4];
					wtr=invGrad[(abs(rawData[i][j-1]-rawData[i][j+1])+abs(rawData[i][j]-rawData[i][j+2])+abs(rawData[i][j+1]-rawData[i][j+3])) /4];

					//store in rgb array the interpolated G value at R/B grid points using directional weighted average
					green[i][j]=(wtu*rawData[i-1][j]+wtd*rawData[i+1][j]+wtl*rawData[i][j-1]+wtr*rawData[i][j+1]) / (wtu+wtd+wtl+wtr);
					//red[i][j] = green[i][j];
					//blue[i][j] = green[i][j];
					
				}
			}
			//progress+=(double)0.33/(H);
			//if(plistener) plistener->setProgress(progress);
		}
		if(plistener) plistener->setProgress(0.4);

#ifdef _OPENMP
#pragma omp for
#endif
		for (int i=bord; i< H-bord; i++) {
			for (int j=bord+(FC(i,2)&1); j < W-bord; j+=2) {
				
				int c=FC(i,j);
				//interpolate B/R colors at R/B sites
				
				if (c==0) {//R site
					red[i][j] = rawData[i][j];
					blue[i][j] = green[i][j] - 0.25f*((green[i-1][j-1]+green[i-1][j+1]+green[i+1][j+1]+green[i+1][j-1]) -
																 min(static_cast<float>(clip_pt),rawData[i-1][j-1]+rawData[i-1][j+1]+rawData[i+1][j+1]+rawData[i+1][j-1]));
				} else {//B site
					red[i][j] = green[i][j] - 0.25f*((green[i-1][j-1]+green[i-1][j+1]+green[i+1][j+1]+green[i+1][j-1]) -
					min(static_cast<float>(clip_pt),rawData[i-1][j-1]+rawData[i-1][j+1]+rawData[i+1][j+1]+rawData[i+1][j-1]));
					blue[i][j] = rawData[i][j];
				}
			}
			//progress+=(double)0.33/(H);
			//if(plistener) plistener->setProgress(progress);
		}
		if(plistener) plistener->setProgress(0.7);

#ifdef _OPENMP
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif

	// interpolate R/B using color differences
	for (int i=bord; i< H-bord; i++) {
		for (int j=bord+1-(FC(i,2)&1); j < W-bord; j+=2) {
			
			//interpolate R and B colors at G sites
			red[i][j] = green[i][j] - 0.25f*((green[i-1][j]-red[i-1][j])+(green[i+1][j]-red[i+1][j])+
													   (green[i][j-1]-red[i][j-1])+(green[i][j+1]-red[i][j+1]));
			blue[i][j] = green[i][j] - 0.25f*((green[i-1][j]-blue[i-1][j])+(green[i+1][j]-blue[i+1][j])+
														(green[i][j-1]-blue[i][j-1])+(green[i][j+1]-blue[i][j+1]));
		}
		progress+=(double)0.33/(H);
		//if(plistener) plistener->setProgress(progress);
	}
	if(plistener) plistener->setProgress(0.99);
	} // End of parallelization
	
#undef bord
	
}
