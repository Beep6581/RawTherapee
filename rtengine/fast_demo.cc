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

#include <slicer.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace rtengine;

void RawImageSource::fast_demo() { 

	if (plistener) {
		plistener->setProgressStr ("Fast demosaicing...");
		plistener->setProgress (0.0);
	}

	//allocate output arrays
	red = new unsigned short*[H];
	for (int i=0; i<H; i++) {
		red[i] = new unsigned short[W];
	}
	green = new unsigned short*[H];
	for (int i=0; i<H; i++) {
		green[i] = new unsigned short[W];
	}	
	blue = new unsigned short*[H];
	for (int i=0; i<H; i++) {
		blue[i] = new unsigned short[W];
	}

#define bord 4

	int i, j, i1, j1, c, sum[6];
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//first, interpolate borders using bilinear
	for (i=0; i<H; i++) {
		for (j=0; j<bord; j++) {//first few columns
			for (i1=i-1; i1<i+2; i1++)
				for (j1=j-1; j1<j+2; j1++) {
					if (i1 > -1 && i1 < H && j1 > -1) {
						c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			c=FC(i,j);
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
		
		for (j=W-bord; j<W; j++) {//last few columns
			for (i1=i-1; i1<i+2; i1++)
				for (j1=j-1; j1<j+2; j1++) {
					if (i1 > -1 && i1 < H && j1 < W) {
						c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			c=FC(i,j);
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
	
	for (j=bord; j<W-bord; j++) {
		for (i=0; i<bord; i++) {//first few rows
			for (i1=i-1; i1<i+2; i1++)
				for (j1=j-1; j1<j+2; j1++) {
					if (j1 > -1 && j1 < W && i1 > -1) {
						c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			c=FC(i,j);
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
		
		for (i=H-bord; i<H; i++) {//last few rows
			for (i1=i-1; i1<i+2; i1++)
				for (j1=j-1; j1<j+2; j1++) {
					if (j1 > -1 && j1 < W && i1 < H) {
						c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			c=FC(i,j);
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
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Block *subRegion = new Block((unsigned int)(4),(unsigned int)(4),(unsigned int)(W-8), (unsigned int)(H-8));
	Slicer *slicedIm = new Slicer((unsigned int)W, (unsigned int)H, subRegion, 100000);

	//printf("\n\nBlock list\n\nDimension image = %dx%d\nDimension sous-region = %dx%d\nNombre demandé de pixel = %d,\nNombre maxi de pixel = %d\n\n", W, H, subRegion->width, subRegion->height, 100000, slicedIm->maxPixelNumber);

#pragma omp parallel private(i,j,c)
{
	int rb;

	float	wtu, wtd, wtl, wtr;

	Block *currBlock = new Block();
	
	#pragma omp for schedule(dynamic) nowait

	// interpolate G using gradient weights
	for (unsigned int blk = 0; blk < slicedIm->blockNumber; blk++) {
		slicedIm->get_block(blk, currBlock);
		//printf("#%d: currBlock->posX=%d, currBlock->posY=%d\n", blk, currBlock->posX, currBlock->posY);

		for (i=(int)currBlock->posY; i< (int)(currBlock->posY+currBlock->height); i++)
			for (j=(int)currBlock->posX; j < (int)(currBlock->posX+currBlock->width); j++) {

				if (FC(i,j)==1) {
					green[i][j] = rawData[i][j];
					//red[i][j] = green[i][j];
					//blue[i][j] = green[i][j];

				} else {
					//compute directional weights using image gradients
					wtu=1/SQR(1.0+fabs((int)rawData[i+1][j]-rawData[i-1][j])+fabs((int)rawData[i][j]-rawData[i-2][j])+fabs((int)rawData[i-1][j]-rawData[i-3][j]));
					wtd=1/SQR(1.0+fabs((int)rawData[i-1][j]-rawData[i+1][j])+fabs((int)rawData[i][j]-rawData[i+2][j])+fabs((int)rawData[i+1][j]-rawData[i+3][j]));
					wtl=1/SQR(1.0+fabs((int)rawData[i][j+1]-rawData[i][j-1])+fabs((int)rawData[i][j]-rawData[i][j-2])+fabs((int)rawData[i][j-1]-rawData[i][j-3]));
					wtr=1/SQR(1.0+fabs((int)rawData[i][j-1]-rawData[i][j+1])+fabs((int)rawData[i][j]-rawData[i][j+2])+fabs((int)rawData[i][j+1]-rawData[i][j+3]));

					//store in rgb array the interpolated G value at R/B grid points using directional weighted average
					green[i][j]=(int)((wtu*rawData[i-1][j]+wtd*rawData[i+1][j]+wtl*rawData[i][j-1]+wtr*rawData[i][j+1])/(wtu+wtd+wtl+wtr));
					//red[i][j] = green[i][j];
					//blue[i][j] = green[i][j];

				}
			}
	}
	
	#pragma omp for schedule(dynamic) nowait

	for (unsigned int blk = 0; blk < slicedIm->blockNumber; blk++) {
		slicedIm->get_block(blk, currBlock);

		for (i=(int)currBlock->posY; i< (int)(currBlock->posY+currBlock->height); i++)
			for (j=(int)currBlock->posX+(FC(i,currBlock->posX)&1); j < (int)(currBlock->posX+currBlock->width); j+=2) {

				c=FC(i,j);
				//interpolate B/R colors at R/B sites
				rb = CLIP((int)(green[i][j] - 0.25*((green[i-1][j-1]-rawData[i-1][j-1])+(green[i-1][j+1]-rawData[i-1][j+1])+  \
												  (green[i+1][j+1]-rawData[i+1][j+1])+(green[i+1][j-1]-rawData[i+1][j-1]))));
				if (c==0) {//R site
					red[i][j] = rawData[i][j];
					blue[i][j] = rb;
				} else {//B site
					red[i][j] = rb;
					blue[i][j] = rawData[i][j];
				}
			}
	}

	#pragma omp for schedule(dynamic) nowait

	// interpolate R/B using color differences
	for (unsigned int blk = 0; blk < slicedIm->blockNumber; blk++) {
		slicedIm->get_block(blk, currBlock);

		for (i=(int)currBlock->posY; i< (int)(currBlock->posY+currBlock->height); i++)
			for (j=(int)currBlock->posX+1-(FC(i,currBlock->posX)&1); j < (int)(currBlock->posX+currBlock->width); j+=2) {

				//interpolate R and B colors at G sites
				red[i][j] = CLIP((int)(green[i][j] - 0.25*((green[i-1][j]-red[i-1][j])+(green[i+1][j]-red[i+1][j])+ \
												 (green[i][j-1]-red[i][j-1])+(green[i][j+1]-red[i][j+1]))));
				blue[i][j] = CLIP((int)(green[i][j] - 0.25*((green[i-1][j]-blue[i-1][j])+(green[i+1][j]-blue[i+1][j])+ \
												 (green[i][j-1]-blue[i][j-1])+(green[i][j+1]-blue[i][j+1]))));
			}
	}

	delete currBlock;
}

	delete subRegion;
	delete slicedIm;

#undef bord
	
}
