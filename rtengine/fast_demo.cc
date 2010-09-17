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



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void RawImageSource::fast_demo() { 
	
	if (plistener) {
		plistener->setProgressStr ("Fast demosaicing...");
		plistener->setProgress (0.0);
	}
	float progress = 0.0;
	
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
			for (c=0; c<6; c++) sum[c]=0;
			for (i1=i-1; i1<i+2; i1++)
				for (j1=j-1; j1<j+2; j1++) {
					if (i1 > -1 && i1 < H && j1 > -1) {
						c = FC(i1,j1);
						sum[c] += ri->data[i1][j1];
						sum[c+3]++;
					}
				}
			c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=ri->data[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=ri->data[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=ri->data[i][j];
				}
			}
		}//j
		
		for (j=W-bord; j<W; j++) {//last few columns
			for (c=0; c<6; c++) sum[c]=0;
			for (i1=i-1; i1<i+2; i1++)
				for (j1=j-1; j1<j+2; j1++) {
					if (i1 > -1 && i1 < H && j1 < W) {
						c = FC(i1,j1);
						sum[c] += ri->data[i1][j1];
						sum[c+3]++;
					}
				}
			c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=ri->data[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=ri->data[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=ri->data[i][j];
				}
			}
		}//j
	}//i
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	for (j=bord; j<W-bord; j++) {
		for (i=0; i<bord; i++) {//first few rows
			for (c=0; c<6; c++) sum[c]=0;
			for (i1=i-1; i1<i+2; i1++)
				for (j1=j-1; j1<j+2; j1++) {
					if (j1 > -1 && j1 < W && i1 > -1) {
						c = FC(i1,j1);
						sum[c] += ri->data[i1][j1];
						sum[c+3]++;
					}
				}
			c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=ri->data[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=ri->data[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=ri->data[i][j];
				}
			}
		}//i
		
		for (i=H-bord; i<H; i++) {//last few rows
			for (c=0; c<6; c++) sum[c]=0;
			for (i1=i-1; i1<i+2; i1++)
				for (j1=j-1; j1<j+2; j1++) {
					if (j1 > -1 && j1 < W && i1 < H) {
						c = FC(i1,j1);
						sum[c] += ri->data[i1][j1];
						sum[c+3]++;
					}
				}
			c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=ri->data[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=ri->data[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=ri->data[i][j];
				}
			}
		}//i
	}//j
	
	if(plistener) plistener->setProgress(0.05);
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
#pragma omp parallel private(i,j,c)
	{
		int rb;
		
		float	wt1, wt2, wt3, wt4;
		
		float * dirwt = new float [0x20000];
		
#pragma omp for schedule(dynamic) nowait

		//set up directional weight function
		for (int i=0; i<0x10000; i++) 
			dirwt[i] = 1.0/SQR(1.0+i); 
		
		
#pragma omp for schedule(dynamic) nowait
		
		// interpolate G using gradient weights
		for (i=bord; i< H-bord; i++) {
			for (j=bord; j < W-bord; j++) {
				
				if (FC(i,j)==1) {
					green[i][j] = ri->data[i][j];
					//red[i][j] = green[i][j];
					//blue[i][j] = green[i][j];
					
				} else {
					//compute directional weights using image gradients
					wt1=dirwt[(abs(ri->data[i+1][j]-ri->data[i-1][j])+abs(ri->data[i][j]-ri->data[i-2][j])+abs(ri->data[i-1][j]-ri->data[i-3][j])) >>4];
					wt2=dirwt[(abs(ri->data[i-1][j]-ri->data[i+1][j])+abs(ri->data[i][j]-ri->data[i+2][j])+abs(ri->data[i+1][j]-ri->data[i+3][j])) >>4];
					wt3=dirwt[(abs(ri->data[i][j+1]-ri->data[i][j-1])+abs(ri->data[i][j]-ri->data[i][j-2])+abs(ri->data[i][j-1]-ri->data[i][j-3])) >>4];
					wt4=dirwt[(abs(ri->data[i][j-1]-ri->data[i][j+1])+abs(ri->data[i][j]-ri->data[i][j+2])+abs(ri->data[i][j+1]-ri->data[i][j+3])) >>4];

					//store in rgb array the interpolated G value at R/B grid points using directional weighted average
					green[i][j]=(int)((wt1*ri->data[i-1][j]+wt2*ri->data[i+1][j]+wt3*ri->data[i][j-1]+wt4*ri->data[i][j+1])/(wt1+wt2+wt3+wt4));
					//red[i][j] = green[i][j];
					//blue[i][j] = green[i][j];
					
				}
			}
			progress+=(double)0.33/(H);
			//if(plistener) plistener->setProgress(progress);
		}
		if(plistener) plistener->setProgress(0.4);

		
#pragma omp for schedule(dynamic) nowait
		
		for (i=bord; i< H-bord; i++) {
			for (j=bord+(FC(i,2)&1); j < W-bord; j+=2) {
				
				c=FC(i,j);
				//interpolate B/R colors at R/B sites
				rb = CLIP((int)(green[i][j] - 0.25*((green[i-1][j-1]-ri->data[i-1][j-1])+(green[i-1][j+1]-ri->data[i-1][j+1])+  \
													(green[i+1][j+1]-ri->data[i+1][j+1])+(green[i+1][j-1]-ri->data[i+1][j-1]))));
				if (c==0) {//R site
					red[i][j] = ri->data[i][j];
					blue[i][j] = rb;
				} else {//B site
					red[i][j] = rb;
					blue[i][j] = ri->data[i][j];
				}
			}
			progress+=(double)0.33/(H);
			//if(plistener) plistener->setProgress(progress);
		}
		if(plistener) plistener->setProgress(0.7);

		
#pragma omp for schedule(dynamic) nowait
		
		// interpolate R/B using color differences
		for (i=bord; i< H-bord; i++) {
			for (j=bord+1-(FC(i,2)&1); j < W-bord; j+=2) {
				
				//interpolate R and B colors at G sites
				red[i][j] = CLIP((int)(green[i][j] - 0.25*((green[i-1][j]-red[i-1][j])+(green[i+1][j]-red[i+1][j])+ \
														   (green[i][j-1]-red[i][j-1])+(green[i][j+1]-red[i][j+1]))));
				blue[i][j] = CLIP((int)(green[i][j] - 0.25*((green[i-1][j]-blue[i-1][j])+(green[i+1][j]-blue[i+1][j])+ \
															(green[i][j-1]-blue[i][j-1])+(green[i][j+1]-blue[i][j+1]))));
			}
			progress+=(double)0.33/(H);
			//if(plistener) plistener->setProgress(progress);
		}
		if(plistener) plistener->setProgress(0.99);

	}
	
	
	
#undef bord
	

	
}//namespace
