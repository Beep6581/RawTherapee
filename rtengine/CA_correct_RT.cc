////////////////////////////////////////////////////////////////
//
//		Chromatic Aberration Auto-correction
//
//		copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: June 14, 2010
//
//	CA_correct_RT.cc is free software: you can redistribute it and/or modify
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
int RawImageSource::LinEqSolve(int nDim, float* pfMatr, float* pfVect, float* pfSolution) 
{
//==============================================================================
// return 1 if system not solving, 0 if system solved
// nDim - system dimension
// pfMatr - matrix with coefficients
// pfVect - vector with free members
// pfSolution - vector with system solution
// pfMatr becames trianglular after function call
// pfVect changes after function call
//
// Developer: Henry Guennadi Levkin
//
//==============================================================================

	float fMaxElem;
	float fAcc;
	
	int i, j, k, m;
	
	for(k=0; k<(nDim-1); k++) {// base row of matrix
		// search of line with max element
		fMaxElem = fabs( pfMatr[k*nDim + k] );
		m = k;
		for (i=k+1; i<nDim; i++) {
			if(fMaxElem < fabs(pfMatr[i*nDim + k]) ) {
				fMaxElem = pfMatr[i*nDim + k];
				m = i;
			}
		}
		
		// permutation of base line (index k) and max element line(index m)
		if(m != k) {
			for(i=k; i<nDim; i++) {
				fAcc               = pfMatr[k*nDim + i];
				pfMatr[k*nDim + i] = pfMatr[m*nDim + i];
				pfMatr[m*nDim + i] = fAcc;
			}
			fAcc = pfVect[k];
			pfVect[k] = pfVect[m];
			pfVect[m] = fAcc;
		}
		
		if( pfMatr[k*nDim + k] == 0.) {
			//linear system has no solution
			return 1; // needs improvement !!!
		}
		
		// triangulation of matrix with coefficients
		for(j=(k+1); j<nDim; j++) {// current row of matrix
			fAcc = - pfMatr[j*nDim + k] / pfMatr[k*nDim + k];
			for(i=k; i<nDim; i++) {
				pfMatr[j*nDim + i] = pfMatr[j*nDim + i] + fAcc*pfMatr[k*nDim + i];
			}
			pfVect[j] = pfVect[j] + fAcc*pfVect[k]; // free member recalculation
		}
	}
	
	for(k=(nDim-1); k>=0; k--) {
		pfSolution[k] = pfVect[k];
		for(i=(k+1); i<nDim; i++) {
			pfSolution[k] -= (pfMatr[k*nDim + i]*pfSolution[i]);
		}
		pfSolution[k] = pfSolution[k] / pfMatr[k*nDim + k];
	}
	
	return 0;
}
//end of linear equation solver
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void RawImageSource::CA_correct_RT() { 

#define TS 256		// Tile size
//#define border 8
//#define border2	16
	
	#define PIX_SORT(a,b) { if ((a)>(b)) {temp=(a);(a)=(b);(b)=temp;} }
	#define SQR(x) ((x)*(x))

	const float clip_pt = ri->defgain; 
		
	// local variables
	int width=W, height=H;
	//temporary array to store simple interpolation of G
	float (*Gtmp);
	Gtmp = (float (*)) calloc ((height)*(width), sizeof *Gtmp);
	
	const int border=8;
	const int border2=16;
	//order of 2d polynomial fit (polyord), and numpar=polyord^2
	int polyord=4, numpar=16;
	//number of blocks used in the fit
	int numblox[3]={0,0,0};
	
	int rrmin, rrmax, ccmin, ccmax;
	int top, left, row, col;
	int rr, cc, c, indx, indx1, i, j, k, m, n, dir;
	//number of pixels in a tile contributing to the CA shift diagnostic
	int areawt[2][3];
	//direction of the CA shift in a tile
	int GRBdir[2][3];
	//offset data of the plaquette where the optical R/B data are sampled
	int offset[2][3];
	int	shifthfloor[3], shiftvfloor[3], shifthceil[3], shiftvceil[3];
	//number of tiles in the image
	int vblsz, hblsz, vblock, hblock, vz1, hz1;
	//int verbose=1;
	//flag indicating success or failure of polynomial fit
	int res;
	//shifts to location of vertical and diagonal neighbors
	const int v1=TS, v2=2*TS, /* v3=3*TS,*/ v4=4*TS;//, p1=-TS+1, p2=-2*TS+2, p3=-3*TS+3, m1=TS+1, m2=2*TS+2, m3=3*TS+3;
	
	float eps=1e-5;			//tolerance to avoid dividing by zero
	
	//adaptive weights for green interpolation
	float	wtu, wtd, wtl, wtr;
	//local quadratic fit to shift data within a tile
	float	coeff[2][3][3];
	//measured CA shift parameters for a tile
	float	CAshift[2][3];
	//polynomial fit coefficients
	float	polymat[3][2][256], shiftmat[3][2][16], fitparams[3][2][16];
	//residual CA shift amount within a plaquette
	float	shifthfrac[3], shiftvfrac[3];
	//temporary storage for median filter
	float	temp, p[9];
	//temporary parameters for tile CA evaluation
	float	gdiff, deltgrb;
	//interpolated G at edge of plaquette
	float	Ginthfloor, Ginthceil, Gint, RBint, gradwt;
	//interpolated color difference at edge of plaquette
	float	grbdiffinthfloor, grbdiffinthceil, grbdiffint, grbdiffold;
	//data for evaluation of block CA shift variance
	float	blockave[2][3]={{0,0,0},{0,0,0}}, blocksqave[2][3]={{0,0,0},{0,0,0}}, blockdenom[2][3]={{0,0,0},{0,0,0}}, blockvar[2][3];
	//low and high pass 1D filters of G in vertical/horizontal directions
	float	glpfh, glpfv;
	
	//max allowed CA shift
	const float bslim = 3.99;
	//gaussians for low pass filtering of G and R/B
	//static const float gaussg[5] = {0.171582, 0.15839, 0.124594, 0.083518, 0.0477063};//sig=2.5
	//static const float gaussrb[3] = {0.332406, 0.241376, 0.0924212};//sig=1.25
	
	//block CA shift values and weight assigned to block
	
	char		*buffer;			// TS*TS*16
	//rgb data in a tile
	float         (*rgb)[3];		// TS*TS*12
	//color differences
	float         (*grbdiff);		// TS*TS*4
	//green interpolated to optical sample points for R/B
	float         (*gshift);		// TS*TS*4
	//high pass filter for R/B in vertical direction
	float         (*rbhpfh);		// TS*TS*4
	//high pass filter for R/B in horizontal direction
	float         (*rbhpfv);		// TS*TS*4
	//low pass filter for R/B in horizontal direction
	float         (*rblpfh);		// TS*TS*4
	//low pass filter for R/B in vertical direction
	float         (*rblpfv);		// TS*TS*4
	//low pass filter for color differences in horizontal direction
	float         (*grblpfh);		// TS*TS*4
	//low pass filter for color differences in vertical direction
	float         (*grblpfv);		// TS*TS*4

	
	/* assign working space; this would not be necessary
	 if the algorithm is part of the larger pre-interpolation processing */
	buffer = (char *) malloc(11*sizeof(float)*TS*TS);
	//merror(buffer,"CA_correct()");
	memset(buffer,0,11*sizeof(float)*TS*TS);
	
	// rgb array
	rgb         = (float (*)[3])		buffer;
	grbdiff		= (float (*))			(buffer +	3*sizeof(float)*TS*TS);
	gshift		= (float (*))			(buffer +	4*sizeof(float)*TS*TS);
	rbhpfh		= (float (*))			(buffer +	5*sizeof(float)*TS*TS);
	rbhpfv		= (float (*))			(buffer +	6*sizeof(float)*TS*TS);
	rblpfh		= (float (*))			(buffer +	7*sizeof(float)*TS*TS);
	rblpfv		= (float (*))			(buffer +	8*sizeof(float)*TS*TS);
	grblpfh		= (float (*))			(buffer +	9*sizeof(float)*TS*TS);
	grblpfv		= (float (*))			(buffer +	10*sizeof(float)*TS*TS);

	
	if((height+border2)%(TS-border2)==0) vz1=1; else vz1=0;
    if((width+border2)%(TS-border2)==0) hz1=1; else hz1=0;
    
    vblsz=ceil((float)(height+border2)/(TS-border2)+2+vz1);
    hblsz=ceil((float)(width+border2)/(TS-border2)+2+hz1);
	
	//block CA shift values and weight assigned to block
	char		*buffer1;				// vblsz*hblsz*(3*2+1)
	float		(*blockwt);				// vblsz*hblsz
	float		(*blockshifts)[3][2];	// vblsz*hblsz*3*2 
	//float blockshifts[1000][3][2]; //fixed memory allocation
	//float blockwt[1000]; //fixed memory allocation
	
	buffer1 = (char *) malloc(vblsz*hblsz*(3*2+1)*sizeof(float));
	//merror(buffer1,"CA_correct()");
	memset(buffer1,0,vblsz*hblsz*(3*2+1)*sizeof(float));
	// block CA shifts
	blockwt		= (float (*))			(buffer1);
	blockshifts	= (float (*)[3][2])		(buffer1+(vblsz*hblsz*sizeof(float)));
		
	int vctr=0, hctr=0;
	
	// Main algorithm: Tile loop
	//#pragma omp parallel for shared(image,height,width) private(top,left,indx,indx1) schedule(dynamic)
	for (top=-border, vblock=1; top < height; top += TS-border2, vblock++) {
		hctr=0;
		vctr++;
		for (left=-border, hblock=1; left < width; left += TS-border2, hblock++) {
			hctr++;
			int bottom = MIN( top+TS,height+border);
			int right  = MIN(left+TS, width+border);
			int rr1 = bottom - top;
			int cc1 = right - left;
			//t1_init = clock();
			// rgb from input CFA data
			// rgb values should be floating point number between 0 and 1 
			// after white balance multipliers are applied 
			if (top<0) {rrmin=border;} else {rrmin=0;}
			if (left<0) {ccmin=border;} else {ccmin=0;}
			if (bottom>height) {rrmax=height-top;} else {rrmax=rr1;}
			if (right>width) {ccmax=width-left;} else {ccmax=cc1;}
			
			for (rr=rrmin; rr < rrmax; rr++)
				for (row=rr+top, cc=ccmin; cc < ccmax; cc++) {
					col = cc+left;
					c = FC(rr,cc);
					indx=row*width+col;
					indx1=rr*TS+cc;
					rgb[indx1][c] = (rawData[row][col])/65535.0f;
					//rgb[indx1][c] = image[indx][c]/65535.0f;//for dcraw implementation
				}
			
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//fill borders
			if (rrmin>0) {
				for (rr=0; rr<border; rr++) 
					for (cc=ccmin; cc<ccmax; cc++) {
						c = FC(rr,cc);
						rgb[rr*TS+cc][c] = rgb[(border2-rr)*TS+cc][c];
					}
			}
			if (rrmax<rr1) {
				for (rr=0; rr<border; rr++) 
					for (cc=ccmin; cc<ccmax; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+cc][c] = (rawData[(height-rr-2)][left+cc])/65535.0f;
						//rgb[(rrmax+rr)*TS+cc][c] = (image[(height-rr-2)*width+left+cc][c])/65535.0f;//for dcraw implementation
					}
			}
			if (ccmin>0) {
				for (rr=rrmin; rr<rrmax; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[rr*TS+cc][c] = rgb[rr*TS+border2-cc][c];
					}
			}
			if (ccmax<cc1) {
				for (rr=rrmin; rr<rrmax; rr++)
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[rr*TS+ccmax+cc][c] = (rawData[(top+rr)][(width-cc-2)])/65535.0f;
						//rgb[rr*TS+ccmax+cc][c] = (image[(top+rr)*width+(width-cc-2)][c])/65535.0f;//for dcraw implementation
					}
			}
			
			//also, fill the image corners
			if (rrmin>0 && ccmin>0) {
				for (rr=0; rr<border; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[(rr)*TS+cc][c] = (rawData[border2-rr][border2-cc])/65535.0f;
						//rgb[(rr)*TS+cc][c] = (rgb[(border2-rr)*TS+(border2-cc)][c]);//for dcraw implementation
					}
			}
			if (rrmax<rr1 && ccmax<cc1) {
				for (rr=0; rr<border; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+ccmax+cc][c] = (rawData[(height-rr-2)][(width-cc-2)])/65535.0f;
						//rgb[(rrmax+rr)*TS+ccmax+cc][c] = (image[(height-rr-2)*width+(width-cc-2)][c])/65535.0f;//for dcraw implementation
					}
			}
			if (rrmin>0 && ccmax<cc1) {
				for (rr=0; rr<border; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[(rr)*TS+ccmax+cc][c] = (rawData[(border2-rr)][(width-cc-2)])/65535.0f;
						//rgb[(rr)*TS+ccmax+cc][c] = (image[(border2-rr)*width+(width-cc-2)][c])/65535.0f;//for dcraw implementation
					}
			}
			if (rrmax<rr1 && ccmin>0) {
				for (rr=0; rr<border; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+cc][c] = (rawData[(height-rr-2)][(border2-cc)])/65535.0f;
						//rgb[(rrmax+rr)*TS+cc][c] = (image[(height-rr-2)*width+(border2-cc)][c])/65535.0f;//for dcraw implementation
					}
			}
			
			//end of border fill
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			
			for (j=0; j<2; j++) 
				for (k=0; k<3; k++)
					for (c=0; c<3; c+=2) {
						coeff[j][k][c]=0;
					}
			//end of initialization
			
			
			for (rr=3; rr < rr1-3; rr++)
				for (row=rr+top, cc=3, indx=rr*TS+cc; cc < cc1-3; cc++, indx++) {
					col = cc+left;
					c = FC(rr,cc);
					
					if (c!=1) {
						//compute directional weights using image gradients
						wtu=1/SQR(eps+fabs(rgb[(rr+1)*TS+cc][1]-rgb[(rr-1)*TS+cc][1])+fabs(rgb[(rr)*TS+cc][c]-rgb[(rr-2)*TS+cc][c])+fabs(rgb[(rr-1)*TS+cc][1]-rgb[(rr-3)*TS+cc][1]));
						wtd=1/SQR(eps+fabs(rgb[(rr-1)*TS+cc][1]-rgb[(rr+1)*TS+cc][1])+fabs(rgb[(rr)*TS+cc][c]-rgb[(rr+2)*TS+cc][c])+fabs(rgb[(rr+1)*TS+cc][1]-rgb[(rr+3)*TS+cc][1]));
						wtl=1/SQR(eps+fabs(rgb[(rr)*TS+cc+1][1]-rgb[(rr)*TS+cc-1][1])+fabs(rgb[(rr)*TS+cc][c]-rgb[(rr)*TS+cc-2][c])+fabs(rgb[(rr)*TS+cc-1][1]-rgb[(rr)*TS+cc-3][1]));
						wtr=1/SQR(eps+fabs(rgb[(rr)*TS+cc-1][1]-rgb[(rr)*TS+cc+1][1])+fabs(rgb[(rr)*TS+cc][c]-rgb[(rr)*TS+cc+2][c])+fabs(rgb[(rr)*TS+cc+1][1]-rgb[(rr)*TS+cc+3][1]));
						
						//store in rgb array the interpolated G value at R/B grid points using directional weighted average
						rgb[indx][1]=(wtu*rgb[indx-v1][1]+wtd*rgb[indx+v1][1]+wtl*rgb[indx-1][1]+wtr*rgb[indx+1][1])/(wtu+wtd+wtl+wtr);
					}
					if (row>-1 && row<height && col>-1 && col<width)
						Gtmp[row*width + col] = rgb[indx][1];
				}

			for (rr=4; rr < rr1-4; rr++)
				for (cc=4+(FC(rr,2)&1), indx=rr*TS+cc, c = FC(rr,cc); cc < cc1-4; cc+=2, indx+=2) {
					

					rbhpfv[indx] = SQR(fabs((rgb[indx][1]-rgb[indx][c])-(rgb[indx+v4][1]-rgb[indx+v4][c])) + \
										fabs((rgb[indx-v4][1]-rgb[indx-v4][c])-(rgb[indx][1]-rgb[indx][c])) - \
										fabs((rgb[indx-v4][1]-rgb[indx-v4][c])-(rgb[indx+v4][1]-rgb[indx+v4][c])));
					rbhpfh[indx] = SQR(fabs((rgb[indx][1]-rgb[indx][c])-(rgb[indx+4][1]-rgb[indx+4][c])) + \
										fabs((rgb[indx-4][1]-rgb[indx-4][c])-(rgb[indx][1]-rgb[indx][c])) - \
										fabs((rgb[indx-4][1]-rgb[indx-4][c])-(rgb[indx+4][1]-rgb[indx+4][c])));
					
					
					glpfv = 0.25*(2*rgb[indx][1]+rgb[indx+v2][1]+rgb[indx-v2][1]);
					glpfh = 0.25*(2*rgb[indx][1]+rgb[indx+2][1]+rgb[indx-2][1]);
					rblpfv[indx] = eps+fabs(glpfv - 0.25*(2*rgb[indx][c]+rgb[indx+v2][c]+rgb[indx-v2][c]));
					rblpfh[indx] = eps+fabs(glpfh - 0.25*(2*rgb[indx][c]+rgb[indx+2][c]+rgb[indx-2][c]));
					grblpfv[indx] = glpfv + 0.25*(2*rgb[indx][c]+rgb[indx+v2][c]+rgb[indx-v2][c]);
					grblpfh[indx] = glpfh + 0.25*(2*rgb[indx][c]+rgb[indx+2][c]+rgb[indx-2][c]);
				}
			
			for (c=0;c<3;c++) {areawt[0][c]=areawt[1][c]=0;}
			
			// along line segments, find the point along each segment that minimizes the color variance
			// averaged over the tile; evaluate for up/down and left/right away from R/B grid point 
			for (rr=rrmin+8; rr < rrmax-8; rr++)
				for (cc=ccmin+8+(FC(rr,2)&1), indx=rr*TS+cc, c = FC(rr,cc); cc < ccmax-8; cc+=2, indx+=2) {
					
					if (rgb[indx][c]>0.8*clip_pt || Gtmp[indx]>0.8*clip_pt) continue;

					//in linear interpolation, color differences are a quadratic function of interpolation position;
					//solve for the interpolation position that minimizes color difference variance over the tile
					
					//vertical
					gdiff=0.3125*(rgb[indx+TS][1]-rgb[indx-TS][1])+0.09375*(rgb[indx+TS+1][1]-rgb[indx-TS+1][1]+rgb[indx+TS-1][1]-rgb[indx-TS-1][1]);
					deltgrb=(rgb[indx][c]-rgb[indx][1])-0.5*((rgb[indx-v4][c]-rgb[indx-v4][1])+(rgb[indx+v4][c]-rgb[indx+v4][1]));
					
					gradwt=fabs(0.25*rbhpfv[indx]+0.125*(rbhpfv[indx+2]+rbhpfv[indx-2]) );//*(grblpfv[indx-v2]+grblpfv[indx+v2])/(eps+0.1*grblpfv[indx-v2]+rblpfv[indx-v2]+0.1*grblpfv[indx+v2]+rblpfv[indx+v2]);
					if (gradwt>eps) {
					coeff[0][0][c] += gradwt*deltgrb*deltgrb;
					coeff[0][1][c] += gradwt*gdiff*deltgrb;
					coeff[0][2][c] += gradwt*gdiff*gdiff;
					areawt[0][c]++;
					}
					
					//horizontal
					gdiff=0.3125*(rgb[indx+1][1]-rgb[indx-1][1])+0.09375*(rgb[indx+1+TS][1]-rgb[indx-1+TS][1]+rgb[indx+1-TS][1]-rgb[indx-1-TS][1]);
					deltgrb=(rgb[indx][c]-rgb[indx][1])-0.5*((rgb[indx-4][c]-rgb[indx-4][1])+(rgb[indx+4][c]-rgb[indx+4][1]));
					
					gradwt=fabs(0.25*rbhpfh[indx]+0.125*(rbhpfh[indx+v2]+rbhpfh[indx-v2]) );//*(grblpfh[indx-2]+grblpfh[indx+2])/(eps+0.1*grblpfh[indx-2]+rblpfh[indx-2]+0.1*grblpfh[indx+2]+rblpfh[indx+2]);
					if (gradwt>eps) {
					coeff[1][0][c] += gradwt*deltgrb*deltgrb;
					coeff[1][1][c] += gradwt*gdiff*deltgrb;
					coeff[1][2][c] += gradwt*gdiff*gdiff;
					areawt[1][c]++;
					}
					
					//	In Mathematica,
					//  f[x_]=Expand[Total[Flatten[
					//  ((1-x) RotateLeft[Gint,shift1]+x RotateLeft[Gint,shift2]-cfapad)^2[[dv;;-1;;2,dh;;-1;;2]]]]];
					//  extremum = -.5Coefficient[f[x],x]/Coefficient[f[x],x^2]	
				}
			for (c=0; c<3; c+=2){
				for (j=0; j<2; j++) {// vert/hor
					//printf("hblock %d vblock %d j %d c %d areawt %d \n",hblock,vblock,j,c,areawt[j][c]);
					//printf("hblock %d vblock %d j %d c %d areawt %d ",hblock,vblock,j,c,areawt[j][c]);

					if (areawt[j][c]>0) {
						CAshift[j][c]=coeff[j][1][c]/coeff[j][2][c];
						blockwt[vblock*hblsz+hblock]= areawt[j][c];//*coeff[j][2][c]/(eps+coeff[j][0][c]) ;
					} else {
						CAshift[j][c]=17.0;
						blockwt[vblock*hblsz+hblock]=0;
					}
					
					//printf("%f  \n",CAshift[j][c]);

					//CAshift[j][c]=coeff[j][1][c]/coeff[j][2][c];
					//blockwt[vblock*hblsz+hblock] = (float)(rr1-8)*(cc1-8)/4 * coeff[j][2][c]/(eps+coeff[j][0][c]) ;
					
					//data structure = CAshift[vert/hor][color]
					//j=0=vert, 1=hor
					
					
					if ((CAshift[j][c])<0) {
						GRBdir[j][c]=-1;
					} else {
						GRBdir[j][c]=1;
					}	
					offset[j][c]=floor(CAshift[j][c]);
					//offset gives NW corner of square containing the min; j=0=vert, 1=hor
					
					if (fabs(CAshift[j][c])<2.0) {
						blockave[j][c] += CAshift[j][c];
						blocksqave[j][c] += SQR(CAshift[j][c]);
						blockdenom[j][c] += 1;
					}
				}//vert/hor
			}//color
			
			
			
			/* CAshift[j][c] are the locations  
			 that minimize color difference variances; 
			 This is the approximate _optical_ location of the R/B pixels */

			for (c=0; c<3; c+=2) {
				//evaluate the shifts to the location that minimizes CA within the tile
				blockshifts[(vblock)*hblsz+hblock][c][0]=(CAshift[0][c]); //vert CA shift for R/B
				blockshifts[(vblock)*hblsz+hblock][c][1]=(CAshift[1][c]); //hor CA shift for R/B
				//data structure: blockshifts[blocknum][R/B][v/h]
			}
			
			if(plistener) plistener->setProgress(0.5*fabs((float)top/height));

		}
	}
	//end of diagnostic pass
	
	for (j=0; j<2; j++)
		for (c=0; c<3; c+=2) {
			if (blockdenom[j][c]) {
				blockvar[j][c] = blocksqave[j][c]/blockdenom[j][c]-SQR(blockave[j][c]/blockdenom[j][c]);
			} else {
				fprintf (stderr,"blockdenom vanishes");
				return;
			}
		}
	
	//printf ("tile variances %f %f %f %f \n",blockvar[0][0],blockvar[1][0],blockvar[0][2],blockvar[1][2] );
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//now prepare for CA correction pass
	//first, fill border blocks of blockshift array
	for (vblock=1; vblock<vblsz-1; vblock++) {//left and right sides
		for (c=0; c<3; c+=2) {
			for (i=0; i<2; i++) {
				blockshifts[vblock*hblsz][c][i]=blockshifts[(vblock)*hblsz+2][c][i];
				blockshifts[vblock*hblsz+hblsz-1][c][i]=blockshifts[(vblock)*hblsz+hblsz-3][c][i];
			}
		}
	}
	for (hblock=0; hblock<hblsz; hblock++) {//top and bottom sides
		for (c=0; c<3; c+=2) {
			for (i=0; i<2; i++) {
				blockshifts[hblock][c][i]=blockshifts[2*hblsz+hblock][c][i];
				blockshifts[(vblsz-1)*hblsz+hblock][c][i]=blockshifts[(vblsz-3)*hblsz+hblock][c][i];
			}
		}
	}
	//end of filling border pixels of blockshift array
	
	//initialize fit arrays
	for (i=0; i<256; i++) {polymat[0][0][i] = polymat[0][1][i] = polymat[2][0][i] = polymat[2][1][i] = 0;}
	for (i=0; i<16; i++) {shiftmat[0][0][i] = shiftmat[0][1][i] = shiftmat[2][0][i] = shiftmat[2][1][i] = 0;}
	
	for (vblock=1; vblock<vblsz-1; vblock++)
		for (hblock=1; hblock<hblsz-1; hblock++) {
			// block 3x3 median of blockshifts for robustness
			for (c=0; c<3; c+=2) {
				for (dir=0; dir<2; dir++) {
					p[0] = blockshifts[(vblock-1)*hblsz+hblock-1][c][dir];
					p[1] = blockshifts[(vblock-1)*hblsz+hblock][c][dir];
					p[2] = blockshifts[(vblock-1)*hblsz+hblock+1][c][dir];
					p[3] = blockshifts[(vblock)*hblsz+hblock-1][c][dir];
					p[4] = blockshifts[(vblock)*hblsz+hblock][c][dir];
					p[5] = blockshifts[(vblock)*hblsz+hblock+1][c][dir];
					p[6] = blockshifts[(vblock+1)*hblsz+hblock-1][c][dir];
					p[7] = blockshifts[(vblock+1)*hblsz+hblock][c][dir];
					p[8] = blockshifts[(vblock+1)*hblsz+hblock+1][c][dir];
					PIX_SORT(p[1],p[2]); PIX_SORT(p[4],p[5]); PIX_SORT(p[7],p[8]);
					PIX_SORT(p[0],p[1]); PIX_SORT(p[3],p[4]); PIX_SORT(p[6],p[7]);
					PIX_SORT(p[1],p[2]); PIX_SORT(p[4],p[5]); PIX_SORT(p[7],p[8]);
					PIX_SORT(p[0],p[3]); PIX_SORT(p[5],p[8]); PIX_SORT(p[4],p[7]);
					PIX_SORT(p[3],p[6]); PIX_SORT(p[1],p[4]); PIX_SORT(p[2],p[5]);
					PIX_SORT(p[4],p[7]); PIX_SORT(p[4],p[2]); PIX_SORT(p[6],p[4]);
					PIX_SORT(p[4],p[2]);
					blockshifts[(vblock)*hblsz+hblock][c][dir] = p[4];
					if (p[4]<0) {GRBdir[dir][c]=-1;} else {GRBdir[dir][c]=1;}
				}
				
				
				//if (verbose) fprintf (stderr,_("tile vshift hshift (%d %d %4f %4f)...\n"),vblock, hblock, blockshifts[(vblock)*hblsz+hblock][c][0], blockshifts[(vblock)*hblsz+hblock][c][1]);
				
				
				//now prepare coefficient matrix; use only data points within two std devs of zero
				if (SQR(blockshifts[(vblock)*hblsz+hblock][c][0])>4.0*blockvar[0][c] || SQR(blockshifts[(vblock)*hblsz+hblock][c][1])>4.0*blockvar[1][c]) continue;
				numblox[c] += 1;
				for (dir=0; dir<2; dir++) {
					for (i=0; i<polyord; i++) {
						for (j=0; j<polyord; j++) {
							for (m=0; m<polyord; m++)
								for (n=0; n<polyord; n++) {
									polymat[c][dir][numpar*(polyord*i+j)+(polyord*m+n)] += (float)pow((float)vblock,i+m)*pow((float)hblock,j+n)*blockwt[vblock*hblsz+hblock];
								}
							shiftmat[c][dir][(polyord*i+j)] += (float)pow((float)vblock,i)*pow((float)hblock,j)*blockshifts[(vblock)*hblsz+hblock][c][dir]*blockwt[vblock*hblsz+hblock];
						}
					}//monomials	
				}//dir 
				
			}//c
		}//blocks
	
	numblox[1]=MIN(numblox[0],numblox[2]);
	//if too few data points, restrict the order of the fit to linear
	if (numblox[1]<32) {
		polyord=2; numpar=4;
		if (numblox[1]< 10) {
			fprintf (stderr,("numblox = %d \n"),numblox[1]);
			return;
		}
	}
	
	//fit parameters to blockshifts
	for (c=0; c<3; c+=2)
		for (dir=0; dir<2; dir++) {
			res = LinEqSolve(numpar, polymat[c][dir], shiftmat[c][dir], fitparams[c][dir]);
			if (res) {
				fprintf (stderr,("CA correction pass failed -- can't solve linear equations for color %d direction %d...\n"),c,dir);
				return;
			}
		}
	//fitparams[polyord*i+j] gives the coefficients of (vblock^i hblock^j) in a polynomial fit for i,j<=4
	
	//end of initialization for CA correction pass
	
	
	// Main algorithm: Tile loop
	//#pragma omp parallel for shared(image,height,width) private(top,left,indx,indx1) schedule(dynamic)
	for (top=-border, vblock=1; top < height; top += TS-border2, vblock++)
		for (left=-border, hblock=1; left < width; left += TS-border2, hblock++) {
			int bottom = MIN( top+TS,height+border);
			int right  = MIN(left+TS, width+border);
			int rr1 = bottom - top;
			int cc1 = right - left;
			//t1_init = clock();
			// rgb from input CFA data
			// rgb values should be floating point number between 0 and 1 
			// after white balance multipliers are applied 
			if (top<0) {rrmin=border;} else {rrmin=0;}
			if (left<0) {ccmin=border;} else {ccmin=0;}
			if (bottom>height) {rrmax=height-top;} else {rrmax=rr1;}
			if (right>width) {ccmax=width-left;} else {ccmax=cc1;}
			
			
			for (rr=rrmin; rr < rrmax; rr++)
				for (row=rr+top, cc=ccmin; cc < ccmax; cc++) {
					col = cc+left;
					c = FC(rr,cc);
					indx=row*width+col;
					indx1=rr*TS+cc;	
					//rgb[indx1][c] = image[indx][c]/65535.0f;
					rgb[indx1][c] = (rawData[row][col])/65535.0f;
					//rgb[indx1][c] = image[indx][c]/65535.0f;//for dcraw implementation

					if ((c&1)==0) rgb[indx1][1] = Gtmp[indx];
				}
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//fill borders
			if (rrmin>0) {
				for (rr=0; rr<border; rr++) 
					for (cc=ccmin; cc<ccmax; cc++) {
						c = FC(rr,cc);
						rgb[rr*TS+cc][c] = rgb[(border2-rr)*TS+cc][c];
						rgb[rr*TS+cc][1] = rgb[(border2-rr)*TS+cc][1];
					}
			}
			if (rrmax<rr1) {
				for (rr=0; rr<border; rr++) 
					for (cc=ccmin; cc<ccmax; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+cc][c] = (rawData[(height-rr-2)][left+cc])/65535.0f;
						//rgb[(rrmax+rr)*TS+cc][c] = (image[(height-rr-2)*width+left+cc][c])/65535.0f;//for dcraw implementation

						rgb[(rrmax+rr)*TS+cc][1] = Gtmp[(height-rr-2)*width+left+cc];
					}
			}
			if (ccmin>0) {
				for (rr=rrmin; rr<rrmax; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[rr*TS+cc][c] = rgb[rr*TS+border2-cc][c];
						rgb[rr*TS+cc][1] = rgb[rr*TS+border2-cc][1];
					}
			}
			if (ccmax<cc1) {
				for (rr=rrmin; rr<rrmax; rr++)
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[rr*TS+ccmax+cc][c] = (rawData[(top+rr)][(width-cc-2)])/65535.0f;
						//rgb[rr*TS+ccmax+cc][c] = (image[(top+rr)*width+(width-cc-2)][c])/65535.0f;//for dcraw implementation

						rgb[rr*TS+ccmax+cc][1] = Gtmp[(top+rr)*width+(width-cc-2)];
					}
			}
			
			//also, fill the image corners
			if (rrmin>0 && ccmin>0) {
				for (rr=0; rr<border; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[(rr)*TS+cc][c] = (rawData[border2-rr][border2-cc])/65535.0f;
						//rgb[(rr)*TS+cc][c] = (rgb[(border2-rr)*TS+(border2-cc)][c]);//for dcraw implementation

						rgb[(rr)*TS+cc][1] = Gtmp[(border2-rr)*width+border2-cc];
					}
			}
			if (rrmax<rr1 && ccmax<cc1) {
				for (rr=0; rr<border; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+ccmax+cc][c] = (rawData[(height-rr-2)][(width-cc-2)])/65535.0f;
						//rgb[(rrmax+rr)*TS+ccmax+cc][c] = (image[(height-rr-2)*width+(width-cc-2)][c])/65535.0f;//for dcraw implementation

						rgb[(rrmax+rr)*TS+ccmax+cc][1] = Gtmp[(height-rr-2)*width+(width-cc-2)];
					}
			}
			if (rrmin>0 && ccmax<cc1) {
				for (rr=0; rr<border; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[(rr)*TS+ccmax+cc][c] = (rawData[(border2-rr)][(width-cc-2)])/65535.0f;
						//rgb[(rr)*TS+ccmax+cc][c] = (image[(border2-rr)*width+(width-cc-2)][c])/65535.0f;//for dcraw implementation

						rgb[(rr)*TS+ccmax+cc][1] = Gtmp[(border2-rr)*width+(width-cc-2)];
					}
			}
			if (rrmax<rr1 && ccmin>0) {
				for (rr=0; rr<border; rr++) 
					for (cc=0; cc<border; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+cc][c] = (rawData[(height-rr-2)][(border2-cc)])/65535.0f;
						//rgb[(rrmax+rr)*TS+cc][c] = (image[(height-rr-2)*width+(border2-cc)][c])/65535.0f;//for dcraw implementation

						rgb[(rrmax+rr)*TS+cc][1] = Gtmp[(height-rr-2)*width+(border2-cc)];
					}
			}
			
			//end of border fill
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			
			blockshifts[(vblock)*hblsz+hblock][0][0] = blockshifts[(vblock)*hblsz+hblock][0][1] = 0;
			blockshifts[(vblock)*hblsz+hblock][2][0] = blockshifts[(vblock)*hblsz+hblock][2][1] = 0;
			for (i=0; i<polyord; i++)
				for (j=0; j<polyord; j++) {
					blockshifts[(vblock)*hblsz+hblock][0][0] += (float)pow((float)vblock,i)*pow((float)hblock,j)*fitparams[0][0][polyord*i+j];
					blockshifts[(vblock)*hblsz+hblock][0][1] += (float)pow((float)vblock,i)*pow((float)hblock,j)*fitparams[0][1][polyord*i+j];
					blockshifts[(vblock)*hblsz+hblock][2][0] += (float)pow((float)vblock,i)*pow((float)hblock,j)*fitparams[2][0][polyord*i+j];
					blockshifts[(vblock)*hblsz+hblock][2][1] += (float)pow((float)vblock,i)*pow((float)hblock,j)*fitparams[2][1][polyord*i+j];
				}
			blockshifts[(vblock)*hblsz+hblock][0][0] = LIM(blockshifts[(vblock)*hblsz+hblock][0][0], -bslim, bslim);
			blockshifts[(vblock)*hblsz+hblock][0][1] = LIM(blockshifts[(vblock)*hblsz+hblock][0][1], -bslim, bslim);
			blockshifts[(vblock)*hblsz+hblock][2][0] = LIM(blockshifts[(vblock)*hblsz+hblock][2][0], -bslim, bslim);
			blockshifts[(vblock)*hblsz+hblock][2][1] = LIM(blockshifts[(vblock)*hblsz+hblock][2][1], -bslim, bslim);
			
			
			for (c=0; c<3; c+=2) {
				
				//some parameters for the bilinear interpolation
				shiftvfloor[c]=floor((float)blockshifts[(vblock)*hblsz+hblock][c][0]);
				shiftvceil[c]=ceil((float)blockshifts[(vblock)*hblsz+hblock][c][0]);
				shiftvfrac[c]=blockshifts[(vblock)*hblsz+hblock][c][0]-shiftvfloor[c];
				
				shifthfloor[c]=floor((float)blockshifts[(vblock)*hblsz+hblock][c][1]);
				shifthceil[c]=ceil((float)blockshifts[(vblock)*hblsz+hblock][c][1]);
				shifthfrac[c]=blockshifts[(vblock)*hblsz+hblock][c][1]-shifthfloor[c];
			}
			
			
			for (rr=4; rr < rr1-4; rr++)
				for (cc=4+(FC(rr,2)&1), c = FC(rr,cc); cc < cc1-4; cc+=2) {
					//perform CA correction using color ratios or color differences
					
					Ginthfloor=(1-shifthfrac[c])*rgb[(rr+shiftvfloor[c])*TS+cc+shifthfloor[c]][1]+(shifthfrac[c])*rgb[(rr+shiftvfloor[c])*TS+cc+shifthceil[c]][1];
					Ginthceil=(1-shifthfrac[c])*rgb[(rr+shiftvceil[c])*TS+cc+shifthfloor[c]][1]+(shifthfrac[c])*rgb[(rr+shiftvceil[c])*TS+cc+shifthceil[c]][1];
					//Gint is blinear interpolation of G at CA shift point
					Gint=(1-shiftvfrac[c])*Ginthfloor+(shiftvfrac[c])*Ginthceil;
					
					//determine R/B at grid points using color differences at shift point plus interpolated G value at grid point
					//but first we need to interpolate G-R/G-B to grid points...
					grbdiff[(rr)*TS+cc]=Gint-rgb[(rr)*TS+cc][c];
					gshift[(rr)*TS+cc]=Gint;
					float tmp1=grbdiff[(rr)*TS+cc];

				}
			
			for (rr=8; rr < rr1-8; rr++)
				for (cc=8+(FC(rr,2)&1), c = FC(rr,cc), indx=rr*TS+cc; cc < cc1-8; cc+=2, indx+=2) {
					
					//if (rgb[indx][c]>clip_pt || Gtmp[indx]>clip_pt) continue;

					grbdiffold = rgb[indx][1]-rgb[indx][c];

					//interpolate color difference from optical R/B locations to grid locations
					grbdiffinthfloor=(1-shifthfrac[c]/2)*grbdiff[indx]+(shifthfrac[c]/2)*grbdiff[indx-2*GRBdir[1][c]];
					grbdiffinthceil=(1-shifthfrac[c]/2)*grbdiff[(rr-2*GRBdir[0][c])*TS+cc]+(shifthfrac[c]/2)*grbdiff[(rr-2*GRBdir[0][c])*TS+cc-2*GRBdir[1][c]];
					//grbdiffint is bilinear interpolation of G-R/G-B at grid point
					grbdiffint=(1-shiftvfrac[c]/2)*grbdiffinthfloor+(shiftvfrac[c]/2)*grbdiffinthceil;
					
					//now determine R/B at grid points using interpolated color differences and interpolated G value at grid point
					RBint=rgb[indx][1]-grbdiffint;
					
					if (fabs(RBint-rgb[indx][c])<0.25*(RBint+rgb[indx][c])) {
						if (fabs(grbdiffold)>fabs(grbdiffint) ) {
							rgb[indx][c]=RBint;
						}
					} else {
						
						//gradient weights using difference from G at CA shift points and G at grid points
						p[0]=1/(eps+fabs(rgb[indx][1]-gshift[indx]));
						p[1]=1/(eps+fabs(rgb[indx][1]-gshift[indx-2*GRBdir[1][c]]));
						p[2]=1/(eps+fabs(rgb[indx][1]-gshift[(rr-2*GRBdir[0][c])*TS+cc]));
						p[3]=1/(eps+fabs(rgb[indx][1]-gshift[(rr-2*GRBdir[0][c])*TS+cc-2*GRBdir[1][c]]));
						
						grbdiffint = (p[0]*grbdiff[indx]+p[1]*grbdiff[indx-2*GRBdir[1][c]]+ \
									  p[2]*grbdiff[(rr-2*GRBdir[0][c])*TS+cc]+p[3]*grbdiff[(rr-2*GRBdir[0][c])*TS+cc-2*GRBdir[1][c]])/(p[0]+p[1]+p[2]+p[3]);
						
						//now determine R/B at grid points using interpolated color differences and interpolated G value at grid point
						if (fabs(grbdiffold)>fabs(grbdiffint) ) {
							rgb[indx][c]=rgb[indx][1]-grbdiffint;
						}
					}
					
					//if color difference interpolation overshot the correction, just desaturate
					if (grbdiffold*grbdiffint<0) {
						rgb[indx][c]=rgb[indx][1]-0.5*(grbdiffold+grbdiffint);
					}
				}
			
			// copy CA corrected results back to image matrix
			for (rr=border; rr < rr1-border; rr++)
				for (row=rr+top, cc=border+(FC(rr,2)&1); cc < cc1-border; cc+=2) {
					col = cc + left;
					indx = row*width + col;
					c = FC(row,col);
					 
					rawData[row][col] = CLIP((int)(65535.0f*rgb[(rr)*TS+cc][c] + 0.5f));
					//image[indx][c] = CLIP((int)(65535.0*rgb[(rr)*TS+cc][c] + 0.5));//for dcraw implementation

				} 
			
			if(plistener) plistener->setProgress(0.5+0.5*fabs((float)top/height));

		}
	
	// clean up
	free(buffer);
	free(Gtmp);
	free(buffer1);
	

	
#undef TS
//#undef border
//#undef border2
#undef PIX_SORT
#undef SQR

}



