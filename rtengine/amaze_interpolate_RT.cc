////////////////////////////////////////////////////////////////
//
//			AMaZE demosaic algorithm
// (Aliasing Minimization and Zipper Elimination)
//
//	copyright (c) 2008-2010  Emil Martinec <ejmartin@uchicago.edu>
//
// incorporating ideas of Luis Sanz Rodrigues and Paul Lee
//
// code dated: May 27, 2010
//
//	amaze_interpolate_RT.cc is free software: you can redistribute it and/or modify
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



void RawImageSource::amaze_demosaic_RT() {  
	
#define SQR(x) ((x)*(x))
	//#define MIN(a,b) ((a) < (b) ? (a) : (b))
	//#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define LIM(x,min,max) MAX(min,MIN(x,max))
#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
	//#define CLIP(x) LIM(x,0,65535)
	
	
	int width=W, height=H;
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
	
#define TS 512	 // Tile size; the image is processed in square tiles to lower memory requirements and facilitate multi-threading
	
	// local variables

	int top, left;
	int ex, ey;
	int rrmin, rrmax, ccmin, ccmax;
	static const int v1=TS, v2=2*TS, v3=3*TS, p1=-TS+1, p2=-2*TS+2, p3=-3*TS+3, m1=TS+1, m2=2*TS+2, m3=3*TS+3;
	int nbr[5] = {-v2,-2,2,v2,0};
	
	static const float eps=1e-10;			//tolerance to avoid dividing by zero
	static const float arthresh=0.75;
	static const float nyqthresh=0.5;//0.5
	static const float pmthresh=0.25;//0.25
	
	static const float gaussodd[4] = {0.14659727707323927f, 0.103592713382435f, 0.0732036125103057f, 0.0365543548389495f};//gaussian on 5x5 quincunx, sigma=1.2
	float gaussgrad[6] = {0.07384411893421103f, 0.06207511968171489f, 0.0521818194747806f, \
	0.03687419286733595f, 0.03099732204057846f, 0.018413194161458882f};//gaussian on 5x5, sigma=1.2
	static const float gauss1[3] = {0.3376688223162362f, 0.12171198028231786f, 0.04387081413862306f};//gaussian on 3x3, sigma =0.7
	static const float gausseven[2] = {0.13719494435797422f, 0.05640252782101291f};//gaussian on 5x5 alt quincunx, sigma=1.5
	static const float gquinc[4] = {0.169917f, 0.108947f, 0.069855f, 0.0287182f};
	
	
	
	//tile vars	
	int bottom, right, row, col;
	int rr, cc, rr1, cc1, c, indx, indx1, dir, i, j, sgn;
	
	
	float cru, crd, crl, crr;
	float vwt, hwt, Gintv, Ginth;
	float guar, gdar, glar, grar, guha, gdha, glha, grha, Ginthar, Ginthha, Gintvar, Gintvha, hcdaltvar, vcdaltvar;
	float Dgrbvvaru, Dgrbvvard, Dgrbhvarl, Dgrbhvarr;
	float sumh, sumv, sumsqh, sumsqv, areawt;
	float nyqtest, vcdvar, hcdvar, hvwtalt, vo, ve, gradp, gradm, gradv, gradh, gradpm, gradhv;
	float vcdvar1, hcdvar1, varwt, diffwt;
	float rbvarp, rbvarm, crp, crm, rbp, rbm;
	float gu, gd, gl, gr;
	float gvarh, gvarv;
	float g[4], f[4];
	
	
	char		*buffer;			// TS*TS*168
	float         (*rgb)[3];		// TS*TS*12
	float         (*delh);			// TS*TS*4
	float         (*delv);			// TS*TS*4
	float         (*delhsq);		// TS*TS*4
	float         (*delvsq);		// TS*TS*4
	float         (*vcd);			// TS*TS*4
	float         (*hcd);			// TS*TS*4
	float         (*vcdalt);		// TS*TS*4
	float         (*hcdalt);		// TS*TS*4
	float         (*vcdsq);			// TS*TS*4
	float         (*hcdsq);			// TS*TS*4
	float         (*cddiffsq);		// TS*TS*4
	float         (*hvwt);			// TS*TS*4
	float         (*Dgrb)[2];		// TS*TS*8
	float         (*delp);			// TS*TS*4
	float         (*delm);			// TS*TS*4
	float         (*rbint);			// TS*TS*4
	float         (*dirwts)[2];		// TS*TS*8
	float         (*Dgrbh1);		// TS*TS*4
	float         (*Dgrbv1);		// TS*TS*4
	float         (*Dgrbhsq1);		// TS*TS*4
	float         (*Dgrbvsq1);		// TS*TS*4
	float         (*Dgrbh2);		// TS*TS*4
	float         (*Dgrbv2);		// TS*TS*4
	float         (*dgintv);		// TS*TS*4
	float         (*dginth);		// TS*TS*4
	float         (*Dgrbp1);		// TS*TS*4
	float         (*Dgrbm1);		// TS*TS*4
	float         (*Dgrbpsq1);		// TS*TS*4
	float         (*Dgrbmsq1);		// TS*TS*4
	float         (*cfa);			// TS*TS*4
	
	int			(*nyquist);		// TS*TS*4
	
	
	// assign working space
	buffer = (char *) malloc(144*TS*TS);
	//merror(buffer,"amaze_interpolate()");
	memset(buffer,0,144*TS*TS);
	// rgb array
	rgb         = (float (*)[3])		buffer; //pointers to array
 	delh		= (float (*))			(buffer +  12*TS*TS);
	delv		= (float (*))			(buffer +  16*TS*TS);
	delhsq		= (float (*))			(buffer +  20*TS*TS);
	delvsq		= (float (*))			(buffer +  24*TS*TS);
	vcd			= (float (*))			(buffer +  28*TS*TS);
	hcd			= (float (*))			(buffer +  32*TS*TS);
	vcdalt		= (float (*))			(buffer +  36*TS*TS);
	hcdalt		= (float (*))			(buffer +  40*TS*TS);
	vcdsq		= (float (*))			(buffer +  44*TS*TS);
	hcdsq		= (float (*))			(buffer +  48*TS*TS);
	cddiffsq	= (float (*))			(buffer +  52*TS*TS);
	hvwt		= (float (*))			(buffer +  56*TS*TS);
	Dgrb		= (float (*)[2])		(buffer +  60*TS*TS);
	delp		= (float (*))			(buffer +  68*TS*TS);
	delm		= (float (*))			(buffer +  72*TS*TS);
	rbint		= (float (*))			(buffer +  76*TS*TS);
	dirwts		= (float (*)[2])		(buffer +  80*TS*TS);
	Dgrbh1		= (float (*))			(buffer +  88*TS*TS);
	Dgrbv1		= (float (*))			(buffer +  92*TS*TS);
	Dgrbhsq1	= (float (*))			(buffer +  96*TS*TS);
	Dgrbvsq1	= (float (*))			(buffer +  100*TS*TS);
	Dgrbh2		= (float (*))			(buffer +  104*TS*TS);
	Dgrbv2		= (float (*))			(buffer +  108*TS*TS);	
	dgintv		= (float (*))			(buffer +  112*TS*TS);
	dginth		= (float (*))			(buffer +  116*TS*TS);
	Dgrbp1		= (float (*))			(buffer +  120*TS*TS);
	Dgrbm1		= (float (*))			(buffer +  124*TS*TS);
	Dgrbpsq1	= (float (*))			(buffer +  128*TS*TS);
	Dgrbmsq1	= (float (*))			(buffer +  132*TS*TS);
	cfa			= (float (*))			(buffer +  136*TS*TS);
	
	nyquist		= (int (*))				(buffer +  140*TS*TS);
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	/*double dt;
	 clock_t t1, t2;
	 
	 clock_t t1_init,       t2_init       = 0;
	 clock_t t1_vcdhcd,      t2_vcdhcd      = 0;
	 clock_t t1_cdvar,		t2_cdvar = 0;
	 clock_t t1_nyqtest,   t2_nyqtest   = 0;
	 clock_t t1_areainterp,  t2_areainterp  = 0;
	 clock_t t1_compare,   t2_compare   = 0;
	 clock_t t1_diag,   t2_diag   = 0;
	 clock_t t1_chroma,    t2_chroma    = 0;*/
	
	
	// start
	//if (verbose) fprintf (stderr,_("AMaZE interpolation ...\n"));
	//t1 = clock();
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if (plistener) {
		plistener->setProgressStr ("AMaZE Demosaicing...");
		plistener->setProgress (0.0);
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	//determine GRBG coset
	if (FC(0,0)==1) {//first pixel is G
		if (FC(0,1)==0) {ey=0; ex=1;} else {ey=1; ex=0;}
	} else {//first pixel is R or B
		if (FC(0,0)==0) {ey=0; ex=0;} else {ey=1; ex=1;}
	}
	
	// Main algorithm: Tile loop
	//#pragma omp parallel for shared(ri->data,height,width,red,green,blue) private(top,left) schedule(dynamic)
	//code is openmp ready; just have to pull local tile variable declarations inside the tile loop
	for (top=-16; top < height; top += TS-32)
		for (left=-16; left < width; left += TS-32) {
			bottom = MIN( top+TS,height+16);
			right  = MIN(left+TS, width+16);
			rr1 = bottom - top;
			cc1 = right - left;
			// rgb from input CFA data
			// rgb values should be floating point number between 0 and 1 
			// after white balance multipliers are applied 
			if (top<0) {rrmin=16;} else {rrmin=0;}
			if (left<0) {ccmin=16;} else {ccmin=0;}
			if (bottom>height) {rrmax=height-top;} else {rrmax=rr1;}
			if (right>width) {ccmax=width-left;} else {ccmax=cc1;}
			
			for (rr=rrmin; rr < rrmax; rr++)
				for (row=rr+top, cc=ccmin; cc < ccmax; cc++) {
					col = cc+left;
					c = FC(rr,cc);
					indx=row*width+col;
					indx1=rr*TS+cc;
					rgb[indx1][c] = (ri->data[row][col])/65535.0f;
					cfa[indx1] = rgb[indx1][c];
				}
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//fill borders
			if (rrmin>0) {
				for (rr=0; rr<16; rr++) 
					for (cc=ccmin; cc<ccmax; cc++) {
						c = FC(rr,cc);
						rgb[rr*TS+cc][c] = rgb[(32-rr)*TS+cc][c];
						cfa[rr*TS+cc] = rgb[rr*TS+cc][c];
					}
			}
			if (rrmax<rr1) {
				for (rr=0; rr<16; rr++) 
					for (cc=ccmin; cc<ccmax; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+cc][c] = (ri->data[(height-rr-2)][left+cc])/65535.0f;
						cfa[(rrmax+rr)*TS+cc] = rgb[(rrmax+rr)*TS+cc][c];
					}
			}
			if (ccmin>0) {
				for (rr=rrmin; rr<rrmax; rr++) 
					for (cc=0; cc<16; cc++) {
						c=FC(rr,cc);
						rgb[rr*TS+cc][c] = rgb[rr*TS+32-cc][c];
						cfa[rr*TS+cc] = rgb[rr*TS+cc][c];
					}
			}
			if (ccmax<cc1) {
				for (rr=rrmin; rr<rrmax; rr++)
					for (cc=0; cc<16; cc++) {
						c=FC(rr,cc);
						rgb[rr*TS+ccmax+cc][c] = (ri->data[(top+rr)][(width-cc-2)])/65535.0f;
						cfa[rr*TS+ccmax+cc] = rgb[rr*TS+ccmax+cc][c];
					}
			}
			
			//also fill the image corners
			if (rrmin>0 && ccmin>0) {
				for (rr=0; rr<16; rr++) 
					for (cc=0; cc<16; cc++) {
						c=FC(rr,cc);
						rgb[(rr)*TS+cc][c] = (ri->data[32-rr][32-cc])/65535.0f;
						cfa[(rr)*TS+cc] = rgb[(rr)*TS+cc][c];
					}
			}
			if (rrmax<rr1 && ccmax<cc1) {
				for (rr=0; rr<16; rr++) 
					for (cc=0; cc<16; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+ccmax+cc][c] = (ri->data[(height-rr-2)][(width-cc-2)])/65535.0f;
						cfa[(rrmax+rr)*TS+ccmax+cc] = rgb[(rrmax+rr)*TS+ccmax+cc][c];
					}
			}
			if (rrmin>0 && ccmax<cc1) {
				for (rr=0; rr<16; rr++) 
					for (cc=0; cc<16; cc++) {
						c=FC(rr,cc);
						rgb[(rr)*TS+ccmax+cc][c] = (ri->data[(rr)][(width-cc-2)])/65535.0f;
						cfa[(rr)*TS+ccmax+cc] = rgb[(rr)*TS+ccmax+cc][c];
					}
			}
			if (rrmax<rr1 && ccmin>0) {
				for (rr=0; rr<16; rr++) 
					for (cc=0; cc<16; cc++) {
						c=FC(rr,cc);
						rgb[(rrmax+rr)*TS+cc][c] = (ri->data[(height-rr-2)][cc])/65535.0f;
						cfa[(rrmax+rr)*TS+cc] = rgb[(rrmax+rr)*TS+cc][c];
					}
			}
			 
			//end of border fill
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			for (rr=1; rr < rr1-1; rr++)
				for (cc=1, indx=(rr)*TS+cc; cc < cc1-1; cc++, indx++) {
					
					delh[indx] = fabs(cfa[indx+1]-cfa[indx-1]);
					delv[indx] = fabs(cfa[indx+v1]-cfa[indx-v1]);
					delhsq[indx] = SQR(delh[indx]);
					delvsq[indx] = SQR(delv[indx]);
					delp[indx] = fabs(cfa[indx+p1]-cfa[indx-p1]);
					delm[indx] = fabs(cfa[indx+m1]-cfa[indx-m1]);
					
				}
			
			for (rr=2; rr < rr1-2; rr++)
				for (cc=2,indx=(rr)*TS+cc; cc < cc1-2; cc++, indx++) {
					
					dirwts[indx][0] = eps+delv[indx+v1]+delv[indx-v1]+delv[indx];//+fabs(cfa[indx+v2]-cfa[indx-v2]);
					//vert directional averaging weights
					dirwts[indx][1] = eps+delh[indx+1]+delh[indx-1]+delh[indx];//+fabs(cfa[indx+2]-cfa[indx-2]);
					//horizontal weights
					
					if (FC(rr,cc)&1) {
						//for later use in diagonal interpolation
						Dgrbp1[indx]=2*cfa[indx]-(cfa[indx-p1]+cfa[indx+p1]);
						Dgrbm1[indx]=2*cfa[indx]-(cfa[indx-m1]+cfa[indx+m1]);
						Dgrbpsq1[indx]=(SQR(cfa[indx]-cfa[indx-p1])+SQR(cfa[indx]-cfa[indx+p1]));
						Dgrbmsq1[indx]=(SQR(cfa[indx]-cfa[indx-m1])+SQR(cfa[indx]-cfa[indx+m1]));
					} 
				}
			
			//t2_init += clock()-t1_init;
			// end of tile initialization
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			//interpolate vertical and horizontal color differences
			//t1_vcdhcd = clock();
			
			for (rr=4; rr<TS-4; rr++)
				//for (cc=4+(FC(rr,2)&1),indx=rr*TS+cc,c=FC(rr,cc); cc<TS-4; cc+=2,indx+=2) {
				for (cc=4,indx=rr*TS+cc; cc<TS-4; cc++,indx++) {
					c=FC(rr,cc);
					if (c&1) {sgn=-1;} else {sgn=1;}
					
					//initialization of nyquist test
					nyquist[indx]=0;
					//preparation for diag interp
					rbint[indx]=0;
					
					//color ratios in each cardinal direction
					cru = cfa[indx-v1]*(dirwts[indx-v2][0]+dirwts[indx][0])/(eps+dirwts[indx-v2][0]*cfa[indx]+dirwts[indx][0]*cfa[indx-v2]);
					crd = cfa[indx+v1]*(dirwts[indx+v2][0]+dirwts[indx][0])/(eps+dirwts[indx+v2][0]*cfa[indx]+dirwts[indx][0]*cfa[indx+v2]);
					crl = cfa[indx-1]*(dirwts[indx-2][1]+dirwts[indx][1])/(eps+dirwts[indx-2][1]*cfa[indx]+dirwts[indx][1]*cfa[indx-2]);
					crr = cfa[indx+1]*(dirwts[indx+2][1]+dirwts[indx][1])/(eps+dirwts[indx+2][1]*cfa[indx]+dirwts[indx][1]*cfa[indx+2]);
					
					guha=cfa[indx-v1]+0.5*(cfa[indx]-cfa[indx-v2]);
					gdha=cfa[indx+v1]+0.5*(cfa[indx]-cfa[indx+v2]);
					glha=cfa[indx-1]+0.5*(cfa[indx]-cfa[indx-2]);
					grha=cfa[indx+1]+0.5*(cfa[indx]-cfa[indx+2]);
					
					if (fabs(1-cru)<arthresh) {guar=cfa[indx]*cru;} else {guar=guha;}
					if (fabs(1-crd)<arthresh) {gdar=cfa[indx]*crd;} else {gdar=gdha;}
					if (fabs(1-crl)<arthresh) {glar=cfa[indx]*crl;} else {glar=glha;}
					if (fabs(1-crr)<arthresh) {grar=cfa[indx]*crr;} else {grar=grha;}
					
					hwt = dirwts[indx-1][1]/(dirwts[indx-1][1]+dirwts[indx+1][1]);
					vwt = dirwts[indx-v1][0]/(dirwts[indx+v1][0]+dirwts[indx-v1][0]);
					
					//interpolated G via adaptive weights of cardinal evaluations
					Gintvar = vwt*gdar+(1-vwt)*guar;
					Ginthar = hwt*grar+(1-hwt)*glar;
					Gintvha = vwt*gdha+(1-vwt)*guha;
					Ginthha = hwt*grha+(1-hwt)*glha;
					//interpolated color differences
					vcd[indx] = sgn*(Gintvar-cfa[indx]);
					hcd[indx] = sgn*(Ginthar-cfa[indx]);
					vcdalt[indx] = sgn*(Gintvha-cfa[indx]);
					hcdalt[indx] = sgn*(Ginthha-cfa[indx]);
					
					//differences of interpolations in opposite directions
					dgintv[indx]=MIN(SQR(guha-gdha),SQR(guar-gdar));
					dginth[indx]=MIN(SQR(glha-grha),SQR(glar-grar));
				}
			//t2_vcdhcd += clock() - t1_vcdhcd;
			
			//t1_cdvar = clock();
			for (rr=4; rr<TS-4; rr++)
				//for (cc=4+(FC(rr,2)&1),indx=rr*TS+cc,c=FC(rr,cc); cc<TS-4; cc+=2,indx+=2) {
				for (cc=4,indx=rr*TS+cc; cc<TS-4; cc++,indx++) {
					c=FC(rr,cc);
					
					hcdvar = 3*(SQR(hcd[indx-2])+SQR(hcd[indx])+SQR(hcd[indx+2]))-SQR(hcd[indx-2]+hcd[indx]+hcd[indx+2]);
					hcdaltvar = 3*(SQR(hcdalt[indx-2])+SQR(hcdalt[indx])+SQR(hcdalt[indx+2]))-SQR(hcdalt[indx-2]+hcdalt[indx]+hcdalt[indx+2]);
					vcdvar = 3*(SQR(vcd[indx-v2])+SQR(vcd[indx])+SQR(vcd[indx+v2]))-SQR(vcd[indx-v2]+vcd[indx]+vcd[indx+v2]);
					vcdaltvar = 3*(SQR(vcdalt[indx-v2])+SQR(vcdalt[indx])+SQR(vcdalt[indx+v2]))-SQR(vcdalt[indx-v2]+vcdalt[indx]+vcdalt[indx+v2]);
					//choose the smallest variance; this yields a smoother interpolation
					if (hcdaltvar<hcdvar) hcd[indx]=hcdalt[indx];
					if (vcdaltvar<vcdvar) vcd[indx]=vcdalt[indx];
					
					//bound the interpolation in regions of high saturation
					if (c&1) {
						Ginth = -hcd[indx]+cfa[indx];//R or B
						Gintv = -vcd[indx]+cfa[indx];//B or R
						if (hcd[indx] < (0.33*(Ginth+cfa[indx]))) hcd[indx]=-ULIM(Ginth,cfa[indx-1],cfa[indx+1])+cfa[indx];
						if (vcd[indx] < (0.33*(Gintv+cfa[indx]))) vcd[indx]=-ULIM(Gintv,cfa[indx-v1],cfa[indx+v1])+cfa[indx];
					} else {
						Ginth = hcd[indx]+cfa[indx];
						Gintv = vcd[indx]+cfa[indx];
						if (hcd[indx] < (-0.33*(Ginth+cfa[indx]))) hcd[indx]=ULIM(Ginth,cfa[indx-1],cfa[indx+1])-cfa[indx];
						if (vcd[indx] < (-0.33*(Gintv+cfa[indx]))) vcd[indx]=ULIM(Gintv,cfa[indx-v1],cfa[indx+v1])-cfa[indx];
					}
					vcdsq[indx] = SQR(vcd[indx]);
					hcdsq[indx] = SQR(hcd[indx]);
					cddiffsq[indx] = SQR(vcd[indx]-hcd[indx]);
				}
			
			for (rr=6; rr<TS-6; rr++)
				for (cc=6+(FC(rr,2)&1),indx=rr*TS+cc; cc<TS-6; cc+=2,indx+=2) {
					
					//compute color difference variances in cardinal directions
					
					Dgrbvvaru = 4*(vcdsq[indx]+vcdsq[indx-v1]+vcdsq[indx-v2]+vcdsq[indx-v3])-SQR(vcd[indx]+vcd[indx-v1]+vcd[indx-v2]+vcd[indx-v3]);
					Dgrbvvard = 4*(vcdsq[indx]+vcdsq[indx+v1]+vcdsq[indx+v2]+vcdsq[indx+v3])-SQR(vcd[indx]+vcd[indx+v1]+vcd[indx+v2]+vcd[indx+v3]);
					Dgrbhvarl = 4*(hcdsq[indx]+hcdsq[indx-1]+hcdsq[indx-2]+hcdsq[indx-3])-SQR(hcd[indx]+hcd[indx-1]+hcd[indx-2]+hcd[indx-3]);
					Dgrbhvarr = 4*(hcdsq[indx]+hcdsq[indx+1]+hcdsq[indx+2]+hcdsq[indx+3])-SQR(hcd[indx]+hcd[indx+1]+hcd[indx+2]+hcd[indx+3]);
					
					hwt = dirwts[indx-1][1]/(dirwts[indx-1][1]+dirwts[indx+1][1]);
					vwt = dirwts[indx-v1][0]/(dirwts[indx+v1][0]+dirwts[indx-v1][0]);
					
					vcdvar = eps+vwt*Dgrbvvard+(1-vwt)*Dgrbvvaru;
					hcdvar = eps+hwt*Dgrbhvarr+(1-hwt)*Dgrbhvarl;
					
					//compute fluctuations in up/down and left/right interpolations of colors
					Dgrbvvaru = (dgintv[indx])+(dgintv[indx-v1])+(dgintv[indx-v2]);
					Dgrbvvard = (dgintv[indx])+(dgintv[indx+v1])+(dgintv[indx+v2]);
					Dgrbhvarl = (dginth[indx])+(dginth[indx-1])+(dginth[indx-2]);
					Dgrbhvarr = (dginth[indx])+(dginth[indx+1])+(dginth[indx+2]);
					
					vcdvar1 = eps+vwt*Dgrbvvard+(1-vwt)*Dgrbvvaru;
					hcdvar1 = eps+hwt*Dgrbhvarr+(1-hwt)*Dgrbhvarl;
					
					//determine adaptive weights for G interpolation
					varwt=hcdvar/(vcdvar+hcdvar);
					diffwt=hcdvar1/(vcdvar1+hcdvar1);
					
					//if both agree on interpolation direction, choose the one with strongest directional discrimination;
					//otherwise, choose the u/d and l/r difference fluctuation weights
					if ((0.5-varwt)*(0.5-diffwt)>0 && fabs(0.5-diffwt)<fabs(0.5-varwt)) {hvwt[indx]=varwt;} else {hvwt[indx]=diffwt;}
				}
			//t2_cdvar += clock() - t1_cdvar;
			
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// Nyquist test				 
			//t1_nyqtest = clock();
			
			for (rr=6; rr<TS-6; rr++)
				for (cc=6+(FC(rr,2)&1),indx=rr*TS+cc; cc<TS-6; cc+=2,indx+=2) {
					
					//nyquist texture test: ask if difference of vcd compared to hcd is larger or smaller than RGGB gradients
					nyqtest = (gaussodd[0]*cddiffsq[indx]+ \
							   gaussodd[1]*(cddiffsq[indx-m1]+cddiffsq[indx+p1]+ \
											cddiffsq[indx-p1]+cddiffsq[indx+m1])+ \
							   gaussodd[2]*(cddiffsq[indx-v2]+cddiffsq[indx-2]+ \
											cddiffsq[indx+2]+cddiffsq[indx+v2])+ \
							   gaussodd[3]*(cddiffsq[indx-m2]+cddiffsq[indx+p2]+ \
											cddiffsq[indx-p2]+cddiffsq[indx+m2]));
					
					nyqtest -= nyqthresh*(gaussgrad[0]*(delhsq[indx]+delvsq[indx])+ \
										  gaussgrad[1]*(delhsq[indx-v1]+delvsq[indx-v1]+delhsq[indx+1]+delvsq[indx+1]+ \
														delhsq[indx-1]+delvsq[indx-1]+delhsq[indx+v1]+delvsq[indx+v1])+ \
										  gaussgrad[2]*(delhsq[indx-m1]+delvsq[indx-m1]+delhsq[indx+p1]+delvsq[indx+p1]+ \
														delhsq[indx-p1]+delvsq[indx-p1]+delhsq[indx+m1]+delvsq[indx+m1])+ \
										  gaussgrad[3]*(delhsq[indx-v2]+delvsq[indx-v2]+delhsq[indx-2]+delvsq[indx-2]+ \
														delhsq[indx+2]+delvsq[indx+2]+delhsq[indx+v2]+delvsq[indx+v2])+ \
										  gaussgrad[4]*(delhsq[indx-2*TS-1]+delvsq[indx-2*TS-1]+delhsq[indx-2*TS+1]+delvsq[indx-2*TS+1]+ \
														delhsq[indx-TS-2]+delvsq[indx-TS-2]+delhsq[indx-TS+2]+delvsq[indx-TS+2]+ \
														delhsq[indx+TS-2]+delvsq[indx+TS-2]+delhsq[indx+TS+2]+delvsq[indx-TS+2]+ \
														delhsq[indx+2*TS-1]+delvsq[indx+2*TS-1]+delhsq[indx+2*TS+1]+delvsq[indx+2*TS+1])+ \
										  gaussgrad[5]*(delhsq[indx-m2]+delvsq[indx-m2]+delhsq[indx+p2]+delvsq[indx+p2]+ \
														delhsq[indx-p2]+delvsq[indx-p2]+delhsq[indx+m2]+delvsq[indx+m2]));
					
					
					if (nyqtest>0) {nyquist[indx]=1;}//nyquist=1 for nyquist region
				}
			
			for (rr=8; rr<TS-8; rr++)
				for (cc=8+(FC(rr,2)&1),indx=rr*TS+cc; cc<TS-8; cc+=2,indx+=2) {
					
					areawt=(nyquist[indx-v2]+nyquist[indx-m1]+nyquist[indx+p1]+ \
							nyquist[indx-2]+nyquist[indx]+nyquist[indx+2]+ \
							nyquist[indx-p1]+nyquist[indx+m1]+nyquist[indx+v2]);
					//if most of your neighbors are named Nyquist, it's likely that you're one too
					if (areawt>4) nyquist[indx]=1;
					//or not
					if (areawt<4) nyquist[indx]=0;
				}
			
			//t2_nyqtest += clock() - t1_nyqtest;
			// end of Nyquist test
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
						
			
			
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// in areas of Nyquist texture, do area interpolation
			//t1_areainterp = clock();
			for (rr=8; rr<TS-8; rr++)
				for (cc=8+(FC(rr,2)&1),indx=rr*TS+cc; cc<TS-8; cc+=2,indx+=2) {
					
					if (nyquist[indx]) {
						// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						// area interpolation
						
						sumh=sumv=sumsqh=sumsqv=areawt=0;
						for (i=-6; i<7; i+=2)
							for (j=-6; j<7; j+=2) {
								indx1=(rr+i)*TS+cc+j;
								if (nyquist[indx1]) {
									sumh += cfa[indx1]-0.5*(cfa[indx1-1]+cfa[indx1+1]);
									sumv += cfa[indx1]-0.5*(cfa[indx1-v1]+cfa[indx1+v1]);
									sumsqh += 0.5*(SQR(cfa[indx1]-cfa[indx1-1])+SQR(cfa[indx1]-cfa[indx1+1]));
									sumsqv += 0.5*(SQR(cfa[indx1]-cfa[indx1-v1])+SQR(cfa[indx1]-cfa[indx1+v1]));
									areawt +=1;
								}
							}
						
						//horizontal and vertical color differences, and adaptive weight
						hcdvar=eps+MAX(0, areawt*sumsqh-sumh*sumh);
						vcdvar=eps+MAX(0, areawt*sumsqv-sumv*sumv);
						hvwt[indx]=hcdvar/(vcdvar+hcdvar);
						
						// end of area interpolation
						// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						
					}
				}
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//t2_areainterp += clock() - t1_areainterp;
			
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			//populate G at R/B sites
			for (rr=8; rr<TS-8; rr++)
				for (cc=8+(FC(rr,2)&1),indx=rr*TS+cc; cc<TS-8; cc+=2,indx+=2) {
					
					//first ask if one gets more directional discrimination from nearby B/R sites
					hvwtalt = 0.25*(hvwt[indx-m1]+hvwt[indx+p1]+hvwt[indx-p1]+hvwt[indx+m1]);
					vo=fabs(0.5-hvwt[indx]);
					ve=fabs(0.5-hvwtalt);
					if (vo<ve) {hvwt[indx]=hvwtalt;}//a better result was obtained from the neighbors
					Dgrb[indx][0] = (hcd[indx]*(1-hvwt[indx]) + vcd[indx]*hvwt[indx]);//evaluate color differences
					rgb[indx][1] = cfa[indx] + Dgrb[indx][0];//evaluate G (finally!)
					
					//local curvature in G (preparation for nyquist refinement step)
					if (nyquist[indx]) {
						Dgrbh2[indx] = SQR(rgb[indx][1] - 0.5*(rgb[indx-1][1]+rgb[indx+1][1]));
						Dgrbv2[indx] = SQR(rgb[indx][1] - 0.5*(rgb[indx-v1][1]+rgb[indx+v1][1]));
					}
				}
			
			//end of standard interpolation
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// refine Nyquist areas using G curvatures
			
			for (rr=8; rr<TS-8; rr++)
				for (cc=8+(FC(rr,2)&1),indx=rr*TS+cc; cc<TS-8; cc+=2,indx+=2) {
					
					if (nyquist[indx]) {
						//local averages (over Nyquist pixels only) of G curvature squared 
						gvarh = eps + (gquinc[0]*Dgrbh2[indx]+ \
									   gquinc[1]*(Dgrbh2[indx-m1]+Dgrbh2[indx+p1]+Dgrbh2[indx-p1]+Dgrbh2[indx+m1])+ \
									   gquinc[2]*(Dgrbh2[indx-v2]+Dgrbh2[indx-2]+Dgrbh2[indx+2]+Dgrbh2[indx+v2])+ \
									   gquinc[3]*(Dgrbh2[indx-m2]+Dgrbh2[indx+p2]+Dgrbh2[indx-p2]+Dgrbh2[indx+m2]));
						gvarv = eps + (gquinc[0]*Dgrbv2[indx]+ \
									   gquinc[1]*(Dgrbv2[indx-m1]+Dgrbv2[indx+p1]+Dgrbv2[indx-p1]+Dgrbv2[indx+m1])+ \
									   gquinc[2]*(Dgrbv2[indx-v2]+Dgrbv2[indx-2]+Dgrbv2[indx+2]+Dgrbv2[indx+v2])+ \
									   gquinc[3]*(Dgrbv2[indx-m2]+Dgrbv2[indx+p2]+Dgrbv2[indx-p2]+Dgrbv2[indx+m2]));
						//use the results as weights for refined G interpolation
						Dgrb[indx][0] = (hcd[indx]*gvarv + vcd[indx]*gvarh)/(gvarv+gvarh);
						rgb[indx][1] = cfa[indx] + Dgrb[indx][0];
					}
				}
			
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						
			//t1_diag = clock();
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// diagonal interpolation correction
			for (rr=8; rr<TS-8; rr++)
				for (cc=8+(FC(rr,2)&1),indx=rr*TS+cc; cc<TS-8; cc+=2,indx+=2) {
					
					
					//evaluate diagonal gradients based on CFA values
					//first, diagonal gradients between G1 and G2
					gradp=0.25*((fabs(cfa[indx-2*TS+1]-cfa[indx-v1])+2*fabs(cfa[indx-v1]-cfa[indx-1])+fabs(cfa[indx-1]-cfa[indx+TS-2]))+ \
								(fabs(cfa[indx-TS+2]-cfa[indx+1])+2*fabs(cfa[indx+1]-cfa[indx+v1])+fabs(cfa[indx+v1]-cfa[indx+2*TS-1])));
					gradm=0.25*((fabs(cfa[indx-2*TS-1]-cfa[indx-v1])+2*fabs(cfa[indx-v1]-cfa[indx+1])+fabs(cfa[indx+1]-cfa[indx+TS+2]))+ \
								(fabs(cfa[indx-TS-2]-cfa[indx-1])+2*fabs(cfa[indx-1]-cfa[indx+v1])+fabs(cfa[indx+v1]-cfa[indx+2*TS+1])));		
					
					//diagonal gradients within RGGB planes
					gradp += eps + (gauss1[0]*delp[indx]+gauss1[1]*(delp[indx-v1]+delp[indx-1]+delp[indx+1]+delp[indx+v1])+ \
									gauss1[2]*(delp[indx-m1]+delp[indx+p1]+delp[indx-p1]+delp[indx+m1]));
					gradm += eps + (gauss1[0]*delm[indx]+gauss1[1]*(delm[indx-v1]+delm[indx-1]+delm[indx+1]+delm[indx+v1])+ \
									gauss1[2]*(delm[indx-m1]+delm[indx+p1]+delm[indx-p1]+delm[indx+m1]));
					
					gradpm = fabs((gradp - gradm)/(gradp + gradm));
					
					//hor/vert gradients within RGGB planes
					gradv = eps + (gauss1[0]*delv[indx]+gauss1[1]*(delv[indx-v1]+delv[indx-1]+delv[indx+1]+delv[indx+v1])+ \
								   gauss1[2]*(delv[indx-m1]+delv[indx+p1]+delv[indx-p1]+delv[indx+m1]));
					gradh = eps + (gauss1[0]*delh[indx]+gauss1[1]*(delh[indx-v1]+delh[indx-1]+delh[indx+1]+delh[indx+v1])+ \
								   gauss1[2]*(delh[indx-m1]+delh[indx+p1]+delh[indx-p1]+delh[indx+m1]));
					gradhv = fabs((gradv - gradh)/(gradv + gradh));
					
					
					// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					
					
					if (gradpm-gradhv<pmthresh) continue;
					
					//otherwise do diagonal interpolation correction
					
					for (dir=0; dir<5; dir++){
						indx1 = indx + nbr[dir];
						if (rbint[indx1]) continue;
						
						rbvarp = eps + (gausseven[0]*(Dgrbpsq1[indx1-v1]+Dgrbpsq1[indx1-1]+Dgrbpsq1[indx1+1]+Dgrbpsq1[indx1+v1]) + \
										gausseven[1]*(Dgrbpsq1[indx1-v2-1]+Dgrbpsq1[indx1-v2+1]+Dgrbpsq1[indx1-2-v1]+Dgrbpsq1[indx1+2-v1]+ \
													  Dgrbpsq1[indx1-2+v1]+Dgrbpsq1[indx1+2+v1]+Dgrbpsq1[indx1+v2-1]+Dgrbpsq1[indx1+v2+1]));
						rbvarp -=  SQR( (gausseven[0]*(Dgrbp1[indx1-v1]+Dgrbp1[indx1-1]+Dgrbp1[indx1+1]+Dgrbp1[indx1+v1]) + \
										 gausseven[1]*(Dgrbp1[indx1-v2-1]+Dgrbp1[indx1-v2+1]+Dgrbp1[indx1-2-v1]+Dgrbp1[indx1+2-v1]+ \
													   Dgrbp1[indx1-2+v1]+Dgrbp1[indx1+2+v1]+Dgrbp1[indx1+v2-1]+Dgrbp1[indx1+v2+1])));
						rbvarm = eps + (gausseven[0]*(Dgrbmsq1[indx1-v1]+Dgrbmsq1[indx1-1]+Dgrbmsq1[indx1+1]+Dgrbmsq1[indx1+v1]) + \
										gausseven[1]*(Dgrbmsq1[indx1-v2-1]+Dgrbmsq1[indx1-v2+1]+Dgrbmsq1[indx1-2-v1]+Dgrbmsq1[indx1+2-v1]+ \
													  Dgrbmsq1[indx1-2+v1]+Dgrbmsq1[indx1+2+v1]+Dgrbmsq1[indx1+v2-1]+Dgrbmsq1[indx1+v2+1]));
						rbvarm -=  SQR( (gausseven[0]*(Dgrbm1[indx1-v1]+Dgrbm1[indx1-1]+Dgrbm1[indx1+1]+Dgrbm1[indx1+v1]) + \
										 gausseven[1]*(Dgrbm1[indx1-v2-1]+Dgrbm1[indx1-v2+1]+Dgrbm1[indx1-2-v1]+Dgrbm1[indx1+2-v1]+ \
													   Dgrbm1[indx1-2+v1]+Dgrbm1[indx1+2+v1]+Dgrbm1[indx1+v2-1]+Dgrbm1[indx1+v2+1])));
						
						//diagonal color ratios
						crp=(cfa[indx1-p1]+cfa[indx1+p1])/(eps+cfa[indx1]+0.5*(cfa[indx1-p2]+cfa[indx1+p2]));
						crm=(cfa[indx1-m1]+cfa[indx1+m1])/(eps+cfa[indx1]+0.5*(cfa[indx1-m2]+cfa[indx1+m2]));
						
						//assign B/R at R/B sites
						if (fabs(1-crp)<arthresh) {rbp=cfa[indx1]*crp;}
						else {rbp=0.5*(cfa[indx1]+cfa[indx1-p1]+cfa[indx1+p1])-0.25*(cfa[indx1-p2]+cfa[indx1+p2]);}
						if (fabs(1-crm)<arthresh) {rbm=cfa[indx1]*crm;}
						else {rbm=0.5*(cfa[indx1]+cfa[indx1-m1]+cfa[indx1+p1])-0.25*(cfa[indx1-m2]+cfa[indx1+m2]);}
						
						
						if (2*rbp < cfa[indx1]) {rbp=ULIM(rbp,cfa[indx1-p1],cfa[indx1+p1]);}
						if (2*rbm < cfa[indx1]) {rbm=ULIM(rbm,cfa[indx1-m1],cfa[indx1+m1]);}
						
						rbint[indx1] = 0.5*(cfa[indx1] + (rbp*rbvarm+rbm*rbvarp)/(rbvarp+rbvarm));//this is R+B, interpolated
						//rbint[indx1] = 0.5*(cfa[indx1] + (rbp*(1-pmwt[indx1])+rbm*pmwt[indx1]));//this is R+B, interpolated
						
					}//end of populating neighbor B/R values
					//now interpolate G vertically/horizontally using R+B values
					//unfortunately, since G interpolation cannot be done diagonally this may lead to color shifts
					//color ratios for G interpolation
					
					cru = cfa[indx-v1]*2/(eps+rbint[indx]+rbint[indx-v2]);
					crd = cfa[indx+v1]*2/(eps+rbint[indx]+rbint[indx+v2]);
					crl = cfa[indx-1]*2/(eps+rbint[indx]+rbint[indx-2]);
					crr = cfa[indx+1]*2/(eps+rbint[indx]+rbint[indx+2]);
					
					//interpolated G via adaptive ratios or Hamilton-Adams in each cardinal direction
					if (fabs(1-cru)<arthresh) {gu=rbint[indx]*cru;} 
					else {gu=cfa[indx-v1]+0.5*(rbint[indx]-rbint[indx-v2]);}
					if (fabs(1-crd)<arthresh) {gd=rbint[indx]*crd;} 
					else {gd=cfa[indx+v1]+0.5*(rbint[indx]-rbint[indx+v2]);}
					if (fabs(1-crl)<arthresh) {gl=rbint[indx]*crl;} 
					else {gl=cfa[indx-1]+0.5*(rbint[indx]-rbint[indx-2]);}
					if (fabs(1-crr)<arthresh) {gr=rbint[indx]*crr;} 
					else {gr=cfa[indx+1]+0.5*(rbint[indx]-rbint[indx+2]);}
					
					//interpolated G via adaptive weights of cardinal evaluations
					Gintv = (dirwts[indx-v1][0]*gd+dirwts[indx+v1][0]*gu)/(dirwts[indx+v1][0]+dirwts[indx-v1][0]);
					Ginth = (dirwts[indx-1][1]*gr+dirwts[indx+1][1]*gl)/(dirwts[indx-1][1]+dirwts[indx+1][1]);
					
					if (rbint[indx]-2*Ginth > (0.33*(2*Ginth+rbint[indx]))) Ginth=ULIM(Ginth,cfa[indx-1],cfa[indx+1]);
					if (rbint[indx]-2*Gintv > (0.33*(2*Gintv+rbint[indx]))) Gintv=ULIM(Gintv,cfa[indx-v1],cfa[indx+v1]);
					
					rgb[indx][1] = Ginth*(1-hvwt[indx]) + Gintv*hvwt[indx];
				}
			//end of diagonal interpolation correction
			//t2_diag += clock() - t1_diag;
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			//t1_chroma = clock();
			//fancy chrominance interpolation
			//(ey,ex) is location of R site
			for (rr=13-ey; rr<TS-12; rr+=2)
				for (cc=13-ex,indx=rr*TS+cc; cc<TS-12; cc+=2,indx+=2) {//B coset
					Dgrb[indx][1]=Dgrb[indx][0];//split out G-B from G-R
					Dgrb[indx][0]=0;
				}
			for (rr=12; rr<TS-12; rr++)
				for (cc=12+(FC(rr,2)&1),indx=rr*TS+cc,c=1-FC(rr,cc)/2; cc<TS-12; cc+=2,indx+=2) {
					f[0]=1.0/(eps+fabs(Dgrb[indx-m1][c]-Dgrb[indx+m1][c])+fabs(Dgrb[indx-m1][c]-Dgrb[indx-m3][c])+fabs(Dgrb[indx+m1][c]-Dgrb[indx-m3][c]));
					f[1]=1.0/(eps+fabs(Dgrb[indx+p1][c]-Dgrb[indx-p1][c])+fabs(Dgrb[indx+p1][c]-Dgrb[indx+p3][c])+fabs(Dgrb[indx-p1][c]-Dgrb[indx+p3][c]));
					f[2]=1.0/(eps+fabs(Dgrb[indx-p1][c]-Dgrb[indx+p1][c])+fabs(Dgrb[indx-p1][c]-Dgrb[indx+m3][c])+fabs(Dgrb[indx+p1][c]-Dgrb[indx-p3][c]));
					f[3]=1.0/(eps+fabs(Dgrb[indx+m1][c]-Dgrb[indx-m1][c])+fabs(Dgrb[indx+m1][c]-Dgrb[indx-p3][c])+fabs(Dgrb[indx-m1][c]-Dgrb[indx+m3][c]));
					
					g[0]=Dgrb[indx-m1][c];
					g[1]=Dgrb[indx+p1][c];
					g[2]=Dgrb[indx-p1][c];
					g[3]=Dgrb[indx+m1][c];
					Dgrb[indx][c]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3])/(f[0]+f[1]+f[2]+f[3]);
				}
			for (rr=12; rr<TS-12; rr++)
				for (cc=12+(FC(rr,1)&1),indx=rr*TS+cc,c=FC(rr,cc+1)/2; cc<TS-12; cc+=2,indx+=2)
					for(c=0;c<2;c++){
						
						Dgrb[indx][c]=((hvwt[indx-v1])*Dgrb[indx-v1][c]+(1-hvwt[indx+1])*Dgrb[indx+1][c]+(1-hvwt[indx-1])*Dgrb[indx-1][c]+(hvwt[indx+v1])*Dgrb[indx+v1][c])/ \
						((hvwt[indx-v1])+(1-hvwt[indx+1])+(1-hvwt[indx-1])+(hvwt[indx+v1]));
						
					}
			for(rr=12; rr<TS-12; rr++)
				for(cc=12,indx=rr*TS+cc; cc<TS-12; cc++,indx++){
					rgb[indx][0]=(rgb[indx][1]-Dgrb[indx][0]);
					rgb[indx][2]=(rgb[indx][1]-Dgrb[indx][1]);
				}
			//t2_chroma += clock() - t1_chroma;
			
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			// copy smoothed results back to image matrix
			for (rr=16; rr < rr1-16; rr++)
				for (row=rr+top, cc=16; cc < cc1-16; cc++) {
					col = cc + left;
					
					indx=rr*TS+cc;
					
					red[row][col] = CLIP((int)(65535.0f*rgb[indx][0] + 0.5f));
					green[row][col] = CLIP((int)(65535.0f*rgb[indx][1] + 0.5f));
					blue[row][col] = CLIP((int)(65535.0f*rgb[indx][2] + 0.5f));
					
				}
			//end of main loop
			
			if(plistener) plistener->setProgress(fabs((float)top/height));
		}	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	
	// clean up
	free(buffer);
	
	// done

#undef TS
	
}
