
////////////////////////////////////////////////////////////////
//
//		//exposure correction before interpolation
//
//  code dated: December 27, 2010
//
//	Expo_before.cc is free software: you can redistribute it and/or modify
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

//		     Jacques Desmis <jdesmis@gmail.com>
//   	     use fast-demo(provisional) from Emil Martinec
//			 inspired from work Guillermo Luijk and Manuel LLorens(Perfectraw)
// I use OMP
  // This function uses parameters:
  //      exposure (lineal): 2^(-8..0..8): currently 0.5 +3
  //      preserve (log)   : 0..8 : currently 0.1 1

//modi : 31/12/2010

#define LIM(x,min,max) MAX(min,MIN(x,max))
#define CLIPF(x) LIM(x,0.0,65535.0)
void RawImageSource::exp_bef(float expos, float preser) {
	double dt,dT2;
	clock_t t1, t2,t3,t4,t5;
	float Yp, exposure2, K, EV;
	//float LUT[65536];
float *LUT = new float[65536];	
	int i;
	int row,col;
	int width=W, height=H;

	// I use with Dcraw FDD interpolate from Luis Sanz , very fast and good, one can change for another : here with Rawtherpee == fastdemo() from Emil Martinec
		float** imgd; 
		imgd = allocArray< float >(W,H);//with memcpy : faster than for (...)
		
		for (int i=0; i<H; i++) {
				memcpy (imgd[i], rawData[i], W*sizeof(**imgd));}//save configuration but perhaps instable...
	 
	 fast_demo (0,0,W,H);//from Emil
// calculate CIE luminosity
float *YY;
YY = (float *)calloc(width*height,sizeof *YY);// for CIE luminosity
#pragma omp parallel default(shared)  
{
#pragma omp for  
	 for(int row=0;row<height;row++)
		for(int col=0;col<width;col++)
			{int i=row*width+col;
			YY[i]=CLIPF(0.299*red[row][col]+0.587*green[row][col]+0.114*blue[row][col]); // CIE luminosity
				}
}	

		for (int i=0; i<H; i++) {
				memcpy (rawData[i], imgd[i], W*sizeof(**imgd));}//restore config
			
freeArray<float>(imgd, H);//free memory imgd

	//exposure correction inspired from G.Luijk
 if(preser==0.0){	// protect highlights 
#pragma omp parallel for  shared(expos)
	 for(int row=0;row<height;row++)
		for(int col=0;col<width;col++)
			{rawData[row][col]=CLIPF(rawData[row][col]*expos);}

  }else{
    // Exposure correction with highlight preservation
    if(expos>1){
      K=65535/expos*exp(-preser*log((double) 2));
      for(int j=0;j<=65535;j++) LUT[(int)j]=CLIPF(((65535-K*expos)/(65535-K)*(j-65535)+65535)/j);

#pragma omp parallel for  shared(expos)
	 for(int row=0;row<height;row++)
		for(int col=0;col<width;col++){
			if(YY[row*width+col]<K){
				rawData[row][col]=CLIPF(rawData[row][col]*expos);}
			else{
				float exposure2=LUT[(int)YY[row*width+col]];
				rawData[row][col]=CLIPF(rawData[row][col]*exposure2);}}
				
			}
	else{
      float EV=log(expos)/log(2.0);                              // Convert exp. linear to EV
      float K=65535.0*exp(-preser*log((double) 2));
      for(int j=0;j<=65535;j++) LUT[(int)j]=CLIPF(exp(EV*(65535.0-j)/(65535.0-K)*log((double) 2)));
#pragma omp parallel for  shared(expos)	  
	 for(int row=0;row<height;row++)
		for(int col=0;col<width;col++){
			if(YY[row*width+col]<K){
				rawData[row][col]=CLIPF(rawData[row][col]*expos);}
				
			else{
				float exposure2=LUT[(int)YY[row*width+col]];
			rawData[row][col]=CLIPF(rawData[row][col]*exposure2);}}
				
    }	
}
free(YY);
delete [] LUT;
}