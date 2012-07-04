#include "EdgePreserveLab.h"
#include "boxblur.h"
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

//#define MAX(a,b) ((a)<(b)?(b):(a))
//#define MIN(a,b) ((a)>(b)?(b):(a))

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EdgePreserveLab::EdgePreserveLab(unsigned int width, unsigned int height){
	w = width;
	h = height;
	n = w*h;

	//Initialize the matrix just once at construction.
	A = new MultiDiagonalSymmetricMatrix(n, 5);
	if(!(
		A->CreateDiagonal(0, 0) &&
		A->CreateDiagonal(1, 1) &&
		A->CreateDiagonal(2, w - 1) &&
		A->CreateDiagonal(3, w) &&
		A->CreateDiagonal(4, w + 1))){
		delete A;
		A = NULL;
		printf("Error in EdgePreserveLab construction: out of memory.\n");
	}else{
		a0    = A->Diagonals[0];
		a_1   = A->Diagonals[1];
		a_w1  = A->Diagonals[2];
		a_w   = A->Diagonals[3];
		a_w_1 = A->Diagonals[4];
	}
}

EdgePreserveLab::~EdgePreserveLab(){
	delete A;
}

float *EdgePreserveLab::CreateBlur(float *Source, float LScale, float abScale, float EdgeStoppingLuma, float EdgeStoppingChroma, unsigned int Iterates, float *Blur, bool UseBlurForEdgeStop){
	if(Blur == NULL)
		UseBlurForEdgeStop = false,	//Use source if there's no supplied Blur.
		Blur = new float[3*n];
	if(LScale == 0.0f){
		memcpy(Blur, Source, 3*n*sizeof(float));
		return Blur;
	}

	//Create the edge stopping function a, rotationally symmetric and just one instead of (ax, ay). Maybe don't need Blur yet, so use its memory.
	float *a, *b, *g;
	if(UseBlurForEdgeStop) a = new float[n], g = Blur;
	else a = Blur, g = Source;
	//b = new float[n];
	
	unsigned int x, y, i;
	unsigned int w1 = w - 1, h1 = h - 1;
	float eps = 0.0001f;
	float scL = powf(100.0f,LScale);
	float scab = powf(200.0f,abScale);
	
	float * var = new float[w*h];
	rtengine::boxvar(g, var, 1, 1, w, h);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(y = 0; y < h1; y++){
		float *rg = &g[w*y];
		for(x = 0; x < w1; x++){
			//Estimate the central difference gradient in the center of a four pixel square. (gx, gy) is actually 2*gradient.
			/*float gx = (fabs((rg[x + 1] - rg[x]) + (rg[x + w + 1] - rg[x + w])));
			float gy = (fabs((rg[x + w] - rg[x]) + (rg[x + w + 1] - rg[x + 1])));
						
			//TODO: combine this with gx, gy if not needing separate quantities
			float hx =  (fabs((rg[x + 1 + n] - rg[x + n]) + (rg[x + w + 1 + n] - rg[x + w + n])) + \
						 fabs((rg[x + 1 + 2*n] - rg[x + 2*n]) + (rg[x + w + 1 + 2*n] - rg[x + w + 2*n]))); 
			float hy = (fabs((rg[x + w + n] - rg[x + n]) + (rg[x + w + 1 + n] - rg[x + 1 + n])) + \
						fabs((rg[x + w + 2*n] - rg[x + 2*n]) + (rg[x + w + 1 + 2*n] - rg[x + 1 + 2*n])));
			*/
			//float gradtot = (gx+gy+hx+hy);
			//gradhisto[MAX(0,MIN(32767,(int)gradtot))] ++;
			
			//Apply power to the magnitude of the gradient to get the edge stopping function.
			//a[x + w*y] = scL*expf(-100.0f*(gx + gy + hx + hy)/(EdgeStoppingLuma));
			//a[x + w*y] = scL*expf(-var[y*w+x]/SQR(0.02*EdgeStoppingLuma));///(0.1+rg[x]);
			a[x + w*y] = scL*expf(-50.0f*sqrt(var[y*w+x])/(EdgeStoppingLuma+eps));///(0.1+rg[x]);

			//b[x + w*y] = (scab)*expf(-20.0f*(gx + gy + Lave*(hx + hy))/(EdgeStoppingChroma));
			//b[x + w*y] = (scab)*expf(-400.0f*SQR(gx + gy + Lave*(hx + hy))/SQR(EdgeStoppingChroma));;

		}
	}

	/* Now setup the linear problem. I use the Maxima CAS, here's code for making an FEM formulation for the smoothness term:
		p(x, y) := (1 - x)*(1 - y);
		P(m, n) := A[m][n]*p(x, y) + A[m + 1][n]*p(1 - x, y) + A[m + 1][n + 1]*p(1 - x, 1 - y) + A[m][n + 1]*p(x, 1 - y);
		Integrate(f) := integrate(integrate(f, x, 0, 1), y, 0, 1);

		Integrate(diff(P(u, v), x)*diff(p(x, y), x) + diff(P(u, v), y)*diff(p(x, y), y));
		Integrate(diff(P(u - 1, v), x)*diff(p(1 - x, y), x) + diff(P(u - 1, v), y)*diff(p(1 - x, y), y));
		Integrate(diff(P(u - 1, v - 1), x)*diff(p(1 - x, 1 - y), x) + diff(P(u - 1, v - 1), y)*diff(p(1 - x, 1 - y), y));
		Integrate(diff(P(u, v - 1), x)*diff(p(x, 1 - y), x) + diff(P(u, v - 1), y)*diff(p(x, 1 - y), y));
	So yeah. Use the numeric results of that to fill the matrix A.*/
	memset(a_1, 0, A->DiagonalLength(1)*sizeof(float));
	memset(a_w1, 0, A->DiagonalLength(w - 1)*sizeof(float));
	memset(a_w, 0, A->DiagonalLength(w)*sizeof(float));
	memset(a_w_1, 0, A->DiagonalLength(w + 1)*sizeof(float));

//TODO: OMP here?
	for(i = y = 0; y != h; y++){
		for(x = 0; x != w; x++, i++){
			float ac;
			a0[i] = 1.0;

			//Remember, only fill the lower triangle. Memory for upper is never made. It's symmetric. Trust.
			if(x > 0 && y > 0)
				ac = a[i - w - 1]/6.0f,
				a_w_1[i - w - 1] -= 2.0f*ac, a_w[i - w] -= ac,
				a_1[i - 1]       -=      ac, a0[i] += 4.0f*ac;

			if(x < w1 && y > 0)
				ac = a[i - w]/6.0f,
				a_w[i - w] -= ac, a_w1[i - w + 1] -= 2.0f*ac,
				a0[i] += 4.0f*ac;

			if(x > 0 && y < h1)
				ac = a[i - 1]/6.0f,
				a_1[i - 1] -= ac, a0[i] += 4.0f*ac;

			if(x < w1 && y < h1)
				a0[i] += 4.0f*a[i]/6.0f;
		}
	}
	if(UseBlurForEdgeStop) delete[] a;
	

	//Solve & return.
	A->CreateIncompleteCholeskyFactorization(1);	//Fill-in of 1 seems to work really good. More doesn't really help and less hurts (slightly).
	if(!UseBlurForEdgeStop) memcpy(Blur, Source, 3*n*sizeof(float));
	// blur L channel
	SparseConjugateGradient(A->PassThroughVectorProduct, Source, n, false, Blur, 0.0f, (void *)A, Iterates, A->PassThroughCholeskyBackSolve);
	
	//reset A for ab channels
	/*memset(a_1, 0, A->DiagonalLength(1)*sizeof(float));
	memset(a_w1, 0, A->DiagonalLength(w - 1)*sizeof(float));
	memset(a_w, 0, A->DiagonalLength(w)*sizeof(float));
	memset(a_w_1, 0, A->DiagonalLength(w + 1)*sizeof(float));
	for(i = y = 0; y != h; y++){
		for(x = 0; x != w; x++, i++){
			float ac;
			a0[i] = 1.0;
			
			//Remember, only fill the lower triangle. Memory for upper is never made. It's symmetric. Trust.
			if(x > 0 && y > 0)
				ac = b[i - w - 1]/6.0f,
				a_w_1[i - w - 1] -= 2.0f*ac, a_w[i - w] -= ac,
				a_1[i - 1]       -=      ac, a0[i] += 4.0f*ac;
			
			if(x < w1 && y > 0)
				ac = b[i - w]/6.0f,
				a_w[i - w] -= ac, a_w1[i - w + 1] -= 2.0f*ac,
				a0[i] += 4.0f*ac;
			
			if(x > 0 && y < h1)
				ac = b[i - 1]/6.0f,
				a_1[i - 1] -= ac, a0[i] += 4.0f*ac;
			
			if(x < w1 && y < h1)
				a0[i] += 4.0f*b[i]/6.0f;
		}
	}*/
	/*if(UseBlurForEdgeStop)*/ //delete[] b;
	
	// blur ab channels
	//A->CreateIncompleteCholeskyFactorization(1);	//Fill-in of 1 seems to work really good. More doesn't really help and less hurts (slightly).
	//SparseConjugateGradient(A->PassThroughVectorProduct, Source+n, n, false, Blur+n, 0.0f, (void *)A, Iterates, A->PassThroughCholeskyBackSolve);
	//SparseConjugateGradient(A->PassThroughVectorProduct, Source+2*n, n, false, Blur+2*n, 0.0f, (void *)A, Iterates, A->PassThroughCholeskyBackSolve);
	
	A->KillIncompleteCholeskyFactorization();
	return Blur;
}

float *EdgePreserveLab::CreateIteratedBlur(float *Source, float LScale, float abScale, float EdgeStoppingLuma, float EdgeStoppingChroma, unsigned int Iterates, unsigned int Reweightings, float *Blur){
	//Simpler outcome?
	if(Reweightings == 0) return CreateBlur(Source, LScale, abScale, EdgeStoppingLuma, EdgeStoppingChroma, Iterates, Blur);

	//Create a blur here, initialize.
	if(Blur == NULL) Blur = new float[3*n];
	memcpy(Blur, Source, 3*n*sizeof(float));

	//Iteratively improve the blur.
	Reweightings++;
	for(unsigned int i = 0; i != Reweightings; i++)
		CreateBlur(Source, LScale, abScale, EdgeStoppingLuma, EdgeStoppingChroma, Iterates, Blur, true);

	return Blur;
}

float *EdgePreserveLab::CompressDynamicRange(float *Source, float LScale, float abScale, float EdgeStoppingLuma, float EdgeStoppingChroma, float CompressionExponent, float DetailBoost, unsigned int Iterates, unsigned int Reweightings, float *Compressed){
	//We're working with luminance, which does better logarithmic.
	unsigned int i;
	//for(i = 0; i != n; i++)
	//	Source[i] = logf(Source[i] + 0.0001f);
	
	//Blur. Also setup memory for Compressed (we can just use u since each element of u is used in one calculation).
	float *u = CreateIteratedBlur(Source, LScale, abScale, EdgeStoppingLuma, EdgeStoppingChroma, Iterates, Reweightings);
	if(Compressed == NULL) Compressed = u;
	
	//Apply compression, detail boost, unlogging. Compression is done on the logged data and detail boost on unlogged.
	for(i = 0; i != n; i++){
		//float ce = expf(Source[i] + u[i]*(CompressionExponent - 1.0f)) - 0.0001f;
		//float ue = expf(u[i]) - 0.0001f;
		//Source[i] = expf(Source[i]) - 0.0001f;
		//Compressed[i] = ce + DetailBoost*(Source[i] - ue);
		Compressed[i] = u[i];//ue;//for testing, to display blur
	}
	
	if(Compressed != u) delete[] u;
	return Compressed;
}




