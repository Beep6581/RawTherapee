#include <cmath>
#include "rt_math.h"

#include "EdgePreservingDecomposition.h"

/* Solves A x = b by the conjugate gradient method, where instead of feeding it the matrix A you feed it a function which
calculates A x where x is some vector. Stops when rms residual < RMSResidual or when maximum iterates is reached.
Stops at n iterates if MaximumIterates = 0 since that many iterates gives exact solution. Applicable to symmetric positive
definite problems only, which is what unconstrained smooth optimization pretty much always is.
Parameter pass can be passed through, containing whatever info you like it to contain (matrix info?).
Takes less memory with OkToModify_b = true, and Preconditioner = NULL. */
float *SparseConjugateGradient(void Ax(float *Product, float *x, void *Pass), float *b, unsigned int n, bool OkToModify_b, 
	float *x, float RMSResidual, void *Pass, unsigned int MaximumIterates, void Preconditioner(float *Product, float *x, void *Pass)){
	unsigned int iterate, i;

	//Start r and x.
	float *r = new float[n];
	if(x == NULL){
		x = new float[n];
		memset(x, 0, sizeof(float)*n);		//Zero initial guess if x == NULL.
		memcpy(r, b, sizeof(float)*n);
	}else{
		Ax(r, x, Pass);
		for(i = 0; i != n; i++) r[i] = b[i] - r[i];		//r = b - A x.
	}

	//s is preconditionment of r. Without, direct to r.
	float *s = r, rs = 0.0f;
	if(Preconditioner != NULL){
		s = new float[n];
		Preconditioner(s, r, Pass);
	}
	for(i = 0; i != n; i++) rs += r[i]*s[i];

	//Search direction d.
	float *d = new float[n];
	memcpy(d, s, sizeof(float)*n);

	//Store calculations of Ax in this.
	float *ax = b;
	if(!OkToModify_b) ax = new float[n];

	//Start iterating!
	if(MaximumIterates == 0) MaximumIterates = n;
	for(iterate = 0; iterate != MaximumIterates; iterate++){
		//Get step size alpha, store ax while at it.
		float ab = 0.0f;
		Ax(ax, d, Pass);
		for(i = 0; i != n; i++) ab += d[i]*ax[i];

		if(ab == 0.0f) break;	//So unlikely. It means perfectly converged or singular, stop either way.
		ab = rs/ab;

		//Update x and r with this step size.
		float rms = 0.0;
		for(i = 0; i != n; i++){
			x[i] += ab*d[i];
			r[i] -= ab*ax[i];	//"Fast recursive formula", use explicit r = b - Ax occasionally?
			rms += r[i]*r[i];
		}
		rms = sqrtf(rms/n);
//printf("%f\n", rms);
		//Quit? This probably isn't the best stopping condition, but ok.
		if(rms < RMSResidual) break;

		if(Preconditioner != NULL) Preconditioner(s, r, Pass);

		//Get beta.
		ab = rs;
		rs = 0.0f;
		for(i = 0; i != n; i++) rs += r[i]*s[i];
		ab = rs/ab;

		//Update search direction p.
		for(i = 0; i != n; i++) d[i] = s[i] + ab*d[i];
	}
	if(iterate == MaximumIterates)
		if(iterate != n && RMSResidual != 0.0f)
			printf("Warning: MaximumIterates (%u) reached in SparseConjugateGradient.\n", MaximumIterates);

	if(ax != b) delete[] ax;
	if(s != r) delete[] s;
	delete[] r;
	delete[] d;
	return x;
}


MultiDiagonalSymmetricMatrix::MultiDiagonalSymmetricMatrix(unsigned int Dimension, unsigned int NumberOfDiagonalsInLowerTriangle){
	n = Dimension;
	m = NumberOfDiagonalsInLowerTriangle;
	IncompleteCholeskyFactorization = NULL;

	Diagonals = new float *[m];
	StartRows = new unsigned int [m];
	memset(Diagonals, 0, sizeof(float *)*m);
	memset(StartRows, 0, sizeof(unsigned int)*m);
}

MultiDiagonalSymmetricMatrix::~MultiDiagonalSymmetricMatrix(){
	for(unsigned int i = 0; i != m; i++) delete[] Diagonals[i];
	delete[] Diagonals;
	delete[] StartRows;
}

bool MultiDiagonalSymmetricMatrix::CreateDiagonal(unsigned int index, unsigned int StartRow){
	if(index >= m){
		printf("Error in MultiDiagonalSymmetricMatrix::CreateDiagonal: invalid index.\n");
		return false;
	}
	if(index > 0)
		if(StartRow <= StartRows[index - 1]){
			printf("Error in MultiDiagonalSymmetricMatrix::CreateDiagonal: each StartRow must exceed the previous.\n");
			return false;
		}

	delete[] Diagonals[index];
	Diagonals[index] = new float[DiagonalLength(StartRow)];
	if(Diagonals[index] == NULL){
		printf("Error in MultiDiagonalSymmetricMatrix::CreateDiagonal: memory allocation failed. Out of memory?\n");
		return false;
	}

	StartRows[index] = StartRow;
	memset(Diagonals[index], 0, sizeof(float)*DiagonalLength(StartRow));
	return true;
}

int MultiDiagonalSymmetricMatrix::FindIndex(unsigned int StartRow){
	//There's GOT to be a better way to do this. "Bidirectional map?"
	for(unsigned int i = 0; i != m; i++)
		if(StartRows[i] == StartRow)
			return i;
	return -1;
}

bool MultiDiagonalSymmetricMatrix::LazySetEntry(float value, unsigned int row, unsigned int column){
	//On the strict upper triangle? Swap, this is ok due to symmetry.
	int i, sr;
	if(column > row)
		i = column,
		column = row,
		row = i;
	if(row >= n) return false;
	sr = row - column;

	//Locate the relevant diagonal.
	i = FindIndex(sr);
	if(i < 0) return false;

	Diagonals[i][column] = value;
	return true;
}

void MultiDiagonalSymmetricMatrix::VectorProduct(float *Product, float *x){
	//Initialize to zero.
	memset(Product, 0, n*sizeof(float));

	//Loop over the stored diagonals.
	for(unsigned int i = 0; i != m; i++){
		unsigned int sr = StartRows[i];
		float *a = Diagonals[i];	//One fewer dereference.
		unsigned int j, l = DiagonalLength(sr);

		if(sr == 0)
			for(j = 0; j != l; j++)
				Product[j] += a[j]*x[j];		//Separate, fairly simple treatment for the main diagonal.
		else
			for(j = 0; j != l; j++)
				Product[j + sr] += a[j]*x[j],	//Contribution from lower...
				Product[j] += a[j]*x[j + sr];	//...and upper triangle.
	}
}

bool MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization(unsigned int MaxFillAbove){
	if(m == 1){
		printf("Error in MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization: just one diagonal? Can you divide?\n");
		return false;
	}
	if(StartRows[0] != 0){
		printf("Error in MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization: main diagonal required to exist for this math.\n");
		return false;
	}

	//How many diagonals in the decomposition?
	MaxFillAbove++;	//Conceptually, now "fill" includes an existing diagonal. Simpler in the math that follows.
	unsigned int i, j, mic;
	for(mic = i = 1; i != m; i++)
		mic += rtengine::min(StartRows[i] - StartRows[i - 1], MaxFillAbove);	//Guarunteed positive since StartRows must be created in increasing order.

	//Initialize the decomposition - setup memory, start rows, etc.
	MultiDiagonalSymmetricMatrix *ic = new MultiDiagonalSymmetricMatrix(n, mic);
	ic->CreateDiagonal(0, 0);	//There's always a main diagonal in this type of decomposition.
	for(mic = i = 1; i != m; i++){
		//Set j to the number of diagonals to be created corresponding to a diagonal on this source matrix...
		j = rtengine::min(StartRows[i] - StartRows[i - 1], MaxFillAbove);

		//...and create those diagonals. I want to take a moment to tell you about how much I love minimalistic loops: very much.
		while(j-- != 0)
			if(!ic->CreateDiagonal(mic++, StartRows[i] - j)){
				//Beware of out of memory, possible for large, sparse problems if you ask for too much fill.
				printf("Error in MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization: Out of memory. Ask for less fill?\n");
				delete ic;
				return false;
			}
	}


	//It's all initialized? Uhkay. Do the actual math then.
	int sss, ss, s;
	unsigned int k, MaxStartRow = StartRows[m - 1];	//Handy number.
	float **l = ic->Diagonals;
	float  *d = ic->Diagonals[0];		//Describes D in LDLt.

	//Loop over the columns.
	for(j = 0; j != n; j++){
		//Calculate d for this column.
		d[j] = Diagonals[0][j];

		//This is a loop over k from 1 to j, inclusive. We'll cover that by looping over the index of the diagonals (s), and get k from it.
		//The first diagonal is d (k = 0), so skip that and have s start at 1. Cover all available s but stop if k exceeds j.
		for(s = 1; s != ic->m; s++){
			k = ic->StartRows[s];
			if(k > j) break;
			d[j] -= l[s][j - k]*l[s][j - k]*d[j - k];
		}

		if(d[j] == 0.0f){
			printf("Error in MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization: division by zero. Matrix not decomposable.\n");
			delete ic;
			return false;
		}
		float id = 1.0f/d[j];

		//Now, calculate l from top down along this column.
		for(s = 1; s != ic->m; s++){
			i = ic->StartRows[s];		//Start row for this entry.
			if(j >= ic->n - i) break;	//Possible values of j are limited.

			//Quicker access for an element of l.
			float *lij = &l[s][j];
			sss = FindIndex(i);	//Find element in same spot in the source matrix. It might be a zero.
			*lij = sss < 0 ? 0.0f : Diagonals[sss][j];

			//Similar to the loop involving d, convoluted by the fact that two l are involved.
			for(ss = 1; ss != ic->m; ss++){
				k = ic->StartRows[ss];
				if(k > j) break;
				if(i + k > MaxStartRow) break;	//Quick exit once k to big.

				int sss = ic->FindIndex(i + k);
				if(sss < 0) continue;			//Asked for diagonal nonexistant. But there may be something later, so don't break.

				/* Let's think about the factors in the term below for a moment.
				j varies from 0 to n - 1, so j - k is bounded inclusive by 0 and j - 1. So d[j - k] is always in the matrix.

				l[sss] and l[ss] are diagonals with corresponding start rows i + k and k.
				For l[sss][j - k] to exist, we must have j - k < n - (i + k) -> j < n - i, which was checked outside this loop and true at this point.
				For l[ ss][j - k] to exist, we must have j - k < n - k -> j < n, which is true straight from definition.

				So, no additional checks, all is good and within bounds at this point.*/
				*lij -= l[sss][j - k]*l[ss][j - k]*d[j - k];
			}

			*lij *= id;
		}
	}

	IncompleteCholeskyFactorization = ic;
	return true;
}

void MultiDiagonalSymmetricMatrix::KillIncompleteCholeskyFactorization(void){
	delete IncompleteCholeskyFactorization;
}

void MultiDiagonalSymmetricMatrix::CholeskyBackSolve(float *x, float *b){
	//We want to solve L D Lt x = b where D is a diagonal matrix described by Diagonals[0] and L is a unit lower triagular matrix described by the rest of the diagonals.
	//Let D Lt x = y. Then, first solve L y = b.
	float *y = new float[n];
	float **d = IncompleteCholeskyFactorization->Diagonals;
	unsigned int *s = IncompleteCholeskyFactorization->StartRows;
	unsigned int M = IncompleteCholeskyFactorization->m, N = IncompleteCholeskyFactorization->n;
	unsigned int i, j;
	for(j = 0; j != N; j++){
		y[j] = b[j];

		for(i = 1; i != M; i++){	//Start at 1 because zero is D.
			int c = (int)j - (int)s[i];
			if(c < 0) break;		//Due to ordering of StartRows, no further contributions.
			y[j] -= d[i][c]*y[c];
		}
	}

	//Now, solve x from D Lt x = y -> Lt x = D^-1 y
	while(j-- != 0){
		x[j] = y[j]/d[0][j];

		for(i = 1; i != M; i++){
			if(j + s[i] >= N) break;
			x[j] -= d[i][j]*x[j + s[i]];
		}
	}

	delete[] y;
}




EdgePreservingDecomposition::EdgePreservingDecomposition(unsigned int width, unsigned int height){
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
		printf("Error in EdgePreservingDecomposition construction: out of memory.\n");
	}else{
		a0    = A->Diagonals[0];
		a_1   = A->Diagonals[1];
		a_w1  = A->Diagonals[2];
		a_w   = A->Diagonals[3];
		a_w_1 = A->Diagonals[4];
	}
}

EdgePreservingDecomposition::~EdgePreservingDecomposition(){
	delete A;
}

float *EdgePreservingDecomposition::CreateBlur(float *Source, float Scale, float EdgeStopping, unsigned int Iterates, float *Blur, bool UseBlurForEdgeStop){
	if(Blur == NULL)
		UseBlurForEdgeStop = false,	//Use source if there's no supplied Blur.
		Blur = new float[n];
	if(Scale == 0.0f){
		memcpy(Blur, Source, n*sizeof(float));
		return Blur;
	}

	//Create the edge stopping function a, rotationally symmetric and just one instead of (ax, ay). Maybe don't need Blur yet, so use its memory.
	float *a, *g;
	if(UseBlurForEdgeStop) a = new float[n], g = Blur;
	else a = Blur, g = Source;

	unsigned int x, y, i;
	unsigned int w1 = w - 1, h1 = h - 1;
	float eps = 0.02f;
	for(y = 0; y != h1; y++){
		float *rg = &g[w*y];
		for(x = 0; x != w1; x++){
			//Estimate the central difference gradient in the center of a four pixel square. (gx, gy) is actually 2*gradient.
			float gx = (rg[x + 1] - rg[x]) + (rg[x + w + 1] - rg[x + w]);
			float gy = (rg[x + w] - rg[x]) + (rg[x + w + 1] - rg[x + 1]);

			//Apply power to the magnitude of the gradient to get the edge stopping function.
			a[x + w*y] = Scale*powf(0.5f*sqrtf(gx*gx + gy*gy + eps*eps), -EdgeStopping);
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
	if(!UseBlurForEdgeStop) memcpy(Blur, Source, n*sizeof(float));
	SparseConjugateGradient(A->PassThroughVectorProduct, Source, n, false, Blur, 0.0f, (void *)A, Iterates, A->PassThroughCholeskyBackSolve);
	A->KillIncompleteCholeskyFactorization();
	return Blur;
}

float *EdgePreservingDecomposition::CreateIteratedBlur(float *Source, float Scale, float EdgeStopping, unsigned int Iterates, unsigned int Reweightings, float *Blur){
	//Simpler outcome?
	if(Reweightings == 0) return CreateBlur(Source, Scale, EdgeStopping, Iterates, Blur);

	//Create a blur here, initialize.
	if(Blur == NULL) Blur = new float[n];
	memcpy(Blur, Source, n*sizeof(float));

	//Iteratively improve the blur.
	Reweightings++;
	for(unsigned int i = 0; i != Reweightings; i++)
		CreateBlur(Source, Scale, EdgeStopping, Iterates, Blur, true);

	return Blur;
}

float *EdgePreservingDecomposition::CompressDynamicRange(float *Source, float Scale, float EdgeStopping, float CompressionExponent, float DetailBoost, unsigned int Iterates, unsigned int Reweightings, float *Compressed){
	//Small number intended to prevent division by zero. This is different from the eps in CreateBlur.
	const float eps = 0.0001f;

	//We're working with luminance, which does better logarithmic.
	unsigned int i;
	for(i = 0; i != n; i++)
		Source[i] = logf(Source[i] + eps);

	//Blur. Also setup memory for Compressed (we can just use u since each element of u is used in one calculation).
	float *u = CreateIteratedBlur(Source, Scale, EdgeStopping, Iterates, Reweightings);
	if(Compressed == NULL) Compressed = u;

	//Apply compression, detail boost, unlogging. Compression is done on the logged data and detail boost on unlogged.
	for(i = 0; i != n; i++){
		float ce = expf(Source[i] + u[i]*(CompressionExponent - 1.0f)) - eps;
		float ue = expf(u[i]) - eps;
		Source[i] = expf(Source[i]) - eps;
		Compressed[i] = ce + DetailBoost*(Source[i] - ue);
	}

	if(Compressed != u) delete[] u;
	return Compressed;
}

