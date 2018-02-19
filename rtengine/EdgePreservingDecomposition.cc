#include <cmath>
#include "rt_math.h"
#include "EdgePreservingDecomposition.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "sleef.c"
#include "opthelper.h"

#define DIAGONALS 5
#define DIAGONALSP1 6

/* Solves A x = b by the conjugate gradient method, where instead of feeding it the matrix A you feed it a function which
calculates A x where x is some vector. Stops when rms residual < RMSResidual or when maximum iterates is reached.
Stops at n iterates if MaximumIterates = 0 since that many iterates gives exact solution. Applicable to symmetric positive
definite problems only, which is what unconstrained smooth optimization pretty much always is.
Parameter pass can be passed through, containing whatever info you like it to contain (matrix info?).
Takes less memory with OkToModify_b = true, and Preconditioner = nullptr. */
float *SparseConjugateGradient(void Ax(float *Product, float *x, void *Pass), float *b, int n, bool OkToModify_b,
                               float *x, float RMSResidual, void *Pass, int MaximumIterates, void Preconditioner(float *Product, float *x, void *Pass))
{
    int iterate;

    float* buffer = (float*)malloc(2 * n * sizeof(float) + 128);
    float *r = (buffer + 16);

    //Start r and x.
    if(x == nullptr) {
        x = new float[n];

        memset(x, 0, sizeof(float)*n);      //Zero initial guess if x == nullptr.
        memcpy(r, b, sizeof(float)*n);
    } else {
        Ax(r, x, Pass);
#ifdef _OPENMP
        #pragma omp parallel for           // removed schedule(dynamic,10)
#endif

        for(int ii = 0; ii < n; ii++) {
            r[ii] = b[ii] - r[ii];    //r = b - A x.
        }
    }

    //s is preconditionment of r. Without, direct to r.
    float *s = r, rs = 0.0f;

    if(Preconditioner != nullptr) {
        s = new float[n];

        Preconditioner(s, r, Pass);
    }

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:rs)  // removed schedule(dynamic,10)
#endif

    for(int ii = 0; ii < n; ii++) {
        rs += r[ii] * s[ii];
    }

    //Search direction d.
    float *d = (buffer + n + 32);

    memcpy(d, s, sizeof(float)*n);

    //Store calculations of Ax in this.
    float *ax = b;

    if(!OkToModify_b) {
        ax = new float[n];
    }

    //Start iterating!
    if(MaximumIterates == 0) {
        MaximumIterates = n;
    }

    for(iterate = 0; iterate < MaximumIterates; iterate++) {
        //Get step size alpha, store ax while at it.
        float ab = 0.0f;
        Ax(ax, d, Pass);
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:ab)
#endif

        for(int ii = 0; ii < n; ii++) {
            ab += d[ii] * ax[ii];
        }

        if(ab == 0.0f) {
            break;    //So unlikely. It means perfectly converged or singular, stop either way.
        }

        ab = rs / ab;

        //Update x and r with this step size.
        float rms = 0.0;
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:rms)
#endif

        for(int ii = 0; ii < n; ii++) {
            x[ii] += ab * d[ii];
            r[ii] -= ab * ax[ii]; //"Fast recursive formula", use explicit r = b - Ax occasionally?
            rms += r[ii] * r[ii];
        }

        rms = sqrtf(rms / n);

        //Quit? This probably isn't the best stopping condition, but ok.
        if(rms < RMSResidual) {
            break;
        }

        if(Preconditioner != nullptr) {
            Preconditioner(s, r, Pass);
        }

        //Get beta.
        ab = rs;
        rs = 0.0f;

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            float c = 0.0f;
#ifdef _OPENMP
            #pragma omp for reduction(+:rs)                            // Summation with error correction
#endif

            for(int ii = 0; ii < n; ii++) {
                float temp = r[ii] * s[ii];
                float t = rs + temp;

                if( fabsf(rs) >= fabsf(temp) ) {
                    c += ((rs - t) + temp);
                } else {
                    c += ((temp - t) + rs);
                }

                rs = t;
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            rs += c;
        }

        ab = rs / ab;

        //Update search direction p.
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for(int ii = 0; ii < n; ii++) {
            d[ii] = s[ii] + ab * d[ii];
        }


    }

    if(iterate == MaximumIterates)
        if(iterate != n && RMSResidual != 0.0f) {
            printf("Warning: MaximumIterates (%u) reached in SparseConjugateGradient.\n", MaximumIterates);
        }

    if(ax != b) {
        delete[] ax;
    }

    if(s != r) {
        delete[] s;
    }

    free(buffer);
    return x;
}

MultiDiagonalSymmetricMatrix::MultiDiagonalSymmetricMatrix(int Dimension, int NumberOfDiagonalsInLowerTriangle) : buffer(nullptr), DiagBuffer(nullptr)
{
    n = Dimension;
    m = NumberOfDiagonalsInLowerTriangle;
    IncompleteCholeskyFactorization = nullptr;

    Diagonals = new float *[m];
    StartRows = new int [m + 1];
    memset(Diagonals, 0, sizeof(float *)*m);
    memset(StartRows, 0, sizeof(int) * (m + 1));
    StartRows[m] = n + 1;
}

MultiDiagonalSymmetricMatrix::~MultiDiagonalSymmetricMatrix()
{
    if(DiagBuffer != nullptr) {
        free(buffer);
    } else
        for(int i = 0; i < m; i++) {
            delete[] Diagonals[i];
        }

    delete[] Diagonals;
    delete[] StartRows;
}

bool MultiDiagonalSymmetricMatrix::CreateDiagonal(int index, int StartRow)
{
    // Changed memory allocation for diagonals to avoid L1 conflict misses
    // Falls back to original version if big block could not be allocated
    int padding = 4096 - ((n * m * sizeof(float)) % 4096);

    if(index == 0) {
        buffer = (char*)calloc( (n + padding) * m * sizeof(float) + (m + 16) * 64 + 63, 1);

        if(buffer == nullptr)
            // no big memory block available => try to allocate smaller blocks
        {
            DiagBuffer = nullptr;
        } else {
            DiagBuffer = (float*)( ( uintptr_t(buffer) + uintptr_t(63)) / 64 * 64);
        }
    }

    if(index >= m) {
        printf("Error in MultiDiagonalSymmetricMatrix::CreateDiagonal: invalid index.\n");
        return false;
    }

    if(index > 0)
        if(StartRow <= StartRows[index - 1]) {
            printf("Error in MultiDiagonalSymmetricMatrix::CreateDiagonal: each StartRow must exceed the previous.\n");
            return false;
        }

    if(DiagBuffer != nullptr) {
        Diagonals[index] = (DiagBuffer + (index * (n + padding)) + ((index + 16) * 16));
    } else {
        Diagonals[index] = new float[DiagonalLength(StartRow)];

        if(Diagonals[index] == nullptr) {
            printf("Error in MultiDiagonalSymmetricMatrix::CreateDiagonal: memory allocation failed. Out of memory?\n");
            return false;
        }

        memset(Diagonals[index], 0, sizeof(float)*DiagonalLength(StartRow));
    }

    StartRows[index] = StartRow;
    return true;
}

inline int MultiDiagonalSymmetricMatrix::FindIndex(int StartRow)
{
    //There's GOT to be a better way to do this. "Bidirectional map?"
    // Issue 1895 : Changed start of loop from zero to one
    // m is small (5 or 6)
    for(int i = 1; i < m; i++)
        if(StartRows[i] == StartRow) {
            return i;
        }

    return -1;
}

bool MultiDiagonalSymmetricMatrix::LazySetEntry(float value, int row, int column)
{
    //On the strict upper triangle? Swap, this is ok due to symmetry.
    int i, sr;

    if(column > row)
        i = column,
        column = row,
        row = i;

    if(row >= n) {
        return false;
    }

    sr = row - column;

    //Locate the relevant diagonal.
    i = FindIndex(sr);

    if(i < 0) {
        return false;
    }

    Diagonals[i][column] = value;
    return true;
}

void MultiDiagonalSymmetricMatrix::VectorProduct(float* RESTRICT Product, float* RESTRICT x)
{

    int srm = StartRows[m - 1];
    int lm = DiagonalLength(srm);
#ifdef _OPENMP
#ifdef __SSE2__
    const int chunkSize = (lm - srm) / (omp_get_num_procs() * 32);
#else
    const int chunkSize = (lm - srm) / (omp_get_num_procs() * 8);
#endif
    #pragma omp parallel
#endif
    {
        // First fill the big part in the middle
        // This can be done without intermediate stores to memory and it can be parallelized too
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,chunkSize) nowait
#endif
#ifdef __SSE2__

        for(int j = srm; j < lm - 3; j += 4) {
            __m128 prodv = LVFU(Diagonals[0][j]) * LVFU(x[j]);

            for(int i = m - 1; i > 0; i--) {
                int s = StartRows[i];
                prodv += (LVFU(Diagonals[i][j - s]) * LVFU(x[j - s])) + (LVFU(Diagonals[i][j]) * LVFU(x[j + s]));
            }

            _mm_storeu_ps(&Product[j], prodv);
        }

#else

        for(int j = srm; j < lm; j++) {
            float prod = Diagonals[0][j] * x[j];

            for(int i = m - 1; i > 0; i--) {
                int s = StartRows[i];
                prod += (Diagonals[i][j - s] * x[j - s]) + (Diagonals[i][j] * x[j + s]);
            }

            Product[j] = prod;
        }

#endif
        #pragma omp single
        {
#ifdef __SSE2__

            for(int j = lm - ((lm - srm) % 4); j < lm; j++) {
                float prod = Diagonals[0][j] * x[j];

                for(int i = m - 1; i > 0; i--) {
                    int s = StartRows[i];
                    prod += (Diagonals[i][j - s] * x[j - s]) + (Diagonals[i][j] * x[j + s]);
                }

                Product[j] = prod;
            }

#endif

            // Fill remaining area
            // Loop over the stored diagonals.
            for(int i = 0; i < m; i++) {
                int sr = StartRows[i];
                float *a = Diagonals[i];    //One fewer dereference.
                int l = DiagonalLength(sr);

                if(sr == 0) {
                    for(int j = 0; j < srm; j++) {
                        Product[j] = a[j] * x[j];    //Separate, fairly simple treatment for the main diagonal.
                    }

                    for(int j = lm; j < l; j++) {
                        Product[j] = a[j] * x[j];    //Separate, fairly simple treatment for the main diagonal.
                    }
                } else {
// Split the loop in 3 parts, so now the big one in the middle can be parallelized without race conditions
                    // updates 0 to sr - 1. Because sr is small (in the range of image-width) no benefit by omp
                    for(int j = 0; j < sr; j++) {
                        Product[j] += a[j] * x[j + sr];     //Contribution from upper triangle
                    }

                    // Updates sr to l - 1. Because sr is small and l is big, this loop is parallelized
                    for(int j = sr; j < srm; j++) {
                        Product[j] += a[j - sr] * x[j - sr] + a[j] * x[j + sr]; //Contribution from lower and upper triangle
                    }

                    for(int j = lm; j < l; j++) {
                        Product[j] += a[j - sr] * x[j - sr] + a[j] * x[j + sr]; //Contribution from lower and upper triangle
                    }

                    // Updates l to l + sr - 1. Because sr is small (in the range of image-width) no benefit by omp
                    for(int j = l; j < l + sr; j++) {
                        Product[j] += a[j - sr] * x[j - sr]; //Contribution from lower triangle
                    }
                }
            }
        }
    }
#ifdef _OPENMP
    static_cast<void>(chunkSize); // to silence cppcheck warning
#endif
}

bool MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization(int MaxFillAbove)
{
    if(m == 1) {
        printf("Error in MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization: just one diagonal? Can you divide?\n");
        return false;
    }

    if(StartRows[0] != 0) {
        printf("Error in MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization: main diagonal required to exist for this math.\n");
        return false;
    }

    //How many diagonals in the decomposition?
    MaxFillAbove++; //Conceptually, now "fill" includes an existing diagonal. Simpler in the math that follows.
    int j, mic, fp;
    mic = 1;
    fp = 1;

    for(int ii = 1; ii < m; ii++) {
        fp = rtengine::min(StartRows[ii] - StartRows[ii - 1], MaxFillAbove);    //Guarunteed positive since StartRows must be created in increasing order.
        mic = mic + fp;
    }

    //Initialize the decomposition - setup memory, start rows, etc.

    MultiDiagonalSymmetricMatrix *ic = new MultiDiagonalSymmetricMatrix(n, mic);
    if(!ic->CreateDiagonal(0, 0)) { //There's always a main diagonal in this type of decomposition.
        delete ic;
        return false;
    }
    mic = 1;

    for(int ii = 1; ii < m; ii++) {
        //Set j to the number of diagonals to be created corresponding to a diagonal on this source matrix...
        j = rtengine::min(StartRows[ii] - StartRows[ii - 1], MaxFillAbove);

        //...and create those diagonals. I want to take a moment to tell you about how much I love minimalistic loops: very much.
        while(j-- != 0)
            if(!ic->CreateDiagonal(mic++, StartRows[ii] - j)) {
                //Beware of out of memory, possible for large, sparse problems if you ask for too much fill.
                printf("Error in MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization: Out of memory. Ask for less fill?\n");
                delete ic;
                return false;
            }
    }

    //It's all initialized? Uhkay. Do the actual math then.
    int sss;
   // int MaxStartRow = StartRows[m - 1];  //Handy number.
    float **l = ic->Diagonals;
    float  *d = ic->Diagonals[0];       //Describes D in LDLt.
    int icm = ic->m;
    int icn = ic->n;
    int* RESTRICT icStartRows = ic->StartRows;

    //Loop over the columns.

    // create array for quicker access to ic->StartRows
    struct s_diagmap {
        int sss;
        int ss;
        int k;
    };


    // Pass one: count number of needed entries
    int entrycount = 0;

    for(int i = 1; i < icm; i++) {
        for(int j = 1; j < icm; j++) {
            if(ic->FindIndex( icStartRows[i] + icStartRows[j]) > 0) {
                entrycount ++;
            }
        }
    }

    // now we can create the array
    struct s_diagmap* RESTRICT DiagMap = new s_diagmap[entrycount];
    // we also need the maxvalues
    int entrynumber = 0;
    int index;
    int* RESTRICT MaxIndizes = new int[icm];

    for(int i = 1; i < icm; i++) {
        for(int j = 1; j < icm; j++) {
            index = ic->FindIndex( icStartRows[i] + icStartRows[j]);

            if(index > 0) {
                DiagMap[entrynumber].ss = j;
                DiagMap[entrynumber].sss = index;
                DiagMap[entrynumber].k = icStartRows[j];
                entrynumber ++;
            }
        }

        MaxIndizes[i] = entrynumber - 1;
    }

    int* RESTRICT findmap = new int[icm];

    for(int j = 0; j < icm; j++) {
        findmap[j] = FindIndex( icStartRows[j]);
    }

    for(j = 0; j < n; j++) {
        //Calculate d for this column.
        d[j] = Diagonals[0][j];

        //This is a loop over k from 1 to j, inclusive. We'll cover that by looping over the index of the diagonals (s), and get k from it.
        //The first diagonal is d (k = 0), so skip that and have s start at 1. Cover all available s but stop if k exceeds j.
        int s = 1;
        int k = icStartRows[s];

        while(k <= j) {
            d[j] -= l[s][j - k] * l[s][j - k] * d[j - k];
            s++;
            k = icStartRows[s];
        }

        if(UNLIKELY(d[j] == 0.0f)) {
            printf("Error in MultiDiagonalSymmetricMatrix::CreateIncompleteCholeskyFactorization: division by zero. Matrix not decomposable.\n");
            delete ic;
            delete[] DiagMap;
            delete[] MaxIndizes;
            delete[] findmap;
            return false;
        }

        float id = 1.0f / d[j];
        //Now, calculate l from top down along this column.

        int mapindex = 0;
        int jMax = icn - j;

        for(s = 1; s < icm; s++) {
            if(icStartRows[s] >= jMax) {
                break;    //Possible values of j are limited
            }

            float temp = 0.0f;

            while(mapindex <= MaxIndizes[s] && ( k = DiagMap[mapindex].k) <= j) {
                temp -= l[DiagMap[mapindex].sss][j - k] * l[DiagMap[mapindex].ss][j - k] * d[j - k];
                mapindex ++;
            }

            sss = findmap[s];
            l[s][j] = id * (sss < 0 ? temp : (Diagonals[sss][j] + temp));
        }
    }

    delete[] DiagMap;
    delete[] MaxIndizes;
    delete[] findmap;
    IncompleteCholeskyFactorization = ic;
    return true;
}

void MultiDiagonalSymmetricMatrix::KillIncompleteCholeskyFactorization()
{
    delete IncompleteCholeskyFactorization;
}

void MultiDiagonalSymmetricMatrix::CholeskyBackSolve(float* RESTRICT x, float* RESTRICT b)
{
    //We want to solve L D Lt x = b where D is a diagonal matrix described by Diagonals[0] and L is a unit lower triagular matrix described by the rest of the diagonals.
    //Let D Lt x = y. Then, first solve L y = b.
    float* RESTRICT  *d = IncompleteCholeskyFactorization->Diagonals;
    int* RESTRICT s = IncompleteCholeskyFactorization->StartRows;
    int M = IncompleteCholeskyFactorization->m, N = IncompleteCholeskyFactorization->n;
    int i, j;

    if(M != DIAGONALSP1) {                  // can happen in theory
        for(j = 0; j < N; j++) {
            float sub = b[j];                   // using local var to reduce memory writes, gave a big speedup
            i = 1;
            int c = j - s[i];

            while(c >= 0) {
                sub -= d[i][c] * x[c];
                i++;
                c = j - s[i];
            }

            x[j] = sub;                         // only one memory-write per j
        }
    } else {                               // that's the case almost every time
        for(j = 0; j <= s[M - 1]; j++) {
            float sub = b[j];                   // using local var to reduce memory writes, gave a big speedup
            i = 1;
            int c = j - s[1];

            while(c >= 0) {
                sub -= d[i][c] * x[c];
                i++;
                c = j - s[i];
            }

            x[j] = sub;                         // only one memory-write per j
        }

        for(j = s[M - 1] + 1; j < N; j++) {
            float sub = b[j];                   // using local var to reduce memory writes, gave a big speedup

            for(int i = DIAGONALSP1 - 1; i > 0; i--) { // using a constant upperbound allows the compiler to unroll this loop (gives a good speedup)
                sub -= d[i][j - s[i]] * x[j - s[i]];
            }

            x[j] = sub;                         // only one memory-write per j
        }
    }

    //Now, solve x from D Lt x = y -> Lt x = D^-1 y
// Took this one out of the while, so it can be parallelized now, which speeds up, because division is expensive
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(j = 0; j < N; j++) {
        x[j] = x[j] / d[0][j];
    }

    if(M != DIAGONALSP1) {                  // can happen in theory
        while(j-- > 0) {
            float sub = x[j];                   // using local var to reduce memory writes, gave a big speedup
            i = 1;
            int c = j + s[1];

            while(c < N) {
                sub -= d[i][j] * x[c];
                i++;
                c = j + s[i];
            }

            x[j] = sub;                         // only one memory-write per j
        }
    } else {                                // that's the case almost every time
        for(j = N - 1; j >= (N - 1) - s[M - 1]; j--) {
            float sub = x[j];                   // using local var to reduce memory writes, gave a big speedup
            i = 1;
            int c = j + s[1];

            while(c < N) {
                sub -= d[i][j] * x[j + s[i]];
                i++;
                c = j + s[i];
            }

            x[j] = sub;                         // only one memory-write per j
        }

        for(j = (N - 2) - s[M - 1]; j >= 0; j--) {
            float sub = x[j];                   // using local var to reduce memory writes, gave a big speedup

            for(int i = DIAGONALSP1 - 1; i > 0; i--) { // using a constant upperbound allows the compiler to unroll this loop (gives a good speedup)
                sub -= d[i][j] * x[j + s[i]];
            }

            x[j] = sub;                         // only one memory-write per j
        }
    }
}

EdgePreservingDecomposition::EdgePreservingDecomposition(int width, int height) : a0(nullptr) , a_1(nullptr), a_w(nullptr), a_w_1(nullptr), a_w1(nullptr)
{
    w = width;
    h = height;
    n = w * h;

    //Initialize the matrix just once at construction.
    A = new MultiDiagonalSymmetricMatrix(n, DIAGONALS);

    if(!(
                A->CreateDiagonal(0, 0) &&
                A->CreateDiagonal(1, 1) &&
                A->CreateDiagonal(2, w - 1) &&
                A->CreateDiagonal(3, w) &&
                A->CreateDiagonal(4, w + 1))) {
        delete A;
        A = nullptr;
        printf("Error in EdgePreservingDecomposition construction: out of memory.\n");
    } else {
        a0    = A->Diagonals[0];
        a_1   = A->Diagonals[1];
        a_w1  = A->Diagonals[2];
        a_w   = A->Diagonals[3];
        a_w_1 = A->Diagonals[4];
    }
}

EdgePreservingDecomposition::~EdgePreservingDecomposition()
{
    delete A;
}

float *EdgePreservingDecomposition::CreateBlur(float *Source, float Scale, float EdgeStopping, int Iterates, float *Blur, bool UseBlurForEdgeStop)
{

    if(Blur == nullptr)
        UseBlurForEdgeStop = false, //Use source if there's no supplied Blur.
        Blur = new float[n];

    if(Scale == 0.0f) {
        memcpy(Blur, Source, n * sizeof(float));
        return Blur;
    }

    //Create the edge stopping function a, rotationally symmetric and just one instead of (ax, ay). Maybe don't need Blur yet, so use its memory.
    float* RESTRICT a;
    float* RESTRICT g;

    if(UseBlurForEdgeStop) {
        a = new float[n], g = Blur;
    } else {
        a = Blur, g = Source;
    }

    int w1 = w - 1, h1 = h - 1;
    const float sqreps = 0.0004f;                           // removed eps*eps from inner loop


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef __SSE2__
        int x;
        __m128 gxv, gyv;
        __m128 Scalev = _mm_set1_ps( Scale );
        __m128 sqrepsv = _mm_set1_ps( sqreps );
        __m128 EdgeStoppingv = _mm_set1_ps( -EdgeStopping );
        __m128 zd5v = _mm_set1_ps( 0.5f );
#endif
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int y = 0; y < h1; y++) {
            float *rg = &g[w * y];
#ifdef __SSE2__

            for(x = 0; x < w1 - 3; x += 4) {
                //Estimate the central difference gradient in the center of a four pixel square. (gx, gy) is actually 2*gradient.
                gxv = (LVFU(rg[x + 1]) -  LVFU(rg[x])) + (LVFU(rg[x + w + 1]) - LVFU(rg[x + w]));
                gyv = (LVFU(rg[x + w]) -  LVFU(rg[x])) + (LVFU(rg[x + w + 1]) - LVFU(rg[x + 1]));
                //Apply power to the magnitude of the gradient to get the edge stopping function.
                _mm_storeu_ps( &a[x + w * y], Scalev * pow_F((zd5v * _mm_sqrt_ps(gxv * gxv + gyv * gyv + sqrepsv)), EdgeStoppingv) );
            }

            for(; x < w1; x++) {
                //Estimate the central difference gradient in the center of a four pixel square. (gx, gy) is actually 2*gradient.
                float gx = (rg[x + 1] - rg[x]) + (rg[x + w + 1] - rg[x + w]);
                float gy = (rg[x + w] - rg[x]) + (rg[x + w + 1] - rg[x + 1]);
                //Apply power to the magnitude of the gradient to get the edge stopping function.
                a[x + w * y] = Scale * pow_F(0.5f * sqrtf(gx * gx + gy * gy + sqreps), -EdgeStopping);
            }

#else

            for(int x = 0; x < w1; x++) {
                //Estimate the central difference gradient in the center of a four pixel square. (gx, gy) is actually 2*gradient.
                float gx = (rg[x + 1] - rg[x]) + (rg[x + w + 1] - rg[x + w]);
                float gy = (rg[x + w] - rg[x]) + (rg[x + w + 1] - rg[x + 1]);

                //Apply power to the magnitude of the gradient to get the edge stopping function.
                a[x + w * y] = Scale * pow_F(0.5f * sqrtf(gx * gx + gy * gy + sqreps), -EdgeStopping);
            }

#endif
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


// checked for race condition here
// a0[] is read and write but addressed by i only
// a[] is read only
// a_w_1 is write only
// a_w is write only
// a_w1 is write only
// a_1 is write only
// So, there should be no race conditions

#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(int y = 0; y < h; y++) {
        int i = y * w;

        for(int x = 0; x < w; x++, i++) {
            float ac, a0temp;
            a0temp = 0.25f;

            //Remember, only fill the lower triangle. Memory for upper is never made. It's symmetric. Trust.
            if(x > 0 && y > 0) {
                ac = a[i - w - 1] / 6.0f;
                a_w_1[i - w - 1] -= 2.0f * ac;
                a_w[i - w] -= ac;
                a_1[i - 1] -= ac;
                a0temp += ac;
            }

            if(x < w1 && y > 0) {
                ac = a[i - w] / 6.0f;
                a_w[i - w] -= ac;
                a_w1[i - w + 1] -= 2.0f * ac;
                a0temp += ac;
            }

            if(x > 0 && y < h1) {
                ac = a[i - 1] / 6.0f;
                a_1[i - 1] -= ac;
                a0temp += ac;
            }

            if(x < w1 && y < h1) {
                a0temp += a[i] / 6.0f;
            }

            a0[i] = 4.0f * a0temp;
        }
    }

    if(UseBlurForEdgeStop) {
        delete[] a;
    }

    //Solve & return.
    bool success = A->CreateIncompleteCholeskyFactorization(1); //Fill-in of 1 seems to work really good. More doesn't really help and less hurts (slightly).

    if(!success) {
        fprintf(stderr, "Error: Tonemapping has failed.\n");
        memset(Blur, 0, sizeof(float)*n);  // On failure, set the blur to zero.  This is subsequently exponentiated in CompressDynamicRange.
        return Blur;
    }

    if(!UseBlurForEdgeStop) {
        memcpy(Blur, Source, n * sizeof(float));
    }

    SparseConjugateGradient(A->PassThroughVectorProduct, Source, n, false, Blur, 0.0f, (void *)A, Iterates, A->PassThroughCholeskyBackSolve);
    A->KillIncompleteCholeskyFactorization();
    return Blur;
}

float *EdgePreservingDecomposition::CreateIteratedBlur(float *Source, float Scale, float EdgeStopping, int Iterates, int Reweightings, float *Blur)
{
    //Simpler outcome?
    if(Reweightings == 0) {
        return CreateBlur(Source, Scale, EdgeStopping, Iterates, Blur);
    }

    //Create a blur here, initialize.
    if(Blur == nullptr) {
        Blur = new float[n];
    }

    memcpy(Blur, Source, n * sizeof(float));

    //Iteratively improve the blur.
    Reweightings++;

    for(int i = 0; i < Reweightings; i++) {
        CreateBlur(Source, Scale, EdgeStopping, Iterates, Blur, true);
    }

    return Blur;
}

void EdgePreservingDecomposition::CompressDynamicRange(float *Source, float Scale, float EdgeStopping, float CompressionExponent, float DetailBoost, int Iterates, int Reweightings)
{
    if(w < 300 && h < 300) { // set number of Reweightings to zero for small images (thumbnails). We could try to find a better solution here.
        Reweightings = 0;
    }

    //Small number intended to prevent division by zero. This is different from the eps in CreateBlur.
    const float eps = 0.0001f;

    //We're working with luminance, which does better logarithmic.
#ifdef __SSE2__
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        __m128 epsv = _mm_set1_ps( eps );
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int ii = 0; ii < n - 3; ii += 4) {
            _mm_storeu_ps( &Source[ii], xlogf(LVFU(Source[ii]) + epsv));
        }
    }

    for(int ii = n - (n % 4); ii < n; ii++) {
        Source[ii] = xlogf(Source[ii] + eps);
    }

#else
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(int ii = 0; ii < n; ii++) {
        Source[ii] = xlogf(Source[ii] + eps);
    }

#endif

    //Blur. Also setup memory for Compressed (we can just use u since each element of u is used in one calculation).
    float *u = CreateIteratedBlur(Source, Scale, EdgeStopping, Iterates, Reweightings);

    //Apply compression, detail boost, unlogging. Compression is done on the logged data and detail boost on unlogged.
    float temp;

    if(DetailBoost > 0.f) {
        float betemp = expf(-(2.f - DetailBoost + 0.694f)) - 1.f; //0.694 = log(2)
        temp = 1.2f * xlogf( -betemp);
    } else {
        temp = CompressionExponent - 1.0f;
    }

#ifdef __SSE2__
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        __m128 cev, uev, sourcev;
        __m128 epsv = _mm_set1_ps( eps );
        __m128 DetailBoostv = _mm_set1_ps( DetailBoost );
        __m128 tempv = _mm_set1_ps( temp );
#ifdef _OPENMP
        #pragma omp for
#endif

        for(int i = 0; i < n - 3; i += 4) {
            cev = xexpf(LVFU(Source[i]) + LVFU(u[i]) * (tempv)) - epsv;
            uev = xexpf(LVFU(u[i])) - epsv;
            sourcev = xexpf(LVFU(Source[i])) - epsv;
            _mm_storeu_ps( &Source[i], cev + DetailBoostv * (sourcev - uev) );
        }
    }

    for(int i = n - (n % 4); i < n; i++) {
        float ce = xexpf(Source[i] + u[i] * (temp)) - eps;
        float ue = xexpf(u[i]) - eps;
        Source[i] = xexpf(Source[i]) - eps;
        Source[i] = ce + DetailBoost * (Source[i] - ue);
    }

#else
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(int i = 0; i < n; i++) {
        float ce = xexpf(Source[i] + u[i] * (temp)) - eps;
        float ue = xexpf(u[i]) - eps;
        Source[i] = xexpf(Source[i]) - eps;
        Source[i] = ce + DetailBoost * (Source[i] - ue);
    }

#endif

    delete[] u;
}

