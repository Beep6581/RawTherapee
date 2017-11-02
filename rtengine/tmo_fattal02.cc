/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Ported from LuminanceHDR by Alberto Griggio <alberto.griggio@gmail.com>
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

/**
 * @file tmo_fattal02.cpp
 * @brief TMO: Gradient Domain High Dynamic Range Compression
 *
 * Implementation of Gradient Domain High Dynamic Range Compression
 * by Raanan Fattal, Dani Lischinski, Michael Werman.
 *
 * @author Grzegorz Krawczyk, <krawczyk@mpi-sb.mpg.de>
 *
 *
 * This file is a part of LuminanceHDR package, based on pfstmo.
 * ----------------------------------------------------------------------
 * Copyright (C) 2003,2004 Grzegorz Krawczyk
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * ----------------------------------------------------------------------
 *
 * $Id: tmo_fattal02.cpp,v 1.3 2008/11/04 23:43:08 rafm Exp $
 */


#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <limits>

#include <math.h>
#include <assert.h>
#include <fftw3.h>

#include "array2D.h"
#include "improcfun.h"
#include "settings.h"

namespace rtengine {

extern const Settings *settings;

using namespace std;

namespace {

class Array2Df: public array2D<float> {
    typedef array2D<float> Super;
public:
    Array2Df(): Super() {}
    Array2Df(int w, int h): Super(w, h) {}

    float &operator()(int w, int h)
    {
        return (*this)[h][w];
    }

    const float &operator()(int w, int h) const
    {
        return (*this)[h][w];
    }

    float &operator()(int i)
    {
        return static_cast<float *>(*this)[i];
    }

    const float &operator()(int i) const
    {
        return const_cast<Array2Df &>(*this).operator()(i);
    }

    int getRows() const
    {
        return const_cast<Array2Df &>(*this).height();
    }

    int getCols() const
    {
        return const_cast<Array2Df &>(*this).width();
    }

    float *data()
    {
        return static_cast<float *>(*this);
    }

    const float *data() const
    {
        return const_cast<Array2Df &>(*this).data();
    }
};


void downSample(const Array2Df& A, Array2Df& B)
{
    const int width = B.getCols();
    const int height = B.getRows();

    // Note, I've uncommented all omp directives. They are all ok but are
    // applied to too small problems and in total don't lead to noticable
    // speed improvements. The main issue is the pde solver and in case of the
    // fft solver uses optimised threaded fftw routines.
    //#pragma omp parallel for
    for ( int y=0 ; y<height ; y++ )
    {
        for ( int x=0 ; x<width ; x++ )
        {
            float p = A(2*x,2*y);
            p += A(2*x+1,2*y);
            p += A(2*x,2*y+1);
            p += A(2*x+1,2*y+1);
            B(x,y) = p * 0.25f; // p / 4.0f;
        }
    }
}

void gaussianBlur(const Array2Df& I, Array2Df& L)
{
    const int width = I.getCols();
    const int height = I.getRows();

    Array2Df T(width,height);

    //--- X blur
    //#pragma omp parallel for shared(I, T)
    for ( int y=0 ; y<height ; y++ )
    {
        for ( int x=1 ; x<width-1 ; x++ )
        {
            float t = 2.f * I(x,y);
            t += I(x-1,y);
            t += I(x+1,y);
            T(x,y) = t * 0.25f; // t / 4.f;
        }
        T(0,y) = ( 3.f * I(0,y)+ I(1,y) ) * 0.25f; // / 4.f;
        T(width-1,y) = ( 3.f * I(width-1,y) + I(width-2,y) ) * 0.25f; // / 4.f;
    }

    //--- Y blur
    //#pragma omp parallel for shared(T, L)
    for ( int x=0 ; x<width ; x++ )
    {
        for ( int y=1 ; y<height-1 ; y++ )
        {
            float t = 2.f * T(x,y);
            t += T(x,y-1);
            t += T(x,y+1);
            L(x,y) = t * 0.25f; // t/4.0f;
        }
        L(x,0) = ( 3.f * T(x,0) + T(x,1) ) * 0.25f; // / 4.0f;
        L(x,height-1) = ( 3.f * T(x,height-1) + T(x,height-2) ) * 0.25f; // / 4.0f;
    }
}

void createGaussianPyramids( Array2Df* H, Array2Df** pyramids, int nlevels)
{
  int width = H->getCols();
  int height = H->getRows();
  const int size = width*height;

  pyramids[0] = new Array2Df(width,height);
//#pragma omp parallel for shared(pyramids, H)
  for( int i=0 ; i<size ; i++ )
    (*pyramids[0])(i) = (*H)(i);

  Array2Df* L = new Array2Df(width,height);
  gaussianBlur( *pyramids[0], *L );

  for ( int k=1 ; k<nlevels ; k++ )
  {
    width /= 2;
    height /= 2;
    pyramids[k] = new Array2Df(width,height);
    downSample(*L, *pyramids[k]);

    delete L;
    L = new Array2Df(width,height);
    gaussianBlur( *pyramids[k], *L );
  }

  delete L;
}

//--------------------------------------------------------------------

float calculateGradients(Array2Df* H, Array2Df* G, int k)
{
  const int width = H->getCols();
  const int height = H->getRows();
  const float divider = pow( 2.0f, k+1 );
  float avgGrad = 0.0f;

//#pragma omp parallel for shared(G,H) reduction(+:avgGrad)
  for( int y=0 ; y<height ; y++ )
  {
    for( int x=0 ; x<width ; x++ )
    {
      float gx, gy;
      int w, n, e, s;
      w = (x == 0 ? 0 : x-1);
      n = (y == 0 ? 0 : y-1);
      s = (y+1 == height ? y : y+1);
      e = (x+1 == width ? x : x+1);

      gx = ((*H)(w,y)-(*H)(e,y)) / divider;

      gy = ((*H)(x,s)-(*H)(x,n)) / divider;
      // note this implicitely assumes that H(-1)=H(0)
      // for the fft-pde slover this would need adjustment as H(-1)=H(1)
      // is assumed, which means gx=0.0, gy=0.0 at the boundaries
      // however, the impact is not visible so we ignore this here

      (*G)(x,y) = sqrt(gx*gx+gy*gy);
      avgGrad += (*G)(x,y);
    }
  }

  return avgGrad / (width*height);
}

//--------------------------------------------------------------------

void upSample(const Array2Df& A, Array2Df& B)
{
    const int width = B.getCols();
    const int height = B.getRows();
    const int awidth = A.getCols();
    const int aheight = A.getRows();

    //#pragma omp parallel for shared(A, B)
    for ( int y=0 ; y<height ; y++ )
    {
        for ( int x=0 ; x<width ; x++ )
        {
            int ax = static_cast<int>(x * 0.5f); //x / 2.f;
            int ay = static_cast<int>(y * 0.5f); //y / 2.f;
            ax = (ax<awidth) ? ax : awidth-1;
            ay = (ay<aheight) ? ay : aheight-1;

            B(x,y) = A(ax,ay);
        }
    }
//--- this code below produces 'use of uninitialized value error'
//   int width = A->getCols();
//   int height = A->getRows();
//   int x,y;

//   for( y=0 ; y<height ; y++ )
//     for( x=0 ; x<width ; x++ )
//     {
//       (*B)(2*x,2*y) = (*A)(x,y);
//       (*B)(2*x+1,2*y) = (*A)(x,y);
//       (*B)(2*x,2*y+1) = (*A)(x,y);
//       (*B)(2*x+1,2*y+1) = (*A)(x,y);
//     }
}


void calculateFiMatrix(Array2Df* FI, Array2Df* gradients[],
                       float avgGrad[], int nlevels, int detail_level,
                       float alfa, float beta, float noise)
{
    const bool newfattal = true;
    int width = gradients[nlevels-1]->getCols();
    int height = gradients[nlevels-1]->getRows();
    Array2Df** fi = new Array2Df*[nlevels];

    fi[nlevels-1] = new Array2Df(width,height);
    if (newfattal)
    {
        //#pragma omp parallel for shared(fi)
        for ( int k = 0 ; k < width*height ; k++ )
        {
            (*fi[nlevels-1])(k) = 1.0f;
        }
    }

    for ( int k = nlevels-1; k >= 0 ; k-- )
    {
        width = gradients[k]->getCols();
        height = gradients[k]->getRows();

        // only apply gradients to levels>=detail_level but at least to the coarsest
        if ( k >= detail_level
             ||k==nlevels-1
             || newfattal == false)
        {
            //DEBUG_STR << "calculateFiMatrix: apply gradient to level " << k << endl;
            //#pragma omp parallel for shared(fi,avgGrad)
            for ( int y = 0; y < height; y++ )
            {
                for ( int x = 0; x < width; x++ )
                {
                    float grad = ((*gradients[k])(x,y) < 1e-4f) ? 1e-4 : (*gradients[k])(x,y);
                    float a = alfa * avgGrad[k];

                    float value = powf((grad+noise)/a, beta - 1.0f);

                    if (newfattal)
                        (*fi[k])(x,y) *= value;
                    else
                        (*fi[k])(x,y) = value;
                }
            }
        }


        // create next level
        if ( k>1 )
        {
            width = gradients[k-1]->getCols();
            height = gradients[k-1]->getRows();
            fi[k-1] = new Array2Df(width,height);
        }
        else
            fi[0] = FI;                         // highest level -> result

        if ( k>0  && newfattal )
        {
            upSample(*fi[k], *fi[k-1]);           // upsample to next level
            gaussianBlur(*fi[k-1], *fi[k-1]);
        }
    }

    for ( int k=1 ; k<nlevels ; k++ )
    {
        delete fi[k];
    }
    delete[] fi;
}

inline
void findMaxMinPercentile(const Array2Df& I,
                                 float minPrct, float& minLum,
                                 float maxPrct, float& maxLum)
{
    const int size = I.getRows() * I.getCols();
    const float* data = I.data();
    std::vector<float> vI;

    std::copy(data, data + size, std::back_inserter(vI));
    std::sort(vI.begin(), vI.end());

    minLum = vI.at( int(minPrct*vI.size()) );
    maxLum = vI.at( int(maxPrct*vI.size()) );
}

void solve_pde_fft(Array2Df *F, Array2Df *U);

void tmo_fattal02(size_t width,
                  size_t height,
                  const Array2Df& Y,
                  Array2Df& L,
                  float alfa,
                  float beta,
                  float noise,
                  int detail_level)
{
// #ifdef TIMER_PROFILING
//     msec_timer stop_watch;
//     stop_watch.start();
// #endif
    static const float black_point = 0.1f;
    static const float white_point = 0.5f;
    static const float gamma = 1.0f; // 0.8f;
    // static const int   detail_level = 3;
    if ( detail_level < 0 ) detail_level = 0;
    if ( detail_level > 3 ) detail_level = 3;

  // ph.setValue(2);
  // if (ph.canceled()) return;

  int MSIZE = 32;         // minimum size of gaussian pyramid
  // I believe a smaller value than 32 results in slightly better overall
  // quality but I'm only applying this if the newly implemented fft solver
  // is used in order not to change behaviour of the old version
  // TODO: best let the user decide this value
  // if (fftsolver)
  {
     MSIZE = 8;
  }

  int size = width*height;
  // unsigned int x,y;
  // int i, k;

  // find max & min values, normalize to range 0..100 and take logarithm
  float minLum = Y(0,0);
  float maxLum = Y(0,0);
  for ( int i=0 ; i<size ; i++ )
  {
      minLum = ( Y(i) < minLum ) ? Y(i) : minLum;
      maxLum = ( Y(i) > maxLum ) ? Y(i) : maxLum;
  }
  Array2Df* H = new Array2Df(width, height);
  //#pragma omp parallel for private(i) shared(H, Y, maxLum)
  for ( int i=0 ; i<size ; i++ )
  {
      (*H)(i) = logf( 100.0f* Y(i)/maxLum + 1e-4 );
  }
  // ph.setValue(4);

  // create gaussian pyramids
  int mins = (width<height) ? width : height;    // smaller dimension
  int nlevels = 0;
  while ( mins >= MSIZE )
  {
    nlevels++;
    mins /= 2;
  }
  // std::cout << "DEBUG: nlevels = " << nlevels << ", mins = " << mins << std::endl;
  // The following lines solves a bug with images particularly small
  if (nlevels == 0) nlevels = 1;

  Array2Df** pyramids = new Array2Df*[nlevels];
  createGaussianPyramids(H, pyramids, nlevels);
  // ph.setValue(8);

  // calculate gradients and its average values on pyramid levels
  Array2Df** gradients = new Array2Df*[nlevels];
  float* avgGrad = new float[nlevels];
  for ( int k=0 ; k<nlevels ; k++ )
  {
    gradients[k] = new Array2Df(pyramids[k]->getCols(), pyramids[k]->getRows());
    avgGrad[k] = calculateGradients(pyramids[k],gradients[k], k);
  }
  // ph.setValue(12);

  // calculate fi matrix
  Array2Df* FI = new Array2Df(width, height);
  calculateFiMatrix(FI, gradients, avgGrad, nlevels, detail_level, alfa, beta, noise);
//  dumpPFS( "FI.pfs", FI, "Y" );
  for ( int i=0 ; i<nlevels ; i++ )
  {
    delete pyramids[i];
    delete gradients[i];
  }
  delete[] pyramids;
  delete[] gradients;
  delete[] avgGrad;
  // ph.setValue(16);
  // if (ph.canceled()){
  //   delete FI;
  //   delete H;
  //   return;
  // }


  // attenuate gradients
  Array2Df* Gx = new Array2Df(width, height);
  Array2Df* Gy = new Array2Df(width, height);

  // the fft solver solves the Poisson pde but with slightly different
  // boundary conditions, so we need to adjust the assembly of the right hand
  // side accordingly (basically fft solver assumes U(-1) = U(1), whereas zero
  // Neumann conditions assume U(-1)=U(0)), see also divergence calculation
  // if (fftsolver)
    for ( size_t y=0 ; y<height ; y++ )
      for ( size_t x=0 ; x<width ; x++ )
      {
        // sets index+1 based on the boundary assumption H(N+1)=H(N-1)
        unsigned int yp1 = (y+1 >= height ? height-2 : y+1);
        unsigned int xp1 = (x+1 >= width ?  width-2  : x+1);
        // forward differences in H, so need to use between-points approx of FI
        (*Gx)(x,y) = ((*H)(xp1,y)-(*H)(x,y)) * 0.5*((*FI)(xp1,y)+(*FI)(x,y));
        (*Gy)(x,y) = ((*H)(x,yp1)-(*H)(x,y)) * 0.5*((*FI)(x,yp1)+(*FI)(x,y));
      }
  // else
  //   for ( size_t y=0 ; y<height ; y++ )
  //     for ( size_t x=0 ; x<width ; x++ )
  //     {
  //       int s, e;
  //       s = (y+1 == height ? y : y+1);
  //       e = (x+1 == width ? x : x+1);

  //       (*Gx)(x,y) = ((*H)(e,y)-(*H)(x,y)) * (*FI)(x,y);
  //       (*Gy)(x,y) = ((*H)(x,s)-(*H)(x,y)) * (*FI)(x,y);
  //     }
  delete H;
  delete FI;
  // ph.setValue(18);


//   dumpPFS( "Gx.pfs", Gx, "Y" );
//   dumpPFS( "Gy.pfs", Gy, "Y" );

  // calculate divergence
  Array2Df DivG(width, height);
  for ( size_t y = 0; y < height; ++y )
  {
      for ( size_t x = 0; x < width; ++x )
      {
          DivG(x,y) = (*Gx)(x,y) + (*Gy)(x,y);
          if ( x > 0 ) DivG(x,y) -= (*Gx)(x-1,y);
          if ( y > 0 ) DivG(x,y) -= (*Gy)(x,y-1);

          // if (fftsolver)
          {
              if (x==0) DivG(x,y) += (*Gx)(x,y);
              if (y==0) DivG(x,y) += (*Gy)(x,y);
          }

      }
  }
  delete Gx;
  delete Gy;
  // ph.setValue(20);
  // if (ph.canceled())
  // {
  //     return;
  // }

//  dumpPFS( "DivG.pfs", DivG, "Y" );

  // solve pde and exponentiate (ie recover compressed image)
  {
  Array2Df U(width, height);
  // if (fftsolver)
  {
      solve_pde_fft(&DivG, &U);//, ph);
  }
  // else
  // {
  //     solve_pde_multigrid(&DivG, &U, ph);
  // }
// #ifndef NDEBUG
//   printf("\npde residual error: %f\n", residual_pde(&U, &DivG));
// #endif
  // ph.setValue(90);
  // if ( ph.canceled() )
  // {
  //     return;
  // }

  for ( size_t idx = 0 ; idx < height*width; ++idx )
  {
      L(idx) = expf( gamma * U(idx) );
  }
  }
  // ph.setValue(95);

  // remove percentile of min and max values and renormalize
  float cut_min = 0.01f * black_point;
  float cut_max = 1.0f - 0.01f * white_point;
  assert(cut_min>=0.0f && (cut_max<=1.0f) && (cut_min<cut_max));
  findMaxMinPercentile(L, cut_min, minLum, cut_max, maxLum);
  for ( size_t idx = 0; idx < height*width; ++idx )
  {
      L(idx) = (L(idx) - minLum) / (maxLum - minLum);
      if ( L(idx) <= 0.0f )
      {
          L(idx) = 0.0;
      }
      // note, we intentionally do not cut off values > 1.0
  }
// #ifdef TIMER_PROFILING
//     stop_watch.stop_and_update();
//     cout << endl;
//     cout << "tmo_fattal02 = " << stop_watch.get_time() << " msec" << endl;
// #endif

  // ph.setValue(96);
}


/**
 *
 * @file pde_fft.cpp
 * @brief Direct Poisson solver using the discrete cosine transform
 *
 * @author Tino Kluge (tino.kluge@hrz.tu-chemnitz.de)
 *
 */

//////////////////////////////////////////////////////////////////////
// Direct Poisson solver using the discrete cosine transform
//////////////////////////////////////////////////////////////////////
// by Tino Kluge (tino.kluge@hrz.tu-chemnitz.de)
//
// let U and F be matrices of order (n1,n2), ie n1=height, n2=width
// and L_x of order (n2,n2) and L_y of order (n1,n1) and both
// representing the 1d Laplace operator with Neumann boundary conditions,
// ie L_x and L_y are tridiagonal matrices of the form
//
//  ( -2  2          )
//  (  1 -2  1       )
//  (     .  .  .    )
//  (        1 -2  1 )
//  (           2 -2 )
//
// then this solver computes U given F based on the equation
//
//  -------------------------
//  L_y U + (L_x U^tr)^tr = F
//  -------------------------
//
// Note, if the first and last row of L_x and L_y contained one's instead of
// two's then this equation would be exactly the 2d Poisson equation with
// Neumann boundary conditions. As a simple rule:
// - Neumann: assume U(-1)=U(0) --> U(i-1) - 2 U(i) + U(i+1) becomes
//        i=0: U(0) - 2 U(0) + U(1) = -U(0) + U(1)
// - our system: assume U(-1)=U(1) --> this becomes
//        i=0: U(1) - 2(0) + U(1) = -2 U(0) + 2 U(1)
//
// The multi grid solver solve_pde_multigrid() solves the 2d Poisson pde
// with the right Neumann boundary conditions, U(-1)=U(0), see function
// atimes(). This means the assembly of the right hand side F is different
// for both solvers.

// #include <iostream>

// #include <boost/math/constants/constants.hpp>

// #include <stdio.h>
// #include <stdlib.h>
// #include "arch/math.h"
// #include <cassert>
// #ifdef _OPENMP
// #include <omp.h>
// #endif
// #include <vector>
// #include <fftw3.h>

// #include "Libpfs/progress.h"
// #include "Libpfs/array2d.h"
// #include "pde.h"

// using namespace std;


// #ifndef SQR
// #define SQR(x) (x)*(x)
// #endif


// returns T = EVy A EVx^tr
// note, modifies input data
void transform_ev2normal(Array2Df *A, Array2Df *T)
{
  int width = A->getCols();
  int height = A->getRows();
  assert((int)T->getCols()==width && (int)T->getRows()==height);

  // the discrete cosine transform is not exactly the transform needed
  // need to scale input values to get the right transformation
  for(int y=1 ; y<height-1 ; y++ )
    for(int x=1 ; x<width-1 ; x++ )
      (*A)(x,y)*=0.25f;

  for(int x=1 ; x<width-1 ; x++ )
  {
    (*A)(x,0)*=0.5f;
    (*A)(x,height-1)*=0.5f;
  }
  for(int y=1 ; y<height-1 ; y++ )
  {
    (*A)(0,y)*=0.5;
    (*A)(width-1,y)*=0.5f;
  }

  // note, fftw provides its own memory allocation routines which
  // ensure that memory is properly 16/32 byte aligned so it can
  // use SSE/AVX operations (2/4 double ops in parallel), if our
  // data is not properly aligned fftw won't use SSE/AVX
  // (I believe new() aligns memory to 16 byte so avoid overhead here)
  //
  // double* in = (double*) fftwf_malloc(sizeof(double) * width*height);
  // fftwf_free(in);

  // executes 2d discrete cosine transform
  fftwf_plan p;
  p=fftwf_plan_r2r_2d(height, width, A->data(), T->data(),
                        FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);
}


// returns T = EVy^-1 * A * (EVx^-1)^tr
void transform_normal2ev(Array2Df *A, Array2Df *T)
{
  int width = A->getCols();
  int height = A->getRows();
  assert((int)T->getCols()==width && (int)T->getRows()==height);

  // executes 2d discrete cosine transform
  fftwf_plan p;
  p=fftwf_plan_r2r_2d(height, width, A->data(), T->data(),
                        FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);

  // need to scale the output matrix to get the right transform
  for(int y=0 ; y<height ; y++ )
    for(int x=0 ; x<width ; x++ )
      (*T)(x,y)*=(1.0f/((height-1)*(width-1)));

  for(int x=0 ; x<width ; x++ )
  {
    (*T)(x,0)*=0.5f;
    (*T)(x,height-1)*=0.5f;
  }
  for(int y=0 ; y<height ; y++ )
  {
    (*T)(0,y)*=0.5f;
    (*T)(width-1,y)*=0.5f;
  }
}

// returns the eigenvalues of the 1d laplace operator
std::vector<double> get_lambda(int n)
{
  assert(n>1);
  std::vector<double> v(n);
  for (int i=0; i<n; i++)
  {
    v[i]=-4.0*SQR(sin((double)i/(2*(n-1))*RT_PI));
  }

  return v;
}

// // makes boundary conditions compatible so that a solution exists
// void make_compatible_boundary(Array2Df *F)
// {
//   int width = F->getCols();
//   int height = F->getRows();

//   double sum=0.0;
//   for(int y=1 ; y<height-1 ; y++ )
//     for(int x=1 ; x<width-1 ; x++ )
//       sum+=(*F)(x,y);

//   for(int x=1 ; x<width-1 ; x++ )
//     sum+=0.5*((*F)(x,0)+(*F)(x,height-1));

//   for(int y=1 ; y<height-1 ; y++ )
//     sum+=0.5*((*F)(0,y)+(*F)(width-1,y));

//   sum+=0.25*((*F)(0,0)+(*F)(0,height-1)+(*F)(width-1,0)+(*F)(width-1,height-1));

//   //DEBUG_STR << "compatible_boundary: int F = " << sum ;
//   //DEBUG_STR << " (should be 0 to be solvable)" << std::endl;

//   double add=-sum/(height+width-3);
//   //DEBUG_STR << "compatible_boundary: adjusting boundary by " << add << std::endl;
//   for(int x=0 ; x<width ; x++ )
//   {
//     (*F)(x,0)+=add;
//     (*F)(x,height-1)+=add;
//   }
//   for(int y=1 ; y<height-1 ; y++ )
//   {
//     (*F)(0,y)+=add;
//     (*F)(width-1,y)+=add;
//   }
// }



// solves Laplace U = F with Neumann boundary conditions
// if adjust_bound is true then boundary values in F are modified so that
// the equation has a solution, if adjust_bound is set to false then F is
// not modified and the equation might not have a solution but an
// approximate solution with a minimum error is then calculated
// double precision version
void solve_pde_fft(Array2Df *F, Array2Df *U)/*, pfs::Progress &ph,
                                              bool adjust_bound)*/
{
   // ph.setValue(20);
  //DEBUG_STR << "solve_pde_fft: solving Laplace U = F ..." << std::endl;
  int width = F->getCols();
  int height = F->getRows();
  assert((int)U->getCols()==width && (int)U->getRows()==height);

  // activate parallel execution of fft routines
#ifdef RT_FFTW3F_OMP
  fftwf_init_threads();
  fftwf_plan_with_nthreads( omp_get_max_threads() );
// #else
//   fftwf_plan_with_nthreads( 2 );
#endif

  // in general there might not be a solution to the Poisson pde
  // with Neumann boundary conditions unless the boundary satisfies
  // an integral condition, this function modifies the boundary so that
  // the condition is exactly satisfied
  // if(adjust_bound)
  // {
  //   //DEBUG_STR << "solve_pde_fft: checking boundary conditions" << std::endl;
  //   make_compatible_boundary(F);
  // }

  // transforms F into eigenvector space: Ftr =
  //DEBUG_STR << "solve_pde_fft: transform F to ev space (fft)" << std::endl;
  Array2Df* F_tr = new Array2Df(width,height);
  transform_normal2ev(F, F_tr);
  // TODO: F no longer needed so could release memory, but as it is an
  // input parameter we won't do that
  // ph.setValue(50);
  // if (ph.canceled())
  // {
  //   delete F_tr;
  //   return;
  // }

  //DEBUG_STR << "solve_pde_fft: F_tr(0,0) = " << (*F_tr)(0,0);
  //DEBUG_STR << " (must be 0 for solution to exist)" << std::endl;

  // in the eigenvector space the solution is very simple
  //DEBUG_STR << "solve_pde_fft: solve in eigenvector space" << std::endl;
  Array2Df* U_tr = new Array2Df(width,height);
  std::vector<double> l1=get_lambda(height);
  std::vector<double> l2=get_lambda(width);
  for(int y=0 ; y<height ; y++ )
  {
    for(int x=0 ; x<width ; x++ )
    {
      if(x==0 && y==0)
        (*U_tr)(x,y)=0.0; // any value ok, only adds a const to the solution
      else
        (*U_tr)(x,y)=(*F_tr)(x,y)/(l1[y]+l2[x]);
    }
  }
  delete F_tr;    // no longer needed so release memory
  // ph.setValue(55);


  // transforms U_tr back to the normal space
  //DEBUG_STR << "solve_pde_fft: transform U_tr to normal space (fft)" << std::endl;
  transform_ev2normal(U_tr, U);
  delete U_tr;    // no longer needed so release memory
  // ph.setValue(85);

  // the solution U as calculated will satisfy something like int U = 0
  // since for any constant c, U-c is also a solution and we are mainly
  // working in the logspace of (0,1) data we prefer to have
  // a solution which has no positive values: U_new(x,y)=U(x,y)-max
  // (not really needed but good for numerics as we later take exp(U))
  //DEBUG_STR << "solve_pde_fft: removing constant from solution" << std::endl;
  double max=0.0;
  for(int i=0; i<width*height; i++)
    if(max<(*U)(i))
      max=(*U)(i);

  for(int i=0; i<width*height; i++)
    (*U)(i)-=max;


  // fft parallel threads cleanup, better handled outside this function?
#ifdef RT_FFTW3F_OMP
  fftwf_cleanup_threads();
#endif

  // ph.setValue(90);
  //DEBUG_STR << "solve_pde_fft: done" << std::endl;
}


// ---------------------------------------------------------------------
// the functions below are only for test purposes to check the accuracy
// of the pde solvers


// // returns the norm of (Laplace U - F) of all interior points
// // useful to compare solvers
// float residual_pde(Array2Df* U, Array2Df* F)
// {
//   int width = U->getCols();
//   int height = U->getRows();
//   assert((int)F->getCols()==width && (int)F->getRows()==height);

//   double res=0.0;
//   for(int y=1;y<height-1;y++)
//     for(int x=1;x<width-1;x++)
//     {
//       double laplace=-4.0*(*U)(x,y)+(*U)(x-1,y)+(*U)(x+1,y)
//                      +(*U)(x,y-1)+(*U)(x,y+1);
//       res += SQR( laplace-(*F)(x,y) );
//     }
//   return static_cast<float>( sqrt(res) );
// }


} // namespace


void ImProcFunctions::ToneMapFattal02(Imagefloat *rgb, int detail_level)
{
    int w = rgb->getWidth();
    int h = rgb->getHeight();
    
    Array2Df Yr(w, h);
    Array2Df L(w, h);

    rgb->normalizeFloatTo1();

    #pragma omp parallel for if (multiThread)
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            Yr(x, y) = Color::rgbLuminance(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x));
        }
    }

    float alpha = params->fattal.alpha;
    float beta = params->fattal.beta;
    float noise = alpha * 0.01f;

    if (settings->verbose) {
        std::cout << "ToneMapFattal02: alpha = " << alpha << ", beta = " << beta
                  << ", detail_level = " << detail_level << std::endl;
    }

    tmo_fattal02(w, h, Yr, L, alpha, beta, noise, detail_level);

    const float epsilon = 1e-4f;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            float Y = std::max(Yr(x, y), epsilon);
            float l = std::max(L(x, y), epsilon);
            rgb->r(y, x) = std::max(rgb->r(y, x)/Y, 0.f) * l;
            rgb->g(y, x) = std::max(rgb->g(y, x)/Y, 0.f) * l;
            rgb->b(y, x) = std::max(rgb->b(y, x)/Y, 0.f) * l;
        }
    }

    rgb->normalizeFloatTo65535();
}

} // namespace rtengine
