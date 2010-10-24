/*
 *  This file is part of RawTherapee.
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
 *
 *  2010 Ilya Popov <ilia_popov@rambler.ru>
 */


#ifndef WAVELET_H_INCLUDED
#define WAVELET_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////
//                        Haar wavelets                                    //
/////////////////////////////////////////////////////////////////////////////

// Bad, strong block effect

template<typename T>
void dwt_haar(T * data, size_t pitch, T * buffer, size_t n)
{
    size_t n2a = (n + 1) / 2;
    size_t n2 = n/2;

    for(size_t i = 0, j = 0; i < n2; i++, j += 2 * pitch)
    {
        T a = data[j];
        T b = data[j + pitch];
        buffer[i] = (a + b) / 2;
        buffer[n2a + i] = (a - b);
    }
    
    if(n2 < n2a)
    {
        buffer[n2] = data[pitch * (n-1)];
    }
    
    for(size_t k = 0, q = 0; k < n; k++, q += pitch)
    {
        data[q] = buffer[k];
    }
}

template<typename T>
void idwt_haar(T * data, size_t pitch, T * buffer, size_t n, int alpha)
{
    size_t n2a = (n + 1) / 2;
    size_t n2 = n/2;

    for(size_t i = 0, j = 0; i < n2; i++, j += 2)
    {
        T p = data[i * pitch];
        T q = (alpha * data[(n2a + i)*pitch]) / 1024;
        buffer[j] = p + q / 2;
        buffer[j + 1] = p - q / 2;
    }
    
    if(n2 < n2a)
    {
        buffer[n-1] = data[pitch * n2];
    }
    
    for(size_t k = 0, q = 0; k < n; k++, q += pitch)
    {
        data[q] = buffer[k];
    }
}

/////////////////////////////////////////////////////////////////////////////
//                        CDF 5/3 wavelets                                 //
/////////////////////////////////////////////////////////////////////////////


// buffer must be of length (n + 4)
template<typename T>
void dwt_53(T * data, size_t pitch, T * buffer, size_t n)
{
    size_t n2 = n/2;
    size_t n2a = (n + 1) / 2;
    T * tmp = buffer + 2;
    
    // copy data
    
    for(size_t i = 0, j = 0; i < n; i++, j += pitch)
    {
        tmp[i] = data[j];
    }
    
    // extend mirror-like
    
    tmp[-1] = tmp[1];
    tmp[-2] = tmp[2];
    
    tmp[n] = tmp[n-2];
    tmp[n+1] = tmp[n-3];

    // calculate coefficients
    
    for(ptrdiff_t i = -1; i < (ptrdiff_t)n + 1; i += 2)
    {
        tmp[i] = tmp[i] - (tmp[i-1] + tmp[i+1]) / 2;
    }
    
    for(ptrdiff_t i = 0; i < (ptrdiff_t)n; i += 2)
    {
        tmp[i] = tmp[i] + (tmp[i-1] + tmp[i+1] + 2) / 4;
    }
    
    // copy with reordering
    
    for(size_t i = 0, j = 0; i < n; i+=2, j += pitch)
    {
        data[j] = tmp[i];
    }

    for(size_t i = 1, j = n2a*pitch; i < n; i+=2, j += pitch)
    {
        data[j] = tmp[i];
    }
}

template<typename T>
void idwt_53(T * data, size_t pitch, T * buffer, size_t n, int alpha)
{
    size_t n2 = n/2;
    size_t n2a = (n + 1) / 2;
    T * tmp = buffer + 2;

    // copy with reordering
    
    for(size_t i = 0, j = 0; i < n; i+=2, j += pitch)
    {
        tmp[i] = data[j];
    }

    for(size_t i = 1, j = n2a*pitch; i < n; i+=2, j += pitch)
    {
        tmp[i] = (alpha * data[j]) / 1024;
    }

    // extend mirror-like
    
    tmp[-1] = tmp[1];
    tmp[-2] = tmp[2];
    
    tmp[n] = tmp[n-2];
    tmp[n+1] = tmp[n-3];

    // calculate coefficients

    for(ptrdiff_t i = 0; i < (ptrdiff_t)n + 1; i += 2)
    {
        tmp[i] = tmp[i] - (tmp[i-1] + tmp[i+1] + 2) / 4;
    }

    for(ptrdiff_t i = 1; i < (ptrdiff_t)n; i += 2)
    {
        tmp[i] = tmp[i] + (tmp[i-1] + tmp[i+1]) / 2;
    }
    
    // copy data

    for(size_t i = 0, j = 0; i < n; i++, j += pitch)
    {
        data[j] = tmp[i];
    }
}

/////////////////////////////////////////////////////////////////////////////
//                        Edge-avoiding wavelets                           //
/////////////////////////////////////////////////////////////////////////////
// based on 
//   Edge-Avoiding Wavelets and their Applications
//     Raanan Fattal <raananf@cs.huji.ac.il>,
//     Hebrew University of Jerusalem, Israel
//
// WCDF variant from this paper is used here
// T must be one of floating-point types 
/////////////////////////////////////////////////////////////////////////////

template<typename T>
inline T wcdf_weight(T a, T b)
{
    static const T eps = 1;
    static const T one = 1.0;
    return one / (fabs(a - b) + eps);
}

// buffer is a temporary storage
// buffer2 must be preserved between dwt and idwt

template<typename T>
void dwt_wcdf(T * data, size_t pitch, T * buffer, size_t n, T * buffer2)
{
    size_t n2 = n/2;
    size_t n2a = (n + 1) / 2;
    T * tmp = buffer + 2;
    T * w = buffer2 + 2;
    
    // copy data
    
    for(size_t i = 0, j = 0; i < n; i++, j += pitch)
    {
        tmp[i] = data[j];
    }
    
    // extend mirror-like
    
    tmp[-1] = tmp[1];
    tmp[-2] = tmp[2];
    
    tmp[n] = tmp[n-2];
    tmp[n+1] = tmp[n-3];
    
    // calculate weights
    
    for(ptrdiff_t i = 0; i < (ptrdiff_t)n - 1; i++)
    {
        w[i] = wcdf_weight(tmp[i], tmp[i+1]);
        //w[i] = 1;
    }
    
    w[-1] = w[-2] = (T)0.0;
    w[n-1] = w[n] = w[n + 1] = (T)0.0;

    // calculate coefficients
    
    for(ptrdiff_t i = 1; i < (ptrdiff_t)n; i += 2)
    {
        tmp[i] = tmp[i] - (w[i-1]*tmp[i-1] + w[i]*tmp[i+1]) / (w[i-1] + w[i]);
    }
    
    for(ptrdiff_t i = 0; i < (ptrdiff_t)n; i += 2)
    {
        tmp[i] = tmp[i] + (T)0.5 * (w[i-1]*tmp[i-1] + w[i]*tmp[i+1]) / (w[i-1] + w[i]);
    }
    
    // copy with reordering
    
    for(size_t i = 0, j = 0; i < n; i+=2, j += pitch)
    {
        data[j] = tmp[i];
    }

    for(size_t i = 1, j = n2a*pitch; i < n; i+=2, j += pitch)
    {
        data[j] = tmp[i];
    }
}

template<typename T>
void idwt_wcdf(T * data, size_t pitch, T * buffer, size_t n, int alpha, T * buffer2)
{
    size_t n2 = n/2;
    size_t n2a = (n + 1) / 2;
    T * tmp = buffer + 2;
    T * w = buffer2 + 2;

    // copy with reordering
    
    for(size_t i = 0, j = 0; i < n; i+=2, j += pitch)
    {
        tmp[i] = data[j];
    }

    for(size_t i = 1, j = n2a*pitch; i < n; i+=2, j += pitch)
    {
        tmp[i] = (alpha * data[j]) / 1024;
    }

    // extend mirror-like
    
    tmp[-1] = tmp[1];
    tmp[-2] = tmp[2];
    
    tmp[n] = tmp[n-2];
    tmp[n+1] = tmp[n-3];

    // calculate coefficients

    for(ptrdiff_t i = 0; i < (ptrdiff_t)n; i += 2)
    {
        tmp[i] = tmp[i] - (T)0.5 * (w[i-1]*tmp[i-1] + w[i]*tmp[i+1]) / (w[i-1] + w[i]);
    }

    for(ptrdiff_t i = 1; i < (ptrdiff_t)n; i += 2)
    {
        tmp[i] = tmp[i] + (w[i-1]*tmp[i-1] + w[i]*tmp[i+1]) / (w[i-1] + w[i]);
    }
   
    // copy data

    for(size_t i = 0, j = 0; i < n; i++, j += pitch)
    {
        data[j] = tmp[i];
    }
}

//////////////////////////////////////////////////////////////////////////////

#endif
