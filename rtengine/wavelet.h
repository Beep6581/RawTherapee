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

template<class T>
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

template<class T>
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


// buffer must be of length (n + 4)
template<class T>
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
    
    for(int j = -1; j < (int)n + 1; j += 2)
    {
        tmp[j]  = tmp[j] - (tmp[j-1] + tmp[j+1]) / 2;
    }
    
    for(int i = 0; i < (int)n; i += 2)
    {
        tmp[i]  = tmp[i] + (tmp[i-1] + tmp[i+1] + 2) / 4;
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

template<class T>
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

    for(int i = 0; i < (int)n + 1; i += 2)
    {
        tmp[i]  = tmp[i] - (tmp[i-1] + tmp[i+1] + 2) / 4;
    }

    for(int j = 1; j < (int)n; j += 2)
    {
        tmp[j]  = tmp[j] + (tmp[j-1] + tmp[j+1]) / 2;
    }
    
    // copy data

    for(size_t i = 0, j = 0; i < n; i++, j += pitch)
    {
        data[j] = tmp[i];
    }
}

#endif
