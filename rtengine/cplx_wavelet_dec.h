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
 *  2012 Emil Martinec <ejmartin@uchicago.edu>
 */

#ifndef CPLX_WAVELET_DEC_H_INCLUDED
#define CPLX_WAVELET_DEC_H_INCLUDED

#include <cstddef>
#include <cmath>

#include "cplx_wavelet_level.h"
#include "cplx_wavelet_filter_coeffs.h"
#include "noncopyable.h"

namespace rtengine
{

class wavelet_decomposition :
    public NonCopyable
{
public:

    typedef float internal_type;
    float *coeff0;
    bool memoryAllocationFailed;

private:

    static const int maxlevels = 10;//should be greater than any conceivable order of decimation

    int lvltot, subsamp;
    int numThreads;
    int m_w, m_h;//dimensions

    int wavfilt_len, wavfilt_offset;
    float *wavfilt_anal;
    float *wavfilt_synth;


    wavelet_level<internal_type> * wavelet_decomp[maxlevels];

public:

    template<typename E>
    wavelet_decomposition(E * src, int width, int height, int maxlvl, int subsampling, int skipcrop = 1, int numThreads = 1, int Daub4Len = 6);

    ~wavelet_decomposition();

    internal_type ** level_coeffs(int level) const
    {
        return wavelet_decomp[level]->subbands();
    }

    int level_W(int level) const
    {
        return wavelet_decomp[level]->width();
    }

    int level_H(int level) const
    {
        return wavelet_decomp[level]->height();
    }

    int level_stride(int level) const
    {
        return wavelet_decomp[level]->stride();
    }

    int maxlevel() const
    {
        return lvltot + 1;
    }

    int subsample() const
    {
        return subsamp;
    }
    template<typename E>
    void reconstruct(E * dst, const float blend = 1.f);
};

template<typename E>
wavelet_decomposition::wavelet_decomposition(E * src, int width, int height, int maxlvl, int subsampling, int skipcrop, int numThreads, int Daub4Len)
    : coeff0(nullptr), memoryAllocationFailed(false), lvltot(0), subsamp(subsampling), numThreads(numThreads), m_w(width), m_h(height)
{

    //initialize wavelet filters
    wavfilt_len = Daub4Len;
    wavfilt_offset = Daub4_offset;
    wavfilt_anal = new float[2 * wavfilt_len];
    wavfilt_synth = new float[2 * wavfilt_len];

    if(wavfilt_len == 6) {
        for (int n = 0; n < 2; n++) {
            for (int i = 0; i < wavfilt_len; i++) {
                wavfilt_anal[wavfilt_len * (n) + i]  = Daub4_anal[n][i];
                wavfilt_synth[wavfilt_len * (n) + i] = Daub4_anal[n][wavfilt_len - 1 - i];
                //n=0 lopass, n=1 hipass
            }
        }
    } else if(wavfilt_len == 8) {
        for (int n = 0; n < 2; n++) {
            for (int i = 0; i < wavfilt_len; i++) {
                wavfilt_anal[wavfilt_len * (n) + i]  = Daub4_anal8[n][i];
                wavfilt_synth[wavfilt_len * (n) + i] = Daub4_anal8[n][wavfilt_len - 1 - i];
                //n=0 lopass, n=1 hipass
            }
        }
    } else if(wavfilt_len == 12) {
        for (int n = 0; n < 2; n++) {
            for (int i = 0; i < wavfilt_len; i++) {
                wavfilt_anal[wavfilt_len * (n) + i]  = Daub4_anal12[n][i];
                wavfilt_synth[wavfilt_len * (n) + i] = Daub4_anal12[n][wavfilt_len - 1 - i];
                //n=0 lopass, n=1 hipass
            }
        }
    } else if(wavfilt_len == 16) {
        for (int n = 0; n < 2; n++) {
            for (int i = 0; i < wavfilt_len; i++) {
                wavfilt_anal[wavfilt_len * (n) + i]  = Daub4_anal16[n][i];
                wavfilt_synth[wavfilt_len * (n) + i] = Daub4_anal16[n][wavfilt_len - 1 - i];
                //n=0 lopass, n=1 hipass
            }
        }
    } else if(wavfilt_len == 4) {
        for (int n = 0; n < 2; n++) {
            for (int i = 0; i < wavfilt_len; i++) {
                wavfilt_anal[wavfilt_len * (n) + i]  = Daub4_anal0[n][i];
                wavfilt_synth[wavfilt_len * (n) + i] = Daub4_anal0[n][wavfilt_len - 1 - i];
                //n=0 lopass, n=1 hipass
            }
        }
    }

    // after coefficient rotation, data structure is:
    // wavelet_decomp[scale][channel={lo,hi1,hi2,hi3}][pixel_array]

    lvltot = 0;
    E *buffer[2];
    buffer[0] = new (std::nothrow) E[(m_w / 2 + 1) * (m_h / 2 + 1)];

    if(buffer[0] == NULL) {
        memoryAllocationFailed = true;
        return;
    }

    buffer[1] = new (std::nothrow) E[(m_w / 2 + 1) * (m_h / 2 + 1)];

    if(buffer[1] == NULL) {
        memoryAllocationFailed = true;
        delete[] buffer[0];
        buffer[0] = NULL;
        return;
    }

    int bufferindex = 0;

    wavelet_decomp[lvltot] = new wavelet_level<internal_type>(src, buffer[bufferindex ^ 1], lvltot/*level*/, subsamp, m_w, m_h, \
            wavfilt_anal, wavfilt_anal, wavfilt_len, wavfilt_offset, skipcrop, numThreads);

    if(wavelet_decomp[lvltot]->memoryAllocationFailed) {
        memoryAllocationFailed = true;
    }

    while(lvltot < maxlvl - 1) {
        lvltot++;
        bufferindex ^= 1;
        wavelet_decomp[lvltot] = new wavelet_level<internal_type>(buffer[bufferindex], buffer[bufferindex ^ 1]/*lopass*/, lvltot/*level*/, subsamp, \
                wavelet_decomp[lvltot - 1]->width(), wavelet_decomp[lvltot - 1]->height(), \
                wavfilt_anal, wavfilt_anal, wavfilt_len, wavfilt_offset, skipcrop, numThreads);

        if(wavelet_decomp[lvltot]->memoryAllocationFailed) {
            memoryAllocationFailed = true;
        }
    }

    coeff0 = buffer[bufferindex ^ 1];
    delete[] buffer[bufferindex];
}

template<typename E>
void wavelet_decomposition::reconstruct(E * dst, const float blend)
{

    if(memoryAllocationFailed) {
        return;
    }

    // data structure is wavcoeffs[scale][channel={lo,hi1,hi2,hi3}][pixel_array]

    if(lvltot >= 1) {
        int width = wavelet_decomp[1]->m_w;
        int height = wavelet_decomp[1]->m_h;

        E *tmpHi = new (std::nothrow) E[width * height];

        if(tmpHi == NULL) {
            memoryAllocationFailed = true;
            return;
        }

        for (int lvl = lvltot; lvl > 0; lvl--) {
            E *tmpLo = wavelet_decomp[lvl]->wavcoeffs[2]; // we can use this as buffer
            wavelet_decomp[lvl]->reconstruct_level(tmpLo, tmpHi, coeff0, coeff0, wavfilt_synth, wavfilt_synth, wavfilt_len, wavfilt_offset);
            delete wavelet_decomp[lvl];
            wavelet_decomp[lvl] = nullptr;
        }

        delete[] tmpHi;
    }

    int width = wavelet_decomp[0]->m_w;
    int height = wavelet_decomp[0]->m_h2;
    E *tmpLo;

    if(wavelet_decomp[0]->bigBlockOfMemoryUsed()) { // bigBlockOfMemoryUsed means that wavcoeffs[2] points to a block of memory big enough to hold the data
        tmpLo = wavelet_decomp[0]->wavcoeffs[2];
    } else {                                      // allocate new block of memory
        tmpLo = new (std::nothrow) E[width * height];

        if(tmpLo == NULL) {
            memoryAllocationFailed = true;
            return;
        }
    }

    E *tmpHi = new (std::nothrow) E[width * height];

    if(tmpHi == NULL) {
        memoryAllocationFailed = true;

        if(!wavelet_decomp[0]->bigBlockOfMemoryUsed()) {
            delete[] tmpLo;
        }

        return;
    }


    wavelet_decomp[0]->reconstruct_level(tmpLo, tmpHi, coeff0, dst, wavfilt_synth, wavfilt_synth, wavfilt_len, wavfilt_offset, blend);

    if(!wavelet_decomp[0]->bigBlockOfMemoryUsed()) {
        delete[] tmpLo;
    }

    delete[] tmpHi;
    delete wavelet_decomp[0];
    wavelet_decomp[0] = nullptr;
    delete[] coeff0;
    coeff0 = nullptr;
}

};

#endif
