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

#include <rtengine.h>
#include <improcfun.h>

#include <wavelet_dec.h>

#include <iostream>

namespace rtengine {

void ImProcFunctions :: waveletEqualizer(Image16 * image) {

    if (!params->equalizer.enabled) {
        return;
    }

    limiter<wavelet_decomposition::internal_type> l(0, 65535);

    wavelet_decomposition r(image->r, image->getW(), image->getH());
    r.reconstruct(image->r, params->equalizer.c, l);

    wavelet_decomposition g(image->g, image->getW(), image->getH());
    g.reconstruct(image->g, params->equalizer.c, l);

    wavelet_decomposition b(image->b, image->getW(), image->getH());
    b.reconstruct(image->b, params->equalizer.c, l);

}

void ImProcFunctions :: waveletEqualizer (LabImage * image, bool luminance, bool chromaticity) {

    if (!params->equalizer.enabled) {
        return;
    }
    
    clock_t start = clock();

    if (luminance) {
        limiter<wavelet_decomposition::internal_type> l1(0, 65535);

        wavelet_decomposition L(image->L, image->W, image->H);
        L.reconstruct(image->L, params->equalizer.c, l1);
    }

    if (chromaticity) {
        limiter<wavelet_decomposition::internal_type> l2(-32768, 32767);

        wavelet_decomposition a(image->a, image->W, image->H);
        a.reconstruct(image->a, params->equalizer.c, l2);

        wavelet_decomposition b(image->b, image->W, image->H);
        b.reconstruct(image->b, params->equalizer.c, l2);
    }
    
    std::cout << "Wavelets done in " << (double)(clock() - start) / CLOCKS_PER_SEC << std::endl;

}

}
