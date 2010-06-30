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

namespace rtengine {

void ImProcFunctions :: waveletEqualizer(Image16 * image, int fw, int fh, const EqualizerParams & params) {

    wavelet_decomposition r(image->r, fw, fh);
    r.reconstruct(image->r, params.c);

    wavelet_decomposition g(image->g, fw, fh);
    g.reconstruct(image->g, params.c);

    wavelet_decomposition b(image->b, fw, fh);
    b.reconstruct(image->b, params.c);

}

}
