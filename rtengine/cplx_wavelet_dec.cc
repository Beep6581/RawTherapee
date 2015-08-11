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

#include "cplx_wavelet_dec.h"

namespace rtengine {
	
	wavelet_decomposition::~wavelet_decomposition()
	{
		for(int i = 0; i <= lvltot; i++) {
			if(wavelet_decomp[i] != NULL)
				delete wavelet_decomp[i];
		}
		delete[] wavfilt_anal;
		delete[] wavfilt_synth;
		if(coeff0)
			delete [] coeff0;
	}

};

