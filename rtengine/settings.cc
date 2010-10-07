/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include "rtengine.h"
#include "iccstore.h"
#include "curves.h"
#include <glibmm.h>

namespace rtengine {

Settings* Settings::settings = NULL;

extern Glib::Mutex* dcrMutex;

Settings::Settings () {

	if (settings)
		delete settings;
	settings = this;

	delete dcrMutex;
	dcrMutex = new Glib::Mutex ();

	colorimetricIntent = 1;
	verbose = true;
    filterList.push_back ("Demosaicing");
    filterList.push_back ("HighlightRecovery");
    filterList.push_back ("--Cache--");
    filterList.push_back ("CoarseTrans");
    filterList.push_back ("WhiteBalance");
    filterList.push_back ("ColorSpaceConversion");
    filterList.push_back ("--Cache--");
    filterList.push_back ("Transform");
    filterList.push_back ("--Cache--");
    filterList.push_back ("ColorMixer");
    filterList.push_back ("ShadowsHighlights");
    filterList.push_back ("ToneCurve");
    filterList.push_back ("--Cache--");
    filterList.push_back ("LuminanceDenoiser");
    filterList.push_back ("ColorDenoiser");
    filterList.push_back ("Sharpener");
    filterList.push_back ("LuminanceCurve");
    filterList.push_back ("ColorCurve");

    previewSkip = 16;

    iccStore.parseDir (iccDirectory);
    CurveFactory::init ();
}


}

