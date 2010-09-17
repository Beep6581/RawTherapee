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
#include "improcfun.h"
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
    filterListStdImage.push_back ("CoarseTrans");
    filterListStdImage.push_back ("WhiteBalance");
    filterListStdImage.push_back ("ColorspaceTrans");
    filterListStdImage.push_back ("--Cache--");
    filterListStdImage.push_back ("Transform");
    filterListStdImage.push_back ("--Cache--");
    filterListStdImage.push_back ("ColorMixer");
    filterListStdImage.push_back ("ShadowsHighlights");
    filterListStdImage.push_back ("ToneCurve");
    filterListStdImage.push_back ("--Cache--");
    filterListStdImage.push_back ("LuminanceDenoiser");
    filterListStdImage.push_back ("ColorDenoiser");
    filterListStdImage.push_back ("Sharpener");
    filterListStdImage.push_back ("LuminanceCurve");
    filterListStdImage.push_back ("ColorCurve");

    filterListRawImage.push_back ("Demosaicing");
    filterListRawImage.push_back ("HighlightRecovery");
    filterListRawImage.push_back ("--Cache--");
    filterListRawImage.push_back ("CoarseTrans");
    filterListRawImage.push_back ("WhiteBalance");
    filterListRawImage.push_back ("ColorspaceTrans");
    filterListRawImage.push_back ("--Cache--");
    filterListRawImage.push_back ("Transform");
    filterListRawImage.push_back ("--Cache--");
    filterListRawImage.push_back ("ColorMixer");
    filterListRawImage.push_back ("ShadowsHighlights");
    filterListRawImage.push_back ("ToneCurve");
    filterListRawImage.push_back ("--Cache--");
    filterListRawImage.push_back ("LuminanceDenoiser");
    filterListRawImage.push_back ("ColorDenoiser");
    filterListRawImage.push_back ("Sharpener");
    filterListRawImage.push_back ("LuminanceCurve");
    filterListRawImage.push_back ("ColorCurve");

    previewSkip = 16;

    iccStore.parseDir (s->iccDirectory);
    CurveFactory::init ();
}


}

