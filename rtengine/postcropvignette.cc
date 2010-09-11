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
#include <postcropvignette.h>
#include <mytime.h>
#include <refreshmap.h>

namespace rtengine {

extern Settings* settings;

PostCropVignette::PostCropVignette ()
{
}

void PostCropVignette::update (ImProcCoordinator* parent) {

    ProcParams& params = parent->params;

    if (!params.postcropvignette.enabled)
        return;

    int amount = params.postcropvignette.amount;

    //Do stuff
    parent->imgsrc->getImage (parent->currWB, tr, origCrop, pp, params.hlrecovery, params.icm);

    cropMutex.unlock ();

    if (!internal)
        parent->mProcessing.unlock ();
}

}

