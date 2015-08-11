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
#include "clipboard.h"

Clipboard clipboard;

Clipboard::Clipboard () : partProfile (false) {}

Clipboard::~Clipboard ()
{
    partProfile.deleteInstance();
}

/*
 * set both the "pparams" and "pedited" field of the PartialProfile; each one can be NULL
 */
void Clipboard::setPartialProfile   (const rtengine::procparams::PartialProfile& pprofile)
{
    if (pprofile.pparams) {
        if (!partProfile.pparams) {
            partProfile.pparams = new rtengine::procparams::ProcParams();
        }

        *partProfile.pparams = *pprofile.pparams;
    } else {
        if (partProfile.pparams) {
            delete partProfile.pparams;
            partProfile.pparams = NULL;
        }
    }

    if (pprofile.pedited) {
        if (!partProfile.pedited) {
            partProfile.pedited = new ParamsEdited();
        }

        *partProfile.pedited = *pprofile.pedited;
    } else {
        if (partProfile.pedited) {
            delete partProfile.pedited;
            partProfile.pedited = NULL;
        }
    }
}

/*
 * this method copy the procparams to "pparams" and delete "pedited"
 */
void Clipboard::setProcParams (const rtengine::procparams::ProcParams& pparams)
{
    // copy procparams
    if (!partProfile.pparams) {
        partProfile.pparams = new rtengine::procparams::ProcParams();
    }

    *partProfile.pparams = pparams;

    // delete pedited
    if (partProfile.pedited) {
        delete partProfile.pedited;
        partProfile.pedited = NULL;
    }
}
