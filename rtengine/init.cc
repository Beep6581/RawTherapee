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
#include <rtengine.h>
#include <iccstore.h>
#include <improcfun.h>
#include <improccoordinator.h>
#include <curves.h>

namespace rtengine {

const Settings* settings;

extern Glib::Mutex* dcrMutex;
Glib::Mutex* lcmsMutex = NULL;

int init (const Settings* s) {

    settings = s;
    iccStore.parseDir (s->iccDirectory);
    CurveFactory::loadCurves ("");
    ImProcFunctions::initCache ();
    delete dcrMutex;
    dcrMutex = new Glib::Mutex;
    delete lcmsMutex;
    lcmsMutex = new Glib::Mutex;
}

StagedImageProcessor* StagedImageProcessor::create (InitialImage* initialImage) {

    ImProcCoordinator* ipc = new ImProcCoordinator ();
    ipc->assign (initialImage->getImageSource ());
    return ipc;
}

void StagedImageProcessor::destroy (StagedImageProcessor* sip) {

    delete sip;
}

Settings* Settings::create  () {
    
    return new Settings;
}

void Settings::destroy (Settings* s) {

    delete s;
}


}

