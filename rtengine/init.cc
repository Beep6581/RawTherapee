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
#include "color.h"
#include "improcfun.h"
#include "improccoordinator.h"
#include "curves.h"
#include "dfmanager.h"
#include "ffmanager.h"
#include "rtthumbnail.h"

namespace rtengine {

const Settings* settings;

Glib::Mutex* lcmsMutex = NULL;

int init (const Settings* s, Glib::ustring baseDir) {

    settings = s;
    iccStore->init (s->iccDirectory, baseDir + "/iccprofiles");
	iccStore->findDefaultMonitorProfile();

    ProcParams::init ();
    Color::init ();
    ImProcFunctions::initMunsell();
    Thumbnail::initGamma ();
    delete lcmsMutex;
    lcmsMutex = new Glib::Mutex;
    dfm.init( s->darkFramesPath );
    ffm.init( s->flatFieldsPath );
	return 0;
}

void cleanup () {

    ProcParams::cleanup ();
	Color::cleanup ();
    Thumbnail::cleanupGamma ();
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

