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
#include "../rtgui/profilestorecombobox.h"
#include "rtengine.h"
#include "iccstore.h"
#include "dcp.h"
#include "camconst.h"
#include "curves.h"
#include "rawimagesource.h"
#include "improcfun.h"
#include "improccoordinator.h"
#include "dfmanager.h"
#include "ffmanager.h"
#include "rtthumbnail.h"
#include "profilestore.h"
#include "../rtgui/threadutils.h"
#include "rtlensfun.h"

namespace rtengine
{

const Settings* settings;

MyMutex* lcmsMutex = nullptr;

int init (const Settings* s, Glib::ustring baseDir, Glib::ustring userSettingsDir, bool loadAll)
{
    settings = s;
    ProfileStore::getInstance()->init (loadAll);
    ICCStore::getInstance()->init (s->iccDirectory, Glib::build_filename (baseDir, "iccprofiles"), loadAll);
    DCPStore::getInstance()->init (Glib::build_filename (baseDir, "dcpprofiles"), loadAll);

    CameraConstantsStore::getInstance ()->init (baseDir, userSettingsDir);
    ProcParams::init ();
    Color::init ();
    PerceptualToneCurve::init ();
    RawImageSource::init ();
    LFDatabase::init(s->lensfunDbDirectory);
    delete lcmsMutex;
    lcmsMutex = new MyMutex;
    dfm.init( s->darkFramesPath );
    ffm.init( s->flatFieldsPath );
    return 0;
}

void cleanup ()
{

    ProcParams::cleanup ();
    Color::cleanup ();
    RawImageSource::cleanup ();
}

StagedImageProcessor* StagedImageProcessor::create (InitialImage* initialImage)
{

    ImProcCoordinator* ipc = new ImProcCoordinator ();
    ipc->assign (initialImage->getImageSource ());
    return ipc;
}

void StagedImageProcessor::destroy (StagedImageProcessor* sip)
{

    delete sip;
}

Settings* Settings::create  ()
{

    return new Settings;
}

void Settings::destroy (Settings* s)
{

    delete s;
}


}

