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
#include <dfmanager.h>
#include <ffmanager.h>
#include <rtthumbnail.h>
#include <exiv2/exiv2.hpp>
#include <imagedata.h>

namespace rtengine {

const Settings* settings;

Glib::Mutex* lcmsMutex = NULL;
Glib::Mutex* exiv2Mutex = NULL;

int init (const Settings* s, Glib::ustring baseDir) {

    settings = s;
    iccStore->init (s->iccDirectory, baseDir + "/iccprofiles");
	iccStore->findDefaultMonitorProfile();

    CurveFactory::init ();
    ImProcFunctions::initCache ();
    Thumbnail::initGamma ();
    delete lcmsMutex;
    lcmsMutex = new Glib::Mutex;
    delete exiv2Mutex;
    exiv2Mutex = new Glib::Mutex;

	//this should be done once
	Exiv2::XmpProperties::registerNs("http://www.rawtherapee.com/1.0/", "rt");
	Exiv2::XmpProperties::registerNs("http://ns.adobe.com/xap/1.0/", "xmp");
	Exiv2::XmpProperties::registerNs("http://ns.adobe.com/exif/1.0/aux/", "aux");
	Exiv2::XmpProperties::registerNs("http://ns.adobe.com/photoshop/1.0/", "photoshop");
	Exiv2::XmpProperties::registerNs("http://purl.org/dc/elements/1.1/", "dc");
	Exiv2::XmpProperties::registerNs("http://iptc.org/std/Iptc4xmpCore/1.0/xmlns/", "Iptc4xmpCore");
	Exiv2::XmpProperties::registerNs("http://ns.useplus.org/ldf/xmp/1.0/", "plus");
	Exiv2::XmpProperties::registerNs("http://iptc.org/std/Iptc4xmpExt/2008-02-29/","Iptc4xmpExt");
	Exiv2::XmpProperties::registerNs("http://ns.adobe.com/xap/1.0/rights/","xmpRights");
	IPTCMeta::initIPTCMeta();

    dfm.init( s->darkFramesPath );
    ffm.init( s->flatFieldsPath );
	return 0;
}

void cleanup () {
	Exiv2::XmpParser::terminate();
    ImProcFunctions::cleanupCache ();
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

