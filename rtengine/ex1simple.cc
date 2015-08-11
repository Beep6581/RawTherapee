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
#include <iostream>
//#include <giomm.h>
#include <helpers.h>

class PListener : public rtengine::ProgressListener
{

public:
    void setProgressStr (Glib::ustring str)
    {
        std::cout << str << std::endl;
    }
    void setProgress (double p)
    {
        std::cout << p << std::endl;
    }
};

int main (int argc, char* argv[])
{

    if (argc < 4) {
        std::cout << "Usage: rtcmd <infile> <paramfile> <outfile>" << std::endl;
        exit(1);
    }

    Glib::thread_init ();

    // create and fill settings
    rtengine::Settings* s = rtengine::Settings::create ();
    s->demosaicMethod = "hphd";
    s->colorCorrectionSteps = 2;
    s->iccDirectory = "";
    s->colorimetricIntent = 1;
    s->monitorProfile = "";
    // init rtengine
    rtengine::init (s);
    // the settings can be modified later through the "s" pointer without calling any api function

    // Create a listener object. Any class is appropriate that inherits from rtengine::ProgressListener
    PListener pl;

    // Load the image given in the first command line parameter
    rtengine::InitialImage* ii;
    int errorCode;
    ii = rtengine::InitialImage::load (argv[1], true, errorCode, &pl);

    if (!ii) {
        ii = rtengine::InitialImage::load (argv[1], false, errorCode, &pl);
    }

    if (!ii) {
        std::cout << "Input file not supported." << std::endl;
        exit(2);
    }

    // create an instance of ProcParams structure that holds the image processing settings. You find the memory map in a separate file and the non-basic types like strings and vectors can be manipulated through helper functions
    rtengine::procparams::ProcParams params;
    params.load (argv[2]);

    /* First, simplest scenario. Develope image and save it in a file */
    // create a processing job with the loaded image and the current processing parameters
    rtengine::ProcessingJob* job = ProcessingJob::create (i, params);
    // process image. The error is given back in errorcode.
    rtengine::IImage16* res = rtengine::processImage (job, errorCode, &pl);
    // save image to disk
    res->saveToFile (argv[3]);
    // through "res" you can access width/height and pixel data, too
}

