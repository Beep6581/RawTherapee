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

class PListener : public rtengine::ProgressListener {
    
    public:
        void setProgressStr (Glib::ustring str) {
            std::cout << str << std::endl;
        }
        void setProgress (double p) {
            std::cout << p << std::endl;
        }
};

int main (int argc, char* argv[]) {

    if (argc<4) {
        std::cout << "Usage: rtcmd <infile> <paramfile> <outfile>" << std::endl;
        exit(1);
    }

    rtengine::Settings s;
    s.demosaicMethod = "hphd";
    s.colorCorrectionSteps = 2;
    s.iccDirectory = "";
    s.colorimetricIntent = 1;
    s.monitorProfile = "";

    Glib::thread_init ();
    rtengine::init (s,"");
    PListener pl;

    rtengine::InitialImage* ii;
    int errorCode;
    ii = rtengine::InitialImage::load (argv[1], true, errorCode, &pl);
    if (!ii)
        ii = rtengine::InitialImage::load (argv[1], false, errorCode, &pl);
    if (!ii) {
        std::cout << "Input file not supported." << std::endl;
        exit(2);
    }

    rtengine::procparams::ProcParams params;
    params.load (argv[2]);

    rtengine::ProcessingJob* job = ProcessingJob::create (ii, params);
    rtengine::IImage16* res = rtengine::processImage (job, errorCode, &pl);
    res->saveToFile (argv[3]);
}

