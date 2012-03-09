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

class MyPrevImgListener : public rtengine::PreviewImageListener {
    
    IImage8* i;

    public:
        // this method is called when the staged image processor creates a new image to store the resulting preview image (this does not happen too often)
        // usually you just have to store it
        void setImage   (IImage8* img, double scale, procparams::CropParams cp) {
            i = img;
        }
        // if the staged image processor wants to delete the image that stores the preview image, it calls this method. You have to destroy the image.
        void delImage   (IImage8* img) {
            if (img) {
                // make sure we dont use this image in an other thread
                IImage8* temp = i;
                i->getMutex().lock ();
                i = NULL;
                temp->getMutex().unlock ();
                // free it
                img->free ();
            }
        }
        // if the preview image changes, this method is called
        void imageReady (procparams::CropParams cp) {
            // initiate a redraw in the background and return as fast as possible
        }
        // a possible redraw function:
        //void redraw () {
        //    if (i) {
        //        i->lock ();
        //        int w = i->getWidth ();
        //        int h = i->getHeigt ();
        //        const char* data = i->getData ();
        //        ... draw it ...
        //        i->unlock ();
        //    }
        // }
};

int main (int argc, char* argv[]) {

    if (argc<4) {
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
    if (!ii)
        ii = rtengine::InitialImage::load (argv[1], false, errorCode, &pl);
    if (!ii) {
        std::cout << "Input file not supported." << std::endl;
        exit(2);
    }

/* Second scenario. Create a stagedimageprocessor with a preview scale of 1:5 and change few things */
    MyPrevImgListener myPrevImgListener;

    StagedImageProcessor* ipc = StagedImageProcessor::create (ii);
    ipc->setProgressListener (&pl);
    ipc->setPreviewImageListener (&myPrevImgListener);
    ipc->setPreviewScale (5); // preview scale = 1:5
    // you can add a histogram listener, too, that is notified when the histogram changes
    // ipc->setHistogramListener (...);
    // you can add autoexplistener that is notified about the exposure settings when the auto exp algorithm finishes
    // ipc->setAutoExpListener (curve);
    // you can add sizelistener if you want to be notified when the size of the final image changes (due to rotation/resize/etc)
    // ipc->setSizeListener (crop);

    // if you want to change the settings you have to ask for the procparams structure of the staged image processor
    // you have to tell it what has changed. At the first time tell it EvPhotoLoaded so a full processing will be performed
    rtengine::procparams::ProcParams* params = ipc->beginUpdateParams ();
    // change this and that...
    params->toneCurve.brightness = 1.0;   
    // you can load it, too, from a file: params->load (argv[2]);
    // finally you have to call this non-blocking method, and the image processing starts in the background. When finished, the preview image listener will be notified
    ipc->endUpdateParams (rtengine::EvPhotoLoaded);
    // you can go on with changing of the settings, following the gui actions
    // now we know that only the brightness has changed compared to the previous settings, to only a part of the processing has to be repeated
    params = ipc->beginUpdateParams ();
    params->toneCurve.brightness = 1.2;   
    ipc->endUpdateParams (rtengine::EvBrightness);
    
    // ... and so on. If you dont need it any more, you can destroy it (make sure that no processing is happening when you destroy it!)
    StagedImageProcessor::destroy (ipc);
}

