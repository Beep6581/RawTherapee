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
#ifndef _THUMBPROCESSINGPARAMETERS_
#define _THUMBPROCESSINGPARAMETERS_

#include <rawmetadatalocation.h>
#include <procparams.h>
#include <glibmm.h>
#include <lcms.h>
#include <image16.h>

namespace rtengine {

    class Thumbnail {

            cmsHPROFILE camProfile;
            double iColorMatrix[3][3];
            double camToD50[3][3];


            void transformPixel (int x, int y, int tran, int& tx, int& ty);

            static bool igammacomputed;
            static unsigned short igammatab[256];
            static unsigned char gammatab[65536];

            Image16* thumbImg; 
            double camwbRed;
            double camwbGreen;
            double camwbBlue;
            double autowbTemp;
            double autowbGreen;
            int* aeHistogram;
            int  aeHistCompression;
            int embProfileLength;
            unsigned char* embProfileData;
            cmsHPROFILE embProfile;
            double redMultiplier;
            double greenMultiplier;
            double blueMultiplier;
            double scale;
            double defGain;
            int scaleForSave;
            bool gammaCorrected;
            double colorMatrix[3][3];

        public:

            bool isRaw;
            
            ~Thumbnail ();
            Thumbnail ();
            
            void init ();
            
            IImage8* processImage   (const procparams::ProcParams& pparams, int rheight, TypeInterpolation interp, double& scale);
            int      getImageWidth  (const procparams::ProcParams& pparams, int rheight);
            void     getFinalSize   (const rtengine::procparams::ProcParams& pparams, int& w, int& h);
            
            static Thumbnail* loadFromRaw (const Glib::ustring& fname, RawMetaDataLocation& rml, int &w, int &h, int fixwh);
            static Thumbnail* loadFromImage (const Glib::ustring& fname, int &w, int &h, int fixwh);           
            
            void getCamWB     (double& temp, double& green);
            void getAutoWB    (double& temp, double& green);
            void getSpotWB    (const procparams::ProcParams& params, int x, int y, int rect, double& temp, double& green);
            void applyAutoExp (procparams::ProcParams& pparams);
            
            bool writeImage (const Glib::ustring& fname, int format);
            bool readImage (const Glib::ustring& fname);
            
            bool readData  (const Glib::ustring& fname);
            bool writeData  (const Glib::ustring& fname);
            
            bool readEmbProfile  (const Glib::ustring& fname);
            bool writeEmbProfile (const Glib::ustring& fname);

            bool readAEHistogram  (const Glib::ustring& fname);
            bool writeAEHistogram (const Glib::ustring& fname);
    };   
}

#endif

