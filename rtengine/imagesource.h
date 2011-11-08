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
#ifndef _IMAGESOURCE_
#define _IMAGESOURCE_

#include <rtengine.h>
#include <colortemp.h>
#include <glibmm.h>
#include <image16.h>
#include <imagefloat.h>
#include <procparams.h>
#include <coord2d.h>
#include <vector>
#include <imagedata.h>
#include "LUT.h"

namespace rtengine {

using namespace procparams;

#define TR_NONE     0
#define TR_R90      1
#define TR_R180     2
#define TR_R270     3
#define TR_VFLIP    4
#define TR_HFLIP    8
#define TR_ROT      3

class PreviewProps {

    public:
    
        int x, y, w, h, skip;
        
        PreviewProps (int _x, int _y, int _w, int _h, int _skip)
            : x(_x), y(_y), w(_w), h(_h), skip(_skip) {}

};

class ImageSource : public InitialImage {

    private:
        int references;

    protected:
		cmsHPROFILE embProfile;
        Glib::ustring fileName;
        ImageMetaData* idata;

    public:
                    ImageSource ( ImageMetaData* meta=NULL) : references (1), embProfile(NULL), idata(meta) {}

        virtual ~ImageSource            () {}
        virtual int         load        (Glib::ustring fname, bool batch = false) =0;
        virtual void        preprocess  (const RAWParams &raw, HRecParams hrp){};
        virtual void        demosaic    (const RAWParams &raw, HRecParams hrp){};
	virtual void		HLRecovery_inpaint (float** red, float** green, float** blue){};

        virtual void        getImage    (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, HRecParams hlp, ColorManagementParams cmp, RAWParams raw) {}
        virtual ColorTemp   getWB       () =0;
        virtual ColorTemp   getAutoWB   () =0;
        virtual ColorTemp   getSpotWB   (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue, int tran) =0;

        virtual double      getDefGain  () { return 1.0; }
        
        virtual double      getGamma    () { return 0.0; }
        
        virtual void        getFullSize (int& w, int& h, int tr = TR_NONE) {}
        virtual void        getSize     (int tran, PreviewProps pp, int& w, int& h) {}

        virtual void        setProgressListener (ProgressListener* pl) {}

                void        increaseRef () { references++; }
                void        decreaseRef () { references--; if (!references) delete this; }

        virtual void        getAutoExpHistogram (LUTu & histogram, int& histcompr)=0;
        virtual void        getRAWHistogram (LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw) { 
            histRedRaw.clear(); histGreenRaw.clear(); histBlueRaw.clear(); // only some sources will supply this
        }  

        // functions inherited from the InitialImage interface
        virtual Glib::ustring getFileName ()        { return fileName; }
        virtual cmsHPROFILE getEmbeddedProfile ()   { return embProfile; }
        virtual ImageMetaData* getMetaData ()       { return idata; }
        virtual ImageSource* getImageSource ()      { return this; }
};
};
#endif
