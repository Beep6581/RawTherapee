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
#ifndef _CROP_H_
#define _CROP_H_

#include <improccoordinator.h>
#include <rtengine.h>
#include <improcfun.h>
#include <image8.h>
#include <image16.h>
#include <imagesource.h>
#include <procevents.h>

namespace rtengine {

using namespace procparams;

class ImProcCoordinator;

class Crop : public DetailedCrop {

    protected:
        Image16* origCrop, *resizeCrop, *transCrop, *baseCrop;
        LabImage *laboCrop, *labnCrop;
        Image8 *cropImg;

        int** cbuffer;
        SHMap* cshmap;

        bool updating, needsNext;
        int skip;
        int cropx, cropy, cropw, croph;         // size of the detail crop image ('skip' taken into account), with border
        int trafx, trafy, trafw, trafh;         // the size and position to get from the imagesource that is transformed to the requested crop area
        int rqcropx, rqcropy, rqcropw, rqcroph; // size of the requested detail crop image (the image might be smaller) (without border)
        int borderRequested, upperBorder, leftBorder;
        

        bool cropAllocated;
        DetailedCropListener* cropImageListener;
        Glib::Mutex cropMutex;
        ImProcCoordinator* parent;

        bool setCropSizes (int cx, int cy, int cw, int ch, int skip, bool internal);
        void freeAll ();
        
    public:
             Crop        (ImProcCoordinator* parent);
            ~Crop        ();
    
        bool hasListener () { return cropImageListener; }
        void update      (int todo, bool internal=false);
        void setWindow   (int cx, int cy, int cw, int ch, int skip) { setCropSizes (cx, cy, cw, ch, skip, false); }
        void fullUpdate  ();
        void setListener (DetailedCropListener* il);
        void destroy     () { delete this; }
};
}
#endif
