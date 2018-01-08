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
#pragma once

#include "improccoordinator.h"
#include "rtengine.h"
#include "improcfun.h"
#include "image8.h"
#include "image16.h"
#include "imagesource.h"
#include "procevents.h"
#include "pipettebuffer.h"
#include "../rtgui/threadutils.h"

namespace rtengine
{

using namespace procparams;

class ImProcCoordinator;

class Crop : public DetailedCrop, public PipetteBuffer
{

protected:
    // --- permanently allocated in RAM and only renewed on size changes
    Imagefloat*  origCrop;   // "one chunk" allocation
    LabImage*    laboCrop;   // "one chunk" allocation
    LabImage*    labnCrop;   // "one chunk" allocation
    Image8*      cropImg;    // "one chunk" allocation ; displayed image in monitor color space, showing the output profile as well (soft-proofing enabled, which then correspond to workimg) or not
    float *      cbuf_real;  // "one chunk" allocation
    SHMap*       cshmap;     // per line allocation

    // --- automatically allocated and deleted when necessary, and only renewed on size changes
    Imagefloat*  transCrop;    // "one chunk" allocation, allocated if necessary
    CieImage*    cieCrop;      // allocating 6 images, each in "one chunk" allocation
    // -----------------------------------------------------------------
    float**      cbuffer;

    bool updating;         /// Flag telling if an updater thread is currently processing
    bool newUpdatePending; /// Flag telling the updater thread that a new update is pending
    int skip;
    int cropx, cropy, cropw, croph;         /// size of the detail crop image ('skip' taken into account), with border
    int trafx, trafy, trafw, trafh;         /// the size and position to get from the imagesource that is transformed to the requested crop area
    int rqcropx, rqcropy, rqcropw, rqcroph; /// size of the requested detail crop image (the image might be smaller) (without border)
    const int borderRequested;              /// requested extra border size for image processing
    int upperBorder, leftBorder;            /// extra border size really allocated for image processing

    bool cropAllocated;
    DetailedCropListener* cropImageListener;

    MyMutex cropMutex;
    ImProcCoordinator* const parent;
    const bool isDetailWindow;
    EditUniqueID getCurrEditID();
    bool setCropSizes (int cropX, int cropY, int cropW, int cropH, int skip, bool internal);
    void freeAll ();

public:
    Crop             (ImProcCoordinator* parent, rtedit::EditDataProvider *editDataProvider, bool isDetailWindow);
    virtual ~Crop    ();

    void setEditSubscriber(rtedit::EditSubscriber* newSubscriber);
    bool hasListener();
    void update      (int todo);
    void setWindow   (int cropX, int cropY, int cropW, int cropH, int skip)
    {
        setCropSizes (cropX, cropY, cropW, cropH, skip, false);
    }

    /** @brief Synchronously look out if a full update is necessary
     * First try, only make fullUpdate if this returns false
     */
    bool tryUpdate   ();
    /** @brief Asynchronously reprocess the detailed crop */
    void fullUpdate  ();  // called via thread

    void setListener    (DetailedCropListener* il);
    void destroy        ();
    int get_skip();
    int getLeftBorder();
    int getUpperBorder();
};
}
