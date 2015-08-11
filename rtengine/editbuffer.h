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
#ifndef _EDITBUFFER_H_
#define _EDITBUFFER_H_

#include "../rtgui/edit.h"
#include "array2D.h"
#include "iimage.h"

namespace rtengine
{

/// @brief Structure that contains information about and pointers to the Edit buffer
class EditBuffer
{
private:

    // Used to draw the objects where the color correspond to the object's ID, in order to find the correct object when hovering
    Cairo::RefPtr<Cairo::ImageSurface> objectMap;
    // If more than 254 objects has to be handled, objectMap2 contains the "upper part" of the 16 bit int value. objectMap2 will be NULL otherwise.
    Cairo::RefPtr<Cairo::ImageSurface> objectMap2;
    ObjectMode objectMode;

protected:

    // To avoid duplicated information, we points to a EditDataProvider that contains the current EditSubscriber
    // instead of pointing to the EditSubscriber directly
    ::EditDataProvider* dataProvider;

    // TODO: Unfortunately, buffer can be of several type, each one representing a floating point image. Maybe we could unify everything one day!?
    // Only one of the following pointers will be allocated at a time, if any; "one chunk" allocation
    Imagefloat* imgFloatBuffer;
    LabImage* LabBuffer;
    PlanarWhateverData<float> singlePlaneBuffer;

    bool ready;  // flag that indicates if the _pipette_ buffer is ready

    void                       createBuffer(int width, int height);
    void                       resize(int newWidth, int newHeight, EditSubscriber* newSubscriber);
    void                       resize(int newWidth, int newHeight);
    void                       flush();

public:
    EditBuffer(::EditDataProvider *dataProvider);
    ~EditBuffer();

    /** @brief Getter to know if the pipette buffer is correctly filled */
    bool                       isReady()
    {
        return ready;
    }

    /** @brief Setter to tell that the pipette buffer is correctly filled
     *  You have to use this method once the pipette is filled, so it can be read. */
    void                       setReady()
    {
        ready = true;
    }

    void                       setObjectMode(ObjectMode newType);
    ::EditDataProvider*        getDataProvider()
    {
        return dataProvider;
    }
    EditUniqueID               getEditID();
    Imagefloat*                getImgFloatBuffer()
    {
        return imgFloatBuffer;
    }
    LabImage*                  getLabBuffer()
    {
        return LabBuffer;
    }
    PlanarWhateverData<float>* getSinglePlaneBuffer()
    {
        return &singlePlaneBuffer;
    }
    ObjectMode                 getObjectMode()
    {
        return objectMode;
    }

    Cairo::RefPtr<Cairo::ImageSurface> &getObjectMap ()
    {
        return objectMap;
    }
    Cairo::RefPtr<Cairo::ImageSurface> &getObjectMap2()
    {
        return objectMap2;
    }

    // return true if the buffer has been allocated
    bool                       bufferCreated();

    int                        getObjectID(const Coord& location);
    // get the pipette values
    void                       getPipetteData(float* v, int x, int y, int squareSize);
};

}

#endif
