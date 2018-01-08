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
#ifndef _PIPETTEBUFFER_H_
#define _PIPETTEBUFFER_H_

#include "../rtgui/edit.h"
#include "array2D.h"
#include "iimage.h"
#include "coord.h"

namespace rtengine
{

/// @brief Structure that contains information about and pointers to the Edit buffer
class PipetteBuffer
{
protected:

    // To avoid duplicated information, we points to a EditDataProvider that contains the current EditSubscriber
    // instead of pointing to the EditSubscriber directly
    rtedit::EditDataProvider* dataProvider;

    // TODO: Unfortunately, buffer can be of several type, each one representing a floating point image. Maybe we could unify everything one day!?
    // Only one of the following pointers will be allocated at a time, if any; "one chunk" allocation
    Imagefloat* imgFloatBuffer;
    LabImage* LabBuffer;
    PlanarWhateverData<float> singlePlaneBuffer;

    bool ready;  // flag that indicates if the _pipette_ buffer is ready

    void createBuffer(int width, int height);
    void resize(int newWidth, int newHeight, rtedit::EditSubscriber* newSubscriber);
    void resize(int newWidth, int newHeight);
    void flush();

public:
    explicit PipetteBuffer(rtedit::EditDataProvider *dataProvider);
    ~PipetteBuffer();

    /** @brief Getter to know if the pipette buffer is correctly filled */
    bool isReady()
    {
        return ready;
    }

    /** @brief Setter to tell that the pipette buffer is correctly filled
     *  You have to use this method once the pipette is filled, so it can be read. */
    void setReady()
    {
        ready = true;
    }

    rtedit::EditDataProvider* getDataProvider()
    {
        return dataProvider;
    }
    EditUniqueID getEditID();
    Imagefloat* getImgFloatBuffer()
    {
        return imgFloatBuffer;
    }
    LabImage* getLabBuffer()
    {
        return LabBuffer;
    }
    PlanarWhateverData<float>* getSinglePlaneBuffer()
    {
        return &singlePlaneBuffer;
    }

    // return true if the buffer has been allocated
    bool bufferCreated();

    // get the pipette values
    void getPipetteData(float* v, int x, int y, int squareSize);
};

}

#endif
