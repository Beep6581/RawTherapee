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

#include "pipettebuffer.h"

namespace rtengine
{

PipetteBuffer::PipetteBuffer(::EditDataProvider *dataProvider) :
    dataProvider(dataProvider), imgFloatBuffer(NULL), LabBuffer(NULL), singlePlaneBuffer(), ready(false) {}

PipetteBuffer::~PipetteBuffer()
{
    flush();
#ifndef NDEBUG
    imgFloatBuffer = (Imagefloat*)(0xbaadf00d);
#endif
#ifndef NDEBUG
    LabBuffer = (LabImage*)(0xbaadf00d);
#endif
}

void PipetteBuffer::createBuffer(int width, int height)
{
    //printf("Appel de createBuffer %d x %d\n", width, height);
    resize (width, height);
}

void PipetteBuffer::flush()
{
    if (imgFloatBuffer) {
        delete imgFloatBuffer;
        imgFloatBuffer = NULL;
    }

    if (LabBuffer) {
        delete LabBuffer;
        LabBuffer = NULL;
    }

    singlePlaneBuffer.flushData();
    ready = false;
}

EditUniqueID PipetteBuffer::getEditID()
{
    if (dataProvider && dataProvider->getCurrSubscriber()) {
        return dataProvider->getCurrSubscriber()->getEditID();
    } else {
        return EUID_None;
    }
}

void PipetteBuffer::resize(int newWidth, int newHeight)
{
    resize(newWidth, newHeight, dataProvider ? dataProvider->getCurrSubscriber() : NULL);
}

// Resize buffers if they already exist
void PipetteBuffer::resize(int newWidth, int newHeight, EditSubscriber* newSubscriber)
{
    if (newSubscriber) {
        if (newSubscriber->getEditingType() == ET_PIPETTE) {
            if (newSubscriber->getPipetteBufferType() == BT_IMAGEFLOAT) {
                if (!imgFloatBuffer) {
                    imgFloatBuffer = new Imagefloat(newWidth, newHeight);
                } else {
                    imgFloatBuffer->allocate(newWidth, newHeight);
                }
            } else if (imgFloatBuffer) {
                delete imgFloatBuffer;
                imgFloatBuffer = NULL;
            }

            if (newSubscriber->getPipetteBufferType() == BT_LABIMAGE) {
                if (LabBuffer && (LabBuffer->W != newWidth && LabBuffer->H != newHeight)) {
                    delete LabBuffer;
                    LabBuffer = NULL;
                }

                if (!LabBuffer) {
                    LabBuffer = new LabImage(newWidth, newHeight);
                }
            } else if (LabBuffer) {
                delete LabBuffer;
                LabBuffer = NULL;
            }

            if (newSubscriber->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
                singlePlaneBuffer.allocate(newWidth, newHeight);
            } else if (singlePlaneBuffer.data) {
                singlePlaneBuffer.allocate(0, 0);
            }
        } else {
            // Should never happen
            flush();
        }
    }

    ready = false;
}

bool PipetteBuffer::bufferCreated()
{
    EditSubscriber* subscriber;

    if (dataProvider && (subscriber = dataProvider->getCurrSubscriber())) {
        if (subscriber->getEditingType() == ET_PIPETTE) {
            switch (dataProvider->getCurrSubscriber()->getPipetteBufferType()) {
            case (BT_IMAGEFLOAT):
                return imgFloatBuffer != NULL;
            case (BT_LABIMAGE):
                return LabBuffer != NULL;
            case (BT_SINGLEPLANE_FLOAT):
                return singlePlaneBuffer.data != NULL;
            }
        } else {
            return false;
        }
    }

    return false;
}

void PipetteBuffer::getPipetteData(float* v, int x, int y, int squareSize)
{
    if (ready && dataProvider && dataProvider->getCurrSubscriber()) {
        switch (dataProvider->getCurrSubscriber()->getPipetteBufferType()) {
        case (BT_IMAGEFLOAT):
            if (imgFloatBuffer) {
                imgFloatBuffer->getPipetteData(v[0], v[1], v[2], x, y, squareSize, 0);
                return;
            }

            break;

        case (BT_LABIMAGE):
            if (LabBuffer) {
                LabBuffer->getPipetteData(v[0], v[1], v[2], x, y, squareSize);
                return;
            }

            break;

        case (BT_SINGLEPLANE_FLOAT):
            if (singlePlaneBuffer.data != NULL) {
                singlePlaneBuffer.getPipetteData(v[0], x, y, squareSize, 0);
                v[1] = v[2] = -1.f;
                return;
            }
        }
    }

    v[0] = v[1] = v[2] = -1.f;
}

}
