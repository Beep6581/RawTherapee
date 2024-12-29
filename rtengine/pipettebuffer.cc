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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pipettebuffer.h"

#include "imagefloat.h"
#include "labimage.h"

#include "rtgui/editcallbacks.h"

namespace rtengine
{

PipetteBuffer::PipetteBuffer(::EditDataProvider *dataProvider) :
    dataProvider(dataProvider), imgFloatBuffer(nullptr), LabBuffer(nullptr), singlePlaneBuffer(), ready(false) {}

PipetteBuffer::~PipetteBuffer()
{
    flush();
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
        imgFloatBuffer = nullptr;
    }

    if (LabBuffer) {
        delete LabBuffer;
        LabBuffer = nullptr;
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
    resize(newWidth, newHeight, dataProvider ? dataProvider->getCurrSubscriber() : nullptr);
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
                imgFloatBuffer = nullptr;
            }

            if (newSubscriber->getPipetteBufferType() == BT_LABIMAGE) {
                if (LabBuffer && (LabBuffer->W != newWidth && LabBuffer->H != newHeight)) {
                    delete LabBuffer;
                    LabBuffer = nullptr;
                }

                if (!LabBuffer) {
                    LabBuffer = new LabImage(newWidth, newHeight);
                }
            } else if (LabBuffer) {
                delete LabBuffer;
                LabBuffer = nullptr;
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
                return imgFloatBuffer != nullptr;
            case (BT_LABIMAGE):
                return LabBuffer != nullptr;
            case (BT_SINGLEPLANE_FLOAT):
                return singlePlaneBuffer.data != nullptr;
            }
        } else {
            return false;
        }
    }

    return false;
}

void PipetteBuffer::getPipetteData(int x, int y, const int squareSize)
{
    if (ready && dataProvider && dataProvider->getCurrSubscriber()) {
        switch (dataProvider->getCurrSubscriber()->getPipetteBufferType()) {
        case (BT_IMAGEFLOAT):
            if (imgFloatBuffer) {
                float v[3];
                imgFloatBuffer->getPipetteData(v[0], v[1], v[2], x, y, squareSize, 0);
                dataProvider->setPipetteVal1(v[0]);
                dataProvider->setPipetteVal2(v[1]);
                dataProvider->setPipetteVal3(v[2]);
                return;
            }

            break;

        case (BT_LABIMAGE):
            if (LabBuffer) {
                float v[3];
                LabBuffer->getPipetteData(v[0], v[1], v[2], x, y, squareSize);
                dataProvider->setPipetteVal1(v[0]);
                dataProvider->setPipetteVal2(v[1]);
                dataProvider->setPipetteVal3(v[2]);
                return;
            }

            break;

        case (BT_SINGLEPLANE_FLOAT):
            if (singlePlaneBuffer.data != nullptr) {
                float v;
                singlePlaneBuffer.getPipetteData(v, x, y, squareSize, 0);
                dataProvider->setPipetteVal1(v);
                dataProvider->setPipetteVal2(-1.f);
                dataProvider->setPipetteVal3(-1.f);
                return;
            }
        }
    }

    if (dataProvider) {
        dataProvider->setPipetteVal1(-1.f);
        dataProvider->setPipetteVal2(-1.f);
        dataProvider->setPipetteVal3(-1.f);
    }
}

}
