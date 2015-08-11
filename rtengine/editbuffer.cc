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

#include "editbuffer.h"

namespace rtengine
{

EditBuffer::EditBuffer(::EditDataProvider *dataProvider) :
    objectMap(NULL), objectMap2(NULL), objectMode(OM_255), dataProvider(dataProvider),
    imgFloatBuffer(NULL), LabBuffer(NULL), singlePlaneBuffer(), ready(false) {}

EditBuffer::~EditBuffer()
{
    flush();
#ifndef NDEBUG
    imgFloatBuffer = (Imagefloat*)(0xbaadf00d);
#endif
#ifndef NDEBUG
    LabBuffer = (LabImage*)(0xbaadf00d);
#endif
}

void EditBuffer::createBuffer(int width, int height)
{
    //printf("Appel de createBuffer %d x %d\n", width, height);
    resize (width, height);
}

void EditBuffer::flush()
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

/* Upgrade or downgrade the objectModeType; we're assuming that objectMap has been allocated */
void EditBuffer::setObjectMode(ObjectMode newType)
{
    switch (newType) {
    case (OM_255):
        if (objectMap2) {
            objectMap2->unreference();
        }

        objectMode = OM_255;
        break;

    case (OM_65535):
        if (!objectMap2) {
            objectMap2 = Cairo::ImageSurface::create(Cairo::FORMAT_A8, objectMap->get_width(), objectMap->get_height());
        }

        objectMode = OM_65535;
        break;
    }
}

EditUniqueID EditBuffer::getEditID()
{
    if (dataProvider && dataProvider->getCurrSubscriber()) {
        return dataProvider->getCurrSubscriber()->getEditID();
    } else {
        return EUID_None;
    }
}

void EditBuffer::resize(int newWidth, int newHeight)
{
    resize(newWidth, newHeight, dataProvider ? dataProvider->getCurrSubscriber() : NULL);
}

// Resize buffers if they already exist
void EditBuffer::resize(int newWidth, int newHeight, EditSubscriber* newSubscriber)
{
    if (newSubscriber) {
        if (newSubscriber->getEditingType() == ET_OBJECTS) {
            if (objectMap && (objectMap->get_width() != newWidth || objectMap->get_height() != newHeight)) {
                objectMap.clear();
            }

            if (!objectMap && newWidth && newHeight) {
                objectMap = Cairo::ImageSurface::create(Cairo::FORMAT_A8, newWidth, newHeight);
            }

            if (objectMode == OM_65535) {
                if (objectMap2) {
                    if (objectMap2->get_width() != newWidth || objectMap2->get_height() != newHeight) {
                        objectMap2.clear();
                    }
                }

                if (!objectMap2 && newWidth && newHeight) {
                    objectMap2 = Cairo::ImageSurface::create(Cairo::FORMAT_A8, newWidth, newHeight);
                }
            }
            // OM_255 -> deleting objectMap2, if any
            else if (objectMap2) {
                objectMap2.clear();
            }

            // Should never happen!
            if (imgFloatBuffer) {
                delete imgFloatBuffer;
                imgFloatBuffer = NULL;
            }

            if (LabBuffer) {
                delete LabBuffer;
                LabBuffer = NULL;
            }

            if (singlePlaneBuffer.data) {
                singlePlaneBuffer.allocate(0, 0);
            }
        }

        if (newSubscriber->getEditingType() == ET_PIPETTE) {
            if (newSubscriber->getEditBufferType() == BT_IMAGEFLOAT) {
                if (!imgFloatBuffer) {
                    imgFloatBuffer = new Imagefloat(newWidth, newHeight);
                } else {
                    imgFloatBuffer->allocate(newWidth, newHeight);
                }
            } else if (imgFloatBuffer) {
                delete imgFloatBuffer;
                imgFloatBuffer = NULL;
            }

            if (newSubscriber->getEditBufferType() == BT_LABIMAGE) {
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

            if (newSubscriber->getEditBufferType() == BT_SINGLEPLANE_FLOAT) {
                singlePlaneBuffer.allocate(newWidth, newHeight);
            } else if (singlePlaneBuffer.data) {
                singlePlaneBuffer.allocate(0, 0);
            }

            // Should never happen!
            if (objectMap ) {
                objectMap.clear();
            }

            if (objectMap2) {
                objectMap2.clear();
            }
        }
    }

    ready = false;
}

bool EditBuffer::bufferCreated()
{
    EditSubscriber* subscriber;

    if (dataProvider && (subscriber = dataProvider->getCurrSubscriber())) {
        switch (subscriber->getEditingType()) {
        case ET_PIPETTE:
            switch (dataProvider->getCurrSubscriber()->getEditBufferType()) {
            case (BT_IMAGEFLOAT):
                return imgFloatBuffer != NULL;

            case (BT_LABIMAGE):
                return LabBuffer != NULL;

            case (BT_SINGLEPLANE_FLOAT):
                return singlePlaneBuffer.data != NULL;
            }

            break;

        case (ET_OBJECTS):
            return bool(objectMap);
        }
    }

    return false;
}

int EditBuffer::getObjectID(const Coord& location)
{
    int id = 0;

    if (objectMap && location.x > 0 && location.y > 0 && location.x < objectMap->get_width() && location.y < objectMap->get_height()) {
        id = (unsigned short)(*( objectMap->get_data() + location.y * objectMap->get_stride() + location.x ));

        if (objectMap2) {
            id |= (unsigned short)(*( objectMap->get_data() + location.y * objectMap->get_stride() + location.x )) << 8;
        }
    }

    return id - 1;
}

void EditBuffer::getPipetteData(float* v, int x, int y, int squareSize)
{
    if (ready && dataProvider && dataProvider->getCurrSubscriber()) {
        switch (dataProvider->getCurrSubscriber()->getEditBufferType()) {
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
