/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Jean-Christophe FRISCH <natureh.510@gmail.com>
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

#include "editbuffer.h"
#include "editcallbacks.h"

ObjectMOBuffer::ObjectMOBuffer(EditDataProvider *dataProvider) : objectMap(nullptr), objectMode(OM_255), dataProvider(dataProvider) {}

ObjectMOBuffer::~ObjectMOBuffer()
{
    flush();
}


/* Upgrade or downgrade the objectModeType */
void ObjectMOBuffer::setObjectMode(ObjectMode newType)
{
    if (!objectMap) {
        objectMode = newType;
        return;
    }

    int w = objectMap->get_width ();
    int h = objectMap->get_height ();
    if (w && h) {
        switch (newType) {
        case (OM_255):
            if (objectMode==OM_65535) {
                objectMap.clear();
                objectMap = Cairo::ImageSurface::create(Cairo::FORMAT_A8, w, h);
            }
            break;

        case (OM_65535):
            if (objectMode==OM_255) {
                objectMap.clear();
                objectMap = Cairo::ImageSurface::create(Cairo::FORMAT_RGB16_565, w, h);
            }
            break;
        }
    }
    objectMode = newType;
}

void ObjectMOBuffer::flush()
{
    if (objectMap ) {
        objectMap.clear();
    }
}

EditSubscriber *ObjectMOBuffer::getEditSubscriber () {
    if (dataProvider) {
        return dataProvider->getCurrSubscriber();
    } else {
        return nullptr;
    }
}


// Resize buffers if they already exist
void ObjectMOBuffer::resize(int newWidth, int newHeight)
{
    if (!dataProvider) {
        return;
    }

    if (const auto currSubscriber = dataProvider->getCurrSubscriber ()) {
        if (currSubscriber->getEditingType() == ET_OBJECTS) {
            if (objectMap && (objectMap->get_width() != newWidth || objectMap->get_height() != newHeight)) {
                objectMap.clear();
            }

            if (!objectMap && newWidth>0 && newHeight>0) {
                objectMap = Cairo::ImageSurface::create(objectMode==OM_255?Cairo::FORMAT_A8:Cairo::FORMAT_RGB16_565, newWidth, newHeight);
            }

        } else {
            flush();
        }
    } else {
        flush();
    }
}

int ObjectMOBuffer::getObjectID(const rtengine::Coord& location)
{
    int id = 0;

    if (!objectMap || location.x < 0 || location.y < 0 || location.x >= objectMap->get_width() || location.y >= objectMap->get_height()) {
        return -1;
    }

    if (objectMode == OM_255) {
        // In OM_255 mode, size of pixel is 1 Byte (i.e. size of uint8_t)
        memcpy(&id, ( objectMap->get_data() + location.y * objectMap->get_stride() + sizeof(std::uint8_t) * location.x ), sizeof(std::uint8_t));
    } else {
        // In OM_65535 mode, size of pixel is 2 Bytes (i.e. size of uint16_t)
        memcpy(&id, ( objectMap->get_data() + location.y * objectMap->get_stride() + sizeof(std::uint16_t) * location.x ), sizeof(std::uint16_t));
    }

    return id - 1;
}

bool ObjectMOBuffer::bufferCreated()
{
    EditSubscriber* subscriber;

    if (dataProvider && (subscriber = dataProvider->getCurrSubscriber())) {
        return subscriber->getEditingType() == ET_OBJECTS ? bool(objectMap) : false;
    }

    return false;
}

