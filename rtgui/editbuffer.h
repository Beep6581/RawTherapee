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
#pragma once

#include "editid.h"
#include <cairomm/cairomm.h>

namespace rtengine {

struct Coord;

}

class EditDataProvider;
class EditSubscriber;

/** @file
 * See editwidgets.h for documentation
 */

class ObjectMOBuffer
{
private:

    // Used to draw the objects where the color correspond to the object's ID, in order to find the correct object when hovering
    Cairo::RefPtr<Cairo::ImageSurface> objectMap;
    ObjectMode objectMode;

protected:

    // To avoid duplicated information, we points to a EditDataProvider that contains the current EditSubscriber
    // instead of pointing to the EditSubscriber directly
    EditDataProvider* dataProvider;

    void createBuffer(int width, int height);
    void resize(int newWidth, int newHeight);
    void flush();
    EditSubscriber *getEditSubscriber ();

public:
    explicit ObjectMOBuffer (EditDataProvider *dataProvider);
    ~ObjectMOBuffer();

    EditDataProvider* getDataProvider ();
    void setObjectMode (ObjectMode newType);
    ObjectMode getObjectMode ();

    Cairo::RefPtr<Cairo::ImageSurface>& getObjectMap ();

    // return true if the buffer has been allocated
    bool bufferCreated();

    int getObjectID(const rtengine::Coord& location);
};

inline EditDataProvider* ObjectMOBuffer::getDataProvider () {
    return dataProvider;
}

inline ObjectMode ObjectMOBuffer::getObjectMode () {
    return objectMode;
}

inline Cairo::RefPtr<Cairo::ImageSurface>& ObjectMOBuffer::getObjectMap () {
    return objectMap;
}


