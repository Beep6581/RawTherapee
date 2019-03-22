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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "editcallbacks.h"

EditSubscriber::EditSubscriber (EditType editType) : ID(EUID_None), editingType(editType), bufferType(BT_SINGLEPLANE_FLOAT), provider(nullptr), action(ES_ACTION_NONE) {}

void EditSubscriber::setEditProvider(EditDataProvider *provider)
{
    this->provider = provider;
}

void EditSubscriber::setEditID(EditUniqueID ID, BufferType buffType)
{
    this->ID = ID;
    bufferType = buffType;
}

bool EditSubscriber::isCurrentSubscriber()
{
    //if (provider && provider->getCurrSubscriber())
    //  return provider->getCurrSubscriber()->getEditID() == ID;

    if (provider) {
        return provider->getCurrSubscriber() == this;
    }

    return false;
}

void EditSubscriber::subscribe()
{
    if (provider) {
        provider->subscribe(this);
    }
}

void EditSubscriber::unsubscribe()
{
    if (provider) {
        provider->unsubscribe();
    }
}

void EditSubscriber::switchOffEditMode()
{
    unsubscribe();
}

EditUniqueID EditSubscriber::getEditID()
{
    return ID;
}

EditType EditSubscriber::getEditingType()
{
    return editingType;
}

BufferType EditSubscriber::getPipetteBufferType()
{
    return bufferType;
}

bool EditSubscriber::isDragging()
{
    return action == ES_ACTION_DRAGGING;
}

bool EditSubscriber::isPicking()
{
    return action == ES_ACTION_PICKING;
}

//--------------------------------------------------------------------------------------------------


EditDataProvider::EditDataProvider() : currSubscriber(nullptr), object(0), posScreen(-1, -1), posImage(-1, -1),
    deltaScreen(0, 0), deltaImage(0, 0), deltaPrevScreen(0, 0), deltaPrevImage(0, 0)
{
    pipetteVal[0] = pipetteVal[1] = pipetteVal[2] = 0.f;
}

void EditDataProvider::subscribe(EditSubscriber *subscriber)
{
    if (currSubscriber) {
        currSubscriber->switchOffEditMode();
    }

    currSubscriber = subscriber;
}

void EditDataProvider::unsubscribe()
{
    currSubscriber = nullptr;
}

void EditDataProvider::switchOffEditMode()
{
    if (currSubscriber) {
        currSubscriber->switchOffEditMode ();
    }
}

CursorShape EditDataProvider::getCursor(int objectID)
{
    if (currSubscriber) {
        currSubscriber->getCursor(objectID);
    }

    return CSHandOpen;
}

EditSubscriber* EditDataProvider::getCurrSubscriber()
{
    return currSubscriber;
}

