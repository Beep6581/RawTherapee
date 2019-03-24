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

EditSubscriber::EditSubscriber (EditType editType) :
    ID(EUID_None),
    editingType(editType),
    bufferType(BT_SINGLEPLANE_FLOAT),
    provider(nullptr),
    action(EditSubscriber::Action::NONE)
{}

void EditSubscriber::setEditProvider(EditDataProvider *provider)
{
    this->provider = provider;
}

void EditSubscriber::setEditID(EditUniqueID ID, BufferType buffType)
{
    this->ID = ID;
    bufferType = buffType;
}

const bool EditSubscriber::isCurrentSubscriber() const
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

const  EditUniqueID EditSubscriber::getEditID() const
{
    return ID;
}

const EditType EditSubscriber::getEditingType() const
{
    return editingType;
}

const  BufferType EditSubscriber::getPipetteBufferType() const
{
    return bufferType;
}

const bool EditSubscriber::isDragging() const
{
    return action == EditSubscriber::Action::DRAGGING;
}

const bool EditSubscriber::isPicking() const
{
    return action == EditSubscriber::Action::PICKING;
}

//--------------------------------------------------------------------------------------------------


EditDataProvider::EditDataProvider() :
    currSubscriber(nullptr),
    object(0),
    pipetteVal1(0.f),
    pipetteVal2(0.f),
    pipetteVal3(0.f),
    posScreen(-1, -1),
    posImage(-1, -1),
    deltaScreen(0, 0),
    deltaImage(0, 0),
    deltaPrevScreen(0, 0),
    deltaPrevImage(0, 0)
{}

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

int EditDataProvider::getObject()
{
    return object;
}

void EditDataProvider::setObject(int newObject)
{
    object = newObject;
}

float EditDataProvider::getPipetteVal1()
{
    return pipetteVal1;
}

float EditDataProvider::getPipetteVal2()
{
    return pipetteVal2;
}

float EditDataProvider::getPipetteVal3()
{
    return pipetteVal3;
}

void EditDataProvider::setPipetteVal1(float newVal)
{
    pipetteVal1 = newVal;
}

void EditDataProvider::setPipetteVal2(float newVal)
{
    pipetteVal2 = newVal;
}

void EditDataProvider::setPipetteVal3(float newVal)
{
    pipetteVal3 = newVal;
}

const CursorShape EditDataProvider::getCursor(int objectID) const
{
    if (currSubscriber) {
        currSubscriber->getCursor(objectID);
    }

    return CSHandOpen;
}

EditSubscriber* EditDataProvider::getCurrSubscriber() const
{
    return currSubscriber;
}

EditDataProvider* EditSubscriber::getEditProvider () {
    return provider;
}

const CursorShape EditSubscriber::getCursor (int objectID) const {
    return CSHandOpen;
}

const bool EditSubscriber::mouseOver (const int modifierKey) {
    return false;
}

bool EditSubscriber::button1Pressed (const int modifierKey) {
    return false;
}

bool EditSubscriber::button1Released () {
    return false;
}

bool EditSubscriber::button2Pressed (const int modifierKey) {
    return false;
}

bool EditSubscriber::button2Released () {
    return false;
}

bool EditSubscriber::button3Pressed (const int modifierKey) {
    return false;
}

bool EditSubscriber::button3Released () {
    return false;
}

bool EditSubscriber::drag1 (const int modifierKey) {
    return false;
}

bool EditSubscriber::drag2 (const int modifierKey) {
    return false;
}

bool EditSubscriber::drag3 (const int modifierKey) {
    return false;
}

bool EditSubscriber::pick1 (const bool picked) {
    return false;
}

bool EditSubscriber::pick2 (const bool picked) {
    return false;
}

bool EditSubscriber::pick3 (const bool picked) {
    return false;
}

const std::vector<Geometry*>& EditSubscriber::getVisibleGeometry () {
    return visibleGeometry;
}

const std::vector<Geometry*>& EditSubscriber::getMouseOverGeometry () {
    return mouseOverGeometry;
}

const int EditDataProvider::getPipetteRectSize () const {
    return 8; // TODO: make a GUI
}
