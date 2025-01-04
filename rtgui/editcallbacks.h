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

#include <vector>

#include "editid.h"
#include "cursormanager.h"
#include "rtengine/coord.h"

class Geometry;
class EditDataProvider;

/** @file
 * See editwidgets.h for documentation
 */

/// @brief Method for client tools needing Edit information
class EditSubscriber
{

public:

private:
    EditUniqueID ID; /// this will be used in improcfun to locate the data that has to be stored in the buffer; it must be unique in RT
    EditType editingType;
    BufferType bufferType;
    EditDataProvider *provider;

protected:
    std::vector<Geometry*> visibleGeometry;    /// displayed geometry
    std::vector<Geometry*> mouseOverGeometry;  /// mouseOver geometry, drawn in a hidden buffer
    enum class Action {
        NONE,      ///
        DRAGGING,  /// set action to this value in the buttonPressed event to start dragging and ask for drag event
        PICKING    /// set action to this value in the buttonPressed event whenever the user is picking something through a single click. In this case, the pickX events will be called INSTEAD of buttonXReleased !
    };

    Action action; /// object mode only, ignored in Pipette mode

public:
    explicit EditSubscriber (EditType editType);
    virtual ~EditSubscriber () = default;

    void               setEditProvider(EditDataProvider *provider);
    EditDataProvider*  getEditProvider () const;
    void               setEditID(EditUniqueID ID, BufferType buffType);
    bool               isCurrentSubscriber() const;
    virtual void       subscribe();
    virtual void       unsubscribe();
    virtual void       switchOffEditMode();           /// Occurs when the user want to stop the editing mode
    EditUniqueID       getEditID() const;
    EditType           getEditingType() const;
    BufferType         getPipetteBufferType() const;
    bool               isDragging() const;            /// Returns true if something is being dragged and drag events has to be sent (object mode only)
    bool               isPicking() const;             /// Returns true if something is being picked

    /** @brief Get the cursor to be displayed when above handles
    @param objectID object currently "hovered"
    @param xPos X cursor position in image space
    @param yPos Y cursor position in image space */
    virtual CursorShape getCursor (int objectID, int xPos, int yPos) const;

    /** @brief Triggered when the mouse is moving over an object
    This method is also triggered when the cursor is moving over the image in ET_PIPETTE mode
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool mouseOver (int modifierKey);

    /** @brief Triggered when mouse button 1 is pressed, together with the CTRL modifier key if the subscriber is of type ET_PIPETTE
    Once the key is pressed, RT will enter in drag1 mode on subsequent mouse movements
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button1Pressed (int modifierKey);

    /** @brief Triggered when mouse button 1 is released
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button1Released ();

    /** @brief Triggered when mouse button 2 is pressed (middle button)
    Once the key is pressed, RT will enter in drag2 mode on subsequent mouse movements
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button2Pressed (int modifierKey);

    /** @brief Triggered when mouse button 2 is released (middle button)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button2Released ();

    /** @brief Triggered when mouse button 3 is pressed (right button)
    Once the key is pressed, RT will enter in drag3 mode on subsequent mouse movements
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button3Pressed (int modifierKey);

    /** @brief Triggered when mouse button 3 is released (right button)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button3Released ();

    /** @brief Triggered when the user is moving while holding down mouse button 1
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool drag1 (int modifierKey);

    /** @brief Triggered when the user is moving while holding down mouse button 2
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool drag2 (int modifierKey);

    /** @brief Triggered when the user is moving while holding down mouse button 3
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool drag3 (int modifierKey);

    /** @brief Triggered when the user is releasing mouse button 1 while in action==ES_ACTION_PICKING mode
    No modifier key is provided, since having a different modifier key than on button press will set picked to false.
    @param picked True if the cursor is still above the the same object than on button pressed and with the same modifier keys.
                  If false, the user moved the cursor away or the modifier key is different, so the element is considered as NOT selected.
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool pick1 (bool picked);

    /** @brief Triggered when the user is releasing mouse button 2 while in action==ES_ACTION_PICKING mode
    @param picked True if the cursor is still above the the same object than on button pressed and with the same modifier keys.
                  If false, the user moved the cursor away or the modifier key is different, so the element is considered as NOT selected.
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool pick2 (bool picked);

    /** @brief Triggered when the user is releasing mouse button 3 while in action==ES_ACTION_PICKING mode
    @param picked True if the cursor is still above the the same object than on button pressed and with the same modifier keys.
                  If false, the user moved the cursor away or the modifier key is different, so the element is considered as NOT selected.
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool pick3 (bool picked);

    /** @brief Get the geometry to be shown to the user */
    const std::vector<Geometry*>& getVisibleGeometry ();

    /** @brief Get the geometry to be drawn in the "mouse over" channel, hidden from the user */
    const std::vector<Geometry*>& getMouseOverGeometry ();
};

/** @brief Class to handle the furniture of data to the subscribers.
 *
 * It is admitted that only one Subscriber can ask data at a time. If the Subscriber is of type ET_PIPETTE, it will have to
 * trigger the usual event so that the image will be reprocessed to compute the buffer of the current subscriber.
 */
class EditDataProvider
{

private:
    EditSubscriber *currSubscriber;
//    int object;            /// ET_OBJECTS mode: Object detected under the cursor, 0 otherwise; ET_PIPETTE mode: 1 if above the image, 0 otherwise
    float pipetteVal1;     /// Current pipette values
    float pipetteVal2;     /// Current pipette values; if bufferType==BT_SINGLEPLANE_FLOAT, will be set to 0
    float pipetteVal3;     /// Current pipette values; if bufferType==BT_SINGLEPLANE_FLOAT, will be set to 0

public:
    int object;            /// ET_OBJECTS mode: Object detected under the cursor, 0 otherwise; ET_PIPETTE mode: 1 if above the image, 0 otherwise

    rtengine::Coord posScreen;       /// Location of the mouse button press, in preview image space
    rtengine::Coord posImage;        /// Location of the mouse button press, in the full image space
    rtengine::Coord deltaScreen;     /// Delta relative to posScreen
    rtengine::Coord deltaImage;      /// Delta relative to posImage
    rtengine::Coord deltaPrevScreen; /// Delta relative to the previous mouse location, in preview image space
    rtengine::Coord deltaPrevImage;  /// Delta relative to the previous mouse location, in the full image space

    EditDataProvider();
    virtual ~EditDataProvider() = default;

    virtual void subscribe(EditSubscriber *subscriber);
    virtual void unsubscribe();         /// Occurs when the subscriber has been switched off first
    virtual void switchOffEditMode ();  /// Occurs when the user want to stop the editing mode
    int getObject() const;
    void setObject(int newObject);
    float getPipetteVal1() const;
    float getPipetteVal2() const;
    float getPipetteVal3() const;
    void setPipetteVal1(float newVal);
    void setPipetteVal2(float newVal);
    void setPipetteVal3(float newVal);
    virtual CursorShape getCursor(int objectID, int xPos, int yPos) const;
    int getPipetteRectSize () const;
    EditSubscriber* getCurrSubscriber() const;
    virtual void getImageSize (int &w, int&h) = 0;
    virtual void getPreviewCenterPos(int &x, int &y) = 0;
    virtual void getPreviewSize(int &w, int &h) = 0;
};
