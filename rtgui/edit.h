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
#ifndef _EDIT_H_
#define _EDIT_H_

#include <gtkmm.h>
#include "../rtengine/imagefloat.h"
#include "editid.h"
#include "cursormanager.h"
#include "../rtengine/rt_math.h"
#include "coord.h"

class EditDataProvider;

namespace rtengine
{
class EditBuffer;
}

/** @file
 *
 * The Edit mechanism is designed to let tools (subscribers) communicate with the preview area (provider).
 * Subscribers will be tools that need to create some graphics in the preview area, to let the user interact
 * with it in a more user friendly way.
 *
 * Do not confuse with a _local_ editing, which is another topic implemented in another class. The Edit feature
 * is also not supported in batch editing from the File Browser.
 *
 * Each Edit tool must have a unique ID, that will identify them, and which let the ImProcCoordinator or other
 * mechanism act as appropriated. They are all defined in rtgui/editid.h.
 *
 * ## How does it works?
 *
 * When a tool is being constructed, unique IDs are affected to the EditSubscribers. Then the EditorPanel class
 * will ask all ToolPanel to register the 'after' preview ImageArea object as data provider. The Subscribers
 * have now to provide a toggle button to click on to start the Edit listening. When toggling on, the Subscriber
 * register itself to the DataProvider, then an event is thrown through the standard ToolPanelListener::panelChanged
 * method to update the preview with new graphics to be displayed, and eventually buffers to create and populate.
 *
 * When the Edit process stops, the Subscriber is removed from the DataProvider, so buffer can be freed up.
 * A new ToolPanelListener::panelChanged event is also thrown to update the preview again, without the tool's
 * graphical objects. The Edit button is also toggled off (by the user or programmatically).
 *
 * It means that each Edit buttons toggled on will start an update of the preview which might or might not create
 * a new History entry, depending on the event.
 *
 */


/** @brief Coordinate system where the widgets will be drawn
 *
 * The EditCoordSystem is used to define a screen and an image coordinate system.
 */
class EditCoordSystem
{
public:
    virtual ~EditCoordSystem() {}

    /// Convert the widget's DrawingArea (i.e. preview area) coords to the edit buffer coords
    virtual void screenCoordToCropBuffer (int phyx, int phyy, int& cropx, int& cropy) = 0;
    /// Convert the widget's DrawingArea (i.e. preview area) coords to the full image coords
    virtual void screenCoordToImage (int phyx, int phyy, int& imgx, int& imgy) = 0;
    /// Convert the image coords to the widget's DrawingArea (i.e. preview area) coords
    virtual void imageCoordToScreen (int imgx, int imgy, int& phyx, int& phyy) = 0;
    /// Convert the image coords to the edit buffer coords
    virtual void imageCoordToCropBuffer (int imgx, int imgy, int& phyx, int& phyy) = 0;
    /// Convert a size value from the preview's scale to the image's scale
    virtual int scaleValueToImage (int value) = 0;
    /// Convert a size value from the preview's scale to the image's scale
    virtual float scaleValueToImage (float value) = 0;
    /// Convert a size value from the preview's scale to the image's scale
    virtual double scaleValueToImage (double value) = 0;
    /// Convert a size value from the image's scale to the preview's scale
    virtual int scaleValueToScreen (int value) = 0;
    /// Convert a size value from the image's scale to the preview's scale
    virtual float scaleValueToScreen (float value) = 0;
    /// Convert a size value from the image's scale to the preview's scale
    virtual double scaleValueToScreen (double value) = 0;
};

class RGBColor
{
    double r;
    double g;
    double b;

public:
    RGBColor () : r(0.), g(0.), b(0.) {}
    explicit RGBColor (double r, double g, double b) : r(r), g(g), b(b) {}
    explicit RGBColor (char r, char g, char b) : r(double(r) / 255.), g(double(g) / 255.), b(double(b) / 255.) {}

    void setColor(double r, double g, double b)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    void setColor(char r, char g, char b)
    {
        this->r = double(r) / 255.;
        this->g = double(g) / 255.;
        this->b = double(b) / 255.;
    }

    double getR()
    {
        return r;
    }
    double getG()
    {
        return g;
    }
    double getB()
    {
        return b;
    }
};

class Geometry
{
public:
    enum State {
        NORMAL,
        PRELIGHT,
        DRAGGED
    };
    enum Datum {
        IMAGE,
        CLICKED_POINT,
        CURSOR
    };
    enum Flags {
        ACTIVE      = 1 << 0, // true if the geometry is active and have to be drawn
        AUTO_COLOR  = 1 << 1, // true if the color depend on the state value, not the color field above
    };

protected:
    RGBColor innerLineColor;
    RGBColor outerLineColor;
    short flags;

public:
    float innerLineWidth;  // ...outerLineWidth = innerLineWidth+2
    Datum datum;
    State state;  // set by the Subscriber

    Geometry () : innerLineColor(char(255), char(255), char(255)), outerLineColor(char(0), char(0), char(0)), flags(ACTIVE | AUTO_COLOR), innerLineWidth(1.f), datum(IMAGE), state(NORMAL) {}
    virtual ~Geometry() {}

    void     setInnerLineColor     (double r, double g, double b)
    {
        innerLineColor.setColor(r, g, b);
        flags &= ~AUTO_COLOR;
    }
    void     setInnerLineColor     (char   r, char   g, char   b)
    {
        innerLineColor.setColor(r, g, b);
        flags &= ~AUTO_COLOR;
    }
    RGBColor getInnerLineColor     ();
    void     setOuterLineColor     (double r, double g, double b)
    {
        outerLineColor.setColor(r, g, b);
        flags &= ~AUTO_COLOR;
    }
    void     setOuterLineColor     (char   r, char   g, char   b)
    {
        outerLineColor.setColor(r, g, b);
        flags &= ~AUTO_COLOR;
    }
    RGBColor getOuterLineColor     ();
    double   getOuterLineWidth     ()
    {
        return double(innerLineWidth) + 2.;
    }
    double   getMouseOverLineWidth ()
    {
        return getOuterLineWidth() + 2.;
    }
    void     setAutoColor          (bool aColor)
    {
        if (aColor) {
            flags |= AUTO_COLOR;
        } else {
            flags &= ~AUTO_COLOR;
        }
    }
    void     setActive             (bool active)
    {
        if (active) {
            flags |= ACTIVE;
        } else {
            flags &= ~ACTIVE;
        }
    }

    virtual void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *parent, EditCoordSystem &coordSystem) = 0;
    virtual void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *parent, EditCoordSystem &coordSystem) = 0;
    virtual void drawToMOChannel   (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) = 0;
};

class Circle : public Geometry
{
public:
    Coord center;
    int radius;
    bool filled;
    bool radiusInImageSpace; /// If true, the radius depend on the image scale; if false, it is a fixed 'screen' size

    Circle () : center(100, 100), radius(10), filled(false), radiusInImageSpace(false) {}
    Circle (Coord &center, int radius, bool filled = false, bool radiusInImageSpace = false) : center(center), radius(radius), filled(filled), radiusInImageSpace(radiusInImageSpace) {}
    Circle (int centerX, int centerY, int radius, bool filled = false, bool radiusInImageSpace = false) : center(centerX, centerY), radius(radius), filled(filled), radiusInImageSpace(radiusInImageSpace) {}

    void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
    void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
    void drawToMOChannel   (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
};

class Line : public Geometry
{
public:
    Coord begin;
    Coord end;

    Line () : begin(10, 10), end(100, 100) {}
    Line (Coord &begin, Coord &end) : begin(begin), end(end) {}
    Line (int beginX, int beginY, int endX, int endY) : begin(beginX, beginY), end(endX, endY) {}

    void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
    void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
    void drawToMOChannel   (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
};

class Polyline : public Geometry
{
public:
    std::vector<Coord> points;
    bool filled;

    Polyline() : filled(false) {}

    void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
    void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
    void drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
};

class Rectangle : public Geometry
{
public:
    Coord topLeft;
    Coord bottomRight;
    bool filled;

    Rectangle() : topLeft(0, 0), bottomRight(10, 10), filled(false) {}

    void setXYWH(int left, int top, int width, int height);
    void setXYXY(int left, int top, int right, int bottom);
    void setXYWH(Coord topLeft, Coord widthHeight);
    void setXYXY(Coord topLeft, Coord bottomRight);
    void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
    void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
    void drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem);
};

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
    std::vector<Geometry*> visibleGeometry;
    std::vector<Geometry*> mouseOverGeometry;

public:
    EditSubscriber (EditType editType);
    virtual ~EditSubscriber () {}

    void              setEditProvider(EditDataProvider *provider);
    EditDataProvider* getEditProvider()
    {
        return provider;
    }
    void              setEditID(EditUniqueID ID, BufferType buffType);
    bool              isCurrentSubscriber();
    virtual void      subscribe();
    virtual void      unsubscribe();
    virtual void      switchOffEditMode ();    /// Occurs when the user want to stop the editing mode
    EditUniqueID      getEditID();
    EditType          getEditingType();
    BufferType        getEditBufferType();

    /** @brief Get the cursor to be displayed when above handles
    @param objectID object currently "hovered" */
    virtual CursorShape getCursor(int objectID)
    {
        return CSOpenHand;
    }

    /** @brief Triggered when the mouse is moving over an object
    This method is also triggered when the cursor is moving over the image in ET_PIPETTE mode
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @param editBuffer buffer to get the pipette values and the from
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool mouseOver(int modifierKey)
    {
        return false;
    }

    /** @brief Triggered when mouse button 1 is pressed, together with the CTRL modifier key if the subscriber is of type ET_PIPETTE
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button1Pressed(int modifierKey)
    {
        return false;
    }

    /** @brief Triggered when mouse button 1 is released
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button1Released()
    {
        return false;
    }

    /** @brief Triggered when the user is moving while holding down mouse button 1
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool drag(int modifierKey)
    {
        return false;
    }

    /** @brief Get the geometry to be shown to the user */
    const std::vector<Geometry*> & getVisibleGeometry()
    {
        return visibleGeometry;
    }

    /** @brief Get the geometry to be drawn in the "mouse over" channel, hidden from the user */
    const std::vector<Geometry*> & getMouseOverGeometry()
    {
        return mouseOverGeometry;
    }
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

public:
    int object;            /// ET_OBJECTS mode: Object detected under the cursor, 0 otherwise; ET_PIPETTE mode: 1 if above the image, 0 otherwise
    float pipetteVal[3];   /// Current pipette values; if bufferType==BT_SINGLEPLANE_FLOAT, #2 & #3 will be set to 0

    Coord posScreen;       /// Location of the mouse button press, in preview image space
    Coord posImage;        /// Location of the mouse button press, in the full image space
    Coord deltaScreen;     /// Delta relative to posScreen
    Coord deltaImage;      /// Delta relative to posImage
    Coord deltaPrevScreen; /// Delta relative to the previous mouse location, in preview image space
    Coord deltaPrevImage;  /// Delta relative to the previous mouse location, in the full image space

    EditDataProvider();
    virtual ~EditDataProvider() {}

    virtual void        subscribe(EditSubscriber *subscriber);
    virtual void        unsubscribe();         /// Occurs when the subscriber has been switched off first
    virtual void        switchOffEditMode ();  /// Occurs when the user want to stop the editing mode
    virtual CursorShape getCursor(int objectID);
    int                 getPipetteRectSize()
    {
        return 8;    // TODO: make a GUI
    }
    EditSubscriber*     getCurrSubscriber();
    virtual void        getImageSize (int &w, int&h) = 0;
};

#endif
