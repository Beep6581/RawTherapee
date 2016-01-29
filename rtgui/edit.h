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
#include "../rtengine/coord.h"
#include "guiutils.h"
#include "options.h"

class EditDataProvider;
class EditSubscriber;

/** @file
 *
 * The Edit mechanism is designed to let tools (subscribers) communicate with the preview area (provider).
 * Subscribers will be tools that need to create some graphics in the preview area, to let the user interact
 * with it in a more user friendly way.
 *
 * Do not confuse with _local_ editing, which is another topic implemented in another class. The Edit feature
 * is also not supported in batch editing from the File Browser.
 *
 * Edit tool can be of 2 types: pipette editing and object editing.
 *
 * ## Pipette edition
 *
 * By using this class, a pipette mechanism can be handled on the preview.
 *
 * Each pipette Edit tool must have a unique ID, that will identify them, and which let the ImProcCoordinator
 * or other mechanism act as appropriated. They are all defined in rtgui/editid.h. A buffer type has to be given
 * too, to know which kind of buffer to allocate (see EditSubscriber::BufferType).
 *
 * Only the first mouse button can be used to manipulate the pipette on the Preview, that's why the developer has
 * to implement at least the following 4 methods:
 *    - mouseOver
 *    - button1Pressed
 *    - drag1
 *    - button1Released
 *
 * Actually, only curves does use this class, and everything is handled for curve implementor (as much as possible).
 * See the curve's class documentation to see how to implement the curve's pipette feature.
 *
 * ### Event handling
 *
 * The mouseOver method is called on each mouse movement, excepted when dragging a point. This method can then access
 * the pipetteVal array values, which contain the mean of the pixel read in the buffer, or -1 of the cursor is outside
 * of the image. In this case, EditDataProvider::object is also set to 0 (and 1 if over the image).
 *
 * When the user will click on the left mouse button while pressing the CTRL key, the button1Pressed will be called.
 * Setting "dragging" to true (or false) is not required for the pipette type editing.
 *
 * The drag1 method will be called on all subsequent mouse move. The pipetteVal[3] array will already be filled with
 * the mean of the read values under the cursor (actually a fixed square of 8px). If the BufferType is BT_SINGLEPLANE_FLOAT,
 * only the first array value will be filled.
 *
 * Then the button1Released will be called to stop the dragging.
 *
 * ## Object edition
 *
 * By using this class, objects can be drawn and manipulated on the preview.
 *
 * The developer has to handle the buttonPress, buttonRelease, drag and mouseOver method that he needs. There
 * are buttonPress, buttonRelease and drag methods dedicated to each mouse button, for better flexibility
 * (e.g.button2Press, button2Release, drag2 will handle event when mouse button 2 is used first). RT actually
 * does not handle multiple mouse button event (e.g. button1 + button2), only one at a time. The first button pressed
 * set the mechanism, all other combined button press are ignored.
 *
 * The developer also have to fill 2 display list with object of the Geometry subclass. Each geometric shape
 * _can_ be used in one or the other, or both list at a time.
 *
 * The first list (visibleGeometry) is used to be displayed on the preview. The developer will have to set their state
 * manually (see Geometry::State), but the display shape, line width and color can be handled automatically, or with
 * specific values. To be displayed, the F_VISIBLE flag has to be set through the setActive or setVisible methods.
 *
 * The second list (mouseOverGeometry) is used in a backbuffer, the color used to draw the shape being the id of the
 * mouseOverGeometry. As an example, you could use a circle line to represent the line to the user, but use another
 * Circle object, filled, to be used as mouseOver detection. The association between both shape (visible and mouseOver)
 * is handled by the developer. To be displayed on this backbuffer, the F_HOVERABLE flag has to be set through the
 * setActive or setHoverable methods.
 *
 *
 * ### Event handling
 *
 * RT will draw in the back buffer all mouseOverGeometry set by the developer once the Edit button is pressed, and handle
 * the events automatically.
 *
 * RT will call the mouseOver method on each mouse movement where no mouse button is pressed.
 *
 * On mouse button press over a mouseOverGeometry, it will call the button press method corresponding to the button
 * (e.g. button1Pressed for mouse button 1), with the modifier key as parameter. Any other mouse button pressed at
 * the same time will be ignored. It's up to the developer to decide whether it leads to a drag movement or not,
 * by setting the "dragging" boolean to true.
 *
 * In this case, RT will then sent drag1 event (to stay with our button 1 pressed example) on each mouse movement. It's up
 * to the developer of the tool to handle the dragging. The EditProvider class will help you in this by handling the actual
 * position in various coordinate system and ways.
 *
 * When the user will release the mouse button, RT will call the button1Release event (in our example). The developer have
 * then to set the "dragging" flag to false.
 *
 * Each of these methods have to returns a boolean value saying that the preview has to be refreshed or not (i.e. the displayed
 * geometry).
 *
 * ## Other general internal implementation notes
 *
 * When a tool is being constructed, unique IDs are affected to the EditSubscribers of the Pipette type.
 * Then the EditorPanel class will ask all ToolPanel to register the 'after' preview ImageArea object as data provider.
 * The Subscribers have now to provide a toggle button to click on to start the Edit listening. When toggling on, the Subscriber
 * register itself to the DataProvider, then an event is thrown through the standard ToolPanelListener::panelChanged
 * method to update the preview with new graphics to be displayed. If a previous Edit button was active, it will be deactivated
 * (the Edit buttons are mutually exclusive). For the Pipette type, a buffer will be created and has to be populated
 * by the developer in rtengine's pipeline. The unique pipette ID will be used to know where to fill the buffer, as each pipette
 * will need different data, corresponding to the state of the image right before the tool that needs pipette values. E.g for
 * the HSV tool, the Hue and Saturation and Value curves are applied on the current state of the image. That's why the pipette
 * of the H, S and V curve will share the same data of this "current state", otherwise the read value would be wrong.
 *
 * When the mouse button 1 is pressed while pressing the CTRL key, the button1Pressed method will be called.
 *
 * When the Edit process stops, the Subscriber is removed from the DataProvider, so buffers can be freed up.
 * A new ToolPanelListener::panelChanged event is also thrown to update the preview again, without the tool's
 * graphical objects. The Edit button is also toggled off (by the user or programmatically).
 *
 * It means that each Edit buttons toggled on will start an update of the preview which might or might not create
 * a new History entry, depending on the ProcEvent used.
 *
 */



class ObjectMOBuffer
{
private:

    // Used to draw the objects where the color correspond to the object's ID, in order to find the correct object when hovering
    Cairo::RefPtr<Cairo::ImageSurface> objectMap;
    // If more than 254 objects has to be handled, objectMap2 contains the "upper part" of the 16 bit int value. objectMap2 will be NULL otherwise.
    Cairo::RefPtr<Cairo::ImageSurface> objectMap2;
    ObjectMode objectMode;

protected:

    // To avoid duplicated information, we points to a EditDataProvider that contains the current EditSubscriber
    // instead of pointing to the EditSubscriber directly
    EditDataProvider* dataProvider;

    void createBuffer(int width, int height);
    void resize(int newWidth, int newHeight, EditSubscriber* newSubscriber);
    void flush();
    EditSubscriber *getEditSubscriber ();

public:
    ObjectMOBuffer(EditDataProvider *dataProvider);
    ~ObjectMOBuffer();

    EditDataProvider* getDataProvider()
    {
        return dataProvider;
    }
    void setObjectMode(ObjectMode newType);
    ObjectMode getObjectMode()
    {
        return objectMode;
    }

    Cairo::RefPtr<Cairo::ImageSurface> &getObjectMap ()
    {
        return objectMap;
    }
    Cairo::RefPtr<Cairo::ImageSurface> &getObjectMap2()
    {
        return objectMap2;
    }

    // return true if the buffer has been allocated
    bool bufferCreated();

    int getObjectID(const rtengine::Coord& location);
};


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
    /// Convert the image coords to the crop's canvas coords
    virtual void imageCoordToCropCanvas (int imgx, int imgy, int& phyx, int& phyy) = 0;
    /// Convert the image coords to the edit buffer coords
    virtual void imageCoordToCropBuffer (int imgx, int imgy, int& phyx, int& phyy) = 0;
    /// Convert a size value from the preview's scale to the image's scale
    virtual int scaleValueToImage (int value) = 0;
    /// Convert a size value from the preview's scale to the image's scale
    virtual float scaleValueToImage (float value) = 0;
    /// Convert a size value from the preview's scale to the image's scale
    virtual double scaleValueToImage (double value) = 0;
    /// Convert a size value from the image's scale to the preview's scale
    virtual int scaleValueToCanvas (int value) = 0;
    /// Convert a size value from the image's scale to the preview's scale
    virtual float scaleValueToCanvas (float value) = 0;
    /// Convert a size value from the image's scale to the preview's scale
    virtual double scaleValueToCanvas (double value) = 0;
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

class RGBAColor : public RGBColor
{
    double a;

public:
    RGBAColor () : RGBColor(0., 0., 0.), a(0.) {}
    explicit RGBAColor (double r, double g, double b, double a) : RGBColor(r, g, b), a(a) {}
    explicit RGBAColor (char r, char g, char b, char a) : RGBColor(r, g, b), a(double(a) / 255.) {}

    void setColor(double r, double g, double b, double a)
    {
        RGBColor::setColor(r, g, b);
        this->a = a;
    }

    void setColor(char r, char g, char b, char a)
    {
        RGBColor::setColor(r, g, b);
        this->a = double(a) / 255.;
    }

    double getA()
    {
        return a;
    }
};

/// @brief Displayable and MouseOver geometry base class
class Geometry
{
public:
    /// @brief Graphical state of the element
    enum State {
        NORMAL,     /// Default state
        ACTIVE,     /// Focused state
        PRELIGHT,   /// Hovered state
        DRAGGED,    /// When being dragged
        INSENSITIVE /// Displayed but insensitive
    };

    /// @brief Coordinate space and origin of the point
    enum Datum {
        IMAGE,          /// Image coordinate system with image's top left corner as origin
        CLICKED_POINT,  /// Screen coordinate system with clicked point as origin
        CURSOR          /// Screen coordinate system with actual cursor position as origin
    };
    enum Flags {
        F_VISIBLE     = 1 << 0, /// true if the geometry have to be drawn on the visible layer
        F_HOVERABLE   = 1 << 1, /// true if the geometry have to be drawn on the "mouse over" layer
        F_AUTO_COLOR  = 1 << 2, /// true if the color depend on the state value, not the color field above
    };

    /// @brief Key point of the image's rectangle that is used to locate the icon copy to the target point:
    enum DrivenPoint {
        DP_CENTERCENTER,
        DP_TOPLEFT,
        DP_TOPCENTER,
        DP_TOPRIGHT,
        DP_CENTERRIGHT,
        DP_BOTTOMRIGHT,
        DP_BOTTOMCENTER,
        DP_BOTTOMLEFT,
        DP_CENTERLEFT
    };

protected:
    RGBColor innerLineColor;
    RGBColor outerLineColor;
    short flags;

public:
    float innerLineWidth;  // ...outerLineWidth = innerLineWidth+2
    Datum datum;
    State state;  // set by the Subscriber

    Geometry () : innerLineColor(char(255), char(255), char(255)), outerLineColor(char(0), char(0), char(0)), flags(F_VISIBLE | F_HOVERABLE | F_AUTO_COLOR), innerLineWidth(1.5f), datum(IMAGE), state(NORMAL) {}
    virtual ~Geometry() {}

    void     setInnerLineColor     (double r, double g, double b)
    {
        innerLineColor.setColor(r, g, b);
        flags &= ~F_AUTO_COLOR;
    }
    void     setInnerLineColor     (char   r, char   g, char   b)
    {
        innerLineColor.setColor(r, g, b);
        flags &= ~F_AUTO_COLOR;
    }
    RGBColor getInnerLineColor     ();
    void     setOuterLineColor     (double r, double g, double b)
    {
        outerLineColor.setColor(r, g, b);
        flags &= ~F_AUTO_COLOR;
    }
    void     setOuterLineColor     (char   r, char   g, char   b)
    {
        outerLineColor.setColor(r, g, b);
        flags &= ~F_AUTO_COLOR;
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
            flags |= F_AUTO_COLOR;
        } else {
            flags &= ~F_AUTO_COLOR;
        }
    }
    bool     isVisible             ()
    {
        return flags & F_VISIBLE;
    }
    void     setVisible            (bool visible)
    {
        if (visible) {
            flags |= F_VISIBLE;
        } else {
            flags &= ~F_VISIBLE;
        }
    }
    bool     isHoverable           ()
    {
        return flags & F_HOVERABLE;
    }
    void     setHoverable          (bool visible)
    {
        if (visible) {
            flags |= F_HOVERABLE;
        } else {
            flags &= ~F_HOVERABLE;
        }
    }

    // setActive will enable/disable the visible and hoverable flags in one shot!
    void     setActive             (bool active)
    {
        if (active) {
            flags |= (F_VISIBLE | F_HOVERABLE);
        } else {
            flags &= ~(F_VISIBLE | F_HOVERABLE);
        }
    }

    virtual void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *parent, EditCoordSystem &coordSystem) = 0;
    virtual void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *parent, EditCoordSystem &coordSystem) = 0;
    virtual void drawToMOChannel   (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem) = 0;
};

class Circle : public Geometry
{
public:
    rtengine::Coord center;
    int radius;
    bool filled;
    bool radiusInImageSpace; /// If true, the radius depend on the image scale; if false, it is a fixed 'screen' size

    Circle () : center(100, 100), radius(10), filled(false), radiusInImageSpace(false) {}
    Circle (rtengine::Coord &center, int radius, bool filled = false, bool radiusInImageSpace = false) : center(center), radius(radius), filled(filled), radiusInImageSpace(radiusInImageSpace) {}
    Circle (int centerX, int centerY, int radius, bool filled = false, bool radiusInImageSpace = false) : center(centerX, centerY), radius(radius), filled(filled), radiusInImageSpace(radiusInImageSpace) {}

    void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
    void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
    void drawToMOChannel   (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
};

class Line : public Geometry
{
public:
    rtengine::Coord begin;
    rtengine::Coord end;

    Line () : begin(10, 10), end(100, 100) {}
    Line (rtengine::Coord &begin, rtengine::Coord &end) : begin(begin), end(end) {}
    Line (int beginX, int beginY, int endX, int endY) : begin(beginX, beginY), end(endX, endY) {}

    void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
    void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
    void drawToMOChannel   (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
};

class Polyline : public Geometry
{
public:
    std::vector<rtengine::Coord> points;
    bool filled;

    Polyline() : filled(false) {}

    void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
    void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
    void drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
};

class Rectangle : public Geometry
{
public:
    rtengine::Coord topLeft;
    rtengine::Coord bottomRight;
    bool filled;

    Rectangle() : topLeft(0, 0), bottomRight(10, 10), filled(false) {}

    void setXYWH(int left, int top, int width, int height);
    void setXYXY(int left, int top, int right, int bottom);
    void setXYWH(rtengine::Coord topLeft, rtengine::Coord widthHeight);
    void setXYXY(rtengine::Coord topLeft, rtengine::Coord bottomRight);
    void drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
    void drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
    void drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, ObjectMOBuffer *pipetteBuffer, EditCoordSystem &coordSystem);
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
    std::vector<Geometry*> visibleGeometry;    /// displayed geometry
    std::vector<Geometry*> mouseOverGeometry;  /// mouseOver geometry, drawn in a hidden buffer
    bool dragging;  /// in object mode, set this to true in buttonPressed events to start dragging and ask for drag event (ignored in pipette mode)

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
    BufferType        getPipetteBufferType();
    bool              isDragging();            /// Returns true if something is being dragged and drag events has to be sent (object mode only)

    /** @brief Get the cursor to be displayed when above handles
    @param objectID object currently "hovered" */
    virtual CursorShape getCursor(int objectID)
    {
        return CSOpenHand;
    }

    /** @brief Triggered when the mouse is moving over an object
    This method is also triggered when the cursor is moving over the image in ET_PIPETTE mode
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool mouseOver(int modifierKey)
    {
        return false;
    }

    /** @brief Triggered when mouse button 1 is pressed, together with the CTRL modifier key if the subscriber is of type ET_PIPETTE
    Once the key is pressed, RT will enter in drag1 mode on subsequent mouse movements
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

    /** @brief Triggered when mouse button 2 is pressed (middle button)
    Once the key is pressed, RT will enter in drag2 mode on subsequent mouse movements
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button2Pressed(int modifierKey)
    {
        return false;
    }

    /** @brief Triggered when mouse button 2 is released (middle button)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button2Released()
    {
        return false;
    }

    /** @brief Triggered when mouse button 3 is pressed (right button)
    Once the key is pressed, RT will enter in drag3 mode on subsequent mouse movements
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button3Pressed(int modifierKey)
    {
        return false;
    }

    /** @brief Triggered when mouse button 3 is released (right button)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool button3Released()
    {
        return false;
    }

    /** @brief Triggered when the user is moving while holding down mouse button 1
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool drag1(int modifierKey)
    {
        return false;
    }

    /** @brief Triggered when the user is moving while holding down mouse button 2
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool drag2(int modifierKey)
    {
        return false;
    }

    /** @brief Triggered when the user is moving while holding down mouse button 3
    @param modifierKey Gtk's event modifier key (GDK_CONTROL_MASK | GDK_SHIFT_MASK | ...)
    @return true if the preview has to be redrawn, false otherwise */
    virtual bool drag3(int modifierKey)
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

    rtengine::Coord posScreen;       /// Location of the mouse button press, in preview image space
    rtengine::Coord posImage;        /// Location of the mouse button press, in the full image space
    rtengine::Coord deltaScreen;     /// Delta relative to posScreen
    rtengine::Coord deltaImage;      /// Delta relative to posImage
    rtengine::Coord deltaPrevScreen; /// Delta relative to the previous mouse location, in preview image space
    rtengine::Coord deltaPrevImage;  /// Delta relative to the previous mouse location, in the full image space

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
