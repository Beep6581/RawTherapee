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
#ifndef __GUI_UTILS_
#define __GUI_UTILS_

#include <gtkmm.h>
#include <cairomm/cairomm.h>
#include "../rtengine/rtengine.h"
#include "rtimage.h"
#include <sstream>
#include <iostream>

Glib::ustring escapeHtmlChars(const Glib::ustring &src);
bool removeIfThere (Gtk::Container* cont, Gtk::Widget* w, bool increference = true);
void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh);
Glib::ustring removeExtension (const Glib::ustring& filename);
Glib::ustring getExtension (const Glib::ustring& filename);
bool confirmOverwrite (Gtk::Window& parent, const std::string& filename);
void writeFailed (Gtk::Window& parent, const std::string& filename);
void drawCrop (Cairo::RefPtr<Cairo::Context> cr, int imx, int imy, int imw, int imh, int startx, int starty, double scale, const rtengine::procparams::CropParams& cparams, bool drawGuide = true, bool useBgColor = true, bool fullImageVisible = true);
gboolean acquireGUI(void* data);
void setExpandAlignProperties(Gtk::Widget *widget, bool hExpand, bool vExpand, enum Gtk::Align hAlign, enum Gtk::Align vAlign);

guint add_idle (GSourceFunc function, gpointer data);

// TODO: The documentation says gdk_threads_enter and gdk_threads_leave should be replaced
// by g_main_context_invoke(), g_idle_add() and related functions, but this will require more extensive changes.
// We silence those warnings until then so that we notice the others.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

/**
 * @brief Lock GTK for critical section.
 *
 * Will unlock on destruction. To use:
 *
 *   <code>
 *     {
 *       GThreadLock lock;
 *       // critical code
 *     }
 *   </code>
 */
class GThreadLock
{
public:
    GThreadLock()
    {
        gdk_threads_enter();
    }
    ~GThreadLock()
    {
        gdk_threads_leave();
    }
};

/**
 * @brief Unlock GTK critical section.
 *
 * Will relock on destruction.
 */
class GThreadUnLock
{
public:
    GThreadUnLock()
    {
        gdk_threads_leave();
    }
    ~GThreadUnLock()
    {
        gdk_threads_enter();
    }
};

#pragma GCC diagnostic pop

class ConnectionBlocker
{
public:
    ConnectionBlocker (sigc::connection& connection) : connection (connection)
    {
        wasBlocked = connection.block();
    }
    ~ConnectionBlocker ()
    {
        connection.block(wasBlocked);
    }
private:
    sigc::connection& connection;
    bool wasBlocked;
};

/**
 * @brief Glue box to control visibility of the MyExpender's content ; also handle the frame around it
 */
class ExpanderBox: public Gtk::EventBox
{
private:
    Gtk::Container *pC;

public:
    ExpanderBox( Gtk::Container *p);
    ~ExpanderBox( )
    {
        delete pC;
    }

    void setLevel(int level);
    void updateStyle();

    void show() {}
    void show_all();
    void hide() {}
    void set_visible(bool isVisible = true) {}

    void showBox();
    void hideBox();

    void on_style_updated ();
//  bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr);
};

/**
 * @brief A custom Expander class, that can handle widgets in the title bar
 *
 * Custom made expander for responsive widgets in the header. It also handle a "enabled/disabled" property that display
 * a different arrow depending on this boolean value.
 *
 * Warning: once you've instantiated this class with a text label or a widget label, you won't be able to revert to the other solution.
 */
class MyExpander : public Gtk::VBox
{
public:
    typedef sigc::signal<void> type_signal_enabled_toggled;
private:
    type_signal_enabled_toggled message;
    static Glib::RefPtr<Gdk::Pixbuf> inconsistentPBuf; /// "inconsistent" image, displayed when useEnabled is true ; in this case, nothing will tell that an expander is opened/closed
    static Glib::RefPtr<Gdk::Pixbuf> enabledPBuf;      ///      "enabled" image, displayed when useEnabled is true ; in this case, nothing will tell that an expander is opened/closed
    static Glib::RefPtr<Gdk::Pixbuf> disabledPBuf;     ///     "disabled" image, displayed when useEnabled is true ; in this case, nothing will tell that an expander is opened/closed
    static Glib::RefPtr<Gdk::Pixbuf> openedPBuf;       ///       "opened" image, displayed when useEnabled is false
    static Glib::RefPtr<Gdk::Pixbuf> closedPBuf;       ///       "closed" image, displayed when useEnabled is false
    bool enabled;               /// Enabled feature (default to true)
    bool inconsistent;          /// True if the enabled button is inconsistent
    Gtk::EventBox *titleEvBox;  /// EventBox of the title, to get a connector from it
    Gtk::HBox *headerHBox;
    bool flushEvent;            /// Flag to control the weird event mechanism of Gtk (please prove me wrong!)
    ExpanderBox* expBox;        /// Frame that includes the child and control its visibility
    Gtk::EventBox *imageEvBox;  /// Enable/Disable or Open/Close arrow event box

    /// Triggered on opened/closed event
    bool on_toggle(GdkEventButton* event);
    /// Triggered on enabled/disabled change -> will emit a toggle event to the connected objects
    bool on_enabled_change(GdkEventButton* event);
    /// Used to handle the colored background for the whole Title
    bool on_enter_leave_title (GdkEventCrossing* event);
    /// Used to handle the colored background for the Enable button
    bool on_enter_leave_enable (GdkEventCrossing* event);
    // The part below can probably be removed from here and the CSS file.
    /// Update the style of this widget, depending in the "slim" option
    void updateStyle();

    void on_style_updated ()
    {
        updateStyle();
    }



protected:
    Gtk::Container* child;      /// Gtk::Contained to display below the expander's title
    Gtk::Widget* headerWidget;  /// Widget to display in the header, next to the arrow image ; can be NULL if the "string" version of the ctor has been used
    Gtk::Image* statusImage;    /// Image to display the opened/closed status (if useEnabled is false) of the enabled/disabled status (if useEnabled is true)
    Gtk::Label* label;          /// Text to display in the header, next to the arrow image ; can be NULL if the "widget" version of the ctor has been used
    bool useEnabled;            /// Set whether to handle an enabled/disabled feature and display the appropriate images

public:

    /** @brief Create a custom expander with a simple header made of a label.
     * @param useEnabled Set whether to handle an enabled/disabled toggle button and display the appropriate image
     * @param titleLabel A string to display in the header. Warning: you won't be able to switch to a widget label.
     */
    MyExpander(bool useEnabled, Glib::ustring titleLabel);

    /** Create a custom expander with a a custom - and responsive - widget
     * @param useEnabled Set whether to handle an enabled/disabled toggle button and display the appropriate image
     * @param titleWidget A widget to display in the header. Warning: you won't be able to switch to a string label.
     */
    MyExpander(bool useEnabled, Gtk::Widget* titleWidget);

    /// Initialize the class by loading the images
    static void init();

    Glib::SignalProxy1< bool, GdkEventButton* > signal_button_release_event()
    {
        return titleEvBox->signal_button_release_event();
    };
    type_signal_enabled_toggled signal_enabled_toggled();

    /// Set the nesting level of the Expander to adapt its style accordingly
    void setLevel(int level);

    /// Set a new label string. If it has been instantiated with a Gtk::Widget, this method will do nothing
    void setLabel (Glib::ustring newLabel);
    /// Set a new label string. If it has been instantiated with a Gtk::Widget, this method will do nothing
    void setLabel (Gtk::Widget *newWidget);

    /// Get whether the enabled option is set (to true or false) or unset (i.e. undefined)
    bool get_inconsistent();
    /// Set whether the enabled option is set (to true or false) or unset (i.e. undefined)
    void set_inconsistent(bool isInconsistent);

    /// Get whether the enabled button is used or not
    bool getUseEnabled();
    /// Get whether the enabled button is on or off
    bool getEnabled();
    /// If not inconsistent, set the enabled button to true or false and emit the message if the state is different
    /// If inconsistent, set the internal value to true or false, but do not update the image and do not emit the message
    void setEnabled(bool isEnabled);
    /// Adds a Tooltip to the Enabled button, if it exist ; do nothing otherwise
    void setEnabledTooltipMarkup(Glib::ustring tooltipMarkup);
    void setEnabledTooltipText(Glib::ustring tooltipText);

    /// Get the header widget. It'll send back the Gtk::Label* if it has been instantiated with a simple text
    Gtk::Widget* getLabelWidget() const
    {
        return headerWidget ? headerWidget : label;
    }

    /// Get the widget shown/hidden by the expander
    Gtk::Container* getChild();

    /// Set the collapsed/expanded state of the expander
    void set_expanded( bool expanded );

    /// Get the collapsed/expanded state of the expander
    bool get_expanded();

    /// Add a Gtk::Container for the content of the expander
    /// Warning: do not manually Show/Hide the widget, because this parameter is handled by the click on the Expander's title
    void add  (Gtk::Container& widget);
};


/**
 * @brief subclass of Gtk::ScrolledWindow in order to handle the scrollwheel
 */
class MyScrolledWindow : public Gtk::ScrolledWindow
{

    bool on_scroll_event (GdkEventScroll* event);
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const;

public:
    MyScrolledWindow();
};

/**
 * @brief subclass of Gtk::ComboBox in order to handle the scrollwheel
 */
class MyComboBox : public Gtk::ComboBox
{

    bool on_scroll_event (GdkEventScroll* event);
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const;
    void get_preferred_width_for_height_vfunc (int width, int &minimum_width, int &natural_width) const;

public:
    MyComboBox ();
};

/**
 * @brief subclass of Gtk::ComboBoxText in order to handle the scrollwheel
 */
class MyComboBoxText : public Gtk::ComboBoxText
{

    bool on_scroll_event (GdkEventScroll* event);
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const;
    void get_preferred_width_for_height_vfunc (int width, int &minimum_width, int &natural_width) const;

public:
    MyComboBoxText ();
};

/**
 * @brief subclass of Gtk::SpinButton in order to handle the scrollwheel
 */
class MySpinButton : public Gtk::SpinButton
{

protected:
    bool on_scroll_event (GdkEventScroll* event);
    bool on_key_press_event (GdkEventKey* event);

public:
    MySpinButton ();
    void updateSize();
};

/**
 * @brief subclass of Gtk::HScale in order to handle the scrollwheel
 */
class MyHScale : public Gtk::HScale
{

    bool on_scroll_event (GdkEventScroll* event);
    bool on_key_press_event (GdkEventKey* event);
};

/**
 * @brief subclass of Gtk::FileChooserButton in order to handle the scrollwheel
 */
class MyFileChooserButton : public Gtk::FileChooserButton
{

protected:
    bool on_scroll_event (GdkEventScroll* event);
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const;
    void get_preferred_width_for_height_vfunc (int width, int &minimum_width, int &natural_width) const;

public:
    MyFileChooserButton (const Glib::ustring& title, Gtk::FileChooserAction action = Gtk::FILE_CHOOSER_ACTION_OPEN);
};

/**
 * A class which maintains the last folder for a FileChooserDialog or Button by
 * caching it in a a variable (which can be persisted externally).
 * Each time the user selects a file or folder, the provided variable is updated
 * with the associated folder. The path found in the variable is set in the
 * dialog instance at constructions time of this object.
 */
class FileChooserLastFolderPersister: public Glib::Object
{
public:

    /**
     * Installs this persister on the provided GtkFileChooser instance and
     * applies the current folder found in @p folderVariable for the dialog.
     *
     * @param chooser file chooser to maintain
     * @param folderVariable variable storage to use for this dialog
     */
    FileChooserLastFolderPersister(Gtk::FileChooser* chooser, Glib::ustring& folderVariable);

    virtual ~FileChooserLastFolderPersister();

private:

    /**
     * Signal handler for the GtkFileChooser selection action.
     */
    void selectionChanged();

    Gtk::FileChooser* chooser;
    Glib::ustring& folderVariable;
    sigc::connection selectionChangedConnetion;

};

typedef enum RTUpdatePolicy {
    RTUP_STATIC,
    RTUP_DYNAMIC
} eUpdatePolicy;

typedef enum RTOrientation {
    RTO_Left2Right,
    RTO_Bottom2Top,
    RTO_Right2Left,
    RTO_Top2Bottom
} eRTOrientation;

enum TOITypes {
    TOI_TEXT,
    TOI_ICON
};

typedef enum RTNav {
    NAV_NONE,
    NAV_NEXT,
    NAV_PREVIOUS
} eRTNav;

/**
 * @brief Handle the switch between text and image to be displayed in the HBox (to be used in a button/toolpanel)
 */
class TextOrIcon : public Gtk::HBox
{

protected:
    Gtk::Image* imgIcon;
    Gtk::Label* label;
    Glib::ustring filename;
    Glib::ustring labelText;
    Glib::ustring tooltipText;

public:
    TextOrIcon (Glib::ustring filename, Glib::ustring labelTx, Glib::ustring tooltipTx, TOITypes type);
    ~TextOrIcon ();

    void switchTo(TOITypes type);
};

class MyImageMenuItem : public Gtk::MenuItem
{
private:
    Gtk::Grid *box;
    RTImage *image;
    Gtk::Label *label;

public:
    MyImageMenuItem (Glib::ustring label, Glib::ustring imageFileName);
    const RTImage *getImage () const;
    const Gtk::Label* getLabel () const;
};

/**
 * @brief Define a gradient milestone
 */
class GradientMilestone
{
public:
    double position;
    double r;
    double g;
    double b;
    double a;

    GradientMilestone(double _p = 0., double _r = 0., double _g = 0., double _b = 0., double _a = 0.)
    {
        position = _p;
        r = _r;
        g = _g;
        b = _b;
        a = _a;
    }
};

/**
 * @brief Handle point coordinates
 */
template <class T>
class Point
{
public:
    T x, y;
    Point()
    {
        x = T(0);
        y = T(0);
    }

    Point(T coordX, T coordY)
    {
        x = coordX;
        y = coordY;
    }

    void setCoords(T coordX, T coordY)
    {
        x = coordX;
        y = coordY;
    }
};

class RefCount
{
private:
    int refCount;
public:
    RefCount() : refCount(1) {}
    virtual ~RefCount() {}

    void reference()
    {
        ++refCount;
    }
    void unreference()
    {
        --refCount;

        if (!refCount) {
            delete this;
        }
    }
};

/**
 * @brief Handle back buffers as automatically as possible, and suitable to be used with Glib::RefPtr
 */
class BackBuffer : public RefCount
{

protected:
    int x, y, w, h;  // Rectangle where the colored bar has to be drawn
    Point<int> offset;  // Offset of the source region to draw, relative to the top left corner
    Cairo::RefPtr<Cairo::ImageSurface> surface;
    bool dirty;  // mean that the Surface has to be (re)allocated

public:
    BackBuffer();
    BackBuffer(int w, int h, Cairo::Format format = Cairo::FORMAT_RGB24);
    BackBuffer(int w, int h, Glib::RefPtr<Gdk::Window> referenceWindow);

    // set the destination drawing rectangle; return true if the dimensions are different
    // Note: newW & newH must be > 0
    bool setDrawRectangle(Glib::RefPtr<Gdk::Window> window, Gdk::Rectangle &rectangle, bool updateBackBufferSize = true);
    bool setDrawRectangle(Glib::RefPtr<Gdk::Window> window, int newX, int newY, int newW, int newH, bool updateBackBufferSize = true);
    bool setDrawRectangle(Cairo::Format format, Gdk::Rectangle &rectangle, bool updateBackBufferSize = true);
    bool setDrawRectangle(Cairo::Format format, int newX, int newY, int newW, int newH, bool updateBackBufferSize = true);
    void setSrcOffset(int x, int y);

    void copyRGBCharData(const unsigned char *srcData, int srcX, int srcY, int srcW, int srcH, int srcRowStride, int dstX, int dstY);
    void copySurface(Glib::RefPtr<Gdk::Window> window, Gdk::Rectangle *rectangle = NULL);
    void copySurface(BackBuffer *destBackBuffer, Gdk::Rectangle *rectangle = NULL);
    void copySurface(Cairo::RefPtr<Cairo::ImageSurface> destSurface, Gdk::Rectangle *rectangle = NULL);
    void copySurface(Cairo::RefPtr<Cairo::Context> crDest, Gdk::Rectangle *destRectangle = NULL);

    void setDirty(bool isDirty)
    {
        dirty = isDirty;

        if (!dirty && !surface) {
            dirty = true;
        }
    }
    bool isDirty()
    {
        return dirty;
    }
    // you have to check if the surface is created thanks to surfaceCreated before starting to draw on it
    bool surfaceCreated()
    {
        return surface;
    }
    Cairo::RefPtr<Cairo::ImageSurface> getSurface()
    {
        return surface;
    }
    void setSurface(Cairo::RefPtr<Cairo::ImageSurface> surf)
    {
        surface = surf;
    }
    void deleteSurface()
    {
        if (surface) {
            surface.clear();
        }

        dirty = true;
    }
    // will let you get a Cairo::Context for Cairo drawing operations
    Cairo::RefPtr<Cairo::Context> getContext()
    {
        return Cairo::Context::create(surface);
    }
    int getWidth()
    {
        return surface ? surface->get_width() : 0;    // sending back the allocated width
    }
    int getHeight()
    {
        return surface ? surface->get_height() : 0;    // sending back the allocated height
    }
};

inline void setActiveTextOrIndex (Gtk::ComboBoxText& comboBox, const Glib::ustring& text, int index)
{
    comboBox.set_active_text (text);

    if (comboBox.get_active_row_number () < 0)
        comboBox.set_active (index);
}

#endif
