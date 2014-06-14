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
#include "../rtengine/rtengine.h"
#include <sstream>
#include <iostream>

Glib::ustring escapeHtmlChars(const Glib::ustring &src);
bool removeIfThere (Gtk::Container* cont, Gtk::Widget* w, bool increference=true);
void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh);
Glib::ustring removeExtension (const Glib::ustring& filename);
Glib::ustring getExtension (const Glib::ustring& filename);
bool confirmOverwrite (Gtk::Window& parent, const std::string& filename);
void writeFailed (Gtk::Window& parent, const std::string& filename);
void drawCrop (Cairo::RefPtr<Cairo::Context> cr, int imx, int imy, int imw, int imh, int startx, int starty, double scale, const rtengine::procparams::CropParams& cparams, bool drawGuide = true);

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

/**
 * @brief subclass of Gtk::ScrolledWindow in order to handle the scrollwheel
 */
class MyScrolledWindow : public Gtk::ScrolledWindow {

	bool on_scroll_event (GdkEventScroll* event);

public:
	MyScrolledWindow();
};

/**
 * @brief subclass of Gtk::ComboBox in order to handle the scrollwheel
 */
class MyComboBox : public Gtk::ComboBox {

	bool on_scroll_event (GdkEventScroll* event);

public:
	MyComboBox ();
};

/**
 * @brief subclass of Gtk::ComboBoxText in order to handle the scrollwheel
 */
class MyComboBoxText : public Gtk::ComboBoxText {

	bool on_scroll_event (GdkEventScroll* event);

public:
	MyComboBoxText ();
};

/**
 * @brief subclass of Gtk::SpinButton in order to handle the scrollwheel
 */
class MySpinButton : public Gtk::SpinButton {

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
class MyHScale : public Gtk::HScale {

	bool on_scroll_event (GdkEventScroll* event);
};

/**
 * @brief subclass of Gtk::FileChooserButton in order to handle the scrollwheel
 */
class MyFileChooserButton : public Gtk::FileChooserButton {

protected:
	bool on_scroll_event (GdkEventScroll* event);

public:
	MyFileChooserButton (const Glib::ustring& title, Gtk::FileChooserAction action=Gtk::FILE_CHOOSER_ACTION_OPEN);
};

/**
 * A class which maintains the last folder for a FileChooserDialog or Button by
 * caching it in a a variable (which can be persisted externally).
 * Each time the user selects a file or folder, the provided variable is updated
 * with the associated folder. The path found in the variable is set in the
 * dialog instance at constructions time of this object.
 */
class FileChooserLastFolderPersister: public Glib::Object {
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
class TextOrIcon : public Gtk::HBox {

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

/**
 * @brief Define a gradient milestone
 */
class GradientMilestone {
public:
	double position;
	double r;
	double g;
	double b;
	double a;

	GradientMilestone(double _p=0., double _r=0., double _g=0., double _b=0., double _a=0.) {
		position = _p; r = _r; g = _g; b = _b; a = _a;
	}
};

/**
 * @brief Handle point coordinates
 */
template <class T>
class Point {
public:
	T x, y;
	Point() {
		x = T(0);
		y = T(0);
	}

	Point(T coordX, T coordY) {
		x = coordX;
		y = coordY;
	}

	void setCoords(T coordX, T coordY) {
		x = coordX;
		y = coordY;
	}
};

/**
 * @brief Handle backbuffers as automatically as possible
 */
class BackBuffer {

protected:
	int x, y, w, h;  // Rectangle where the colored bar has to be drawn
	Cairo::RefPtr<Cairo::Surface> surface;
	bool dirty;  // mean that the Surface has to be (re)allocated

public:
	BackBuffer();

	// set the destination drawing rectangle; return true if the dimensions are different
	bool setDrawRectangle(Glib::RefPtr<Gdk::Window> window, int newX, int newY, int newW, int newH);

	void copySurface(Glib::RefPtr<Gdk::Window> window, GdkRectangle *rectangle=NULL);
	void copySurface(BackBuffer *destBackBuffer, GdkRectangle *rectangle=NULL);
	void copySurface(Cairo::RefPtr<Cairo::Surface> destSurface, GdkRectangle *rectangle=NULL);

	void setDirty(bool isDirty) { dirty = isDirty; if (!dirty && !surface) dirty = true; }
	bool isDirty() { return dirty; }
	// you have to check if the surface is created thanks to surfaceCreated before starting to draw on it
	bool surfaceCreated() { return surface; }
	Cairo::RefPtr<Cairo::Surface> getSurface() { return surface; }
	void deleteSurface() { surface.clear(); dirty=true; }
	// will let you get a Cairo::Context for Cairo drawing operations
	Cairo::RefPtr<Cairo::Context> getContext() { return Cairo::Context::create(surface); }
	int getWidth() { return w; }
	int getHeight() { return h; }
};


#endif
