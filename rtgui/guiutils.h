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

bool removeIfThere (Gtk::Container* cont, Gtk::Widget* w, bool increference=true);
void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh);
Glib::ustring removeExtension (const Glib::ustring& filename);
Glib::ustring getExtension (const Glib::ustring& filename);
void drawCrop (Cairo::RefPtr<Cairo::Context> cr, int imx, int imy, int imw, int imh, int startx, int starty, double scale, const rtengine::procparams::CropParams& cparams);

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

enum TOITypes {
	TOI_TEXT,
	TOI_ICON
};

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

#endif
