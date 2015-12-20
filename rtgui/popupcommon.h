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
 *
 *  Class created by Jean-Christophe FRISCH, aka 'Hombre'
 */
#ifndef _POPUPCOMMON_
#define _POPUPCOMMON_

#include <vector>
#include <glibmm/ustring.h>
#include <sigc++/signal.h>

namespace Gtk
{
class HBox;
class Menu;
class Button;
class ImageMenuItem;
}

typedef struct _GdkEventButton GdkEventButton;

class RTImage;

class PopUpCommon
{

public:
    typedef sigc::signal<void, int> type_signal_changed;
    type_signal_changed signal_changed();
    Gtk::HBox* buttonGroup;     // this is the widget to be packed

    PopUpCommon (Gtk::Button* button, const Glib::ustring& label = "");
    virtual ~PopUpCommon ();
    bool addEntry (const Glib::ustring& fileName, const Glib::ustring& label);
    int getEntryCount () const;
    bool setSelected (int entryNum);
    int  getSelected () const;
    void setButtonHint();
    void show ();
    void set_tooltip_text (const Glib::ustring &text);

private:
    type_signal_changed message;

    std::vector<Glib::ustring> imageFilenames;
    std::vector<RTImage*> images;
    Glib::ustring buttonHint;
    RTImage* buttonImage;
    Gtk::HBox* imageContainer;
    Gtk::Menu* menu;
    Gtk::Button* button;
    int selected;
    bool hasMenu;

    void showMenu(GdkEventButton* event);
    void setItemSensitivity (int i, bool isSensitive);

protected:
    void entrySelected (int i);

};

inline PopUpCommon::type_signal_changed PopUpCommon::signal_changed ()
{
    return message;
}

inline int PopUpCommon::getEntryCount () const
{
    return images.size();
}

inline int PopUpCommon::getSelected () const
{
    return selected;
}

#endif
