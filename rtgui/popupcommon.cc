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

#include <gtkmm.h>
#include "multilangmgr.h"
#include "popupcommon.h"
#include "rtimage.h"
#include "guiutils.h"

PopUpCommon::PopUpCommon (Gtk::Button* thisButton, const Glib::ustring& label)
    : selected (-1) // -1 means that the button is invalid
    , menu (nullptr)
    , buttonImage (nullptr)
{
    button = thisButton;
    hasMenu = false;
    imageContainer = Gtk::manage( new Gtk::Grid());
    setExpandAlignProperties(imageContainer, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    button->set_relief (Gtk::RELIEF_NORMAL);
    button->add(*imageContainer);

    if (!label.empty()) {
        Gtk::Label* buttonLabel = Gtk::manage ( new Gtk::Label(label + " ") );
        setExpandAlignProperties(buttonLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
        imageContainer->attach(*buttonLabel, 0, 0, 1, 1);
    }

    // Create the global container and put the button in it
    buttonGroup = Gtk::manage( new Gtk::Grid());
    setExpandAlignProperties(buttonGroup, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    buttonGroup->attach(*button, 0, 0, 1, 1);
}

PopUpCommon::~PopUpCommon ()
{
    delete menu;
    delete buttonImage;
}

bool PopUpCommon::addEntry (const Glib::ustring& fileName, const Glib::ustring& label)
{
    if (label.empty ())
         return false;

    // Create the menu item and image
    MyImageMenuItem* newItem = Gtk::manage (new MyImageMenuItem (label, fileName));
    imageFilenames.push_back (fileName);
    images.push_back (newItem->getImage ());

    if (selected == -1) {
        // Create the menu on the first item
        menu = new Gtk::Menu ();
        // Create the image for the button
        buttonImage = new RTImage(fileName);
        setExpandAlignProperties(buttonImage, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
        // Use the first image by default
        imageContainer->attach_next_to(*buttonImage, Gtk::POS_RIGHT, 1, 1);
        selected = 0;
    }

    // When there is at least 1 choice, we add the arrow button
    if (images.size() == 1) {
        Gtk::Button* arrowButton = Gtk::manage( new Gtk::Button() );
        RTImage* arrowImage = Gtk::manage( new RTImage("popuparrow.png") );
        setExpandAlignProperties(arrowButton, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
        arrowButton->add(*arrowImage); //menuSymbol);
        buttonGroup->attach_next_to(*arrowButton, *button, Gtk::POS_RIGHT, 1, 1);
        arrowButton->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &PopUpCommon::showMenu) );
        button->get_style_context()->add_class("Left");
        arrowButton->get_style_context()->add_class("Right");
        hasMenu = true;
    }

    newItem->signal_activate ().connect (sigc::bind (sigc::mem_fun (*this, &PopUpCommon::entrySelected), images.size () - 1));
    menu->append (*newItem);

    return true;
}

// TODO: 'PopUpCommon::removeEntry' method to be created...

void PopUpCommon::entrySelected (int i)
{
    // Emit a a signal if the selected item has changed
    if (setSelected (i))
        message (selected);
}

void PopUpCommon::setItemSensitivity (int index, bool isSensitive) {
    const auto items = menu->get_children ();
    if (index < items.size ()) {
        items[index]->set_sensitive (isSensitive);
    }
}


/*
 * Set the button image with the selected item
 */
bool PopUpCommon::setSelected (int entryNum)
{
    if (entryNum < 0 || entryNum > ((int)images.size() - 1) || (int)entryNum == selected) {
        return false;
    } else {
        // Maybe we could do something better than loading the image file each time the selection is changed !?
        buttonImage->changeImage(imageFilenames.at(entryNum));
        selected = entryNum;
        setButtonHint();
        return true;
    }
}

void PopUpCommon::show()
{
    menu->reposition();
    setButtonHint();
    menu->show_all();
    buttonGroup->show_all();
}

void PopUpCommon::setButtonHint()
{
    Glib::ustring hint;

    if (!buttonHint.empty()) {
        hint = buttonHint;

        if (selected > -1) {
            hint += " ";
        }
    }

    if (selected > -1) {
        auto widget = menu->get_children ()[selected];
        auto item = dynamic_cast<MyImageMenuItem*>(widget);

        if (item) {
            hint += item->getLabel ()->get_text ();
        }
    }

    button->set_tooltip_markup(hint);
}

void PopUpCommon::showMenu(GdkEventButton* event)
{
    if (event->button == 1) {
        menu->popup(event->button, event->time);
    }
}

void PopUpCommon::set_tooltip_text (const Glib::ustring &text)
{
    buttonHint = text;
    setButtonHint();
}
