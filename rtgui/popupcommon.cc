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
#include "../rtengine/safegtk.h"
#include "rtimage.h"

PopUpCommon::PopUpCommon (Gtk::Button* thisButton, const Glib::ustring& label)
    : selected (-1) // -1 means that the button is invalid
    , menu (0)
    , buttonImage (0)
{
    button = thisButton;
    hasMenu = false;
    imageContainer = Gtk::manage( new Gtk::HBox(false, 0));
    button->set_relief (Gtk::RELIEF_NORMAL);
    button->set_border_width (0);
    button->add(*imageContainer);

    if (label.size()) {
        Gtk::Label* buttonLabel = Gtk::manage ( new Gtk::Label(label + " ") );
        imageContainer->pack_start(*buttonLabel, Gtk::PACK_SHRINK, 0);
    }

    // Create the global container and put the button in it
    buttonGroup = Gtk::manage( new Gtk::HBox(false, 0));
    buttonGroup->pack_start(*button, Gtk::PACK_EXPAND_WIDGET, 0);
}

PopUpCommon::~PopUpCommon ()
{
    for (std::vector<RTImage*>::iterator i = images.begin(); i != images.end(); ++i) {
        delete *i;
    }

    delete menu;
    delete buttonImage;
    delete buttonGroup;
}

bool PopUpCommon::addEntry (const Glib::ustring& fileName, const Glib::ustring& label)
{
    if (label.empty ())
         return false;

    // Create the image
    RTImage* newImage = new RTImage(fileName);
    images.push_back(newImage);
    imageFilenames.push_back(fileName);
    int currPos = (int)images.size();
    // Create the menu item
    Gtk::ImageMenuItem* newItem = Gtk::manage(new Gtk::ImageMenuItem (*newImage, label));

    if (selected == -1) {
        // Create the menu on the first item
        menu = new Gtk::Menu ();
        // Create the image for the button
        buttonImage = new RTImage(fileName);
        // Use the first image by default
        imageContainer->pack_start(*buttonImage, Gtk::PACK_EXPAND_WIDGET);
        selected = 0;
    }

    // When there is at least 1 choice, we add the arrow button
    if (images.size() == 1) {
        Gtk::Button* arrowButton = Gtk::manage( new Gtk::Button() );
        RTImage* arrowImage = Gtk::manage( new RTImage("popuparrow.png") );
        arrowButton->add(*arrowImage); //menuSymbol);
        arrowButton->set_relief (Gtk::RELIEF_NONE);
        arrowButton->set_border_width (0);
        buttonGroup->pack_start(*arrowButton, Gtk::PACK_SHRINK, 0);
        arrowButton->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &PopUpCommon::showMenu) );
        hasMenu = true;
    }

    newItem->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &PopUpCommon::entrySelected), currPos - 1));
    menu->attach (*newItem, 0, 1, currPos - 1, currPos);

    return true;
}

// TODO: 'PopUpCommon::removeEntry' method to be created...

void PopUpCommon::entrySelected (int i)
{
    // Emit a a signal if the selected item has changed
    if (setSelected (i))
        message (selected);
}

void PopUpCommon::setItemSensitivity (int i, bool isSensitive) {
    Gtk::Menu_Helpers::MenuList items = menu->items();
    if (i < items.size()) {
        items[i].set_sensitive(isSensitive);
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
        // HACK: Gtk::MenuItem::get_label does not seem to work reliably.
        Gtk::MenuItem& item = menu->items ()[selected];
        Gtk::Label* label = dynamic_cast<Gtk::Label*>(item.get_child ());

        if (label)
            hint += label->get_text ();
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
