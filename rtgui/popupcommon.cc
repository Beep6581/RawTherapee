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

#include <popupcommon.h>
#include <safegtk.h>

PopUpCommon::PopUpCommon (Gtk::Button* thisButton, const Glib::ustring& label, bool imgRight) {
	button = thisButton;
	if (label.size()) {
		hasText = true;
		button->set_label(label + " ");
	}
	else
		hasText = false;
	// Create the list entry
	imagePaths.clear();
	images.clear();
	sItems.clear();
	items.clear();
	selected = -1;		// -1 : means that the button is invalid
	menu = 0;
	buttonImage = 0;
	buttonHint = "";
	imageRight = imgRight; // By default, image is on the left in the menu
}

PopUpCommon::~PopUpCommon () {
    for (std::vector<Gtk::Image*>::iterator i = images.begin(); i != images.end(); ++i)
    {
        delete *i;
    }
    for (std::vector<Gtk::ImageMenuItem*>::iterator i = items.begin(); i != items.end(); ++i)
    {
        delete *i;
    }
    if (menu) delete menu;
    if (buttonImage) delete buttonImage;
}

PopUpCommon::type_signal_changed PopUpCommon::signal_changed() {
	return message;
}

bool PopUpCommon::addEntry (Glib::ustring imagePath, Glib::ustring label) {
	bool added = false;
	if ( Glib::file_test(safe_locale_from_utf8(imagePath), Glib::FILE_TEST_EXISTS) && label.size() ) {
		imagePaths.push_back(imagePath);
		sItems.push_back(label);
		// Create the image
		Gtk::Image* newImage = new Gtk::Image(imagePath);
		images.push_back(newImage);
		int currPos = (int)images.size();
		// Create the menu item
		Gtk::ImageMenuItem* newItem = new Gtk::ImageMenuItem (*newImage, label);
		items.push_back(newItem);
		if (selected == -1) {
			// Create the menu on the first item
			menu = new Gtk::Menu ();
			// Create the image for the button
			buttonImage = new Gtk::Image(imagePath);
			// Use the first image by default
			if (hasText) {
				button->set_image_position(imageRight ? Gtk::POS_RIGHT : Gtk::POS_LEFT);
				button->set_image(*buttonImage);
			}
			else {
				button->add(*buttonImage);
			}
			selected = 0;
		}
		newItem->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &PopUpCommon::entrySelected), currPos-1));
		menu->attach (*newItem, 0, 1, currPos-1, currPos);
		// When there is at least 2 choice, we add the RMB connector
		if (images.size() == 2) {
			button->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &PopUpCommon::showMenu) );
		}
		// The item has been created
		added = true;
	}
	return added;
}

// TODO: 'PopUpCommon::removeEntry' method to be created...

void PopUpCommon::entrySelected (int i) {
	if (setSelected((unsigned int)i))
		// Emit a a signal if the selected item has changed
		message.emit(selected);
}

/*
 * Set the button image with the selected item
 */
bool PopUpCommon::setSelected (int entryNum) {
	if (entryNum < 0 || entryNum > (int)images.size()-1 || (int)entryNum == selected)
		return false;
	else {
		// Maybe we could do something better than loading the image file each time the selection is changed !?
		buttonImage->set(imagePaths.at(entryNum));
		selected = entryNum;
		setButtonHint();
		return true;
	}
}

void PopUpCommon::show() {
	menu->reposition();
	setButtonHint();
	menu->show_all();
	button->show();
}

void PopUpCommon::setButtonHint() {
	button->set_tooltip_text(buttonHint.size() ? buttonHint + " : " + sItems.at(selected) : sItems.at(selected));
}

void PopUpCommon::showMenu(GdkEventButton* event) {
	if (event->button == 3)	menu->popup(event->button, event->time);
}

void PopUpCommon::set_tooltip_text (const Glib::ustring &text) {
	buttonHint = text;
}
