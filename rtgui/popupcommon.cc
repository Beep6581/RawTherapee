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

#include "multilangmgr.h"
#include "popupcommon.h"
#include "../rtengine/safegtk.h"
#include "rtimage.h"

PopUpCommon::PopUpCommon (Gtk::Button* thisButton, const Glib::ustring& label) {
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
	// Create the list entry
	imageFilenames.clear();
	images.clear();
	sItems.clear();
	items.clear();
	selected = -1;		// -1 : means that the button is invalid
	menu = 0;
	buttonImage = 0;
	buttonHint = "";
}

PopUpCommon::~PopUpCommon () {
    for (std::vector<RTImage*>::iterator i = images.begin(); i != images.end(); ++i)
    {
        delete *i;
    }
    for (std::vector<Gtk::ImageMenuItem*>::iterator i = items.begin(); i != items.end(); ++i)
    {
        delete *i;
    }
    if (menu) delete menu;
    if (buttonImage) delete buttonImage;
    delete buttonGroup;
}

PopUpCommon::type_signal_changed PopUpCommon::signal_changed() {
	return message;
}

bool PopUpCommon::addEntry (Glib::ustring fileName, Glib::ustring label) {
	bool added = false;
	if ( label.size() ) {
		imageFilenames.push_back(fileName);
		sItems.push_back(label);
		// Create the image
		RTImage* newImage = new RTImage(fileName);
		images.push_back(newImage);
		int currPos = (int)images.size();
		// Create the menu item
		Gtk::ImageMenuItem* newItem = new Gtk::ImageMenuItem (*newImage, label);
		items.push_back(newItem);
		if (selected == -1) {
			// Create the menu on the first item
			menu = new Gtk::Menu ();
			// Create the image for the button
			buttonImage = new RTImage(fileName);
			// Use the first image by default
			imageContainer->pack_start(*buttonImage,Gtk::PACK_EXPAND_WIDGET);
			selected = 0;
		}
		// When there is at least 1 choice, we add the arrow button
		if (images.size() == 1) {
			Gtk::Button* arrowButton = Gtk::manage( new Gtk::Button() );
			RTImage* arrowImage = Gtk::manage( new RTImage("popuparrow.png") );
			arrowButton->add(*arrowImage); //menuSymbol);
			arrowButton->set_relief (Gtk::RELIEF_NONE);
			arrowButton->set_border_width (0);
			buttonGroup->pack_start(*arrowButton,Gtk::PACK_SHRINK, 0);
			arrowButton->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &PopUpCommon::showMenu) );
			hasMenu = true;
		}
		newItem->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &PopUpCommon::entrySelected), currPos-1));
		menu->attach (*newItem, 0, 1, currPos-1, currPos);
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
	if (entryNum < 0 || entryNum > ((int)images.size()-1) || (int)entryNum == selected)
		return false;
	else {
		// Maybe we could do something better than loading the image file each time the selection is changed !?
		buttonImage->changeImage(imageFilenames.at(entryNum));
		selected = entryNum;
		setButtonHint();
		return true;
	}
}

void PopUpCommon::show() {
	menu->reposition();
	setButtonHint();
	menu->show_all();
	buttonGroup->show_all();
}

void PopUpCommon::setButtonHint() {
	Glib::ustring hint;
	if (!buttonHint.empty()) {
		hint = buttonHint;
		if (selected > -1)
			hint += " ";
	}
	if (selected > -1)
		hint += sItems.at(selected);
	button->set_tooltip_markup(hint);
}

void PopUpCommon::showMenu(GdkEventButton* event) {
	if (event->button == 1)	menu->popup(event->button, event->time);
}

void PopUpCommon::set_tooltip_text (const Glib::ustring &text) {
	buttonHint = text;
	setButtonHint();
}
