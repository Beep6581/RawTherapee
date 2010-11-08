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

#include <popupbutton.h>

/*
 * PopUpButton::PopUpButton (const Glib::ustring& label, bool imgRight)
 *
 * Creates a button with a contextual menu where you can select an item that the button content will reflect
 *
 * Parameters:
 * 		label = label displayed in the button
 * 		imRight = 0: the image is displayed at the left of the label (default)
 * 				  1: the image is displayed at the right of the label
 */
PopUpButton::PopUpButton (const Glib::ustring& label, bool imgRight) : Gtk::Button(), PopUpCommon(this, label, imgRight) { }

void PopUpButton::show() {
	PopUpCommon::show();
}
void PopUpButton::set_tooltip_text (const Glib::ustring &text) {
	PopUpCommon::set_tooltip_text (text);
}
