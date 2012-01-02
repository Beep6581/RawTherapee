/*
 *  This file is part of RawTherapee.
 *
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
#ifndef _PREVIEWMODEPANEL_
#define _PREVIEWMODEPANEL_

#include <gtkmm.h>
#include "adjuster.h"//dev

class ImageArea;
class PreviewModePanel : public Gtk::HBox {

	protected:
        Gtk::ToggleButton* previewR;
        Gtk::ToggleButton* previewG;
        Gtk::ToggleButton* previewB;
        Gtk::ToggleButton* previewL;
        Gtk::ToggleButton* previewFocusMask;
		ImageArea* imageArea;

	public:
		PreviewModePanel (ImageArea* ia);

		void toggleR ();
		void toggleG ();
		void toggleB ();
		void toggleL ();
		void toggleFocusMask ();

		sigc::connection connR, connB, connG, connL, connFocusMask;

		void buttonToggled(Gtk::ToggleButton* tbpreview);

		bool showR         () { return previewR->get_active (); }
		bool showG         () { return previewG->get_active (); }
		bool showB         () { return previewB->get_active (); }
		bool showL         () { return previewL->get_active (); }
		bool showFocusMask () { return previewFocusMask->get_active (); }

};

#endif
