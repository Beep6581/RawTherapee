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
#ifndef _INDCLIPPEDPANEL_
#define _INDCLIPPEDPANEL_

#include <gtkmm.h>

class ImageArea;
class IndicateClippedPanel : public Gtk::HBox
{

protected:
    Gtk::ToggleButton* indclippedh;
    Gtk::ToggleButton* indclippeds;
    ImageArea* imageArea;

public:
    IndicateClippedPanel (ImageArea* ia);

    void buttonToggled   ();

    void toggleClipped (bool highlights);  // inverts a toggle programmatically

    bool showClippedShadows    ()
    {
        return indclippeds->get_active ();
    }
    bool showClippedHighlights ()
    {
        return indclippedh->get_active ();
    }
};

#endif
