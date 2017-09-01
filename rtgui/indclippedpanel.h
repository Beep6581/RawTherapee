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
    Gtk::Image* iFon, *iFoff;
    Gtk::ToggleButton* previewFocusMask;
    Gtk::ToggleButton* indClippedH;
    Gtk::ToggleButton* indClippedS;
    ImageArea* imageArea;

public:
    explicit IndicateClippedPanel(ImageArea* ia);
    ~IndicateClippedPanel();

    void buttonToggled(Gtk::ToggleButton* tb);
    void toggleClipped(bool highlights);  // inverts a toggle programmatically
    void toggleFocusMask();

    sigc::connection connFocusMask, connClippedS, connClippedH;


    bool showFocusMask ()
    {
        return previewFocusMask->get_active ();
    }
    bool showClippedShadows()
    {
        return indClippedS->get_active();
    }
    bool showClippedHighlights()
    {
        return indClippedH->get_active();
    }
};

#endif
