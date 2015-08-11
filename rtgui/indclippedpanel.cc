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
#include "indclippedpanel.h"
#include "options.h"
#include "multilangmgr.h"
#include "imagearea.h"
#include "rtimage.h"

IndicateClippedPanel::IndicateClippedPanel (ImageArea* ia) : imageArea(ia)
{

    Glib::ustring tt;

    indclippedh = Gtk::manage (new Gtk::ToggleButton ());
    indclippedh->set_relief(Gtk::RELIEF_NONE);
    indclippedh->add (*Gtk::manage (new RTImage ("warnhl.png")));
    tt = Glib::ustring::compose("%1\n%2 = %3", M("MAIN_TOOLTIP_INDCLIPPEDH"), M("MAIN_TOOLTIP_THRESHOLD"), options.highlightThreshold);

    if (tt.find("&lt;") == Glib::ustring::npos && tt.find("&gt;") == Glib::ustring::npos) {
        indclippedh->set_tooltip_text (tt);
    } else {
        indclippedh->set_tooltip_markup (tt);
    }

    indclippeds = Gtk::manage (new Gtk::ToggleButton ());
    indclippeds->set_relief(Gtk::RELIEF_NONE);
    indclippeds->add (*Gtk::manage (new RTImage ("warnsh.png")));
    tt = Glib::ustring::compose("%1\n%2 = %3", M("MAIN_TOOLTIP_INDCLIPPEDS"), M("MAIN_TOOLTIP_THRESHOLD"), options.shadowThreshold);

    if (tt.find("&lt;") == Glib::ustring::npos && tt.find("&gt;") == Glib::ustring::npos) {
        indclippeds->set_tooltip_text (tt);
    } else {
        indclippeds->set_tooltip_markup (tt);
    }

    indclippedh->set_active (options.showClippedHighlights);
    indclippeds->set_active (options.showClippedShadows);

    pack_start (*indclippedh, Gtk::PACK_SHRINK, 0);
    pack_start (*indclippeds, Gtk::PACK_SHRINK, 0);

    indclippedh->signal_toggled().connect( sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled) );
    indclippeds->signal_toggled().connect( sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled) );

    show_all ();
}

// inverts a toggle programmatically
void IndicateClippedPanel::toggleClipped (bool highlights)
{
    if (highlights) {
        indclippedh->set_active(!indclippedh->get_active());
    } else {
        indclippeds->set_active(!indclippeds->get_active());
    }
}

void IndicateClippedPanel::buttonToggled ()
{
    imageArea->queue_draw ();

    // this will redraw the linked Before image area
    // which is set when before/after view is enabled
    if (imageArea->iLinkedImageArea != NULL) {
        imageArea->iLinkedImageArea->queue_draw ();
    }
}
