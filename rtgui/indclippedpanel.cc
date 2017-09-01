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

    iFon  = new RTImage ("previewmodeF-focusScreen-on.png");
    iFoff = new RTImage ("previewmodeF-focusScreen-off.png");

    previewFocusMask = Gtk::manage (new Gtk::ToggleButton ());
    previewFocusMask->set_relief(Gtk::RELIEF_NONE);
    previewFocusMask->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWFOCUSMASK"));
    previewFocusMask->set_image(*iFoff);

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

    previewFocusMask->set_active (false);
    indclippedh->set_active (options.showClippedHighlights);
    indclippeds->set_active (options.showClippedShadows);

    pack_start (*previewFocusMask, Gtk::PACK_SHRINK, 0);
    pack_start (*indclippeds, Gtk::PACK_SHRINK, 0);
    pack_start (*indclippedh, Gtk::PACK_SHRINK, 0);

    connFocusMask = previewFocusMask->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled), previewFocusMask) );
    connClippedS = indclippeds->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled), indclippeds) );
    connClippedH = indclippedh->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled), indclippedh) );

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

void IndicateClippedPanel::toggleFocusMask ()
{
    previewFocusMask->set_active(!previewFocusMask->get_active());
}

void IndicateClippedPanel::buttonToggled (Gtk::ToggleButton* tb)
{

    connFocusMask.block(true);
    connClippedS.block(true);
    connClippedH.block(true);

    if (tb != previewFocusMask) {
        previewFocusMask->set_active(false);
    } else {
        if (indclippeds->get_active()) {
            indclippeds->set_active(false);
        }
        if (indclippedh->get_active()) {
            indclippedh->set_active(false);
        }
    }

    previewFocusMask->set_image(previewFocusMask->get_active() ? *iFon : *iFoff);

    connFocusMask.block(false);
    connClippedS.block(false);
    connClippedH.block(false);

    imageArea->queue_draw ();

    // this will redraw the linked Before image area
    // which is set when before/after view is enabled
    if (imageArea->iLinkedImageArea != nullptr) {
        imageArea->iLinkedImageArea->queue_draw ();
    }
}

IndicateClippedPanel::~IndicateClippedPanel ()
{
    delete iFon;
    delete iFoff;
}
