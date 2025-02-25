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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "indclippedpanel.h"
#include "options.h"
#include "multilangmgr.h"
#include "imagearea.h"
#include "rtimage.h"

IndicateClippedPanel::IndicateClippedPanel (ImageArea* ia) :
    Fon("focusscreen-on"),
    Foff("focusscreen-off"),
    Son("contrastmask-on"),
    Soff("contrastmask-off"),
    iF(Gtk::manage(new RTImage(Foff, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    iS(Gtk::manage(new RTImage(Soff, Gtk::ICON_SIZE_LARGE_TOOLBAR))),
    imageArea(ia)
{
    previewFocusMask = Gtk::manage (new Gtk::ToggleButton ());
    previewFocusMask->set_relief(Gtk::RELIEF_NONE);
    previewFocusMask->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWFOCUSMASK"));
    previewFocusMask->set_image(*iF);

    previewSharpMask = Gtk::manage (new Gtk::ToggleButton ());
    previewSharpMask->set_relief(Gtk::RELIEF_NONE);
    previewSharpMask->set_tooltip_markup (M("MAIN_TOOLTIP_PREVIEWSHARPMASK"));
    previewSharpMask->set_image(*iS);

    Glib::ustring tt;

    indClippedH = Gtk::manage (new Gtk::ToggleButton ());
    indClippedH->set_relief(Gtk::RELIEF_NONE);
    indClippedH->add (*Gtk::manage (new RTImage ("warning-highlights", Gtk::ICON_SIZE_LARGE_TOOLBAR)));
    tt = Glib::ustring::compose("%1\n%2 = %3", M("MAIN_TOOLTIP_INDCLIPPEDH"), M("MAIN_TOOLTIP_THRESHOLD"), options.highlightThreshold);

    if (tt.find("&lt;") == Glib::ustring::npos && tt.find("&gt;") == Glib::ustring::npos) {
        indClippedH->set_tooltip_text (tt);
    } else {
        indClippedH->set_tooltip_markup (tt);
    }

    indClippedS = Gtk::manage (new Gtk::ToggleButton ());
    indClippedS->set_relief(Gtk::RELIEF_NONE);
    indClippedS->add (*Gtk::manage (new RTImage ("warning-shadows", Gtk::ICON_SIZE_LARGE_TOOLBAR)));
    tt = Glib::ustring::compose("%1\n%2 = %3", M("MAIN_TOOLTIP_INDCLIPPEDS"), M("MAIN_TOOLTIP_THRESHOLD"), options.shadowThreshold);

    if (tt.find("&lt;") == Glib::ustring::npos && tt.find("&gt;") == Glib::ustring::npos) {
        indClippedS->set_tooltip_text (tt);
    } else {
        indClippedS->set_tooltip_markup (tt);
    }

    previewFocusMask->set_active (false);
    previewSharpMask->set_active (false);
    indClippedH->set_active (options.showClippedHighlights);
    indClippedS->set_active (options.showClippedShadows);

    pack_start (*previewFocusMask, Gtk::PACK_SHRINK, 0);
    pack_start (*previewSharpMask, Gtk::PACK_SHRINK, 0);
    pack_start (*indClippedS, Gtk::PACK_SHRINK, 0);
    pack_start (*indClippedH, Gtk::PACK_SHRINK, 0);

    connSharpMask = previewSharpMask->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled), previewSharpMask) );
    connFocusMask = previewFocusMask->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled), previewFocusMask) );
    connClippedS = indClippedS->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled), indClippedS) );
    connClippedH = indClippedH->signal_toggled().connect( sigc::bind(sigc::mem_fun(*this, &IndicateClippedPanel::buttonToggled), indClippedH) );

    show_all ();
}

// inverts a toggle programmatically
void IndicateClippedPanel::toggleClipped (bool highlights)
{
    if (highlights) {
        indClippedH->set_active(!indClippedH->get_active());
    } else {
        indClippedS->set_active(!indClippedS->get_active());
    }
}

void IndicateClippedPanel::toggleFocusMask ()
{
    previewFocusMask->set_active(!previewFocusMask->get_active());
}

void IndicateClippedPanel::silentlyDisableSharpMask ()
{
    ConnectionBlocker conBlocker(connSharpMask);
    previewSharpMask->set_active(false);
    iS->set_from_icon_name(Soff);

}

void IndicateClippedPanel::toggleSharpMask ()
{
    previewSharpMask->set_active(!previewSharpMask->get_active());
}

void IndicateClippedPanel::buttonToggled (Gtk::ToggleButton* tb)
{

    connFocusMask.block(true);
    connSharpMask.block(true);
    connClippedS.block(true);
    connClippedH.block(true);

    if (tb == previewFocusMask) {
        if (indClippedS->get_active()) {
            indClippedS->set_active(false);
        }
        if (indClippedH->get_active()) {
            indClippedH->set_active(false);
        }
        previewSharpMask->set_active(false);
    } else if (tb == previewSharpMask) {
        if (indClippedS->get_active()) {
            indClippedS->set_active(false);
        }
        if (indClippedH->get_active()) {
            indClippedH->set_active(false);
        }
        previewFocusMask->set_active(false);
    } else {
        previewFocusMask->set_active(false);
        previewSharpMask->set_active(false);
    }

    imageArea->sharpMaskSelected(previewSharpMask->get_active());
    iF->set_from_icon_name(previewFocusMask->get_active() ? Fon : Foff);
    iS->set_from_icon_name(previewSharpMask->get_active() ? Son : Soff);

    connFocusMask.block(false);
    connSharpMask.block(false);
    connClippedS.block(false);
    connClippedH.block(false);

    imageArea->queue_draw ();

    // this will redraw the linked Before image area
    // which is set when before/after view is enabled
    if (imageArea->iLinkedImageArea != nullptr) {
        imageArea->iLinkedImageArea->queue_draw ();
    }
}

IndicateClippedPanel::~IndicateClippedPanel () {}
