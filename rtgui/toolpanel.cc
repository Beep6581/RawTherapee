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
 */
#include "toolpanel.h"
#include "toolpanelcoord.h"
#include "guiutils.h"

using namespace rtengine::procparams;


ToolVBox::ToolVBox()
{
    updateStyle();
}

void ToolVBox::updateStyle()
{
    set_spacing(1);       // Vertical space between tools
    set_border_width(3);  // Space separating the tab's frame and the tools
}

void ToolVBox::on_style_updated ()
{
    updateStyle();
}

ToolParamBlock::ToolParamBlock()
{
    updateStyle();
}

void ToolParamBlock::updateStyle()
{
    set_spacing(2);       // Vertical space between parameters in a single tool
    set_border_width(5);  // Space separating the parameters of a tool and its surrounding frame
}

void ToolParamBlock::on_style_updated ()
{
    updateStyle();
}

FoldableToolPanel::FoldableToolPanel(Gtk::Box* content, Glib::ustring toolName, Glib::ustring UILabel, bool need11, bool useEnabled) : ToolPanel(toolName, need11), parentContainer(nullptr), exp(nullptr), lastEnabled(true)
{
    if (!content) {
        return;
    }

//  exp->set_use_markup (true);
    if (need11) {
        Gtk::HBox *titleHBox = Gtk::manage(new Gtk::HBox());

        Gtk::Label *label = Gtk::manage(new Gtk::Label());
        label->set_markup(Glib::ustring("<b>") + escapeHtmlChars(UILabel) + Glib::ustring("</b>"));
        label->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
        titleHBox->pack_start(*label, Gtk::PACK_EXPAND_WIDGET, 0);

        RTImage *image = Gtk::manage (new RTImage("zoom-100-identifier.png"));
        image->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
        titleHBox->pack_end(*image, Gtk::PACK_SHRINK, 0);

        exp = Gtk::manage (new MyExpander (useEnabled, titleHBox));
    } else {
        exp = Gtk::manage (new MyExpander (useEnabled, UILabel));
    }

    exp->signal_button_release_event().connect_notify( sigc::mem_fun(this, &FoldableToolPanel::foldThemAll) );
    enaConn = signal_enabled_toggled().connect( sigc::mem_fun(*this, &FoldableToolPanel::enabled_toggled) );

    exp->add (*content);
    exp->show ();
}

void FoldableToolPanel::foldThemAll (GdkEventButton* event)
{
    if (event->button == 3) {
        if (listener) {
            (static_cast<ToolPanelCoordinator*>(listener))->foldAllButOne( parentContainer, this);
        } else {
            (static_cast<ToolPanelCoordinator*>(tmp))->foldAllButOne( parentContainer, this);
        }
    }
}

void FoldableToolPanel::enabled_toggled()
{
    if (multiImage) {
        if (exp->get_inconsistent()) {
            exp->set_inconsistent (false);
            enaConn.block (true);
            exp->setEnabled (false);
            enaConn.block (false);
        } else if (lastEnabled) {
            exp->set_inconsistent (true);
        }

        lastEnabled = exp->getEnabled();
    }

    enabledChanged();
}

bool FoldableToolPanel::get_inconsistent()
{
    return exp->get_inconsistent();
}

void FoldableToolPanel::set_inconsistent(bool isInconsistent)
{
    exp->set_inconsistent(isInconsistent);
}

void FoldableToolPanel::setLevel (int level)
{
    if (exp) {
        exp->setLevel(level);
    }
}

bool FoldableToolPanel::getEnabled()
{
    return exp->getEnabled();
}

// do not emit the enabled_toggled event
void FoldableToolPanel::setEnabled(bool isEnabled)
{
    enaConn.block (true);
    exp->setEnabled(isEnabled);
    lastEnabled = isEnabled;
    enaConn.block (false);
}

void FoldableToolPanel::setEnabledTooltipMarkup(Glib::ustring tooltipMarkup)
{
    if (exp) {
        exp->set_tooltip_markup(tooltipMarkup);
    }
}

void FoldableToolPanel::setEnabledTooltipText(Glib::ustring tooltipText)
{
    if (exp) {
        exp->set_tooltip_text(tooltipText);
    }
}



