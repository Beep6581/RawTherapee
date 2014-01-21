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

using namespace rtengine::procparams;


class Frame2: public Gtk::Frame
{
	Gtk::Container *pC;
public:
	Frame2( Gtk::Container *p):pC(p){}
	~Frame2( ){ delete pC;}
};

FoldableToolPanel::FoldableToolPanel(Gtk::Box* content) : ToolPanel(), parentContainer(NULL), exp(NULL) {

	exp = Gtk::manage (new Gtk::Expander ());
	exp->set_border_width (4);
	exp->set_use_markup (true);
	exp->signal_button_release_event().connect_notify( sigc::mem_fun(this, &FoldableToolPanel::foldThemAll) );

	Frame2* pframe = Gtk::manage (new Frame2 (content));

	pframe->set_name ("ToolPanel");

	pframe->add (*content);

	exp->add (*pframe);
	pframe->set_shadow_type (Gtk::SHADOW_ETCHED_IN);
	pframe->show ();
	exp->show ();
}

void FoldableToolPanel::foldThemAll (GdkEventButton* event) {
	if (event->button == 3)	{
		if (listener)
		  (static_cast<ToolPanelCoordinator*>(listener))->foldAllButOne( parentContainer, this);
		else
		  (static_cast<ToolPanelCoordinator*>(tmp))->foldAllButOne( parentContainer, this);
	}
}

void FoldableToolPanel::setLabel (Glib::ustring label, bool need100Percent) {
	if (!need100Percent)
		exp->set_label(Glib::ustring("<b>") + label + Glib::ustring("</b>"));
	else {
		Gtk::Label *labelWidget = Gtk::manage (new Gtk::Label(Glib::ustring("<b>") + label + Glib::ustring("</b>")));
		labelWidget->set_use_markup();
		RTImage *image = Gtk::manage (new RTImage("zoom-100-identifier.png"));
		image->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
		Gtk::HBox *hbox = Gtk::manage (new Gtk::HBox());

		hbox->set_spacing(4);
		hbox->pack_start(*labelWidget, false, false, 0);
		hbox->pack_end(*image, false, false, 0);
		exp->set_label_widget(*hbox);
		exp->set_label_fill();
	}
}
