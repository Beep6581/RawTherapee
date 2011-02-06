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
#include <toolpanel.h>
#include <toolpanelcoord.h>

using namespace rtengine::procparams;

FoldableToolPanel::FoldableToolPanel(Gtk::Box* content) : ToolPanel(), parentContainer(NULL), exp(NULL) {

	// Gtk::Expander* exp = new Gtk::Expander ();
	// exp->set_label_widget (*(new ILabel (Glib::ustring("<b>") + label + "</b>")));
	exp = Gtk::manage (new Gtk::Expander ());
	exp->set_border_width (4);
	//Glib::RefPtr<Gtk::Style> *style = new Gtk::Style();
	//exp->set_style()
	exp->set_use_markup (true);
	exp->signal_button_release_event().connect_notify( sigc::mem_fun(this, &FoldableToolPanel::foldThemAll) );

	Gtk::Frame* pframe = Gtk::manage (new Gtk::Frame ());

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
			((ToolPanelCoordinator*)listener)->foldAllButOne( parentContainer, this);
		else
			((ToolPanelCoordinator*)tmp)->foldAllButOne( parentContainer, this);
	}
}
