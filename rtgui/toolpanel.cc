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


class Frame2: public Gtk::EventBox
{
	Gtk::Container *pC;

public:
	Frame2( Gtk::Container *p):pC(p){ updateStyle(); }
	~Frame2( ){ delete pC;}

	void updateStyle() {
		set_border_width(options.slimUI ? 2 : 8);  // Outer space around the tool's frame 2:7
	}

	void on_style_changed (const Glib::RefPtr<Gtk::Style>& style) { updateStyle(); }
	bool on_expose_event(GdkEventExpose* event);
};

bool Frame2::on_expose_event(GdkEventExpose* event) {
	bool retVal = Gtk::EventBox::on_expose_event(event);

	if (!options.useSystemTheme) {
		Glib::RefPtr<Gdk::Window> window = get_window();
		Glib::RefPtr<Gtk::Style> style = get_style ();
		Cairo::RefPtr<Cairo::Context> cr = window->create_cairo_context();

		int x_, y_, w_, h_, foo;
		window->get_geometry(x_, y_, w_, h_, foo);
		double x = 0.;
		double y = 0.;
		double w = double(w_);
		double h = double(h_);

		cr->set_antialias (Cairo::ANTIALIAS_NONE);

		// draw a frame
		cr->set_line_width (1.0);
		Gdk::Color c = style->get_fg (Gtk::STATE_NORMAL);
		cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
		cr->move_to(x+0.5, y+0.5);
		cr->line_to(x+w, y+0.5);
		cr->line_to(x+w, y+h);
		cr->line_to(x+0.5, y+h);
		cr->line_to(x+0.5, y+0.5);
		cr->stroke ();
	}
	return retVal;
}

ToolVBox::ToolVBox() {
	updateStyle();
}

void ToolVBox::updateStyle() {
	if (options.slimUI) {
		set_spacing(1);       // Vertical space between tools
		set_border_width(1);  // Space separating the tab's frame and the tools
	}
	else {
		set_spacing(2);       // Vertical space between tools
		set_border_width(1);  // Space separating the tab's frame and the tools  3
	}
}

void ToolVBox::on_style_changed (const Glib::RefPtr<Gtk::Style>& style) {
	updateStyle();
}

ToolParamBlock::ToolParamBlock() {
	updateStyle();
}

void ToolParamBlock::updateStyle() {
	if (options.slimUI) {
		set_spacing(2);       // Vertical space between parameters in a single tool
		set_border_width(6);  // Space separating the parameters of a tool and its surrounding frame  6
	}
	else {
		set_spacing(4);       // Vertical space between parameters in a single tool
		set_border_width(8);  // Space separating the parameters of a tool and its surrounding frame  8
	}
}

void ToolParamBlock::on_style_changed (const Glib::RefPtr<Gtk::Style>& style) {
	updateStyle();
}

FoldableToolPanel::FoldableToolPanel(Gtk::Box* content, Glib::ustring toolName, Glib::ustring UILabel, bool need11, bool useEnabled) : ToolPanel(toolName, need11), parentContainer(NULL), exp(NULL), lastEnabled(true)
{
	if (!content)
		return;

//	exp->set_border_width (5);
//	exp->set_use_markup (true);
	if (need11) {
		Gtk::HBox *titleHBox = Gtk::manage(new Gtk::HBox());

		Gtk::Label *label = Gtk::manage(new Gtk::Label());
		label->set_markup(Glib::ustring("<b>") + escapeHtmlChars(UILabel) + Glib::ustring("</b>"));
		label->set_alignment(Gtk::ALIGN_LEFT, Gtk::ALIGN_CENTER);
		titleHBox->pack_start(*label, Gtk::PACK_EXPAND_WIDGET, 0);

		RTImage *image = Gtk::manage (new RTImage("zoom-100-identifier.png"));
		image->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
		titleHBox->pack_end(*image, Gtk::PACK_SHRINK, 0);

		exp = Gtk::manage (new MyExpander (useEnabled, titleHBox));
	}
	else {
		exp = Gtk::manage (new MyExpander (useEnabled, UILabel));
	}
	exp->signal_button_release_event().connect_notify( sigc::mem_fun(this, &FoldableToolPanel::foldThemAll) );
	enaConn = signal_enabled_toggled().connect( sigc::mem_fun(*this, &FoldableToolPanel::enabled_toggled) );

	Frame2* pframe = Gtk::manage (new Frame2 (content));

	pframe->set_name ("ToolPanel");

	pframe->add (*content);

	exp->add (*pframe);
	pframe->show ();
	exp->show ();
}

void FoldableToolPanel::foldThemAll (GdkEventButton* event) {
	if (event->button == 3) {
		if (listener)
		  (static_cast<ToolPanelCoordinator*>(listener))->foldAllButOne( parentContainer, this);
		else
		  (static_cast<ToolPanelCoordinator*>(tmp))->foldAllButOne( parentContainer, this);
	}
}

void FoldableToolPanel::enabled_toggled() {
	if (multiImage) {
		if (exp->get_inconsistent()) {
			exp->set_inconsistent (false);
			enaConn.block (true);
			exp->setEnabled (false);
			enaConn.block (false);
		}
		else if (lastEnabled)
			exp->set_inconsistent (true);

		lastEnabled = exp->getEnabled();
	}
	enabledChanged();
}

bool FoldableToolPanel::get_inconsistent() {
	return exp->get_inconsistent();
}

void FoldableToolPanel::set_inconsistent(bool isInconsistent) {
	exp->set_inconsistent(isInconsistent);
}

bool FoldableToolPanel::getEnabled() {
	return exp->getEnabled();
}

// do not emit the enabled_toggled event
void FoldableToolPanel::setEnabled(bool isEnabled) {
	enaConn.block (true);
	exp->setEnabled(isEnabled);
	lastEnabled = isEnabled;
	enaConn.block (false);
}

void FoldableToolPanel::setEnabledTooltipMarkup(Glib::ustring tooltipMarkup) {
	if (exp)
		exp->set_tooltip_markup(tooltipMarkup);
}

void FoldableToolPanel::setEnabledTooltipText(Glib::ustring tooltipText) {
	if (exp)
		exp->set_tooltip_text(tooltipText);
}



