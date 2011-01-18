/*
 *  This file is part of RawTherapee.
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
 *
 *  2010 Ilya Popov <ilia_popov@rambler.ru>
 */
 
#ifndef HSVEQUALIZER_H_INCLUDED
#define HSVEQUALIZER_H_INCLUDED

#include <gtkmm.h>
#include <adjuster.h>
#include <toolpanel.h>
#include <guiutils.h>


class HSVEqualizer : public Gtk::VBox, public AdjusterListener, public ToolPanel 
{

protected:

    Gtk::CheckButton * enabled;
	Gtk::ComboBoxText* hsvchannel;
	
	Gtk::VBox* satbox;
    Gtk::VBox* valbox;
    Gtk::VBox* huebox;

    Adjuster* sat[8]; 
    Adjuster* val[8]; 
    Adjuster* hue[8]; 

    sigc::connection enaConn;
	sigc::connection neutralPressedConn;

    bool lastEnabled;

public:

    HSVEqualizer ();
    virtual ~HSVEqualizer ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);
   
    void adjusterChanged (Adjuster* a, double newval);
    void enabledToggled ();
	void hsvchannelChanged ();
	
	void neutralPressed ();

};

#endif
