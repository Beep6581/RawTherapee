/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2011 Oliver Duis <oduis@oliverduis.de>
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
#include <glibmm.h>
#include "lensprofile.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include "../rtengine/lcp.h"
#include <sstream>
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

LensProfilePanel::LensProfilePanel () : Gtk::VBox(), FoldableToolPanel(this)
{
    hbLCPFile = Gtk::manage(new Gtk::HBox());

    lLCPFileHead = Gtk::manage(new Gtk::Label(M("GENERAL_FILE")));
    hbLCPFile->pack_start(*lLCPFileHead, Gtk::PACK_SHRINK, 4);

    fcbLCPFile = Gtk::manage(new MyFileChooserButton(M("TP_LENSPROFILE_LABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
    
    Gtk::FileFilter filterLCP;
    filterLCP.set_name(M("TP_LENSPROFILE_FILEDLGFILTERLCP"));
    filterLCP.add_pattern("*.lcp");
    filterLCP.add_pattern("*.LCP");
    fcbLCPFile->add_filter(filterLCP);
    
    Glib::ustring defDir=lcpStore->getDefaultCommonDirectory();
    if (!defDir.empty()) {
#ifdef WIN32
        fcbLCPFile->set_show_hidden(true);  // ProgramData is hidden on Windows
#endif
        fcbLCPFile->set_current_folder(defDir);
    }

    hbLCPFile->pack_start(*fcbLCPFile);

    btnReset = Gtk::manage(new Gtk::Button());
    btnReset->set_image (*Gtk::manage(new RTImage ("gtk-cancel.png")));
    hbLCPFile->pack_start(*btnReset, Gtk::PACK_SHRINK, 4);

    pack_start(*hbLCPFile, Gtk::PACK_SHRINK, 4);

    conLCPFile = fcbLCPFile->signal_file_set().connect( sigc::mem_fun(*this, &LensProfilePanel::onLCPFileChanged), true);
    btnReset->signal_clicked().connect( sigc::mem_fun(*this, &LensProfilePanel::onLCPFileReset), true);
}

void LensProfilePanel::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if (pp->lensProf.lcpFile.length()>0 && lcpStore->isValidLCPFileName(pp->lensProf.lcpFile))
        fcbLCPFile->set_filename (pp->lensProf.lcpFile);
    else
        fcbLCPFile->unselect_filename(fcbLCPFile->get_filename());

    lcpFileChanged = false;

    enableListener ();
}

void LensProfilePanel::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    if (lcpStore->isValidLCPFileName(fcbLCPFile->get_filename()))
        pp->lensProf.lcpFile = fcbLCPFile->get_filename();
    else
        pp->lensProf.lcpFile = "";

    if (pedited) pedited->lensProf.lcpFile = lcpFileChanged;
}

void LensProfilePanel::onLCPFileChanged()
{
    lcpFileChanged=true;
    if (listener)
        listener->panelChanged (EvLCPFile, Glib::path_get_basename(fcbLCPFile->get_filename()));
}

void LensProfilePanel::onLCPFileReset()
{
    lcpFileChanged=true;

    fcbLCPFile->unselect_filename(fcbLCPFile->get_filename());

    if (listener)
        listener->panelChanged (EvLCPFile, M("GENERAL_NONE"));
}
