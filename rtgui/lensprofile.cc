/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2012 Oliver Duis <oduis@oliverduis.de>
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
#include "../rtengine/lcp.h"
#include <sstream>
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

LensProfilePanel::LensProfilePanel () : FoldableToolPanel(this, "lensprof", M("TP_LENSPROFILE_LABEL")), isRaw(true)
{
    hbLCPFile = Gtk::manage(new Gtk::HBox());

    lLCPFileHead = Gtk::manage(new Gtk::Label(M("GENERAL_FILE")));
    hbLCPFile->pack_start(*lLCPFileHead, Gtk::PACK_SHRINK, 4);

    fcbLCPFile = Gtk::manage(new MyFileChooserButton(M("TP_LENSPROFILE_LABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));

    Glib::RefPtr<Gtk::FileFilter> filterLCP = Gtk::FileFilter::create();
    filterLCP->set_name(M("FILECHOOSER_FILTER_LCP"));
    filterLCP->add_pattern("*.lcp");
    filterLCP->add_pattern("*.LCP");
    fcbLCPFile->add_filter(filterLCP);

    Glib::ustring defDir = lcpStore->getDefaultCommonDirectory();

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

    ckbUseDist = Gtk::manage (new Gtk::CheckButton (M("TP_LENSPROFILE_USEDIST")));
    pack_start (*ckbUseDist, Gtk::PACK_SHRINK, 4);
    ckbUseVign = Gtk::manage (new Gtk::CheckButton (M("TP_LENSPROFILE_USEVIGN")));
    pack_start (*ckbUseVign, Gtk::PACK_SHRINK, 4);
    ckbUseCA = Gtk::manage (new Gtk::CheckButton (M("TP_LENSPROFILE_USECA")));
    pack_start (*ckbUseCA, Gtk::PACK_SHRINK, 4);

    conLCPFile = fcbLCPFile->signal_file_set().connect( sigc::mem_fun(*this, &LensProfilePanel::onLCPFileChanged), true);
    btnReset->signal_clicked().connect( sigc::mem_fun(*this, &LensProfilePanel::onLCPFileReset), true);
    conUseDist = ckbUseDist->signal_toggled().connect( sigc::mem_fun(*this, &LensProfilePanel::onUseDistChanged) );
    ckbUseVign->signal_toggled().connect( sigc::mem_fun(*this, &LensProfilePanel::onUseVignChanged) );
    ckbUseCA->signal_toggled().connect( sigc::mem_fun(*this, &LensProfilePanel::onUseCAChanged) );

    allowFocusDep = true;
}

void LensProfilePanel::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    conUseDist.block(true);

    if (!pp->lensProf.lcpFile.empty() && lcpStore->isValidLCPFileName(pp->lensProf.lcpFile)) {
        fcbLCPFile->set_filename (pp->lensProf.lcpFile);
        updateDisabled(true);
    } else {
        Glib::ustring fname = fcbLCPFile->get_filename();

        if (!pp->lensProf.lcpFile.empty()) {
            fcbLCPFile->unselect_filename(fname);
        } else {
            Glib::ustring lastFolder = fcbLCPFile->get_current_folder();
            fcbLCPFile->set_filename("");
            fcbLCPFile->set_current_folder(lastFolder);
        }

        updateDisabled(false);
    }

    ckbUseDist->set_active (pp->lensProf.useDist);
    ckbUseVign->set_active (pp->lensProf.useVign && isRaw);
    ckbUseCA->set_active   (pp->lensProf.useCA && isRaw);

    lcpFileChanged = useDistChanged = useVignChanged = useCAChanged = false;

    enableListener ();
    conUseDist.block(false);
}

void LensProfilePanel::setRawMeta(bool raw, const rtengine::ImageMetaData* pMeta)
{
    if (!raw || pMeta->getFocusDist() <= 0) {
        disableListener();

        // CA is very focus layer dependend, otherwise it might even worsen things
        allowFocusDep = false;
        ckbUseCA->set_active(false);
        ckbUseCA->set_sensitive(false);

        enableListener();
    }

    isRaw = raw;
}

void LensProfilePanel::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    if (lcpStore->isValidLCPFileName(fcbLCPFile->get_filename())) {
        pp->lensProf.lcpFile = fcbLCPFile->get_filename();
    } else {
        pp->lensProf.lcpFile = "";
    }

    pp->lensProf.useDist = ckbUseDist->get_active();
    pp->lensProf.useVign = ckbUseVign->get_active();
    pp->lensProf.useCA   = ckbUseCA->get_active();

    if (pedited) {
        pedited->lensProf.lcpFile = lcpFileChanged;
        pedited->lensProf.useDist = useDistChanged;
        pedited->lensProf.useVign = useVignChanged;
        pedited->lensProf.useCA   = useCAChanged;
    }
}

void LensProfilePanel::onLCPFileChanged()
{

    // Disable Auto-Fill when enabling LCP Distortion Correction, #1791
    lensgeomLcpFill->disableAutoFillIfActive();

    lcpFileChanged = true;
    updateDisabled(lcpStore->isValidLCPFileName(fcbLCPFile->get_filename()));

    if (listener) {
        listener->panelChanged (EvLCPFile, Glib::path_get_basename(fcbLCPFile->get_filename()));
    }
}

void LensProfilePanel::onLCPFileReset()
{
    lcpFileChanged = true;

    fcbLCPFile->unselect_filename(fcbLCPFile->get_filename());
    updateDisabled(false);

    if (listener) {
        listener->panelChanged (EvLCPFile, M("GENERAL_NONE"));
    }
}

void LensProfilePanel::onUseDistChanged()
{

    // Disable Auto-Fill when enabling LCP Distortion Correction, #1791
    if (ckbUseDist->get_active()) {
        lensgeomLcpFill->disableAutoFillIfActive();
    }

    useDistChanged = true;

    if (listener) {
        listener->panelChanged (EvLCPUseDist, ckbUseDist->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
void LensProfilePanel::onUseVignChanged()
{
    useVignChanged = true;

    if (listener) {
        listener->panelChanged (EvLCPUseVign, ckbUseVign->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
void LensProfilePanel::onUseCAChanged()
{
    useCAChanged = true;

    if (listener) {
        listener->panelChanged (EvLCPUseCA, ckbUseCA->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void LensProfilePanel::updateDisabled(bool enable)
{
    ckbUseDist->set_sensitive(enable);
    ckbUseVign->set_sensitive(enable && isRaw);
    ckbUseCA->set_sensitive(enable && allowFocusDep);
}
