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
#include "../rtengine/rtlensfun.h"
#include <map>
#include <set>

using namespace rtengine;
using namespace rtengine::procparams;

LensProfilePanel::LFDbHelper *LensProfilePanel::lf(nullptr);

LensProfilePanel::LensProfilePanel () :
    FoldableToolPanel(this, "lensprof", M("TP_LENSPROFILE_LABEL")),
    lcpFileChanged(false),
    useDistChanged(false),
    useVignChanged(false),
    useCAChanged(false),
    isRaw(true),
    metadata(nullptr),
    lensgeomLcpFill(nullptr),
    useLensfunChanged(false),
    lensfunAutoChanged(false),
    lensfunCameraChanged(false),
    lensfunLensChanged(false)
{
    if (!lf) {
        lf = new LFDbHelper();
    }
    
    corrOff = Gtk::manage(new Gtk::RadioButton(M("LENSPROFILE_CORRECTION_OFF")));
    pack_start(*corrOff);

    corrGroup = corrOff->get_group();
    
    corrLensfunAuto = Gtk::manage(new Gtk::RadioButton(corrGroup, M("LENSPROFILE_CORRECTION_AUTOMATCH")));
    pack_start(*corrLensfunAuto);
    
    corrLensfunManual = Gtk::manage(new Gtk::RadioButton(corrGroup, M("LENSPROFILE_CORRECTION_MANUAL")));
    pack_start(*corrLensfunManual);

    lensfunCameras = Gtk::manage(new MyComboBox());
    lensfunCameras->set_model(lf->lensfunCameraModel);
    lensfunCameras->pack_start(lf->lensfunModelCam.model);
    lensfunLenses = Gtk::manage(new MyComboBox());
    lensfunLenses->set_model(lf->lensfunLensModel);
    lensfunLenses->pack_start(lf->lensfunModelLens.prettylens);
    
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("EXIFFILTER_CAMERA"))), Gtk::PACK_SHRINK, 4);
    hb->pack_start(*lensfunCameras);
    pack_start(*hb);

    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("EXIFFILTER_LENS"))), Gtk::PACK_SHRINK, 4);
    hb->pack_start(*lensfunLenses);
    pack_start(*hb);

    corrLcpFile = Gtk::manage(new Gtk::RadioButton(corrGroup));
    hbLCPFile = Gtk::manage(new Gtk::HBox());
    hbLCPFile->pack_start(*corrLcpFile, Gtk::PACK_SHRINK);

    lLCPFileHead = Gtk::manage(new Gtk::Label(M("LENSPROFILE_CORRECTION_LCPFILE")));
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
    } else if (!options.lastLensProfileDir.empty()) {
        fcbLCPFile->set_current_folder(options.lastLensProfileDir);
    }
    bindCurrentFolder(*fcbLCPFile, options.lastLensProfileDir);

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

    lensfunCameras->signal_changed().connect(sigc::mem_fun(*this, &LensProfilePanel::onLensfunCameraChanged));
    lensfunLenses->signal_changed().connect(sigc::mem_fun(*this, &LensProfilePanel::onLensfunLensChanged));
    corrOff->signal_toggled().connect(sigc::mem_fun(*this, &LensProfilePanel::onCorrModeChanged));
    corrLensfunAuto->signal_toggled().connect(sigc::mem_fun(*this, &LensProfilePanel::onCorrModeChanged));
    corrLensfunManual->signal_toggled().connect(sigc::mem_fun(*this, &LensProfilePanel::onCorrModeChanged));
    corrLcpFile->signal_toggled().connect(sigc::mem_fun(*this, &LensProfilePanel::onCorrModeChanged));

    allowFocusDep = true;
}

void LensProfilePanel::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();
    conUseDist.block(true);

    corrLensfunAuto->set_sensitive(true);
    
    if (pp->lensProf.useLensfun) {
        if (pp->lensProf.lfAutoMatch) {
            corrLensfunAuto->set_active(true);
        } else {
            corrLensfunManual->set_active(true);
        }
    } else if (!pp->lensProf.lcpFile.empty() && lcpStore->isValidLCPFileName(pp->lensProf.lcpFile)) {
        corrLcpFile->set_active(true);
        fcbLCPFile->set_filename (pp->lensProf.lcpFile);
        updateDisabled(true);
    } else {
        Glib::ustring fname = fcbLCPFile->get_filename();

        if (!pp->lensProf.lcpFile.empty()) {
            fcbLCPFile->unselect_filename(fname);
        } else {
            Glib::ustring lastFolder = fcbLCPFile->get_current_folder();
            fcbLCPFile->set_current_folder(lastFolder);
            fcbLCPFile->set_filename(lastFolder + "/.");
            bindCurrentFolder(*fcbLCPFile, options.lastLensProfileDir);
        }

        updateDisabled(false);

        corrOff->set_active(true);
    }

    ckbUseDist->set_active (pp->lensProf.useDist);
    ckbUseVign->set_active (pp->lensProf.useVign && isRaw);
    ckbUseCA->set_active   (pp->lensProf.useCA && isRaw);

    const LFDatabase *db = LFDatabase::getInstance();
    LFCamera c;
    LFLens l;
    if (metadata) {
        c = db->findCamera(metadata->getMake(), metadata->getModel());
        l = db->findLens(c, metadata->getLens());
    }
    
    if (!setLensfunCamera(pp->lensProf.lfCameraMake, pp->lensProf.lfCameraModel) && pp->lensProf.lfAutoMatch) {
        setLensfunCamera(c.getMake(), c.getModel());
    }
    if (!setLensfunLens(pp->lensProf.lfLens) && pp->lensProf.lfAutoMatch) {
        setLensfunLens(l.getLens());
    }
    
    lcpFileChanged = useDistChanged = useVignChanged = useCAChanged = false;
    useLensfunChanged = lensfunAutoChanged = lensfunCameraChanged = lensfunLensChanged = false;

    if (!checkLensfunCanCorrect(true)) {
        if (corrLensfunAuto->get_active()) {
            corrOff->set_active(true);
        }
        corrLensfunAuto->set_sensitive(false);
    }

    if (corrLensfunManual->get_active() && !checkLensfunCanCorrect(false)) {
        corrOff->set_active(true);
    }
    // if (metadata) {
    //     std::unique_ptr<LFModifier> mod(LFDatabase::findModifier(pp->lensProf, metadata, 100, 100, pp->coarse, -1));
    //     if (!mod) {
    //         if (pp->lensProf.useLensfun) {
    //             corrOff->set_active(true);
    //         }
    //         corrLensfunAuto->set_sensitive(false);
    //     }
    // }

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
    metadata = pMeta;
}

void LensProfilePanel::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    if (corrLcpFile->get_active() && lcpStore->isValidLCPFileName(fcbLCPFile->get_filename())) {
        pp->lensProf.lcpFile = fcbLCPFile->get_filename();
    } else {
        pp->lensProf.lcpFile = "";
    }

    pp->lensProf.useDist = ckbUseDist->get_active();
    pp->lensProf.useVign = ckbUseVign->get_active();
    pp->lensProf.useCA   = ckbUseCA->get_active();

    pp->lensProf.useLensfun = corrLensfunAuto->get_active() || corrLensfunManual->get_active();
    pp->lensProf.lfAutoMatch = corrLensfunAuto->get_active();
    auto itc = lensfunCameras->get_active();
    if (itc) {
        pp->lensProf.lfCameraMake = (*itc)[lf->lensfunModelCam.make];
        pp->lensProf.lfCameraModel = (*itc)[lf->lensfunModelCam.model];
    } else {
        pp->lensProf.lfCameraMake = "";
        pp->lensProf.lfCameraModel = "";
    }
    auto itl = lensfunLenses->get_active();
    if (itl) {
        pp->lensProf.lfLens = (*itl)[lf->lensfunModelLens.lens];
    } else {
        pp->lensProf.lfLens = "";
    }

    if (pedited) {
        pedited->lensProf.lcpFile = lcpFileChanged;
        pedited->lensProf.useDist = useDistChanged;
        pedited->lensProf.useVign = useVignChanged;
        pedited->lensProf.useCA   = useCAChanged;
        pedited->lensProf.useLensfun = useLensfunChanged;
        pedited->lensProf.lfAutoMatch = lensfunAutoChanged;
        pedited->lensProf.lfCameraMake = lensfunCameraChanged;
        pedited->lensProf.lfCameraModel = lensfunCameraChanged;
        pedited->lensProf.lfLens = lensfunLensChanged;
    }
}

void LensProfilePanel::onLCPFileChanged()
{
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

void LensProfilePanel::setBatchMode(bool yes)
{
    FoldableToolPanel::setBatchMode(yes);
}


bool LensProfilePanel::setLensfunCamera(const Glib::ustring &make, const Glib::ustring &model)
{
    if (!make.empty() && !model.empty()) {
        auto it = lensfunCameras->get_active();
        if (it && (*it)[lf->lensfunModelCam.make] == make && (*it)[lf->lensfunModelCam.model] == model) {
            return true;
        }
        
        // search for the active row
        for (auto row : lf->lensfunCameraModel->children()) {
            if (row[lf->lensfunModelCam.make] == make) {
                auto &c = row.children();
                for (auto it = c.begin(), end = c.end(); it != end; ++it) {
                    auto &childrow = *it;
                    if (childrow[lf->lensfunModelCam.model] == model) {
                        lensfunCameras->set_active(it);
                        return true;
                    }
                }
                break;
            }
        }
    }
    lensfunCameras->set_active(-1);
    return false;
}


bool LensProfilePanel::setLensfunLens(const Glib::ustring &lens)
{
    if (!lens.empty()) {
        auto it = lensfunLenses->get_active();
        if (it && (*it)[lf->lensfunModelLens.lens] == lens) {
            return true;
        }
        
        for (auto row : lf->lensfunLensModel->children()) {
            if (lens.find(row[lf->lensfunModelLens.lens]) == 0) {
                auto &c = row.children();
                for (auto it = c.begin(), end = c.end(); it != end; ++it) {
                    auto &childrow = *it;
                    if (childrow[lf->lensfunModelLens.lens] == lens) {
                        lensfunLenses->set_active(it);
                        return true;
                    }
                }
                break;
            }
        }
    }
    lensfunLenses->set_active(-1);
    return false;
}



void LensProfilePanel::onLensfunCameraChanged()
{
    auto iter = lensfunCameras->get_active();

    if (iter) {
        lensfunCameraChanged = true;

        if (listener) {
            Glib::ustring name = (*iter)[lf->lensfunModelCam.model];
            listener->panelChanged(EvLensCorrLensfunCamera, name);
        }
    }
}


void LensProfilePanel::onLensfunLensChanged()
{
    auto iter = lensfunLenses->get_active();

    if (iter) {
        lensfunLensChanged = true;

        if (listener) {
            Glib::ustring name = (*iter)[lf->lensfunModelLens.lens];
            listener->panelChanged(EvLensCorrLensfunLens, name);
        }
    }
}


void LensProfilePanel::onCorrModeChanged()
{
    Glib::ustring mode;

    if (corrOff->get_active()) {
        useLensfunChanged = true;
        lensfunAutoChanged = true;
        lcpFileChanged = true;
        
        lensfunCameras->set_sensitive(false);
        lensfunLenses->set_sensitive(false);
        ckbUseDist->set_sensitive(false);
        ckbUseVign->set_sensitive(false);
        ckbUseCA->set_sensitive(false);
        
        mode = M("LENSPROFILE_CORRECTION_OFF");
    } else if (corrLensfunAuto->get_active()) {
        useLensfunChanged = true;
        lensfunAutoChanged = true;
        lcpFileChanged = true;
        useDistChanged = true;
        useVignChanged = true;

        lensfunCameras->set_sensitive(false);
        lensfunLenses->set_sensitive(false);

        ckbUseDist->set_sensitive(true);
        ckbUseVign->set_sensitive(true);
        ckbUseCA->set_sensitive(false);

        if (metadata) {
            bool b = disableListener();
            const LFDatabase *db = LFDatabase::getInstance();
            LFCamera c = db->findCamera(metadata->getMake(), metadata->getModel());
            LFLens l = db->findLens(c, metadata->getLens());
            setLensfunCamera(c.getMake(), c.getModel());
            setLensfunLens(l.getLens());
            if (b) {
                enableListener();
            }
        }

        mode = M("LENSPROFILE_CORRECTION_AUTOMATCH");
    } else if (corrLensfunManual->get_active()) {
        useLensfunChanged = true;
        lensfunAutoChanged = true;
        lcpFileChanged = true;
        useDistChanged = true;
        useVignChanged = true;

        lensfunCameras->set_sensitive(true);
        lensfunLenses->set_sensitive(true);

        ckbUseDist->set_sensitive(true);
        ckbUseVign->set_sensitive(true);
        ckbUseCA->set_sensitive(false);

        mode = M("LENSPROFILE_CORRECTION_MANUAL");
    } else if (corrLcpFile->get_active()) {
        useLensfunChanged = true;
        lensfunAutoChanged = true;
        lcpFileChanged = true;
        useDistChanged = true;
        useVignChanged = true;

        lensfunCameras->set_sensitive(false);
        lensfunLenses->set_sensitive(false);
        updateDisabled(true);

        mode = M("LENSPROFILE_CORRECTION_LCPFILE");
    }
    
    if (listener) {
        listener->panelChanged(EvLensCorrMode, mode);
    }
}


bool LensProfilePanel::checkLensfunCanCorrect(bool automatch)
{
    if (!metadata) {
        return false;
    }
    rtengine::procparams::ProcParams lpp;
    write(&lpp);
    lpp.lensProf.lfAutoMatch = automatch;
    std::unique_ptr<LFModifier> mod(LFDatabase::findModifier(lpp.lensProf, metadata, 100, 100, lpp.coarse, -1));
    return mod.get() != nullptr;
}


//-----------------------------------------------------------------------------
// LFDbHelper
//-----------------------------------------------------------------------------

LensProfilePanel::LFDbHelper::LFDbHelper()
{
    lensfunCameraModel = Gtk::TreeStore::create(lensfunModelCam);
    lensfunLensModel = Gtk::TreeStore::create(lensfunModelLens);

    fillLensfunCameras();
    fillLensfunLenses();
}

void LensProfilePanel::LFDbHelper::fillLensfunCameras()
{
    std::map<Glib::ustring, std::set<Glib::ustring>> camnames;
    auto camlist = LFDatabase::getInstance()->getCameras();
    for (auto &c : camlist) {
        camnames[c.getMake()].insert(c.getModel());
    }
    for (auto &p : camnames) {
        Gtk::TreeModel::Row row = *(lensfunCameraModel->append());
        row[lensfunModelCam.make] = p.first;
        row[lensfunModelCam.model] = p.first;
        for (auto &c : p.second) {
            Gtk::TreeModel::Row child = *(lensfunCameraModel->append(row.children()));
            child[lensfunModelCam.make] = p.first;
            child[lensfunModelCam.model] = c;
        }
    }
}


void LensProfilePanel::LFDbHelper::fillLensfunLenses()
{
    std::map<Glib::ustring, std::set<Glib::ustring>> lenses;
    auto lenslist = LFDatabase::getInstance()->getLenses();
    for (auto &l : lenslist) {
        auto name = l.getLens();
        auto make = l.getMake();
        lenses[make].insert(name);
    }
    for (auto &p : lenses) {
        Gtk::TreeModel::Row row = *(lensfunLensModel->append());
        row[lensfunModelLens.lens] = p.first;
        row[lensfunModelLens.prettylens] = p.first;
        for (auto &c : p.second) {
            Gtk::TreeModel::Row child = *(lensfunLensModel->append(row.children()));
            child[lensfunModelLens.lens] = c;
            if (c.find(p.first, p.first.size()+1) == p.first.size()+1) {
                child[lensfunModelLens.prettylens] = c.substr(p.first.size()+1);
            } else {
                child[lensfunModelLens.prettylens] = c;
            }
        }
    }
}
