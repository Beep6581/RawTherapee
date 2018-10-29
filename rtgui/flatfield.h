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
#ifndef _FLATFIELD_H_
#define _FLATFIELD_H_

#include <memory>
#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "../rtengine/rawimage.h"
#include "guiutils.h"

class FFProvider
{
public:
    virtual ~FFProvider() {}
    virtual rtengine::RawImage* getFF() = 0;
    virtual Glib::ustring GetCurrentImageFilePath() = 0;
    // add other info here
};

class FlatField : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public rtengine::FlatFieldAutoClipListener
{

protected:

    MyFileChooserButton *flatFieldFile;
    Gtk::Label *ffLabel;
    Gtk::Label *ffInfo;
    Gtk::Button *flatFieldFileReset;
    Gtk::CheckButton* flatFieldAutoSelect;
    Adjuster* flatFieldClipControl;
    Adjuster* flatFieldBlurRadius;
    MyComboBoxText* flatFieldBlurType;
    Gtk::HBox *hbff;
    bool ffChanged;
    bool lastFFAutoSelect;
    bool lastFFAutoClipCtrl;
    FFProvider *ffp;
    sigc::connection flatFieldFileconn, flatFieldAutoSelectconn, flatFieldBlurTypeconn;
    Glib::ustring lastShortcutPath;
    bool b_filter_asCurrent;
    bool israw;

    IdleRegister idle_register;
public:

    FlatField ();
    ~FlatField ();

    void read                (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write               (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setBatchMode        (bool batchMode);
    void setAdjusterBehavior (bool clipctrladd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void setDefaults         (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void adjusterChanged            (Adjuster* a, double newval);
    void adjusterAutoToggled        (Adjuster* a, bool newval);
    void flatFieldFileChanged       ();
    void flatFieldFile_Reset        ();
    void flatFieldAutoSelectChanged ();
    void flatFieldBlurTypeChanged   ();
    void setShortcutPath (const Glib::ustring& path);
    void setFFProvider              (FFProvider* p)
    {
        ffp = p;
    };
    void flatFieldAutoClipValueChanged(int n = 0);
};

#endif
