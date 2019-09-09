/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 *  2017 Jacques Desmis <jdesmis@gmail.com>
 *  2019 Pierre Cabrera <pierre.cab@gmail.com>
 */
#ifndef _LOCALLAB_H_
#define _LOCALLAB_H_

#include "controlspotpanel.h"
#include "locallabtools.h"

/* ==== LocallabToolListListener ==== */
class LocallabToolList;
class LocallabToolListListener
{
public:
    LocallabToolListListener() {};
    virtual ~LocallabToolListListener() {};

    virtual void locallabToolToAdd(const Glib::ustring &toolname) = 0;
};

/* ==== LocallabToolList ==== */
class LocallabToolList:
    public Gtk::VBox
{
private:
    // Tree model to manage ComboBox rows
    class ToolRow:
        public Gtk::TreeModel::ColumnRecord
    {
    public:
        Gtk::TreeModelColumn<int> id;
        Gtk::TreeModelColumn<Glib::ustring> name;

        ToolRow()
        {
            add(id);
            add(name);
        }
    };

    // Tool list GUI widgets
    MyComboBox* const list;
    sigc::connection listConn;
    ToolRow toolRow;
    Glib::RefPtr<Gtk::ListStore> listTreeModel;

    // Tool list listener
    LocallabToolListListener* listListener;

public:
    LocallabToolList();

    // Setter for tool list listener
    void setLocallabToolListListener(LocallabToolListListener* ltll)
    {
        listListener = ltll;
    }

    // Tool list management function
    void addToolRow(const Glib::ustring &toolname, const int id);
    void removeToolRow(const Glib::ustring &toolname);
    void removeAllTool();

private:
    // Tool list event management function
    void toolRowSelected();
};

/* ==== Locallab ==== */
class Locallab :
    public ToolParamBlock,
    public FoldableToolPanel,
    public rtengine::LocallabListener,
    public LocallabToolListener,
    public LocallabToolListListener
{
private:
    // Spot control panel widget
    ControlSpotPanel* const expsettings;

    // Tool list widget
    LocallabToolList* const toollist;

    // Locallab tool widgets
    LocallabColor* const expcolor;
    LocallabExposure* const expexpose;
    LocallabShadow* const expshadhigh;
    LocallabVibrance* const expvibrance;
    LocallabSoft* const expsoft;
    LocallabBlur* const expblur;
    LocallabTone* const exptonemap;
    LocallabRetinex* const expreti;
    LocallabSharp* const expsharp;
    LocallabContrast* const expcontrast;
    LocallabCBDL* const expcbdl;
    LocallabDenoise* const expdenoi;

    std::vector<LocallabTool*> locallabTools;

    // Locallab tools mask background management data
    std::vector<locallabRef> maskBackRef;

public:
    Locallab();

    // FoldableToolPanel management functions
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void setListener(ToolPanelListener* tpl);

    // Locallab tools mask background management function
    void refChanged(const std::vector<locallabRef> &ref, int selspot);

    // Mask visibility management functions
    struct llMaskVisibility {
        int colorMask;
        int expMask;
        int SHMask;
        int cbMask;
        int retiMask;
        int softMask;
        int tmMask;
    };

    void resetMaskVisibility();
    llMaskVisibility getMaskVisibility() const;

    // EditProvider management function
    void setEditProvider(EditDataProvider* provider);
    void subscribe();
    void unsubscribe();

    // FoldableToolPanel event function
    void enabledChanged();

    // Curve management function
    void autoOpenCurve();

    // Locallab tools expanders management functions
    void foldAllButOne(LocallabTool* except);

private:
    // Locallab tools management functions
    void addTool(Gtk::Box* where, LocallabTool* tool);

    // Locallab tools management functions
    void openAllTools();

    // Locallab GUI management function
    void setParamEditable(bool cond);

    // LocallabToolListener function
    void resetOtherMaskView(LocallabTool* current);
    void toolRemoved(LocallabTool* current);

    // LocallabToolListListener function
    void locallabToolToAdd(const Glib::ustring &toolname);
};

#endif
