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
#ifndef _WB_H_
#define _WB_H_

#include <gtkmm.h>
#include "toolpanel.h"
#include "adjuster.h"
#include "guiutils.h"
#include "wbprovider.h"
#include "../rtengine/procparams.h"

class SpotWBListener {

    public: 
        virtual void spotWBRequested (int size) {}
};

class WhiteBalance : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel {

    enum WB_LabelType {
        WBLT_GUI,
        WBLT_PP
    };

  protected:
    class MethodColumns : public Gtk::TreeModel::ColumnRecord {
      public:
        Gtk::TreeModelColumn< Glib::RefPtr<Gdk::Pixbuf> > colIcon;
        Gtk::TreeModelColumn<Glib::ustring> colLabel;
        Gtk::TreeModelColumn<int> colId;
        MethodColumns() { add(colIcon); add(colLabel); add(colId); }
    };

    static Glib::RefPtr<Gdk::Pixbuf> wbPixbufs[rtengine::procparams::WBT_CUSTOM+1];
    Glib::RefPtr<Gtk::TreeStore> refTreeModel;
    MethodColumns methodColumns;
    MyComboBox* method;
    MyComboBoxText* spotsize;
    Adjuster* temp;
    Adjuster* green;
    Gtk::Button* spotbutton;
    int opt;
    double nextTemp;
    double nextGreen;
    WBProvider *wbp;
    SpotWBListener* wblistener;
    sigc::connection methconn;
    int custom_temp;
    double* custom_green;
    void cache_customWB    (int temp, double green); //cache custom WB setting to allow its recall
    void cache_customTemp  (int temp);               //cache Temperature only to allow its recall
    void cache_customGreen (double green);           //cache Green only to allow its recall
    int  setActiveMethod   (Glib::ustring label);
    int _setActiveMethod   (Glib::ustring &label, Gtk::TreeModel::Children &children);

    Gtk::TreeModel::Row            getActiveMethod ();
    int                            findWBEntryId   (Glib::ustring label, enum WB_LabelType lblType=WBLT_GUI);
    rtengine::procparams::WBEntry* findWBEntry     (Glib::ustring label, enum WB_LabelType lblType=WBLT_GUI);

  public:

    WhiteBalance ();
    ~WhiteBalance ();

    static void init    ();
    static void cleanup ();
    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);
    
    
    void optChanged ();
    void spotPressed ();
    void spotSizeChanged ();
    void adjusterChanged (Adjuster* a, double newval);
    int  getSize (); 
    void setWBProvider (WBProvider* p) { wbp = p; }
    void setSpotWBListener (SpotWBListener* l) { wblistener = l; }
    void setWB (int temp, double green);
    
    void setAdjusterBehavior (bool tempadd, bool greenadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
};

#endif
