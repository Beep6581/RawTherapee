/*
 *  This file is part of RawTherapee.
 */

#ifndef _CONTROLSPOTPANEL_H_
#define _CONTROLSPOTPANEL_H_

#include "../rtengine/coord.h"
#include "adjuster.h"
#include "edit.h"
#include "guiutils.h"
#include "toolpanel.h"
#include <gtkmm.h>
#include <string>

class ControlSpotPanel:
    public ToolParamBlock,
    public AdjusterListener,
    public EditSubscriber
{
public:
    ControlSpotPanel();
    void setEditProvider(EditDataProvider* provider);

private:
    // cell renderer
    void render_id (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_name (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_isvisible (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);

    void on_button_add();
    void on_button_delete();
    void on_button_rename();
    // TODO Add visibility button
    // TODO Add duplication button

    void save_ControlSpot_param();
    void load_ControlSpot_param();

    void controlspotChanged();

    void shapeChanged();
    void shapeMethodeChanged();
    void updateParamVisibility();
    void adjusterChanged(Adjuster* a, double newval);
    void disableParamlistener(bool cond);
    void setParamEditable(bool cond);

    void addControlSpotCurve(Gtk::TreeModel::Row row);
    void updateControlSpotCurve(Gtk::TreeModel::Row row);
    void deleteControlSpotCurve(Gtk::TreeModel::Row row);
    CursorShape getCursor(int objectID);
    bool mouseOver(int modifierKey);
    bool button1Pressed(int modifierKey);
    bool button1Released();
    bool drag1(int modifierKey);

    class ControlSpots:
        public Gtk::TreeModel::ColumnRecord
    {
    public:
        ControlSpots();

        Gtk::TreeModelColumn<int> id; // Control spot id
        Gtk::TreeModelColumn<Glib::ustring> name;
        Gtk::TreeModelColumn<bool> isvisible;
        Gtk::TreeModelColumn<int> curveid; // Associated curve id
        Gtk::TreeModelColumn<int> shape; // 0 = Ellipse, 1 = Rectangle
        Gtk::TreeModelColumn<int> spotMethod; // 0 = Normal, 1 = Excluding
        Gtk::TreeModelColumn<int> shapeMethod; // 0 = Independent (mouse), 1 = Symmetrical (mouse), 2 = Independent (mouse + sliders), 3 = Symmetrical (mouse + sliders)
        Gtk::TreeModelColumn<int> locX;
        Gtk::TreeModelColumn<int> locXL;
        Gtk::TreeModelColumn<int> locY;
        Gtk::TreeModelColumn<int> locYT;
        Gtk::TreeModelColumn<int> centerX;
        Gtk::TreeModelColumn<int> centerY;
        Gtk::TreeModelColumn<int> circrad;
        Gtk::TreeModelColumn<int> qualityMethod; // 0 = Standard, 1 = Enhanced, 2 = Enhanced + chroma denoise
        Gtk::TreeModelColumn<int> transit;
        Gtk::TreeModelColumn<int> thresh;
        Gtk::TreeModelColumn<int> iter;
    };

    class RenameDialog:
        public Gtk::Dialog
    {
    public:
        RenameDialog (const Glib::ustring &actualname, Gtk::Window &parent);
        Glib::ustring get_new_name();

    private:
        Gtk::Entry newname_;
    };

    ControlSpots spots_;

    // Child widgets
    Gtk::ScrolledWindow scrolledwindow_;
    Gtk::TreeView treeview_;
    Glib::RefPtr<Gtk::ListStore> treemodel_;

    Gtk::ButtonBox buttonbox_;
    Gtk::Button button_add_;
    Gtk::Button button_delete_;
    Gtk::Button button_rename_;

    MyComboBoxText* const shape_;
    sigc::connection shapeconn_;
    MyComboBoxText* const spotMethod_;
    sigc::connection spotMethodconn_;
    MyComboBoxText* const shapeMethod_;
    sigc::connection shapeMethodconn_;
    MyComboBoxText* const qualityMethod_;
    sigc::connection qualityMethodconn_;

    Adjuster* const locX_;
    Adjuster* const locXL_;
    Adjuster* const locY_;
    Adjuster* const locYT_;
    Adjuster* const centerX_;
    Adjuster* const centerY_;
    Adjuster* const circrad_;
    Adjuster* const transit_;
    Adjuster* const thresh_;
    Adjuster* const iter_;

    int lastObject_;
    rtengine::Coord* lastCoord_;
};

#endif // _CONTROLSPOTPANEL_H_
