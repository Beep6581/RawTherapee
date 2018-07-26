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
    public EditSubscriber,
    public FoldableToolPanel
{
public:
    /** A SpotRow structure allows exchanges from and to ControlSpotClass */
    struct SpotRow {
        int id; // Control spot id
        Glib::ustring name;
        bool isvisible;
        int shape; // 0 = Ellipse, 1 = Rectangle
        int spotMethod; // 0 = Normal, 1 = Excluding
        int shapeMethod; // 0 = Independent (mouse), 1 = Symmetrical (mouse), 2 = Independent (mouse + sliders), 3 = Symmetrical (mouse + sliders)
        int locX;
        int locXL;
        int locY;
        int locYT;
        int centerX;
        int centerY;
        int circrad;
        int qualityMethod; // 0 = Standard, 1 = Enhanced, 2 = Enhanced + chroma denoise
        int transit;
        int thresh;
        int iter;
    };

    /** A SpotEdited structure allows exchanges of spot panel widgets edited states from and to ControlSpotClass */
    struct SpotEdited {
        bool addbutton;
        bool deletebutton;
        bool treeview;
        bool name;
        bool isvisible;
        bool shape;
        bool spotMethod;
        bool shapeMethod;
        bool locX;
        bool locXL;
        bool locY;
        bool locYT;
        bool centerX;
        bool centerY;
        bool circrad;
        bool qualityMethod;
        bool transit;
        bool thresh;
        bool iter;
    };

    // Constructor and management functions
    ControlSpotPanel();
    void setEditProvider(EditDataProvider* provider);
    int getEventType();
    SpotRow* getSpot(int id);
    std::vector<int>* getSpotIdList();
    int getSelectedSpot();
    void setSelectedSpot(int id);

    // Control spot creation functions
    int getNewId();
    void addControlSpot(SpotRow* newSpot);

    // Control spot update function
    int updateControlSpot(SpotRow* spot);

    // Control spot delete function
    void deleteControlSpot(int id);

    // Panel widgets edited states management functions
    SpotEdited* getEditedStates();
    void setEditedStates(SpotEdited* se);

private:
    // cell renderer
    void render_id(Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_name(Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_isvisible(Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);

    void on_button_add();
    void on_button_delete();
    void on_button_rename();
    // TODO Add visibility button
    // TODO Add duplication button

    void save_ControlSpot_param();
    void load_ControlSpot_param();

    void controlspotChanged();

    void shapeChanged();
    void spotMethodChanged();
    void shapeMethodChanged();
    void qualityMethodChanged();
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
        RenameDialog(const Glib::ustring &actualname, Gtk::Window &parent);
        Glib::ustring get_new_name();

    private:
        Gtk::Entry newname_;
    };

    ControlSpots spots_;

    // Child widgets
    Gtk::ScrolledWindow scrolledwindow_;
    Gtk::TreeView treeview_;
    sigc::connection treeviewconn_;
    Glib::RefPtr<Gtk::ListStore> treemodel_;

    Gtk::ButtonBox buttonbox_;
    Gtk::Button button_add_;
    sigc::connection buttonaddconn_;
    Gtk::Button button_delete_;
    sigc::connection buttondeleteconn_;
    Gtk::Button button_rename_;
    sigc::connection buttonrenameconn_;

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

    int eventType; // 0 = No event, 1 = Spot creation event, 2 = Spot deletion event, 3 = Spot selection event
};

#endif // _CONTROLSPOTPANEL_H_
