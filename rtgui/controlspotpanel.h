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
 *  2018 Pierre Cabrera <pierre.cab@gmail.com>
 */

#ifndef _CONTROLSPOTPANEL_H_
#define _CONTROLSPOTPANEL_H_

#include "../rtengine/coord.h"
#include "adjuster.h"
#include "editcallbacks.h"
#include "guiutils.h"
#include "threadutils.h"
#include "toolpanel.h"
#include <gtkmm.h>
#include <string>
#include "thresholdadjuster.h"

class ControlSpotPanel:
    public ToolParamBlock,
    public AdjusterListener,
    public EditSubscriber,
    public FoldableToolPanel
{
public:
    /**
     * A SpotRow structure allows exchanges from and to ControlSpotClass
     */
    struct SpotRow {
        int id; // Control spot id
        Glib::ustring name;
        bool isvisible;
        int shape; // 0 = Ellipse, 1 = Rectangle
        int spotMethod; // 0 = Normal, 1 = Excluding
        int sensiexclu;
        int structexclu;
        double struc;
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
        double thresh;
        double iter;
        double balan;
        double transitweak;
        double transitgrad;
        int scopemask;
        bool avoid;
        bool laplac;
        bool deltae;
    };

    /**
     * A SpotEdited structure allows exchanges of spot panel widgets edited states from and to ControlSpotClass
     */
    struct SpotEdited {
        bool nbspot;
        bool selspot;
        bool name;
        bool isvisible;
        bool shape;
        bool spotMethod;
        bool sensiexclu;
        bool structexclu;
        bool struc;
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
        bool balan;
        bool transitweak;
        bool transitgrad;
        bool scopemask;
        bool avoid;
        bool laplac;
        bool deltae;
    };

    /**
     * An event type enumeration allows exchanges of spot panel event type from and to ControlSpotClass
     */
    enum eventType {
        None = 0,
        SpotCreation = 1,
        SpotDeletion = 2,
        SpotSelection = 3,
        SpotDuplication = 4,
        SpotAllVisibilityChanged = 5
    };

    // Constructor and management functions
    /**
     * Default constructor of ControlSpotPanel class
     */
    ControlSpotPanel();
    /**
     * Destructor of ControlSpotPanel class
     */
    ~ControlSpotPanel();
    /**
     * Implementation of setEditProvider function of toolpanel.h
     *
     * @param provider The EditDataProvider to be linked to the panel to manage curves
     */
    void setEditProvider(EditDataProvider* provider);
    /**
     * Getter of the event type raised by this panel
     *
     * @return The raised event type (refer to eventType enumeration)
     */
    int getEventType();
    /**
     * Getter of params of associated spot
     *
     * @param id The spot id to get params
     * @return A SpotRow structure containing params of associated spot
     */
    SpotRow* getSpot(const int id);
    /**
     * Get of spot id list
     *
     * @return A vector contening the list of spot id
     */
    std::vector<int>* getSpotIdList();
    /**
     * Getter of selected spot id
     *
     * @return The id of selected spot in treeview (return 0 if no selected spot)
     */
    int getSelectedSpot();
    /**
     * Setter of selected spot
     *
     * @param id The id of spot to be selected
     * @return True if a spot corresponding to the id has been selected
     */
    bool setSelectedSpot(const int id);

    // Control spot creation functions
    /**
     * Getter of available id for new spot creation
     *
     * @return An available id (i.e. max existing ones + 1)
     */
    int getNewId();
    /**
     * Add a new spot (and its associated curve)
     *
     * @param newSpot A SpotRow structure containing new spot params
     */
    void addControlSpot(SpotRow* newSpot);

    // Control spot update function
    /**
     * Update a spot (and its associated curve)
     *
     * @param spot A SpotRow structure containing spot params to update
     */
    int updateControlSpot(SpotRow* spot);

    // Control spot delete function
    /**
     * Delete a spot (and its associated curve)
     *
     * @param id The id of the spot to be deleted
     */
    void deleteControlSpot(const int id);

    // Panel widgets management functions
    /**
     * Getter of panel widgets edited states
     *
     * @return A SpotEdited structure containing the widgets edited states
     */
    SpotEdited* getEditedStates();
    /**
     * Setter of panel widgets edited states
     *
     * @param se A SpotEdited structure containing the widgets edited states to update
     */
    void setEditedStates(SpotEdited* se);
    /**
     * Implementation of setDefaults function of toolpanel.h
     *
     * @param defParams ProcParams containing default values to set to the adjusters
     * @param pedited   ParamsEdited containing default state values to set to the adjusters
     * @param id Spot id to consider to update adjusters default values and default state values
     */
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr, int id = 0);
    /**
     * Variant of setDefaults function which only update adjuster default states
     *
     * @param pedited ParamsEdited containing default states to set to the adjusters
     * @param id Spot id to consider to update adjusters default states
     */
    // void updateDefaultsStates(const ParamsEdited* pedited, int id = 0);
    /**
     * Enable or disable the interactions with panel widgets
     *
     * @param cond Condition to enable interactions
     */
    void setParamEditable(bool cond);

    // Batch mode management
    /**
     * Implementation of setBatchMode function of toolpanel.h
     *
     * @param batchMode Condition to enable batch mode
     */
    void setBatchMode(bool batchMode);

private:
    // Cell renderer
    void render_id(Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_name(Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_isvisible(Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);

    void on_button_add();
    void on_button_delete();
    void on_button_duplicate();
    void on_button_rename();
    bool on_button_visibility(GdkEventButton* event);

    bool blockTreeviewSearch(GdkEventKey* event);

    void load_ControlSpot_param();

    void controlspotChanged();

    void shapeChanged();
    void spotMethodChanged();
    void shapeMethodChanged();
    void qualityMethodChanged();
    void updateParamVisibility();
    void adjusterChanged(Adjuster* a, double newval);
    void adjusterAutoToggled(Adjuster* a, bool newval);
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop);
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight);
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop);
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight);
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR);
    void avoidChanged();
    void laplacChanged();
    void deltaeChanged();

    void disableParamlistener(bool cond);

    void addControlSpotCurve(Gtk::TreeModel::Row& row);
    void updateControlSpotCurve(const Gtk::TreeModel::Row& row);
    void deleteControlSpotCurve(Gtk::TreeModel::Row& row);
    void updateCurveOpacity(const Gtk::TreeModel::Row& selectedRow);
    CursorShape getCursor(int objectID) const;
    bool mouseOver(int modifierKey);
    bool button1Pressed(int modifierKey);
    bool button1Released();
    bool drag1(int modifierKey);

    using ToolPanel::setDefaults;

    class ControlSpots:
        public Gtk::TreeModel::ColumnRecord
    {
    public:
        ControlSpots();

        Gtk::TreeModelColumn<bool> mouseover; // Used to manage spot enlightening when mouse over
        Gtk::TreeModelColumn<int> id; // Control spot id
        Gtk::TreeModelColumn<Glib::ustring> name;
        Gtk::TreeModelColumn<bool> isvisible;
        Gtk::TreeModelColumn<int> curveid; // Associated curve id
        Gtk::TreeModelColumn<int> shape; // 0 = Ellipse, 1 = Rectangle
        Gtk::TreeModelColumn<int> spotMethod; // 0 = Normal, 1 = Excluding
        Gtk::TreeModelColumn<int> sensiexclu;
        Gtk::TreeModelColumn<int> structexclu;
        Gtk::TreeModelColumn<double> struc;
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
        Gtk::TreeModelColumn<double> thresh;
        Gtk::TreeModelColumn<double> iter;
        Gtk::TreeModelColumn<double> balan;
        Gtk::TreeModelColumn<double> transitweak;
        Gtk::TreeModelColumn<double> transitgrad;
        Gtk::TreeModelColumn<int> scopemask;
        Gtk::TreeModelColumn<bool> avoid;
        Gtk::TreeModelColumn<bool> laplac;
        Gtk::TreeModelColumn<bool> deltae;
    };

    class RenameDialog:
        public Gtk::Dialog
    {
    public:
        enum DialogButton {
            OkButton = 1,
            CancelButton = 2
        };

        RenameDialog(const Glib::ustring &actualname, Gtk::Window &parent);
        Glib::ustring get_new_name();

    private:
        Gtk::Entry* const newname_;
    };

    ControlSpots spots_;

    // Child widgets
    Gtk::ScrolledWindow* const scrolledwindow_;
    Gtk::TreeView* const treeview_;
    sigc::connection treeviewconn_;
    Glib::RefPtr<Gtk::ListStore> treemodel_;

    Gtk::Button* const button_add_;
    sigc::connection buttonaddconn_;
    Gtk::Button* const button_delete_;
    sigc::connection buttondeleteconn_;
    Gtk::Button* const button_duplicate_;
    sigc::connection buttonduplicateconn_;

    Gtk::Button* const button_rename_;
    sigc::connection buttonrenameconn_;
    Gtk::Button* const button_visibility_;
    sigc::connection buttonvisibilityconn_;

    MyComboBoxText* const shape_;
    sigc::connection shapeconn_;
    MyComboBoxText* const spotMethod_;
    sigc::connection spotMethodconn_;
    MyComboBoxText* const shapeMethod_;
    sigc::connection shapeMethodconn_;
    MyComboBoxText* const qualityMethod_;
    sigc::connection qualityMethodconn_;

    Adjuster* const sensiexclu_;
    Adjuster* const structexclu_;
    Adjuster* const struc_;
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
    Adjuster* const balan_;
    Adjuster* const transitweak_;
    Adjuster* const transitgrad_;
    Adjuster* const scopemask_;

    Gtk::CheckButton* const avoid_;
    sigc::connection avoidConn_;
    Gtk::CheckButton* const laplac_;
    sigc::connection laplacConn_;
    Gtk::CheckButton* const deltae_;
    sigc::connection deltaeConn_;

    // Internal variables
    int lastObject_;
    rtengine::Coord lastCoord_;
    bool nbSpotChanged_;
    bool selSpotChanged_;
    bool nameChanged_;
    bool visibilityChanged_;
    int eventType; // 0 = No event, 1 = Spot creation event, 2 = Spot deletion event, 3 = Spot selection event, 4 = Spot duplication event
    Gtk::Frame* const excluFrame;

    // Row background color
    Gdk::RGBA colorMouseover, colorNominal, colorMouseovertext;

    // Treeview mutex
    MyMutex mTreeview;
};

#endif // _CONTROLSPOTPANEL_H_
