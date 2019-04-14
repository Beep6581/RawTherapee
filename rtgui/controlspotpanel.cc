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

#include "../rtengine/rt_math.h"
#include "controlspotpanel.h"
#include "multilangmgr.h"
#include <iomanip>
#include "editwidgets.h"
#include "options.h"

using namespace rtengine;
extern Options options;

//-----------------------------------------------------------------------------
// ControlSpotPanel
//-----------------------------------------------------------------------------

ControlSpotPanel::ControlSpotPanel():
    EditSubscriber(ET_OBJECTS),
    FoldableToolPanel(this, "controlspotpanel", M("TP_LOCALLAB_SETTINGS")),

    button_add_(M("TP_LOCALLAB_BUTTON_ADD")),
    button_delete_(M("TP_LOCALLAB_BUTTON_DEL")),
    button_duplicate_(M("TP_LOCALLAB_BUTTON_DUPL")),

    button_rename_(M("TP_LOCALLAB_BUTTON_REN")),
    button_visibility_(M("TP_LOCALLAB_BUTTON_VIS")),

    shape_(Gtk::manage(new MyComboBoxText())),
    spotMethod_(Gtk::manage(new MyComboBoxText())),
    shapeMethod_(Gtk::manage(new MyComboBoxText())),
    qualityMethod_(Gtk::manage(new MyComboBoxText())),

    sensiexclu_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIEXCLU"), 0, 100, 1, 12))),
    structexclu_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUCCOL"), 0, 100, 1, 0))),
    
    struc_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_THRES"), 1.0, 12.0, 0.1, 4.0))),
    locX_(Gtk::manage(new Adjuster(M("TP_LOCAL_WIDTH"), 2, 2250, 1, 250))),
    locXL_(Gtk::manage(new Adjuster(M("TP_LOCAL_WIDTH_L"), 2, 2250, 1, 250))),
    locY_(Gtk::manage(new Adjuster(M("TP_LOCAL_HEIGHT"), 2, 2250, 1, 250))),
    locYT_(Gtk::manage(new Adjuster(M("TP_LOCAL_HEIGHT_T"), 2, 2250, 1, 250))),
    centerX_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTER_X"), -1000, 1000, 1, 0))),
    centerY_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTER_Y"), -1000, 1000, 1, 0))),
    circrad_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CIRCRADIUS"), 2, 150, 1, 18))),
    transit_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_TRANSITVALUE"), 5, 95, 1, 60))),
    thresh_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_THRESDELTAE"), 0.0, 10.0, 0.1, 2.0))),
    iter_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_PROXI"), 0.2, 10.0, 0.1, 2.0))),
    balan_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BALAN"), 0.2, 2.5, 0.1, 1.0, Gtk::manage(new RTImage("rawtherapee-logo-16.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    transitweak_(Gtk::manage(new Adjuster(M("TP_LOCALLAB_TRANSITWEAK"), 0.5, 8.0, 0.1, 1.0))),

    avoid_(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_AVOID")))),

    lastObject_(-1),
    lastCoord_(new Coord()),
    nbSpotChanged_(false),
    selSpotChanged_(false),
    nameChanged_(false),
    visibilityChanged_(false),
    eventType(0),
    excluFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_EXCLUF"))))
{
    bool showtooltip = options.showtooltip;

    Gtk::HBox* const hbox1_ = Gtk::manage(new Gtk::HBox(true, 4));
    hbox1_->pack_start(button_add_);
    hbox1_->pack_start(button_delete_);
    hbox1_->pack_start(button_duplicate_);
    pack_start(*hbox1_);

    Gtk::HBox* const hbox2_ = Gtk::manage(new Gtk::HBox(true, 4));
    hbox2_->pack_start(button_rename_);
    hbox2_->pack_start(button_visibility_);
    pack_start(*hbox2_);

    buttonaddconn_ = button_add_.signal_clicked().connect(
                         sigc::mem_fun(*this, &ControlSpotPanel::on_button_add));
    buttondeleteconn_ = button_delete_.signal_clicked().connect(
                            sigc::mem_fun(*this, &ControlSpotPanel::on_button_delete));
    buttonduplicateconn_ = button_duplicate_.signal_clicked().connect(
                            sigc::mem_fun(*this, &ControlSpotPanel::on_button_duplicate));


    buttonrenameconn_ = button_rename_.signal_clicked().connect(
                            sigc::mem_fun(*this, &ControlSpotPanel::on_button_rename));
    buttonvisibilityconn_ = button_visibility_.signal_clicked().connect(
                                sigc::mem_fun(*this, &ControlSpotPanel::on_button_visibility));

    treeview_.set_grid_lines(Gtk::TREE_VIEW_GRID_LINES_VERTICAL);

    scrolledwindow_.add(treeview_);
    scrolledwindow_.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    scrolledwindow_.set_min_content_height(150);
    pack_start(scrolledwindow_);

    treemodel_ = Gtk::ListStore::create(spots_);

    treeview_.set_model(treemodel_);
    treeviewconn_ = treeview_.get_selection()->signal_changed().connect(
                        sigc::mem_fun(
                            *this, &ControlSpotPanel::controlspotChanged));

    auto cell = Gtk::manage(new Gtk::CellRendererText());
    int cols_count = treeview_.append_column("ID", *cell);
    auto col = treeview_.get_column(cols_count - 1);

    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &ControlSpotPanel::render_id));
    }

    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(M("TP_LOCALLAB_COL_NAME"), *cell);
    col = treeview_.get_column(cols_count - 1);

    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &ControlSpotPanel::render_name));
    }

    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(M("TP_LOCALLAB_COL_VIS"), *cell);
    col = treeview_.get_column(cols_count - 1);

    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &ControlSpotPanel::render_isvisible));
    }

    Gtk::HBox* const ctboxshape = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const labelshape = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_SHAPETYPE") + ":"));
    ctboxshape->pack_start(*labelshape, Gtk::PACK_SHRINK, 4);
    shape_->append(M("TP_LOCALLAB_ELI"));
    shape_->append(M("TP_LOCALLAB_RECT"));
    shape_->set_active(0);
    shapeconn_ = shape_->signal_changed().connect(
                     sigc::mem_fun(
                         *this, &ControlSpotPanel::shapeChanged));
    ctboxshape->pack_start(*shape_);
    pack_start(*ctboxshape);

    Gtk::HBox* const ctboxspotmethod = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const labelspotmethod = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_EXCLUTYPE") + ":"));
    ctboxspotmethod->pack_start(*labelspotmethod, Gtk::PACK_SHRINK, 4);
    if(showtooltip) ctboxspotmethod->set_tooltip_markup(M("TP_LOCALLAB_EXCLUTYPE_TOOLTIP"));
    spotMethod_->append(M("TP_LOCALLAB_EXNORM"));
    spotMethod_->append(M("TP_LOCALLAB_EXECLU"));
    spotMethod_->set_active(0);
    spotMethodconn_ = spotMethod_->signal_changed().connect(
                          sigc::mem_fun(
                              *this, &ControlSpotPanel::spotMethodChanged));
    ctboxspotmethod->pack_start(*spotMethod_);
    pack_start(*ctboxspotmethod);

    excluFrame->set_label_align(0.025, 0.5);
    if(showtooltip) excluFrame->set_tooltip_text(M("TP_LOCALLAB_EXCLUF_TOOLTIP"));
    ToolParamBlock* const excluBox = Gtk::manage(new ToolParamBlock());
    if(showtooltip) sensiexclu_->set_tooltip_text(M("TP_LOCALLAB_SENSIEXCLU_TOOLTIP"));
    sensiexclu_->setAdjusterListener(this);
    structexclu_->setAdjusterListener(this);
    excluBox->pack_start(*sensiexclu_);
    excluBox->pack_start(*structexclu_);
    excluFrame->add(*excluBox);
    pack_start(*excluFrame);

    Gtk::HBox* const ctboxshapemethod = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const labelshapemethod = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_STYPE") + ":"));
    ctboxshapemethod->pack_start(*labelshapemethod, Gtk::PACK_SHRINK, 4);
    if(showtooltip) ctboxshapemethod->set_tooltip_markup(M("TP_LOCALLAB_STYPE_TOOLTIP"));
    shapeMethod_->append(M("TP_LOCALLAB_IND"));
    shapeMethod_->append(M("TP_LOCALLAB_SYM"));
    shapeMethod_->append(M("TP_LOCALLAB_INDSL"));
    shapeMethod_->append(M("TP_LOCALLAB_SYMSL"));
    shapeMethod_->set_active(0);
    shapeMethodconn_ = shapeMethod_->signal_changed().connect(
                           sigc::mem_fun(
                               *this, &ControlSpotPanel::shapeMethodChanged));
    ctboxshapemethod->pack_start(*shapeMethod_);
    pack_start(*ctboxshapemethod);

    pack_start(*locX_);
    locX_->setAdjusterListener(this);

    pack_start(*locXL_);
    locXL_->setAdjusterListener(this);

    pack_start(*locY_);
    locY_->setAdjusterListener(this);

    pack_start(*locYT_);
    locYT_->setAdjusterListener(this);

    pack_start(*centerX_);
    centerX_->setAdjusterListener(this);

    pack_start(*centerY_);
    centerY_->setAdjusterListener(this);

    pack_start(*circrad_);
    circrad_->setAdjusterListener(this);

    Gtk::HBox* const ctboxqualitymethod = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const labelqualitymethod = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_QUAL_METHOD") + ":"));
    ctboxqualitymethod->pack_start(*labelqualitymethod, Gtk::PACK_SHRINK, 4);
    if(showtooltip) ctboxqualitymethod->set_tooltip_markup(M("TP_LOCALLAB_METHOD_TOOLTIP"));
//    qualityMethod_->append(M("TP_LOCALLAB_STD"));
    qualityMethod_->append(M("TP_LOCALLAB_ENH"));
    qualityMethod_->append(M("TP_LOCALLAB_ENHDEN"));
    qualityMethod_->set_active(1);
    qualityMethodconn_ = qualityMethod_->signal_changed().connect(
                             sigc::mem_fun(
                                 *this, &ControlSpotPanel::qualityMethodChanged));
    ctboxqualitymethod->pack_start(*qualityMethod_);
    pack_start(*ctboxqualitymethod);

    Gtk::Frame* const transitFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_TRANSIT")));
    transitFrame->set_label_align(0.025, 0.5);
    if(showtooltip) transitFrame->set_tooltip_text(M("TP_LOCALLAB_TRANSIT_TOOLTIP"));
    ToolParamBlock* const transitBox = Gtk::manage(new ToolParamBlock());
    if(showtooltip) transit_->set_tooltip_text(M("TP_LOCALLAB_TRANSIT_TOOLTIP"));
    if(showtooltip) transitweak_->set_tooltip_text(M("TP_LOCALLAB_TRANSITWEAK_TOOLTIP"));
    transit_->setAdjusterListener(this);
    transitweak_->setAdjusterListener(this);

    transitBox->pack_start(*transit_);
    transitBox->pack_start(*transitweak_);
    transitFrame->add(*transitBox);
    pack_start(*transitFrame);
    

    Gtk::Frame* const artifFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_ARTIF")));
    artifFrame->set_label_align(0.025, 0.5);
    if(showtooltip) artifFrame->set_tooltip_text(M("TP_LOCALLAB_ARTIF_TOOLTIP"));
    ToolParamBlock* const artifBox = Gtk::manage(new ToolParamBlock());
    thresh_->setAdjusterListener(this);
    struc_->setAdjusterListener(this);
    artifBox->pack_start(*struc_);
    artifBox->pack_start(*thresh_);
    artifBox->pack_start(*iter_);
    artifBox->pack_start(*balan_);
    iter_->setAdjusterListener(this);
    balan_->setAdjusterListener(this);
    artifFrame->add(*artifBox);
    pack_start(*artifFrame);

    avoidConn_  = avoid_->signal_toggled().connect(
            sigc::mem_fun(*this, &ControlSpotPanel::avoidChanged));
    pack_start(*avoid_);

    show_all();
}

void ControlSpotPanel::setEditProvider(EditDataProvider* provider)
{
    EditSubscriber::setEditProvider(provider);
}

void ControlSpotPanel::render_id(
    Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter)
{
    auto row = *iter;
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *>(cell);
    int value = row[spots_.id];
    ct->property_text() = std::to_string(value);
    // Render background color
    Gdk::RGBA color;
    if (row[spots_.mouseover]) { // Orange
        color.set_red(1.);
        color.set_green(100. / 255.);
        color.set_blue(0.);
        color.set_alpha(1.);
    } else { // Transparent black
        color.set_red(0.);
        color.set_green(0.);
        color.set_blue(0.);
        color.set_alpha(0.);
    }
    ct->property_background_rgba() = color;
}

void ControlSpotPanel::render_name(
    Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter)
{
    auto row = *iter;
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *>(cell);
    auto value = row[spots_.name];
    ct->property_text() = value;
    Gdk::RGBA color;
    // Render background color
    if (row[spots_.mouseover]) { // Orange
        color.set_red(1.);
        color.set_green(100. / 255.);
        color.set_blue(0.);
        color.set_alpha(1.);
    } else { // Transparent black
        color.set_red(0.);
        color.set_green(0.);
        color.set_blue(0.);
        color.set_alpha(0.);
    }
    ct->property_background_rgba() = color;
}

void ControlSpotPanel::render_isvisible(
    Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter)
{
    auto row = *iter;
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *>(cell);
    auto value = row[spots_.isvisible];

    if (value) {
        ct->property_text() = M("TP_LOCALLAB_ROW_VIS");
    } else {
        ct->property_text() = M("TP_LOCALLAB_ROW_NVIS");
    }

    // Render background color
    Gdk::RGBA color;
    if (row[spots_.mouseover]) { // Orange
        color.set_red(1.);
        color.set_green(100. / 255.);
        color.set_blue(0.);
        color.set_alpha(1.);
    } else { // Transparent black
        color.set_red(0.);
        color.set_green(0.);
        color.set_blue(0.);
        color.set_alpha(0.);
    }
    ct->property_background_rgba() = color;
}

void ControlSpotPanel::on_button_add()
{
    // printf("on_button_add\n");

    if (!listener) {
        return;
    }

    // Raise event
    nbSpotChanged_ = true;
    selSpotChanged_ = true;
    eventType = 1; // 1 = Spot creation event
    const int newId = getNewId();
    listener->panelChanged(EvLocallabSpotCreated, "ID#" + std::to_string(newId));
}

void ControlSpotPanel::on_button_delete()
{
    // printf("on_button_delete\n");

    if (!listener) {
        return;
    }

    // Raise event
    nbSpotChanged_ = true;
    selSpotChanged_ = true;
    eventType = 2; // 2 = Spot deletion event
    const int delId = getSelectedSpot();
    listener->panelChanged(EvLocallabSpotDeleted, "ID#" + std::to_string(delId));
}

void ControlSpotPanel::on_button_duplicate()
{
    // printf("on_button_duplicate\n");

    if (!listener) {
        return;
    }

    // Raise event
    const int selId = getSelectedSpot();
    if (selId == 0) { // No selected spot to duplicate
        return;
    }
    nbSpotChanged_ = true;
    selSpotChanged_ = true;
    eventType = 4; // 4 = Spot duplication event
    const int newId = getNewId();
    listener->panelChanged(EvLocallabSpotCreated, "ID#" + std::to_string(newId)
                                + " (" + M("TP_LOCALLAB_EV_DUPL") + " ID#"
                                + std::to_string(selId) + ")");
}

void ControlSpotPanel::on_button_rename()
{
    // printf("on_button_rename\n");

    if (!listener) {
        return;
    }

    // Get actual control spot name
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    const Gtk::TreeModel::Row row = *iter;
    const Glib::ustring actualname = row[spots_.name];

    // Launch windows to update spot name
    RenameDialog d(actualname,
                   static_cast<Gtk::Window &>(*get_toplevel()));
    int status = d.run();

    // Update actual name and raise event
    if (status == 1) {
        nameChanged_ = true;
        const auto newname = d.get_new_name();
        row[spots_.name] = newname;
        treeview_.columns_autosize();
        listener->panelChanged(EvLocallabSpotName, newname);
    }
}

void ControlSpotPanel::on_button_visibility()
{
    // printf("on_button_visibility\n");

    if (!listener) {
        return;
    }

    // Get selected control spot
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    const Gtk::TreeModel::Row row = *iter;

    // Update visibility
    row[spots_.isvisible] = !(bool)row[spots_.isvisible];
    updateControlSpotCurve(row);

    // Raise event
    visibilityChanged_ = true;
    const int id = getSelectedSpot();

    if ((bool)row[spots_.isvisible]) {
        listener->panelChanged(EvLocallabSpotVisibility, M("TP_LOCALLAB_EV_VIS") + " ID#" + std::to_string(id));
    } else {
        listener->panelChanged(EvLocallabSpotVisibility, M("TP_LOCALLAB_EV_NVIS") + " ID#" + std::to_string(id));
    }
}

void ControlSpotPanel::load_ControlSpot_param()
{
    // printf("load_ControlSpot_param\n");

    // Get selected control spot
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    const Gtk::TreeModel::Row row = *iter;

    // Load param in selected control spot
    shape_->set_active(row[spots_.shape]);
    spotMethod_->set_active(row[spots_.spotMethod]);
    sensiexclu_->setValue(static_cast<double>(row[spots_.sensiexclu]));
    structexclu_->setValue(static_cast<double>(row[spots_.structexclu]));
    struc_->setValue(static_cast<double>(row[spots_.struc]));
    shapeMethod_->set_active(row[spots_.shapeMethod]);
    locX_->setValue(static_cast<double>(row[spots_.locX]));
    locXL_->setValue(static_cast<double>(row[spots_.locXL]));
    locY_->setValue(static_cast<double>(row[spots_.locY]));
    locYT_->setValue(static_cast<double>(row[spots_.locYT]));
    centerX_->setValue(static_cast<double>(row[spots_.centerX]));
    centerY_->setValue(static_cast<double>(row[spots_.centerY]));
    circrad_->setValue(static_cast<double>(row[spots_.circrad]));
    qualityMethod_->set_active(row[spots_.qualityMethod]);
    transit_->setValue(static_cast<double>(row[spots_.transit]));
    thresh_->setValue(static_cast<double>(row[spots_.thresh]));
    iter_->setValue(static_cast<double>(row[spots_.iter]));
    balan_->setValue(static_cast<double>(row[spots_.balan]));
    transitweak_->setValue(static_cast<double>(row[spots_.transitweak]));
    avoid_->set_active(row[spots_.avoid]);
}

void ControlSpotPanel::controlspotChanged()
{
    // printf("controlspotChanged\n");

    if (!listener) {
        return;
    }

    // Raise event
    selSpotChanged_ = true;
    eventType = 3; // 3 = Spot selection event
    const int selId = getSelectedSpot();
    listener->panelChanged(EvLocallabSpotSelected, "ID#" + std::to_string(selId));
}

void ControlSpotPanel::shapeChanged()
{
    // printf("shapeChanged\n");

    // Get selected control spot
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;

    row[spots_.shape] = shape_->get_active_row_number();
    updateControlSpotCurve(row);

    // Raise event
    if (listener) {
        listener->panelChanged(EvLocallabSpotShape, shape_->get_active_text());
    }
}

void ControlSpotPanel::spotMethodChanged()
{
    // printf("spotMethodChanged\n");

    // Get selected control spot
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;

    row[spots_.spotMethod] = spotMethod_->get_active_row_number();

    // Update Control Spot GUI according to spotMethod_ combobox state (to be compliant with updateParamVisibility function)
    if (multiImage && spotMethod_->get_active_text() == M("GENERAL_UNCHANGED")) {
        excluFrame->show();
    } else if (spotMethod_->get_active_row_number() == 0) { // Normal case
        excluFrame->hide();
    } else { // Excluding case
        excluFrame->show();
    }

    // Raise event
    if (listener) {
        listener->panelChanged(EvLocallabSpotSpotMethod, spotMethod_->get_active_text());
    }
}

void ControlSpotPanel::shapeMethodChanged()
{
    // printf("shapeMethodChanged\n");

    const int method = shapeMethod_->get_active_row_number();

    // Get selected control spot
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;

    if (!batchMode && (method == 1 || method == 3)) { // Symmetrical cases
        disableParamlistener(true);
        locXL_->setValue(locX_->getValue());
        locYT_->setValue(locY_->getValue());
        disableParamlistener(false);

        row[spots_.shapeMethod] = shapeMethod_->get_active_row_number();
        row[spots_.locXL] = static_cast<int>(locX_->getValue());
        row[spots_.locYT] = static_cast<int>(locY_->getValue());

        updateControlSpotCurve(row);
    } else { // In batch mode, sliders are always independent
        row[spots_.shapeMethod] = shapeMethod_->get_active_row_number();
    }

    // Update Control Spot GUI according to shapeMethod_ combobox state (to be compliant with updateParamVisibility function)
    if (!batchMode) {
        if (method == 1 || method == 3) { // Symmetrical cases
            locXL_->hide();
            locYT_->hide();

            if (method == 1) { // 1 = Symmetrical (mouse)
                locX_->hide();
                locY_->hide();
                centerX_->hide();
                centerY_->hide();
            } else { // 3 = Symmetrical (mouse + sliders)
                locX_->show();
                locY_->show();
                centerX_->show();
                centerY_->show();
            }
        } else { // Independent cases
            if (method == 0) { // 0 = Independent (mouse)
                locX_->hide();
                locXL_->hide();
                locY_->hide();
                locYT_->hide();
                centerX_->hide();
                centerY_->hide();
            } else { // 2 = Independent (mouse + sliders)
                locX_->show();
                locXL_->show();
                locY_->show();
                locYT_->show();
                centerX_->show();
                centerY_->show();
            }
        }
    } else { // In batch mode, sliders are necessary shown
        locX_->show();
        locXL_->show();
        locY_->show();
        locYT_->show();
        centerX_->show();
        centerY_->show();
    }

    // Raise event
    if (listener) {
        listener->panelChanged(EvLocallabSpotShapeMethod, shapeMethod_->get_active_text());
    }
}

void ControlSpotPanel::qualityMethodChanged()
{
    // printf("qualityMethodChanged\n");

    // Get selected control spot
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;

    row[spots_.qualityMethod] = qualityMethod_->get_active_row_number();

    // Raise event
    if (listener) {
        listener->panelChanged(EvLocallabSpotQualityMethod, qualityMethod_->get_active_text());
    }
}

void ControlSpotPanel::updateParamVisibility()
{
    // printf("updateParamVisibility\n");

    // Update Control Spot GUI according to shapeMethod_ combobox state (to be compliant with shapeMethodChanged function)
    const int method = shapeMethod_->get_active_row_number();
    if (!batchMode) {
        if (method == 1 || method == 3) { // Symmetrical cases
            locXL_->hide();
            locYT_->hide();

            if (method == 1) { // 1 = Symmetrical (mouse)
                locX_->hide();
                locY_->hide();
                centerX_->hide();
                centerY_->hide();
            } else { // 3 = Symmetrical (mouse + sliders)
                locX_->show();
                locY_->show();
                centerX_->show();
                centerY_->show();
            }
        } else { // Independent cases
            if (method == 0) { // 0 = Independent (mouse)
                locX_->hide();
                locXL_->hide();
                locY_->hide();
                locYT_->hide();
                centerX_->hide();
                centerY_->hide();
            } else { // 2 = Independent (mouse + sliders)
                locX_->show();
                locXL_->show();
                locY_->show();
                locYT_->show();
                centerX_->show();
                centerY_->show();
            }
        }
    } else { // In batch mode, sliders are necessary shown
        locX_->show();
        locXL_->show();
        locY_->show();
        locYT_->show();
        centerX_->show();
        centerY_->show();
    }

    // Update Control Spot GUI according to spotMethod_ combobox state (to be compliant with spotMethodChanged function)
    if (multiImage && spotMethod_->get_active_text() == M("GENERAL_UNCHANGED")) {
        excluFrame->show();
    } else if (spotMethod_->get_active_row_number() == 0) { // Normal case
        excluFrame->hide();
    } else { // Excluding case
        excluFrame->show();
    }
}

void ControlSpotPanel::adjusterAutoToggled(Adjuster* a, bool newval)
{
}
void ControlSpotPanel::adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop)
{
}
void ControlSpotPanel::adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight)
{
}
void ControlSpotPanel::adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop)
{
}
void ControlSpotPanel::adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight)
{
}
void ControlSpotPanel::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
}

void ControlSpotPanel::adjusterChanged(Adjuster* a, double newval)
{
    // printf("adjusterChanged\n");

    const int method = shapeMethod_->get_active_row_number();

    // Get selected control spot
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;

    if (a == sensiexclu_) {
        row[spots_.sensiexclu] = (int) sensiexclu_->getValue();

        if (listener) {
            listener->panelChanged(EvLocallabSpotSensiexclu, sensiexclu_->getTextValue());
        }
    }

    if (a == structexclu_) {
        row[spots_.structexclu] = (int) structexclu_->getValue();

        if (listener) {
            listener->panelChanged(Evlocallabstructexlu, structexclu_->getTextValue());
        }
    }
    
    if (a == struc_) {
        row[spots_.struc] = struc_->getValue();

        if (listener) {
            listener->panelChanged(EvLocallabSpotStruc, struc_->getTextValue());
        }
    }

    if (a == locX_) {
        row[spots_.locX] = (int) locX_->getValue();

        if (!batchMode && (method == 1 || method == 3)) { // Symmetrical cases (in batch mode, sliders are always independent)
            disableParamlistener(true);
            locXL_->setValue(locX_->getValue());
            disableParamlistener(false);
            row[spots_.locXL] = (int) locXL_->getValue();
        }

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotLocX, locX_->getTextValue());
        }
    }

    if (a == locXL_) {
        row[spots_.locXL] = (int) locXL_->getValue();

        if (!batchMode && (method == 1 || method == 3)) { // Symmetrical cases (in batch mode, sliders are always independent)
            disableParamlistener(true);
            locX_->setValue(locXL_->getValue());
            disableParamlistener(false);
            row[spots_.locX] = (int) locX_->getValue();
        }

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotLocXL, locXL_->getTextValue());
        }
    }

    if (a == locY_) {
        row[spots_.locY] = (int) locY_->getValue();

        if (!batchMode && (method == 1 || method == 3)) { // Symmetrical cases (in batch mode, sliders are always independent)
            disableParamlistener(true);
            locYT_->setValue(locY_->getValue());
            disableParamlistener(false);
            row[spots_.locYT] = (int) locYT_->getValue();
        }

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotLocY, locY_->getTextValue());
        }
    }

    if (a == locYT_) {
        row[spots_.locYT] = (int) locYT_->getValue();

        if (!batchMode && (method == 1 || method == 3)) { // Symmetrical cases (in batch mode, sliders are always independent)
            disableParamlistener(true);
            locY_->setValue(locYT_->getValue());
            disableParamlistener(false);
            row[spots_.locY] = (int) locY_->getValue();
        }

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotLocYT, locYT_->getTextValue());
        }
    }

    if (a == centerX_ || a == centerY_) {
        row[spots_.centerX] = (int) centerX_->getValue();
        row[spots_.centerY] = (int) centerY_->getValue();

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotCenter, "X=" + centerX_->getTextValue() + ", Y=" + centerY_->getTextValue());
        }
    }

    if (a == circrad_) {
        row[spots_.circrad] = (int) circrad_->getValue();

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotCircrad, circrad_->getTextValue());
        }
    }

    if (a == transit_) {
        row[spots_.transit] = (int) transit_->getValue();

        if (listener) {
            listener->panelChanged(EvLocallabSpotTransit, transit_->getTextValue());
        }
    }


    if (a == thresh_) {
        row[spots_.thresh] = thresh_->getValue();

        if (listener) {
            listener->panelChanged(EvLocallabSpotThresh, thresh_->getTextValue());
        }
    }

    if (a == iter_) {
      //  row[spots_.iter] = (int) iter_->getValue();
        row[spots_.iter] = iter_->getValue();

        if (listener) {
            listener->panelChanged(EvLocallabSpotIter, iter_->getTextValue());
        }
    }

    if (a == balan_) {
        row[spots_.balan] = balan_->getValue();

        if (listener) {
            listener->panelChanged(EvLocallabSpotbalan, balan_->getTextValue());
        }
    }

    if (a == transitweak_) {
        row[spots_.transitweak] = transitweak_->getValue();

        if (listener) {
            listener->panelChanged(EvLocallabSpotTransitweak, transitweak_->getTextValue());
        }
    }

}

void ControlSpotPanel::avoidChanged()
{
    // printf("avoidChanged\n");

    // Get selected control spot
    const auto s = treeview_.get_selection();

    if (!s->count_selected_rows()) {
        return;
    }

    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;

    if (multiImage) {
        if (avoid_->get_inconsistent()) {
            avoid_->set_inconsistent(false);
            avoidConn_.block(true);
            avoid_->set_active(false);
            avoidConn_.block(false);
        }
    }

    row[spots_.avoid] = avoid_->get_active();

    // Raise event
    if (listener) {
        if (avoid_->get_active()) {
            listener->panelChanged(Evlocallabavoid, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabavoid, M("GENERAL_DISABLED"));
        }
    }
}

void ControlSpotPanel::disableParamlistener(bool cond)
{
    // printf("disableParamlistener: %d\n", cond);

    treeviewconn_.block(cond);
    buttonaddconn_.block(cond);
    buttondeleteconn_.block(cond);
    buttonduplicateconn_.block(cond);
    buttonrenameconn_.block(cond);
    buttonvisibilityconn_.block(cond);
    shapeconn_.block(cond);
    spotMethodconn_.block(cond);
    sensiexclu_->block(cond);
    structexclu_->block(cond);
    struc_->block(cond);
    shapeMethodconn_.block(cond);
    locX_->block(cond);
    locXL_->block(cond);
    locY_->block(cond);
    locYT_->block(cond);
    centerX_->block(cond);
    centerY_->block(cond);
    circrad_->block(cond);
    qualityMethodconn_.block(cond);
    transit_->block(cond);
    thresh_->block(cond);
    iter_->block(cond);
    balan_->block(cond);
    transitweak_->block(cond);
    avoidConn_.block(cond);
}

void ControlSpotPanel::setParamEditable(bool cond)
{
    // printf("setParamEditable: %d\n", cond);

    shape_->set_sensitive(cond);
    spotMethod_->set_sensitive(cond);
    sensiexclu_->set_sensitive(cond);
    structexclu_->set_sensitive(cond);
    struc_->set_sensitive(cond);
    shapeMethod_->set_sensitive(cond);
    locX_->set_sensitive(cond);
    locXL_->set_sensitive(cond);
    locY_->set_sensitive(cond);
    locYT_->set_sensitive(cond);
    centerX_->set_sensitive(cond);
    centerY_->set_sensitive(cond);
    circrad_->set_sensitive(cond);
    qualityMethod_->set_sensitive(cond);
    transit_->set_sensitive(cond);
    thresh_->set_sensitive(cond);
    iter_->set_sensitive(cond);
    balan_->set_sensitive(cond);
    transitweak_->set_sensitive(cond);
    avoid_->set_sensitive(cond);
}

void ControlSpotPanel::addControlSpotCurve(Gtk::TreeModel::Row row)
{
    // printf("addControlSpotCurve\n");

    if (row[spots_.curveid] > 0) { // Row has already an associated curve
        return;
    }

    // Creation of visibleGeometry
    Circle* cirX;
    cirX = new Circle();
    cirX->radius = 4.;
    cirX->filled = true;
    cirX->datum = Geometry::IMAGE;
    Circle* cirXL;
    cirXL = new Circle();
    cirXL->radius = 4.;
    cirXL->filled = true;
    cirXL->datum = Geometry::IMAGE;
    Circle* cirY;
    cirY = new Circle();
    cirY->radius = 4.;
    cirY->filled = true;
    cirY->datum = Geometry::IMAGE;
    Circle* cirYT;
    cirYT = new Circle();
    cirYT->radius = 4.;
    cirYT->filled = true;
    cirYT->datum = Geometry::IMAGE;
    Circle* centerCircle;
    centerCircle = new Circle();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = true;
    Ellipse* shape_ellipse;
    shape_ellipse = new Ellipse();
    shape_ellipse->datum = Geometry::IMAGE;
    shape_ellipse->radiusInImageSpace = true;
    Rectangle* shape_rectangle;
    shape_rectangle = new Rectangle();
    shape_rectangle->datum = Geometry::IMAGE;
    EditSubscriber::visibleGeometry.push_back(centerCircle); // (curveid - 1) * 7
    EditSubscriber::visibleGeometry.push_back(shape_ellipse); // (curveid - 1) * 7 + 1
    EditSubscriber::visibleGeometry.push_back(shape_rectangle); // (curveid - 1) * 7 + 2
    EditSubscriber::visibleGeometry.push_back(cirX); // (curveid - 1) * 7 + 3
    EditSubscriber::visibleGeometry.push_back(cirXL); // (curveid - 1) * 7 + 4
    EditSubscriber::visibleGeometry.push_back(cirY); // (curveid - 1) * 7 + 5
    EditSubscriber::visibleGeometry.push_back(cirYT); // (curveid - 1) * 7 + 6

    // Creation of mouseOverGeometry
    cirX = new Circle();
    cirX->radius = 4.;
    cirX->filled = true;
    cirX->datum = Geometry::IMAGE;
    cirXL = new Circle();
    cirXL->radius = 4.;
    cirXL->filled = true;
    cirXL->datum = Geometry::IMAGE;
    cirY = new Circle();
    cirY->radius = 4.;
    cirY->filled = true;
    cirY->datum = Geometry::IMAGE;
    cirYT = new Circle();
    cirYT->radius = 4.;
    cirYT->filled = true;
    cirYT->datum = Geometry::IMAGE;
    centerCircle = new Circle();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = true;
    shape_ellipse = new Ellipse();
    shape_ellipse->datum = Geometry::IMAGE;
    shape_ellipse->radiusInImageSpace = true;
    shape_rectangle = new Rectangle();
    shape_rectangle->datum = Geometry::IMAGE;
    EditSubscriber::mouseOverGeometry.push_back(centerCircle);  // (curveid - 1) * 7
    EditSubscriber::mouseOverGeometry.push_back(shape_ellipse);  // (curveid - 1) * 7 + 1
    EditSubscriber::mouseOverGeometry.push_back(shape_rectangle);  // (curveid - 1) * 7 + 2
    EditSubscriber::mouseOverGeometry.push_back(cirX);  // (curveid - 1) * 7 + 3
    EditSubscriber::mouseOverGeometry.push_back(cirXL);  // (curveid - 1) * 7 + 4
    EditSubscriber::mouseOverGeometry.push_back(cirY);  // (curveid - 1) * 7 + 5
    EditSubscriber::mouseOverGeometry.push_back(cirYT);  // (curveid - 1) * 7 + 6

    row[spots_.curveid] = EditSubscriber::visibleGeometry.size() / 7;
}

void ControlSpotPanel::updateControlSpotCurve(Gtk::TreeModel::Row row)
{
    const int curveid_ = row[spots_.curveid];
    EditDataProvider* const dataProvider = getEditProvider();

    // printf("updateControlSpotCurve: %d\n", curveid_);

    if (curveid_ == 0 || !dataProvider) { // Row has no associated curve or there is no EditProvider
        return;
    }

    int imW = 0;
    int imH = 0;
    dataProvider->getImageSize(imW, imH);

    if (!imW || !imH) { // No image loaded
        return;
    }

    const int centerX_ = row[spots_.centerX];
    const int centerY_ = row[spots_.centerY];
    const int circrad_ = row[spots_.circrad];
    const int locX_ = row[spots_.locX];
    const int locXL_ = row[spots_.locXL];
    const int locY_ = row[spots_.locY];
    const int locYT_ = row[spots_.locYT];
    const int shape_ = row[spots_.shape];
    const bool isvisible_ = row[spots_.isvisible];

    const int decayX = (double)locX_ * (double)imW / 2000.;
    const int decayXL = (double)locXL_ * (double)imW / 2000.;
    const int decayY = (double)locY_ * (double)imH / 2000.;
    const int decayYT = (double)locYT_ * (double)imH / 2000.;
    const rtengine::Coord origin((double)imW / 2. + (double)centerX_ * (double)imW / 2000., (double)imH / 2. + (double)centerY_ * (double)imH / 2000.);

    const auto updateSelectionCircle = [&](Geometry * geometry, const int offsetX, const int offsetY) {
        const auto cir = static_cast<Circle*>(geometry);
        cir->center.x = origin.x + offsetX;
        cir->center.y = origin.y + offsetY;
    };

    const auto updateCenterCircle = [&](Geometry * geometry) {
        const auto circle = static_cast<Circle*>(geometry);
        circle->center = origin;
        circle->radius = circrad_;
    };

    const auto updateEllipse = [&](Geometry * geometry) {
        const auto ellipse = static_cast<Ellipse*>(geometry);
        ellipse->center = origin;
        ellipse->radX = decayX;
        ellipse->radXL = decayXL;
        ellipse->radY = decayY;
        ellipse->radYT = decayYT;
    };

    const auto updateRectangle = [&](Geometry * geometry) {
        const auto rectangle = static_cast<Rectangle*>(geometry);
        rectangle->bottomRight.x = origin.x + decayX;
        rectangle->bottomRight.y = origin.y + decayY;
        rectangle->topLeft.x = origin.x - decayXL;
        rectangle->topLeft.y = origin.y - decayYT;
    };

    updateCenterCircle(visibleGeometry.at((curveid_ - 1) * 7));
    updateCenterCircle(mouseOverGeometry.at((curveid_ - 1) * 7));

    updateEllipse(visibleGeometry.at((curveid_ - 1) * 7 + 1));
    updateEllipse(mouseOverGeometry.at((curveid_ - 1) * 7 + 1));

    updateRectangle(visibleGeometry.at((curveid_ - 1) * 7 + 2));
    updateRectangle(mouseOverGeometry.at((curveid_ - 1) * 7 + 2));

    updateSelectionCircle(visibleGeometry.at((curveid_ - 1) * 7 + 3), decayX, 0.);
    updateSelectionCircle(mouseOverGeometry.at((curveid_ - 1) * 7 + 3), decayX, 0.);

    updateSelectionCircle(visibleGeometry.at((curveid_ - 1) * 7 + 4), -decayXL, 0.);
    updateSelectionCircle(mouseOverGeometry.at((curveid_ - 1) * 7 + 4), -decayXL, 0.);

    updateSelectionCircle(visibleGeometry.at((curveid_ - 1) * 7 + 5), 0., decayY);
    updateSelectionCircle(mouseOverGeometry.at((curveid_ - 1) * 7 + 5), 0., decayY);

    updateSelectionCircle(visibleGeometry.at((curveid_ - 1) * 7 + 6), 0., -decayYT);
    updateSelectionCircle(mouseOverGeometry.at((curveid_ - 1) * 7 + 6), 0., -decayYT);

    // Update Arcellipse/Rectangle visibility according to shape and visibility
    if (isvisible_) {
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7)->setActive(true); // centerCircle
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 3)->setActive(true); // cirX
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 4)->setActive(true); // cirXL
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 5)->setActive(true); // cirY
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 6)->setActive(true); // cirYT

        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7)->setActive(true); // centerCircle
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 3)->setActive(true); // cirX
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 4)->setActive(true); // cirXL
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 5)->setActive(true); // cirY
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 6)->setActive(true); // cirYT

        if (shape_ == 0) { // 0 = Ellipse
            EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 1)->setActive(true); // shape_ellipse
            EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 2)->setActive(false); // shape_rectangle

            EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 1)->setActive(true); // shape_ellipse
            EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 2)->setActive(false); // shape_rectangle
        } else { // 1 = Rectangle
            EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 1)->setActive(false); // shape_ellipse
            EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 2)->setActive(true); // shape_rectangle

            EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 1)->setActive(false); // shape_ellipse
            EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 2)->setActive(true); // shape_rectangle
        }
    } else {
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7)->setActive(false); // centerCircle
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 1)->setActive(false); // shape_ellipse
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 2)->setActive(false); // shape_rectangle
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 3)->setActive(false); // cirX
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 4)->setActive(false); // cirXL
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 5)->setActive(false); // cirY
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 7 + 6)->setActive(false); // cirYT

        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7)->setActive(false); // centerCircle
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 1)->setActive(false); // shape_ellipse
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 2)->setActive(false); // shape_rectangle
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 3)->setActive(false); // cirX
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 4)->setActive(false); // cirXL
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 5)->setActive(false); // cirY
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 7 + 6)->setActive(false); // cirYT
    }
}

void ControlSpotPanel::deleteControlSpotCurve(Gtk::TreeModel::Row row)
{
    const int curveid_ = row[spots_.curveid];

    // printf("deleteControlSpotCurve: %d\n", curveid_);

    if (curveid_ == 0) { // Row has no associated curve
        return;
    }

    // visibleGeometry
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin() + (curveid_ - 1) * 7 + 6);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin() + (curveid_ - 1) * 7 + 5);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin() + (curveid_ - 1) * 7 + 4);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin() + (curveid_ - 1) * 7 + 3);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin() + (curveid_ - 1) * 7 + 2);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin() + (curveid_ - 1) * 7 + 1);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin() + (curveid_ - 1) * 7);

    // mouseOverGeometry
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin() + (curveid_ - 1) * 7 + 6);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin() + (curveid_ - 1) * 7 + 5);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin() + (curveid_ - 1) * 7 + 4);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin() + (curveid_ - 1) * 7 + 3);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin() + (curveid_ - 1) * 7 + 2);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin() + (curveid_ - 1) * 7 + 1);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin() + (curveid_ - 1) * 7);

    row[spots_.curveid] = 0; // Reset associated curve id

    // Reordering curve id
    Gtk::TreeModel::Children children = treemodel_->children();

    for (auto iter = children.begin(); iter != children.end(); iter++) {
        Gtk::TreeModel::Row r = *iter;

        if (r[spots_.curveid] > curveid_) {
            r[spots_.curveid] = r[spots_.curveid] - 1;
        }
    }
}

void ControlSpotPanel::updateCurveOpacity(Gtk::TreeModel::Row selectedRow)
{
    const int curveid_ = selectedRow[spots_.curveid];

    // printf("updateCurveOpacity: %d\n", curveid_);

    if (curveid_ == 0) { // Row has no associated curve
        return;
    }

    for (int it_ = 0; it_ < (int) EditSubscriber::visibleGeometry.size(); it_++) {
        if ((it_ < ((curveid_ - 1) * 7)) || (it_ > ((curveid_ - 1) * 7) + 6)) { // it_ does not belong to selected curve
            EditSubscriber::visibleGeometry.at(it_)->opacity = 25.;
        } else {
            EditSubscriber::visibleGeometry.at(it_)->opacity = 75.;
        }
    }
}

CursorShape ControlSpotPanel::getCursor(int objectID) const
{
    // printf("Object ID: %d\n", objectID);

    // When there is no control spot (i.e. no selected row), objectID can unexpectedly be different from -1 and produced not desired behavior
    const auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
         return CSHandOpen;
     }

    int rem_ = objectID % 7;

    switch (rem_) {
        case (0): // centerCircle: (curveid_ - 1) * 7
            return CSMove2D;

        case (1): // shape_ellipse: (curveid_ - 1) * 7 + 1
            return CSMove2D;

        case (2): // shape_rectangle: (curveid_ - 1) * 7 + 2
            return CSMove2D;

        case (3): // cirX: (curveid_ - 1) * 7 + 3
            return CSMove1DH;

        case (4): // cirXL: (curveid_ - 1) * 7 + 4
            return CSMove1DH;

        case (5): // cirY: (curveid_ - 1) * 7 + 5
            return CSMove1DV;

        case (6): // cirYT: (curveid_ - 1) * 7 + 6
            return CSMove1DV;

        default:
            return CSHandOpen;
    }
}

bool ControlSpotPanel::mouseOver(int modifierKey)
{
    EditDataProvider* editProvider_ = getEditProvider();
    const auto s = treeview_.get_selection();

    if (!editProvider_ || !s->count_selected_rows()) { // When there is no control spot (i.e. no selected row), objectID can unexpectedly be different from -1 and produced not desired behavior
        return false;
    }

    // Get selected row
    const auto selIter = s->get_selected();
    Gtk::TreeModel::Row selRow = *selIter;

    int object_ = editProvider_->object;

    if (object_ != lastObject_) {
        if (object_ == -1) {
            // Reset mouseOver preview for visibleGeometry
            for (int it_ = 0; it_ < (int) EditSubscriber::visibleGeometry.size(); it_++) {
                EditSubscriber::visibleGeometry.at(it_)->state = Geometry::NORMAL;
            }

            // Reset mouseOver preview for TreeView
            Gtk::TreeModel::Children children = treemodel_->children();
            Gtk::TreeModel::Children::iterator iter;

            for (iter = children.begin(); iter != children.end(); iter++) {
                Gtk::TreeModel::Row row = *iter;
                row[spots_.mouseover] = false;
            }

            // Actualize lastObject_
            lastObject_ = object_;
            return false;
        }

        int curveId_ = object_ / 7 + 1;
        int rem = object_ % 7;

        // Manage mouseOver preview for TreeView
        Gtk::TreeModel::Children children = treemodel_->children();
        Gtk::TreeModel::Children::iterator iter;
        for (iter = children.begin(); iter != children.end(); iter++) {
            Gtk::TreeModel::Row row = *iter;
            if (row[spots_.curveid] == curveId_ && *row != *selRow) {
                row[spots_.mouseover] = true;
            } else {
                row[spots_.mouseover] = false;
            }
        }

        for (int it_ = 0; it_ < (int) EditSubscriber::visibleGeometry.size(); it_++) {
            if ((it_ < ((curveId_ - 1) * 7)) || (it_ > ((curveId_ - 1) * 7) + 6)) { // it_ does not belong to cursor pointed curve
                EditSubscriber::visibleGeometry.at(it_)->state = Geometry::NORMAL;
            }
        }

        const int method = shapeMethod_->get_active_row_number();

        // Circle, Arcellipses and Rectangle
        if (rem >= 0 && rem < 3) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 1)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 2)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 3)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 4)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 5)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 6)->state = Geometry::PRELIGHT;
        } else {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 2)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 3)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 4)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 4)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 5)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 6)->state = Geometry::NORMAL;
        }

        // cirX
        if (rem == 3) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 3)->state = Geometry::PRELIGHT;

            if (method == 1 || method == 3) { // Symmetrical cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 4)->state = Geometry::PRELIGHT;
            }
        }

        // cirXL
        if (rem == 4) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 4)->state = Geometry::PRELIGHT;

            if (method == 1 || method == 3) { // Symmetrical cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 3)->state = Geometry::PRELIGHT;
            }
        }

        // cirY
        if (rem == 5) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 5)->state = Geometry::PRELIGHT;

            if (method == 1 || method == 3) { // Symmetrical cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 6)->state = Geometry::PRELIGHT;
            }
        }

        // cirYT
        if (rem == 6) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 6)->state = Geometry::PRELIGHT;

            if (method == 1 || method == 3) { // Symmetrical cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 7 + 5)->state = Geometry::PRELIGHT;
            }
        }

        lastObject_ = object_;
        return true;
    }

    return false;
}

bool ControlSpotPanel::button1Pressed(int modifierKey)
{
    // printf("button1Pressed\n");

    EditDataProvider *provider = getEditProvider();
    const auto s = treeview_.get_selection();

    if (!provider || lastObject_ == -1 || !s->count_selected_rows()) { // When there is no control spot (i.e. no selected row), objectID can unexpectedly be different from -1 and produced not desired behavior
        return false;
    }

    // Select associated control spot
    int curveId_ = lastObject_ / 7 + 1;
    Gtk::TreeModel::Children children = treemodel_->children();

    for (auto iter = children.begin(); iter != children.end(); iter++) {
        Gtk::TreeModel::Row r = *iter;

        if (r[spots_.curveid] == curveId_) {
            treeview_.set_cursor(treemodel_->get_path(r));
            break;
        }
    }

    lastCoord_->set(provider->posImage.x + provider->deltaImage.x, provider->posImage.y + provider->deltaImage.y);
    EditSubscriber::action = EditSubscriber::Action::DRAGGING;    
    return true;
}

bool ControlSpotPanel::button1Released()
{
    // printf("button1Released\n");
    EditSubscriber::action = EditSubscriber::Action::NONE;
    return true;
}

bool ControlSpotPanel::drag1(int modifierKey)
{
    // printf("drag1\n");

    EditDataProvider *provider = getEditProvider();
    const auto s = treeview_.get_selection();

    if (!provider || lastObject_ == -1 || !s->count_selected_rows()) { // When there is no control spot (i.e. no selected row), objectID can unexpectedly be different from -1 and produced not desired behavior
        return false;
    }

    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;

    int imW, imH;
    provider->getImageSize(imW, imH);
    int rem = lastObject_ % 7;
    int method = shapeMethod_->get_active_row_number();
    Coord* newCoord = new Coord(provider->posImage.x + provider->deltaImage.x, provider->posImage.y + provider->deltaImage.y);

    // Circle, Ellipses and Rectangle
    if (rem >= 0 && rem < 3) {
        double deltaX = (double (newCoord->x) - double (lastCoord_->x)) * 2000. / double (imW);
        double deltaY = (double (newCoord->y) - double (lastCoord_->y)) * 2000. / double (imH);
        centerX_->setValue(centerX_->getValue() + deltaX);
        centerY_->setValue(centerY_->getValue() + deltaY);
        row[spots_.centerX] = (int) centerX_->getValue();
        row[spots_.centerY] = (int) centerY_->getValue();

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotCenter, "X=" + centerX_->getTextValue() + ", Y=" + centerY_->getTextValue());
        }
    }

    // cirX
    if (rem == 3) {
        double deltaX = (double (newCoord->x) - double (lastCoord_->x)) * 2000. / double (imW);
        locX_->setValue(locX_->getValue() + deltaX);
        row[spots_.locX] = (int) locX_->getValue();

        if (method == 1 || method == 3) { // Symmetrical cases
            disableParamlistener(true);
            locXL_->setValue(locX_->getValue());
            disableParamlistener(false);
            row[spots_.locXL] = (int) locXL_->getValue();
        }

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotLocX, locX_->getTextValue());
        }
    }

    // cirXL
    if (rem == 4) {
        double deltaXL = (double (lastCoord_->x) - double (newCoord->x)) * 2000. / double (imW);
        locXL_->setValue(locXL_->getValue() + deltaXL);
        row[spots_.locXL] = (int) locXL_->getValue();

        if (method == 1 || method == 3) { // Symmetrical cases
            disableParamlistener(true);
            locX_->setValue(locXL_->getValue());
            disableParamlistener(false);
            row[spots_.locX] = (int) locX_->getValue();
        }

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotLocXL, locXL_->getTextValue());
        }
    }

    // cirY
    if (rem == 5) {
        double deltaY = (double (newCoord->y) - double (lastCoord_->y)) * 2000. / double (imH);
        locY_->setValue(locY_->getValue() + deltaY);
        row[spots_.locY] = (int) locY_->getValue();

        if (method == 1 || method == 3) { // Symmetrical cases
            disableParamlistener(true);
            locYT_->setValue(locY_->getValue());
            disableParamlistener(false);
            row[spots_.locYT] = (int) locYT_->getValue();
        }

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotLocY, locY_->getTextValue());
        }
    }

    // cirYT
    if (rem == 6) {
        double deltaYT = (double (lastCoord_->y) - double (newCoord->y)) * 2000. / double (imH);
        locYT_->setValue(locYT_->getValue() + deltaYT);
        row[spots_.locYT] = (int) locYT_->getValue();

        if (method == 1 || method == 3) { // Symmetrical cases
            disableParamlistener(true);
            locY_->setValue(locYT_->getValue());
            disableParamlistener(false);
            row[spots_.locY] = (int) locY_->getValue();
        }

        updateControlSpotCurve(row);

        if (listener) {
            listener->panelChanged(EvLocallabSpotLocYT, locYT_->getTextValue());
        }
    }

    lastCoord_->set(newCoord->x, newCoord->y);
    return true;
}

int ControlSpotPanel::getEventType()
{
    const int tmp = eventType;
    eventType = 0; // Re-initialization at 0 if event type gotten
    return tmp;
}

ControlSpotPanel::SpotRow* ControlSpotPanel::getSpot(int id)
{
    // printf("getSpot: %d\n", id);

    MyMutex::MyLock lock(mTreeview);

    SpotRow* r = new SpotRow();

    Gtk::TreeModel::Children children = treemodel_->children();
    Gtk::TreeModel::Children::iterator iter;

    for (iter = children.begin(); iter != children.end(); iter++) {
        Gtk::TreeModel::Row row = *iter;

        if (row[spots_.id] == id) {
            r->id = row[spots_.id];
            r->name = row[spots_.name];
            r->isvisible = row[spots_.isvisible];
            r->shape = row[spots_.shape];
            r->spotMethod = row[spots_.spotMethod];
            r->sensiexclu = row[spots_.sensiexclu];
            r->structexclu = row[spots_.structexclu];
            r->struc = row[spots_.struc];
            r->shapeMethod = row[spots_.shapeMethod];
            r->locX = row[spots_.locX];
            r->locXL = row[spots_.locXL];
            r->locY = row[spots_.locY];
            r->locYT = row[spots_.locYT];
            r->centerX = row[spots_.centerX];
            r->centerY = row[spots_.centerY];
            r->circrad = row[spots_.circrad];
            r->qualityMethod = row[spots_.qualityMethod];
            r->transit = row[spots_.transit];
            r->thresh = row[spots_.thresh];
            r->iter = row[spots_.iter];
            r->balan = row[spots_.balan];
            r->transitweak = row[spots_.transitweak];
            r->avoid = row[spots_.avoid];

            return r;
        }
    }

    return nullptr;
}

std::vector<int>* ControlSpotPanel::getSpotIdList()
{
    MyMutex::MyLock lock(mTreeview);

    std::vector<int>* r = new std::vector<int>();

    Gtk::TreeModel::Children children = treemodel_->children();
    Gtk::TreeModel::Children::iterator iter;

    for (iter = children.begin(); iter != children.end(); iter++) {
        Gtk::TreeModel::Row row = *iter;
        r->push_back(row[spots_.id]);
    }

    return r;
}

int ControlSpotPanel::getSelectedSpot()
{
    // printf("getSelectedSpot\n");

    MyMutex::MyLock lock(mTreeview);

    const auto s = treeview_.get_selection();

    // Check if treeview has row, otherwise return 0
    if (!s->count_selected_rows()) {
        return 0;
    }

    auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;
    int id = row[spots_.id];

    return id;
}

void ControlSpotPanel::setSelectedSpot(int id)
{
    // printf("setSelectedSpot: %d\n", id);

    MyMutex::MyLock lock(mTreeview);

    disableParamlistener(true);

    Gtk::TreeModel::Children children = treemodel_->children();
    Gtk::TreeModel::Children::iterator iter;

    for (iter = children.begin(); iter != children.end(); iter++) {
        Gtk::TreeModel::Row row = *iter;

        if (row[spots_.id] == id) {
            treeview_.set_cursor(treemodel_->get_path(row));
            load_ControlSpot_param();
            updateParamVisibility();
            updateCurveOpacity(row);
        }
    }

    disableParamlistener(false);
}

int ControlSpotPanel::getNewId()
{
    MyMutex::MyLock lock(mTreeview);

    // Looking for maximum used id
    int max_row_id = 0;
    Gtk::TreeModel::Children children = treemodel_->children();
    Gtk::TreeModel::Children::iterator iter;

    for (iter = children.begin(); iter != children.end(); iter++) {
        Gtk::TreeModel::Row row = *iter;
        int iter_id = row[spots_.id];
        max_row_id = std::max(max_row_id, iter_id);
    }

    max_row_id++;

    return max_row_id;
}


void ControlSpotPanel::addControlSpot(SpotRow* newSpot)
{
    // printf("addControlSpot: %d\n", newSpot->id);

    MyMutex::MyLock lock(mTreeview);

    disableParamlistener(true);
    Gtk::TreeModel::Row row = * (treemodel_->append());
    row[spots_.mouseover] = false;
    row[spots_.id] = newSpot->id;
    row[spots_.name] = newSpot->name;
    row[spots_.isvisible] = newSpot->isvisible;
    row[spots_.curveid] = 0; // No associated curve
    row[spots_.shape] = newSpot->shape;
    row[spots_.spotMethod] = newSpot->spotMethod;
    row[spots_.sensiexclu] = newSpot->sensiexclu;
    row[spots_.structexclu] = newSpot->structexclu;
    row[spots_.struc] = newSpot->struc;
    row[spots_.shapeMethod] = newSpot->shapeMethod;
    row[spots_.locX] = newSpot->locX;
    row[spots_.locXL] = newSpot->locXL;
    row[spots_.locY] = newSpot->locY;
    row[spots_.locYT] = newSpot->locYT;
    row[spots_.centerX] = newSpot->centerX;
    row[spots_.centerY] = newSpot->centerY;
    row[spots_.circrad] = newSpot->circrad;
    row[spots_.qualityMethod] = newSpot->qualityMethod;
    row[spots_.transit] = newSpot->transit;
    row[spots_.thresh] = newSpot->thresh;
    row[spots_.iter] = newSpot->iter;
    row[spots_.balan] = newSpot->balan;
    row[spots_.transitweak] = newSpot->transitweak;
    row[spots_.avoid] = newSpot->avoid;
    updateParamVisibility();
    disableParamlistener(false);

    // Add associated control spot curve
    addControlSpotCurve(row);
    updateControlSpotCurve(row);
}

int ControlSpotPanel::updateControlSpot(SpotRow* spot)
{
    // printf("updateControlSpot: %d\n", spot->id);

    MyMutex::MyLock lock(mTreeview);

    disableParamlistener(true);

    Gtk::TreeModel::Children children = treemodel_->children();
    Gtk::TreeModel::Children::iterator iter;

    for (iter = children.begin(); iter != children.end(); iter++) {
        Gtk::TreeModel::Row row = *iter;

        if (row[spots_.id] == spot->id) {
            row[spots_.name] = spot->name;
            row[spots_.isvisible] = spot->isvisible;
            row[spots_.shape] = spot->shape;
            row[spots_.spotMethod] = spot->spotMethod;
            row[spots_.sensiexclu] = spot->sensiexclu;
            row[spots_.structexclu] = spot->structexclu;
            row[spots_.struc] = spot->struc;
            row[spots_.shapeMethod] = spot->shapeMethod;
            row[spots_.locX] = spot->locX;
            row[spots_.locXL] = spot->locXL;
            row[spots_.locY] = spot->locY;
            row[spots_.locYT] = spot->locYT;
            row[spots_.centerX] = spot->centerX;
            row[spots_.centerY] = spot->centerY;
            row[spots_.circrad] = spot->circrad;
            row[spots_.qualityMethod] = spot->qualityMethod;
            row[spots_.transit] = spot->transit;
            row[spots_.thresh] = spot->thresh;
            row[spots_.iter] = spot->iter;
            row[spots_.balan] = spot->balan;
            row[spots_.transitweak] = spot->transitweak;
            row[spots_.avoid] = spot->avoid;

            updateControlSpotCurve(row);
            updateParamVisibility();
            disableParamlistener(false);

            return 1;
        }
    }

    disableParamlistener(false);
    return 0;
}

void ControlSpotPanel::deleteControlSpot(int id)
{
    // printf("deleteControlSpot: %d\n", id);

    MyMutex::MyLock lock(mTreeview);

    disableParamlistener(true);

    Gtk::TreeModel::Children children = treemodel_->children();
    Gtk::TreeModel::Children::iterator iter;

    for (iter = children.begin(); iter != children.end(); iter++) {
        Gtk::TreeModel::Row row = *iter;

        if (row[spots_.id] == id) {
            deleteControlSpotCurve(row);
            treemodel_->erase(iter);
            break;
        }
    }

    disableParamlistener(false);
}

ControlSpotPanel::SpotEdited* ControlSpotPanel::getEditedStates()
{
    // printf("getEditedStates\n");

    SpotEdited* se = new SpotEdited();

    if (nbSpotChanged_) {
        se->nbspot = true;
        // nbSpotChanged_ = false;
    } else {
        se->nbspot = false;
    }

    if (selSpotChanged_) {
        se->selspot = true;
        // selSpotChanged_ = false;
    } else {
        se->selspot = false;
    }

    if (nameChanged_) {
        se->name = true;
        // nameChanged_ = false;
    } else {
        se->name = false;
    }

    if (visibilityChanged_) {
        se->isvisible = true;
        // visibilityChanged_ = false;
    } else {
        se->isvisible = false;
    }

    se->shape = shape_->get_active_text() != M("GENERAL_UNCHANGED");
    se->spotMethod = spotMethod_->get_active_text() != M("GENERAL_UNCHANGED");
    se->sensiexclu = sensiexclu_->getEditedState();
    se->structexclu = structexclu_->getEditedState();
    se->struc = struc_->getEditedState();
    se->shapeMethod = shapeMethod_->get_active_text() != M("GENERAL_UNCHANGED");
    se->locX = locX_->getEditedState();
    se->locXL = locXL_->getEditedState();
    se->locY = locY_->getEditedState();
    se->locYT = locYT_->getEditedState();
    se->centerX = centerX_->getEditedState();
    se-> centerY = centerY_->getEditedState();
    se->circrad = circrad_->getEditedState();
    se->qualityMethod = qualityMethod_->get_active_text() != M("GENERAL_UNCHANGED");
    se->transit = transit_->getEditedState();
    se->thresh = thresh_->getEditedState();
    se->iter = iter_->getEditedState();
    se->balan = balan_->getEditedState();
    se->transitweak = transitweak_->getEditedState();
    se->avoid = !avoid_->get_inconsistent();

    return se;
}

void ControlSpotPanel::setEditedStates(SpotEdited* se)
{
    // printf("setEditedStates\n");

    // Reset treeview edited states
    nbSpotChanged_ = false;
    selSpotChanged_ = false;
    nameChanged_ = false;
    visibilityChanged_ = false;

    // Disable params listeners
    disableParamlistener(true);

    // Set widgets edited states
    if (!se->nbspot || !se->selspot) {
        treeview_.set_sensitive(false);
        button_add_.set_sensitive(false);
        button_delete_.set_sensitive(false);
        button_duplicate_.set_sensitive(false);
        button_rename_.set_sensitive(false);
        button_visibility_.set_sensitive(false);
    } else {
        treeview_.set_sensitive(true);
        button_add_.set_sensitive(true);
        button_delete_.set_sensitive(true);
        button_duplicate_.set_sensitive(true);
        button_rename_.set_sensitive(se->name);
        button_visibility_.set_sensitive(se->isvisible);
    }

    if (!se->shape) {
        shape_->set_active_text(M("GENERAL_UNCHANGED"));
    }

    if (!se->spotMethod) {
        spotMethod_->set_active_text(M("GENERAL_UNCHANGED"));
    }

    sensiexclu_->setEditedState(se->sensiexclu ? Edited : UnEdited);
    structexclu_->setEditedState(se->structexclu ? Edited : UnEdited);
    struc_->setEditedState(se->struc ? Edited : UnEdited);

    if (!se->shapeMethod) {
        shapeMethod_->set_active_text(M("GENERAL_UNCHANGED"));
    }

    locX_->setEditedState(se->locX ? Edited : UnEdited);
    locXL_->setEditedState(se->locXL ? Edited : UnEdited);
    locY_->setEditedState(se->locY ? Edited : UnEdited);
    locYT_->setEditedState(se->locYT ? Edited : UnEdited);
    centerX_->setEditedState(se->centerX ? Edited : UnEdited);
    centerY_->setEditedState(se->centerY ? Edited : UnEdited);
    circrad_->setEditedState(se->circrad ? Edited : UnEdited);

    if (!se->qualityMethod) {
        qualityMethod_->set_active_text(M("GENERAL_UNCHANGED"));
    }

    transit_->setEditedState(se->transit ? Edited : UnEdited);
    thresh_->setEditedState(se->thresh ? Edited : UnEdited);
    iter_->setEditedState(se->iter ? Edited : UnEdited);
    balan_->setEditedState(se->balan ? Edited : UnEdited);
    transitweak_->setEditedState(se->transitweak ? Edited : UnEdited);
    avoid_->set_inconsistent(multiImage && !se->avoid);

    // Update Control Spot GUI according to widgets edited states
    updateParamVisibility();

    // Enable params listeners
    disableParamlistener(false);
}

void ControlSpotPanel::setDefaults(const rtengine::procparams::ProcParams * defParams, const ParamsEdited * pedited, int id)
{
    // Find vector index of given spot id (index = -1 if not found)
    int index = -1;

    for (int i = 0; i < (int)defParams->locallab.spots.size(); i++) {
        if (defParams->locallab.spots.at(i).id == id) {
            index = i;

            break;
        }
    }

    // Set default values for adjusters
    const rtengine::procparams::LocallabParams::LocallabSpot* defSpot = new rtengine::procparams::LocallabParams::LocallabSpot();

    if (index != -1 && index < (int)defParams->locallab.spots.size()) {
        defSpot = &defParams->locallab.spots.at(index);
    }

    sensiexclu_->setDefault((double)defSpot->sensiexclu);
    structexclu_->setDefault((double)defSpot->structexclu);
    struc_->setDefault(defSpot->struc);
    locX_->setDefault((double)defSpot->locX);
    locXL_->setDefault((double)defSpot->locXL);
    locY_->setDefault((double)defSpot->locY);
    locYT_->setDefault((double)defSpot->locYT);
    centerX_->setDefault((double)defSpot->centerX);
    centerY_->setDefault((double)defSpot->centerY);
    circrad_->setDefault((double)defSpot->circrad);
    transit_->setDefault((double)defSpot->transit);
    thresh_->setDefault(defSpot->thresh);
    iter_->setDefault(defSpot->iter);
    balan_->setDefault(defSpot->balan);
    transitweak_->setDefault(defSpot->transitweak);

    // Set default edited states for adjusters
    if (!pedited) {
        sensiexclu_->setDefaultEditedState(Irrelevant);
        structexclu_->setDefaultEditedState(Irrelevant);
        struc_->setDefaultEditedState(Irrelevant);
        locX_->setDefaultEditedState(Irrelevant);
        locXL_->setDefaultEditedState(Irrelevant);
        locY_->setDefaultEditedState(Irrelevant);
        locYT_->setDefaultEditedState(Irrelevant);
        centerX_->setDefaultEditedState(Irrelevant);
        centerY_->setDefaultEditedState(Irrelevant);
        circrad_->setDefaultEditedState(Irrelevant);
        transit_->setDefaultEditedState(Irrelevant);
        thresh_->setDefaultEditedState(Irrelevant);
        iter_->setDefaultEditedState(Irrelevant);
        balan_->setDefaultEditedState(Irrelevant);
        transitweak_->setDefaultEditedState(Irrelevant);
    } else {
        const LocallabParamsEdited::LocallabSpotEdited* defSpotState = new LocallabParamsEdited::LocallabSpotEdited(true);

        if (index != 1 && index < (int)pedited->locallab.spots.size()) {
            defSpotState = &pedited->locallab.spots.at(index);
        }

        sensiexclu_->setDefaultEditedState(defSpotState->sensiexclu ? Edited : UnEdited);
        structexclu_->setDefaultEditedState(defSpotState->structexclu ? Edited : UnEdited);
        struc_->setDefaultEditedState(defSpotState->struc ? Edited : UnEdited);
        locX_->setDefaultEditedState(defSpotState->locX ? Edited : UnEdited);
        locXL_->setDefaultEditedState(defSpotState->locXL ? Edited : UnEdited);
        locY_->setDefaultEditedState(defSpotState->locY ? Edited : UnEdited);
        locYT_->setDefaultEditedState(defSpotState->locYT ? Edited : UnEdited);
        centerX_->setDefaultEditedState(defSpotState->centerX ? Edited : UnEdited);
        centerY_->setDefaultEditedState(defSpotState->centerY ? Edited : UnEdited);
        circrad_->setDefaultEditedState(defSpotState->circrad ? Edited : UnEdited);
        transit_->setDefaultEditedState(defSpotState->transit ? Edited : UnEdited);
        thresh_->setDefaultEditedState(defSpotState->thresh ? Edited : UnEdited);
        iter_->setDefaultEditedState(defSpotState->iter ? Edited : UnEdited);
        balan_->setDefaultEditedState(defSpotState->balan ? Edited : UnEdited);
        transitweak_->setDefaultEditedState(defSpotState->transitweak ? Edited : UnEdited);
    }
}

void ControlSpotPanel::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    // Set batch mode for adjusters
    sensiexclu_->showEditedCB();
    structexclu_->showEditedCB();
    struc_->showEditedCB();
    locX_->showEditedCB();
    locXL_->showEditedCB();
    locY_->showEditedCB();
    locYT_->showEditedCB();
    centerX_->showEditedCB();
    centerY_->showEditedCB();
    circrad_->showEditedCB();
    transit_->showEditedCB();
    thresh_->showEditedCB();
    iter_->showEditedCB();
    balan_->showEditedCB();
    transitweak_->showEditedCB();

    // Set batch mode for comboBoxText
    shape_->append(M("GENERAL_UNCHANGED"));
    spotMethod_->append(M("GENERAL_UNCHANGED"));
    shapeMethod_->append(M("GENERAL_UNCHANGED"));
    qualityMethod_->append(M("GENERAL_UNCHANGED"));
}

//-----------------------------------------------------------------------------
// ControlSpots
//-----------------------------------------------------------------------------

ControlSpotPanel::ControlSpots::ControlSpots()
{
    add(mouseover);
    add(id);
    add(name);
    add(isvisible);
    add(curveid);
    add(shape);
    add(spotMethod);
    add(sensiexclu);
    add(structexclu);
    add(struc);
    add(shapeMethod);
    add(locX);
    add(locXL);
    add(locYT);
    add(locY);
    add(centerX);
    add(centerY);
    add(circrad);
    add(qualityMethod);
    add(transit);
    add(thresh);
    add(iter);
    add(balan);
    add(transitweak);
    add(avoid);
}

//-----------------------------------------------------------------------------
// RenameDialog
//-----------------------------------------------------------------------------

ControlSpotPanel::RenameDialog::RenameDialog(const Glib::ustring &actualname, Gtk::Window &parent):
    Gtk::Dialog(M("TP_LOCALLAB_REN_DIALOG_NAME"), parent)
{
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_REN_DIALOG_LAB"))), false, false, 4);

    newname_.set_text(actualname);
    hb->pack_start(newname_);

    get_content_area()->pack_start(*hb, Gtk::PACK_SHRINK, 4);

    add_button(M("GENERAL_OK"), 1);
    add_button(M("GENERAL_CANCEL"), 2);

    show_all_children();
}

Glib::ustring ControlSpotPanel::RenameDialog::get_new_name()
{
    return newname_.get_text();
}
