/*
 *  This file is part of RawTherapee.
 */

#include "../rtengine/rt_math.h"
#include "controlspotpanel.h"
#include "multilangmgr.h"
#include <iomanip>

using namespace rtengine;

//-----------------------------------------------------------------------------
// ControlSpotPanel
//-----------------------------------------------------------------------------

ControlSpotPanel::ControlSpotPanel():
    EditSubscriber(ET_OBJECTS),

    button_add_ ("Add"),
    button_delete_ ("Delete"),
    button_rename_ ("Rename"),

    shape_ (Gtk::manage (new MyComboBoxText ())),
    spotMethod_ (Gtk::manage (new MyComboBoxText ())),
    shapeMethod_ (Gtk::manage (new MyComboBoxText ())),
    qualityMethod_ (Gtk::manage (new MyComboBoxText ())),

    locX_ (Gtk::manage (new Adjuster (M ("TP_LOCAL_WIDTH"), 0, 2250, 1, 250))),
    locXL_ (Gtk::manage (new Adjuster (M ("TP_LOCAL_WIDTH_L"), 0, 2250, 1, 250))),
    locY_ (Gtk::manage (new Adjuster (M ("TP_LOCAL_HEIGHT"), 0, 2250, 1, 250))),
    locYT_ (Gtk::manage (new Adjuster (M ("TP_LOCAL_HEIGHT_T"), 0, 2250, 1, 250))),
    centerX_ (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CENTER_X"), -1000, 1000, 1, 0))),
    centerY_ (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CENTER_Y"), -1000, 1000, 1, 0))),
    circrad_ (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_CIRCRADIUS"), 2, 150, 1, 18))),
    transit_ (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_TRANSIT"), 5, 95, 1, 60))),
    thresh_ (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_THRES"), 1, 35, 1, 18))),
    iter_ (Gtk::manage (new Adjuster (M ("TP_LOCALLAB_PROXI"), 0, 60, 1, 0))),

    lastObject_ (-1),
    lastCoord_ (new Coord ())
{
    treeview_.set_grid_lines (Gtk::TREE_VIEW_GRID_LINES_VERTICAL);

    scrolledwindow_.add (treeview_);
    scrolledwindow_.set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    scrolledwindow_.set_min_content_height(150);

    pack_start (buttonbox_);
    pack_start (scrolledwindow_);

    buttonbox_.pack_start(button_add_, Gtk::PACK_SHRINK, 4);
    buttonbox_.pack_start(button_delete_, Gtk::PACK_SHRINK, 4);
    buttonbox_.pack_start(button_rename_);
    buttonbox_.set_layout(Gtk::BUTTONBOX_START);

    button_add_.signal_clicked().connect (
        sigc::mem_fun (*this, &ControlSpotPanel::on_button_add));
    button_delete_.signal_clicked().connect (
        sigc::mem_fun (*this, &ControlSpotPanel::on_button_delete));
    button_rename_.signal_clicked().connect (
        sigc::mem_fun (*this, &ControlSpotPanel::on_button_rename));

    treemodel_ = Gtk::ListStore::create (spots_);

    treeview_.set_model (treemodel_);
    treeview_.get_selection()->signal_changed().connect(
        sigc::mem_fun (
            *this, &ControlSpotPanel::controlspotChanged));

    auto cell = Gtk::manage (new Gtk::CellRendererText());
    int cols_count = treeview_.append_column ("ID", *cell);
    auto col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func (
            *cell, sigc::mem_fun (
                *this, &ControlSpotPanel::render_id));
    }

    cell = Gtk::manage (new Gtk::CellRendererText());
    cols_count = treeview_.append_column("Name", *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func (
            *cell, sigc::mem_fun (
                *this, &ControlSpotPanel::render_name));
    }

    cell = Gtk::manage (new Gtk::CellRendererText());
    cols_count = treeview_.append_column("Status", *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func (
            *cell, sigc::mem_fun (
                *this, &ControlSpotPanel::render_isvisible));
    }

    // TODO Reload saved control spots (don't forget autosize)

    // TODO Rectangle

    Gtk::HBox* const ctboxshape = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* const labelshape = Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_SHAPETYPE") + ":"));
    ctboxshape->pack_start (*labelshape, Gtk::PACK_SHRINK, 4);
    shape_->append (M ("TP_LOCALLAB_ELI"));
    shape_->append (M ("TP_LOCALLAB_RECT"));
    shape_->set_active (0);
    shapeconn_ = shape_->signal_changed ().connect (
                     sigc::mem_fun (
                         *this, &ControlSpotPanel::shapeChanged));
    ctboxshape->pack_start (*shape_);
    pack_start (*ctboxshape);

    Gtk::HBox* const ctboxspotmethod = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* const labelspotmethod = Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_EXCLUTYPE") + ":"));
    ctboxspotmethod->pack_start (*labelspotmethod, Gtk::PACK_SHRINK, 4);
    ctboxspotmethod->set_tooltip_markup (M ("TP_LOCALLAB_EXCLUTYPE_TOOLTIP"));
    spotMethod_->append (M ("TP_LOCALLAB_EXNORM"));
    spotMethod_->append (M ("TP_LOCALLAB_EXECLU"));
    spotMethod_->set_active (0);
    spotMethodconn_ = spotMethod_->signal_changed ().connect (
                          sigc::mem_fun (
                              *this, &ControlSpotPanel::save_ControlSpot_param));
    ctboxspotmethod->pack_start(*spotMethod_);
    pack_start(*ctboxspotmethod);

    Gtk::HBox* const ctboxshapemethod = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* const labelshapemethod = Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_STYPE") + ":"));
    ctboxshapemethod->pack_start (*labelshapemethod, Gtk::PACK_SHRINK, 4);
    ctboxshapemethod->set_tooltip_markup (M ("TP_LOCALLAB_STYPE_TOOLTIP"));
    shapeMethod_->append (M ("TP_LOCALLAB_IND"));
    shapeMethod_->append (M ("TP_LOCALLAB_SYM"));
    shapeMethod_->append (M ("TP_LOCALLAB_INDSL"));
    shapeMethod_->append (M ("TP_LOCALLAB_SYMSL"));
    shapeMethod_->set_active (0);
    shapeMethodconn_ = shapeMethod_->signal_changed ().connect (
                           sigc::mem_fun (
                               *this, &ControlSpotPanel::shapeMethodeChanged));
    ctboxshapemethod->pack_start (*shapeMethod_);
    pack_start (*ctboxshapemethod);

    pack_start (*locX_);
    locX_->setAdjusterListener (this);

    pack_start (*locXL_);
    locXL_->setAdjusterListener (this);

    pack_start (*locY_);
    locY_->setAdjusterListener (this);

    pack_start (*locYT_);
    locYT_->setAdjusterListener (this);

    pack_start (*centerX_);
    centerX_->setAdjusterListener (this);

    pack_start (*centerY_);
    centerY_->setAdjusterListener (this);

    pack_start (*circrad_);
    circrad_->setAdjusterListener (this);

    Gtk::HBox* const ctboxqualitymethod = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* const labelqualitymethod = Gtk::manage (new Gtk::Label (M ("TP_LOCALLAB_QUAL_METHOD") + ":"));
    ctboxqualitymethod->pack_start (*labelqualitymethod, Gtk::PACK_SHRINK, 4);
    ctboxqualitymethod->set_tooltip_markup(M("TP_LOCALLAB_METHOD_TOOLTIP"));
    qualityMethod_->append (M ("TP_LOCALLAB_STD"));
    qualityMethod_->append (M ("TP_LOCALLAB_ENH"));
    qualityMethod_->append (M ("TP_LOCALLAB_ENHDEN"));
    qualityMethod_->set_active(0);
    qualityMethodconn_ = qualityMethod_->signal_changed ().connect (
                             sigc::mem_fun (
                                 *this, &ControlSpotPanel::save_ControlSpot_param));
    ctboxqualitymethod->pack_start (*qualityMethod_);
    pack_start (*ctboxqualitymethod);

    pack_start (*transit_);
    transit_->set_tooltip_text (M ("TP_LOCALLAB_TRANSIT_TOOLTIP"));
    transit_->setAdjusterListener (this);

    Gtk::Frame* const artifFrame = Gtk::manage (new Gtk::Frame (M ("TP_LOCALLAB_ARTIF")));
    artifFrame->set_label_align (0.025, 0.5);
    artifFrame->set_tooltip_text (M ("TP_LOCALLAB_ARTIF_TOOLTIP"));
    ToolParamBlock* const artifBox = Gtk::manage (new ToolParamBlock ());
    artifBox->pack_start (*thresh_);
    thresh_->setAdjusterListener (this);
    artifBox->pack_start (*iter_);
    iter_->setAdjusterListener (this);
    artifFrame->add (*artifBox);
    pack_start (*artifFrame);

    // Set param widgets sensitive if there is at least one control spot
    auto s = treeview_.get_selection();
    if (!s->count_selected_rows ()) {
        setParamEditable (false);
    } else {
        setParamEditable (true);
    }

    show_all ();
}

void ControlSpotPanel::setEditProvider (EditDataProvider* provider)
{
    EditSubscriber::setEditProvider (provider);
}

namespace
{

template <class V>
Glib::ustring to_str (V n, int precision = 1)
{
    std::ostringstream buf;
    buf << std::setprecision (precision) << std::fixed << n;
    return buf.str ();
}

} // namespace

void ControlSpotPanel::render_id (
    Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter)
{
    auto row = *iter;
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *> (cell);
    int value = row[spots_.id];
    ct->property_text() = to_str (value);
}

void ControlSpotPanel::render_name (
    Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter)
{
    auto row = *iter;
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *> (cell);
    auto value = row[spots_.name];
    ct->property_text() = value;
}

void ControlSpotPanel::render_isvisible (
    Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter)
{
    auto row = *iter;
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *> (cell);
    auto value = row[spots_.isvisible];
    if (value) {
        ct->property_text () = "Visible";
    } else {
        ct->property_text () = "Not visible";
    }
}

void ControlSpotPanel::on_button_add ()
{
    printf("on_button_add\n");
    // Looking for maximum used id
    int max_row_id = 0;
    Gtk::TreeModel::Children children = treemodel_->children ();
    Gtk::TreeModel::Children::iterator iter;
    for (iter = children.begin (); iter != children.end (); iter++)
    {
        Gtk::TreeModel::Row row = *iter;
        int iter_id = row[spots_.id];
        max_row_id = std::max (max_row_id, iter_id);
    }

    // Adding row
    Gtk::TreeModel::Row row = * (treemodel_->append ());
    row[spots_.id] = max_row_id + 1;
    row[spots_.name] = "Control Spot #" + to_str(row[spots_.id]);
    row[spots_.isvisible] = true;
    row[spots_.curveid] = 0; // No associated curve
    row[spots_.shape] = 0;
    row[spots_.spotMethod] = 0;
    row[spots_.shapeMethod] = 2;
    row[spots_.locX] = 250;
    row[spots_.locXL] = 250;
    row[spots_.locY] = 250;
    row[spots_.locYT] = 250;
    row[spots_.centerX] = 0;
    row[spots_.centerY] = 0;
    row[spots_.circrad] = 18;
    row[spots_.qualityMethod] = 0;
    row[spots_.transit] = 60;
    row[spots_.thresh] = 18;
    row[spots_.iter] = 0;
    setParamEditable (true);

    // Select newly added row
    treeview_.set_cursor (treemodel_->get_path (row));

    // Add associated control spot curve
    addControlSpotCurve (row);
    updateControlSpotCurve (row);
    subscribe ();
}

void ControlSpotPanel::on_button_delete ()
{
    auto s = treeview_.get_selection ();
    if (!s->count_selected_rows ()) {
        return;
    }
    auto iter = s->get_selected ();
    Gtk::TreeModel::Row row = *iter;
    deleteControlSpotCurve (row);
    treemodel_->erase (iter);

    // Set param widgets unsensitive and unsubscribe if there is no more control spot
    s = treeview_.get_selection ();
    if (!s->count_selected_rows ()) {
        unsubscribe ();
        setParamEditable (false);
    }
}

void ControlSpotPanel::on_button_rename ()
{
    // Get actual control spot name
    const auto s = treeview_.get_selection ();
    if (!s->count_selected_rows ()) {
        return;
    }
    const auto iter = s->get_selected ();
    const Gtk::TreeModel::Row row = *iter;
    const Glib::ustring actualname = row[spots_.name];

    RenameDialog d (actualname,
                    static_cast<Gtk::Window &> (*get_toplevel ()));
    int status = d.run ();

    if (status == 1) {
        const auto newname = d.get_new_name ();
        row[spots_.name] = newname;
    }

    treeview_.columns_autosize ();
}

void ControlSpotPanel::save_ControlSpot_param ()
{
    printf("save_ControlSpot_param\n");
    // Get selected control spot
    const auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return;
    }
    const auto iter = s->get_selected();
    const Gtk::TreeModel::Row row = *iter;

    // Save param in selected control spot
    row[spots_.shape] = shape_->get_active_row_number ();
    row[spots_.spotMethod] = spotMethod_->get_active_row_number ();
    row[spots_.shapeMethod] = shapeMethod_->get_active_row_number ();
    row[spots_.locX] = static_cast<int> (locX_->getValue ());
    row[spots_.locXL] = static_cast<int> (locXL_->getValue ());
    row[spots_.locY] = static_cast<int> (locY_->getValue ());
    row[spots_.locYT] = static_cast<int> (locYT_->getValue ());
    row[spots_.centerX] = static_cast<int> (centerX_->getValue ());
    row[spots_.centerY] = static_cast<int> (centerY_->getValue ());
    row[spots_.circrad] = static_cast<int> (circrad_->getValue ());
    row[spots_.qualityMethod] = qualityMethod_->get_active_row_number ();
    row[spots_.transit] = static_cast<int> (transit_->getValue ());
    row[spots_.thresh] = static_cast<int> (thresh_->getValue ());
    row[spots_.iter] = static_cast<int> (iter_->getValue ());
}

void ControlSpotPanel::load_ControlSpot_param()
{
    printf("load_ControlSpot_param\n");
    // Get selected control spot
    const auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return;
    }
    const auto iter = s->get_selected();
    const Gtk::TreeModel::Row row = *iter;

    // Listener are deactivated to avoid unexpected even during param load
    disableParamlistener (true);

    // Load param in selected control spot
    shape_->set_active (row[spots_.shape]);
    spotMethod_->set_active (row[spots_.spotMethod]);
    shapeMethod_->set_active (row[spots_.shapeMethod]);
    locX_->setValue (static_cast<double> (row[spots_.locX]));
    locXL_->setValue (static_cast<double> (row[spots_.locXL]));
    locY_->setValue (static_cast<double> (row[spots_.locY]));
    locYT_->setValue (static_cast<double> (row[spots_.locYT]));
    centerX_->setValue (static_cast<double> (row[spots_.centerX]));
    centerY_->setValue (static_cast<double> (row[spots_.centerY]));
    circrad_->setValue (static_cast<double> (row[spots_.circrad]));
    qualityMethod_->set_active (row[spots_.qualityMethod]);
    transit_->setValue (static_cast<double> (row[spots_.transit]));
    thresh_->setValue (static_cast<double> (row[spots_.thresh]));
    iter_->setValue (static_cast<double> (row[spots_.iter]));

    // Listener are reactivated
    disableParamlistener (false);

    updateParamVisibility ();
}

void ControlSpotPanel::controlspotChanged ()
{
    printf("controlspotChanged\n");
    load_ControlSpot_param();
}

void ControlSpotPanel::shapeChanged ()
{
    save_ControlSpot_param();

    printf("shapeChanged\n");
    const auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return;
    }
    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;
    updateControlSpotCurve (row);
}

void ControlSpotPanel::shapeMethodeChanged ()
{
    printf("shapeMethodeChanged\n");
    const int method = shapeMethod_->get_active_row_number ();

    if (method == 1 || method == 3) { // Symmetrical cases
        locXL_->setValue (locX_->getValue ());
        locYT_->setValue (locY_->getValue ());

        // Update associated control spot curve
        const auto s = treeview_.get_selection();
        if (!s->count_selected_rows()) {
            return;
        }
        const auto iter = s->get_selected();
        Gtk::TreeModel::Row row = *iter;
        save_ControlSpot_param ();
        updateControlSpotCurve (row);
    } else {
        save_ControlSpot_param ();
    }

    updateParamVisibility ();
}

void ControlSpotPanel::updateParamVisibility ()
{
    printf("updateParamVisibility\n");
    const int method = shapeMethod_->get_active_row_number ();

    if (method == 1 || method == 3) { // Symmetrical cases
        locXL_->hide ();
        locYT_->hide ();
        if (method == 1) { // 1 = Symmetrical (mouse)
            locX_->hide ();
            locY_->hide ();
            centerX_->hide ();
            centerY_->hide ();
        } else { // 3 = Symmetrical (mouse + sliders)
            locX_->show ();
            locY_->show ();
            centerX_->show ();
            centerY_->show ();
        }
    } else { // Independent cases
        if (method == 0) { // 0 = Independent (mouse)
            locX_->hide ();
            locXL_->hide ();
            locY_->hide ();
            locYT_->hide ();
            centerX_->hide ();
            centerY_->hide ();
        } else { // 2 = Independent (mouse + sliders)
            locX_->show ();
            locXL_->show ();
            locY_->show ();
            locYT_->show ();
            centerX_->show ();
            centerY_->show ();
        }
    }
}

void ControlSpotPanel::adjusterChanged(Adjuster* a, double newval)
{
    printf("adjusterChanged\n");
    const int method = shapeMethod_->get_active_row_number ();
    if (method == 1 || method == 3) { // Symmetrical cases
        locXL_->setValue (locX_->getValue ());
        locYT_->setValue (locY_->getValue ());
    }

    save_ControlSpot_param();

    // Update associated control spot curve
    const auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return;
    }
    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;
    updateControlSpotCurve (row);
}

void ControlSpotPanel::disableParamlistener(bool cond)
{
    printf("disableParamlistener: %d\n", cond);
    shapeconn_.block(cond);
    spotMethodconn_.block(cond);
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
}

void ControlSpotPanel::setParamEditable(bool cond)
{
    printf("setParamEditable: %d\n", cond);
    shape_->set_sensitive(cond);
    spotMethod_->set_sensitive(cond);
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
}

void ControlSpotPanel::addControlSpotCurve (Gtk::TreeModel::Row row)
{
    printf("addControlSpotCurve\n");
    if (row[spots_.curveid] > 0) { // Row has already an associated curve
        return;
    }

    // Creation of visibleGeometry
    Line* lineX;
    lineX = new Line ();
    lineX->innerLineWidth = 2.5;
    lineX->datum = Geometry::IMAGE;
    Line* lineXL;
    lineXL = new Line ();
    lineXL->innerLineWidth = 2.5;
    lineXL->datum = Geometry::IMAGE;
    Line* lineY;
    lineY = new Line ();
    lineY->innerLineWidth = 2.5;
    lineY->datum = Geometry::IMAGE;
    Line* lineYT;
    lineYT = new Line ();
    lineYT->innerLineWidth = 2.5;
    lineYT->datum = Geometry::IMAGE;
    Circle* centerCircle;
    centerCircle = new Circle ();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = true;
    Arcellipse* arc1;
    arc1 = new Arcellipse ();
    arc1->innerLineWidth = 0.7;
    arc1->datum = Geometry::IMAGE;
    arc1->radiusInImageSpace = true;
    Arcellipse* arc2;
    arc2 = new Arcellipse ();
    arc2->innerLineWidth = 0.7;
    arc2->datum = Geometry::IMAGE;
    arc2->radiusInImageSpace = true;
    Arcellipse* arc3;
    arc3 = new Arcellipse ();
    arc3->innerLineWidth = 0.7;
    arc3->datum = Geometry::IMAGE;
    arc3->radiusInImageSpace = true;
    Arcellipse* arc4;
    arc4 = new Arcellipse ();
    arc4->innerLineWidth = 0.7;
    arc4->datum = Geometry::IMAGE;
    arc4->radiusInImageSpace = true;
    Rectangle* rec;
    rec = new Rectangle ();
    rec->innerLineWidth = 0.7;
    rec->datum = Geometry::IMAGE;
    EditSubscriber::visibleGeometry.push_back(lineX); // (curveid - 1) * 10
    EditSubscriber::visibleGeometry.push_back(lineXL); // (curveid - 1) * 10 + 1
    EditSubscriber::visibleGeometry.push_back(lineY); // (curveid - 1) * 10 + 2
    EditSubscriber::visibleGeometry.push_back(lineYT); // (curveid - 1) * 10 + 3
    EditSubscriber::visibleGeometry.push_back(centerCircle); // (curveid - 1) * 10 + 4
    EditSubscriber::visibleGeometry.push_back(arc1); // (curveid - 1) * 10 + 5
    EditSubscriber::visibleGeometry.push_back(arc2); // (curveid - 1) * 10 + 6
    EditSubscriber::visibleGeometry.push_back(arc3); // (curveid - 1) * 10 + 7
    EditSubscriber::visibleGeometry.push_back(arc4); // (curveid - 1) * 10 + 8
    EditSubscriber::visibleGeometry.push_back(rec); // (curveid - 1) * 10 + 9

    // Creation of mouseOverGeometry
    lineX = new Line ();
    lineX->innerLineWidth = 2.5;
    lineX->datum = Geometry::IMAGE;
    lineXL = new Line ();
    lineXL->innerLineWidth = 2.5;
    lineXL->datum = Geometry::IMAGE;
    lineY = new Line ();
    lineY->innerLineWidth = 2.5;
    lineY->datum = Geometry::IMAGE;
    lineYT = new Line ();
    lineYT->innerLineWidth = 2.5;
    lineYT->datum = Geometry::IMAGE;
    centerCircle = new Circle ();
    centerCircle->datum = Geometry::IMAGE;
    centerCircle->radiusInImageSpace = true;
    arc1 = new Arcellipse ();
    arc1->innerLineWidth = 0.7;
    arc1->datum = Geometry::IMAGE;
    arc1->radiusInImageSpace = true;
    arc2 = new Arcellipse ();
    arc2->innerLineWidth = 0.7;
    arc2->datum = Geometry::IMAGE;
    arc2->radiusInImageSpace = true;
    arc3 = new Arcellipse ();
    arc3->innerLineWidth = 0.7;
    arc3->datum = Geometry::IMAGE;
    arc3->radiusInImageSpace = true;
    arc4 = new Arcellipse ();
    arc4->innerLineWidth = 0.7;
    arc4->datum = Geometry::IMAGE;
    arc4->radiusInImageSpace = true;
    rec = new Rectangle ();
    rec->innerLineWidth = 0.7;
    rec->datum = Geometry::IMAGE;
    EditSubscriber::mouseOverGeometry.push_back(lineX);  // (curveid - 1) * 10
    EditSubscriber::mouseOverGeometry.push_back(lineXL);  // (curveid - 1) * 10 + 1
    EditSubscriber::mouseOverGeometry.push_back(lineY);  // (curveid - 1) * 10 + 2
    EditSubscriber::mouseOverGeometry.push_back(lineYT);  // (curveid - 1) * 10 + 3
    EditSubscriber::mouseOverGeometry.push_back(centerCircle);  // (curveid - 1) * 10 + 4
    EditSubscriber::mouseOverGeometry.push_back(arc1);  // (curveid - 1) * 10 + 5
    EditSubscriber::mouseOverGeometry.push_back(arc2);  // (curveid - 1) * 10 + 6
    EditSubscriber::mouseOverGeometry.push_back(arc3);  // (curveid - 1) * 10 + 7
    EditSubscriber::mouseOverGeometry.push_back(arc4);  // (curveid - 1) * 10 + 8
    EditSubscriber::mouseOverGeometry.push_back(rec);  // (curveid - 1) * 10 + 9

    row[spots_.curveid] = EditSubscriber::visibleGeometry.size () / 10;
}

void ControlSpotPanel::updateControlSpotCurve (Gtk::TreeModel::Row row)
{
    const int curveid_ = static_cast<int> (row[spots_.curveid]);
    if (curveid_ == 0) { // Row has no associated curve
        return;
    }
    const int centerX_ = static_cast<int> (row[spots_.centerX]);
    const int centerY_ = static_cast<int> (row[spots_.centerY]);
    const int circrad_ = static_cast<int> (row[spots_.circrad]);
    const int locX_ = static_cast<int> (row[spots_.locX]);
    const int locXL_ = static_cast<int> (row[spots_.locXL]);
    const int locY_ = static_cast<int> (row[spots_.locY]);
    const int locYT_ = static_cast<int> (row[spots_.locYT]);
    const int shape_ = static_cast<int> (row[spots_.shape]);

    printf("updateControlSpotCurve: %d\n", curveid_);

    EditDataProvider* dataProvider = getEditProvider();

    if (!dataProvider) {
        return;
    }
    int imW = 0;
    int imH = 0;
    dataProvider->getImageSize(imW, imH);
    if (!imW || !imH) {
        return;
    }

    const double decayX = (locX_) * (double (imW)) / 2000.;
    const double decayXL = (locXL_) * (double (imW)) / 2000.;
    const double decayY = (locY_) * double (imH) / 2000.;
    const double decayYT = (locYT_) * double (imH) / 2000.;
    rtengine::Coord origin (imW / 2 + centerX_ * imW / 2000.f, imH / 2 + centerY_ * imH / 2000.f);

    const auto updateLineWithDecay = [&](Geometry * geometry, const float radius, const float decal, const float offSetAngle, const double decay) {
        const auto line = static_cast<Line*>(geometry); // 180
        line->begin = PolarCoord (radius, decal) + PolarCoord (decay, offSetAngle);
        line->begin += origin; // 0
        line->end = PolarCoord(radius, decal - 180) + PolarCoord(decay, offSetAngle);
        line->end += origin;
    };

    const auto updateCircle = [&](Geometry * geometry) {
        const auto circle = static_cast<Circle*> (geometry);
        circle->center = origin;
        circle->radius = circrad_;
    };

    const auto updateArcellipse = [&](Geometry * geometry, const double dRad_, const double dRad2_, const double begang_, const double endang_) {
        const auto arcellipse = static_cast<Arcellipse*> (geometry);
        arcellipse->center = origin;
        arcellipse->begang = begang_;
        arcellipse->endang = endang_;
        arcellipse->radius = dRad_;
        arcellipse->radius2 = dRad2_;
    };

    const auto updateRectangle = [&](Geometry * geometry) {
        const auto rectangle = static_cast<Rectangle*> (geometry);
        rectangle->bottomRight.x = origin.x + (int) decayX;
        rectangle->bottomRight.y = origin.y + (int) decayY;
        rectangle->topLeft.x = origin.x - (int) decayXL;
        rectangle->topLeft.y = origin.y - (int) decayYT;
    };

    updateLineWithDecay(visibleGeometry.at((curveid_ - 1) * 10), 500., 90., 0., decayX);
    updateLineWithDecay(mouseOverGeometry.at((curveid_ - 1) * 10), 500., 90., 0., decayX);

    updateLineWithDecay(visibleGeometry.at((curveid_ - 1) * 10 + 1), 500., 90., 180., decayXL);
    updateLineWithDecay(mouseOverGeometry.at((curveid_ - 1) * 10 + 1), 500., 90., 180., decayXL);

    updateLineWithDecay(visibleGeometry.at((curveid_ - 1) * 10 + 2), 500., 180., 90., decayY);
    updateLineWithDecay(mouseOverGeometry.at((curveid_ - 1) * 10 + 2), 500., 180., 90., decayY);

    updateLineWithDecay(visibleGeometry.at((curveid_ - 1) * 10 + 3), 500., 180., 270., decayYT);
    updateLineWithDecay(mouseOverGeometry.at((curveid_ - 1) * 10 + 3), 500., 180., 270., decayYT);

    updateCircle(visibleGeometry.at((curveid_ - 1) * 10 + 4));
    updateCircle(mouseOverGeometry.at((curveid_ - 1) * 10 + 4));

    updateArcellipse(visibleGeometry.at((curveid_ - 1) * 10 + 5), decayX, decayYT, 3*RT_PI_2, 2 * RT_PI);
    updateArcellipse(visibleGeometry.at((curveid_ - 1) * 10 + 6), decayXL, decayYT, RT_PI, 3*RT_PI_2);
    updateArcellipse(visibleGeometry.at((curveid_ - 1) * 10 + 7), decayXL, decayY, RT_PI_2, RT_PI);
    updateArcellipse(visibleGeometry.at((curveid_ - 1) * 10 + 8), decayX, decayY, 0., RT_PI_2);
    updateArcellipse(mouseOverGeometry.at((curveid_ - 1) * 10 + 5), decayX, decayYT, 3*RT_PI_2, 2 * RT_PI);
    updateArcellipse(mouseOverGeometry.at((curveid_ - 1) * 10 + 6), decayXL, decayYT, RT_PI, 3*RT_PI_2);
    updateArcellipse(mouseOverGeometry.at((curveid_ - 1) * 10 + 7), decayXL, decayY, RT_PI_2, RT_PI);
    updateArcellipse(mouseOverGeometry.at((curveid_ - 1) * 10 + 8), decayX, decayY, 0., RT_PI_2);

    updateRectangle(visibleGeometry.at((curveid_ - 1) * 10 + 9));
    updateRectangle(mouseOverGeometry.at((curveid_ - 1) * 10 + 9));

    // Update Arcellipse/Rectangle visibility according to shape
    if (shape_ == 0) { // 0 = Ellipse
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 5)->setActive(true); // arc1
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 6)->setActive(true); // arc2
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 7)->setActive(true); // arc3
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 8)->setActive(true); // arc4
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 9)->setActive(false); // rec

        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 5)->setActive(true); // arc1
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 6)->setActive(true); // arc2
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 7)->setActive(true); // arc3
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 8)->setActive(true); // arc4
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 9)->setActive(false); // rec
    } else { // 1 = Rectangle
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 5)->setActive(false); // arc1
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 6)->setActive(false); // arc2
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 7)->setActive(false); // arc3
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 8)->setActive(false); // arc4
        EditSubscriber::visibleGeometry.at((curveid_ - 1) * 10 + 9)->setActive(true); // rec

        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 5)->setActive(false); // arc1
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 6)->setActive(false); // arc2
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 7)->setActive(false); // arc3
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 8)->setActive(false); // arc4
        EditSubscriber::mouseOverGeometry.at((curveid_ - 1) * 10 + 9)->setActive(true); // rec
    }
}

void ControlSpotPanel::deleteControlSpotCurve (Gtk::TreeModel::Row row)
{
    const int curveid_ = static_cast<int> (row[spots_.curveid]);
    if (curveid_ == 0) { // Row has no associated curve
        return;
    }
    printf("deleteControlSpotCurve: %d\n", curveid_);

    // visibleGeometry
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 9);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 8);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 7);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 6);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 5);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 4);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 3);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 2);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10 + 1);
    EditSubscriber::visibleGeometry.erase(EditSubscriber::visibleGeometry.begin () + (curveid_ - 1) * 10);

    // mouseOverGeometry
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 9);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 8);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 7);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 6);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 5);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 4);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 3);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 2);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10 + 1);
    EditSubscriber::mouseOverGeometry.erase(EditSubscriber::mouseOverGeometry.begin () + (curveid_ - 1) * 10);

    row[spots_.curveid] = 0; // Reset associated curve id

    // Reordering curve id
    Gtk::TreeModel::Children children = treemodel_->children ();
    for (auto iter = children.begin (); iter != children.end (); iter++) {
        Gtk::TreeModel::Row r = *iter;
        if (r[spots_.curveid] > curveid_) {
            r[spots_.curveid] = r[spots_.curveid] - 1;
        }
    }
}

CursorShape ControlSpotPanel::getCursor (int objectID)
{
    printf ("Object ID: %d\n", objectID);
    int rem_ = objectID % 10;

    switch (rem_) {
    case (0): // LocX: (curveid_ - 1) * 10
        return CSMove1DH;
    case (1): // LocXL: (curveid_ - 1) * 10 + 1
        return CSMove1DH;
    case (2): // LocY: (curveid_ - 1) * 10 + 2
        return CSMove1DV;
    case (3): // LocYT: (curveid_ - 1) * 10 + 3
        return CSMove1DV;
    case (4): // centerCircle: (curveid_ - 1) * 10 + 4
        return CSMove2D;
    case (5): // arc1: (curveid_ - 1) * 10 + 5
        return CSMove2D;
    case (6): // arc2: (curveid_ - 1) * 10 + 6
        return CSMove2D;
    case (7): // arc3: (curveid_ - 1) * 10 + 7
        return CSMove2D;
    case (8): // arc4: (curveid_ - 1) * 10 + 8
        return CSMove2D;
    case (9): // rec: (curveid_ - 1) * 10 + 9
        return CSMove2D;
    default:
        return CSOpenHand;
    }
}

bool ControlSpotPanel::mouseOver (int modifierKey)
{
    EditDataProvider* editProvider_ = getEditProvider();
    if (!editProvider_) {
        return false;
    }

    int object_ = editProvider_->object;

    if (object_ != lastObject_) {
        if (object_ == -1) {
            for (int it_ = 0; it_ < (int) EditSubscriber::visibleGeometry.size (); it_++) {
                EditSubscriber::visibleGeometry.at(it_)->state = Geometry::NORMAL;
            }
            lastObject_ = object_;
            return false;
        }

        int curveId_ = object_ / 10 + 1;
        int rem = object_ % 10;
        for (int it_ = 0; it_ < (int) EditSubscriber::visibleGeometry.size (); it_++) {
            if ((it_ < ((curveId_ - 1) * 10)) || (it_ > ((curveId_ - 1) * 10) + 9)) { // it_ does not belong to cursor pointed curve
                EditSubscriber::visibleGeometry.at(it_)->state = Geometry::NORMAL;
            }
        }

        const int method = shapeMethod_->get_active_row_number ();

        // LocX
        if (rem == 0) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10)->state = Geometry::PRELIGHT;

            if (method == 1 || method == 3) { // Symmetrical cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 1)->state = Geometry::PRELIGHT;
            }
        } else {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10)->state = Geometry::NORMAL;
        }

        // LocXL
        if (rem == 1) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 1)->state = Geometry::PRELIGHT;

            if (method == 1 || method == 3) { // Symmetrical cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10)->state = Geometry::PRELIGHT;
            }
        } else {
            if (method == 0 || method == 2) { // Independent cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 1)->state = Geometry::NORMAL;
            }
        }

        // LocY
        if (rem == 2) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 2)->state = Geometry::PRELIGHT;

            if (method == 1 || method == 3) { // Symmetrical cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 3)->state = Geometry::PRELIGHT;
            }
        } else {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 2)->state = Geometry::NORMAL;
        }

        // LocYT
        if (rem == 3) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 3)->state = Geometry::PRELIGHT;

            if (method == 1 || method == 3) { // Symmetrical cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 2)->state = Geometry::PRELIGHT;
            }
        } else {
            if (method == 0 || method == 2) { // Independent cases
                EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 3)->state = Geometry::NORMAL;
            }
        }

        // Circle, Arcellipses and Rectangle
        if (rem >= 4) {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 4)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 5)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 6)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 7)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 8)->state = Geometry::PRELIGHT;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 9)->state = Geometry::PRELIGHT;
        } else {
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 4)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 5)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 6)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 7)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 8)->state = Geometry::NORMAL;
            EditSubscriber::visibleGeometry.at((curveId_ - 1) * 10 + 9)->state = Geometry::NORMAL;
        }

        lastObject_ = object_;
        return true;
    }

    return false;
}

bool ControlSpotPanel::button1Pressed (int modifierKey)
{
    printf("button1Pressed\n");
    EditDataProvider *provider = getEditProvider ();
    if (!provider || lastObject_ == -1) {
        return false;
    }

    // Select associated control spot
    int curveId_ = lastObject_ / 10 + 1;
    Gtk::TreeModel::Children children = treemodel_->children ();
    for (auto iter = children.begin (); iter != children.end (); iter++) {
        Gtk::TreeModel::Row r = *iter;
        if (r[spots_.curveid] == curveId_) {
            treeview_.set_cursor (treemodel_->get_path (r));
            break;
        }
    }

    lastCoord_->set (provider->posImage.x + provider->deltaImage.x, provider->posImage.y + provider->deltaImage.y);
    EditSubscriber::action = ES_ACTION_DRAGGING;
    return true;
}

bool ControlSpotPanel::button1Released ()
{
    printf("button1Released\n");
    EditSubscriber::action = ES_ACTION_NONE;
    return true;
}

bool ControlSpotPanel::drag1 (int modifierKey)
{
    printf("drag1\n");

    EditDataProvider *provider = getEditProvider ();
    if (!provider || lastObject_ == -1) {
        return false;
    }

    // Get associated control spot
    const auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return false;
    }
    const auto iter = s->get_selected();
    Gtk::TreeModel::Row row = *iter;

    int imW, imH;
    provider->getImageSize (imW, imH);
    int rem = lastObject_ % 10;
    int method = shapeMethod_->get_active_row_number ();
    Coord* newCoord = new Coord (provider->posImage.x + provider->deltaImage.x, provider->posImage.y + provider->deltaImage.y);

    // LocX
    if (rem == 0) {
        double deltaX = (double (newCoord->x) - double (lastCoord_->x)) * 2000. / double (imW);
        locX_->setValue (locX_->getValue () + deltaX);
        if (method == 1 || method == 3) { // Symmetrical cases
            locXL_->setValue (locX_->getValue ());
        }
        save_ControlSpot_param ();
        updateControlSpotCurve (row);
    }

    // LocXL
    if (rem == 1) {
        double deltaXL = (double (lastCoord_->x) - double (newCoord->x)) * 2000. / double (imW);
        locXL_->setValue (locXL_->getValue () + deltaXL);
        if (method == 1 || method == 3) { // Symmetrical cases
            locX_->setValue (locXL_->getValue ());
        }
        save_ControlSpot_param ();
        updateControlSpotCurve (row);
    }

    // LocY
    if (rem == 2) {
        double deltaY = (double (newCoord->y) - double (lastCoord_->y)) * 2000. / double (imH);
        locY_->setValue (locY_->getValue () + deltaY);
        if (method == 1 || method == 3) { // Symmetrical cases
            locYT_->setValue (locY_->getValue ());
        }
        save_ControlSpot_param ();
        updateControlSpotCurve (row);
    }

    // LocYT
    if (rem == 3) {
        double deltaYT = (double (lastCoord_->y) - double (newCoord->y)) * 2000. / double (imH);
        locYT_->setValue (locYT_->getValue () + deltaYT);
        if (method == 1 || method == 3) { // Symmetrical cases
            locY_->setValue (locYT_->getValue ());
        }
        save_ControlSpot_param ();
        updateControlSpotCurve (row);
    }

    // Circle, Arcellipses and Rectangle
    if (rem >= 4) {
        double deltaX = (double (newCoord->x) - double (lastCoord_->x)) * 2000. / double (imW);
        double deltaY = (double (newCoord->y) - double (lastCoord_->y)) * 2000. / double (imH);
        centerX_->setValue (centerX_->getValue () + deltaX);
        centerY_->setValue (centerY_->getValue () + deltaY);
        save_ControlSpot_param ();
        updateControlSpotCurve (row);
    }

    lastCoord_->set (newCoord->x, newCoord->y);
    return true;
}

//-----------------------------------------------------------------------------
// ControlSpots
//-----------------------------------------------------------------------------

ControlSpotPanel::ControlSpots::ControlSpots()
{
    add (id);
    add (name);
    add (isvisible);
    add (curveid);
    add (shape);
    add (spotMethod);
    add (shapeMethod);
    add (locX);
    add (locXL);
    add (locYT);
    add (locY);
    add (centerX);
    add (centerY);
    add (circrad);
    add (qualityMethod);
    add (transit);
    add (thresh);
    add (iter);
}

//-----------------------------------------------------------------------------
// RenameDialog
//-----------------------------------------------------------------------------

ControlSpotPanel::RenameDialog::RenameDialog(const Glib::ustring &actualname, Gtk::Window &parent):
    Gtk::Dialog ("Renaming Control Spot", parent)
{
    Gtk::HBox *hb = Gtk::manage (new Gtk::HBox());
    hb->pack_start (*Gtk::manage (new Gtk::Label ("Enter the new Control Spot name")), false, false, 4);

    newname_.set_text(actualname);
    hb->pack_start(newname_);

    get_content_area()->pack_start (*hb, Gtk::PACK_SHRINK, 4);

    add_button (M ("GENERAL_OK"), 1);
    add_button (M ("GENERAL_CANCEL"), 2);

    show_all_children();
}

Glib::ustring ControlSpotPanel::RenameDialog::get_new_name()
{
    return newname_.get_text();
}
