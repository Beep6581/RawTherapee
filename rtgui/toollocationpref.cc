/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2021 Lawrence Lee
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <unordered_set>

#include "guiutils.h"
#include "options.h"
#include "rtscalable.h"
#include "toollocationpref.h"
#include "toolpanelcoord.h"

using Tool = ToolPanelCoordinator::Tool;
using Favorites = std::unordered_set<Tool, ScopedEnumHash>;

namespace
{

/**
 * Returns the language key for the panel's title.
 */
Glib::ustring getToolPanelTitleKey(ToolPanelCoordinator::Panel panel)
{
    switch (panel) {
        case ToolPanelCoordinator::Panel::FAVORITE:
            return "MAIN_TAB_FAVORITES";
        case ToolPanelCoordinator::Panel::EXPOSURE:
            return "MAIN_TAB_EXPOSURE";
        case ToolPanelCoordinator::Panel::DETAILS:
            return "MAIN_TAB_DETAIL";
        case ToolPanelCoordinator::Panel::COLOR:
            return "MAIN_TAB_COLOR";
        case ToolPanelCoordinator::Panel::ADVANCED:
            return "MAIN_TAB_ADVANCED";
        case ToolPanelCoordinator::Panel::LOCALLAB:
            return "MAIN_TAB_LOCALLAB";
        case ToolPanelCoordinator::Panel::TRANSFORM_PANEL:
            return "MAIN_TAB_TRANSFORM";
        case ToolPanelCoordinator::Panel::RAW:
            return "MAIN_TAB_RAW";
    }
    assert(false);
    return "";
}

/**
 * Returns the language key for the tool's title.
 */
Glib::ustring getToolTitleKey(Tool tool)
{
    using Tool = Tool;
    switch (tool) {
        case Tool::TONE_CURVE:
            return "TP_EXPOSURE_LABEL";
        case Tool::SHADOWS_HIGHLIGHTS:
            return "TP_SHADOWSHLIGHTS_LABEL";
        case Tool::TONE_EQUALIZER:
            return "TP_TONE_EQUALIZER_LABEL";
        case Tool::IMPULSE_DENOISE:
            return "TP_IMPULSEDENOISE_LABEL";
        case Tool::DEFRINGE_TOOL:
            return "TP_DEFRINGE_LABEL";
        case Tool::COMPRESSGAMUT_TOOL:
            return "TP_COMPRESSGAMUT_LABEL";
        case Tool::SPOT:
            return "TP_SPOT_LABEL";
        case Tool::DIR_PYR_DENOISE:
            return "TP_DIRPYRDENOISE_LABEL";
        case Tool::EPD:
            return "TP_EPD_LABEL";
        case Tool::SHARPENING_TOOL:
            return "TP_SHARPENING_LABEL";
        case Tool::LOCAL_CONTRAST:
            return "TP_LOCALCONTRAST_LABEL";
        case Tool::SHARPEN_EDGE:
            return "TP_SHARPENEDGE_LABEL";
        case Tool::SHARPEN_MICRO:
            return "TP_SHARPENMICRO_LABEL";
        case Tool::L_CURVE:
            return "TP_LABCURVE_LABEL";
        case Tool::RGB_CURVES:
            return "TP_RGBCURVES_LABEL";
        case Tool::COLOR_TONING:
            return "TP_COLORTONING_LABEL";
        case Tool::LENS_GEOM:
            return "TP_LENSGEOM_LABEL";
        case Tool::LENS_PROF:
            return "TP_LENSPROFILE_LABEL";
        case Tool::DISTORTION:
            return "TP_DISTORTION_LABEL";
        case Tool::ROTATE:
            return "TP_ROTATE_LABEL";
        case Tool::VIBRANCE:
            return "TP_VIBRANCE_LABEL";
        case Tool::COLOR_APPEARANCE:
            return "TP_COLORAPP_LABEL";
        case Tool::WHITE_BALANCE:
            return "TP_WBALANCE_LABEL";
        case Tool::VIGNETTING:
            return "TP_VIGNETTING_LABEL";
        case Tool::RETINEX_TOOL:
            return "TP_RETINEX_LABEL";
        case Tool::GRADIENT:
            return "TP_GRADIENT_LABEL";
        case Tool::LOCALLAB:
            return "TP_LOCALLAB_LABEL";
        case Tool::PC_VIGNETTE:
            return "TP_PCVIGNETTE_LABEL";
        case Tool::PERSPECTIVE:
            return "TP_PERSPECTIVE_LABEL";
        case Tool::CA_CORRECTION:
            return "TP_CACORRECTION_LABEL";
        case Tool::CH_MIXER:
            return "TP_CHMIXER_LABEL";
        case Tool::BLACK_WHITE:
            return "TP_BWMIX_LABEL";
        case Tool::RESIZE_TOOL:
            return "TP_RESIZE_LABEL";
        case Tool::PR_SHARPENING:
            return "TP_PRSHARPENING_LABEL";
        case Tool::FRAMING:
            return "TP_FRAMING_LABEL";
        case Tool::CROP_TOOL:
            return "TP_CROP_LABEL";
        case Tool::ICM:
            return "TP_ICM_LABEL";
        case Tool::WAVELET:
            return "TP_WAVELET_LABEL";
        case Tool::DIR_PYR_EQUALIZER:
            return "TP_DIRPYREQUALIZER_LABEL";
        case Tool::HSV_EQUALIZER:
            return "TP_HSVEQUALIZER_LABEL";
        case Tool::FILM_SIMULATION:
            return "TP_FILMSIMULATION_LABEL";
        case Tool::SOFT_LIGHT:
            return "TP_SOFTLIGHT_LABEL";
        case Tool::DEHAZE:
            return "TP_DEHAZE_LABEL";
        case Tool::SENSOR_BAYER:
            return "TP_RAW_SENSOR_BAYER_LABEL";
        case Tool::SENSOR_XTRANS:
            return "TP_RAW_SENSOR_XTRANS_LABEL";
        case Tool::BAYER_PROCESS:
            return "TP_RAW_LABEL";
        case Tool::XTRANS_PROCESS:
            return "TP_RAW_LABEL";
        case Tool::BAYER_PREPROCESS:
            return "TP_PREPROCESS_LABEL";
        case Tool::PREPROCESS:
            return "TP_PREPROCESS_LABEL";
        case Tool::DARKFRAME_TOOL:
            return "TP_DARKFRAME_LABEL";
        case Tool::FLATFIELD_TOOL:
            return "TP_FLATFIELD_LABEL";
        case Tool::RAW_CA_CORRECTION:
            return "TP_RAWCACORR_LABEL";
        case Tool::RAW_EXPOSURE:
            return "TP_EXPOS_WHITEPOINT_LABEL";
        case Tool::PREPROCESS_WB:
            return "TP_PREPROCWB_LABEL";
        case Tool::BAYER_RAW_EXPOSURE:
            return "TP_EXPOS_BLACKPOINT_LABEL";
        case Tool::XTRANS_RAW_EXPOSURE:
            return "TP_EXPOS_BLACKPOINT_LABEL";
        case Tool::FATTAL:
            return "TP_TM_FATTAL_LABEL";
        case Tool::FILM_NEGATIVE:
            return "TP_FILMNEGATIVE_LABEL";
        case Tool::PD_SHARPENING:
            return "TP_PDSHARPENING_LABEL";
    };
    assert(false);
    return "";
}

/**
 * A widget with buttons (packed vertically) for modifying a list store with a
 * tree view.
 *
 * Includes buttons for moving single rows up or down and a button for removing
 * selected rows.
 */
class ListEditButtons : public Gtk::Box
{
private:
    Gtk::TreeView &list;
    Glib::RefPtr<Gtk::ListStore> listStore;
    Gtk::Button buttonUp;
    Gtk::Button buttonDown;
    Gtk::Button buttonRemove;

    sigc::signal<void, const std::vector<Gtk::TreeModel::Path> &> signalRowsPreErase;

    void onButtonDownPressed();
    void onButtonRemovePressed();
    void onButtonUpPressed();
    void onListSelectionChanged();
    void updateButtonSensitivity();

public:
    /**
     * Constructs an edit buttons widget for modifying the provided list store.
     *
     * @param list The tree view for which selections are made in. The tree
     * view's model MUST be the list store.
     * @param listStore The list store that the widget will modify.
     */
    explicit ListEditButtons(Gtk::TreeView &list, Glib::RefPtr<Gtk::ListStore> listStore);

    /**
     * Returns the signal that gets emitted right before this widget removes
     * rows from its list store.
     *
     * The signal contains a vector of tree model paths of the rows that will be
     * erased.
     */
    sigc::signal<void, const std::vector<Gtk::TreeModel::Path> &> getSignalRowsPreErase() const;
};

/**
 * Model columns for the favorites list.
 */
class FavoritesColumns : public Gtk::TreeModelColumnRecord
{
public:
    /** The tool's display name. */
    Gtk::TreeModelColumn<Glib::ustring> toolName;
    /** The tool. */
    Gtk::TreeModelColumn<Tool> tool;

    FavoritesColumns()
    {
        add(toolName);
        add(tool);
    }
};

/**
 * Model columns for the available tools list.
 */
class ToolListColumns : public Gtk::TreeModelColumnRecord
{
public:
    /** The tool's display name. */
    Gtk::TreeModelColumn<Glib::ustring> toolName;
    /** The tool. */
    Gtk::TreeModelColumn<Tool> tool;
    /** Is the tool added to favorites. */
    Gtk::TreeModelColumn<bool> isFavorite;
    /** Can the tool be added to favorites. */
    Gtk::TreeModelColumn<bool> isEditable;

    ToolListColumns()
    {
        add(toolName);
        add(tool);
        add(isFavorite);
        add(isEditable);
    }
};

ListEditButtons::ListEditButtons(Gtk::TreeView &list, Glib::RefPtr<Gtk::ListStore> listStore) :
    Gtk::Box(Gtk::Orientation::ORIENTATION_VERTICAL),
    list(list),
    listStore(listStore)
{
    assert(list.get_model() == listStore);

    // Set button images.
    buttonUp.set_image_from_icon_name("arrow-up-small");
    buttonDown.set_image_from_icon_name("arrow-down-small");
    buttonRemove.set_image_from_icon_name("remove-small");

    // Connect signals for changing button sensitivity.
    const auto on_list_sel_changed_fun = sigc::mem_fun(
        *this, &ListEditButtons::onListSelectionChanged);
    const auto on_row_deleted_fun = sigc::hide(on_list_sel_changed_fun);
    const auto on_row_inserted_fun = sigc::hide(on_row_deleted_fun);
    list.get_selection()->signal_changed().connect(on_list_sel_changed_fun);
    listStore->signal_row_deleted().connect(on_row_deleted_fun);
    listStore->signal_row_inserted().connect(on_row_inserted_fun);

    // Connect signals for buttons.
    buttonUp.signal_pressed().connect(sigc::mem_fun(
        *this, &ListEditButtons::onButtonUpPressed));
    buttonDown.signal_pressed().connect(sigc::mem_fun(
        *this, &ListEditButtons::onButtonDownPressed));
    buttonRemove.signal_pressed().connect(sigc::mem_fun(
        *this, &ListEditButtons::onButtonRemovePressed));

    updateButtonSensitivity();

    add(buttonUp);
    add(buttonDown);
    add(buttonRemove);
}

void ListEditButtons::onButtonDownPressed()
{
    const auto list_children = listStore->children();
    const std::vector<Gtk::TreeModel::Path> selected =
        list.get_selection()->get_selected_rows();

    if (selected.size() != 1) { // Only one can be selected.
        return;
    }

    // Get the selected row and next row.
    const auto selected_row_iter = listStore->get_iter(selected[0]);
    auto next_row_iter = selected_row_iter;
    next_row_iter++;

    if (next_row_iter == list_children.end()) { // Can't be last row.
        return;
    }

    // Move the selected row down and update the buttons.
    listStore->iter_swap(selected_row_iter, next_row_iter);
    updateButtonSensitivity();
}

void ListEditButtons::onButtonRemovePressed()
{
    const std::vector<Gtk::TreeModel::Path> selected_paths =
        list.get_selection()->get_selected_rows();
    std::vector<Gtk::TreeModel::RowReference> selected(selected_paths.size());

    // Get row references, which are valid until the row is removed.
    std::transform(
        selected_paths.begin(),
        selected_paths.end(),
        selected.begin(),
        [this](const Gtk::TreeModel::Path &row_path) {
            return Gtk::TreeModel::RowReference(listStore, row_path);
        });

    signalRowsPreErase.emit(selected_paths);

    // Remove the selected rows.
    for (const auto &row_ref : selected) {
        const auto row_path = row_ref.get_path();
        if (row_path) {
            listStore->erase(listStore->get_iter(row_path));
        } else if (rtengine::settings->verbose) {
            std::cout << "Unable to remove row because it does not exist anymore." << std::endl;
        }
    }

    updateButtonSensitivity();
}

void ListEditButtons::onButtonUpPressed()
{
    const auto list_children = listStore->children();
    const std::vector<Gtk::TreeModel::Path> selected =
        list.get_selection()->get_selected_rows();

    if (selected.size() != 1) { // Only one can be selected.
        return;
    }

    const auto selected_row_iter = listStore->get_iter(selected[0]);

    if (selected_row_iter == list_children.begin()) { // Can't be first row.
        return;
    }

    // Swap selected row with the previous row.
    auto prev_row_iter = selected_row_iter;
    prev_row_iter--;
    listStore->iter_swap(selected_row_iter, prev_row_iter);
    updateButtonSensitivity();
}

void ListEditButtons::onListSelectionChanged()
{
    updateButtonSensitivity();
}

void ListEditButtons::updateButtonSensitivity()
{
    assert(list.get_model() == listStore);

    const std::vector<Gtk::TreeModel::Path> selected =
        list.get_selection()->get_selected_rows();

    // Update sensitivity of the move up/down buttons.
    if (selected.size() != 1) {
        // Items can only be moved if one row is selected.
        buttonDown.set_sensitive(false);
        buttonUp.set_sensitive(false);
    } else {
        // Up button cannot be used on the first row. Down button cannot be used
        // on the last row.
        auto selected_row_iter = list.get_model()->get_iter(selected[0]);
        const auto list_children = listStore->children();
        buttonUp.set_sensitive(!selected_row_iter->equal(list_children.begin()));
        buttonDown.set_sensitive(!(++selected_row_iter)->equal(list_children.end()));
    }

    // Update sensitivity of the remove button.
    buttonRemove.set_sensitive(selected.size() > 0);
}

sigc::signal<void, const std::vector<Gtk::TreeModel::Path> &>
ListEditButtons::getSignalRowsPreErase() const
{
    return signalRowsPreErase;
}

}

struct ToolLocationPreference::Impl {
    Options &options;

    // General options.
    Gtk::CheckButton *cloneFavoriteToolsToggleWidget;

    // Tool list.
    ToolListColumns toolListColumns;
    Glib::RefPtr<Gtk::TreeStore> toolListModelPtr;
    Gtk::CellRendererToggle toolListCellRendererFavorite;
    Gtk::CellRendererText toolListCellRendererToolName;
    Gtk::TreeViewColumn toolListViewColumnFavorite;
    Gtk::TreeViewColumn toolListViewColumnToolName;
    Gtk::TreeView *toolListViewPtr;
    std::unordered_map<Tool, Gtk::TreeModel::iterator, ScopedEnumHash>
        toolListToolToRowIterMap;

    // Favorites list.
    FavoritesColumns favoritesColumns;
    Glib::RefPtr<Gtk::ListStore> favoritesModelPtr;
    Gtk::CellRendererText favoritesCellRendererToolName;
    Gtk::TreeViewColumn favoritesViewColumnToolName;
    Gtk::TreeView *favoritesViewPtr;
    ListEditButtons favoritesListEditButtons;

    /**
     * Constructs an implementation that gets values from and updates the
     * provided options object.
     */
    explicit Impl(Options &options);

    /**
     * Adds the tools in the tool tree as a child in the provided row.
     *
     * @param tools The tool tree.
     * @param parentRowIter An iterator to the row under which to add the tools.
     * @param favorites The tools which are currently marked as favorites.
     */
    void addToolListRowGroup(
        const std::vector<ToolPanelCoordinator::ToolTree> &tools,
        const Gtk::TreeIter &parentRowIter,
        const Favorites &favorites);
    /**
     * Toggles the tool list favorite column and updates the favorites list.
     *
     * @param row_path Path to the tool list model row.
     */
    void favoriteToggled(const Glib::ustring &row_path);
    /**
     * Initializes the favorites list.
     *
     * @param favorites Tools that are currently marked as favorites.
     */
    void initFavoritesRows(const std::vector<Tool> &favorites);
    /**
     * Initializes the available tools list.
     *
     * @param favorites Tools that are currently marked as favorites.
     */
    void initToolListRows(const std::vector<Tool> &favorites);
    /**
     * Updates the favorites column of the available tools list when tools are
     * about to be removed from the favorites list.
     *
     * @param paths Paths in the favorites list pointing to the rows that are
     * about to be removed.
     */
    void onFavoritesRowsPreRemove(const std::vector<Gtk::TreeModel::Path> &paths);
    /**
     * Converts tool names to their corresponding tools.
     *
     * @param tool_names The tool names that need to be converted.
     * @return The tools.
     */
    std::vector<Tool> toolNamesToTools(
        const std::vector<Glib::ustring> &tool_names) const;
    /**
     * Updates the options object associated with this object with the current
     * favorites preferences.
     */
    void updateOptions();
};

ToolLocationPreference::Impl::Impl(Options &options) :
    options(options),

    // General options.
    cloneFavoriteToolsToggleWidget(Gtk::manage(
        new Gtk::CheckButton(M("PREFERENCES_TOOLPANEL_CLONE_FAVORITES")))),

    // Tool list.
    toolListModelPtr(Gtk::TreeStore::create(toolListColumns)),
    toolListViewColumnFavorite(
        Gtk::TreeViewColumn(M("PREFERENCES_TOOLPANEL_FAVORITE"))),
    toolListViewColumnToolName(
        Gtk::TreeViewColumn(M("PREFERENCES_TOOLPANEL_TOOL"))),
    toolListViewPtr(Gtk::manage(new Gtk::TreeView(toolListModelPtr))),

    // Favorites list.
    favoritesModelPtr(Gtk::ListStore::create(favoritesColumns)),
    favoritesViewColumnToolName(
        Gtk::TreeViewColumn(M("PREFERENCES_TOOLPANEL_TOOL"))),
    favoritesViewPtr(Gtk::manage(new Gtk::TreeView(favoritesModelPtr))),
    favoritesListEditButtons(*favoritesViewPtr, favoritesModelPtr)
{
    const std::vector<Tool> favorites = toolNamesToTools(options.favorites);

    // General options.
    cloneFavoriteToolsToggleWidget->set_active(options.cloneFavoriteTools);
    cloneFavoriteToolsToggleWidget->set_tooltip_text(
        M("PREFERENCES_TOOLPANEL_CLONE_FAVORITES_TOOLTIP"));

    // Tool list.
    toolListViewPtr->append_column(toolListViewColumnToolName);
    toolListViewColumnToolName.pack_start(toolListCellRendererToolName);
    toolListViewColumnToolName.set_renderer(
        toolListCellRendererToolName, toolListColumns.toolName);
    toolListViewPtr->append_column(toolListViewColumnFavorite);
    toolListViewColumnFavorite.pack_start(toolListCellRendererFavorite, false);
    toolListViewColumnFavorite.set_renderer(
        toolListCellRendererFavorite, toolListColumns.isFavorite);
    toolListViewColumnFavorite.add_attribute(
        toolListCellRendererFavorite.property_visible(), toolListColumns.isEditable);
    toolListCellRendererFavorite.signal_toggled().connect(
        sigc::mem_fun(*this, &ToolLocationPreference::Impl::favoriteToggled));
    initToolListRows(favorites);
    toolListViewPtr->expand_all();

    // Favorites list.
    favoritesViewPtr->append_column(favoritesViewColumnToolName);
    favoritesViewPtr->set_reorderable(true);
    favoritesViewColumnToolName.pack_start(favoritesCellRendererToolName);
    favoritesViewColumnToolName.set_renderer(
        favoritesCellRendererToolName, favoritesColumns.toolName);
    favoritesListEditButtons.getSignalRowsPreErase().connect(sigc::mem_fun(
        *this, &ToolLocationPreference::Impl::onFavoritesRowsPreRemove));
    favoritesViewPtr->get_selection()->set_mode(
        Gtk::SelectionMode::SELECTION_MULTIPLE);
    initFavoritesRows(favorites);
}

void ToolLocationPreference::Impl::favoriteToggled(const Glib::ustring &row_path)
{
    auto row_iter = toolListModelPtr->get_iter(row_path);
    const bool is_favorite = !row_iter->get_value(toolListColumns.isFavorite);
    const Tool tool = row_iter->get_value(toolListColumns.tool);

    // Update favorite column in the tool list.
    row_iter->set_value(toolListColumns.isFavorite, is_favorite);

    // Update the favorites list.
    if (is_favorite) {
        // Add to favorites list.
        const auto new_favorite_row_iter = favoritesModelPtr->append();
        new_favorite_row_iter->set_value(
            favoritesColumns.toolName,
            M(getToolTitleKey(tool)));
        new_favorite_row_iter->set_value(favoritesColumns.tool, tool);
    } else {
        // Remove from favorites list.
        const auto favorites_rows = favoritesModelPtr->children();
        auto row = favorites_rows.begin();
        while (
            row != favorites_rows.end() &&
            row->get_value(favoritesColumns.tool) != tool) {
            row++;
        }
        if (row != favorites_rows.end()) {
            favoritesModelPtr->erase(row);
        }
    }
}

void ToolLocationPreference::Impl::initFavoritesRows(
    const std::vector<Tool> &favorites)
{
    // Add the favorites to the favorites list store.
    for (const auto tool : favorites) {
        const auto favorite_row_iter = favoritesModelPtr->append();
        favorite_row_iter->set_value(
            favoritesColumns.toolName,
            M(getToolTitleKey(tool)));
        favorite_row_iter->set_value(favoritesColumns.tool, tool);
    }
}

void ToolLocationPreference::Impl::addToolListRowGroup(
    const std::vector<ToolPanelCoordinator::ToolTree> &tools,
    const Gtk::TreeIter &parentRowIter,
    const Favorites &favorites)
{
    // Recursively add the tool and its children to the tool list tree store.
    for (const ToolPanelCoordinator::ToolTree &tool : tools) {
        const auto tool_row_iter =
            toolListModelPtr->append(parentRowIter->children());
        tool_row_iter->set_value(
            toolListColumns.toolName,
            M(getToolTitleKey(tool.id)));
        tool_row_iter->set_value(toolListColumns.tool, tool.id);
        tool_row_iter->set_value(
            toolListColumns.isFavorite,
            favorites.count(tool.id) > 0);
        tool_row_iter->set_value(
            toolListColumns.isEditable,
            ToolPanelCoordinator::isFavoritable(tool.id));
        toolListToolToRowIterMap[tool.id] = tool_row_iter;
        addToolListRowGroup(tool.children, tool_row_iter, favorites);
    }
};

void ToolLocationPreference::Impl::initToolListRows(const std::vector<Tool> &favorites)
{
    const auto panel_tools = ToolPanelCoordinator::getDefaultToolLayout();
    Favorites favorites_set;

    // Convert the favorites vector into a set for fast lookup.
    for (const auto &tool : favorites) {
        favorites_set.insert(tool);
    }

    // Add each panel and their children to the tool list.
    for (const auto panel : {
             ToolPanelCoordinator::Panel::EXPOSURE,
             ToolPanelCoordinator::Panel::DETAILS,
             ToolPanelCoordinator::Panel::COLOR,
             ToolPanelCoordinator::Panel::ADVANCED,
             ToolPanelCoordinator::Panel::LOCALLAB,
             ToolPanelCoordinator::Panel::TRANSFORM_PANEL,
             ToolPanelCoordinator::Panel::RAW,
         }) {
        const auto tool_group_iter = toolListModelPtr->append();
        tool_group_iter->set_value(
            toolListColumns.toolName,
            M(getToolPanelTitleKey(panel)));
        addToolListRowGroup(panel_tools.at(panel), tool_group_iter, favorites_set);
    }
}

void ToolLocationPreference::Impl::onFavoritesRowsPreRemove(
    const std::vector<Gtk::TreeModel::Path> &paths)
{
    // Unset the favorite column in the tools list for tools about to be removed
    // from the favorites list.
    for (const auto &path : paths) {
        const auto &row_iter = toolListToolToRowIterMap.at(
            favoritesModelPtr->get_iter(path)->get_value(favoritesColumns.tool));
        row_iter->set_value(toolListColumns.isFavorite, false);
    }
}

std::vector<Tool> ToolLocationPreference::Impl::toolNamesToTools(
    const std::vector<Glib::ustring> &tool_names) const
{
    std::vector<Tool> tool_set;

    for (const auto &tool_name : tool_names) {
        Tool tool;
        try {
            tool = ToolPanelCoordinator::getToolFromName(tool_name);
        } catch (const std::out_of_range &e) {
            if (rtengine::settings->verbose) {
                std::cerr
                    << "Unrecognized tool name \"" << tool_name << "\"."
                    << std::endl;
            }
            assert(false);
            continue;
        }
        tool_set.push_back(tool);
    }

    return tool_set;
}

void ToolLocationPreference::Impl::updateOptions()
{
    options.cloneFavoriteTools = cloneFavoriteToolsToggleWidget->get_active();

    const auto favorites_rows = favoritesModelPtr->children();
    options.favorites.resize(favorites_rows.size());
    for (Gtk::TreeNodeChildren::size_type i = 0; i < favorites_rows.size(); i++) {
        const Tool tool = favorites_rows[i].get_value(favoritesColumns.tool);
        options.favorites[i] = ToolPanelCoordinator::getToolName(tool);
    }
}

ToolLocationPreference::ToolLocationPreference(Options &options) :
    impl(new Impl(options))
{
    // Layout grid.
    Gtk::Grid *layout_grid = Gtk::manage(new Gtk::Grid());
    layout_grid->set_column_spacing(4);
    layout_grid->set_row_spacing(4);
    layout_grid->set_column_homogeneous(true);
    pack_start(*layout_grid);

    // Tool list.
    Gtk::Frame *tool_list_frame = Gtk::manage(new Gtk::Frame(
        M("PREFERENCES_TOOLPANEL_AVAILABLETOOLS")));
    Gtk::ScrolledWindow *tool_list_scrolled_window =
        Gtk::manage(new Gtk::ScrolledWindow());
    tool_list_scrolled_window->set_min_content_width(RTScalable::scalePixelSize(400));
    layout_grid->attach_next_to(*tool_list_frame, Gtk::PositionType::POS_RIGHT, 1, 1);
    tool_list_frame->add(*tool_list_scrolled_window);
    tool_list_scrolled_window->add(*impl->toolListViewPtr);
    setExpandAlignProperties(
        tool_list_frame, true, true, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    // Favorites list.
    Gtk::Frame *favorites_frame = Gtk::manage(new Gtk::Frame(
        M("PREFERENCES_TOOLPANEL_FAVORITESPANEL")));
    Gtk::Box *favorites_box = Gtk::manage(new Gtk::Box());
    Gtk::ScrolledWindow *favorites_list_scrolled_window =
        Gtk::manage(new Gtk::ScrolledWindow());
    favorites_list_scrolled_window->set_min_content_width(RTScalable::scalePixelSize(400));
    layout_grid->attach_next_to(*favorites_frame, Gtk::PositionType::POS_RIGHT, 1, 1);
    favorites_box->pack_start(impl->favoritesListEditButtons, false, false);
    favorites_box->pack_start(*favorites_list_scrolled_window, true, true);
    favorites_frame->add(*favorites_box);
    favorites_list_scrolled_window->add(*impl->favoritesViewPtr);
    setExpandAlignProperties(
        favorites_frame, true, true, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    // General options.
    layout_grid->attach_next_to(
        *impl->cloneFavoriteToolsToggleWidget, Gtk::PositionType::POS_BOTTOM, 2, 1);
}

void ToolLocationPreference::updateOptions()
{
    impl->updateOptions();
}
