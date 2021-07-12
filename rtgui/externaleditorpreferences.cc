/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2021 Lawrence Lee <billee@ucdavis.edu>
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
#include "externaleditorpreferences.h"
#include "multilangmgr.h"
#include "rtimage.h"


ExternalEditorPreferences::ExternalEditorPreferences():
    Box(Gtk::Orientation::ORIENTATION_VERTICAL),
    list_model(Gtk::ListStore::create(model_columns)),
    toolbar(Gtk::Orientation::ORIENTATION_HORIZONTAL)
{
    // List view.
    list_view = Gtk::make_managed<Gtk::TreeView>();
    list_view->set_model(list_model);
    list_view->append_column(*Gtk::manage(makeAppColumn()));
    list_view->append_column(*Gtk::manage(makeCommandColumn()));

    for (auto &&column : list_view->get_columns()) {
        column->set_sizing(Gtk::TreeViewColumnSizing::TREE_VIEW_COLUMN_FIXED);
    }

    list_view->set_grid_lines(Gtk::TREE_VIEW_GRID_LINES_VERTICAL);
    list_view->set_reorderable();

    // List scroll area.
    list_scroll_area.set_hexpand();
    list_scroll_area.set_vexpand();
    list_scroll_area.add(*list_view);

    // Toolbar buttons.
    auto add_image = Gtk::make_managed<RTImage>("add-small.png");
    auto remove_image = Gtk::make_managed<RTImage>("remove-small.png");
    button_add = Gtk::make_managed<Gtk::Button>();
    button_remove = Gtk::make_managed<Gtk::Button>();
    button_add->set_image(*add_image);
    button_remove->set_image(*remove_image);
    button_app_chooser = Gtk::make_managed<Gtk::Button>(M("PREFERENCES_EXTERNALEDITOR_CHANGE"));

    button_app_chooser->signal_pressed().connect(sigc::mem_fun(
                *this, &ExternalEditorPreferences::openAppChooserDialog));
    button_add->signal_pressed().connect(sigc::mem_fun(
            *this, &ExternalEditorPreferences::addEditor));
    button_remove->signal_pressed().connect(sigc::mem_fun(
            *this, &ExternalEditorPreferences::removeSelectedEditors));

    list_view->get_selection()->signal_changed().connect(sigc::mem_fun(
                *this, &ExternalEditorPreferences::updateToolbarSensitivity));
    updateToolbarSensitivity();

    // Toolbar.
    toolbar.set_halign(Gtk::Align::ALIGN_END);
    toolbar.add(*button_app_chooser);
    toolbar.add(*button_add);
    toolbar.add(*button_remove);

    // This widget's children.
    add(list_scroll_area);
    add(toolbar);
    show_all();
}

std::vector<ExternalEditorPreferences::EditorInfo>
ExternalEditorPreferences::getEditors() const
{
    std::vector<ExternalEditorPreferences::EditorInfo> editors;

    auto children = list_model->children();

    for (auto rowIter = children.begin(); rowIter != children.end(); rowIter++) {
        const Gio::Icon *const icon = rowIter->get_value(model_columns.icon).get();
        const auto &icon_name = icon == nullptr ? "" : icon->to_string();
        editors.push_back(ExternalEditorPreferences::EditorInfo(
                              rowIter->get_value(model_columns.name),
                              rowIter->get_value(model_columns.command),
                              icon_name,
                              rowIter->get_value(model_columns.other_data)
                          ));
    }

    return editors;
}

void ExternalEditorPreferences::setEditors(
    const std::vector<ExternalEditorPreferences::EditorInfo> &editors)
{
    list_model->clear();

    for (const ExternalEditorPreferences::EditorInfo & editor : editors) {
        auto row = *list_model->append();
        row[model_columns.name] = editor.name;
        row[model_columns.icon] = editor.icon_name.empty() ? Glib::RefPtr<Gio::Icon>() : Gio::Icon::create(editor.icon_name);
        row[model_columns.command] = editor.command;
        row[model_columns.other_data] = editor.other_data;
    }
}

void ExternalEditorPreferences::addEditor()
{
    Gtk::TreeModel::Row row;
    auto selected = list_view->get_selection()->get_selected_rows();

    if (selected.size()) {
        row = *list_model->insert_after(list_model->get_iter(selected.back()));
    } else {
        row = *list_model->append();
    }

    row[model_columns.name] = "-";
    list_view->get_selection()->select(row);
}

Gtk::TreeViewColumn *ExternalEditorPreferences::makeAppColumn()
{
    auto name_renderer = Gtk::make_managed<Gtk::CellRendererText>();
    auto icon_renderer = Gtk::make_managed<Gtk::CellRendererPixbuf>();
    auto col = Gtk::make_managed<Gtk::TreeViewColumn>();

    col->set_title(M("PREFERENCES_EXTERNALEDITOR_COLUMN_NAME"));
    col->set_resizable();
    col->pack_start(*icon_renderer, false);
    col->pack_start(*name_renderer);
    col->add_attribute(*icon_renderer, "gicon", model_columns.icon);
    col->add_attribute(*name_renderer, "text", model_columns.name);
    col->set_min_width(20);

    name_renderer->property_editable() = true;
    name_renderer->signal_edited().connect(
        sigc::mem_fun(*this, &ExternalEditorPreferences::setAppName));

    return col;
}

Gtk::TreeViewColumn *ExternalEditorPreferences::makeCommandColumn()
{
    auto command_renderer = Gtk::make_managed<Gtk::CellRendererText>();
    auto col = Gtk::make_managed<Gtk::TreeViewColumn>();

    col->set_title(M("PREFERENCES_EXTERNALEDITOR_COLUMN_COMMAND"));
    col->pack_start(*command_renderer);
    col->add_attribute(*command_renderer, "text", model_columns.command);

    command_renderer->property_editable() = true;
    command_renderer->signal_edited().connect(
        sigc::mem_fun(*this, &ExternalEditorPreferences::setAppCommand));

    return col;
}

void ExternalEditorPreferences::onAppChooserDialogResponse(
    int response_id, Gtk::AppChooserDialog *dialog)
{
    switch (response_id) {
        case Gtk::RESPONSE_OK:
            dialog->close();
            setApp(dialog->get_app_info());
            break;

        case Gtk::RESPONSE_CANCEL:
        case Gtk::RESPONSE_CLOSE:
            dialog->close();
            break;

        default:
            break;
    }
}

void ExternalEditorPreferences::openAppChooserDialog()
{
    if (app_chooser_dialog.get()) {
        app_chooser_dialog->refresh();
        app_chooser_dialog->show();
        return;
    }

    app_chooser_dialog.reset(new Gtk::AppChooserDialog("image/tiff"));
    app_chooser_dialog->signal_response().connect(sigc::bind(
                sigc::mem_fun(*this, &ExternalEditorPreferences::onAppChooserDialogResponse),
                app_chooser_dialog.get()
            ));
    app_chooser_dialog->set_modal();
    app_chooser_dialog->show();
}

void ExternalEditorPreferences::removeSelectedEditors()
{
    auto selection = list_view->get_selection()->get_selected_rows();

    for (const auto &selected : selection) {
        list_model->erase(list_model->get_iter(selected));
    }
}

void ExternalEditorPreferences::setApp(const Glib::RefPtr<Gio::AppInfo> app_info)
{
    auto selection = list_view->get_selection()->get_selected_rows();

    for (const auto &selected : selection) {
        auto row = *list_model->get_iter(selected);
        row[model_columns.icon] = app_info->get_icon();
        row[model_columns.name] = app_info->get_name();
        row[model_columns.command] = app_info->get_commandline();
    }
}

void ExternalEditorPreferences::setAppCommand(
    const Glib::ustring & path, const Glib::ustring & new_text)
{
    auto row_iter = list_model->get_iter(path);

    if (!row_iter->get_value(model_columns.command).compare(new_text)) {
        return;
    }

    row_iter->set_value(model_columns.command, new_text);
    row_iter->set_value(model_columns.icon, Glib::RefPtr<Gio::Icon>(nullptr));
}

void ExternalEditorPreferences::setAppName(
    const Glib::ustring & path, const Glib::ustring & new_text)
{
    list_model->get_iter(path)->set_value(model_columns.name, new_text);
}

void ExternalEditorPreferences::updateToolbarSensitivity()
{
    bool selected = list_view->get_selection()->count_selected_rows();
    button_app_chooser->set_sensitive(selected);
    button_remove->set_sensitive(selected);
}

ExternalEditorPreferences::EditorInfo::EditorInfo(
    Glib::ustring name, Glib::ustring command, Glib::ustring icon_name, void *other_data
) : name(name), icon_name(icon_name), command(command), other_data(other_data)
{
}

ExternalEditorPreferences::ModelColumns::ModelColumns()
{
    add(name);
    add(icon);
    add(command);
    add(other_data);
}
