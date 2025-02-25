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
#include <iostream>

#include <giomm/contenttype.h>
#include <glibmm/shell.h>
#include <gtkmm/filechooserdialog.h>
#include <gtkmm/stock.h>

#include "externaleditorpreferences.h"
#include "multilangmgr.h"


ExternalEditorPreferences::ExternalEditorPreferences():
    Box(Gtk::Orientation::ORIENTATION_VERTICAL),
    list_model(Gtk::ListStore::create(model_columns)),
    toolbar(Gtk::Orientation::ORIENTATION_HORIZONTAL)
{
    // List view.
    list_view = Gtk::manage(new Gtk::TreeView());
    list_view->set_model(list_model);
    list_view->append_column(*Gtk::manage(makeAppColumn()));
#ifndef __APPLE__
    list_view->append_column(*Gtk::manage(makeNativeCommandColumn()));
#endif
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
    button_add = Gtk::manage(new Gtk::Button());
    button_remove = Gtk::manage(new Gtk::Button());
    button_add->set_image_from_icon_name("add-small");
    button_remove->set_image_from_icon_name("remove-small");
    button_app_chooser =
#ifdef __APPLE__
        nullptr;
#else
        Gtk::manage(new Gtk::Button(M("PREFERENCES_EXTERNALEDITOR_CHANGE")));
#endif
    button_file_chooser = Gtk::manage(new Gtk::Button(M("PREFERENCES_EXTERNALEDITOR_CHANGE_FILE")));

    if (button_app_chooser) {
        button_app_chooser->signal_pressed().connect(sigc::mem_fun(
                    *this, &ExternalEditorPreferences::openAppChooserDialog));
    }
    button_add->signal_pressed().connect(sigc::mem_fun(
            *this, &ExternalEditorPreferences::addEditor));
    button_file_chooser->signal_pressed().connect(sigc::mem_fun(
        *this, &ExternalEditorPreferences::openFileChooserDialog));
    button_remove->signal_pressed().connect(sigc::mem_fun(
            *this, &ExternalEditorPreferences::removeSelectedEditors));

    list_view->get_selection()->signal_changed().connect(sigc::mem_fun(
                *this, &ExternalEditorPreferences::updateToolbarSensitivity));
    updateToolbarSensitivity();

    // Toolbar.
    toolbar.set_halign(Gtk::Align::ALIGN_END);
    if (button_app_chooser) {
        toolbar.add(*button_app_chooser);
    }
    toolbar.add(*button_file_chooser);
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
    std::vector<EditorInfo> editors;

    auto children = list_model->children();

    for (auto rowIter = children.begin(); rowIter != children.end(); rowIter++) {
        const auto icon = rowIter->get_value(model_columns.icon);
        const auto &icon_serialized = !icon ? "" : icon->serialize().print();
        editors.emplace_back(
            rowIter->get_value(model_columns.name),
            rowIter->get_value(model_columns.command),
            icon_serialized,
            rowIter->get_value(model_columns.native_command),
            rowIter->get_value(model_columns.other_data));
    }

    return editors;
}

void ExternalEditorPreferences::setEditors(
    const std::vector<ExternalEditorPreferences::EditorInfo> &editors)
{
    list_model->clear();

    for (const EditorInfo & editor : editors) {
        auto row = *list_model->append();
        Glib::RefPtr<Gio::Icon> icon;

        // Get icon.
        if (editor.icon_serialized.empty()) {
            icon = Glib::RefPtr<Gio::Icon>();
        } else {
            GError *e = nullptr;
            GVariant *icon_variant = g_variant_parse(
                nullptr, editor.icon_serialized.c_str(), nullptr, nullptr, &e);
            if (e) {
                std::cerr
                    << "Error loading external editor icon from \""
                    << editor.icon_serialized << "\": " << e->message
                    << std::endl;
                icon = Glib::RefPtr<Gio::Icon>();
            } else {
                icon = Gio::Icon::deserialize(Glib::VariantBase(icon_variant));
            }
        }

        row[model_columns.name] = editor.name;
        row[model_columns.icon] = icon;
        row[model_columns.command] = editor.command;
        row[model_columns.native_command] = editor.native_command;
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
#ifdef __APPLE__
    row[model_columns.native_command] = true;
#endif
    list_view->get_selection()->select(row);
}

Gtk::TreeViewColumn *ExternalEditorPreferences::makeAppColumn()
{
    auto name_renderer = Gtk::manage(new Gtk::CellRendererText());
    auto icon_renderer = Gtk::manage(new Gtk::CellRendererPixbuf());
    auto col = Gtk::manage(new Gtk::TreeViewColumn());

    col->set_title(M("PREFERENCES_EXTERNALEDITOR_COLUMN_NAME"));
    col->set_resizable();
    col->pack_start(*icon_renderer, false);
    col->pack_start(*name_renderer);
    col->add_attribute(icon_renderer->property_gicon(), model_columns.icon);
    col->add_attribute(name_renderer->property_text(), model_columns.name);
    col->set_min_width(20);

    name_renderer->property_editable() = true;
    name_renderer->signal_edited().connect(
        sigc::mem_fun(*this, &ExternalEditorPreferences::setAppName));

    return col;
}

Gtk::TreeViewColumn *ExternalEditorPreferences::makeCommandColumn()
{
    auto command_renderer = Gtk::manage(new Gtk::CellRendererText());
    auto col = Gtk::manage(new Gtk::TreeViewColumn());

    col->set_title(M("PREFERENCES_EXTERNALEDITOR_COLUMN_COMMAND"));
    col->pack_start(*command_renderer);
    col->add_attribute(command_renderer->property_text(), model_columns.command);

    command_renderer->property_editable() = true;
    command_renderer->signal_edited().connect(
        sigc::mem_fun(*this, &ExternalEditorPreferences::setAppCommand));

    return col;
}

Gtk::TreeViewColumn *ExternalEditorPreferences::makeNativeCommandColumn()
{
    auto toggle_renderer = Gtk::manage(new Gtk::CellRendererToggle());
    auto col = Gtk::manage(new Gtk::TreeViewColumn());

    col->set_title(M("PREFERENCES_EXTERNALEDITOR_COLUMN_NATIVE_COMMAND"));
    col->pack_start(*toggle_renderer);
    col->add_attribute(toggle_renderer->property_active(), model_columns.native_command);

    toggle_renderer->signal_toggled().connect([this](const Glib::ustring &path) {
        const auto row_iter = list_model->get_iter(path);
        bool new_value = !row_iter->get_value(model_columns.native_command);
        row_iter->set_value(model_columns.native_command, new_value);
    });

    return col;
}

void ExternalEditorPreferences::onAppChooserDialogResponse(
    int response_id, RTAppChooserDialog *dialog)
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

void ExternalEditorPreferences::onFileChooserDialogResponse(
        int response_id, Gtk::FileChooserDialog *dialog)
{
    switch (response_id) {
        case Gtk::RESPONSE_OK: {
            dialog->close();

            auto selection = list_view->get_selection()->get_selected_rows();
            for (const auto &selected : selection) {
                auto row = *list_model->get_iter(selected);
                row[model_columns.icon] = Glib::RefPtr<Gio::Icon>(nullptr);
                row[model_columns.native_command] =
#ifdef __APPLE__
                    true;
#else
                    false;
#endif
                row[model_columns.command] =
#ifdef _WIN32
                    '"' + dialog->get_filename() + '"';
#else
                    Glib::shell_quote(dialog->get_filename());
#endif
            }

            break;
        }

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

    app_chooser_dialog.reset(new RTAppChooserDialog("image/tiff"));
    app_chooser_dialog->signal_response().connect(sigc::bind(
                sigc::mem_fun(*this, &ExternalEditorPreferences::onAppChooserDialogResponse),
                app_chooser_dialog.get()
            ));
    app_chooser_dialog->set_modal();
    app_chooser_dialog->show();
}

void ExternalEditorPreferences::openFileChooserDialog()
{
    if (file_chooser_dialog.get()) {
        file_chooser_dialog->show();
        return;
    }

    file_chooser_dialog.reset(new Gtk::FileChooserDialog(M("PREFERENCES_EXTERNALEDITOR_CHANGE_FILE")));

    const auto exe_filter = Gtk::FileFilter::create();
    exe_filter->set_name(M("FILECHOOSER_FILTER_EXECUTABLE"));
    exe_filter->add_custom(Gtk::FILE_FILTER_MIME_TYPE, [](const Gtk::FileFilter::Info &info) {
#ifdef _WIN32
        return info.mime_type == "application/x-msdownload";
#else
        return Gio::content_type_can_be_executable(info.mime_type);
#endif
    });
    const auto all_filter = Gtk::FileFilter::create();
    all_filter->set_name(M("FILECHOOSER_FILTER_ANY"));
    all_filter->add_pattern("*");
    file_chooser_dialog->add_filter(exe_filter);
    file_chooser_dialog->add_filter(all_filter);

    file_chooser_dialog->signal_response().connect(sigc::bind(
        sigc::mem_fun(*this, &ExternalEditorPreferences::onFileChooserDialogResponse),
        file_chooser_dialog.get()));
    file_chooser_dialog->set_modal();
    file_chooser_dialog->add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    file_chooser_dialog->add_button(M("GENERAL_OPEN"), Gtk::RESPONSE_OK);
    file_chooser_dialog->show();
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
        row[model_columns.native_command] = false;
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
    if (button_app_chooser) {
        button_app_chooser->set_sensitive(selected);
    }
    button_file_chooser->set_sensitive(selected);
    button_remove->set_sensitive(selected);
}

ExternalEditorPreferences::EditorInfo::EditorInfo(
    const Glib::ustring &name,
    const Glib::ustring &command,
    const Glib::ustring &icon_serialized,
    bool native_command,
    EditorTag other_data) :
    name(name),
    icon_serialized(icon_serialized),
    command(command),
    native_command(native_command),
    other_data(other_data)
{
}

ExternalEditorPreferences::ModelColumns::ModelColumns()
{
    add(name);
    add(icon);
    add(command);
    add(native_command);
    add(other_data);
}
