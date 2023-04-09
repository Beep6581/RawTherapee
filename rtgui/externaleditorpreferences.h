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
#pragma once

#include <gtkmm/button.h>
#include <gtkmm/box.h>
#include <gtkmm/liststore.h>
#include <gtkmm/scrolledwindow.h>
#include <gtkmm/treemodelcolumn.h>
#include <gtkmm/treeview.h>

#include "rtappchooserdialog.h"


namespace Gtk
{

class FileChooserDialog;

}


/**
 * Widget for editing the external editors options.
 */
class ExternalEditorPreferences : public Gtk::Box
{
public:
    struct EditorTag {
        bool selected;
        EditorTag(): selected(false) {}
        explicit EditorTag(bool selected): selected(selected) {}
    };

    /**
     * Data struct containing information about an external editor.
     */
    struct EditorInfo {
        explicit EditorInfo(
            const Glib::ustring &name = Glib::ustring(),
            const Glib::ustring &command = Glib::ustring(),
            const Glib::ustring &icon_serialized = Glib::ustring(),
            bool native_command = false,
            EditorTag other_data = EditorTag()
        );
        /**
         * Name of the external editor.
         */
        Glib::ustring name;
        /**
         * The string representation of the icon. See Gio::Icon::serialize().
         */
        Glib::ustring icon_serialized;
        /**
         * The commandline for running the program. See
         * Gio::AppInfo::get_commandline()
         */
        Glib::ustring command;
        /**
         * Use the OS native launcher instead of Gio.
         */
        bool native_command;
        /**
         * Holds any other data associated with the editor. For example, it can
         * be used as a tag to uniquely identify the editor.
         */
        EditorTag other_data;
    };

    ExternalEditorPreferences();

    /**
     * Creates and returns a vector representing the external editors shown in
     * this widget.
     */
    std::vector<EditorInfo> getEditors() const;
    /**
     * Populates this widget with the external editors described in the
     * argument.
     */
    void setEditors(const std::vector<EditorInfo> &editors);

private:
    /**
     * Model representing the data fields each external editor entry has.
     */
    class ModelColumns : public Gtk::TreeModelColumnRecord
    {
    public:
        ModelColumns();
        Gtk::TreeModelColumn<Glib::ustring> name;
        Gtk::TreeModelColumn<Glib::RefPtr<Gio::Icon>> icon;
        Gtk::TreeModelColumn<Glib::ustring> command;
        Gtk::TreeModelColumn<bool> native_command;
        Gtk::TreeModelColumn<EditorTag> other_data;
    };

    ModelColumns model_columns;
    Glib::RefPtr<Gtk::ListStore> list_model; // The list of editors.
    Gtk::ScrolledWindow list_scroll_area; // Allows the list to be scrolled.
    Gtk::TreeView *list_view; // Widget for displaying the list.
    Gtk::Box toolbar; // Contains buttons for editing the list.
    Gtk::Button *button_app_chooser;
    Gtk::Button *button_add;
    Gtk::Button *button_file_chooser;
    Gtk::Button *button_remove;
    std::unique_ptr<RTAppChooserDialog> app_chooser_dialog;
    std::unique_ptr<Gtk::FileChooserDialog> file_chooser_dialog;

    /**
     * Inserts a new editor entry after the current selection, or at the end if
     * no editor is selected.
     */
    void addEditor();
    /**
     * Constructs the column for displaying the external editor name (and icon).
     */
    Gtk::TreeViewColumn *makeAppColumn();
    /**
     * Constructs the column for displaying an editable commandline.
     */
    Gtk::TreeViewColumn *makeCommandColumn();
    /**
     * Constructs the column for displaying the native command toggle.
     */
    Gtk::TreeViewColumn *makeNativeCommandColumn();
    /**
     * Called when the user is done interacting with the app chooser dialog.
     * Closes the dialog and updates the selected entry if an app was chosen.
     */
    void onAppChooserDialogResponse(int responseId, RTAppChooserDialog *dialog);
    /**
     * Called when the user is done interacting with the file chooser dialog.
     * Closes the dialog and updates the selected entry if a file was chosen.
     */
    void onFileChooserDialogResponse(int responseId, Gtk::FileChooserDialog *dialog);
    /**
     * Shows the app chooser dialog.
     */
    void openAppChooserDialog();
    /**
     * Shows the file chooser dialog for picking an executable.
     */
    void openFileChooserDialog();
    /**
     * Removes all selected editors.
     */
    void removeSelectedEditors();
    /**
     * Sets the selected entries with the provided information.
     */
    void setApp(const Glib::RefPtr<Gio::AppInfo> app_info);
    /**
     * Updates the application command and removes the icon for the given row.
     */
    void setAppCommand(const Glib::ustring & path, const Glib::ustring & new_text);
    /**
     * Updates the application name for the given row.
     */
    void setAppName(const Glib::ustring & path, const Glib::ustring & new_text);
    /**
     * Sets the sensitivity of the widgets in the toolbar to reflect the current
     * state of the list. For example, makes the remove button insensitive if no
     * entries are selected.
     */
    void updateToolbarSensitivity();
};
