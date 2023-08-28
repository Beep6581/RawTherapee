/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2012 Oliver Duis <www.oliverduis.de>
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
#include "extprog.h"

#include <cstring>
#include <iostream>

#ifdef _WIN32
#include <shlobj.h>
#include <shellapi.h>
#endif

#include <gtkmm.h>

#include "options.h"
#include "multilangmgr.h"

Glib::ustring ExtProgAction::getFullName () const
{
    return name + " [" + M(Glib::ustring::compose("EXTPROGTARGET_%1", target)) + "]";
}

bool ExtProgAction::execute (const std::vector<Glib::ustring>& fileNames) const
{
    if (fileNames.empty ()) {
        return false;
    }

    // Check if they all exists as they may not be processed yet.
    for (const auto& fileName : fileNames) {

        if (Glib::file_test (fileName, Glib::FILE_TEST_EXISTS)) {
            continue;
        }

        Gtk::MessageDialog (M("MAIN_MSG_IMAGEUNPROCESSED"), true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true).run ();
        return false;
    }

    Glib::ustring cmdLine = "\"" + filePathEXE + "\"";

    if (!preparams.empty()) {
        cmdLine += " " + preparams;
    }

    for (const auto& fileName : fileNames) {
        cmdLine += " " + Glib::shell_quote(fileName);
    }

    return ExtProgStore::spawnCommandAsync (cmdLine);
}

ExtProgStore* ExtProgStore::getInstance()
{
    static ExtProgStore instance_;
    return &instance_;
}

// Reads all profiles from the given profiles dir
void ExtProgStore::init ()
{
    MyMutex::MyLock lock(mtx);

    actions.clear ();

#ifdef _WIN32

    // Please do not add obscure little tools here, only widely used programs.
    // They should also have a proper setup program and therefore a standard path.

    searchProgram ("Photoshop", "Adobe\\Adobe Photoshop CS%1 (64 Bit)\\Photoshop.exe", "Adobe\\Adobe Photoshop CS%1\\Photoshop.exe", 9, false, true);
    searchProgram ("Photomatix Pro", "PhotomatixPro%1\\PhotomatixPro.exe", "", 9, true, true);
    searchProgram ("Paint.NET", "Paint.NET\\PaintDotNet.exe", "", 0, false, true);
    searchProgram ("MS Image Composition Editor", "Microsoft Research\\Image Composite Editor\\ICE.exe", "", 0, false, true);
    searchProgram ("PTGui", "PTGui\\PTGui.exe", "", 0, false, true);
    searchProgram ("GeoSetter", "GeoSetter\\GeoSetter.exe", "", 0, true, true);
    searchProgram ("FastStone Image Viewer", "FastStone Image Viewer\\FSViewer.exe", "", 0, true, true);
    searchProgram ("FastPictureViewer", "FastPictureViewer\\FastPictureViewer.exe", "", 0, true, true);

    if (!searchProgram ("Autopano Giga 3", "Kolor\\Autopano Giga 3.%1\\AutopanoGiga_x64.exe", "Kolor\\Autopano Giga 3.%1\\AutopanoGiga.exe", 15, true, true)) {
        if (!searchProgram ("Autopano Pro 3", "Kolor\\Autopano Pro 3.%1\\AutopanoPro_x64.exe", "Kolor\\Autopano Pro 3.%1\\AutopanoPro.exe", 15, true, true))   {
            if (!searchProgram ("Autopano Giga 2", "Kolor\\Autopano Giga 2.%1\\AutopanoGiga_x64.exe", "Kolor\\Autopano Giga 2.%1\\AutopanoGiga.exe", 6, true, true)) {
                searchProgram ("Autopano Pro 2", "Kolor\\Autopano Pro 2.%1\\AutopanoPro_x64.exe", "Kolor\\Autopano Pro 2.%1\\AutopanoPro.exe", 6, true, true);
            }
        }
    }

#endif

}
#ifdef _WIN32
bool ExtProgStore::searchProgram (const Glib::ustring& name,
                                  const Glib::ustring& exePath,
                                  const Glib::ustring& exePath86,
                                  int maxVer,
                                  bool allowRaw,
                                  bool allowQueueProcess)
{

    // get_user_special_dir crashes on some Windows configurations.
    static Glib::ustring progFilesDir, progFilesDirx86;

    if (progFilesDir.empty ()) {
        WCHAR pathW[MAX_PATH];
        char pathA[MAX_PATH];

        if (SHGetSpecialFolderPathW (NULL, pathW, CSIDL_PROGRAM_FILES, false)) {
            if (WideCharToMultiByte (CP_UTF8, 0, pathW, -1, pathA, MAX_PATH, 0, 0)) {
                progFilesDir = pathA;
            }
        }

        if (SHGetSpecialFolderPathW (NULL, pathW, CSIDL_PROGRAM_FILESX86, false)) {
            if (WideCharToMultiByte (CP_UTF8, 0, pathW, -1, pathA, MAX_PATH, 0, 0)) {
                progFilesDirx86 = pathA;
            }
        }
    }

    ExtProgAction action;
    action.name = name;
    action.target = (allowRaw ? 1 : 2);

    auto& filePath = action.filePathEXE;

    if (maxVer > 0) {

        for (auto ver = maxVer; ver >= 0; ver--) {

            filePath = progFilesDir + "\\" + Glib::ustring::compose(exePath, ver);

            if (Glib::file_test(filePath, Glib::FILE_TEST_EXISTS)) {
                break;
            }

            if (!exePath86.empty()) {

                filePath = progFilesDirx86 + "\\" + Glib::ustring::compose(exePath86, ver);

                if (Glib::file_test(filePath, Glib::FILE_TEST_EXISTS)) {
                    break;
                }
            }
            filePath.clear();
        }
    } else {

        do {

            filePath = progFilesDir + "\\" + exePath;

            if (Glib::file_test(filePath, Glib::FILE_TEST_EXISTS)) {
                break;
            }

            if (!exePath86.empty()) {

                filePath = progFilesDirx86 + "\\" + exePath86;

                if (Glib::file_test(filePath, Glib::FILE_TEST_EXISTS)) {
                    break;
                }
            }
            filePath.clear();
        } while (false);
    }

    if (!action.filePathEXE.empty ()) {

        actions.push_back (action);

        if (allowRaw && allowQueueProcess) {

            action.target = 2;
            actions.push_back (action);
        }

        return true;
    }


    return false;
}
#endif

bool ExtProgStore::spawnCommandAsync (const Glib::ustring& cmd)
{
    try {

        const auto encodedCmd = Glib::filename_from_utf8 (cmd);
        Glib::spawn_command_line_async (encodedCmd.c_str ());

        return true;

    } catch (const Glib::Exception& exception) {

        if (rtengine::settings->verbose) {
            std::cerr << "Failed to execute \"" << cmd << "\": " << exception.what() << std::endl;
        }

        return false;

    }
}

bool ExtProgStore::spawnCommandSync (const Glib::ustring& cmd)
{
    auto exitStatus = -1;

    try {

        Glib::spawn_command_line_sync (cmd, nullptr, nullptr, &exitStatus);

    } catch (const Glib::Exception& exception) {

        if (rtengine::settings->verbose) {
            std::cerr << "Failed to execute \"" << cmd << "\": " << exception.what() << std::endl;
        }

    }

    return exitStatus == 0;
}

bool ExtProgStore::openInGimp (const Glib::ustring& fileName)
{
#if defined _WIN32

    auto executable = Glib::build_filename (options.gimpDir, "bin", "gimp-win-remote");
    auto success = ShellExecute( NULL, "open", executable.c_str(), fileName.c_str(), NULL, SW_SHOWNORMAL );

#elif defined __APPLE__

    // Apps should be opened using the simplest, case-insensitive form, "open -a NameOfProgram"
    // Calling the executable directly is said to often cause trouble,
    // https://discuss.pixls.us/t/affinity-photo-as-external-editor-how-to/1756/18
    auto cmdLine = Glib::ustring("open -a GIMP \'") + fileName + Glib::ustring("\'");
    auto success = spawnCommandAsync (cmdLine);

#else

    auto cmdLine = Glib::ustring("gimp ") + Glib::shell_quote(fileName);
    auto success = spawnCommandAsync (cmdLine);

#endif

#ifdef _WIN32
    if (reinterpret_cast<uintptr_t>(success) > 32) {
        return true;
    }
#else
    if (success) {
        return true;
    }

#endif

#ifdef _WIN32

    for (auto ver = 12; ver >= 0; --ver) {

        executable = Glib::build_filename (options.gimpDir, "bin", Glib::ustring::compose (Glib::ustring("gimp-2.%1.exe"), ver));
        Glib::ustring escapedFileName = Glib::ustring::compose ("\"%1\"", fileName);
        auto lsuccess = ShellExecute( NULL, "open", executable.c_str(), escapedFileName.c_str(), NULL, SW_SHOWNORMAL );
        if (reinterpret_cast<uintptr_t>(lsuccess) > 32) {
            return true;
        }
    }

#elif defined __APPLE__

    cmdLine = Glib::ustring("open -a GIMP-dev \'") + fileName + Glib::ustring("\'");
    success = ExtProgStore::spawnCommandAsync (cmdLine);

#else

    cmdLine = Glib::ustring("gimp-remote ") + Glib::shell_quote(fileName);
    success = ExtProgStore::spawnCommandAsync (cmdLine);

#endif

    return success;
}

bool ExtProgStore::openInPhotoshop (const Glib::ustring& fileName)
{
#if defined _WIN32

    const auto executable = Glib::build_filename(options.psDir, "Photoshop.exe");
    const auto cmdLine = Glib::ustring("\"") + executable + Glib::ustring("\" \"") + fileName + Glib::ustring("\"");

#elif defined __APPLE__

    const auto cmdLine = Glib::ustring("open -a Photoshop \'") + fileName + Glib::ustring("\'");

#else

    const auto cmdLine = Glib::ustring("\"") + Glib::build_filename(options.psDir, "Photoshop.exe") + "\" " + Glib::shell_quote(fileName);

#endif

    return spawnCommandAsync (cmdLine);
}

bool ExtProgStore::openInCustomEditor (const Glib::ustring& fileName, const Glib::ustring* command)
{
    if (!command) {
        command = &(options.customEditorProg);
    }

#if defined _WIN32

    const auto cmdLine = Glib::ustring("\"") + *command + Glib::ustring("\"");
    auto success = ShellExecute( NULL, "open", cmdLine.c_str(), ('"' + fileName + '"').c_str(), NULL, SW_SHOWNORMAL );
    return (uintptr_t)success > 32;

#elif defined __APPLE__

    const auto cmdLine = *command + Glib::ustring(" \"") + fileName + Glib::ustring("\"");
    return spawnCommandAsync (cmdLine);

#else

    const auto cmdLine = *command + Glib::ustring(" ") + Glib::shell_quote(fileName);
    return spawnCommandAsync (cmdLine);

#endif

}

bool ExtProgStore::openInExternalEditor(const Glib::ustring &fileName, const EditorInfo &editorInfo)
{
    if (editorInfo.isNativeCommand) {
        if (rtengine::settings->verbose) {
            std::cout << "Launching external editor as native command." << std::endl;
        }
        const Glib::ustring command = editorInfo.commandline;
        return openInCustomEditor(fileName, &command);
    }

    if (rtengine::settings->verbose) {
        std::cout << "Launching external editor with Gio." << std::endl;
    }

    bool success = false;

    try {
        Glib::RefPtr<Gio::AppInfo> appInfo =
            Gio::AppInfo::create_from_commandline(
                editorInfo.commandline, editorInfo.name, Gio::APP_INFO_CREATE_NONE);
        success = appInfo->launch(Gio::File::create_for_path(fileName));
    } catch (const Glib::Error &e) {
        std::cerr
            << "Error launching external editor.\n"
            << "Error code #" << e.code() << ": " << e.what()
            << std::endl;
        success = false;
    }

    if (success) {
        return true;
    }

    if (rtengine::settings->verbose) {
        std::cout << "Unable to launch external editor with Gio. Trying custom launcher." << std::endl;
    }
    Glib::ustring command = editorInfo.commandline;
#if defined _WIN32
    if (command.length() > 2 && command[0] == '"' && command[command.length() - 1] == '"') {
        command = command.substr(1, command.length() - 2);
    }
#endif
    return openInCustomEditor(fileName, &command);
}
