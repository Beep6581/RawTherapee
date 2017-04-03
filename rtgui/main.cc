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
 */

#ifdef __GNUC__
#if defined(__FAST_MATH__)
#error Using the -ffast-math CFLAG is known to lead to problems. Disable it to compile RawTherapee.
#endif
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "config.h"
#include <gtkmm.h>
#include <giomm.h>
#include <iostream>
#include <tiffio.h>
#include "../rtengine/icons.h"
#include "rtwindow.h"
#include <cstring>
#include <cstdlib>
#include <locale.h>
#include "options.h"
#include "soundman.h"
#include "rtimage.h"
#include "version.h"
#include "extprog.h"
#include "dynamicprofile.h"

#ifndef WIN32
#include <glibmm/fileutils.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <glibmm/threads.h>
#else
#include <glibmm/thread.h>
#include "conio.h"
#endif

extern Options options;

// stores path to data files
Glib::ustring argv0;
Glib::ustring creditsPath;
Glib::ustring licensePath;
Glib::ustring argv1;
bool simpleEditor;
Glib::RefPtr<Gtk::CssProvider> cssForced;
Glib::RefPtr<Gtk::CssProvider> cssRT;
//Glib::Threads::Thread* mainThread;

namespace
{

// For an unknown reason, Glib::filename_to_utf8 doesn't work on reliably Windows,
// so we're using Glib::filename_to_utf8 for Linux/Apple and Glib::locale_to_utf8 for Windows.
Glib::ustring fname_to_utf8 (const char* fname)
{
#ifdef WIN32

    try {
        return Glib::locale_to_utf8 (fname);
    } catch (Glib::Error&) {
        return Glib::convert_with_fallback (fname, "UTF-8", "ISO-8859-1", "?");
    }

#else

    return Glib::filename_to_utf8 (fname);

#endif
}

}

// This recursive mutex will be used by gdk_threads_enter/leave instead of a simple mutex
static Glib::Threads::RecMutex myGdkRecMutex;

static void myGdkLockEnter()
{
    myGdkRecMutex.lock();
}
static void myGdkLockLeave()
{
    // Automatic gdk_flush for non main tread
#if AUTO_GDK_FLUSH
    //if (Glib::Thread::self() != mainThread) {
    //    gdk_flush();
    //}

#endif
    myGdkRecMutex.unlock();
}

/* Process line command options
 * Returns
 *  0 if process in batch has executed
 *  1 to start GUI (with a dir or file option)
 *  2 to start GUI because no files found
 *  -1 if there is an error in parameters
 *  -2 if an error occurred during processing
 *  -3 if at least one required procparam file was not found */
int processLineParams( int argc, char **argv );

int main(int argc, char **argv)
{
    setlocale(LC_ALL, "");
    setlocale(LC_NUMERIC, "C"); // to set decimal point to "."
    gtk_init (&argc, &argv);  // use the "--g-fatal-warnings" command line flag to make warnings fatal

    Glib::init();  // called by Gtk::Main, but this may be important for thread handling, so we call it ourselves now
    gdk_threads_set_lock_functions(G_CALLBACK(myGdkLockEnter), (G_CALLBACK(myGdkLockLeave)));
    gdk_threads_init();
    Gio::init ();

    //mainThread = Glib::Threads::Thread::self();

#ifdef BUILD_BUNDLE
    char exname[512] = {0};
    Glib::ustring exePath;
    // get the path where the rawtherapee executable is stored
#ifdef WIN32
    WCHAR exnameU[512] = {0};
    GetModuleFileNameW (NULL, exnameU, 512);
    WideCharToMultiByte(CP_UTF8, 0, exnameU, -1, exname, 512, 0, 0 );
#else

    if (readlink("/proc/self/exe", exname, 512) < 0) {
        strncpy(exname, argv[0], 512);
    }

#endif
    exePath = Glib::path_get_dirname(exname);

    // set paths
    if (Glib::path_is_absolute(DATA_SEARCH_PATH)) {
        argv0 = DATA_SEARCH_PATH;
    } else {
        argv0 = Glib::build_filename(exePath, DATA_SEARCH_PATH);
    }

    if (Glib::path_is_absolute(CREDITS_SEARCH_PATH)) {
        creditsPath = CREDITS_SEARCH_PATH;
    } else {
        creditsPath = Glib::build_filename(exePath, CREDITS_SEARCH_PATH);
    }

    if (Glib::path_is_absolute(LICENCE_SEARCH_PATH)) {
        licensePath = LICENCE_SEARCH_PATH;
    } else {
        licensePath = Glib::build_filename(exePath, LICENCE_SEARCH_PATH);
    }

#else
    argv0 = DATA_SEARCH_PATH;
    creditsPath = CREDITS_SEARCH_PATH;
    licensePath = LICENCE_SEARCH_PATH;
#endif

    if (!Options::load ()) {
        Gtk::Main m(&argc, &argv);
        Gtk::MessageDialog msgd ("Fatal error!\nThe RT_SETTINGS and/or RT_PATH environment variables are set, but use a relative path. The path must be absolute!", true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();
        return -2;
    }

    profileStore.init ();
    extProgStore->init();
    SoundManager::init();

#ifdef WIN32
    bool consoleOpened = false;

    // suppression of annoying error boxes
    SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);

    if (argc > 1 || options.rtSettings.verbose) {
        if (options.rtSettings.verbose || ( !Glib::file_test (fname_to_utf8 (argv[1]), Glib::FILE_TEST_EXISTS ) && !Glib::file_test (fname_to_utf8 (argv[1]), Glib::FILE_TEST_IS_DIR))) {
            bool stdoutRedirectedtoFile = (GetFileType(GetStdHandle(STD_OUTPUT_HANDLE)) == 0x0001);
            bool stderrRedirectedtoFile = (GetFileType(GetStdHandle(STD_ERROR_HANDLE)) == 0x0001);

            // no console, if stdout and stderr both are redirected to file
            if( !(stdoutRedirectedtoFile && stderrRedirectedtoFile)) {
                // check if parameter -w was passed.
                // We have to do that in this step, because it controls whether to open a console to show the output of following steps
                bool Console = true;

                for(int i = 1; i < argc; i++)
                    if(!strcmp(argv[i], "-w")) {
                        Console = false;
                        break;
                    }

                if(Console && AllocConsole()) {
                    AttachConsole( GetCurrentProcessId() ) ;
                    // Don't allow CTRL-C in console to terminate RT
                    SetConsoleCtrlHandler( NULL, true );
                    // Set title of console
                    char consoletitle[128];
                    sprintf(consoletitle, "RawTherapee %s Console", RTVERSION);
                    SetConsoleTitle(consoletitle);
                    // increase size of screen buffer
                    COORD c;
                    c.X = 200;
                    c.Y = 1000;
                    SetConsoleScreenBufferSize( GetStdHandle( STD_OUTPUT_HANDLE ), c );
                    // Disable console-Cursor
                    CONSOLE_CURSOR_INFO cursorInfo;
                    cursorInfo.dwSize = 100;
                    cursorInfo.bVisible = false;
                    SetConsoleCursorInfo( GetStdHandle( STD_OUTPUT_HANDLE ), &cursorInfo );

                    if(!stdoutRedirectedtoFile) {
                        freopen( "CON", "w", stdout ) ;
                    }

                    if(!stderrRedirectedtoFile) {
                        freopen( "CON", "w", stderr ) ;
                    }

                    freopen( "CON", "r", stdin ) ;

                    consoleOpened = true;

                    // printing RT's version in every case, particularly useful for the 'verbose' mode, but also for the batch processing
                    std::cout << "RawTherapee, version " << RTVERSION << std::endl;
                    std::cout << "WARNING: closing this window will close RawTherapee!" << std::endl << std::endl;
                }
            }
        }

        if(argc > 1) {
            int ret = processLineParams( argc, argv);

            if( ret <= 0 ) {
                if(consoleOpened) {
                    printf("Press any key to exit RawTherapee\n");
                    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
                    getch();
                }

                return ret;
            }
        }
    }

#else

    if (argc > 1 || options.rtSettings.verbose) {
        // printing RT's version in all case, particularly useful for the 'verbose' mode, but also for the batch processing
        std::cout << "RawTherapee, version " << RTVERSION << std::endl;
#ifdef WIN32
        std::cout << "WARNING: closing this window will close RawTherapee!" << std::endl << std::endl;
#endif

        if (argc > 1) {
            int ret = processLineParams( argc, argv);

            if( ret <= 0 ) {
                return ret;
            }
        }
    }

#endif

    if( !options.rtSettings.verbose ) {
        TIFFSetWarningHandler(nullptr);    // avoid annoying message boxes
    }

#ifndef WIN32

    // Move the old path to the new one if the new does not exist
    if (Glib::file_test(Glib::build_filename(options.rtdir, "cache"), Glib::FILE_TEST_IS_DIR) && !Glib::file_test(options.cacheBaseDir, Glib::FILE_TEST_IS_DIR)) {
        g_rename(Glib::build_filename (options.rtdir, "cache").c_str (), options.cacheBaseDir.c_str ());
    }

#endif

    simpleEditor = false;

    if( !argv1.empty() )
        if( Glib::file_test(argv1, Glib::FILE_TEST_EXISTS) && !Glib::file_test(argv1, Glib::FILE_TEST_IS_DIR)) {
            simpleEditor = true;
        }

    Gtk::Main m(&argc, &argv);

    Glib::ustring icon_path = Glib::build_filename(argv0, "images");
    Glib::RefPtr<Gtk::IconTheme> defaultIconTheme = Gtk::IconTheme::get_default();
    defaultIconTheme->append_search_path(icon_path);

    rtengine::setPaths(options);
    MyExpander::init();  // has to stay AFTER rtengine::setPaths

    // ------- loading theme files

    Glib::RefPtr<Gdk::Screen> screen = Gdk::Screen::get_default();

    if (screen) {
        Gtk::Settings::get_for_screen(screen)->property_gtk_theme_name() = "Adwaita";
        Gtk::Settings::get_for_screen(screen)->property_gtk_application_prefer_dark_theme() = true;

        Glib::RefPtr<Glib::Regex> regex = Glib::Regex::create(THEMEREGEXSTR, Glib::RegexCompileFlags::REGEX_CASELESS);
        Glib::ustring filename = Glib::build_filename(argv0, "themes", options.theme + ".css");
        if (!regex->match(options.theme + ".css") || !Glib::file_test(filename, Glib::FILE_TEST_EXISTS)) {
            options.theme = "RawTherapee-GTK";
            // We're not testing GTK_MAJOR_VERSION == 3 here, since this branch requires Gtk3 only
            if (GTK_MINOR_VERSION < 20) {
                options.theme = options.theme + "3-_19";
            } else {
                options.theme = options.theme + "3-20_";
            }
            filename = Glib::build_filename(argv0, "themes", options.theme + ".css");
        }
        cssRT = Gtk::CssProvider::create();

        try {
            cssRT->load_from_path (filename);
            Gtk::StyleContext::add_provider_for_screen(screen, cssRT, GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
        } catch (Glib::Error &err) {
            printf("Error: Can't load css file \"%s\"\nMessage: %s\n", filename.c_str(), err.what().c_str());
        } catch (...) {
            printf("Error: Can't load css file \"%s\"\n", filename.c_str());
        }

        // Set the font face and size
        if (options.fontFamily != "default") {
            try {
                cssForced = Gtk::CssProvider::create();
                //GTK318
                #if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
                cssForced->load_from_data (Glib::ustring::compose("* { font-family: %1; font-size: %2px }", options.fontFamily, options.fontSize));
                #else
                cssForced->load_from_data (Glib::ustring::compose("* { font-family: %1; font-size: %2pt }", options.fontFamily, options.fontSize));
                #endif
                //GTK318
                Gtk::StyleContext::add_provider_for_screen(screen, cssForced, GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
            } catch (Glib::Error &err) {
                printf("Error: \"%s\"\n", err.what().c_str());
            } catch (...) {
                printf("Error: Can't find the font named \"%s\"\n", options.fontFamily.c_str());
            }
        }
    }

#ifndef NDEBUG
    else if (!screen) {
        printf("ERROR: Can't get default screen!\n");
    }

#endif

    // ------- end loading theme files

    gdk_threads_enter ();
    RTWindow *rtWindow = new class RTWindow();

    // alerting users if the default raw and image profiles are missing
    if (options.is_defProfRawMissing()) {
        Gtk::MessageDialog msgd (Glib::ustring::compose(M("OPTIONS_DEFRAW_MISSING"), options.defProfRaw), true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();
    }

    if (options.is_defProfImgMissing()) {
        Gtk::MessageDialog msgd (Glib::ustring::compose(M("OPTIONS_DEFIMG_MISSING"), options.defProfImg), true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();
    }

    // opening the main window
    m.run(*rtWindow);

    gdk_threads_leave ();
    delete rtWindow;
    rtengine::cleanup();

#ifdef WIN32

    if (consoleOpened) {
        printf("Press any key to exit RawTherapee\n");
        FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
        getch();
    }

#endif

    return 0;
}

void deleteProcParams(std::vector<rtengine::procparams::PartialProfile*> &pparams)
{
    for (unsigned int i = 0; i < pparams.size(); i++) {
        pparams[i]->deleteInstance();
        delete pparams[i];
        pparams[i] = NULL;
    }

    return;
}

int processLineParams( int argc, char **argv )
{
    unsigned errors = 0;

    for( int iArg = 1; iArg < argc; iArg++) {
        if( argv[iArg][0] == '-' ) {
            switch( argv[iArg][1] ) {
#ifdef WIN32

            case 'w': // This case is handled outside this function
                break;
#endif

            case 'h':
            case '?':
            default: {
                Glib::ustring pparamsExt = paramFileExtension.substr(1);
                std::cout << "  An advanced, cross-platform program for developing raw photos." << std::endl;
                std::cout << std::endl;
                std::cout << "  Website: http://www.rawtherapee.com/" << std::endl;
                std::cout << "  Documentation: http://rawpedia.rawtherapee.com/" << std::endl;
                std::cout << "  Forum: https://discuss.pixls.us/c/software/rawtherapee" << std::endl;
                std::cout << "  Code and bug reports: https://github.com/Beep6581/RawTherapee" << std::endl;
                std::cout << std::endl;
                std::cout << "Symbols:" << std::endl;
                std::cout << "  <Chevrons> indicate parameters you can change." << std::endl;
                std::cout << "  [Square brackets] mean the parameter is optional." << std::endl;
                std::cout << "  The pipe symbol | indicates a choice of one or the other." << std::endl;
                std::cout << "  The dash symbol - denotes a range of possible values from one to the other." << std::endl;
                std::cout << std::endl;
                std::cout << "Usage:" << std::endl;
                std::cout << "  " << Glib::path_get_basename(argv[0]) << " <folder>           Start File Browser inside folder." << std::endl;
                std::cout << "  " << Glib::path_get_basename(argv[0]) << " <file>             Start Image Editor with file." << std::endl;
                std::cout << std::endl;
                std::cout << "Options:" << std::endl;
#ifdef WIN32
                std::cout << "  -w Do not open the Windows console" << std::endl;
#endif
                std::cout << "  -h -? Display this help message" << std::endl;
                return -1;
            }
            }
        } else {
            argv1 = fname_to_utf8 (argv[iArg]);
            break;
        }
    }

    if( !argv1.empty() ) {
        return 1;
    }

    return errors > 0 ? -2 : 0;
}
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif // __GNUC__