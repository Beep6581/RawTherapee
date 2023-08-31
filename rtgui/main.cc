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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef __GNUC__
#if defined(__FAST_MATH__)
#error Using the -ffast-math CFLAG is known to lead to problems. Disable it to compile RawTherapee.
#endif
#endif

#include "config.h"
#include <gtkmm.h>
#include <giomm.h>
#include <iostream>
#include <tiffio.h>
#include "rtwindow.h"
#include <cstring>
#include <cstdlib>
#include <locale.h>
#include <lensfun.h>
#include "cachemanager.h"
#include "editorpanel.h"
#include "filecatalog.h"
#include "filepanel.h"
#include "options.h"
#include "soundman.h"
#include "rtimage.h"
#include "version.h"
#include "extprog.h"
#include "../rtengine/dynamicprofile.h"
#include "../rtengine/procparams.h"
#include "pathutils.h"

#ifndef _WIN32
#include <glibmm/fileutils.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <glibmm/threads.h>
#else
#include <glibmm/thread.h>
#include "conio.h"
#include "windows.h"
#endif

// Set this to 1 to make RT work when started with Eclipse and arguments, at least on Windows platform
#define ECLIPSE_ARGS 0

// stores path to data files
Glib::ustring argv0;
Glib::ustring creditsPath;
Glib::ustring licensePath;
Glib::ustring argv1;
Glib::ustring argv2;
bool simpleEditor = false;
bool gimpPlugin = false;
bool remote = false;
unsigned char initialGdkScale = 1;
//Glib::Threads::Thread* mainThread;

namespace {

// This recursive mutex will be used by gdk_threads_enter/leave instead of a simple mutex
static Glib::Threads::RecMutex myGdkRecMutex;

static void myGdkLockEnter()
{
    myGdkRecMutex.lock();
}
static void myGdkLockLeave()
{
    // Automatic gdk_flush for non main thread
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
//int processLineParams ( int argc, char **argv );
int processLineParams ( int argc, char **argv )
{
    int ret = 1;
    for ( int iArg = 1; iArg < argc; iArg++) {
        Glib::ustring currParam (argv[iArg]);
        if ( currParam.empty() ) {
            continue;
        }
#if ECLIPSE_ARGS
        currParam = currParam.substr (1, currParam.length() - 2);
#endif

        if ( currParam.at (0) == '-' && currParam.size() > 1 ) {
            switch ( currParam.at (1) ) {
                case '-':
                    // GTK --argument, we're skipping it
                    break;

#ifdef _WIN32

                case 'w': // This case is handled outside this function
                    break;
#endif

                case 'v':
                    printf("RawTherapee, version %s\n", RTVERSION);
                    ret = 0;
                    break;

#ifndef __APPLE__ // TODO agriggio - there seems to be already some "single instance app" support for OSX in rtwindow. Disabling it here until I understand how to merge the two

                case 'R':
                    if (!gimpPlugin) {
                        remote = true;
                    }

                    break;
#endif

                case 'g':
                    if (currParam == "-gimp") {
                        gimpPlugin = true;
                        simpleEditor = true;
                        remote = false;
                        break;
                    }

                // no break here on purpose

                case 'h':
                case '?':
                default: {
                    printf("  An advanced, cross-platform program for developing raw photos.\n\n");
                    printf("  Website: http://www.rawtherapee.com/\n");
                    printf("  Documentation: http://rawpedia.rawtherapee.com/\n");
                    printf("  Forum: https://discuss.pixls.us/c/software/rawtherapee\n");
                    printf("  Code and bug reports: https://github.com/Beep6581/RawTherapee\n\n");
                    printf("Symbols:\n");
                    printf("  <Chevrons> indicate parameters you can change.\n\n");
                    printf("Usage:\n");
                    printf("  %s <folder>           Start File Browser inside folder.\n",Glib::path_get_basename (argv[0]).c_str());
                    printf("  %s <file>             Start Image Editor with file.\n\n",Glib::path_get_basename (argv[0]).c_str());
                    std::cout << std::endl;
                    printf("Options:\n");
#ifdef _WIN32
                    printf("  -w Do not open the Windows console\n");
#endif
                    printf("  -v Print RawTherapee version number and exit\n");
#ifndef __APPLE__
                    printf("  -R Raise an already running RawTherapee instance (if available)\n");
#endif
                    printf("  -h -? Display this help message\n");

                    ret = -1;
                    break;
                }
            }
        } else {
            if (argv1.empty()) {
                argv1 = Glib::ustring (fname_to_utf8 (argv[iArg]));
#if ECLIPSE_ARGS
                argv1 = argv1.substr (1, argv1.length() - 2);
#endif
            } else if (gimpPlugin) {
                argv2 = Glib::ustring (fname_to_utf8 (argv[iArg]));
                break;
            }

            if (!gimpPlugin) {
                break;
            }
        }
    }

    return ret;
}


bool init_rt()
{
    extProgStore->init();
    SoundManager::init();

    if (!rtengine::settings->verbose) {
        TIFFSetWarningHandler (nullptr);   // avoid annoying message boxes
    }

#ifndef _WIN32

    // Move the old path to the new one if the new does not exist
    if (Glib::file_test (Glib::build_filename (options.rtdir, "cache"), Glib::FILE_TEST_IS_DIR) && !Glib::file_test (options.cacheBaseDir, Glib::FILE_TEST_IS_DIR)) {
        g_rename (Glib::build_filename (options.rtdir, "cache").c_str (), options.cacheBaseDir.c_str ());
    }

#endif

    return true;
}


void cleanup_rt()
{
    rtengine::cleanup();
}


RTWindow *create_rt_window()
{
    Glib::ustring icon_path = Glib::build_filename (argv0, "images");
    Glib::RefPtr<Gtk::IconTheme> defaultIconTheme = Gtk::IconTheme::get_default();
    defaultIconTheme->append_search_path (icon_path);

    //gdk_threads_enter ();
    RTWindow *rtWindow = new RTWindow();
    rtWindow->setWindowSize(); // Need to be called after RTWindow creation to work with all OS Windows Manager
    return rtWindow;
}


class RTApplication: public Gtk::Application
{
public:
    RTApplication():
        Gtk::Application ("com.rawtherapee.application",
                          Gio::APPLICATION_HANDLES_OPEN),
        rtWindow (nullptr)
    {
    }

    ~RTApplication() override
    {
        if (rtWindow) {
            delete rtWindow;
        }

        cleanup_rt();
    }

private:
    bool create_window()
    {
        if (rtWindow) {
            return true;
        }

        if (!init_rt()) {
            Gtk::MessageDialog msgd ("Fatal error!\nThe RT_SETTINGS and/or RT_PATH environment variables are set, but use a relative path. The path must be absolute!", true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
            add_window (msgd);
            msgd.run ();
            return false;
        } else {
            rtWindow = create_rt_window();
            add_window (*rtWindow);
            return true;
        }
    }

    // Override default signal handlers:
    void on_activate() override
    {
        if (create_window()) {
            rtWindow->present();
        }
    }

    void on_open (const Gio::Application::type_vec_files& files,
                  const Glib::ustring& hint) override
    {
        if (create_window()) {
            struct Data {
                std::vector<Thumbnail *> entries;
                Glib::ustring lastfilename;
                FileCatalog *filecatalog;
            };
            Data *d = new Data;
            d->filecatalog = rtWindow->fpanel->fileCatalog;

            for (const auto &f : files) {
                Thumbnail *thm = cacheMgr->getEntry (f->get_path());

                if (thm) {
                    d->entries.push_back (thm);
                    d->lastfilename = f->get_path();
                }
            }

            if (!d->entries.empty()) {
                const auto doit =
                [] (gpointer data) -> gboolean {
                    Data *d = static_cast<Data *> (data);
                    d->filecatalog->openRequested (d->entries);
                    d->filecatalog->selectImage (d->lastfilename, true);
                    delete d;
                    return FALSE;
                };
                gdk_threads_add_idle (doit, d);
            } else {
                delete d;
            }

            rtWindow->present();
        }
    }

private:
    RTWindow *rtWindow;
};

void show_gimp_plugin_info_dialog(Gtk::Window *parent)
{
    if (options.gimpPluginShowInfoDialog) {
        Gtk::MessageDialog info(*parent, M("GIMP_PLUGIN_INFO"), false, Gtk::MESSAGE_INFO, Gtk::BUTTONS_OK, true);
        Gtk::Box *box = info.get_message_area();
        Gtk::CheckButton dontshowagain(M("DONT_SHOW_AGAIN"));
        dontshowagain.show();
        box->pack_start(dontshowagain);
        info.run();
        options.gimpPluginShowInfoDialog = !dontshowagain.get_active();
    }
}

} // namespace


int main (int argc, char **argv)
{
    setlocale (LC_ALL, "");
    setlocale (LC_NUMERIC, "C"); // to set decimal point to "."

    simpleEditor = false;
    gimpPlugin = false;
    remote = false;
    argv0 = "";
    argv1 = "";
    argv2 = "";

    Glib::init();  // called by Gtk::Main, but this may be important for thread handling, so we call it ourselves now
    Gio::init ();

#ifdef _WIN32
    if (GetFileType (GetStdHandle (STD_OUTPUT_HANDLE)) == 0x0003) {
        // started from msys2 console => do not buffer stdout
        setbuf(stdout, NULL);
    }
#endif

#ifdef BUILD_BUNDLE
    char exname[512] = {0};
    Glib::ustring exePath;
    // get the path where the rawtherapee executable is stored
#ifdef _WIN32
    WCHAR exnameU[512] = {0};
    GetModuleFileNameW (NULL, exnameU, 511);
    WideCharToMultiByte (CP_UTF8, 0, exnameU, -1, exname, 511, 0, 0 );
#else

    if (readlink ("/proc/self/exe", exname, 511) < 0) {
        strncpy (exname, argv[0], 511);
    }

#endif
    exePath = Glib::path_get_dirname (exname);

    // set paths
    if (Glib::path_is_absolute (DATA_SEARCH_PATH)) {
        argv0 = DATA_SEARCH_PATH;
    } else {
        argv0 = Glib::build_filename (exePath, DATA_SEARCH_PATH);
    }

    if (Glib::path_is_absolute (CREDITS_SEARCH_PATH)) {
        creditsPath = CREDITS_SEARCH_PATH;
    } else {
        creditsPath = Glib::build_filename (exePath, CREDITS_SEARCH_PATH);
    }

    if (Glib::path_is_absolute (LICENCE_SEARCH_PATH)) {
        licensePath = LICENCE_SEARCH_PATH;
    } else {
        licensePath = Glib::build_filename (exePath, LICENCE_SEARCH_PATH);
    }

    options.rtSettings.lensfunDbDirectory = LENSFUN_DB_PATH;
    options.rtSettings.lensfunDbBundleDirectory = LENSFUN_DB_PATH;

#else
    argv0 = DATA_SEARCH_PATH;
    creditsPath = CREDITS_SEARCH_PATH;
    licensePath = LICENCE_SEARCH_PATH;
    options.rtSettings.lensfunDbDirectory = LENSFUN_DB_PATH;
    options.rtSettings.lensfunDbBundleDirectory = LENSFUN_DB_PATH;
#endif

#ifdef _WIN32
    bool consoleOpened = false;

    // suppression of annoying error boxes
    SetErrorMode (SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);

    if (argc > 1) {
        if (!remote && !Glib::file_test (argv1, Glib::FILE_TEST_EXISTS ) && !Glib::file_test (argv1, Glib::FILE_TEST_IS_DIR)) {
            const bool stdoutRedirecttoConsole = (GetFileType (GetStdHandle (STD_OUTPUT_HANDLE)) == 0x0000);
            // open console, if stdout is invalid
            if (stdoutRedirecttoConsole) {
                // check if parameter -w was passed.
                // We have to do that in this step, because it controls whether to open a console to show the output of following steps
                bool Console = true;

                for (int i = 1; i < argc; i++)
                    if (!strcmp (argv[i], "-w") || !strcmp (argv[i], "-R") || !strcmp (argv[i], "-gimp")) {
                        Console = false;
                        break;
                    }

                if (Console && AllocConsole()) {
                    AttachConsole ( GetCurrentProcessId() ) ;
                    // Don't allow CTRL-C in console to terminate RT
                    SetConsoleCtrlHandler ( NULL, true );
                    // Set title of console
                    char consoletitle[128];
                    snprintf(consoletitle, sizeof(consoletitle), "RawTherapee %s Console", RTVERSION);
                    SetConsoleTitle (consoletitle);
                    // increase size of screen buffer
                    COORD c;
                    c.X = 200;
                    c.Y = 1000;
                    SetConsoleScreenBufferSize ( GetStdHandle ( STD_OUTPUT_HANDLE ), c );
                    // Disable console-Cursor
                    CONSOLE_CURSOR_INFO cursorInfo;
                    cursorInfo.dwSize = 100;
                    cursorInfo.bVisible = false;
                    SetConsoleCursorInfo ( GetStdHandle ( STD_OUTPUT_HANDLE ), &cursorInfo );

                    // we also redirect stderr to console
                    freopen ( "CON", "w", stdout ) ;
                    freopen ( "CON", "w", stderr ) ;

                    freopen ( "CON", "r", stdin ) ;

                    consoleOpened = true;
                }
            }
        }
        int ret = processLineParams ( argc, argv);

        if ( ret <= 0 ) {
            fflush(stdout);
            if (consoleOpened) {
                printf ("Press any key to exit RawTherapee\n");
                FlushConsoleInputBuffer (GetStdHandle (STD_INPUT_HANDLE));
                getch();
            }

            return ret;
        }
    }

#else

    if (argc > 1) {
        int ret = processLineParams ( argc, argv);

        if ( ret <= 0 ) {
            return ret;
        }
    }

#endif

    Glib::ustring fatalError;

    try {
        Options::load();
    } catch (Options::Error &e) {
        fatalError = e.get_msg();
    }

    if (gimpPlugin) {
        if (!Glib::file_test (argv1, Glib::FILE_TEST_EXISTS) || Glib::file_test (argv1, Glib::FILE_TEST_IS_DIR)) {
            printf ("Error: argv1 doesn't exist\n");
            return 1;
        }

        if (argv2.empty()) {
            printf ("Error: -gimp requires two arguments\n");
            return 1;
        }
    } else if (!remote && Glib::file_test(argv1, Glib::FILE_TEST_EXISTS) && !Glib::file_test(argv1, Glib::FILE_TEST_IS_DIR)) {
        simpleEditor = true;
    }

    int ret = 0;

    if (options.pseudoHiDPISupport) {
        // Reading/updating GDK_SCALE early if it exists
        const gchar *gscale = g_getenv("GDK_SCALE");
        if (gscale && gscale[0] == '2') {
            initialGdkScale = 2;
        }
        // HOMBRE: On Windows, if resolution is set to 200%, Gtk internal variables are SCALE=2 and DPI=96
        g_setenv("GDK_SCALE", "1", true);
    }

    gdk_threads_set_lock_functions (G_CALLBACK (myGdkLockEnter), (G_CALLBACK (myGdkLockLeave)));
    gdk_threads_init();
    gtk_init (&argc, &argv);  // use the "--g-fatal-warnings" command line flag to make warnings fatal

    if (fatalError.empty() && remote) {
        char *app_argv[2] = { const_cast<char *> (argv0.c_str()) };
        int app_argc = 1;

        if (!argv1.empty()) {
            app_argc = 2;
            app_argv[1] = const_cast<char *> (argv1.c_str());
        }

        RTApplication app;
        ret = app.run (app_argc, app_argv);
    } else {
        if (fatalError.empty() && init_rt()) {
            Gtk::Main m (&argc, &argv);
            gdk_threads_enter();
            const std::unique_ptr<RTWindow> rtWindow (create_rt_window());
            if (gimpPlugin) {
                show_gimp_plugin_info_dialog(rtWindow.get());
            }
            m.run (*rtWindow);
            gdk_threads_leave();

            if (gimpPlugin && rtWindow->epanel && rtWindow->epanel->isRealized()) {
                if (!rtWindow->epanel->saveImmediately(argv2, SaveFormat())) {
                    ret = -2;
                }
            }

            cleanup_rt();
        } else {
            Gtk::Main m (&argc, &argv);
            Gtk::MessageDialog msgd (Glib::ustring::compose("FATAL ERROR!\n\n%1", fatalError), true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
            msgd.run ();
            ret = -2;
        }
    }

#ifdef _WIN32

    if (consoleOpened) {
        printf ("Press any key to exit RawTherapee\n");
        fflush(stdout);
        FlushConsoleInputBuffer (GetStdHandle (STD_INPUT_HANDLE));
        getch();
    }

#endif

    return ret;
}
