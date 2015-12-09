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
// generated 2004/6/3 19:15:32 CEST by gabor@darkstar.(none)
// using glademm V2.5.0
//
// newer (non customized) versions of this file go to raw.cc_new

// This file is for your program, I won't touch it again!

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
#include "options.h"
#include "soundman.h"
#include "rtimage.h"
#include "version.h"
#include "extprog.h"

#ifndef WIN32
#include <glibmm/fileutils.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <glibmm/threads.h>
#else
#include <glibmm/thread.h>
#include "conio.h"
#endif

#include "../rtengine/safegtk.h"

extern Options options;

// stores path to data files
Glib::ustring argv0;
Glib::ustring creditsPath;
Glib::ustring licensePath;
Glib::ustring argv1;
bool simpleEditor;
Glib::Thread* mainThread;


// This recursive mutex will be used by gdk_threads_enter/leave instead of a simple mutex
#ifdef WIN32
static Glib::RecMutex myGdkRecMutex;
#else
static Glib::Threads::RecMutex myGdkRecMutex;
#endif

static void myGdkLockEnter()
{
    myGdkRecMutex.lock();
}
static void myGdkLockLeave()
{
    // Automatic gdk_flush for non main tread
#if AUTO_GDK_FLUSH
    if (Glib::Thread::self() != mainThread) {
        gdk_flush();
    }

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
    // Uncomment the following line if you want to use the "--g-fatal-warnings" command line flag
    //gtk_init (&argc, &argv);

    Glib::thread_init();
    gdk_threads_set_lock_functions(G_CALLBACK(myGdkLockEnter), (G_CALLBACK(myGdkLockLeave)));
    gdk_threads_init();
    Gio::init ();

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

    mainThread = Glib::Thread::self();

    if (!Options::load ()) {
        Gtk::Main m(&argc, &argv);
        Gtk::MessageDialog msgd ("Fatal error!\nThe RT_SETTINGS and/or RT_PATH environment variables are set, but use a relative path. The path must be absolute!", true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();
        return -2;
    }

    extProgStore->init();
    SoundManager::init();

#ifdef WIN32
    bool consoleOpened = false;

    if (argc > 1 || options.rtSettings.verbose) {
        if(options.rtSettings.verbose || ( !safe_file_test( safe_filename_to_utf8(argv[1]), Glib::FILE_TEST_EXISTS ) && !safe_file_test( safe_filename_to_utf8(argv[1]), Glib::FILE_TEST_IS_DIR ))) {
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

                if(Console) {
                    AllocConsole();
                    AttachConsole( GetCurrentProcessId() ) ;
                    // Don't allow CTRL-C in console to terminate RT
                    SetConsoleCtrlHandler( NULL, true );
                    // Set title of console
                    char consoletitle[128];
                    sprintf(consoletitle, "RawTherapee %s Console", VERSION);
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
                    std::cout << "RawTherapee, version " << VERSION << std::endl;
                    std::cout << "WARNING: closing this window will close RawTherapee!" << std::endl << std::endl;
                }
            }
        }

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

#else

    if (argc > 1 || options.rtSettings.verbose) {
        // printing RT's version in all case, particularly useful for the 'verbose' mode, but also for the batch processing
        std::cout << "RawTherapee, version " << VERSION << std::endl;
        std::cout << "WARNING: closing this window will close RawTherapee!" << std::endl << std::endl;

        if (argc > 1) {
            int ret = processLineParams( argc, argv);

            if( ret <= 0 ) {
                return ret;
            }
        }
    }

#endif

    if( !options.rtSettings.verbose ) {
        TIFFSetWarningHandler(NULL);    // avoid annoying message boxes
    }

#ifndef WIN32

    // Move the old path to the new one if the new does not exist
    if (safe_file_test(Glib::build_filename(options.rtdir, "cache"), Glib::FILE_TEST_IS_DIR) && !safe_file_test(options.cacheBaseDir, Glib::FILE_TEST_IS_DIR)) {
        safe_g_rename(Glib::build_filename(options.rtdir, "cache"), options.cacheBaseDir);
    }

#endif

    simpleEditor = false;

    if( !argv1.empty() )
        if( safe_file_test(argv1, Glib::FILE_TEST_EXISTS) && !safe_file_test(argv1, Glib::FILE_TEST_IS_DIR)) {
            simpleEditor = true;
        }

    if (options.theme.empty()) {
        options.theme = "21-Gray-Gray";
    } else {
        std::string themeFile = argv0 + "/themes/" + options.theme + ".gtkrc";
        if (!std::ifstream(themeFile)) {
            printf ("Current theme in options file is invalid:  %s\nChanging to 21-Gray-Gray\n", options.theme.c_str());
            options.theme = "21-Gray-Gray";
        }
    }

    if (!options.useSystemTheme) {
        std::vector<Glib::ustring> rcfiles;
        rcfiles.push_back (argv0 + "/themes/" + options.theme + ".gtkrc");

        if (options.slimUI) {
            rcfiles.push_back (argv0 + "/themes/slim");
        }

        // Set the font face and size
        Gtk::RC::parse_string (Glib::ustring::compose(
                                   "style \"clearlooks-default\" { font_name = \"%1\" }", options.font));
        Gtk::RC::set_default_files (rcfiles);
    }

    Gtk::Main m(&argc, &argv);

    Glib::ustring icon_path = Glib::build_filename(argv0, "images");
    Glib::RefPtr<Gtk::IconTheme> defaultIconTheme = Gtk::IconTheme::get_default();
    defaultIconTheme->append_search_path(icon_path);

    RTImage::setPaths(options);
    MyExpander::init();  // has to stay AFTER RTImage::setPaths

#ifndef WIN32
    // For an unknown reason, gtkmm 2.22 don't know the gtk-button-images property, while it exists in the documentation...
    // Anyway, the problem was Linux only
    static Glib::RefPtr<Gtk::Settings> settings = Gtk::Settings::get_default();

    if (settings) {
        settings->property_gtk_button_images().set_value(true);
    } else {
        printf("Error: no default settings to update!\n");
    }

#endif

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
    rtengine::procparams::PartialProfile *rawParams = NULL, *imgParams = NULL;
    std::vector<Glib::ustring> inputFiles;
    Glib::ustring outputPath = "";
    std::vector<rtengine::procparams::PartialProfile*> processingParams;
    bool outputDirectory = false;
    bool overwriteFiles = false;
    bool sideProcParams = false;
    bool copyParamsFile = false;
    bool skipIfNoSidecar = false;
    bool useDefault = false;
    unsigned int sideCarFilePos = 0;
    int compression = 92;
    int subsampling = 3;
    int bits = -1;
    std::string outputType = "";
    unsigned errors = 0;

    for( int iArg = 1; iArg < argc; iArg++) {
        if( argv[iArg][0] == '-' ) {
            switch( argv[iArg][1] ) {
            case 'O':
                copyParamsFile = true;

            case 'o': // outputfile or dir
                if( iArg + 1 < argc ) {
                    iArg++;
                    outputPath = safe_filename_to_utf8 (argv[iArg]);

                    if( safe_file_test (outputPath, Glib::FILE_TEST_IS_DIR)) {
                        outputDirectory = true;
                    }
                }

                break;

            case 'p': // processing parameters for all inputs; all set procparams are required, so

                // RT stop if any of them can't be loaded for any reason.
                if( iArg + 1 < argc ) {
                    iArg++;
                    Glib::ustring fname = safe_filename_to_utf8 ( argv[iArg] );

                    if (fname.at(0) == '-') {
                        std::cerr << "Error: filename missing next to the -p switch" << std::endl;
                        deleteProcParams(processingParams);
                        return -3;
                    }

                    rtengine::procparams::PartialProfile* currentParams = new rtengine::procparams::PartialProfile(true);

                    if (!(currentParams->load ( fname ))) {
                        processingParams.push_back(currentParams);
                    } else {
                        std::cerr << "Error: \"" << fname << "\" not found" << std::endl;
                        deleteProcParams(processingParams);
                        return -3;
                    }
                }

                break;

            case 'S':
                skipIfNoSidecar = true;

            case 's': // Processing params next to file (file extension appended)
                sideProcParams = true;
                sideCarFilePos = processingParams.size();
                break;

            case 'd':
                useDefault = true;
                break;

            case 'Y':
                overwriteFiles = true;
                break;

            case 'j':
                if (strlen(argv[iArg]) > 2 && argv[iArg][2] == 's') {
                    if (strlen(argv[iArg]) == 3) {
                        std::cerr << "Error: the -js switch requires a mandatory value!" << std::endl;
                        deleteProcParams(processingParams);
                        return -3;
                    }

                    // looking for the subsampling parameter
                    sscanf(&argv[iArg][3], "%d", &subsampling);

                    if (subsampling < 1 || subsampling > 3) {
                        std::cerr << "Error: the value accompanying the -js switch has to be in the [1-3] range!" << std::endl;
                        deleteProcParams(processingParams);
                        return -3;
                    }
                } else {
                    outputType = "jpg";
                    sscanf(&argv[iArg][2], "%d", &compression);

                    if (compression < 0 || compression > 100) {
                        std::cerr << "Error: the value accompanying the -j switch has to be in the [0-100] range!" << std::endl;
                        deleteProcParams(processingParams);
                        return -3;
                    }
                }

                break;

            case 'b':
                sscanf(&argv[iArg][2], "%d", &bits);

                if (bits != 8 && bits != 16) {
                    std::cerr << "Error: specify -b8 for 8-bit or -b16 for 16-bit output." << std::endl;
                    deleteProcParams(processingParams);
                    return -3;
                }

                break;

            case 't':
                outputType = "tif";
                compression = ((argv[iArg][2] != 'z') ? 0 : 1);
                break;

            case 'n':
                outputType = "png";
                compression = -1;
                break;

            case 'c': // MUST be last option
                while( iArg + 1 < argc ) {
                    iArg++;

                    if( !safe_file_test( safe_filename_to_utf8(argv[iArg]), Glib::FILE_TEST_EXISTS )) {
                        std::cerr << argv[iArg] << " doesn't exist." << std::endl;
                        continue;
                    }

                    if( safe_file_test( safe_filename_to_utf8(argv[iArg]), Glib::FILE_TEST_IS_DIR )) {
                        std::vector<Glib::ustring> names;
                        Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path ( argv[iArg] );
                        safe_build_file_list (dir, names, argv[iArg] );

                        for(size_t iFile = 0; iFile < names.size(); iFile++ ) {
                            if( !safe_file_test( names[iFile] , Glib::FILE_TEST_IS_DIR)) {
                                // skip files without extension and without sidecar files
                                Glib::ustring s(names[iFile]);
                                Glib::ustring::size_type ext = s.find_last_of('.');

                                if( Glib::ustring::npos == ext ) {
                                    continue;
                                }

                                if( ! s.substr(ext).compare( paramFileExtension )) {
                                    continue;
                                }

                                inputFiles.push_back( names[iFile] );
                            }
                        }
                    } else {
                        inputFiles.push_back( safe_filename_to_utf8 (argv[iArg]) );
                    }
                }

                break;
#ifdef WIN32

            case 'w': // This case is handled outside this function
                break;
#endif

            case 'h':
            case '?':
            default: {
                Glib::ustring pparamsExt = paramFileExtension.substr(1);
                std::cout << "<Chevrons> indicate parameters you can change." << std::endl;
                std::cout << "[Square brackets] mean the parameter is not mandatory." << std::endl;
                std::cout << "The pipe symbol | indicates a choice of one or the other." << std::endl;
                std::cout << "The dash symbol - denotes a range of possible values from one to the other." << std::endl;
                cout << std::endl;
                std::cout << "Usage:" << std::endl;
                std::cout << "  " << Glib::path_get_basename(argv[0]) << " <selected dir>     Start File Browser inside directory." << std::endl;
                std::cout << "  " << Glib::path_get_basename(argv[0]) << " <file>             Start Image Editor with file." << std::endl;
                std::cout << "  " << Glib::path_get_basename(argv[0]) << " -c <dir>|<files>   Convert files in batch with default parameters." << std::endl << std::endl;
#ifdef WIN32
                std::cout << "  -w Do not open the Windows console" << std::endl;
#endif
                std::cout << "Other options used with -c (-c must be the last option):" << std::endl;
                std::cout << Glib::path_get_basename(argv[0]) << " [-o <output>|-O <output>] [-s|-S] [-p <files>] [-d] [-j[1-100] [-js<1-3>]|[-b<8|16>] <[-t[z] | [-n]]] [-Y] -c <input>" << std::endl;
                std::cout << "  -o <file>|<dir>  Select output file or directory." << std::endl;
                std::cout << "  -O <file>|<dir>  Select output file or directory and copy " << pparamsExt << " file into it." << std::endl;
                std::cout << "  -s               Include the " << pparamsExt << " file next to the input file (with the same" << std::endl;
                std::cout << "                   name) to build the image parameters," << std::endl;
                std::cout << "                   e.g. for photo.raw there should be a photo.raw." << pparamsExt << " file in" << std::endl;
                std::cout << "                   the same directory. If the file does not exist, internal" << std::endl;
                std::cout << "                   default (neutral) values (not those in Default." << pparamsExt << ") will be" << std::endl;
                std::cout << "                   used." << std::endl;
                std::cout << "  -S               Like -s but skip if the " << pparamsExt << " file does not exist." << std::endl;
                std::cout << "  -p <file.pp3>    Specify " << pparamsExt << " file to be used for all conversions." << std::endl;
                std::cout << "                   You can specify as many -p options as you like (see" << std::endl;
                std::cout << "                   description below)." << std::endl;
                std::cout << "  -d               Use the default raw or non-raw " << pparamsExt << " file as set in" << std::endl;
                std::cout << "                   Preferences > Image Processing > Default Processing Profile" << std::endl;
                std::cout << "  -j[1-100]        Specify output to be JPEG (on by default). Optionally add" << std::endl;
                std::cout << "                   compression 1-100 (default value: 92)." << std::endl;
                std::cout << "  -js<1-3>         Specify the JPEG chroma subsampling parameter, where:" << std::endl;
                std::cout << "                   1 = Best compression: 2x2, 1x1, 1x1 (4:2:0)" << std::endl;
                std::cout << "                       Chroma halved vertically and horizontally." << std::endl;
                std::cout << "                   2 = Balanced:         2x1, 1x1, 1x1 (4:2:2)" << std::endl;
                std::cout << "                       Chroma halved horizontally." << std::endl;
                std::cout << "                   3 = Best quality:     1x1, 1x1, 1x1 (4:4:4)" << std::endl;
                std::cout << "                       No chroma subsampling." << std::endl;
                std::cout << "  -b<8|16>         Specify bit depth per channel (only applies to TIFF and PNG output)." << std::endl;
                std::cout << "  -t[z]            Specify output to be TIFF (16-bit if -b8 is not set)." << std::endl;
                std::cout << "                   Uncompressed by default, or ZIP compression with 'z'" << std::endl;
                std::cout << "  -n               Specify output to be compressed PNG (16-bit if -b8 is not set)." << std::endl;
                std::cout << "  -Y               Overwrite output if present." << std::endl << std::endl;
                std::cout << "Your " << pparamsExt << " files can be incomplete, RawTherapee will set the values as follows:" << std::endl;
                std::cout << "  1- A new profile is created using internal default (neutral) values" << std::endl;
                std::cout << "     (hard-coded into RawTherapee)," << std::endl;
                std::cout << "  2- then overridden by those found in the default raw or non-raw " << pparamsExt << " file" << std::endl;
                std::cout << "     (if -d has been set)," << std::endl;
                std::cout << "  3- then overridden by those found in the " << pparamsExt << " files provided by -p, each one" << std::endl;
                std::cout << "     overriding the previous values," << std::endl;
                std::cout << "  4- then overridden by the sidecar file if -s is set and if the file exists;" << std::endl;
                std::cout << "     the time where the sidecar file is used depends on the position of the -s" << std::endl;
                std::cout << "     switch in the command line relative to the -p parameters," << std::endl;
                std::cout << "     e.g. -p first." << pparamsExt << " -p second." << pparamsExt << " -s -p fourth." << pparamsExt << std::endl;
                return -1;
            }
            }
        } else {
            argv1 = safe_filename_to_utf8 ( argv[iArg] );

            if( outputDirectory ) {
                options.savePathFolder = outputPath;
                options.saveUsePathTemplate = false;
            } else {
                options.saveUsePathTemplate = true;

                if (options.savePathTemplate.empty())
                    // If the save path template is empty, we use its default value
                {
                    options.savePathTemplate = "%p1/converted/%f";
                }
            }

            if (outputType == "jpg") {
                options.saveFormat.format = outputType;
                options.saveFormat.jpegQuality = compression;
                options.saveFormat.jpegSubSamp = subsampling;
            } else if (outputType == "tif") {
                options.saveFormat.format = outputType;
            } else if (outputType == "png") {
                options.saveFormat.format = outputType;
            }

            break;
        }
    }

    if( !argv1.empty() ) {
        return 1;
    }

    if( inputFiles.empty() ) {
        return 2;
    }

    if (useDefault) {
        rawParams = new rtengine::procparams::PartialProfile(true, true);
        Glib::ustring profPath = options.findProfilePath(options.defProfRaw);

        if (options.is_defProfRawMissing() || profPath.empty() || rawParams->load(profPath == DEFPROFILE_INTERNAL ? DEFPROFILE_INTERNAL : Glib::build_filename(profPath, options.defProfRaw.substr(5) + paramFileExtension))) {
            std::cerr << "Error: default raw processing profile not found" << std::endl;
            rawParams->deleteInstance();
            delete rawParams;
            deleteProcParams(processingParams);
            return -3;
        }

        imgParams = new rtengine::procparams::PartialProfile(true);
        profPath = options.findProfilePath(options.defProfImg);

        if (options.is_defProfImgMissing() || profPath.empty() || imgParams->load(profPath == DEFPROFILE_INTERNAL ? DEFPROFILE_INTERNAL : Glib::build_filename(profPath, options.defProfImg.substr(5) + paramFileExtension))) {
            std::cerr << "Error: default non-raw processing profile not found" << std::endl;
            imgParams->deleteInstance();
            delete imgParams;
            rawParams->deleteInstance();
            delete rawParams;
            deleteProcParams(processingParams);
            return -3;
        }
    }

    for( size_t iFile = 0; iFile < inputFiles.size(); iFile++) {

        // Has to be reinstanciated at each profile to have a ProcParams object with default values
        rtengine::procparams::ProcParams currentParams;

        Glib::ustring inputFile = inputFiles[iFile];
        std::cout << "Processing: " << inputFile << std::endl;

        rtengine::InitialImage* ii = NULL;
        rtengine::ProcessingJob* job = NULL;
        int errorCode;
        bool isRaw = false;

        Glib::ustring outputFile;

        if( outputType.empty() ) {
            outputType = "jpg";
        }

        if( outputPath.empty() ) {
            Glib::ustring s = inputFile;
            Glib::ustring::size_type ext = s.find_last_of('.');
            outputFile = s.substr(0, ext) + "." + outputType;
        } else if( outputDirectory ) {
            Glib::ustring s = Glib::path_get_basename( inputFile );
            Glib::ustring::size_type ext = s.find_last_of('.');
            outputFile = outputPath + "/" + s.substr(0, ext) + "." + outputType;
        } else {
            Glib::ustring s = outputPath;
            Glib::ustring::size_type ext = s.find_last_of('.');
            outputFile =  s.substr(0, ext) + "." + outputType;
        }

        if( inputFile == outputFile) {
            std::cerr << "Cannot overwrite: " << inputFile << std::endl;
            continue;
        }

        if( !overwriteFiles && safe_file_test( outputFile , Glib::FILE_TEST_EXISTS ) ) {
            std::cerr << outputFile  << " already exists: use -Y option to overwrite. This image has been skipped." << std::endl;
            continue;
        }

        // Load the image
        isRaw = true;
        Glib::ustring ext = getExtension (inputFile);

        if (ext.lowercase() == "jpg" || ext.lowercase() == "jpeg" || ext.lowercase() == "tif" || ext.lowercase() == "tiff" || ext.lowercase() == "png") {
            isRaw = false;
        }

        ii = rtengine::InitialImage::load ( inputFile, isRaw, &errorCode, NULL );

        if (!ii) {
            errors++;
            std::cerr << "Error loading file: " << inputFile << std::endl;
            continue;
        }

        if (useDefault) {
            if (isRaw) {
                std::cout << "  Merging default raw processing profile" << std::endl;
                rawParams->applyTo(&currentParams);
            } else {
                std::cout << "  Merging default non-raw processing profile" << std::endl;
                imgParams->applyTo(&currentParams);
            }
        }

        bool sideCarFound = false;
        unsigned int i = 0;

        // Iterate the procparams file list in order to build the final ProcParams
        do {
            if (sideProcParams && i == sideCarFilePos) {
                // using the sidecar file
                Glib::ustring sideProcessingParams = inputFile + paramFileExtension;

                // the "load" method don't reset the procparams values anymore, so values found in the procparam file override the one of currentParams
                if( !safe_file_test( sideProcessingParams, Glib::FILE_TEST_EXISTS ) || currentParams.load ( sideProcessingParams )) {
                    std::cerr << "Warning: sidecar file requested but not found for: " << sideProcessingParams << std::endl;
                } else {
                    sideCarFound = true;
                    std::cout << "  Merging sidecar procparams" << std::endl;
                }
            }

            if( processingParams.size() > i  ) {
                std::cout << "  Merging procparams #" << i << std::endl;
                processingParams[i]->applyTo(&currentParams);
            }

            i++;
        } while (i < processingParams.size() + (sideProcParams ? 1 : 0));

        if( sideProcParams && !sideCarFound && skipIfNoSidecar ) {
            delete ii;
            errors++;
            std::cerr << "Error: no sidecar procparams found for: " << inputFile << std::endl;
            continue;
        }

        job = rtengine::ProcessingJob::create (ii, currentParams);

        if( !job ) {
            errors++;
            std::cerr << "Error creating processing for: " << inputFile << std::endl;
            ii->decreaseRef();
            continue;
        }

        // Process image
        rtengine::IImage16* resultImage = rtengine::processImage (job, errorCode, NULL, options.tunnelMetaData);

        if( !resultImage ) {
            errors++;
            std::cerr << "Error processing: " << inputFile << std::endl;
            rtengine::ProcessingJob::destroy( job );
            continue;
        }

        // save image to disk
        if( outputType == "jpg" ) {
            errorCode = resultImage->saveAsJPEG( outputFile, compression, subsampling );
        } else if( outputType == "tif" ) {
            errorCode = resultImage->saveAsTIFF( outputFile, bits, compression == 0  );
        } else if( outputType == "png" ) {
            errorCode = resultImage->saveAsPNG( outputFile, compression, bits );
        } else {
            errorCode = resultImage->saveToFile (outputFile);
        }

        if(errorCode) {
            errors++;
            std::cerr << "Error saving to: " << outputFile << std::endl;
        } else {
            if( copyParamsFile ) {
                Glib::ustring outputProcessingParams = outputFile + paramFileExtension;
                currentParams.save( outputProcessingParams );
            }
        }

        ii->decreaseRef();
        resultImage->free();
    }

    if (imgParams) {
        imgParams->deleteInstance();
        delete imgParams;
    }

    if (rawParams) {
        rawParams->deleteInstance();
        delete rawParams;
    }

    deleteProcParams(processingParams);

    return errors > 0 ? -2 : 0;
}

