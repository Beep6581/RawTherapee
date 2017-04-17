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
#include "../rtengine/icons.h"
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

extern Options options;

// stores path to data files
Glib::ustring argv0;
Glib::ustring creditsPath;
Glib::ustring licensePath;
Glib::ustring argv1;
bool simpleEditor;
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

bool fast_export = false;

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

bool dontLoadCache( int argc, char **argv );

int main(int argc, char **argv)
{
    setlocale(LC_ALL, "");
    setlocale(LC_NUMERIC, "C"); // to set decimal point to "."

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

    bool quickstart = dontLoadCache(argc, argv);

    if (!Options::load (quickstart)) {
        printf("Fatal error!\nThe RT_SETTINGS and/or RT_PATH environment variables are set, but use a relative path. The path must be absolute!\n");
        return -2;
    }

    rtengine::setPaths(options);

    TIFFSetWarningHandler(nullptr);    // avoid annoying message boxes

#ifndef WIN32

    // Move the old path to the new one if the new does not exist
    if (Glib::file_test(Glib::build_filename(options.rtdir, "cache"), Glib::FILE_TEST_IS_DIR) && !Glib::file_test(options.cacheBaseDir, Glib::FILE_TEST_IS_DIR)) {
        g_rename(Glib::build_filename (options.rtdir, "cache").c_str (), options.cacheBaseDir.c_str ());
    }

#endif

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

                if(Console) {
                    AllocConsole();
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
                    std::cout << "RawTherapee, version " << RTVERSION << ", command line" << std::endl;
                    std::cout << "WARNING: closing this window will close RawTherapee!" << std::endl << std::endl;
                }
            }
        }
    }
#endif

    int ret = 0;

    // printing RT's version in all case, particularly useful for the 'verbose' mode, but also for the batch processing
    std::cout << "RawTherapee, version " << RTVERSION << ", command line" << std::endl;
    if (argc > 1) {
        ret = processLineParams(argc, argv);
    }
    else {
        std::cout << "Terminating without anything to do." << std::endl;
    }

#ifdef WIN32
    if(consoleOpened) {
        printf("Press any key to exit RawTherapee\n");
        FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
        getch();
    }
#endif

    return ret;
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


bool dontLoadCache( int argc, char **argv )
{
    for( int iArg = 1; iArg < argc; iArg++) {
        if( argv[iArg][0] == '-' && argv[iArg][1] == 'q' ) {
            return true;
        }
    }

    return false;
}

int processLineParams( int argc, char **argv )
{
    rtengine::procparams::PartialProfile *rawParams = nullptr, *imgParams = nullptr;
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
                    outputPath = fname_to_utf8 (argv[iArg]);

                    if( Glib::file_test (outputPath, Glib::FILE_TEST_IS_DIR)) {
                        outputDirectory = true;
                    }
                }

                break;

            case 'p': // processing parameters for all inputs; all set procparams are required, so

                // RT stop if any of them can't be loaded for any reason.
                if( iArg + 1 < argc ) {
                    iArg++;
                    Glib::ustring fname = fname_to_utf8 (argv[iArg]);

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

            case 'q':
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

            case 'f':
                fast_export = true;
                break;
                
            case 'c': // MUST be last option
                while (iArg + 1 < argc) {
                    iArg++;

                    const auto argument = fname_to_utf8 (argv[iArg]);

                    if (Glib::file_test (argument, Glib::FILE_TEST_IS_REGULAR)) {
                        inputFiles.emplace_back (argument);
                        continue;
                    }

                    if (Glib::file_test (argument, Glib::FILE_TEST_IS_DIR)) {

                        auto dir = Gio::File::create_for_path (argument);
                        if (!dir || !dir->query_exists()) {
                            continue;
                        }

                        try {

                            auto enumerator = dir->enumerate_children ();

                            while (auto file = enumerator->next_file ()) {

                                const auto fileName = Glib::build_filename (argument, file->get_name ());

                                if (Glib::file_test (fileName, Glib::FILE_TEST_IS_DIR)) {
                                    continue;
                                }

                                // skip files without extension and sidecar files
                                auto lastdot = fileName.find_last_of('.');
                                if (lastdot == Glib::ustring::npos) {
                                    continue;
                                }

                                if (fileName.substr (lastdot).compare (paramFileExtension) == 0) {
                                    continue;
                                }

                                inputFiles.emplace_back (fileName);
                            }

                        } catch (Glib::Exception&) {}

                        continue;
                    }

                    std::cerr << "\"" << argument << "\" is neither a regular file nor a directory." << std::endl;
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
                std::cout << "  " << Glib::path_get_basename(argv[0]) << " -c <dir>|<files>   Convert files in batch with default parameters." << std::endl;
                std::cout << "  " << Glib::path_get_basename(argv[0]) << " <other options> -c <dir>|<files>   Convert files in batch with your own settings." << std::endl;
                std::cout << std::endl;
#ifdef WIN32
                std::cout << "  -w Do not open the Windows console" << std::endl;
                std::cout << std::endl;
#endif
                std::cout << "Options:" << std::endl;
                std::cout << "  " << Glib::path_get_basename(argv[0]) << " [-o <output>|-O <output>] [-s|-S] [-p <one.pp3> [-p <two.pp3> ...] ] [-d] [ -j[1-100] [-js<1-3>] | [-b<8|16>] [-t[z] | [-n]] ] [-Y] [-f] -c <input>" << std::endl;
                std::cout << std::endl;
                std::cout << "  -q               Quick Start mode : do not load cached files to speedup start time." << std::endl;
                std::cout << "  -c <files>       Specify one or more input files." << std::endl;
                std::cout << "                   -c must be the last option." << std::endl;
                std::cout << "  -o <file>|<dir>  Set output file or folder." << std::endl;
                std::cout << "                   Saves output file alongside input file if -o is not specified." << std::endl;
                std::cout << "  -O <file>|<dir>  Set output file or folder and copy " << pparamsExt << " file into it." << std::endl;
                std::cout << "                   Saves output file alongside input file if -O is not specified." << std::endl;
                std::cout << "  -s               Use the existing sidecar file to build the processing parameters," << std::endl;
                std::cout << "                   e.g. for photo.raw there should be a photo.raw." << pparamsExt << " file in the same folder." << std::endl;
                std::cout << "                   If the sidecar file does not exist, neutral values will be used." << std::endl;
                std::cout << "  -S               Like -s but skip if the sidecar file does not exist." << std::endl;
                std::cout << "  -p <file.pp3>    Specify processing profile to be used for all conversions." << std::endl;
                std::cout << "                   You can specify as many sets of \"-p <file.pp3>\" options as you like," << std::endl;
                std::cout << "                   each will be built on top of the previous one, as explained below." << std::endl;
                std::cout << "  -d               Use the default raw or non-raw processing profile as set in" << std::endl;
                std::cout << "                   Preferences > Image Processing > Default Processing Profile" << std::endl;
                std::cout << "  -j[1-100]        Specify output to be JPEG (default, if -t and -n are not set)." << std::endl;
                std::cout << "                   Optionally, specify compression 1-100 (default value: 92)." << std::endl;
                std::cout << "  -js<1-3>         Specify the JPEG chroma subsampling parameter, where:" << std::endl;
                std::cout << "                   1 = Best compression:   2x2, 1x1, 1x1 (4:2:0)" << std::endl;
                std::cout << "                       Chroma halved vertically and horizontally." << std::endl;
                std::cout << "                   2 = Balanced (default): 2x1, 1x1, 1x1 (4:2:2)" << std::endl;
                std::cout << "                       Chroma halved horizontally." << std::endl;
                std::cout << "                   3 = Best quality:       1x1, 1x1, 1x1 (4:4:4)" << std::endl;
                std::cout << "                       No chroma subsampling." << std::endl;
                std::cout << "  -b<8|16>         Specify bit depth per channel (default value: 16 for TIFF, 8 for PNG)." << std::endl;
                std::cout << "                   Only applies to TIFF and PNG output, JPEG is always 8." << std::endl;
                std::cout << "  -t[z]            Specify output to be TIFF." << std::endl;
                std::cout << "                   Uncompressed by default, or deflate compression with 'z'." << std::endl;
                std::cout << "  -n               Specify output to be compressed PNG." << std::endl;
                std::cout << "                   Compression is hard-coded to 6." << std::endl;
                std::cout << "  -Y               Overwrite output if present." << std::endl;
                std::cout << "  -f               Use the custom fast-export processing pipeline." << std::endl;
                std::cout << std::endl;
                std::cout << "Your " << pparamsExt << " files can be incomplete, RawTherapee will build the final values as follows:" << std::endl;
                std::cout << "  1- A new processing profile is created using neutral values," << std::endl;
                std::cout << "  2- If the \"-d\" option is set, the values are overridden by those found in" << std::endl;
                std::cout << "     the default raw or non-raw processing profile." << std::endl;
                std::cout << "  3- If one or more \"-p\" options are set, the values are overridden by those" << std::endl;
                std::cout << "     found in these processing profiles." << std::endl;
                std::cout << "  4- If the \"-s\" or \"-S\" options are set, the values are finally overridden by those" << std::endl;
                std::cout << "     found in the sidecar files." << std::endl;
                std::cout << "  The processing profiles are processed in the order specified on the command line." << std::endl;
                return -1;
            }
            }
        } else {
            argv1 = fname_to_utf8 (argv[iArg]);

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

        if (options.is_defProfRawMissing() || profPath.empty() || (profPath != DEFPROFILE_DYNAMIC && rawParams->load(profPath == DEFPROFILE_INTERNAL ? DEFPROFILE_INTERNAL : Glib::build_filename(profPath, options.defProfRaw.substr(5) + paramFileExtension)))) {
            std::cerr << "Error: default raw processing profile not found" << std::endl;
            rawParams->deleteInstance();
            delete rawParams;
            deleteProcParams(processingParams);
            return -3;
        }

        imgParams = new rtengine::procparams::PartialProfile(true);
        profPath = options.findProfilePath(options.defProfImg);

        if (options.is_defProfImgMissing() || profPath.empty() || (profPath != DEFPROFILE_DYNAMIC && imgParams->load(profPath == DEFPROFILE_INTERNAL ? DEFPROFILE_INTERNAL : Glib::build_filename(profPath, options.defProfImg.substr(5) + paramFileExtension)))) {
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

        rtengine::InitialImage* ii = nullptr;
        rtengine::ProcessingJob* job = nullptr;
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

        if( !overwriteFiles && Glib::file_test( outputFile , Glib::FILE_TEST_EXISTS ) ) {
            std::cerr << outputFile  << " already exists: use -Y option to overwrite. This image has been skipped." << std::endl;
            continue;
        }

        // Load the image
        isRaw = true;
        Glib::ustring ext = getExtension (inputFile);

        if (ext.lowercase() == "jpg" || ext.lowercase() == "jpeg" || ext.lowercase() == "tif" || ext.lowercase() == "tiff" || ext.lowercase() == "png") {
            isRaw = false;
        }

        ii = rtengine::InitialImage::load ( inputFile, isRaw, &errorCode, nullptr );

        if (!ii) {
            errors++;
            std::cerr << "Error loading file: " << inputFile << std::endl;
            continue;
        }

        if (useDefault) {
            if (isRaw) {
                if (options.defProfRaw == DEFPROFILE_DYNAMIC) {
                    rawParams->deleteInstance();
                    delete rawParams;
                    rawParams = loadDynamicProfile(ii->getMetaData());
                }
                std::cout << "  Merging default raw processing profile" << std::endl;
                rawParams->applyTo(&currentParams);
             } else {
                if (options.defProfImg == DEFPROFILE_DYNAMIC) {
                    imgParams->deleteInstance();
                    delete imgParams;
                    imgParams = loadDynamicProfile(ii->getMetaData());
                }
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
                if( !Glib::file_test( sideProcessingParams, Glib::FILE_TEST_EXISTS ) || currentParams.load ( sideProcessingParams )) {
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

        job = rtengine::ProcessingJob::create (ii, currentParams, fast_export);

        if( !job ) {
            errors++;
            std::cerr << "Error creating processing for: " << inputFile << std::endl;
            ii->decreaseRef();
            continue;
        }

        // Process image
        rtengine::IImage16* resultImage = rtengine::processImage (job, errorCode, nullptr, options.tunnelMetaData);

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
