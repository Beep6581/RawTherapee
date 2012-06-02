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

#include "config.h"
#include <gtkmm.h>
#include <giomm.h>
#include <iostream>
#include "rtwindow.h"
#include <cstring>
#include <cstdlib>
#include "options.h"
#include "soundman.h"
#include "rtimage.h"
#include "version.h"
#include "extprog.h"

#ifndef WIN32
#include <glibmm/fileutils.h>
#include <glib.h>
#include <glib/gstdio.h>
#endif

#include "../rtengine/safegtk.h"

extern Options options;

// stores path to data files
Glib::ustring argv0;
Glib::ustring creditsPath;
Glib::ustring licensePath;
Glib::ustring argv1;
bool simpleEditor;

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
#ifdef BUILD_BUNDLE
    char exname[512] = {0};
    Glib::ustring exePath;
    // get the path where the rawtherapee executable is stored
    #ifdef WIN32
        WCHAR exnameU[512] = {0};
        GetModuleFileNameW (NULL, exnameU, 512);
        WideCharToMultiByte(CP_UTF8,0,exnameU,-1,exname,512,0,0 );
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
  
   Glib::thread_init();
   gdk_threads_init();
   Gio::init ();

   Options::load ();
   extProgStore->init();
   SoundManager::init();

   if (argc>1){
	   int ret = processLineParams( argc, argv);
	   if( ret <= 0 )
		   return ret;
   }

#ifndef WIN32
   // Move the old path to the new one if the new does not exist
   if (safe_file_test(Glib::build_filename(options.rtdir,"cache"), Glib::FILE_TEST_IS_DIR) && !safe_file_test(options.cacheBaseDir, Glib::FILE_TEST_IS_DIR))
       safe_g_rename(Glib::build_filename(options.rtdir,"cache"), options.cacheBaseDir);
#endif

   simpleEditor=false;
   if( !argv1.empty() )
      if( safe_file_test(argv1, Glib::FILE_TEST_EXISTS) && !safe_file_test(argv1, Glib::FILE_TEST_IS_DIR))
         simpleEditor = true;

   if (!options.useSystemTheme)
   {
       std::vector<Glib::ustring> rcfiles;
       rcfiles.push_back (argv0+"/themes/"+options.theme+".gtkrc");
   	   if (options.slimUI)
           rcfiles.push_back (argv0+"/themes/slim");
       // Set the font face and size
       Gtk::RC::parse_string (Glib::ustring::compose(
          "style \"clearlooks-default\" { font_name = \"%1\" }", options.font));
       Gtk::RC::set_default_files (rcfiles);
   }
   Gtk::Main m(&argc, &argv);

   Glib::ustring icon_path = Glib::build_filename(argv0,"images");
   Glib::RefPtr<Gtk::IconTheme> defaultIconTheme = Gtk::IconTheme::get_default();
   defaultIconTheme->append_search_path(icon_path);

   RTImage::setPaths(options);

#ifndef WIN32
   // For an unknown reason, gtkmm 2.22 don't know the gtk-button-images property, while it exists in the documentation...
   // Anyway, the problem was Linux only
   static Glib::RefPtr<Gtk::Settings> settings = Gtk::Settings::get_default();
   if (settings)
      settings->property_gtk_button_images().set_value(true);
   else
      printf("Error: no default settings to update!\n");
#endif


   RTWindow *rtWindow = new class RTWindow();
   gdk_threads_enter ();

   // alerting users if the default raw and image profiles are missing
   if (options.is_defProfRawMissing()) {
       Gtk::MessageDialog msgd (Glib::ustring::compose(M("OPTIONS_DEFRAW_MISSING"), options.defProfRaw), false, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
       msgd.run ();
   }
   if (options.is_defProfImgMissing()) {
       Gtk::MessageDialog msgd (Glib::ustring::compose(M("OPTIONS_DEFIMG_MISSING"), options.defProfImg), false, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
       msgd.run ();
   }

   // opening the main window
   m.run(*rtWindow);

   gdk_threads_leave ();
   delete rtWindow;
   rtengine::cleanup();
   return 0;
}

void deleteProcParams(std::vector<rtengine::procparams::PartialProfile*> &pparams) {
	for (unsigned int i=0; i<pparams.size(); i++) {
		pparams[i]->deleteInstance();
		delete pparams[i];
		pparams[i] = NULL;
	}
	return;
}

int processLineParams( int argc, char **argv )
{
	rtengine::procparams::PartialProfile *rawParams=NULL, *imgParams=NULL;
	std::vector<Glib::ustring> inputFiles;
	Glib::ustring outputPath = "";
	std::vector<rtengine::procparams::PartialProfile*> processingParams;
	bool isDirectory=false;
	bool outputDirectory=false;
	bool overwriteFiles=false;
	bool sideProcParams=false;
	bool copyParamsFile=false;
	bool skipIfNoSidecar=false;
	bool useDefault=false;
	unsigned int sideCarFilePos = 0;
	int compression=100;
	int bits=-1;
	std::string outputType = "";
	unsigned errors=0;
	for( int iArg=1; iArg<argc; iArg++){
		if( argv[iArg][0]=='-' ){
			switch( argv[iArg][1] ){
			case 'O':
				copyParamsFile = true;
			case 'o': // outputfile or dir
				if( iArg+1 <argc ){
					iArg++;
					outputPath = safe_filename_to_utf8 (argv[iArg]);
					if( safe_file_test (outputPath, Glib::FILE_TEST_IS_DIR))
						outputDirectory=true;
				}
				break;
			case 'p': // processing parameters for all inputs; all set procparams are required, so
				      // RT stop if any of them can't be loaded for any reason.
				if( iArg+1 <argc ){
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
					}
					else {
						std::cerr << "Error: \""<< fname <<"\" not found" << std::endl;
						deleteProcParams(processingParams);
						return -3;
					}
				}
				break;
			case 'S':
				skipIfNoSidecar=true;
			case 's': // Processing params next to file (file extension appended)
				sideProcParams = true;
				sideCarFilePos = processingParams.size();
				break;
			case 'd':
				useDefault = true;
				break;
			case 'Y':
				overwriteFiles =true;
				break;
			case 'j':
				outputType = "jpg";
				sscanf(&argv[iArg][2],"%d",&compression);
				break;
			case 't':
				outputType = "tif";
				compression = ((argv[iArg][2]!='1')?0:1);
				break;
			case 'n':
				outputType = "png";
				compression = -1;
				break;
			case 'c': // MUST be last option
				while( iArg+1 <argc ){
					iArg++;
					if( !safe_file_test( safe_filename_to_utf8(argv[iArg]), Glib::FILE_TEST_EXISTS )){
						std::cerr << argv[iArg] << " doesn't exist."<< std::endl;
						continue;
					}
					if( safe_file_test( safe_filename_to_utf8(argv[iArg]), Glib::FILE_TEST_IS_DIR )){
						isDirectory = true;
						std::vector<Glib::ustring> names;
						Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path ( argv[iArg] );
						safe_build_file_list (dir, names, argv[iArg] );
						for(size_t iFile=0; iFile< names.size(); iFile++ ){
							if( !safe_file_test( names[iFile] , Glib::FILE_TEST_IS_DIR)){
								// skip files without extension and without sidecar files
								Glib::ustring s(names[iFile]);
								Glib::ustring::size_type ext= s.find_last_of('.');
								if( Glib::ustring::npos == ext )
									continue;
								if( ! s.substr(ext).compare( paramFileExtension ))
									continue;
								inputFiles.push_back( names[iFile] );
							}
						}
					}else{
						inputFiles.push_back( safe_filename_to_utf8 (argv[iArg]) );
					}
				}
				break;
			case 'h':
			case '?':
			default:
			{
				Glib::ustring pparamsExt = paramFileExtension.substr(1);
				std::cerr << "RawTherapee, V" << VERSION << std::endl;
				std::cerr << "Copyright (c)2004-2012 Gabor Horvath <hgabor@rawtherapee.com>"<< std::endl << std::endl;
				std::cerr << "Usage:"<< std::endl;
				std::cerr << Glib::path_get_basename(argv[0]) << " [<selected dir>] : start RT GUI browser inside dir."<< std::endl;
				std::cerr << Glib::path_get_basename(argv[0]) << " <file> : start GUI editor with file."<< std::endl;
				std::cerr << Glib::path_get_basename(argv[0]) << " -c <inputDir>|<file list> : convert files in batch with default parameters."<< std::endl<< std::endl;
				std::cerr << "Other options used with -c (that must be last option) "<< std::endl;
				std::cerr << Glib::path_get_basename(argv[0]) <<" [-o <output> | -O <output>] [-s|-S] [-p <file>] [-d] [-j[1-100]|-t|-n] -Y -c <input>"<< std::endl;
				std::cerr << " -o <outputFile>|<outputDir> : select output directory."<< std::endl;
				std::cerr << " -O <outputFile>|<outputDir> : select output dir and copy " << pparamsExt << " file into it"<< std::endl;
				std::cerr << " -s : include the " << pparamsExt << " file next to the input file (with same name) to build the image's parameters"<< std::endl;
				std::cerr << "      ex: for IMG001.NEF there should be IMG001.NEF." << pparamsExt << " in the same dir" << std::endl;
				std::cerr << "      if absent use default" << std::endl;
				std::cerr << " -S : like -s but skip if " << pparamsExt << " file not found." << std::endl;
				std::cerr << " -p <file." << pparamsExt << "> : specify " << pparamsExt << " file to be used for all conversions."<< std::endl;
				std::cerr << "                 you can specify as much -p option as you want (see the note below)."<< std::endl;
				std::cerr << " -d : use the default Raw or Image " << pparamsExt << " file to build the image's parameters"<< std::endl;
				std::cerr << " -j[compression] : specify output to be jpeg.(default)"<< std::endl;
				std::cerr << " -t : specify output to be uncompressed tiff."<< std::endl;
				std::cerr << " -t1: specify output to be compressed tiff."<< std::endl;
				std::cerr << " -n : specify output to be png."<< std::endl;
				std::cerr << " -Y : overwrite output if present."<< std::endl<<std::endl;
				std::cerr << "NOTE: Your " << pparamsExt << " files can be incomplete, RawTherapee will set the values like that:"<< std::endl;
				std::cerr << "   - the " << pparamsExt << " file is built with internal default values;"<< std::endl;
				std::cerr << "   - then overrided by those found in the default Raw or Image " << pparamsExt << " file (if -d has been set);"<< std::endl;
				std::cerr << "   - then overrided by those found in the " << pparamsExt << " files provided by -p,"<< std::endl;
				std::cerr << "     each one overriding the previous values;" << std::endl;
				std::cerr << "   * then overrided by the sidecar file if -s is set and if the file exists;"<< std::endl;
				std::cerr << "     the time where the sidecar file is used depend of the position of the -s switch"<< std::endl;
				std::cerr << "     in the command line regarding to the -p parameters (e.g. \"-p first." << pparamsExt << " -p second." << pparamsExt << " -s -p fourth." << pparamsExt << "\")"<< std::endl;
				return -1;
			}
			}
		}else{
			argv1 = safe_filename_to_utf8 ( argv[iArg] );
			if( outputDirectory ){
				options.savePathFolder = outputPath;
				options.saveUsePathTemplate = false;
			}
			else {
				options.saveUsePathTemplate = true;
				if (options.savePathTemplate.empty())
					// If the save path template is empty, we use its default value
					options.savePathTemplate = "%p1/converted/%f";
			}
			if (outputType == "jpg") {
				options.saveFormat.format = outputType;
				options.saveFormat.jpegQuality = compression;
			} else if (outputType == "tif") {
				options.saveFormat.format = outputType;
			} else if (outputType == "png") {
				options.saveFormat.format = outputType;
			}
			break;
		}
	}
	if( !argv1.empty() )
		return 1;
	if( inputFiles.empty() )
		return 2;

	if (useDefault) {
		rawParams = new rtengine::procparams::PartialProfile(true);
		Glib::ustring profPath = options.findProfilePath(options.defProfRaw);
		if (options.is_defProfRawMissing() || profPath.empty() || rawParams->load(Glib::build_filename(profPath, options.defProfRaw + paramFileExtension))) {
			std::cerr << "Error: default Raw procparams file not found" << std::endl;
			rawParams->deleteInstance();
			delete rawParams;
			deleteProcParams(processingParams);
			return -3;
		}
		imgParams = new rtengine::procparams::PartialProfile(true);
		profPath = options.findProfilePath(options.defProfImg);
		if (options.is_defProfImgMissing() || profPath.empty() || imgParams->load(Glib::build_filename(profPath, options.defProfImg + paramFileExtension))) {
			std::cerr << "Error: default Image procparams file not found" << std::endl;
			imgParams->deleteInstance();
			delete imgParams;
			rawParams->deleteInstance();
			delete rawParams;
			deleteProcParams(processingParams);
			return -3;
		}
	}

	ParamsEdited paramsEdited;
	for( size_t iFile=0; iFile< inputFiles.size(); iFile++){

		// Has to be reinstanciated at each profile to have a ProcParams object with default values
		rtengine::procparams::ProcParams currentParams;

		Glib::ustring inputFile = inputFiles[iFile];
		std::cout << "Processing: " << inputFile << std::endl;

		rtengine::InitialImage* ii=NULL;
		rtengine::ProcessingJob* job =NULL;
		int errorCode;
		bool isRaw=false;

		Glib::ustring outputFile;
		if( outputType.empty() )
			outputType = "jpg";
		if( outputPath.empty() ){
			Glib::ustring s = inputFile;
			Glib::ustring::size_type ext= s.find_last_of('.');
			outputFile = s.substr(0,ext)+ "." + outputType;
		}else if( outputDirectory ){
			Glib::ustring s = Glib::path_get_basename( inputFile );
			Glib::ustring::size_type ext= s.find_last_of('.');
			outputFile = outputPath + "/" + s.substr(0,ext) + "." + outputType;
		}else{
			Glib::ustring s = outputPath;
			Glib::ustring::size_type ext= s.find_last_of('.');
			outputFile =  s.substr(0,ext) + "." + outputType;
		}
		if( inputFile == outputFile){
			std::cerr << "Cannot overwrite: " << inputFile << std::endl;
			continue;
		}
		if( !overwriteFiles && safe_file_test( outputFile , Glib::FILE_TEST_EXISTS ) ){
			std::cerr << outputFile  <<" already exists: use -Y option to overwrite. This image has been skipped." << std::endl;
			continue;
		}

		// Load the image
		ii = rtengine::InitialImage::load ( inputFile, true, &errorCode, NULL );
		if (ii)
			isRaw=true;
		else
			ii = rtengine::InitialImage::load ( inputFile , false, &errorCode, NULL );
		if (!ii) {
			errors++;
			std::cerr << "Error loading file: "<< inputFile << std::endl;
			continue;
		}

		if (useDefault) {
			if (isRaw) {
				std::cout << "  Merging default Raw profile" << std::endl;
				rawParams->applyTo(&currentParams);
			}
			else {
				std::cout << "  Merging default Image profile" << std::endl;
				imgParams->applyTo(&currentParams);
			}
		}

		bool sideCarFound = false;
		unsigned int i=0;
		// Iterate the procparams file list in order to build the final ProcParams
		do {
			if (sideProcParams && i==sideCarFilePos) {
				// using the sidecar file
				Glib::ustring sideProcessingParams = inputFile + paramFileExtension;
				// the "load" method don't reset the procparams values anymore, so values found in the procparam file override the one of currentParams
				if( !safe_file_test( sideProcessingParams, Glib::FILE_TEST_EXISTS ) || currentParams.load ( sideProcessingParams ))
					std::cerr << "Warning: sidecar file requested but not found for: "<< sideProcessingParams << std::endl;
				else {
					sideCarFound = true;
					std::cout << "  Merging sidecar procparams" << std::endl;
				}
			}
			if( processingParams.size()>i  ) {
				std::cout << "  Merging procparams #" << i << std::endl;
				processingParams[i]->applyTo(&currentParams);
			}
			i++;
		} while (i < processingParams.size()+(sideProcParams?1:0));

		if( sideProcParams && !sideCarFound && skipIfNoSidecar ){
			delete ii;
			errors++;
			std::cerr << "Error: no sidecar procparams found for: "<< inputFile << std::endl;
			continue;
		}

		job = rtengine::ProcessingJob::create (ii, currentParams);
		if( !job ){
			errors++;
			std::cerr << "Error creating processing for: "<< inputFile << std::endl;
			ii->decreaseRef();
			continue;
		}

		// Process image
		rtengine::IImage16* resultImage = rtengine::processImage (job, errorCode, NULL, options.tunnelMetaData);
        if( !resultImage ){
        	errors++;
        	std::cerr << "Error processing: "<< inputFile << std::endl;
        	rtengine::ProcessingJob::destroy( job );
    		continue;
        }
		// save image to disk
		if( outputType=="jpg" )
			errorCode = resultImage->saveAsJPEG( outputFile, compression );
		else if( outputType=="tif" )
			errorCode = resultImage->saveAsTIFF( outputFile, bits, compression==0  );
		else if( outputType=="png" )
			errorCode = resultImage->saveAsPNG( outputFile,compression, bits );
		else
			errorCode = resultImage->saveToFile (outputFile);

		if(errorCode){
			errors++;
			std::cerr << "Error saving to: "<< outputFile << std::endl;
		}else{
			if( copyParamsFile ){
			   Glib::ustring outputProcessingParams = outputFile + paramFileExtension;
			   currentParams.save( outputProcessingParams );
			}
		}

		ii->decreaseRef();
		resultImage->free();
	}

	if (imgParams) { imgParams->deleteInstance(); delete imgParams; }
	if (rawParams) { rawParams->deleteInstance(); delete rawParams; }
	deleteProcParams(processingParams);

	return errors>0?-2:0;
}
