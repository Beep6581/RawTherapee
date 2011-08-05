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
#include <options.h>
#include <stdio.h>
#include <glib/gstdio.h>
#include <sstream>
#include <multilangmgr.h>
#include <safekeyfile.h>
#include <addsetids.h>
#include <guiutils.h>
#include <safegtk.h>
#include "version.h"

#ifdef WIN32
#include <windows.h>
#include <Shlobj.h>
#endif

Options options;
Glib::ustring versionString      = VERSION;
Glib::ustring paramFileExtension = ".pp3";

Options::Options () {

    setDefaults ();
}

const char *DefaultLanguage = "English (US)";

void Options::setDefaults () {

	font = "sans, 10";
    windowWidth = 900;
    windowHeight = 560;
    windowMaximized = false;
    firstRun = true;
    savesParamsAtExit = true;
    saveFormat.format = "jpg";
    saveFormat.jpegQuality = 100;
    saveFormat.pngCompression = 6;
    saveFormat.pngBits = 8;
    saveFormat.tiffBits = 8;
    saveFormat.tiffUncompressed = true;
    saveFormat.saveParams = true;		// was false

    saveFormatBatch.format = "jpg";
    saveFormatBatch.jpegQuality = 100;
    saveFormatBatch.pngCompression = 6;
    saveFormatBatch.pngBits = 8;
    saveFormatBatch.tiffBits = 8;
    saveFormatBatch.tiffUncompressed = true;
    saveFormatBatch.saveParams = true;		// was false

    savePathTemplate = "%p1/converted/%f";
    savePathFolder = "";
    saveUsePathTemplate = true;
    defProfRaw = "default";
    defProfImg = "neutral";
    dateFormat = "%y-%m-%d";
    adjusterDelay = 0;
    startupDir = STARTUPDIR_LAST;		// was STARTUPDIR_HOME ; an empty startupPath is now correctly handled (open in the Home dir)
    startupPath = "";
    profilePath = "profiles";
    dirBrowserWidth = 200;
    dirBrowserHeight = 150;
	preferencesWidth = 0;
	preferencesHeight = 0;
    toolPanelWidth = 300;
    browserToolPanelWidth = 300;
    browserToolPanelHeight = 300;
    historyPanelWidth = 230;			// was 150
    lastScale = 5;						// was 4
    lastCropSize = 1;
    fbOnlyRaw = false;
    fbShowDateTime = true;
    fbShowBasicExif = true;
    fbShowHidden = false;
    fbArrangement = 2;					// was 0
    multiUser = true;
    version = VERSION;
    thumbSize = 240;					// was 80
    thumbSizeTab = 80;
    showHistory = true;
    showFilePanelState = 0;				// Not used anymore ; was the thumb strip state
    showInfo = true;					// was false
    cropPPI = 600;						// was 300
    showClippedHighlights = false;
    showClippedShadows = false;
    highlightThreshold = 253;			// was 254
    shadowThreshold = 8;				// was 0
    bgcolor = 0;
    blinkClipped = false;				// was true
    language = DefaultLanguage;
    languageAutoDetect= langMgr.isOSLanguageDetectSupported();
    lastSaveAsPath = "";
    overwriteOutputFile = false;		// if TRUE, existing output JPGs/PNGs are overwritten, instead of adding ..-1.jpg, -2.jpg etc.
    theme = "17-Gray-Red";
    slimUI = false;		// TODO: Should this be TRUE for worst case screen resolution or FALSE for nicer interface by default ???
    useSystemTheme = true;
    maxThumbnailHeight = 400;
    maxCacheEntries = 20000;			// was 10000
    thumbnailFormat = FT_Custom;		// was FT_Custom16
    thumbInterp = 1;
    autoSuffix = false;
    saveParamsFile = true;				// was false, but saving the procparams files next to the file make more sense when reorganizing file tree than in a cache
    saveParamsCache = false;			// there's no need to save the procparams files in a cache if saveParamsFile is true
    paramsLoadLocation = PLL_Input;		// was PLL_Cache
    procQueueEnabled = false;			// was true
    gimpDir = "";
    psDir = "";
    customEditorProg = "";
    editorToSendTo = 1;
    liveThumbnails = true;
    favoriteDirs.clear();
    tpOpen.clear ();
    //crvOpen.clear ();
    parseExtensions.clear ();
    parseExtensionsEnabled.clear ();
    parsedExtensions.clear ();
    renameUseTemplates = false;
    renameTemplates.clear ();
    thumbnailZoomRatios.clear ();
    thumbnailZoomRatios.push_back (0.2);
    thumbnailZoomRatios.push_back (0.3);
    thumbnailZoomRatios.push_back (0.45);
    thumbnailZoomRatios.push_back (0.6);
    thumbnailZoomRatios.push_back (0.8);
    thumbnailZoomRatios.push_back (1.0);
    overlayedFileNames = true;
    internalThumbIfUntouched = true; 	// if TRUE, only fast, internal preview images are taken if the image is not edited yet
    showFileNames = true;
    tabbedUI = true;					// was false;
    multiDisplayMode = 0;
    tunnelMetaData = false;
    histogramPosition = 2;
    showProfileSelector = true;
    FileBrowserToolbarSingleRow = true;
    menuGroupRank = true;
    menuGroupLabel = true;
    menuGroupFileOperations = true;
    menuGroupProfileOperations = true;

    cutOverlayBrush = std::vector<double> (4);
    cutOverlayBrush[3] = 0.667;  // :-p

    sndEnable=true;
    sndLngEditProcDoneSecs=3.0;

    int babehav[] = {
			1,  // ADDSET_TC_EXPCOMP
			1,  // ADDSET_TC_BRIGHTNESS
			1,  // ADDSET_TC_BLACKLEVEL
			1,  // ADDSET_TC_CONTRAST
			1,  // ADDSET_SH_HIGHLIGHTS
			1,  // ADDSET_SH_SHADOWS
			1,  // ADDSET_SH_LOCALCONTRAST
			1,  // ADDSET_LC_BRIGHTNESS
			1,  // ADDSET_LC_CONTRAST
			0,  // ADDSET_SHARP_AMOUNT
			0,  // ADDSET_LD_EDGETOLERANCE
			1,  // ADDSET_WB_TEMPERATURE
			1,  // ADDSET_WB_GREEN
			0,  // ADDSET_CBOOST_AMOUNT
			0,  // ADDSET_CS_BLUEYELLOW
			0,  // ADDSET_CS_GREENMAGENTA
			0,  // ADDSET_ROTATE_DEGREE
			0,  // ADDSET_DIST_AMOUNT
			0,  // ADDSET_PERSPECTIVE
			0,  // ADDSET_CA
			1,  // ADDSET_VIGN_AMOUNT
			1,  // ADDSET_LC_SATURATION
			1,  // ADDSET_TC_SATURATION
			1,  // ADDSET_TC_HLCOMPAMOUNT
			0,  // ADDSET_TC_HLCOMPTHRESH
			1,  // ADDSET_TC_SHCOMP
			1,  // ADDSET_DIRPYREQ
			0,  // ADDSET_DIRPYRDN_CHLUM
			0,  // ADDSET_DIRPYRDN_GAMMA
			1,  // ADDSET_CHMIXER
			1,  // ADDSET_PREPROCESS_GREENEQUIL
			0,  // ADDSET_PREPROCESS_LINEDENOISE
			0,  // ADDSET_RAWCACORR
			1,  // ADDSET_RAWEXPOS_LINEAR
			1,  // ADDSET_RAWEXPOS_PRESER
			1,  // ADDSET_RAWEXPOS_BLACKS
			0,  // ADDSET_CLAR_STREN
			0,  // ADDSET_CLAR_MLSTREN
			0,  // ADDSET_CLAR_PASS
			0   // ADDSET_CLAR_RELAT
	};
    baBehav = std::vector<int> (babehav, babehav+ADDSET_PARAM_NUM);
    
    rtSettings.dualThreadEnabled = true;
    rtSettings.darkFramesPath = "";
	rtSettings.flatFieldsPath = "";
#ifdef WIN32
	const gchar* sysRoot = g_getenv("SystemRoot");  // Returns e.g. "c:\Windows"
	if (sysRoot!=NULL) 
		rtSettings.iccDirectory = Glib::ustring(sysRoot) + Glib::ustring("\\System32\\spool\\drivers\\color");
	else 
		rtSettings.iccDirectory = "C:\\WINDOWS\\System32\\spool\\drivers\\color";
#else
    rtSettings.iccDirectory = "/usr/share/color/icc";
#endif
    rtSettings.colorimetricIntent = 1;
    rtSettings.monitorProfile = "";
	rtSettings.autoMonitorProfile = false;
    rtSettings.LCMSSafeMode = true;
    rtSettings.adobe = "AdobeRGB1998"; // put the name of yours profiles (here windows)
    rtSettings.prophoto = "ProPhoto"; // these names appear in the menu "output profile"
    rtSettings.widegamut = "WideGamutRGB";
    rtSettings.srgb = "sRGB Color Space Profile";
    rtSettings.bruce = "Bruce";
    rtSettings.beta = "BetaRGB";
    rtSettings.best = "BestRGB";

    rtSettings.verbose = false;
}

Options* Options::copyFrom (Options* other) {

    *this = *other;
	return this;
}

void Options::filterOutParsedExtensions () {
	parsedExtensions.clear();
	for (unsigned int i=0; i<parseExtensions.size(); i++)
		if (parseExtensionsEnabled[i]) parsedExtensions.push_back(parseExtensions[i].lowercase());
}

int Options::readFromFile (Glib::ustring fname) {

    rtengine::SafeKeyFile keyFile;

    try {
    	if( !safe_file_test(fname,Glib::FILE_TEST_EXISTS))
    		return 1;
        if (!keyFile.load_from_file (fname)) 
            return 1;
    }
    catch (Glib::FileError) {
        return 1;
    }

    setDefaults ();
    
if (keyFile.has_group ("General")) {
    if (keyFile.has_key ("General", "TabbedEditor"))    tabbedUI= keyFile.get_boolean ("General", "TabbedEditor");    
    if (keyFile.has_key ("General", "StartupDirectory")){
    	if( keyFile.get_string ("General", "StartupDirectory") == "home")          startupDir = STARTUPDIR_HOME;
        else if ( keyFile.get_string ("General", "StartupDirectory") == "current") startupDir = STARTUPDIR_CURRENT;
        else if ( keyFile.get_string ("General", "StartupDirectory") == "last")    startupDir = STARTUPDIR_LAST;
        else if ( keyFile.get_string ("General", "StartupDirectory") == "custom")  startupDir = STARTUPDIR_CUSTOM;
    }
        
    if (keyFile.has_key ("General", "StartupPath"))      startupPath     = keyFile.get_string ("General", "StartupPath");
    if (keyFile.has_key ("General", "DateFormat"))       dateFormat      = keyFile.get_string ("General", "DateFormat");
    if (keyFile.has_key ("General", "AdjusterDelay"))    adjusterDelay   = keyFile.get_integer ("General", "AdjusterDelay");
    if (keyFile.has_key ("General", "StoreLastProfile")) savesParamsAtExit = keyFile.get_boolean ("General", "StoreLastProfile");
    if (keyFile.has_key ("General", "DualProcSupport"))  rtSettings.dualThreadEnabled = keyFile.get_boolean ("General", "DualProcSupport");
    if (keyFile.has_key ("General", "MultiUser"))        multiUser       = keyFile.get_boolean ("General", "MultiUser");
    if (keyFile.has_key ("General", "Version"))          version         = keyFile.get_string ("General", "Version");
    if (keyFile.has_key ("General", "Language"))         language        = keyFile.get_string ("General", "Language");
    if (keyFile.has_key ("General", "LanguageAutoDetect")) languageAutoDetect = keyFile.get_boolean ("General", "LanguageAutoDetect");
    if (keyFile.has_key ("General", "Theme"))            theme           = keyFile.get_string ("General", "Theme");
    if (keyFile.has_key ("General", "SlimUI"))           slimUI          = keyFile.get_boolean ("General", "SlimUI");
    if (keyFile.has_key ("General", "UseSystemTheme"))   useSystemTheme  = keyFile.get_boolean ("General", "UseSystemTheme");
    if (keyFile.has_key ("General", "FirstRun"))         firstRun        = keyFile.get_boolean ("General", "FirstRun");
    if( keyFile.has_key ("General", "DarkFramesPath"))   rtSettings.darkFramesPath = keyFile.get_string("General", "DarkFramesPath");
    if( keyFile.has_key ("General", "FlatFieldsPath"))   rtSettings.flatFieldsPath = keyFile.get_string("General", "FlatFieldsPath");
    if( keyFile.has_key ("General", "Verbose"))          rtSettings.verbose = keyFile.get_boolean ( "General", "Verbose");
}

if (keyFile.has_group ("External Editor")) { 
    if (keyFile.has_key ("External Editor", "EditorKind"))      editorToSendTo   = keyFile.get_integer ("External Editor", "EditorKind");
    if (keyFile.has_key ("External Editor", "GimpDir"))         gimpDir          = keyFile.get_string  ("External Editor", "GimpDir");
    if (keyFile.has_key ("External Editor", "PhotoshopDir"))    psDir            = keyFile.get_string  ("External Editor", "PhotoshopDir");
    if (keyFile.has_key ("External Editor", "CustomEditor"))    customEditorProg = keyFile.get_string  ("External Editor", "CustomEditor");
}

if (keyFile.has_group ("Output")) { 
    if (keyFile.has_key ("Output", "Format"))           saveFormat.format          = keyFile.get_string ("Output", "Format");
    if (keyFile.has_key ("Output", "JpegQuality"))      saveFormat.jpegQuality     = keyFile.get_integer ("Output", "JpegQuality");
    if (keyFile.has_key ("Output", "PngCompression"))   saveFormat.pngCompression  = keyFile.get_integer ("Output", "PngCompression");
    if (keyFile.has_key ("Output", "PngBps"))           saveFormat.pngBits         = keyFile.get_integer ("Output", "PngBps");
    if (keyFile.has_key ("Output", "TiffBps"))          saveFormat.tiffBits        = keyFile.get_integer ("Output", "TiffBps");
    if (keyFile.has_key ("Output", "TiffUncompressed")) saveFormat.tiffUncompressed= keyFile.get_boolean ("Output", "TiffUncompressed");
    if (keyFile.has_key ("Output", "SaveProcParams"))   saveFormat.saveParams      = keyFile.get_boolean ("Output", "SaveProcParams");


    if (keyFile.has_key ("Output", "FormatBatch"))           saveFormatBatch.format          = keyFile.get_string ("Output", "FormatBatch");
    if (keyFile.has_key ("Output", "JpegQualityBatch"))      saveFormatBatch.jpegQuality     = keyFile.get_integer ("Output", "JpegQualityBatch");
    if (keyFile.has_key ("Output", "PngCompressionBatch"))   saveFormatBatch.pngCompression  = keyFile.get_integer ("Output", "PngCompressionBatch");
    if (keyFile.has_key ("Output", "PngBpsBatch"))           saveFormatBatch.pngBits         = keyFile.get_integer ("Output", "PngBpsBatch");
    if (keyFile.has_key ("Output", "TiffBpsBatch"))          saveFormatBatch.tiffBits        = keyFile.get_integer ("Output", "TiffBpsBatch");
    if (keyFile.has_key ("Output", "TiffUncompressedBatch")) saveFormatBatch.tiffUncompressed= keyFile.get_boolean ("Output", "TiffUncompressedBatch");
    if (keyFile.has_key ("Output", "SaveProcParamsBatch"))   saveFormatBatch.saveParams      = keyFile.get_boolean ("Output", "SaveProcParamsBatch");

    if (keyFile.has_key ("Output", "Path"))             savePathTemplate           = keyFile.get_string ("Output", "Path");
    if (keyFile.has_key ("Output", "PathTemplate"))     savePathTemplate           = keyFile.get_string ("Output", "PathTemplate");
    if (keyFile.has_key ("Output", "PathFolder"))       savePathFolder             = keyFile.get_string ("Output", "PathFolder");
    if (keyFile.has_key ("Output", "AutoSuffix"))       autoSuffix                 = keyFile.get_boolean("Output", "AutoSuffix");
    if (keyFile.has_key ("Output", "UsePathTemplate"))  saveUsePathTemplate        = keyFile.get_boolean("Output", "UsePathTemplate");
    if (keyFile.has_key ("Output", "LastSaveAsPath"))   lastSaveAsPath             = keyFile.get_string ("Output", "LastSaveAsPath");
	if (keyFile.has_key ("Output", "OverwriteOutputFile"))  overwriteOutputFile    = keyFile.get_boolean("Output", "OverwriteOutputFile");
	if (keyFile.has_key ("Output", "TunnelMetaData"))   tunnelMetaData             = keyFile.get_boolean("Output", "TunnelMetaData");
}

if (keyFile.has_group ("Profiles")) { 
    if (keyFile.has_key ("Profiles", "Directory"))      profilePath     = keyFile.get_string ("Profiles", "Directory");
    if (keyFile.has_key ("Profiles", "RawDefault"))     defProfRaw      = keyFile.get_string ("Profiles", "RawDefault");
    if (keyFile.has_key ("Profiles", "ImgDefault"))     defProfImg      = keyFile.get_string ("Profiles", "ImgDefault");
    if (keyFile.has_key ("Profiles", "SaveParamsWithFile")) saveParamsFile  = keyFile.get_boolean ("Profiles", "SaveParamsWithFile");
    if (keyFile.has_key ("Profiles", "SaveParamsToCache"))  saveParamsCache = keyFile.get_boolean ("Profiles", "SaveParamsToCache");
    if (keyFile.has_key ("Profiles", "LoadParamsFromLocation")) paramsLoadLocation = (PPLoadLocation)keyFile.get_integer ("Profiles", "LoadParamsFromLocation");
    if (keyFile.has_key ("Profiles", "CustomProfileBuilder"))   customProfileBuilder = keyFile.get_string  ("Profiles", "CustomProfileBuilder");
}

if (keyFile.has_group ("File Browser")) { 
    if (keyFile.has_key ("File Browser", "ThumbnailSize"))      thumbSize          = keyFile.get_integer ("File Browser", "ThumbnailSize");
    if (keyFile.has_key ("File Browser", "ThumbnailSizeTab"))      thumbSizeTab    = keyFile.get_integer ("File Browser", "ThumbnailSizeTab");
    if (keyFile.has_key ("File Browser", "BrowseOnlyRaw"))      fbOnlyRaw          = keyFile.get_boolean ("File Browser", "BrowseOnlyRaw");
    if (keyFile.has_key ("File Browser", "BrowserShowsDate"))   fbShowDateTime     = keyFile.get_boolean ("File Browser", "BrowserShowsDate");
    if (keyFile.has_key ("File Browser", "BrowserShowsExif"))   fbShowBasicExif    = keyFile.get_boolean ("File Browser", "BrowserShowsExif");
    if (keyFile.has_key ("File Browser", "BrowserShowsHidden")) fbShowHidden       = keyFile.get_boolean ("File Browser", "BrowserShowsHidden");
    if (keyFile.has_key ("File Browser", "MaxPreviewHeight"))   maxThumbnailHeight = keyFile.get_integer ("File Browser", "MaxPreviewHeight");
    if (keyFile.has_key ("File Browser", "MaxCacheEntries"))    maxCacheEntries    = keyFile.get_integer ("File Browser", "MaxCacheEntries");
    if (keyFile.has_key ("File Browser", "ThumbnailFormat"))    thumbnailFormat    = (ThFileType)keyFile.get_integer ("File Browser", "ThumbnailFormat");
    if (keyFile.has_key ("File Browser", "ParseExtensions"))    parseExtensions    = keyFile.get_string_list ("File Browser", "ParseExtensions");
    if (keyFile.has_key ("File Browser", "ParseExtensionsEnabled")) parseExtensionsEnabled    = keyFile.get_integer_list ("File Browser", "ParseExtensionsEnabled");
    if (keyFile.has_key ("File Browser", "ThumbnailArrangement")) fbArrangement    = keyFile.get_integer ("File Browser", "ThumbnailArrangement");
    if (keyFile.has_key ("File Browser", "ThumbnailInterpolation")) thumbInterp    = keyFile.get_integer ("File Browser", "ThumbnailInterpolation");
    if (keyFile.has_key ("File Browser", "LiveThumbnails"))     liveThumbnails     = keyFile.get_boolean ("File Browser", "LiveThumbnails");
    if (keyFile.has_key ("File Browser", "FavoriteDirs"))       favoriteDirs       = keyFile.get_string_list ("File Browser", "FavoriteDirs");
    if (keyFile.has_key ("File Browser", "RenameTemplates"))    renameTemplates    = keyFile.get_string_list ("File Browser", "RenameTemplates");
    if (keyFile.has_key ("File Browser", "RenameUseTemplates")) renameUseTemplates = keyFile.get_boolean ("File Browser", "RenameUseTemplates");
    if (keyFile.has_key ("File Browser", "ThumbnailZoomRatios"))thumbnailZoomRatios= keyFile.get_double_list ("File Browser", "ThumbnailZoomRatios");
    if (keyFile.has_key ("File Browser", "OverlayedFileNames")) overlayedFileNames = keyFile.get_boolean ("File Browser", "OverlayedFileNames");
    if (keyFile.has_key ("File Browser", "ShowFileNames"))      showFileNames = keyFile.get_boolean ("File Browser", "ShowFileNames");
    if (keyFile.has_key ("File Browser", "InternalThumbIfUntouched")) internalThumbIfUntouched = keyFile.get_boolean ("File Browser", "InternalThumbIfUntouched");
    if (keyFile.has_key ("File Browser", "menuGroupRank")) menuGroupRank = keyFile.get_boolean ("File Browser", "menuGroupRank");
    if (keyFile.has_key ("File Browser", "menuGroupLabel")) menuGroupLabel = keyFile.get_boolean ("File Browser", "menuGroupLabel");
    if (keyFile.has_key ("File Browser", "menuGroupFileOperations")) menuGroupFileOperations = keyFile.get_boolean ("File Browser", "menuGroupFileOperations");
    if (keyFile.has_key ("File Browser", "menuGroupProfileOperations")) menuGroupProfileOperations = keyFile.get_boolean ("File Browser", "menuGroupProfileOperations");
}

if (keyFile.has_group ("Clipping Indication")) { 
    if (keyFile.has_key ("Clipping Indication", "HighlightThreshold"))  highlightThreshold= keyFile.get_integer ("Clipping Indication", "HighlightThreshold");
    if (keyFile.has_key ("Clipping Indication", "ShadowThreshold"))     shadowThreshold   = keyFile.get_integer ("Clipping Indication", "ShadowThreshold");
    if (keyFile.has_key ("Clipping Indication", "BlinkClipped"))        blinkClipped      = keyFile.get_boolean ("Clipping Indication", "BlinkClipped");
}

if (keyFile.has_group ("GUI")) { 
    if (keyFile.has_key ("GUI", "Font"))            font            = keyFile.get_string  ("GUI", "Font");
    if (keyFile.has_key ("GUI", "WindowWidth"))     windowWidth     = keyFile.get_integer ("GUI", "WindowWidth");
    if (keyFile.has_key ("GUI", "WindowHeight"))    windowHeight    = keyFile.get_integer ("GUI", "WindowHeight");
    if (keyFile.has_key ("GUI", "WindowMaximized")) windowMaximized = keyFile.get_boolean ("GUI", "WindowMaximized");
    if (keyFile.has_key ("GUI", "DirBrowserWidth"))     dirBrowserWidth          = keyFile.get_integer ("GUI", "DirBrowserWidth");
    if (keyFile.has_key ("GUI", "DirBrowserHeight"))    dirBrowserHeight         = keyFile.get_integer ("GUI", "DirBrowserHeight");
	if (keyFile.has_key ("GUI", "PreferencesWidth"))    preferencesWidth         = keyFile.get_integer ("GUI", "PreferencesWidth");
	if (keyFile.has_key ("GUI", "PreferencesHeight"))   preferencesHeight        = keyFile.get_integer ("GUI", "PreferencesHeight"); 
    if (keyFile.has_key ("GUI", "SaveAsDialogWidth"))   saveAsDialogWidth        = keyFile.get_integer ("GUI", "SaveAsDialogWidth");
    if (keyFile.has_key ("GUI", "SaveAsDialogHeight"))  saveAsDialogHeight       = keyFile.get_integer ("GUI", "SaveAsDialogHeight");
    if (keyFile.has_key ("GUI", "ToolPanelWidth"))      toolPanelWidth           = keyFile.get_integer ("GUI", "ToolPanelWidth");
    if (keyFile.has_key ("GUI", "BrowserToolPanelWidth"))browserToolPanelWidth   = keyFile.get_integer ("GUI", "BrowserToolPanelWidth");
    if (keyFile.has_key ("GUI", "BrowserToolPanelHeight"))browserToolPanelHeight = keyFile.get_integer ("GUI", "BrowserToolPanelHeight");
    if (keyFile.has_key ("GUI", "HistoryPanelWidth"))   historyPanelWidth = keyFile.get_integer ("GUI", "HistoryPanelWidth");
    if (keyFile.has_key ("GUI", "LastPreviewScale"))    lastScale         = keyFile.get_integer ("GUI", "LastPreviewScale");
    if (keyFile.has_key ("GUI", "LastCropSize"))        lastCropSize      = keyFile.get_integer ("GUI", "LastCropSize");
    if (keyFile.has_key ("GUI", "ShowHistory"))         showHistory       = keyFile.get_boolean ("GUI", "ShowHistory");
    if (keyFile.has_key ("GUI", "ShowFilePanelState"))  showFilePanelState= keyFile.get_integer ("GUI", "ShowFilePanelState");
    if (keyFile.has_key ("GUI", "ShowInfo"))            showInfo          = keyFile.get_boolean ("GUI", "ShowInfo");
    if (keyFile.has_key ("GUI", "MainNBVertical"))      mainNBVertical    = keyFile.get_boolean ("GUI", "MainNBVertical");
    if (keyFile.has_key ("GUI", "ShowClippedHighlights"))showClippedHighlights = keyFile.get_boolean ("GUI", "ShowClippedHighlights");
    if (keyFile.has_key ("GUI", "ShowClippedShadows"))  showClippedShadows= keyFile.get_boolean ("GUI", "ShowClippedShadows");
    if (keyFile.has_key ("GUI", "FrameColor"))          bgcolor           = keyFile.get_integer ("GUI", "FrameColor");
    if (keyFile.has_key ("GUI", "ProcessingQueueEnbled"))procQueueEnabled = keyFile.get_boolean ("GUI", "ProcessingQueueEnbled");
    if (keyFile.has_key ("GUI", "ToolPanelsExpanded"))  tpOpen            = keyFile.get_integer_list ("GUI", "ToolPanelsExpanded");
    if (keyFile.has_key ("GUI", "MultiDisplayMode"))    multiDisplayMode  = keyFile.get_integer ("GUI", "MultiDisplayMode");
    //if (keyFile.has_key ("GUI", "CurvePanelsExpanded")) crvOpen           = keyFile.get_integer_list ("GUI", "CurvePanelsExpanded");
    if (keyFile.has_key ("GUI", "CutOverlayBrush"))     cutOverlayBrush   = keyFile.get_double_list ("GUI", "CutOverlayBrush");
    if (keyFile.has_key ("GUI", "HistogramPosition"))   histogramPosition   = keyFile.get_integer ("GUI", "HistogramPosition");
    if (keyFile.has_key ("GUI", "ShowProfileSelector")) showProfileSelector = keyFile.get_boolean ("GUI", "ShowProfileSelector");
    if (keyFile.has_key ("GUI", "FileBrowserToolbarSingleRow")) FileBrowserToolbarSingleRow = keyFile.get_boolean ("GUI", "FileBrowserToolbarSingleRow");
}



if (keyFile.has_group ("Crop Settings")) { 
    if (keyFile.has_key ("Crop Settings", "PPI"))       cropPPI      = keyFile.get_integer ("Crop Settings", "PPI");
}

if (keyFile.has_group ("Color Management")) { 
    if (keyFile.has_key ("Color Management", "ICCDirectory"))   rtSettings.iccDirectory         = keyFile.get_string ("Color Management", "ICCDirectory");
    if (keyFile.has_key ("Color Management", "MonitorProfile")) rtSettings.monitorProfile       = keyFile.get_string ("Color Management", "MonitorProfile");
    if (keyFile.has_key ("Color Management", "AutoMonitorProfile")) rtSettings.autoMonitorProfile = keyFile.get_boolean ("Color Management", "AutoMonitorProfile");

    if (keyFile.has_key ("Color Management", "Intent"))         rtSettings.colorimetricIntent   = keyFile.get_integer("Color Management", "Intent");

    // Disabled (default is true) till issues are sorted out
    //if (keyFile.has_key ("Color Management", "LCMSSafeMode")) rtSettings.LCMSSafeMode = keyFile.get_boolean ("Color Management", "LCMSSafeMode");
}

if (keyFile.has_group ("Batch Processing")) { 
    if (keyFile.has_key ("Batch Processing", "AdjusterBehavior")) baBehav = keyFile.get_integer_list ("Batch Processing", "AdjusterBehavior");
}

if (keyFile.has_group ("Sounds")) { 
    if (keyFile.has_key ("Sounds", "Enable"))              sndEnable              = keyFile.get_boolean ("Sounds", "Enable");
    if (keyFile.has_key ("Sounds", "BatchQueueDone")) sndBatchQueueDone = keyFile.get_string ("Sounds", "BatchQueueDone");
    if (keyFile.has_key ("Sounds", "LngEditProcDone"))     sndLngEditProcDone     = keyFile.get_string ("Sounds", "LngEditProcDone");
    if (keyFile.has_key ("Sounds", "LngEditProcDoneSecs")) sndLngEditProcDoneSecs = keyFile.get_double ("Sounds", "LngEditProcDoneSecs");
}

        filterOutParsedExtensions ();

        return 0;
}

int Options::saveToFile (Glib::ustring fname) {

    rtengine::SafeKeyFile keyFile;
    keyFile.set_boolean ("General", "TabbedEditor", tabbedUI);

    keyFile.set_boolean ("General", "StoreLastProfile", savesParamsAtExit);
    if (startupDir==STARTUPDIR_HOME)
        keyFile.set_string ("General", "StartupDirectory", "home");
    else if (startupDir==STARTUPDIR_CURRENT)
        keyFile.set_string ("General", "StartupDirectory", "current");
    else if (startupDir==STARTUPDIR_CUSTOM)
        keyFile.set_string ("General", "StartupDirectory", "custom");
    else if (startupDir==STARTUPDIR_LAST)
        keyFile.set_string ("General", "StartupDirectory", "last");
    keyFile.set_string  ("General", "StartupPath", startupPath);
    keyFile.set_string  ("General", "DateFormat", dateFormat);
    keyFile.set_integer ("General", "AdjusterDelay", adjusterDelay);
    keyFile.set_boolean ("General", "DualProcSupport", rtSettings.dualThreadEnabled);
    keyFile.set_boolean ("General", "MultiUser", multiUser);
    keyFile.set_string  ("General", "Language", language);
    keyFile.set_boolean ("General", "LanguageAutoDetect", languageAutoDetect);
    keyFile.set_string  ("General", "Theme", theme);
    keyFile.set_boolean ("General", "SlimUI", slimUI);
    keyFile.set_boolean ("General", "UseSystemTheme", useSystemTheme);
    keyFile.set_string  ("General", "Version", VERSION);
    keyFile.set_boolean ("General", "FirstRun", false);
    keyFile.set_string  ("General", "DarkFramesPath", rtSettings.darkFramesPath);
    keyFile.set_string  ("General", "FlatFieldsPath", rtSettings.flatFieldsPath);
    keyFile.set_boolean ("General", "Verbose", rtSettings.verbose);

    keyFile.set_integer ("External Editor", "EditorKind", editorToSendTo);
    keyFile.set_string  ("External Editor", "GimpDir", gimpDir);
    keyFile.set_string  ("External Editor", "PhotoshopDir", psDir);
    keyFile.set_string  ("External Editor", "CustomEditor", customEditorProg);
    
    
    keyFile.set_boolean ("File Browser", "BrowseOnlyRaw", fbOnlyRaw);
    keyFile.set_boolean ("File Browser", "BrowserShowsDate", fbShowDateTime);
    keyFile.set_boolean ("File Browser", "BrowserShowsExif", fbShowBasicExif);
    keyFile.set_boolean ("File Browser", "BrowserShowsHidden", fbShowHidden);
    keyFile.set_integer ("File Browser", "ThumbnailSize", thumbSize);
    keyFile.set_integer ("File Browser", "ThumbnailSizeTab", thumbSizeTab);
    keyFile.set_integer ("File Browser", "MaxPreviewHeight", maxThumbnailHeight);
    keyFile.set_integer ("File Browser", "MaxCacheEntries", maxCacheEntries);
    keyFile.set_integer ("File Browser", "ThumbnailFormat", (int)thumbnailFormat);
    Glib::ArrayHandle<Glib::ustring> pext = parseExtensions;
    keyFile.set_string_list ("File Browser", "ParseExtensions", pext);
    Glib::ArrayHandle<int> pextena = parseExtensionsEnabled;
    keyFile.set_integer_list ("File Browser", "ParseExtensionsEnabled", pextena);
    keyFile.set_integer ("File Browser", "ThumbnailArrangement", fbArrangement);
    keyFile.set_integer ("File Browser", "ThumbnailInterpolation", thumbInterp);
    keyFile.set_boolean ("File Browser", "LiveThumbnails", liveThumbnails);
    Glib::ArrayHandle<Glib::ustring> pfav = favoriteDirs;
    keyFile.set_string_list ("File Browser", "FavoriteDirs", pfav);
    Glib::ArrayHandle<Glib::ustring> pren = renameTemplates;
    keyFile.set_string_list ("File Browser", "RenameTemplates", pren);
    keyFile.set_boolean ("File Browser", "RenameUseTemplates", renameUseTemplates);
    Glib::ArrayHandle<double> ptzoom = thumbnailZoomRatios;
    keyFile.set_double_list ("File Browser", "ThumbnailZoomRatios", ptzoom);
    keyFile.set_boolean ("File Browser", "OverlayedFileNames", overlayedFileNames);
    keyFile.set_boolean ("File Browser", "ShowFileNames", showFileNames );
    keyFile.set_boolean ("File Browser", "InternalThumbIfUntouched", internalThumbIfUntouched );
    keyFile.set_boolean ("File Browser", "menuGroupRank", menuGroupRank);
    keyFile.set_boolean ("File Browser", "menuGroupLabel", menuGroupLabel);
    keyFile.set_boolean ("File Browser", "menuGroupFileOperations", menuGroupFileOperations);
    keyFile.set_boolean ("File Browser", "menuGroupProfileOperations", menuGroupProfileOperations);
   
    keyFile.set_integer ("Clipping Indication", "HighlightThreshold", highlightThreshold);
    keyFile.set_integer ("Clipping Indication", "ShadowThreshold", shadowThreshold);
    keyFile.set_boolean ("Clipping Indication", "BlinkClipped", blinkClipped);

    keyFile.set_string  ("Output", "Format", saveFormat.format);
    keyFile.set_integer ("Output", "JpegQuality", saveFormat.jpegQuality);
    keyFile.set_integer ("Output", "PngCompression", saveFormat.pngCompression);
    keyFile.set_integer ("Output", "PngBps", saveFormat.pngBits);
    keyFile.set_integer ("Output", "TiffBps", saveFormat.tiffBits);
    keyFile.set_boolean ("Output", "TiffUncompressed", saveFormat.tiffUncompressed);
    keyFile.set_boolean ("Output", "SaveProcParams", saveFormat.saveParams);

    keyFile.set_string  ("Output", "FormatBatch", saveFormatBatch.format);
    keyFile.set_integer ("Output", "JpegQualityBatch", saveFormatBatch.jpegQuality);
    keyFile.set_integer ("Output", "PngCompressionBatch", saveFormatBatch.pngCompression);
    keyFile.set_integer ("Output", "PngBpsBatch", saveFormatBatch.pngBits);
    keyFile.set_integer ("Output", "TiffBpsBatch", saveFormatBatch.tiffBits);
    keyFile.set_boolean ("Output", "TiffUncompressedBatch", saveFormatBatch.tiffUncompressed);
    keyFile.set_boolean ("Output", "SaveProcParamsBatch", saveFormatBatch.saveParams);

    keyFile.set_string  ("Output", "PathTemplate", savePathTemplate);
    keyFile.set_string  ("Output", "PathFolder", savePathFolder);
    keyFile.set_boolean ("Output", "AutoSuffix", autoSuffix);
    keyFile.set_boolean ("Output", "UsePathTemplate", saveUsePathTemplate);
    keyFile.set_string  ("Output", "LastSaveAsPath", lastSaveAsPath);
	keyFile.set_boolean ("Output", "OverwriteOutputFile", overwriteOutputFile);
    keyFile.set_boolean ("Output", "TunnelMetaData", tunnelMetaData);

    keyFile.set_string  ("Profiles", "Directory", profilePath);
    keyFile.set_string  ("Profiles", "RawDefault", defProfRaw);
    keyFile.set_string  ("Profiles", "ImgDefault", defProfImg);
    keyFile.set_boolean ("Profiles", "SaveParamsWithFile", saveParamsFile);
    keyFile.set_boolean ("Profiles", "SaveParamsToCache", saveParamsCache);
    keyFile.set_integer ("Profiles", "LoadParamsFromLocation", paramsLoadLocation);
    keyFile.set_string  ("Profiles", "CustomProfileBuilder", customProfileBuilder);
    
    keyFile.set_string  ("GUI", "Font", font);
    keyFile.set_integer ("GUI", "WindowWidth", windowWidth);
    keyFile.set_integer ("GUI", "WindowHeight", windowHeight);
    keyFile.set_boolean ("GUI", "WindowMaximized", windowMaximized);
    keyFile.set_integer ("GUI", "DirBrowserWidth", dirBrowserWidth);
    keyFile.set_integer ("GUI", "DirBrowserHeight", dirBrowserHeight);
	keyFile.set_integer ("GUI", "PreferencesWidth", preferencesWidth);
	keyFile.set_integer ("GUI", "PreferencesHeight", preferencesHeight); 
    keyFile.set_integer ("GUI", "SaveAsDialogWidth", saveAsDialogWidth);
    keyFile.set_integer ("GUI", "SaveAsDialogHeight", saveAsDialogHeight);
    keyFile.set_integer ("GUI", "ToolPanelWidth", toolPanelWidth);
    keyFile.set_integer ("GUI", "BrowserToolPanelWidth", browserToolPanelWidth);
    keyFile.set_integer ("GUI", "BrowserToolPanelHeight", browserToolPanelHeight);
    keyFile.set_integer ("GUI", "HistoryPanelWidth", historyPanelWidth);
    keyFile.set_integer ("GUI", "LastPreviewScale", lastScale);
    keyFile.set_integer ("GUI", "LastCropSize", lastCropSize);
    keyFile.set_boolean ("GUI", "ShowHistory", showHistory);
    keyFile.set_integer ("GUI", "ShowFilePanelState", showFilePanelState);
    keyFile.set_boolean ("GUI", "ShowInfo", showInfo);
    keyFile.set_boolean ("GUI", "MainNBVertical", mainNBVertical);
    keyFile.set_boolean ("GUI", "ShowClippedHighlights", showClippedHighlights);
    keyFile.set_boolean ("GUI", "ShowClippedShadows", showClippedShadows);
    keyFile.set_integer ("GUI", "FrameColor", bgcolor);
    keyFile.set_boolean ("GUI", "ProcessingQueueEnbled", procQueueEnabled);
    Glib::ArrayHandle<int> tpopen = tpOpen;
    keyFile.set_integer_list ("GUI", "ToolPanelsExpanded", tpopen);
    keyFile.set_integer ("GUI", "MultiDisplayMode", multiDisplayMode);
    keyFile.set_double_list ("GUI", "CutOverlayBrush", cutOverlayBrush);
    keyFile.set_integer ("GUI", "HistogramPosition", histogramPosition);
    keyFile.set_boolean ("GUI", "ShowProfileSelector", showProfileSelector);
    keyFile.set_boolean ("GUI", "FileBrowserToolbarSingleRow", FileBrowserToolbarSingleRow);

    //Glib::ArrayHandle<int> crvopen = crvOpen;
    //keyFile.set_integer_list ("GUI", "CurvePanelsExpanded", crvopen);

    keyFile.set_integer ("Crop Settings", "PPI", cropPPI);

    keyFile.set_string  ("Color Management", "ICCDirectory",   rtSettings.iccDirectory);
    keyFile.set_string  ("Color Management", "MonitorProfile", rtSettings.monitorProfile);
	keyFile.set_boolean ("Color Management", "AutoMonitorProfile", rtSettings.autoMonitorProfile);
    keyFile.set_integer ("Color Management", "Intent",         rtSettings.colorimetricIntent);
    keyFile.set_boolean ("Color Management", "LCMSSafeMode", rtSettings.LCMSSafeMode);
    keyFile.set_string  ("Color Management", "Adobe_RGB", rtSettings.adobe);
    keyFile.set_string  ("Color Management", "Pro_Photo", rtSettings.prophoto);
    keyFile.set_string  ("Color Management", "Wide_Gamut", rtSettings.widegamut);
    keyFile.set_string  ("Color Management", "S_rgb", rtSettings.srgb);
    keyFile.set_string  ("Color Management", "B_eta", rtSettings.beta);
    keyFile.set_string  ("Color Management", "B_est", rtSettings.best);
    keyFile.set_string  ("Color Management", "B_ruce", rtSettings.bruce);

    Glib::ArrayHandle<int> bab = baBehav;
    keyFile.set_integer_list ("Batch Processing", "AdjusterBehavior", bab);

    keyFile.set_boolean ("Sounds", "Enable", sndEnable);
    keyFile.set_string  ("Sounds", "BatchQueueDone", sndBatchQueueDone);
    keyFile.set_string  ("Sounds", "LngEditProcDone", sndLngEditProcDone);
    keyFile.set_double  ("Sounds", "LngEditProcDoneSecs", sndLngEditProcDoneSecs);


    FILE *f = safe_g_fopen (fname, "wt");
    if (f==NULL)
        return 1;
    else {
        fprintf (f, "%s", keyFile.to_data().c_str());
        fclose (f);
        return 0;
    }
}

Glib::ustring Options::rtdir;
Glib::ustring Options::cacheBaseDir;

void Options::load () {

	// Find the application data path

#ifdef WIN32
	/*
	 * If LOCALAPPDATA exists, RT run on a WinVista/7 system, so we use LOCALAPPDATA as is
	 * otherwise RT run on a Win2000/XP system, so we rebuild the path like this: %USERPROFILE%\Local Settings\Application Data
	 *
	 * Folder redirection is then fully supported on WinVista/7, but not on Win2000/XP
	 */
	const gchar* dataPath;
	Glib::ustring dPath;

	// ->ODUIS: How to make that commented out code work ?

	/*WCHAR path[MAX_PATH] = {0};
	if (SHGetSpecialFolderPathW(NULL, path, CSIDL_LOCAL_APPDATA, false)) {
		dPath = path;
		printf("SHGetSpecialFolderPathW: \"%s\"\n", dPath.c_str());
	}
	else {
		printf("SHGetSpecialFolderPathW: Fail!\n");
	}*/

	dataPath = g_getenv("LOCALAPPDATA");
	if (dataPath != NULL)
		rtdir = Glib::ustring(dataPath) + Glib::ustring("\\") + Glib::ustring(CACHEFOLDERNAME);
	else {
		dataPath = g_getenv("USERPROFILE");
		if (dataPath != NULL)
			rtdir = Glib::ustring(dataPath) + Glib::ustring("\\Local Settings\\Application Data\\") + Glib::ustring(CACHEFOLDERNAME);
	}
#else
    rtdir = Glib::ustring(g_get_user_config_dir ()) + Glib::ustring("/") + Glib::ustring(CACHEFOLDERNAME);
#endif

    // Set the cache folder in RT's base folder
    cacheBaseDir = argv0 + "/cache";

    // Read the global option file (the one located in the application's base folder)
    options.readFromFile (argv0+"/options");

    // Check if RT is installed in Multi-User mode
    if (options.multiUser) {
        // Read the user option file (the one located somewhere in the user's home folder)
    	// Those values supersets those of the global option file
        int r = options.readFromFile (rtdir + "/options");
        // If the local option file does not exist or is broken, and the local cache folder does not exist, recreate it
        if (r && !safe_g_mkdir_with_parents (rtdir, 511)) {
        	// Recreate the user's profile folder
            Glib::ustring profdir = rtdir + "/profiles";
            safe_g_mkdir_with_parents (profdir, 511);
            // Save the option file
            options.saveToFile (rtdir + "/options");
        }
        // Modify the path of the cache folder to the user's personal folder
#ifdef WIN32
        cacheBaseDir = rtdir + "/cache";
#else
        cacheBaseDir = Glib::ustring(g_get_user_cache_dir()) + Glib::ustring("/") + Glib::ustring(CACHEFOLDERNAME);
#endif
    }

	//We handle languages using a hierarchy of translations.  The top of the hierarchy is default.  This includes a default translation for all items
	// (most likely using simple English).  The next level is the language: for instance, English, French, Chinese, etc.  This file should contain a 
	// generic translation for all items which differ from default.  Finally there is the locale.  This is region-specific items which differ from the
	// language file.  These files must be name in the format <Language> (<LC>), where Language is the name of the language which it inherits from,
	// and LC is the local code.  Some examples of this would be English (US) (American English), French (FR) (Franch French), French (CA) (Canadian
	// French), etc.
	//
	// Each level will only contain the differences between itself and its parent translation.  For instance, English (UK) or English (CA) may 
	// include the translation "HISTORY_MSG_34;Avoid Colour Clipping" where English would translate it as "HISTORY_MSG_34;Avoid Color Clipping" (note
	// the difference in the spelling of 'colour').
	//
	// It is important that when naming the translation files, that you stick to the format <Language> or <Language> (<LC>).  We depend on that to figure 
	// out which are the parent translations.  Furthermore, there must be a file <Language> for each locale <Language> (<LC>) -- you cannot have 
	// 'French (CA)' unless there is a file 'French'.

    Glib::ustring defaultTranslation = argv0 + "/languages/default";
	Glib::ustring languageTranslation = "";
	Glib::ustring localeTranslation = "";

    if (options.languageAutoDetect) options.language=langMgr.getOSUserLanguage();

	if (!options.language.empty()){
		std::vector<Glib::ustring> langPortions = Glib::Regex::split_simple(" ", options.language);
		if (langPortions.size() >= 1){
			languageTranslation = argv0 + "/languages/" + langPortions.at(0);
		}
		if (langPortions.size() >= 2){
			localeTranslation = argv0 + "/languages/" + options.language;
		}
	}

	langMgr.load(localeTranslation, new MultiLangMgr(languageTranslation, new MultiLangMgr(defaultTranslation)));

	rtengine::init (&options.rtSettings);
}

void Options::save () {

    if (options.multiUser==false) {
        options.saveToFile (argv0+"/options");
    }
    else {
        options.saveToFile (rtdir + "/options");
    }
}

/*
 * return true if fname ends with one of the retained image file extensions
 */
bool Options::has_retained_extention (Glib::ustring fname) {

	Glib::ustring ext = getExtension(fname).lowercase();

	if (ext.length()) {
		// there is an extension to the filename

		// look out if it has one of the retained extensions
		for (unsigned int i=0; i<parsedExtensions.size(); i++) {
			if (ext == parsedExtensions[i]) {
				return true;
			}
		}
	}
	return false;
}

/*
 * return true if ext is an enabled extension
 */
bool Options::is_extention_enabled (Glib::ustring ext) {
		for (int j=0; j<(int)parseExtensions.size(); j++)
      if (parseExtensions[j].casefold() == ext.casefold())
				return j>=(int)parseExtensionsEnabled.size() || parseExtensionsEnabled[j];
		return false;
}
