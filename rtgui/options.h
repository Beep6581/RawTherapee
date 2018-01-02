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
#ifndef _OPTIONS_
#define _OPTIONS_

#include <gtkmm.h>
#include "../rtengine/rtengine.h"
#include <exception>

#define STARTUPDIR_CURRENT 0
#define STARTUPDIR_HOME    1
#define STARTUPDIR_CUSTOM  2
#define STARTUPDIR_LAST    3

#define THEMEREGEXSTR      "^(.+)-GTK3-(\\d{1,2})?_(\\d{1,2})?\\.css$"

// Default bundled profile name to use for Raw images
#ifdef WIN32
#define DEFPROFILE_RAW      "${G}\\Default"
#else
#define DEFPROFILE_RAW      "${G}/Default"
#endif
// Default bundled profile name to use for Standard images
#define DEFPROFILE_IMG      "Neutral"
// Profile name to use for internal values' profile
#define DEFPROFILE_INTERNAL "Neutral"
// Special name for the Dynamic profile
#define DEFPROFILE_DYNAMIC  "Dynamic"

struct SaveFormat {
    SaveFormat() :
        format ("jpg"),
        pngBits (8),
        jpegQuality (90),
        jpegSubSamp (2),
        tiffBits (8),
        tiffUncompressed (true),
        saveParams (true)
    {
    }

    Glib::ustring format;
    int pngBits;
    int jpegQuality;
    int jpegSubSamp;  // 1=best compression, 3=best quality
    int tiffBits;
    bool tiffUncompressed;
    bool saveParams;
};

enum ThFileType {FT_Invalid = -1, FT_None = 0, FT_Raw = 1, FT_Jpeg = 2, FT_Tiff = 3, FT_Png = 4, FT_Custom = 5, FT_Tiff16 = 6, FT_Png16 = 7, FT_Custom16 = 8};
enum PPLoadLocation {PLL_Cache = 0, PLL_Input = 1};
enum CPBKeyType {CPBKT_TID = 0, CPBKT_NAME = 1, CPBKT_TID_NAME = 2};
enum prevdemo_t {PD_Sidecar = 1, PD_Fast = 0};
enum mip_t {MI_prev = 0, MI_opt = 1};
enum locaaju_t {lo_std = 0, lo_enh = 1, lo_enhde = 2};

class Options
{
public:
    class Error: public std::exception
    {
    public:
        Error (const Glib::ustring &msg): msg_ (msg) {}
        const char *what() const throw()
        {
            return msg_.c_str();
        }
        const Glib::ustring &get_msg() const throw()
        {
            return msg_;
        }

    private:
        Glib::ustring msg_;
    };

private:
    bool defProfRawMissing;
    bool defProfImgMissing;
    Glib::ustring userProfilePath;
    Glib::ustring globalProfilePath;
    bool checkProfilePath (Glib::ustring &path);
    bool checkDirPath (Glib::ustring &path, Glib::ustring errString);
    void updatePaths();
    int getString (const char* src, char* dst);
    void error (int line);
    /**
     * Safely reads a directory from the configuration file and only applies it
     * to the provided destination variable if there is a non-empty string in
     * the configuration.
     *
     * @param keyFile file to read configuration from
     * @param section name of the section in the configuration file
     * @param entryName name of the entry in the configuration file
     * @param destination destination variable to store to
     * @return @c true if @p destination was changed
     */
    bool safeDirGet (const Glib::KeyFile& keyFile, const Glib::ustring& section,
                     const Glib::ustring& entryName, Glib::ustring& destination);

public:

    enum class NavigatorUnit {
        PERCENT,
        R0_255,
        R0_1,
        _COUNT
    };
    bool savesParamsAtExit;
    SaveFormat saveFormat, saveFormatBatch;
    Glib::ustring savePathTemplate;
    Glib::ustring savePathFolder;
    bool saveUsePathTemplate;
    Glib::ustring defProfRaw;
    Glib::ustring defProfImg;
    Glib::ustring dateFormat;
    int adjusterMinDelay;
    int adjusterMaxDelay;
    int  startupDir;
    Gtk::SortType dirBrowserSortType;
    Glib::ustring startupPath;
    Glib::ustring profilePath; // can be an absolute or relative path; depending on this value, bundled profiles may not be found
    bool useBundledProfiles;   // only used if multiUser == true
    Glib::ustring loadSaveProfilePath;
    Glib::ustring lastSaveAsPath;
    int saveAsDialogWidth;
    int saveAsDialogHeight;
    int toolPanelWidth;
    int browserToolPanelWidth;
    int browserToolPanelHeight;
    bool browserToolPanelOpened;
    bool browserDirPanelOpened;
    bool editorFilmStripOpened;
    int historyPanelWidth;
    int windowX;
    int windowY;
    int windowWidth;
    int windowHeight;
    bool windowMaximized;
    int windowMonitor;
    int meowMonitor;
    bool meowFullScreen;
    bool meowMaximized;
    int meowWidth;
    int meowHeight;
    int meowX;
    int meowY;
    int detailWindowWidth;
    int detailWindowHeight;
    int dirBrowserWidth;
    int dirBrowserHeight;
    int preferencesWidth;
    int preferencesHeight;
    bool lastShowAllExif;
    int lastScale;
    int panAccelFactor;
    int lastCropSize;
    Glib::ustring fontFamily;    // RT's main font family
    int fontSize;                // RT's main font size (units: pt)
    Glib::ustring CPFontFamily;  // ColorPicker font family
    int CPFontSize;              // ColorPicker font size (units: pt)
    bool fbOnlyRaw;
    bool fbShowDateTime;
    bool fbShowBasicExif;
    bool fbShowExpComp;
    bool fbShowHidden;
    int  fbArrangement;
    NavigatorUnit navRGBUnit;
    NavigatorUnit navHSVUnit;
    bool multiUser;
    static Glib::ustring rtdir;
    Glib::ustring version;
    int thumbSize, thumbSizeTab, thumbSizeQueue;
    bool sameThumbSize;     // Will use only one thumb size for the file browser and the single editor tab, and avoid recomputing them
    bool showHistory;
    int showFilePanelState; // 0: normal, 1: maximized, 2: normal, 3: hidden
    bool showInfo;
    bool mainNBVertical;  // main notebook vertical tabs?
    int cropPPI;
    bool showClippedHighlights;
    bool showClippedShadows;
    int highlightThreshold;
    int shadowThreshold;
    bool blinkClipped;
    int bgcolor;
    Glib::ustring language;
    bool languageAutoDetect;
    Glib::ustring theme;
    static Glib::ustring cacheBaseDir;
    bool autoSuffix;
    bool forceFormatOpts;
    int saveMethodNum;
    bool saveParamsFile;
    bool saveParamsCache;
    PPLoadLocation paramsLoadLocation;
    bool procQueueEnabled;
    Glib::ustring gimpDir;
    Glib::ustring psDir;
    Glib::ustring customEditorProg;
    Glib::ustring CPBPath; // Custom Profile Builder's path
    CPBKeyType CPBKeys; // Custom Profile Builder's key type
    int editorToSendTo;
    int maxThumbnailHeight;
    std::size_t maxCacheEntries;
    int thumbInterp; // 0: nearest, 1: bilinear
    std::vector<Glib::ustring> parseExtensions;   // List containing all extensions type
    std::vector<int> parseExtensionsEnabled;      // List of bool to retain extension or not
    std::vector<Glib::ustring> parsedExtensions;  // List containing all retained extensions (lowercase)
    std::vector<int> tpOpen;
    bool autoSaveTpOpen;
    //std::vector<int> crvOpen;
    std::vector<int> baBehav;
    rtengine::Settings rtSettings;

    std::vector<Glib::ustring> favoriteDirs;
    std::vector<Glib::ustring> renameTemplates;
    bool renameUseTemplates;
    bool internalThumbIfUntouched;
    bool overwriteOutputFile;

    std::vector<double> thumbnailZoomRatios;
    bool overlayedFileNames;
    bool filmStripOverlayedFileNames;
    bool showFileNames;
    bool filmStripShowFileNames;
    bool tabbedUI;
    bool rememberZoomAndPan;
    int multiDisplayMode;  // 0=none, 1=Edit panels on other display
    std::vector<double> cutOverlayBrush;  // Red;Green;Blue;Alpha , all ranging 0..1
    std::vector<double> navGuideBrush;  // Red;Green;Blue;Alpha , all ranging 0..1

    Glib::ustring sndBatchQueueDone;
    Glib::ustring sndLngEditProcDone;
    double sndLngEditProcDoneSecs;  // Minimum processing time seconds till the sound is played
    bool sndEnable;

    int histogramPosition;  // 0=disabled, 1=left pane, 2=right pane
    //int histogramWorking;  // 0=disabled, 1=left pane, 2=right pane
    bool histogramBar;
    bool histogramFullMode;
    bool FileBrowserToolbarSingleRow;
    bool hideTPVScrollbar;
    bool UseIconNoText;
    int whiteBalanceSpotSize;
    int curvebboxpos; // 0=above, 1=right, 2=below, 3=left

    bool showFilmStripToolBar;
    bool showdelimspot;

    // Performance options
    Glib::ustring clutsDir;
    int rgbDenoiseThreadLimit; // maximum number of threads for the denoising tool ; 0 = use the maximum available
    int maxInspectorBuffers;   // maximum number of buffers (i.e. images) for the Inspector feature
    int clutCacheSize;
    bool filledProfile;  // Used as reminder for the ProfilePanel "mode"
    prevdemo_t prevdemo; // Demosaicing method used for the <100% preview
    bool serializeTiffRead;
    mip_t mip; // MIP
    locaaju_t locaaju;

    bool menuGroupRank;
    bool menuGroupLabel;
    bool menuGroupFileOperations;
    bool menuGroupProfileOperations;
    bool menuGroupExtProg;

    // fast export options
    bool fastexport_bypass_sharpening;
    bool fastexport_bypass_sharpenEdge;
    bool fastexport_bypass_sharpenMicro;
    //bool fastexport_bypass_lumaDenoise;
    //bool fastexport_bypass_colorDenoise;
    bool fastexport_bypass_defringe;
    bool fastexport_bypass_dirpyrDenoise;
    bool fastexport_bypass_sh_hq;
    bool fastexport_bypass_dirpyrequalizer;
    bool fastexport_bypass_wavelet;
    Glib::ustring fastexport_raw_bayer_method;
    //bool fastexport_bypass_raw_bayer_all_enhance;
    bool fastexport_bypass_raw_bayer_dcb_iterations;
    bool fastexport_bypass_raw_bayer_dcb_enhance;
    bool fastexport_bypass_raw_bayer_lmmse_iterations;
    bool fastexport_bypass_raw_bayer_linenoise;
    bool fastexport_bypass_raw_bayer_greenthresh;
    Glib::ustring fastexport_raw_xtrans_method;
    bool fastexport_bypass_raw_ccSteps;
    bool fastexport_bypass_raw_ca;
    bool fastexport_bypass_raw_df;
    bool fastexport_bypass_raw_ff;
    Glib::ustring fastexport_icm_input;
    Glib::ustring fastexport_icm_working;
    Glib::ustring fastexport_icm_output;
    rtengine::RenderingIntent fastexport_icm_outputIntent;
    bool          fastexport_icm_outputBPC;
    Glib::ustring fastexport_icm_gamma;
    bool          fastexport_resize_enabled;
    double        fastexport_resize_scale;
    Glib::ustring fastexport_resize_appliesTo;
    Glib::ustring fastexport_resize_method;
    int           fastexport_resize_dataspec;
    int           fastexport_resize_width;
    int           fastexport_resize_height;
    bool fastexport_use_fast_pipeline;

    // Dialog settings
    Glib::ustring lastIccDir;
    Glib::ustring lastDarkframeDir;
    Glib::ustring lastFlatfieldDir;
    Glib::ustring lastRgbCurvesDir;
    Glib::ustring lastLabCurvesDir;
    Glib::ustring lastRetinexDir;
    Glib::ustring lastDenoiseCurvesDir;
    Glib::ustring lastWaveletCurvesDir;
    Glib::ustring lastlocalCurvesDir;
    Glib::ustring lastPFCurvesDir;
    Glib::ustring lastHsvCurvesDir;
    Glib::ustring lastToneCurvesDir;
    Glib::ustring lastColorToningCurvesDir;
    Glib::ustring lastVibranceCurvesDir;
    Glib::ustring lastProfilingReferenceDir;
    Glib::ustring lastBWCurvesDir;
    Glib::ustring lastLensProfileDir;
    bool gimpPluginShowInfoDialog;

    size_t maxRecentFolders;                   // max. number of recent folders stored in options file
    std::vector<Glib::ustring> recentFolders;  // List containing all recent folders


    Options ();

    Options*    copyFrom        (Options* other);
    void        filterOutParsedExtensions ();
    void        setDefaults     ();
    void readFromFile (Glib::ustring fname);
    void saveToFile (Glib::ustring fname);
    static void load (bool lightweight = false);
    static void save();

    // if multiUser=false, send back the global profile path
    Glib::ustring getPreferredProfilePath();
    Glib::ustring getUserProfilePath()
    {
        return userProfilePath;
    }
    Glib::ustring getGlobalProfilePath()
    {
        return globalProfilePath;
    }
    Glib::ustring findProfilePath (Glib::ustring &profName);
    bool        is_parse_extention (Glib::ustring fname);
    bool        has_retained_extention (Glib::ustring fname);
    bool        is_extention_enabled (Glib::ustring ext);
    bool        is_defProfRawMissing()
    {
        return defProfRawMissing;
    }
    bool        is_defProfImgMissing()
    {
        return defProfImgMissing;
    }
    void        setDefProfRawMissing (bool value)
    {
        defProfRawMissing = value;
    }
    void        setDefProfImgMissing (bool value)
    {
        defProfImgMissing = value;
    }
};

extern Options options;
extern Glib::ustring argv0;
extern Glib::ustring argv1;
extern bool simpleEditor;
extern bool gimpPlugin;
extern bool remote;
extern Glib::ustring versionString;
extern Glib::ustring paramFileExtension;

#endif
