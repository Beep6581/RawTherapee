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

#include <glibmm.h>
#include "../rtengine/rtengine.h"

#define STARTUPDIR_CURRENT 0
#define STARTUPDIR_HOME    1
#define STARTUPDIR_CUSTOM  2
#define STARTUPDIR_LAST    3

class SaveFormat {

    public:
        Glib::ustring format;
        int pngBits;
        int pngCompression;
        int jpegQuality;
        int tiffBits;
        bool tiffUncompressed;
        bool saveParams;
};

enum ThFileType {FT_Invalid=-1, FT_None=0, FT_Raw=1, FT_Jpeg=2, FT_Tiff=3, FT_Png=4, FT_Custom=5, FT_Tiff16=6, FT_Png16=7, FT_Custom16=8}; 
enum PPLoadLocation {PLL_Cache=0, PLL_Input=1};

class Options {

  private:
    int getString (const char* src, char* dst);
    void error (int line); 

  public:
    bool savesParamsAtExit;
    SaveFormat saveFormat,saveFormatBatch;
    Glib::ustring savePathTemplate;
    Glib::ustring savePathFolder;
    bool saveUsePathTemplate;
    Glib::ustring defProfRaw;
    Glib::ustring defProfImg;
    Glib::ustring dateFormat;
    int adjusterDelay;
    int  startupDir;
    Glib::ustring startupPath;
    Glib::ustring profilePath;
    Glib::ustring loadSaveProfilePath;
    Glib::ustring lastSaveAsPath;
    int saveAsDialogWidth;
    int saveAsDialogHeight;
    int toolPanelWidth;
    int browserToolPanelWidth;
    int browserToolPanelHeight;
    int historyPanelWidth;
    Glib::ustring font;
    int windowWidth;
    int windowHeight;
    bool windowMaximized;
    int dirBrowserWidth;
    int dirBrowserHeight;
	int preferencesWidth;
	int preferencesHeight; 
    int lastScale;
    int panAccelFactor;
    int lastCropSize;
    bool fbOnlyRaw;
    bool fbShowDateTime;
    bool fbShowBasicExif;
    bool fbShowExpComp;
    bool fbShowHidden;
    int  fbArrangement;
    bool multiUser;
    static Glib::ustring rtdir;
    Glib::ustring version;
    int thumbSize,thumbSizeTab;
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
    bool slimUI;
    bool useSystemTheme;
    static Glib::ustring cacheBaseDir;
    bool autoSuffix;
    bool saveParamsFile;
    bool saveParamsCache;
    PPLoadLocation paramsLoadLocation;
    bool procQueueEnabled;
    Glib::ustring gimpDir;
    Glib::ustring psDir;
    Glib::ustring customEditorProg;
    Glib::ustring customProfileBuilder;
    int editorToSendTo;   
    int maxThumbnailHeight;
    int maxCacheEntries;
    ThFileType thumbnailFormat;
    int thumbInterp; // 0: nearest, 1: bilinear
    bool liveThumbnails;
    std::vector<Glib::ustring> parseExtensions;   // List containing all extensions type
    std::vector<int> parseExtensionsEnabled;      // List of bool to retain extension or not
    std::vector<Glib::ustring> parsedExtensions;  // List containing all retained extensions (lowercase)
    std::vector<int> tpOpen;
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
    bool showFileNames;
    bool tabbedUI;
    int previewSizeTab,previewSizeBrowser;
    int multiDisplayMode;  // 0=none, 1=Edit panels on other display
    std::vector<double> cutOverlayBrush;  // Red;Green;Blue;Alpha , all ranging 0..1
    
    Glib::ustring sndBatchQueueDone;
    Glib::ustring sndLngEditProcDone;
    double sndLngEditProcDoneSecs;  // Minimum processing time seconds till the sound is played
    bool sndEnable;

    bool tunnelMetaData;    // Pass through IPTC and XMP unchanged
    int histogramPosition;  // 0=disabled, 1=left pane, 2=right pane
    bool histogramBar;
    bool showProfileSelector;
    bool squareDetailWindow;
    bool FileBrowserToolbarSingleRow;
    bool hideTPVScrollbar;
    bool UseIconNoText;
    int whiteBalanceSpotSize;

    bool menuGroupRank;
    bool menuGroupLabel;
    bool menuGroupFileOperations;
    bool menuGroupProfileOperations;

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
    bool fastexport_bypass_raw_all_enhance;
    bool fastexport_bypass_raw_ccSteps;
    bool fastexport_bypass_raw_dcb_iterations;
    bool fastexport_bypass_raw_dcb_enhance;
    bool fastexport_bypass_raw_ca;
    bool fastexport_bypass_raw_linenoise;
    bool fastexport_bypass_raw_greenthresh;
    bool fastexport_bypass_raw_df;
    bool fastexport_bypass_raw_ff;
    Glib::ustring fastexport_raw_dmethod;
    Glib::ustring fastexport_icm_input;
    Glib::ustring fastexport_icm_working;
    Glib::ustring fastexport_icm_output;
    Glib::ustring fastexport_icm_gamma;
    bool          fastexport_resize_enabled;
    double        fastexport_resize_scale;
    Glib::ustring fastexport_resize_appliesTo;
    Glib::ustring fastexport_resize_method;
    int           fastexport_resize_dataspec;
    int           fastexport_resize_width;
    int           fastexport_resize_height;

    Options ();

    Options*    copyFrom        (Options* other);
    void        filterOutParsedExtensions ();
    void        setDefaults     ();
    int         readFromFile    (Glib::ustring fname);
    int         saveToFile      (Glib::ustring fname);
    static void load            ();
    static void save            ();

    bool        has_retained_extention (Glib::ustring fname);
    bool        is_extention_enabled(Glib::ustring ext);
};

extern Options options;
extern Glib::ustring argv0;
extern Glib::ustring argv1;
extern bool simpleEditor;
extern Glib::ustring versionString;
extern Glib::ustring paramFileExtension;

#endif
