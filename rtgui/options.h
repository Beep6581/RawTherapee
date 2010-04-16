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
#include <rtengine.h>

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
    bool firstRun;
    bool savesParamsAtExit;
    SaveFormat saveFormat;
    Glib::ustring savePathTemplate;
    Glib::ustring savePathFolder;
    bool saveUsePathTemplate;
    Glib::ustring defProfRaw;
    Glib::ustring defProfImg;
    Glib::ustring dateFormat;
    int  startupDir;
    Glib::ustring startupPath;
    Glib::ustring profilePath;
    Glib::ustring lastSaveAsPath;
    int toolPanelWidth;
    int browserToolPanelWidth;
    int historyPanelWidth;
    int windowWidth;
    int windowHeight;
    int dirBrowserWidth;
    int dirBrowserHeight;
    int lastScale;
    int lastCropSize;
    bool fbOnlyRaw;
    bool fbShowDateTime;
    bool fbShowBasicExif;
    bool fbShowHidden;
    int  fbArrangement;
    bool multiUser;
    static Glib::ustring rtdir;
    int version;
    int thumbSize;
    bool showHistory;
    int showFilePanelState; // 0: normal, 1: maximized, 2: normal, 3: hidden
    bool showInfo;
    int cropDPI;
    bool showClippedHighlights;
    bool showClippedShadows;
    int highlightThreshold;
    int shadowThreshold;
    bool blinkClipped;
    int bgcolor;
    Glib::ustring language;
    Glib::ustring theme;
    static Glib::ustring cacheBaseDir;
    bool saveParamsFile;
    bool saveParamsCache;
    PPLoadLocation paramsLoadLocation;
    bool procQueueEnabled;
    Glib::ustring gimpDir;
    Glib::ustring psDir;
    Glib::ustring customEditorProg;
    int editorToSendTo;   
    int maxThumbnailHeight;
    int maxCacheEntries;
    ThFileType thumbnailFormat;
    int thumbInterp; // 0: nearest, 1: bilinear
    bool liveThumbnails;
    std::vector<Glib::ustring> parseExtensions;
    std::vector<int> parseExtensionsEnabled;
    std::vector<int> tpOpen;
    std::vector<int> crvOpen;
    std::vector<int> baBehav;
    rtengine::Settings rtSettings;
    
    std::vector<Glib::ustring> favoriteDirs;
    std::vector<Glib::ustring> renameTemplates;
    bool renameUseTemplates;
    
    std::vector<double> thumbnailZoomRatios;
    bool overlayedFileNames;
    
    
                Options         ();

    Options*    copyFrom        (Options* other);
    void        setDefaults     ();
    int         readFromFile    (Glib::ustring fname);
    int         saveToFile      (Glib::ustring fname);
    static void load            ();
    static void save            ();

    bool        is_extention_enabled(Glib::ustring ext);
};

extern Options options;
extern Glib::ustring argv0;
extern Glib::ustring argv1;
extern Glib::ustring versionString;

#endif
