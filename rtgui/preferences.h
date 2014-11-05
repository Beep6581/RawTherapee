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
#ifndef __PREFERENCES_H__
#define __PREFERENCES_H__

#include <gtkmm.h>
#include "adjuster.h"
#include "options.h"
#include <vector>
#include "rtwindow.h"

class Preferences : public Gtk::Dialog, public ProfileStoreListener {

        class ExtensionColumns : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<bool>  enabled;
                Gtk::TreeModelColumn<Glib::ustring>  ext;
                ExtensionColumns() { add(enabled); add(ext); }
        };
        ExtensionColumns extensionColumns;
        Glib::RefPtr<Gtk::ListStore> extensionModel;


        class BehavColumns : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<Glib::ustring> label;
                Gtk::TreeModelColumn<bool>          badd;
                Gtk::TreeModelColumn<bool>          bset;
                Gtk::TreeModelColumn<bool>          visible;
                Gtk::TreeModelColumn<int>           addsetid;
                BehavColumns() { add(label); add(badd); add(bset); add(visible); add(addsetid); }
        };
        Glib::RefPtr<Gtk::TreeStore> behModel;
        BehavColumns behavColumns;


  protected:
    Splash* splash;
    ProfileStoreComboBox* rprofiles;
    Gtk::TreeIter currRawRow; // :)
    ProfileStoreComboBox* iprofiles;
    Gtk::TreeIter currImgRow;
    Gtk::ComboBoxText* languages;
    Gtk::CheckButton* ckbLangAutoDetect;
    Gtk::Entry* dateformat;
    Gtk::Entry* startupdir;
    Gtk::RadioButton* sdcurrent;
    Gtk::RadioButton* sdlast;
    Gtk::RadioButton* sdhome;
    Gtk::RadioButton* sdother;
    Gtk::FileChooserButton* gimpDir;
    Gtk::FileChooserButton* psDir;
    Gtk::Entry* editorToSendTo;
    Gtk::RadioButton* edGimp;
    Gtk::RadioButton* edPS;
    Gtk::RadioButton* edOther;
    Gtk::FileChooserButton* darkFrameDir;
    Gtk::FileChooserButton* flatFieldDir;
    Gtk::FileChooserButton* clutsDir;
    Gtk::Label *dfLabel;
    Gtk::Label *ffLabel;

    Gtk::CheckButton* showDateTime;
    Gtk::CheckButton* showBasicExif;
    Gtk::CheckButton* showExpComp;

    Gtk::FileChooserButton* iccDir;
    Gtk::FileChooserButton* monProfile;
	Gtk::CheckButton* cbAutoMonProfile;
	//Gtk::CheckButton* cbAutocielab;
	Gtk::CheckButton* cbciecamfloat;

	Gtk::SpinButton*  hlThresh;
	Gtk::SpinButton*  shThresh;

	Gtk::SpinButton*  panFactor;

    Gtk::ComboBoxText* intent;
    Gtk::ComboBoxText* view;
    Gtk::ComboBoxText* grey;
    Gtk::ComboBoxText* greySc;

    Gtk::ComboBoxText* theme;
    Gtk::CheckButton* slimUI;
    Gtk::HBox* hbtheme;
    Gtk::CheckButton* chUseSystemTheme;
    Gtk::FontButton* fontbutton;
    Gtk::ColorButton* butCropCol;
    Gtk::ColorButton* butNavGuideCol;
	
    Gtk::SpinButton*   maxThumbSize;
    Gtk::SpinButton*   maxCacheEntries;
    Gtk::Button*       clearThumbnails;
    Gtk::Button*       clearProfiles;
    Gtk::Button*       clearAll;
    Gtk::Entry*     extension;
    Gtk::TreeView*  extensions;
    Gtk::Button*    addExt;
    Gtk::Button*    delExt;
    Gtk::CheckButton* overlayedFileNames;
    Gtk::CheckButton* filmStripOverlayedFileNames;
    Gtk::CheckButton* sameThumbSize;

    Gtk::SpinButton*  rgbDenoiseTreadLimitSB;

    Gtk::CheckButton* ckbmenuGroupRank;
    Gtk::CheckButton* ckbmenuGroupLabel;
    Gtk::CheckButton* ckbmenuGroupFileOperations;
    Gtk::CheckButton* ckbmenuGroupProfileOperations;
    Gtk::CheckButton* ckbmenuGroupExtProg;

    Gtk::Button*      behAddAll;
    Gtk::Button*      behSetAll;
    Gtk::CheckButton* chOverwriteOutputFile;

    Gtk::CheckButton* saveParamsFile;
    Gtk::CheckButton* saveParamsCache;
    Gtk::CheckButton* useBundledProfiles;
    Gtk::ComboBoxText* loadParamsPreference;
    Gtk::ComboBoxText* editorLayout;
    RTWindow* parent;
	
    Gtk::CheckButton* ckbSndEnable;
    Gtk::Entry* txtSndBatchQueueDone;
    Gtk::Entry* txtSndLngEditProcDone;
	Gtk::SpinButton* spbSndLngEditProcDoneSecs;
	
    Gtk::CheckButton* ckbTunnelMetaData;
	Gtk::CheckButton* ckbInternalThumbIfUntouched;

    Gtk::Entry* txtCustProfBuilderPath;
    Gtk::ComboBoxText* custProfBuilderLabelType;

    Gtk::CheckButton* ckbHistogramPositionLeft;
    Gtk::CheckButton* ckbHistogramWorking;
    Gtk::CheckButton* ckbShowProfileSelector;
    Gtk::CheckButton* ckbFileBrowserToolbarSingleRow;
    Gtk::CheckButton* ckbHideTPVScrollbar;
    Gtk::CheckButton* ckbUseIconNoText;

    Glib::ustring storedValueRaw;
    Glib::ustring storedValueImg;

    Options moptions;
    sigc::connection tconn, sconn, fconn, usethcon, addc, setc, dfconn, ffconn, bpconn, rpconn, ipconn;
    sigc::connection autoMonProfileConn, sndEnableConn, langAutoDetectConn, autocielabConn;
    Glib::ustring initialTheme;
    Glib::ustring initialFont;

    bool oldSlimUI;

    void fillPreferences ();
    void storePreferences ();
    void parseDir       (Glib::ustring dirname, std::vector<Glib::ustring>& items, Glib::ustring ext);
    void updateDFinfos ();
    void updateFFinfos ();
    void workflowUpdate();
    void themeChanged  ();
    void useThemeChanged();
    void fontChanged   ();
    void forRAWComboChanged ();
    void forImageComboChanged ();
    void layoutComboChanged ();
    void bundledProfilesChanged();
    void switchThemeTo (Glib::ustring newTheme, bool slimInterface);
    void switchFontTo  (Glib::ustring newFont);
    bool splashClosed(GdkEventAny* event);

    void appendBehavList (Gtk::TreeModel::iterator& parent, Glib::ustring label, int id, bool set);

    Gtk::Widget* getProcParamsPanel ();
    Gtk::Widget* getColorManagementPanel ();
    Gtk::Widget* getFileBrowserPanel ();
    Gtk::Widget* getGeneralPanel ();
    Gtk::Widget* getBatchProcPanel ();
    Gtk::Widget* getPerformancePanel ();
    Gtk::Widget* getSoundPanel ();
    
  public:
	Preferences (RTWindow *rtwindow);
	~Preferences ();
    
    void savePressed ();
    void loadPressed ();
    void okPressed ();
    void cancelPressed ();
    void aboutPressed ();
	void autoMonProfileToggled ();
    void sndEnableToggled ();
    void langAutoDetectToggled ();
    void autocielabToggled ();

    void selectStartupDir ();
    void addExtPressed ();
    void delExtPressed ();
    void darkFrameChanged ();
    void flatFieldChanged ();
    void clearProfilesPressed ();
    void clearThumbImagesPressed ();
    void clearAllPressed ();

    void behAddRadioToggled (const Glib::ustring& path);
    void behSetRadioToggled (const Glib::ustring& path);
    void behAddAllPressed ();
    void behSetAllPressed ();

    virtual void storeCurrentValue();
    virtual void updateProfileList();
    virtual void restoreValue();

//    void selectICCProfileDir ();
//    void selectMonitorProfile ();
};

#endif
