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
#ifndef _IPTCPANEL_
#define _IPTCPANEL_

#include <gtkmm.h>
#include <toolpanel.h>
#include <guiutils.h>
#include <iptcmeta.h>

class XRTWidget{
protected:
	std::string key;
	rtengine::MetadataList &list;
	sigc::connection connChange;

	static void attachToTable( Gtk::Table &table, int row, Gtk::Label &label, Gtk::Widget &widg);
	virtual int readValue(){};
	virtual int writeValue(){};

public:
	XRTWidget(const std::string &k, rtengine::MetadataList &l ):key(k),list(l){}

	// Write gui control value into MetadataList
	void updateList();

	// Read value from MetadataList and update control
	void updateFromList();
};

class XRTLabel: public XRTWidget
{
	Gtk::Label *control;
	int readValue( );
	int writeValue( ){return 0;};

public:
	XRTLabel( Gtk::Table* table, int row, const std::string &key, rtengine::MetadataList &l );

};

class XRTEntry: public XRTWidget
{
	Gtk::Entry *control;
	int readValue( );
	int writeValue( );

public:
	XRTEntry( Gtk::Table* table, int row, const std::string &key, rtengine::MetadataList &l );

};

class XRTEntryMultiline: public XRTWidget
{
	Glib::RefPtr<Gtk::TextBuffer> control;
	int readValue( );
	int writeValue( );

public:
	XRTEntryMultiline( Gtk::Table* table, int row, const std::string &key, rtengine::MetadataList &l );

};

class XRTEntryMultivalue: public XRTWidget
{
    Gtk::Button *addKW;
    Gtk::Button *delKW;
    Gtk::Image* addKWImg;
    Gtk::Image* delKWImg;
    Gtk::ComboBoxEntryText *control;
    Gtk::ListViewText *controlList;

    void add();
    void del();
	int readValue( );
	int writeValue( );

public:
	XRTEntryMultivalue( Gtk::Table* table, int row, const std::string &key, rtengine::MetadataList &l );

};

class XRTCombo: public XRTWidget
{
	Gtk::ComboBoxEntryText *control;
	rtengine::IPTCPairList_t &predefValues;
	int readValue( );
	int writeValue( );

public:
	XRTCombo( Gtk::Table* table, int row, const std::string &key, rtengine::IPTCPairList_t &info, rtengine::MetadataList &l );

};



class IPTCPanel : public Gtk::VBox, public ToolPanel {

    private: 
        rtengine::MetadataList chgList;
        
/*
        Gtk::TextView*  captionView;
        Glib::RefPtr<Gtk::TextBuffer> captionText;
        Gtk::Entry*     captionWriter;
        Gtk::Entry*     headline;
        Gtk::Entry*     instructions;
        Gtk::Entry*     usage;
        Gtk::ComboBoxEntryText* keyword;
        Gtk::ListViewText*  keywords;

        Gtk::ComboBoxEntryText* category;
        Gtk::ComboBoxEntryText* suppCategory;
		Gtk::ComboBoxEntryText* genre;
		Gtk::ComboBoxEntryText* scene;
		Gtk::ComboBoxEntryText* subjectCode;
        Gtk::ListViewText*      suppCategories;
        Gtk::Entry*     eventPhoto;
        Gtk::Entry*     person;
        Gtk::Entry*     modelAge;
        Gtk::Entry*     modelInfo;

        Gtk::Button*    addKW;
        Gtk::Button*    delKW;
        Gtk::Button*    addSC;
        Gtk::Button*    delSC;

        Gtk::Entry*     author;
        Gtk::Entry*     authorPos;
        Gtk::Entry*     creatorCity;
        Gtk::Entry*     creatorCountry;
        Gtk::Entry*     creatorAddr;
        Gtk::Entry*     creatorPCode;
        Gtk::Entry*     creatorRegion;
        Gtk::Entry*     creatorEmail;
        Gtk::Entry*     creatorTel;
        Gtk::Entry*     creatorUrl;
        Gtk::Entry*     credit;
        Gtk::Entry*     source;
        Gtk::Entry*     copyright;
        Gtk::Entry*     rightsOwner;
        Gtk::Entry*     imageCreator;
        Gtk::Entry*     licensor;
        Gtk::Entry*     imageSupplier;
        Gtk::Entry*     guidSupplier;
        Gtk::Entry*     minorDisclosure;
        Gtk::Entry*     modeReleaseID;
        Gtk::ComboBoxEntryText* modeleReleaseStatus;
        Gtk::Entry*     propReleaseID;
        Gtk::ComboBoxEntryText* propReleaseStatus;
        Gtk::Entry*     registryID;
        Gtk::Entry*     registrOrgID;
		Gtk::Entry*     location;
        Gtk::Entry*     city;
        Gtk::Entry*     province;
        Gtk::Entry*     country;
		Gtk::ComboBoxEntryText* countryCode;
		Gtk::Entry*     locationCity;
		Gtk::ComboBoxEntryText* locationCode;
		Gtk::Entry*     locationCountry;
		Gtk::Entry*     locationState;
		Gtk::Entry*     locationSubloc;
		Gtk::ComboBoxEntryText* locationRegion;
		Gtk::Entry*     locCreateCity;
		Gtk::ComboBoxEntryText* locCreateCode;
		Gtk::Entry*     locCreateCountry;
		Gtk::Entry*     locCreateState;
		Gtk::Entry*     locCreateSubloc;
		Gtk::ComboBoxEntryText* locCreateRegion;
		Gtk::Entry*     artworkRights;
		Gtk::Entry*     artworkCreator;
		Gtk::Entry*     artworkDate;
		Gtk::Entry*     artworkSource;
		Gtk::Entry*     artworkNumber;
		Gtk::Entry*     artworkTitle;
        Gtk::Entry*     title;
        Gtk::Label*     dateCreated;
        Gtk::Label*     imageGuid;
        Gtk::Label*     imageHeight;
        Gtk::Label*     imageWidth;
        Gtk::Entry*     transReference;
		Gtk::Entry*     urgency;
*/
        std::vector< XRTWidget* > wdgt;

        Gtk::Button*    reset;
        Gtk::Button*    file;
        Gtk::Button*    copy;
        Gtk::Button*    paste;

        std::vector< sigc::connection > conns;
        
        void applyChangeList ();
        void updateChangeList ();
        
    public:
        IPTCPanel ();
        ~IPTCPanel ();
        
		void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
        
        void setImageData   (const rtengine::ImageMetaData* id);       
		void writeImageData ( rtengine::ImageMetaData* id );

        void notifyListener ();

        void addKeyWord     ();
        void delKeyWord     ();
        void addSuppCategory ();
        void delSuppCategory ();
        
        void resetClicked   ();
        void fileClicked    ();
        void copyClicked    ();
        void pasteClicked   ();
};

#endif
