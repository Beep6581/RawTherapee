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
#include <iptcpanel.h>
#include <clipboard.h>
#include <rtimage.h>
#include <partialpastedlg.h>

extern Glib::ustring argv0;

using namespace rtengine;
using namespace rtengine::procparams;

IPTCPanel::IPTCPanel () {

    set_border_width (2);

    Gtk::Table* iptc = Gtk::manage( new Gtk::Table (72, 2) );
    
    int row = 0;

    Gtk::Label* H0 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_DESCRIPTION")) );
    H0->set_use_markup (true);
    H0->set_alignment(Gtk::ALIGN_LEFT);
    iptc->attach (*H0, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCHeadline, chgList ) );
    wdgt.push_back( new XRTEntryMultiline( iptc, row++, kIPTCDescription, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCWriter, chgList ) );

    Gtk::VBox* vbox1 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H1 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_LOCATIONSHOWNCORE")) );
    H1->set_use_markup (true);
    H1->set_alignment(Gtk::ALIGN_LEFT);
    vbox1->pack_start(*hsept, Gtk::PACK_SHRINK, 4);
    vbox1->pack_start(*H1);
    iptc->attach (*vbox1, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocation, chgList ) ); //Iptc4xmpCore:Location
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCity, chgList ) );     //photoshop:City
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCState, chgList ) );    //photoshop:State
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCountry, chgList ) );  //photoshop:Country
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCCountryCode, rtengine::IPTCMeta::IPTCISO3166, chgList) ); //Iptc4xmpCore:CountryCode

    Gtk::VBox* vbox2 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept2 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H2 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_LOCATIONSHOWNEXT")) );
    H2->set_use_markup (true);
    H2->set_alignment(Gtk::ALIGN_LEFT);
    vbox2->pack_start(*hsept2, Gtk::PACK_SHRINK, 4);
    vbox2->pack_start(*H2);
    iptc->attach (*vbox2, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    //Location based on Iptc4xmpExt
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocationSubloc, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocationCity, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocationState, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocationCtry, chgList ) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCLocationCode, rtengine::IPTCMeta::IPTCISO3166, chgList) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCLocationRegion, rtengine::IPTCMeta::IPTCWorldRegion, chgList) );

    Gtk::VBox* vbox3 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept3 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H3 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_LOCATIONCREATED")) );
    H3->set_use_markup (true);
    H3->set_alignment(Gtk::ALIGN_LEFT);
    vbox3->pack_start(*hsept3, Gtk::PACK_SHRINK, 4);
    vbox3->pack_start(*H3);
    iptc->attach (*vbox3, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocCreateSubloc, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocCreateCity, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocCreateState, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLocCreateCtry, chgList ) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCLocCreateCode, rtengine::IPTCMeta::IPTCISO3166, chgList) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCLocCreateRegion, rtengine::IPTCMeta::IPTCWorldRegion, chgList) );

    Gtk::VBox* vbox4 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept4 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H4 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_ARTWORKOROBJECT")) );
    H4->set_use_markup (true);
    H4->set_alignment(Gtk::ALIGN_LEFT);
    vbox4->pack_start(*hsept4, Gtk::PACK_SHRINK, 4);
    vbox4->pack_start(*H4);
    iptc->attach (*vbox4, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCArtworkTitle, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCArtworkCreator, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCArtworkRights, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCArtworkDate, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCArtworkSource, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCArtworkNumber, chgList ) );

    Gtk::VBox* vbox5 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept5 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H5 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_CLASSIFICATION")) );
    H5->set_use_markup (true);
    H5->set_alignment(Gtk::ALIGN_LEFT);
    vbox5->pack_start(*hsept5, Gtk::PACK_SHRINK, 4);
    vbox5->pack_start(*H5);
    iptc->attach (*vbox5, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntryMultivalue( iptc, row++, kIPTCKeywords, chgList ) );row++;
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCCategory, rtengine::IPTCMeta::IPTCSubject, chgList) );
    wdgt.push_back( new XRTEntryMultivalue( iptc, row++, kIPTCSuppCateg, chgList ) );row++;
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCGenre, rtengine::IPTCMeta::IPTCGenre, chgList) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCScene, rtengine::IPTCMeta::IPTCScene, chgList) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCSubjCode, rtengine::IPTCMeta::IPTCSubject, chgList) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCEvent, chgList ) );

    Gtk::VBox* vbox6 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept6 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H6 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_ABOUTMODELS")) );
    H6->set_use_markup (true);
    H6->set_alignment(Gtk::ALIGN_LEFT);
    vbox6->pack_start(*hsept6, Gtk::PACK_SHRINK, 4);
    vbox6->pack_start(*H6);
    iptc->attach (*vbox6, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCPerson, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCModelAge, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCModelInfo, chgList ) );

    Gtk::VBox* vbox7 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept7 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H7 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_CREATORINFO")) );
    H7->set_use_markup (true);
    H7->set_alignment(Gtk::ALIGN_LEFT);
    vbox7->pack_start(*hsept7, Gtk::PACK_SHRINK, 4);
    vbox7->pack_start(*H7);
    iptc->attach (*vbox7, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreator, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCAuthorPos, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreatorExtadr, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreatorPcode, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreatorAdrCity, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreatorRegion, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreatorAdrCtry, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreatorEmail, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreatorTel, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCreatorUrl, chgList ) );

    Gtk::VBox* vbox8 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept8 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H8 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_RIGHTSINFO")) );
    H8->set_use_markup (true);
    H8->set_alignment(Gtk::ALIGN_LEFT);
    vbox8->pack_start(*hsept8, Gtk::PACK_SHRINK, 4);
    vbox8->pack_start(*H8);
    iptc->attach (*vbox8, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCCredit, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCSource, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCRights, chgList ) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCRightsStatus, rtengine::IPTCMeta::IPTCCopyrightStatus, chgList ) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCRightsMarked, rtengine::IPTCMeta::IPTCRightsMarked, chgList) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCRightsCertif, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCRightsStatement, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCInstruct, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCUsageTerms, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCRightsOwner, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCImageCreator, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCLicensor, chgList ) );

    Gtk::VBox* vbox9 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept9 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H9 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_RELEASEINFO")) );
    H9->set_use_markup (true);
    H9->set_alignment(Gtk::ALIGN_LEFT);
    vbox9->pack_start(*hsept9, Gtk::PACK_SHRINK, 4);
    vbox9->pack_start(*H9);
    iptc->attach (*vbox9, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCMinorModelAgeDisclosure, rtengine::IPTCMeta::IPTCMinorModelAgeDisclosure, chgList) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCModelReleaseID, chgList ) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCModelReleaseSt, rtengine::IPTCMeta::IPTCModelReleaseStatus, chgList) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCPropertyRelID, chgList ) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCPropertyRelSt, rtengine::IPTCMeta::IPTCPropertyReleaseStatus, chgList) );


    Gtk::VBox* vbox13 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept13 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H13 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_FEATUREDORGANISATION")) );
    H13->set_use_markup (true);
    H13->set_alignment(Gtk::ALIGN_LEFT);
    vbox13->pack_start(*hsept13, Gtk::PACK_SHRINK, 4);
    vbox13->pack_start(*H13);
    iptc->attach (*vbox13, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCOrganisationName, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCOrganisationCode, chgList ) );


    Gtk::VBox* vbox10 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept10 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H10 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_REFERENCEADMINISTRATION")) );
    H10->set_use_markup (true);
    H10->set_alignment(Gtk::ALIGN_LEFT);
    vbox10->pack_start(*hsept10, Gtk::PACK_SHRINK, 4);
    vbox10->pack_start(*H10);
    iptc->attach (*vbox10, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCTitle, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCReference, chgList ) );
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCUrgency, rtengine::IPTCMeta::IPTCUrgency, chgList) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCRegistryID, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCRegistryOrgID, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCImageSupplier, chgList ) );
    wdgt.push_back( new XRTEntry( iptc, row++, kIPTCGUIDSupplier, chgList ) );
    //TODO!!! kIPTCDigitalSourceType contains long values, may be display it in 2 rows?
    wdgt.push_back( new XRTCombo( iptc, row++, kIPTCDigitalSourceType, rtengine::IPTCMeta::IPTCDigitalSourceType, chgList) );

    Gtk::VBox* vbox11 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept11 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* H11 = Gtk::manage( new Gtk::Label (M("IPTC_GROUPNAME_READONLYINFORMATION")) );
    H11->set_use_markup (true);
    H11->set_alignment(Gtk::ALIGN_LEFT);
    vbox11->pack_start(*hsept11, Gtk::PACK_SHRINK, 4);
    vbox11->pack_start(*H11);
    iptc->attach (*vbox11, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
    row=row+1;

    wdgt.push_back( new XRTLabel( iptc, row++, kIPTCDate, chgList ) );
    wdgt.push_back( new XRTLabel( iptc, row++, kIPTCGUID, chgList ) );
    wdgt.push_back( new XRTLabel( iptc, row++, kIPTCMaxHeight, chgList ) );
    wdgt.push_back( new XRTLabel( iptc, row++, kIPTCMaxWidth, chgList ) );

    // add panel ending
    Gtk::VBox* vbox12 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept12 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Image* peImg = Gtk::manage (new RTImage("PanelEnding.png"));
    vbox12->pack_start(*hsept12, Gtk::PACK_SHRINK, 4);
    vbox12->pack_start(*peImg);
    iptc->attach (*vbox12, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);

    Gtk::ScrolledWindow* scrolledWindow = Gtk::manage( new Gtk::ScrolledWindow() );
    scrolledWindow->set_border_width(2);
    scrolledWindow->set_shadow_type(Gtk::SHADOW_NONE);
    scrolledWindow->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    scrolledWindow->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);
    scrolledWindow->add(*iptc);
    
    pack_start (*scrolledWindow);
    
    Gtk::HBox* bbox = Gtk::manage( new Gtk::HBox () );
    
    reset = Gtk::manage( new Gtk::Button (M("IPTCPANEL_RESET")) );
    reset->set_image (*Gtk::manage(new RTImage ("gtk-undo-ltr.png", "gtk-undo-rtl.png")));
    bbox->pack_start (*reset, Gtk::PACK_SHRINK, 0);

    fileOpen = Gtk::manage( new Gtk::Button () );
    fileOpen->set_image (*Gtk::manage(new RTImage ("gtk-open.png")));

    fileSave = Gtk::manage( new Gtk::Button () );
    fileSave->set_image (*Gtk::manage(new RTImage ("gtk-save-large.png")));
    
    copy = Gtk::manage( new Gtk::Button () );
    copy->set_image (*Gtk::manage(new RTImage ("edit-copy.png")));

    paste = Gtk::manage( new Gtk::Button () );
    paste->set_image (*Gtk::manage(new RTImage ("edit-paste.png")));
    bbox->pack_end (*paste, Gtk::PACK_SHRINK, 0);
    bbox->pack_end (*copy, Gtk::PACK_SHRINK, 0);
    bbox->pack_end (*fileSave, Gtk::PACK_SHRINK, 0);
    bbox->pack_end (*fileOpen, Gtk::PACK_SHRINK, 0);
    
    pack_end (*bbox, Gtk::PACK_SHRINK, 2);

    Gtk::Tooltips* toolTip = Gtk::manage( new Gtk::Tooltips () );
    toolTip->set_tip (*reset, M("IPTCPANEL_RESETHINT"));
    toolTip->set_tip (*fileOpen, M("IPTCPANEL_OPENHINT"));
    toolTip->set_tip (*fileSave, M("IPTCPANEL_SAVEHINT"));
    toolTip->set_tip (*copy, M("IPTCPANEL_COPYHINT"));
    toolTip->set_tip (*paste, M("IPTCPANEL_PASTEHINT"));
    
    reset->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::resetClicked) );
    fileOpen->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::fileOpenClicked) );
    fileSave->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::fileSaveClicked) );
    copy->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::copyClicked) );
    paste->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::pasteClicked) );

    show_all ();
}

IPTCPanel::~IPTCPanel ()
{
	for( unsigned i=0; i< wdgt.size(); i++)
		delete wdgt[i];
}


void IPTCPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    //defChangeList = defParams->iptc;
}
        
void IPTCPanel::setImageData (const ImageMetaData* id) {
    
	if( id ){
		idata = id;
		chgList = idata->getIPTCData();
		applyChangeList();
	}
}

void IPTCPanel::writeImageData ( rtengine::ImageMetaData* id )
{
	updateChangeList ();
	if( id )
		id->setIPTCData( chgList );
}

void IPTCPanel::notifyListener () {
    
    if (listener)
        listener->panelChanged (EvIPTC, M("HISTORY_CHANGED"));
}

void IPTCPanel::applyChangeList () {

	for( unsigned i=0; i< wdgt.size(); i++)
		wdgt[i]->updateFromList();

}

void IPTCPanel::updateChangeList ()
{

	for( unsigned i=0; i< wdgt.size(); i++)
		wdgt[i]->updateList();

}

void IPTCPanel::resetClicked ()
{
    if(idata){
    	disableListener ();
		chgList = idata->getIPTCData();
		applyChangeList();
		enableListener ();
		notifyListener ();
    }
}

void IPTCPanel::fileOpenClicked ()
{

    Gtk::FileChooserDialog dialog(M("IPTCPANEL_LOADDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN);
    if (options.multiUser)
       dialog.set_current_folder (Options::rtdir + "/iptc" );
    else
       dialog.set_current_folder (argv0 + "/iptc" );

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-open"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:

    Gtk::FileFilter filter_xmp;
    filter_xmp.set_name(M("IPTCPANEL_FILEDLGFILTERXMP"));
    filter_xmp.add_pattern("*"+paramFileExtension);
    dialog.add_filter(filter_xmp);

    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK)
    {
    	disableListener ();
        Glib::ustring filename( dialog.get_filename() );
        rtengine::ImageMetaData *id = rtengine::ImageMetaData::fromFile("",filename,"",false );
        if( id ){
        	rtengine::MetadataList loaded = id->getIPTCData();
        	for( rtengine::MetadataList::iterator iter = loaded.begin(); iter != loaded.end();iter++)
       			chgList[ iter->first ] = iter->second;
        	delete id;
        }
        applyChangeList ();
    }

    enableListener ();
    notifyListener ();
}

void IPTCPanel::fileSaveClicked ()
{
	updateChangeList ();

    Gtk::FileChooserDialog dialog(M("IPTCPANEL_SAVEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    if (options.multiUser)
       dialog.set_current_folder (Options::rtdir + "/iptc" );
    else
       dialog.set_current_folder (argv0 + "/iptc" );

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:
    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("IPTCPANEL_FILEDLGFILTERXMP"));
    filter_pp.add_pattern("*"+paramFileExtension);
    dialog.add_filter(filter_pp);

    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK)
    {
    	PartialPasteIPTCDlg selectDlg(chgList);
    	selectDlg.set_title (M("PARTIALPASTE_DIALOGIPTCSELECT"));
        if ( selectDlg.run ()) {
        	rtengine::MetadataList iptc = selectDlg.getIPTC();
        	Glib::ustring filename( dialog.get_filename() );
        	filename = removeExtension(filename)+".xmp";
        	rtengine::ImageMetaData *id = new rtengine::ImageMetaData("",filename,"");
        	if( id ){
				id->setIPTCData( iptc );
				id->saveXMP();
				delete id;
        	}
        }
    }
}

void IPTCPanel::copyClicked () {
	updateChangeList ();
    clipboard.setIPTC (chgList);
}

void IPTCPanel::pasteClicked () {

    disableListener ();
    chgList = clipboard.getIPTC ();
    applyChangeList ();
    enableListener ();
    notifyListener ();
}

/************************************ New widgets ************************************************/

void XRTWidget::attachToTable( Gtk::Table &table, int row, Gtk::Label &label, Gtk::Widget &widg )
{
	label.set_alignment(Gtk::ALIGN_RIGHT);
    table.attach (label, 0, 1, row, row+1, Gtk::FILL, Gtk::FILL, 1, 1);
    table.attach (widg, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 1, 1);
}

void XRTWidget::updateList()
{
	writeValue();
}

void XRTWidget::updateFromList()
{
	connChange.block(true);
	readValue();
	connChange.block(false);
}

XRTLabel::XRTLabel( Gtk::Table* table, int row, const std::string &key, rtengine::MetadataList &l ):XRTWidget( key,l )
{
	IPTCMeta  meta =IPTCMeta::IPTCtags[ key ];

	Gtk::Label* lab = Gtk::manage( new Gtk::Label (M(meta.guiName)+":") );
	control = Gtk::manage( new Gtk::Label () );
	lab->set_tooltip_text (M(meta.description));
	control->set_alignment(Gtk::ALIGN_LEFT);

	XRTWidget::attachToTable( *table, row, *lab, *control );
}

int XRTLabel::readValue( )
{
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
    	control->set_text ( iter->second[0] );
    }else{
    	control->set_text ("");
    }
	return 0;
}


XRTEntry::XRTEntry( Gtk::Table* table, int row, const std::string &key, rtengine::MetadataList &l ):XRTWidget( key,l )
{
	IPTCMeta  meta =IPTCMeta::IPTCtags[ key ];

	Gtk::Label* lab = Gtk::manage( new Gtk::Label (M(meta.guiName)+":") );
	lab->set_tooltip_text (M(meta.description));
	control = Gtk::manage( new Gtk::Entry () );
	connChange = control->signal_changed().connect( sigc::mem_fun(*this, &XRTWidget::updateList) );

	XRTWidget::attachToTable( *table, row, *lab, *control );
}

int XRTEntry::readValue( )
{
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
    	control->set_text ( iter->second[0] );
    }else{
    	control->set_text ("");
    }
	return 0;
}

int XRTEntry::writeValue( )
{
	Glib::ustring newValue(control->get_text ());
	Glib::ustring oldValue;
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
    	oldValue = iter->second[0];
    }
    if( !oldValue.empty() || !newValue.empty() ){
    	std::vector<Glib::ustring> v;
    	v.push_back(  newValue );
    	list[ key ] = v; // if newValue empty, but old was not: value="" means delete it!
	}else if( iter != list.end())
    	list.erase(iter);
	return 0;
}

XRTEntryMultiline::XRTEntryMultiline( Gtk::Table* table, int row, const std::string &key, rtengine::MetadataList &l ):XRTWidget( key,l )
{
	IPTCMeta  &meta =IPTCMeta::IPTCtags[ key ];

    Gtk::Label* lab = Gtk::manage( new Gtk::Label (M(meta.guiName)+":") );
    lab->set_tooltip_text (M(meta.description));

    control = Gtk::TextBuffer::create ();
    Gtk::TextView*  captionView = Gtk::manage( new Gtk::TextView (control) );
    Gtk::ScrolledWindow* scrolledWindowc = Gtk::manage( new Gtk::ScrolledWindow() );
    scrolledWindowc->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    scrolledWindowc->add(*captionView);
    connChange = control->signal_changed().connect( sigc::mem_fun(*this, &XRTWidget::updateList) );

    XRTWidget::attachToTable( *table, row, *lab, *scrolledWindowc );
}
int XRTEntryMultiline::readValue( )
{
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
    	control->set_text ( iter->second[0] );
    }else{
    	control->set_text ("");
    }
	return 0;
}

int XRTEntryMultiline::writeValue( )
{
	Glib::ustring newValue(control->get_text ());
	Glib::ustring oldValue;
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
    	oldValue = iter->second[0];
    }
    if( !oldValue.empty() || !newValue.empty() ){
    	std::vector<Glib::ustring> v;
    	v.push_back(  newValue );
    	list[ key ] = v; // if newValue empty, but old was not: value="" means delete it!
    }else if( iter != list.end())
    	list.erase(iter);
	return 0;
}

XRTEntryMultivalue::XRTEntryMultivalue( Gtk::Table* table, int row, const std::string &key, rtengine::MetadataList &l ):XRTWidget( key,l )
{
	IPTCMeta  &meta =IPTCMeta::IPTCtags[ key ];

    Gtk::Label* keyl = Gtk::manage( new Gtk::Label (M(meta.guiName)+":"));
    keyl->set_alignment(Gtk::ALIGN_RIGHT);
    keyl->set_tooltip_text (M(meta.description));
    controlList = Gtk::manage( new Gtk::ListViewText (1, false, Gtk::SELECTION_MULTIPLE) );
    controlList->set_headers_visible (false);
    Gtk::ScrolledWindow* scrolledWindowkw = Gtk::manage( new Gtk::ScrolledWindow() );
    scrolledWindowkw->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    scrolledWindowkw->add(*controlList);

    control  = Gtk::manage( new MyComboBoxEntryText() );
    control->set_size_request (32, -1);
    addKW = Gtk::manage( new Gtk::Button () );
    delKW = Gtk::manage( new Gtk::Button () );
    addKWImg = Gtk::manage (new RTImage ("list-add-small.png"));
    delKWImg = Gtk::manage (new RTImage ("list-remove-red-small.png"));
    addKW->add (*addKWImg);
    delKW->add (*delKWImg);
    addKW->signal_clicked().connect( sigc::mem_fun(*this, &XRTEntryMultivalue::add) );
    delKW->signal_clicked().connect( sigc::mem_fun(*this, &XRTEntryMultivalue::del) );
    Gtk::HBox* kwhb = Gtk::manage( new Gtk::HBox () );
    kwhb->pack_start (*control);
    kwhb->pack_start (*addKW, Gtk::PACK_SHRINK, 2);
    kwhb->pack_start (*delKW, Gtk::PACK_SHRINK, 2);

    XRTWidget::attachToTable( *table, row, *keyl, *kwhb );
    // another row here ...
    table->attach (*scrolledWindowkw, 1, 2, row+1, row+2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
}

void XRTEntryMultivalue::add()
{
	control->get_entry()->select_region (0, control->get_entry()->get_text().size());

    for (int i=0; i<controlList->size(); i++)
        if (controlList->get_text (i) == control->get_entry()->get_text())
            return;

    controlList->append_text (control->get_entry()->get_text());
    control->prepend_text (control->get_entry()->get_text());
    std::vector<Glib::ustring> items;
    for (Gtk::TreeModel::iterator i = control->get_model()->children().begin(); i!=control->get_model()->children().end(); i++) {
        Glib::ustring s;
        i->get_value (0, s);
        items.push_back (s);
    }
    control->clear_items ();
    for (int i=0; i<10 && i<items.size(); i++)
    	control->append_text (items[i]);
    controlList->scroll_to_row (controlList->get_model()->get_path(--controlList->get_model()->children().end()));

    updateList();
}

void XRTEntryMultivalue::del()
{
    std::vector<int> selection = controlList->get_selected ();
    if (selection.size()>0) {
        std::vector<Glib::ustring> keep;
        for (int i=0; i<controlList->size(); i++)
            if (std::find (selection.begin(), selection.end(), i) == selection.end())
                keep.push_back (controlList->get_text (i));
        controlList->clear_items ();
        for (int i=0; i<keep.size(); i++)
        	controlList->append_text (keep[i]);
    }
    updateList();
}

int XRTEntryMultivalue::readValue( )
{
	controlList->clear_items ();
	control->get_entry()->set_text ("");
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
        for (int j=0; j< iter->second.size(); j++)
        	controlList->append_text ( iter->second[j] );
    }
	return 0;
}

int XRTEntryMultivalue::writeValue( )
{
	int newValue(controlList->size());
	int oldValue(0);
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
    	oldValue = iter->second.size();
    }
    if( oldValue || newValue ){
    	std::vector<Glib::ustring> v;
    	for( int i=0; i< newValue; i++ )
    		v.push_back( controlList->get_text(i));
    	list[ key ] = v;
    }else if( iter != list.end())
    	list.erase(iter);
	return 0;
}

XRTCombo::XRTCombo( Gtk::Table* table, int row, const std::string &key, rtengine::IPTCPairList_t &info, rtengine::MetadataList &l ):XRTWidget( key,l ),predefValues(info)
{
	IPTCMeta  meta =IPTCMeta::IPTCtags[ key ];

    Gtk::Label* lab = Gtk::manage( new Gtk::Label (M(meta.guiName)+":") );
    lab->set_tooltip_text (M(meta.description));
    control = Gtk::manage( new MyComboBoxEntryText() );
    control->set_size_request (32, -1);

    Glib::ustring s;
    for( IPTCPairList_t::iterator iter = info.begin(); iter != info.end(); iter ++){
    	Glib::ustring::size_type iPos = iter->second.find_first_of(':',0);
    	if( iPos == Glib::ustring::npos )
    		control->append_text( iter->first + ": " + iter->second.substr(0,32-iter->first.size()));
    	else
    		control->append_text( iter->first + ": " + iter->second.substr(0,iPos));
    }
    //control->set_tooltip_text (s);
    connChange = control->signal_changed().connect( sigc::mem_fun(*this, &XRTWidget::updateList) );
    XRTWidget::attachToTable( *table, row, *lab, *control );

}
int XRTCombo::readValue( )
{
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
    	std::string val( iter->second[0] );
    	rtengine::IPTCPairList_t::iterator vIter = predefValues.find( val );
    	if( vIter != predefValues.end()) {
    		Glib::ustring::size_type iPos = vIter->second.find_first_of(':',0);
    		if( iPos == Glib::ustring::npos ){
    			val = vIter->first +": " + vIter->second.substr(0,32);
    		}else
    			val = vIter->first +": " + vIter->second.substr(0,iPos);
    		control->set_tooltip_text ( vIter->second );
    	}
    	control->get_entry()->set_text (  val );
    }else{
    	control->get_entry()->set_text ("");
    }
	return 0;
}

int XRTCombo::writeValue( )
{
	Glib::ustring newValue(control->get_entry()->get_text ());
	Glib::ustring::size_type iPos = newValue.find_first_of(':',0);
	if( iPos != Glib::ustring::npos )
		newValue = newValue.substr(0,iPos );
	rtengine::IPTCPairList_t::iterator vIter = predefValues.find( newValue );
	if( vIter != predefValues.end())
		control->set_tooltip_text ( vIter->second );
	Glib::ustring oldValue;
	rtengine::MetadataList::iterator iter = list.find( key );
    if( iter != list.end() && !iter->second.empty() ){
    	oldValue = iter->second[0];
    }
    if( !oldValue.empty() || !newValue.empty() ){
    	std::vector<Glib::ustring> v;
    	v.push_back(  newValue );
    	list[ key ] = v; // if newValue empty, but old was not: value="" means delete it!
    }else if( iter != list.end())
    	list.erase(iter);
	return 0;
}

/*IPTCSectionTitle::IPTCSectionTitle(Gtk::Table* table, int row, Glib::ustring &label)
{
    Gtk::VBox* vbox1 = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hsept1 = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Label* Title = Gtk::manage( new Gtk::Label (label) );
    Title->set_use_markup (true);
    Title->set_alignment(Gtk::ALIGN_LEFT);
    vbox1->pack_start(*hsept1);
    vbox1->pack_start(*Title);
    table->attach (*vbox1, 0, 2, row, row+1, Gtk::FILL, Gtk::FILL, 4, 4);
}*/
