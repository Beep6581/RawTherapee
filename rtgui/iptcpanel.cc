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
#include "iptcpanel.h"
#include "clipboard.h"
#include "rtimage.h"

extern Glib::ustring argv0;

using namespace rtengine;
using namespace rtengine::procparams;

IPTCPanel::IPTCPanel () {

    set_border_width (2);

    Gtk::Table* iptc = Gtk::manage( new Gtk::Table (27, 2) );
    
    int row = 0;

    Gtk::Label* capl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_CAPTION")+":") );
    captionText = Gtk::TextBuffer::create ();
    captionView = Gtk::manage( new Gtk::TextView (captionText) );
    Gtk::ScrolledWindow* scrolledWindowc = Gtk::manage( new Gtk::ScrolledWindow() );
    scrolledWindowc->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    scrolledWindowc->add(*captionView);
    capl->set_tooltip_text (M("IPTCPANEL_CAPTIONHINT"));
    captionView->set_tooltip_text (M("IPTCPANEL_CAPTIONHINT"));
    captionView->set_size_request(32, 80);
    iptc->attach (*capl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*scrolledWindowc, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

    Gtk::Label* capwl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_CAPTIONWRITER")+":") );
    captionWriter = Gtk::manage( new Gtk::Entry () );
    capwl->set_tooltip_text (M("IPTCPANEL_CAPTIONWRITERHINT"));
    captionWriter->set_tooltip_text (M("IPTCPANEL_CAPTIONWRITERHINT"));
    iptc->attach (*capwl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*captionWriter, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    
    row++;

    Gtk::Label* headl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_HEADLINE")+":") );
    headline = Gtk::manage( new Gtk::Entry () );
    headl->set_tooltip_text (M("IPTCPANEL_HEADLINEHINT"));
    headline->set_tooltip_text (M("IPTCPANEL_HEADLINEHINT"));
    iptc->attach (*headl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*headline, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

    Gtk::Label* instl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_INSTRUCTIONS")+":") );
    instructions = Gtk::manage( new Gtk::Entry () );
    instl->set_tooltip_text (M("IPTCPANEL_INSTRUCTIONSHINT"));
    instructions->set_tooltip_text (M("IPTCPANEL_INSTRUCTIONSHINT"));
    iptc->attach (*instl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*instructions, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    
    row++;

    Gtk::HSeparator* hsep1 = Gtk::manage( new Gtk::HSeparator () );
    iptc->attach (*hsep1, 0, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    
    row++;

    Gtk::Label* keyl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_KEYWORDS")+":"));
    keywords = Gtk::manage( new Gtk::ListViewText (1, false, Gtk::SELECTION_MULTIPLE) );
    keywords->set_headers_visible (false);
    keywords->set_size_request (50, 80);
    Gtk::ScrolledWindow* scrolledWindowkw = Gtk::manage( new Gtk::ScrolledWindow() );
    scrolledWindowkw->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    scrolledWindowkw->add(*keywords);
    keyword  = Gtk::manage( new Gtk::ComboBoxEntryText () );
    keyword->set_size_request (32, -1);
    keyl->set_tooltip_text (M("IPTCPANEL_KEYWORDSHINT"));
    keywords->set_tooltip_text (M("IPTCPANEL_KEYWORDSHINT"));
    keyword->set_tooltip_text (M("IPTCPANEL_KEYWORDSHINT"));
    addKW = Gtk::manage( new Gtk::Button () );
    delKW = Gtk::manage( new Gtk::Button () );
    Gtk::Image* addKWImg = Gtk::manage( new RTImage ("list-add-small.png") );
    Gtk::Image* delKWImg = Gtk::manage( new RTImage ("list-remove-red-small.png") );
    addKW->add (*addKWImg);
    delKW->add (*delKWImg);
    Gtk::HBox* kwhb = Gtk::manage( new Gtk::HBox () );
    kwhb->pack_start (*keyword);
    kwhb->pack_start (*addKW, Gtk::PACK_SHRINK, 2);
    kwhb->pack_start (*delKW, Gtk::PACK_SHRINK, 2);
    iptc->attach (*keyl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*kwhb, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    row++;
    iptc->attach (*scrolledWindowkw, 0, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    row++;

    Gtk::HSeparator* hsep2 = Gtk::manage( new Gtk::HSeparator () );
    iptc->attach (*hsep2, 0, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    row++;

    Gtk::Label* catl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_CATEGORY")+":") );
    category = Gtk::manage( new Gtk::ComboBoxEntryText () );
    category->set_size_request (32, -1);
    catl->set_tooltip_text (M("IPTCPANEL_CATEGORYHINT"));
    category->set_tooltip_text (M("IPTCPANEL_CATEGORYHINT"));
    Gtk::Label* scl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_SUPPCATEGORIES")+":") );
    suppCategories = Gtk::manage( new Gtk::ListViewText (1, false, Gtk::SELECTION_MULTIPLE) );
    suppCategories->set_headers_visible (false);
    suppCategories->set_size_request(50, 80);
    Gtk::ScrolledWindow* scrolledWindowsc = Gtk::manage( new Gtk::ScrolledWindow() );
    scrolledWindowsc->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    scrolledWindowsc->add(*suppCategories);
    suppCategory  = Gtk::manage( new Gtk::ComboBoxEntryText () );
    suppCategory->set_size_request (32, -1);
    scl->set_tooltip_text (M("IPTCPANEL_SUPPCATEGORIESHINT"));
    suppCategories->set_tooltip_text (M("IPTCPANEL_SUPPCATEGORIESHINT"));
    suppCategory->set_tooltip_text (M("IPTCPANEL_SUPPCATEGORIESHINT"));
    addSC = Gtk::manage( new Gtk::Button () );
    delSC = Gtk::manage( new Gtk::Button () );
    Gtk::Image* addSCImg = Gtk::manage( new RTImage ("list-add-small.png") );
    Gtk::Image* delSCImg = Gtk::manage( new RTImage ("list-remove-red-small.png") );
    addSC->add (*addSCImg);
    delSC->add (*delSCImg);
    Gtk::HBox* schb = Gtk::manage( new Gtk::HBox () );
    schb->pack_start (*suppCategory);
    schb->pack_start (*addSC, Gtk::PACK_SHRINK, 2);
    schb->pack_start (*delSC, Gtk::PACK_SHRINK, 2);
    iptc->attach (*catl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*category, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    row++;
    iptc->attach (*scl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*schb, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    row++;
    iptc->attach (*scrolledWindowsc, 0, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    row++;
    
    Gtk::HSeparator* hsep3 = Gtk::manage( new Gtk::HSeparator () );
    iptc->attach (*hsep3, 0, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);  
    row++;
    
    Gtk::Label* authl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_AUTHOR")+":") );
    author = Gtk::manage( new Gtk::Entry () );
    authl->set_tooltip_text (M("IPTCPANEL_CREDITHINT"));
    author->set_tooltip_text (M("IPTCPANEL_CREDITHINT"));
    iptc->attach (*authl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*author, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

    Gtk::Label* aupl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_AUTHORSPOSITION")+":") );
    authorPos = Gtk::manage(  new Gtk::Entry () );
    aupl->set_tooltip_text (M("IPTCPANEL_AUTHORSPOSITIONHINT"));
    authorPos->set_tooltip_text (M("IPTCPANEL_AUTHORSPOSITIONHINT"));
    iptc->attach (*aupl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*authorPos, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

    Gtk::Label* credl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_CREDIT")+":") );
    credit = Gtk::manage( new Gtk::Entry () );
    credl->set_tooltip_text (M("IPTCPANEL_CREDITHINT"));
    credit->set_tooltip_text (M("IPTCPANEL_CREDITHINT"));
    iptc->attach (*credl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*credit, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;
    
    Gtk::Label* sourl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_SOURCE")+":") );
    source = Gtk::manage( new Gtk::Entry () );
    sourl->set_tooltip_text (M("IPTCPANEL_SOURCEHINT"));
    source->set_tooltip_text (M("IPTCPANEL_SOURCEHINT"));
    iptc->attach (*sourl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*source, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

    Gtk::Label* cprl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_COPYRIGHT")+":") );
    copyright = Gtk::manage( new Gtk::Entry () );
    cprl->set_tooltip_text (M("IPTCPANEL_COPYRIGHTHINT"));
    copyright->set_tooltip_text (M("IPTCPANEL_COPYRIGHTHINT"));
    iptc->attach (*cprl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*copyright, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

    Gtk::HSeparator* hsep4 = Gtk::manage( new Gtk::HSeparator () );
    iptc->attach (*hsep4, 0, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    
    row++;

    Gtk::Label* cityl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_CITY")+":") );
    city = Gtk::manage( new Gtk::Entry () );
    cityl->set_tooltip_text (M("IPTCPANEL_CITYHINT"));
    city->set_tooltip_text (M("IPTCPANEL_CITYHINT"));
    iptc->attach (*cityl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*city, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;
    
    Gtk::Label* provl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_PROVINCE")+":") );
    province = Gtk::manage( new Gtk::Entry () );
    provl->set_tooltip_text (M("IPTCPANEL_PROVINCEHINT"));
    province->set_tooltip_text (M("IPTCPANEL_PROVINCEHINT"));
    iptc->attach (*provl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*province, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;
    
    Gtk::Label* ctrl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_COUNTRY")+":") );
    country = Gtk::manage( new Gtk::Entry () );
    ctrl->set_tooltip_text (M("IPTCPANEL_COUNTRYHINT"));
    country->set_tooltip_text (M("IPTCPANEL_COUNTRYHINT"));
    iptc->attach (*ctrl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*country, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;
    
    Gtk::Label* titll = Gtk::manage( new Gtk::Label (M("IPTCPANEL_TITLE")+":") );
    title = Gtk::manage( new Gtk::Entry () );
    titll->set_tooltip_text (M("IPTCPANEL_TITLEHINT"));
    title->set_tooltip_text (M("IPTCPANEL_TITLEHINT"));
    iptc->attach (*titll, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*title, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

    Gtk::Label* dcl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_DATECREATED")+":") );
    dateCreated = Gtk::manage(  new Gtk::Entry () );
    dcl->set_tooltip_text (M("IPTCPANEL_DATECREATEDHINT"));
    dateCreated->set_tooltip_text (M("IPTCPANEL_DATECREATEDHINT"));
    iptc->attach (*dcl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*dateCreated, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

    Gtk::Label* trl = Gtk::manage( new Gtk::Label (M("IPTCPANEL_TRANSREFERENCE")+":") );
    transReference = Gtk::manage( new Gtk::Entry () );
    trl->set_tooltip_text (M("IPTCPANEL_TRANSREFERENCEHINT"));
    transReference->set_tooltip_text (M("IPTCPANEL_TRANSREFERENCEHINT"));
    iptc->attach (*trl, 0, 1, row, row+1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    iptc->attach (*transReference, 1, 2, row, row+1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    row++;

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
    bbox->pack_start (*reset);

    file = Gtk::manage( new Gtk::Button (M("IPTCPANEL_EMBEDDED")) );
    file->set_image (*Gtk::manage(new RTImage ("gtk-open.png")));
    bbox->pack_start (*file);
    
    copy = Gtk::manage( new Gtk::Button () );
    copy->set_image (*Gtk::manage(new RTImage ("edit-copy.png")));
    bbox->pack_start (*copy, Gtk::PACK_SHRINK, 0);

    paste = Gtk::manage( new Gtk::Button () );
    paste->set_image (*Gtk::manage(new RTImage ("edit-paste.png")));
    bbox->pack_start (*paste, Gtk::PACK_SHRINK, 0);
    
    pack_end (*bbox, Gtk::PACK_SHRINK, 2);

    Gtk::Tooltips* toolTip = Gtk::manage( new Gtk::Tooltips () );
    toolTip->set_tip (*reset, M("IPTCPANEL_RESETHINT"));
    toolTip->set_tip (*file, M("IPTCPANEL_EMBEDDEDHINT"));
    toolTip->set_tip (*copy, M("IPTCPANEL_COPYHINT"));
    toolTip->set_tip (*paste, M("IPTCPANEL_PASTEHINT"));
    
    reset->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::resetClicked) );
    file->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::fileClicked) );
    copy->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::copyClicked) );
    paste->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::pasteClicked) );


   addKW->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::addKeyWord) );
   delKW->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::delKeyWord) );
   addSC->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::addSuppCategory) );
   delSC->signal_clicked().connect( sigc::mem_fun(*this, &IPTCPanel::delSuppCategory) );
   keyword->get_entry()->signal_activate().connect( sigc::mem_fun(*this, &IPTCPanel::addKeyWord) );
   suppCategory->get_entry()->signal_activate().connect( sigc::mem_fun(*this, &IPTCPanel::addSuppCategory) );

   conns[0] = captionText->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[1] = captionWriter->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[2] = headline->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[3] = instructions->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[4] = category->get_entry()->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[5] = author->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[6] = authorPos->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[7] = credit->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[8] = source->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[9] = copyright->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[10] = city->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[11] = province->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[12] = country->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[13] = title->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[14] = dateCreated->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   conns[15] = transReference->signal_changed().connect( sigc::mem_fun(*this, &IPTCPanel::updateChangeList) );
   
   category->get_entry()->set_max_length (3);
   keyword->get_entry()->set_max_length (64);
   captionWriter->set_max_length (32);
   instructions->set_max_length (256);
   author->set_max_length (32);
   authorPos->set_max_length (32);
   credit->set_max_length (32);
   source->set_max_length (32);
   copyright->set_max_length (128);
   city->set_max_length (32);
   province->set_max_length (32);
   country->set_max_length (64);
   title->set_max_length (64);
   dateCreated->set_max_length (8);
   transReference->set_max_length (32);
   
    show_all ();
}

void IPTCPanel::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();
    changeList.clear();
    if (!pp->iptc.empty())
        changeList = pp->iptc;
    else
        changeList = embeddedData;
    applyChangeList ();
    enableListener ();
}

void IPTCPanel::write (ProcParams* pp, ParamsEdited* pedited) {
    
    pp->iptc = changeList;
}

void IPTCPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    defChangeList = defParams->iptc;
}
        
void IPTCPanel::setImageData (const ImageMetaData* id) {
    
    if (id) 
        embeddedData = id->getIPTCData ();
    else
        embeddedData.clear ();

    file->set_sensitive (!embeddedData.empty());
}

void IPTCPanel::notifyListener () {
    
    if (listener)
        listener->panelChanged (EvIPTC, M("HISTORY_CHANGED"));
}

void IPTCPanel::addKeyWord () {

    keyword->get_entry()->select_region (0, keyword->get_entry()->get_text().size());

    for (unsigned int i=0; i<keywords->size(); i++)
        if (keywords->get_text (i) == keyword->get_entry()->get_text()) 
            return;
        
    keywords->append_text (keyword->get_entry()->get_text());
    keyword->prepend_text (keyword->get_entry()->get_text());
    std::vector<Glib::ustring> items;
    for (Gtk::TreeModel::iterator i = keyword->get_model()->children().begin(); i!=keyword->get_model()->children().end(); i++) {
        Glib::ustring s;
        i->get_value (0, s);
        items.push_back (s);
    }
    keyword->clear_items ();
    for (unsigned int i=0; i<10 && i<items.size(); i++)
        keyword->append_text (items[i]);
    keywords->scroll_to_row (keywords->get_model()->get_path(--keywords->get_model()->children().end()));
    
    updateChangeList ();
}

void IPTCPanel::delKeyWord () {

    std::vector<int> selection = keywords->get_selected ();
    if (!selection.empty()) {
        std::vector<Glib::ustring> keep;
        for (unsigned int i=0; i<keywords->size(); i++)
            if (std::find (selection.begin(), selection.end(), i) == selection.end())
                keep.push_back (keywords->get_text (i));
        keywords->clear_items ();
        for (unsigned int i=0; i<keep.size(); i++)
            keywords->append_text (keep[i]);
    }

    updateChangeList ();
}

void IPTCPanel::addSuppCategory () {

    for (unsigned int i=0; i<suppCategories->size(); i++)
        if (suppCategories->get_text (i) == suppCategory->get_entry()->get_text())
            return;

    suppCategories->append_text (suppCategory->get_entry()->get_text());
    suppCategory->prepend_text (suppCategory->get_entry()->get_text());
    std::vector<Glib::ustring> items;
    for (Gtk::TreeModel::iterator i = suppCategory->get_model()->children().begin(); i!=suppCategory->get_model()->children().end(); i++) {
        Glib::ustring s;
        i->get_value (0, s);
        items.push_back (s);
    }
    suppCategory->clear_items ();
    for (unsigned int i=0; i<10 && i<items.size(); i++)
        suppCategory->append_text (items[i]);
    suppCategories->scroll_to_row (suppCategories->get_model()->get_path(--suppCategories->get_model()->children().end()));
    suppCategory->get_entry()->select_region (0, suppCategory->get_entry()->get_text().size());

    updateChangeList ();    
}

void IPTCPanel::delSuppCategory () {

    std::vector<int> selection = suppCategories->get_selected ();
    if (!selection.empty()) {
        std::vector<Glib::ustring> keep;
        for (unsigned int i=0; i<suppCategories->size(); i++)
            if (std::find (selection.begin(), selection.end(), i) == selection.end())
                keep.push_back (suppCategories->get_text (i));
        suppCategories->clear_items ();
        for (unsigned int i=0; i<keep.size(); i++)
            suppCategories->append_text (keep[i]);
    }

    updateChangeList ();
}

void IPTCPanel::updateChangeList () {

    changeList.clear ();
    changeList["Caption"        ].push_back (captionText->get_text ());
    changeList["CaptionWriter"  ].push_back (captionWriter->get_text ());
    changeList["Headline"       ].push_back (headline->get_text ());
    changeList["Instructions"   ].push_back (instructions->get_text ());

    for (unsigned int i=0; i<keywords->size(); i++)
    changeList["Keywords"       ].push_back (keywords->get_text (i));

    changeList["Category"       ].push_back (category->get_entry()->get_text ());

    for (unsigned int i=0; i<suppCategories->size(); i++)
    changeList["SupplementalCategories"].push_back (suppCategories->get_text (i));

    changeList["Author"         ].push_back (author->get_text ());
    changeList["AuthorsPosition"].push_back (authorPos->get_text ());
    changeList["Credit"         ].push_back (credit->get_text ());
    changeList["Source"         ].push_back (source->get_text ());
    changeList["Copyright"      ].push_back (copyright->get_text ());
    changeList["City"           ].push_back (city->get_text ());
    changeList["Province"       ].push_back (province->get_text ());
    changeList["Country"        ].push_back (country->get_text ());
    changeList["Title"          ].push_back (title->get_text ());
    changeList["DateCreated"    ].push_back (dateCreated->get_text ());
    changeList["TransReference" ].push_back (transReference->get_text ());

    notifyListener ();
}

void IPTCPanel::applyChangeList () {

    for (int i=0; i<16; i++)
        conns[i].block (true);

    captionText->set_text ("");
    captionWriter->set_text ("");
    headline->set_text ("");
    instructions->set_text ("");
    keywords->clear_items ();
    category->get_entry()->set_text ("");
    suppCategories->clear_items ();
    author->set_text ("");
    authorPos->set_text ("");
    credit->set_text ("");
    source->set_text ("");
    copyright->set_text ("");
    city->set_text ("");
    province->set_text ("");
    country->set_text ("");
    title->set_text ("");
    dateCreated->set_text ("");
    transReference->set_text ("");
    keyword->get_entry()->set_text ("");
    suppCategory->get_entry()->set_text ("");
    
    for (rtengine::procparams::IPTCPairs::iterator i=changeList.begin(); i!=changeList.end(); i++) {
        if (i->first == "Caption" && !i->second.empty())
            captionText->set_text (i->second.at(0));
        else if (i->first == "CaptionWriter" && !i->second.empty())
            captionWriter->set_text (i->second.at(0));
        else if (i->first == "Headline" && !i->second.empty())
            headline->set_text (i->second.at(0));
        else if (i->first == "Instructions" && !i->second.empty())
            instructions->set_text (i->second.at(0));
        else if (i->first == "Keywords")
            for (unsigned int j=0; j<i->second.size(); j++)
                keywords->append_text (i->second.at(j));
        else if (i->first == "Category" && !i->second.empty())
            category->get_entry()->set_text (i->second.at(0));
        else if (i->first == "SupplementalCategories")
            for (unsigned int j=0; j<i->second.size(); j++)
                suppCategories->append_text (i->second.at(j));
        else if (i->first == "Author" && !i->second.empty())
            author->set_text (i->second.at(0));
        else if (i->first == "AuthorsPosition" && !i->second.empty())
            authorPos->set_text (i->second.at(0));
        else if (i->first == "Credit" && !i->second.empty())
            credit->set_text (i->second.at(0));
        else if (i->first == "Source" && !i->second.empty())
            source->set_text (i->second.at(0));
        else if (i->first == "Copyright" && !i->second.empty())
            copyright->set_text (i->second.at(0));
        else if (i->first == "City" && !i->second.empty())
            city->set_text (i->second.at(0));
        else if (i->first == "Province" && !i->second.empty())
            province->set_text (i->second.at(0));
        else if (i->first == "Country" && !i->second.empty())
            country->set_text (i->second.at(0));
        else if (i->first == "Title" && !i->second.empty())
            title->set_text (i->second.at(0));
        else if (i->first == "DateCreated" && !i->second.empty())
            dateCreated->set_text (i->second.at(0));
        else if (i->first == "TransReference" && !i->second.empty())
            transReference->set_text (i->second.at(0));
}

    for (int i=0; i<16; i++)
        conns[i].block (false);
}

void IPTCPanel::resetClicked () {

    disableListener ();
    changeList = defChangeList;
    applyChangeList ();
    enableListener ();
    notifyListener ();
}

void IPTCPanel::fileClicked () {

    disableListener ();
    changeList = embeddedData;
    applyChangeList ();
    enableListener ();
    notifyListener ();
}

void IPTCPanel::copyClicked () {

    clipboard.setIPTC (changeList);
}

void IPTCPanel::pasteClicked () {

    disableListener ();
    changeList = clipboard.getIPTC ();
    applyChangeList ();
    enableListener ();
    notifyListener ();
}
