/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2012 Michael Ezra <michael@michaelezra.com>
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
#include "exportpanel.h"
#include "multilangmgr.h"
#include "options.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

ExportPanel::ExportPanel () : listener (NULL) {

	set_border_width (4);

	/*enabled = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_ENABLE")) );
	pack_start(*enabled, Gtk::PACK_SHRINK, 4);
	pack_start (*Gtk::manage(new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 2);*/

	Gtk::Label* labExportTitle = Gtk::manage ( new Gtk::Label (M("EXPORT_FASTEXPORTOPTIONS")) );
	labExportTitle->set_use_markup (true);
	labExportTitle->set_alignment(Gtk::ALIGN_LEFT);
	pack_start(*labExportTitle, Gtk::PACK_SHRINK, 4);

	Gtk::Label* labInstructions = Gtk::manage ( new Gtk::Label (M("EXPORT_INSTRUCTIONS")) );
	labInstructions->set_use_markup (true);
	labInstructions->set_line_wrap (true);
	labInstructions->set_alignment(Gtk::ALIGN_LEFT);
	Gtk::HBox* hbInstructions = Gtk::manage (new Gtk::HBox ());
	hbInstructions->pack_start(*labInstructions, Gtk::PACK_SHRINK, 4);
	pack_start(*hbInstructions, Gtk::PACK_SHRINK, 4);

	bypass_ALL              = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_ALL")));
	bypass_sharpening       = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_SHARPENING")));
	bypass_sharpenEdge      = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_SHARPENEDGE")));
	bypass_sharpenMicro     = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_SHARPENMICRO")));
	bypass_lumaDenoise      = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_LUMADENOISE")));
	bypass_colorDenoise     = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_COLORDENOISE")));
	bypass_defringe         = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_DEFRINGE")));
	bypass_dirpyrDenoise    = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_DIRPYRDENOISE")));
	bypass_sh_hq            = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_SH_HQ")));
	bypass_dirpyrequalizer  = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_DIRPYREQUALIZER")));
	bypass_raw_all_enhance  = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_ALL_ENHANCE")));
	bypass_raw_ccSteps      = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_CCSTEPS")));
	bypass_raw_linenoise    = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_LINENOISE")));
	bypass_raw_greenthresh  = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_GREENTHRESH")));
	bypass_raw_ca           = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_CA")));
	bypass_raw_df           = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_DF")));
	bypass_raw_ff           = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_FF")));
	bypass_raw_dcb_iterations = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_DCB_ITERATIONS")));
	bypass_raw_dcb_enhance    = Gtk::manage ( new Gtk::CheckButton (M("EXPORT_BYPASS_RAW_DCB_ENHANCE")));

	Gtk::HBox* hb_raw_dmethod = Gtk::manage (new Gtk::HBox ());
	hb_raw_dmethod->pack_start (*Gtk::manage (new Gtk::Label ( M("EXPORT_RAW_DMETHOD") +": ")),Gtk::PACK_SHRINK, 4);
	raw_dmethod = Gtk::manage (new MyComboBoxText ());
	for( size_t i=0; i< procparams::RAWParams::numMethods;i++)
		raw_dmethod->append_text(procparams::RAWParams::methodstring[i]);

	raw_dmethod->set_active(0);
	hb_raw_dmethod->pack_end (*raw_dmethod, Gtk::PACK_EXPAND_WIDGET, 4);

	// start packing
	pack_start(*bypass_ALL          , Gtk::PACK_SHRINK, 4);
	pack_start(*Gtk::manage(new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_sharpening   , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_sharpenEdge  , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_sharpenMicro , Gtk::PACK_SHRINK, 4);
	//pack_start(*bypass_lumaDenoise  , Gtk::PACK_SHRINK, 4);
	//pack_start(*bypass_colorDenoise , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_defringe     , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_dirpyrDenoise, Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_sh_hq        , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_dirpyrequalizer , Gtk::PACK_SHRINK, 4);

	pack_start(*hb_raw_dmethod         , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_raw_all_enhance , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_raw_ccSteps     , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_raw_dcb_iterations, Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_raw_dcb_enhance   , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_raw_ca         , Gtk::PACK_SHRINK, 4);

	pack_start(*bypass_raw_linenoise   , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_raw_greenthresh , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_raw_df         , Gtk::PACK_SHRINK, 4);
	pack_start(*bypass_raw_ff         , Gtk::PACK_SHRINK, 4);

	pack_start (*Gtk::manage(new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 2);

	// Resize options

	resize_method = Gtk::manage (new MyComboBoxText ());
	resize_method->append_text (M("TP_RESIZE_NEAREST"));
	resize_method->append_text (M("TP_RESIZE_BILINEAR"));
	resize_method->append_text (M("TP_RESIZE_BICUBIC"));
	resize_method->append_text (M("TP_RESIZE_BICUBICSF"));
	resize_method->append_text (M("TP_RESIZE_BICUBICSH"));
	resize_method->append_text (M("TP_RESIZE_LANCZOS"));
	resize_method->set_active (5);

	Gtk::HBox* rmbox = Gtk::manage (new Gtk::HBox ());
	rmbox->pack_start (*Gtk::manage (new Gtk::Label (M("EXPORT_RESIZEMETHOD"))), Gtk::PACK_SHRINK, 4);
	rmbox->pack_start (*resize_method);
	pack_start (*rmbox, Gtk::PACK_SHRINK, 4);

	Gtk::HBox* wbox = Gtk::manage (new Gtk::HBox ());
	Gtk::HBox* hbox = Gtk::manage (new Gtk::HBox ());
	MaxWidth = Gtk::manage (new MySpinButton ());
	MaxHeight = Gtk::manage (new MySpinButton ());
	wbox->pack_start (*Gtk::manage (new Gtk::Label (M("EXPORT_MAXWIDTH"))), Gtk::PACK_SHRINK, 4);
	wbox->pack_start (*MaxWidth);
	hbox->pack_start (*Gtk::manage (new Gtk::Label (M("EXPORT_MAXHEIGHT"))), Gtk::PACK_SHRINK, 4);
	hbox->pack_start (*MaxHeight);
	pack_start (*wbox, Gtk::PACK_SHRINK, 4);
	pack_start (*hbox, Gtk::PACK_SHRINK, 4);

	MaxWidth->set_digits (0);
	MaxWidth->set_increments (1,100);
	MaxWidth->set_value (options.fastexport_resize_width);
	MaxWidth->set_range (32, 3000);

	MaxHeight->set_digits (0);
	MaxHeight->set_increments (1,100);
	MaxHeight->set_value (options.fastexport_resize_height);
	MaxHeight->set_range (32, 3000);

	// Buttons
	btnFastExport =  Gtk::manage ( new Gtk::Button (M("EXPORT_PUTTOQUEUEFAST")) );
	btnFastExport->set_image (*Gtk::manage (new RTImage ("processing.png")));
	pack_start(*btnFastExport, Gtk::PACK_SHRINK, 4);


    // add panel ending
    Gtk::VBox* vboxpe = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hseptpe = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Image* peImg = Gtk::manage (new RTImage("PanelEnding.png"));
    vboxpe->pack_start(*hseptpe, Gtk::PACK_SHRINK, 4);
    vboxpe->pack_start(*peImg);
    pack_start(*vboxpe, Gtk::PACK_SHRINK, 0);


	btnFastExport->signal_clicked().connect( sigc::mem_fun(*this, &ExportPanel::FastExportPressed) );
	//btnExportLoadSettings->signal_clicked().connect( sigc::mem_fun(*this, &ExportPanel::LoadSettings) );
	//btnExportSaveSettings->signal_clicked().connect( sigc::mem_fun(*this, &ExportPanel::SaveSettings) );
	bypass_ALLconn = bypass_ALL->signal_toggled().connect (sigc::mem_fun(*this, &ExportPanel::bypassALL_Toggled));

	bypass_sharpeningConn         = bypass_sharpening->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_sharpenEdgeConn        = bypass_sharpenEdge->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_sharpenMicroConn       = bypass_sharpenMicro->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_lumaDenoiseConn        = bypass_lumaDenoise->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_colorDenoiseConn       = bypass_colorDenoise->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_defringeConn           = bypass_defringe->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_dirpyrDenoiseConn      = bypass_dirpyrDenoise->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_sh_hqConn              = bypass_sh_hq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_dirpyrequalizerConn    = bypass_dirpyrequalizer->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_all_enhanceConn    = bypass_raw_all_enhance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_ccStepsConn        = bypass_raw_ccSteps->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_dcb_iterationsConn = bypass_raw_dcb_iterations->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_dcb_enhanceConn    = bypass_raw_dcb_enhance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_caConn             = bypass_raw_ca->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_linenoiseConn      = bypass_raw_linenoise->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_greenthreshConn    = bypass_raw_greenthresh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_dfConn             = bypass_raw_df->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));
	bypass_raw_ffConn             = bypass_raw_ff->signal_toggled().connect (sigc::bind (sigc::mem_fun(*bypass_ALL, &Gtk::CheckButton::set_inconsistent), true));

	LoadDefaultSettings();
}

/*bool ExportPanel::isEnabled () {

	return enabled->get_active () && is_sensitive();
}*/


void ExportPanel::FastExportPressed (){
	// options is the container for making these settings available globally.
	// Therefore, settings must be saved to options before they are used further
	SaveSettingsAsDefault();

	if (listener)
		listener->exportRequested ();
}

void ExportPanel::SaveSettingsAsDefault(){
	// Save fast export settings to options
	options.fastexport_bypass_sharpening         = bypass_sharpening->get_active        ();
	options.fastexport_bypass_sharpenEdge        = bypass_sharpenEdge->get_active       ();
	options.fastexport_bypass_sharpenMicro       = bypass_sharpenMicro->get_active      ();
	options.fastexport_bypass_lumaDenoise        = bypass_lumaDenoise->get_active       ();
	options.fastexport_bypass_colorDenoise       = bypass_colorDenoise->get_active      ();
	options.fastexport_bypass_defringe           = bypass_defringe->get_active          ();
	options.fastexport_bypass_dirpyrDenoise      = bypass_dirpyrDenoise->get_active     ();
	options.fastexport_bypass_sh_hq              = bypass_sh_hq->get_active             ();
	options.fastexport_bypass_dirpyrequalizer    = bypass_dirpyrequalizer->get_active   ();
	options.fastexport_bypass_raw_all_enhance    = bypass_raw_all_enhance->get_active   ();
	options.fastexport_bypass_raw_ccSteps        = bypass_raw_ccSteps->get_active       ();
	options.fastexport_bypass_raw_dcb_iterations = bypass_raw_dcb_iterations->get_active();
	options.fastexport_bypass_raw_dcb_enhance    = bypass_raw_dcb_enhance->get_active   ();
	options.fastexport_bypass_raw_ca             = bypass_raw_ca->get_active            ();
	options.fastexport_bypass_raw_linenoise      = bypass_raw_linenoise->get_active     ();
	options.fastexport_bypass_raw_greenthresh    = bypass_raw_greenthresh->get_active   ();
	options.fastexport_bypass_raw_df             = bypass_raw_df->get_active            ();
	options.fastexport_bypass_raw_ff             = bypass_raw_ff->get_active            ();

	//saving demosaic_method
	int currentRow = raw_dmethod->get_active_row_number();
	if( currentRow>=0 && currentRow < procparams::RAWParams::numMethods)
		options.fastexport_raw_dmethod = procparams::RAWParams::methodstring[currentRow];

//	options.fastexport_icm_input        = icm_input       ;
//	options.fastexport_icm_working      = icm_working     ;
//	options.fastexport_icm_output       = icm_output      ;
//	options.fastexport_icm_gamma        = icm_gamma       ;
//	options.fastexport_resize_enabled   = resize_enabled  ;
//	options.fastexport_resize_scale     = resize_scale    ;
//	options.fastexport_resize_appliesTo = resize_appliesTo;
//	options.fastexport_resize_dataspec  = resize_dataspec ;

	options.fastexport_resize_method = "Lanczos";
    if (resize_method->get_active_row_number() == 0)
        options.fastexport_resize_method = "Nearest";
    else if (resize_method->get_active_row_number() == 1)
        options.fastexport_resize_method = "Bilinear";
    else if (resize_method->get_active_row_number() == 2)
        options.fastexport_resize_method = "Bicubic";
    else if (resize_method->get_active_row_number() == 3)
        options.fastexport_resize_method = "Bicubic (Softer)";
    else if (resize_method->get_active_row_number() == 4)
        options.fastexport_resize_method = "Bicubic (Sharper)";
    else if (resize_method->get_active_row_number() == 5)
        options.fastexport_resize_method = "Lanczos";

	options.fastexport_resize_width     = MaxWidth->get_value_as_int ();
	options.fastexport_resize_height    = MaxHeight->get_value_as_int ();
}

void ExportPanel::LoadDefaultSettings(){
	// Load fast export settings from options
	bypass_sharpening->set_active        (options.fastexport_bypass_sharpening         );
	bypass_sharpenEdge->set_active       (options.fastexport_bypass_sharpenEdge        );
	bypass_sharpenMicro->set_active      (options.fastexport_bypass_sharpenMicro       );
	bypass_lumaDenoise->set_active       (options.fastexport_bypass_lumaDenoise        );
	bypass_colorDenoise->set_active      (options.fastexport_bypass_colorDenoise       );
	bypass_defringe->set_active          (options.fastexport_bypass_defringe           );
	bypass_dirpyrDenoise->set_active     (options.fastexport_bypass_dirpyrDenoise      );
	bypass_sh_hq->set_active             (options.fastexport_bypass_sh_hq              );
	bypass_dirpyrequalizer->set_active   (options.fastexport_bypass_dirpyrequalizer    );
	bypass_raw_all_enhance->set_active   (options.fastexport_bypass_raw_all_enhance    );
	bypass_raw_ccSteps->set_active       (options.fastexport_bypass_raw_ccSteps        );
	bypass_raw_dcb_iterations->set_active(options.fastexport_bypass_raw_dcb_iterations );
	bypass_raw_dcb_enhance->set_active   (options.fastexport_bypass_raw_dcb_enhance    );
	bypass_raw_ca->set_active            (options.fastexport_bypass_raw_ca             );
	bypass_raw_linenoise->set_active     (options.fastexport_bypass_raw_linenoise      );
	bypass_raw_greenthresh->set_active   (options.fastexport_bypass_raw_greenthresh    );
	bypass_raw_df->set_active            (options.fastexport_bypass_raw_df             );
	bypass_raw_ff->set_active            (options.fastexport_bypass_raw_ff             );

	//demosaic method
	raw_dmethod->set_active(procparams::RAWParams::numMethods);
    for( size_t i=0; i< procparams::RAWParams::numMethods;i++)
	   if( options.fastexport_raw_dmethod == procparams::RAWParams::methodstring[i]){
		   raw_dmethod->set_active(i);
		   break;
	   }

//    icm_input        = options.fastexport_icm_input       ;
//    icm_working      = options.fastexport_icm_working     ;
//    icm_output       = options.fastexport_icm_output      ;
//    icm_gamma        = options.fastexport_icm_gamma       ;
//    resize_enabled   = options.fastexport_resize_enabled  ;
//    resize_scale     = options.fastexport_resize_scale    ;
//    resize_appliesTo = options.fastexport_resize_appliesTo;
//    resize_dataspec  = options.fastexport_resize_dataspec ;

    resize_method->set_active (2);
	if (options.fastexport_resize_method == "Nearest")
	  resize_method->set_active (0);
	else if (options.fastexport_resize_method == "Bilinear")
	  resize_method->set_active (1);
	else if (options.fastexport_resize_method == "Bicubic")
	  resize_method->set_active (2);
	else if (options.fastexport_resize_method == "Bicubic (Softer)")
	  resize_method->set_active (3);
	else if (options.fastexport_resize_method == "Bicubic (Sharper)")
	  resize_method->set_active (4);
	else if (options.fastexport_resize_method == "Lanczos")
	  resize_method->set_active (5);
	else if (options.fastexport_resize_method == "Downscale (Better)" ||
		   options.fastexport_resize_method == "Downscale (Faster)")
	{
	  resize_method->set_active (5);
	}

    MaxWidth->set_value(options.fastexport_resize_width);
    MaxHeight->set_value(options.fastexport_resize_height);
}

void ExportPanel::LoadSettings(){
	// load settings from file
}

void ExportPanel::SaveSettings(){
	// save settings to file

}

void ExportPanel::bypassALL_Toggled(){
	bypass_sharpeningConn.block         (true);
	bypass_sharpenEdgeConn.block        (true);
	bypass_sharpenMicroConn.block       (true);
	bypass_lumaDenoiseConn.block        (true);
	bypass_colorDenoiseConn.block       (true);
	bypass_defringeConn.block           (true);
	bypass_dirpyrDenoiseConn.block      (true);
	bypass_sh_hqConn.block              (true);
	bypass_dirpyrequalizerConn.block    (true);
	bypass_raw_all_enhanceConn.block    (true);
	bypass_raw_ccStepsConn.block        (true);
	bypass_raw_dcb_iterationsConn.block (true);
	bypass_raw_dcb_enhanceConn.block    (true);
	bypass_raw_caConn.block             (true);
	bypass_raw_linenoiseConn.block      (true);
	bypass_raw_greenthreshConn.block    (true);
	bypass_raw_dfConn.block             (true);
	bypass_raw_ffConn.block             (true);

	bypass_ALL->set_inconsistent (false);

	bypass_sharpening->set_active(bypass_ALL->get_active());
	bypass_sharpenEdge->set_active(bypass_ALL->get_active());
	bypass_sharpenMicro->set_active(bypass_ALL->get_active());
	bypass_lumaDenoise->set_active(bypass_ALL->get_active());
	bypass_colorDenoise->set_active(bypass_ALL->get_active());
	bypass_defringe->set_active(bypass_ALL->get_active());
	bypass_dirpyrDenoise->set_active(bypass_ALL->get_active());
	bypass_sh_hq->set_active(bypass_ALL->get_active());
	bypass_dirpyrequalizer->set_active(bypass_ALL->get_active());
	bypass_raw_all_enhance->set_active(bypass_ALL->get_active());
	bypass_raw_ccSteps->set_active(bypass_ALL->get_active());
	bypass_raw_dcb_iterations->set_active(bypass_ALL->get_active());
	bypass_raw_dcb_enhance->set_active(bypass_ALL->get_active());
	bypass_raw_ca->set_active(bypass_ALL->get_active());
	bypass_raw_linenoise->set_active(bypass_ALL->get_active());
	bypass_raw_greenthresh->set_active(bypass_ALL->get_active());
	bypass_raw_df->set_active(bypass_ALL->get_active());
	bypass_raw_ff->set_active(bypass_ALL->get_active());

	bypass_sharpeningConn.block         (false);
	bypass_sharpenEdgeConn.block        (false);
	bypass_sharpenMicroConn.block       (false);
	bypass_lumaDenoiseConn.block        (false);
	bypass_colorDenoiseConn.block       (false);
	bypass_defringeConn.block           (false);
	bypass_dirpyrDenoiseConn.block      (false);
	bypass_sh_hqConn.block              (false);
	bypass_dirpyrequalizerConn.block    (false);
	bypass_raw_all_enhanceConn.block    (false);
	bypass_raw_ccStepsConn.block        (false);
	bypass_raw_dcb_iterationsConn.block (false);
	bypass_raw_dcb_enhanceConn.block    (false);
	bypass_raw_caConn.block             (false);
	bypass_raw_linenoiseConn.block      (false);
	bypass_raw_greenthreshConn.block    (false);
	bypass_raw_dfConn.block             (false);
	bypass_raw_ffConn.block             (false);
}

/*
fastexport_bypass_sharpening
fastexport_bypass_sharpenEdge
fastexport_bypass_sharpenMicro
fastexport_bypass_lumaDenoise
fastexport_bypass_colorDenoise
fastexport_bypass_defringe
fastexport_bypass_dirpyrDenoise
fastexport_bypass_sh_hq
fastexport_bypass_dirpyrequalizer
fastexport_bypass_raw_all_enhance
fastexport_bypass_raw_ccSteps
fastexport_bypass_raw_dcb_iterations
fastexport_bypass_raw_dcb_enhance
fastexport_bypass_raw_ca
fastexport_bypass_raw_linenoise
fastexport_bypass_raw_greenthresh
fastexport_bypass_raw_df
fastexport_bypass_raw_ff
fastexport_raw_dmethod
fastexport_icm_input
fastexport_icm_working
fastexport_icm_output
fastexport_icm_gamma
fastexport_resize_enabled
fastexport_resize_scale
fastexport_resize_appliesTo
fastexport_resize_method
fastexport_resize_dataspec
fastexport_resize_width
fastexport_resize_height
*/
