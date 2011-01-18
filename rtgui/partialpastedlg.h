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
#ifndef _PARTIALPASTEDLG_
#define _PARTIALPASTEDLG_

#include <gtkmm.h>
#include <rtengine.h>

class PartialPasteDlg : public Gtk::Dialog {

    public:
    		Gtk::CheckButton* everything;
    		
        // main groups:
        Gtk::CheckButton* basic;
        Gtk::CheckButton* luminance;
        Gtk::CheckButton* color;
        Gtk::CheckButton* lens;
        Gtk::CheckButton* composition;
        Gtk::CheckButton* metaicm;
        Gtk::CheckButton* raw;

        // options in basic:
        Gtk::CheckButton* wb;
        Gtk::CheckButton* exposure;
        Gtk::CheckButton* hlrec;

        // options in luminance:
        Gtk::CheckButton* sharpen;
        Gtk::CheckButton* impden;
        Gtk::CheckButton* lumaden;
        Gtk::CheckButton* labcurve;
        Gtk::CheckButton* sh;
        Gtk::CheckButton* dirpyreq;
        Gtk::CheckButton* waveq;

        // options in color:
        Gtk::CheckButton* colormixer;
        Gtk::CheckButton* colorshift;
        Gtk::CheckButton* colorboost;
        Gtk::CheckButton* hsveq;
        Gtk::CheckButton* colorden;
        Gtk::CheckButton* dirpyrden;


        // options in lens:
        Gtk::CheckButton* distortion;
        Gtk::CheckButton* cacorr;
        Gtk::CheckButton* vignetting;

        // options in composition:
        Gtk::CheckButton* coarserot;
        Gtk::CheckButton* finerot;
        Gtk::CheckButton* crop;
        Gtk::CheckButton* resize;

        // options in metaicm:
        Gtk::CheckButton* exifch;
        Gtk::CheckButton* iptc;
        Gtk::CheckButton* icm;

        // options in raw:
        Gtk::CheckButton* df_file;
        Gtk::CheckButton* df_AutoSelect;                                                                             
        Gtk::CheckButton* ff_file;
        Gtk::CheckButton* ff_AutoSelect;
        Gtk::CheckButton* ff_BlurRadius;
        Gtk::CheckButton* ff_BlurType;
                                                                                                      
        sigc::connection everythingConn, basicConn, luminanceConn, colorConn, lensConn, compositionConn, metaicmConn, rawConn;

        sigc::connection wbConn, exposureConn, hlrecConn;
        sigc::connection sharpenConn, impdenConn, lumadenConn, labcurveConn, shConn, dirpyreqConn, waveqConn, hsveqConn;
        sigc::connection colormixerConn, colorshiftConn, colorboostConn, colordenConn, dirpyrdenConn;
        sigc::connection distortionConn, cacorrConn, vignettingConn;
        sigc::connection coarserotConn, finerotConn, cropConn, resizeConn;
        sigc::connection exifchConn, iptcConn, icmConn;
        sigc::connection df_fileConn, df_AutoSelectConn,ff_fileConn, ff_AutoSelectConn, ff_BlurRadiusConn, ff_BlurTypeConn;


        public:
            PartialPasteDlg ();

            void applyPaste (rtengine::procparams::ProcParams* dst, const rtengine::procparams::ProcParams* src);

            void everythingToggled ();
            void basicToggled ();
            void luminanceToggled ();
            void colorToggled ();
            void lensToggled ();
            void compositionToggled ();
            void metaicmToggled ();
            void rawToggled ();
};

#endif

