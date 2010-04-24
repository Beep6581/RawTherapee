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
        // main groups:
        Gtk::CheckButton* basic;
        Gtk::CheckButton* luminance;
        Gtk::CheckButton* color;
        Gtk::CheckButton* lens;
        Gtk::CheckButton* composition;
        Gtk::CheckButton* metaicm;

        // options in basic:
        Gtk::CheckButton* wb;
        Gtk::CheckButton* exposure;
        Gtk::CheckButton* hlrec;

        // options in luminance:
        Gtk::CheckButton* sharpen;
        Gtk::CheckButton* lumaden;
        Gtk::CheckButton* lumacurve;
        Gtk::CheckButton* sh;

        // options in color:
        Gtk::CheckButton* colormixer;
        Gtk::CheckButton* colorshift;
        Gtk::CheckButton* colorboost;
        Gtk::CheckButton* colorden;

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

        sigc::connection basicConn, luminanceConn, colorConn, lensConn, compositionConn, metaicmConn;
        sigc::connection wbConn, exposureConn, hlrecConn;
        sigc::connection sharpenConn, lumadenConn, lumacurveConn, shConn;
        sigc::connection colormixerConn, colorshiftConn, colorboostConn, colordenConn;
        sigc::connection distortionConn, cacorrConn, vignettingConn;
        sigc::connection coarserotConn, finerotConn, cropConn, resizeConn;
        sigc::connection exifchConn, iptcConn, icmConn;


        public:
            PartialPasteDlg ();

            void applyPaste (rtengine::procparams::ProcParams* dst, const rtengine::procparams::ProcParams* src);

            void basicToggled ();
            void luminanceToggled ();
            void colorToggled ();
            void lensToggled ();
            void compositionToggled ();
            void metaicmToggled ();
};

#endif

