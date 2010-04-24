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
#include <guiutils.h>
#include <options.h>
#include <utils.h>

bool removeIfThere (Gtk::Container* cont, Gtk::Widget* w, bool increference) {

    Glib::ListHandle<Gtk::Widget*> list = cont->get_children ();
    Glib::ListHandle<Gtk::Widget*>::iterator i = list.begin ();
    for (; i!=list.end() && *i!=w; i++);
    if (i!=list.end()) {
        if (increference)
            w->reference ();
        cont->remove (*w);
        return true;
    }
    else
        return false;
}

void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh) {

    if (options.thumbInterp==0)
        rtengine::nearestInterp (src, sw, sh, dst, dw, dh);
    else if (options.thumbInterp==1)
        rtengine::bilinearInterp (src, sw, sh, dst, dw, dh);
}

Glib::ustring removeExtension (const Glib::ustring& filename) {

    Glib::ustring bname = Glib::path_get_basename(filename);
    int lastdot = bname.find_last_of ('.');
    if (lastdot!=bname.npos)
        return filename.substr (0, filename.size()-(bname.size()-lastdot));
    else
        return filename;
}

Glib::ustring getExtension (const Glib::ustring& filename) {

    Glib::ustring bname = Glib::path_get_basename(filename);
    int lastdot = bname.find_last_of ('.');
    if (lastdot!=bname.npos)
        return filename.substr (filename.size()-(bname.size()-lastdot)+1, filename.npos);
    else
        return "";
}

void drawCrop (Cairo::RefPtr<Cairo::Context> cr, int imx, int imy, int imw, int imh, int startx, int starty, double scale, const rtengine::procparams::CropParams& cparams) {

    cr->set_line_width (1.0);
    cr->rectangle (imx+0.5, imy+0.5, imw, imh);
    cr->clip ();

    double c1x = (cparams.x-startx)*scale;
    double c1y = (cparams.y-starty)*scale;
    double c2x = (cparams.x+cparams.w-1-startx)*scale;
    double c2y = (cparams.y+cparams.h-1-starty)*scale;

    cr->set_source_rgba (0, 0, 0, 2.0/3.0);
    cr->rectangle (imx+0.5, imy+0.5, imw, c1y);
    cr->rectangle (imx+0.5, imy+0.5+c2y, imw, imh-c2y);
    cr->rectangle (imx+0.5, imy+0.5+c1y, c1x, c2y-c1y+1);
    cr->rectangle (imx+0.5+c2x, imy+0.5+c1y, imw-c2x, c2y-c1y+1);
    cr->fill ();

    // rectangle around the cropped area and guides
    if (cparams.guide!="None") {
        double rectx1 = c1x + imx + 0.5;
        double recty1 = c1y + imy + 0.5;
        double rectx2 = c2x + imx + 0.5;
        double recty2 = c2y + imy + 0.5;
        cr->set_line_width (1.0);
        cr->set_source_rgb (1.0, 1.0, 1.0);
        cr->move_to (rectx1, recty1);
        cr->line_to (rectx2, recty1);
        cr->line_to (rectx2, recty2);
        cr->line_to (rectx1, recty2);
        cr->line_to (rectx1, recty1);
        cr->stroke ();
        cr->set_source_rgb (0.0, 0.0, 0.0);
        std::valarray<double> ds (1);
        ds[0] = 4;
        cr->set_dash (ds, 0);
        cr->move_to (rectx1, recty1);
        cr->line_to (rectx2, recty1);
        cr->line_to (rectx2, recty2);
        cr->line_to (rectx1, recty2);
        cr->line_to (rectx1, recty1);
        cr->stroke ();
        ds.resize (0);
        cr->set_dash (ds, 0);
        
        if (cparams.guide!="Rule of diagonals") {
            // draw guide lines
            std::vector<double> horiz_ratios;
            std::vector<double> vert_ratios;
            
            if (cparams.guide=="Rule of thirds") {
                horiz_ratios.push_back (1.0/3.0);
                horiz_ratios.push_back (2.0/3.0);
                vert_ratios.push_back (1.0/3.0);
                vert_ratios.push_back (2.0/3.0);
            }
            else if (cparams.guide=="Harmonic means 1") {
                horiz_ratios.push_back (1.0-0.618);
                vert_ratios.push_back (1.0-0.618);
            }
            else if (cparams.guide=="Harmonic means 2") {
                horiz_ratios.push_back (0.618);
                vert_ratios.push_back (1.0-0.618);
            }
            else if (cparams.guide=="Harmonic means 3") {
                horiz_ratios.push_back (1.0-0.618);
                vert_ratios.push_back (0.618);
            }
            else if (cparams.guide=="Harmonic means 4") {
                horiz_ratios.push_back (0.618);
                vert_ratios.push_back (0.618);
            }
            for (int i=0; i<vert_ratios.size(); i++) {
                cr->set_source_rgb (1.0, 1.0, 1.0);
                cr->move_to (rectx1 + (rectx2-rectx1) * vert_ratios[i], recty1);
                cr->line_to (rectx1 + (rectx2-rectx1) * vert_ratios[i], recty2);
                cr->move_to (rectx1, recty1 + (recty2-recty1) * horiz_ratios[i]);
                cr->line_to (rectx2, recty1 + (recty2-recty1) * horiz_ratios[i]);
                cr->stroke ();
                cr->set_source_rgb (0.0, 0.0, 0.0);
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_dash (ds, 0);
                cr->move_to (rectx1 + (rectx2-rectx1) * vert_ratios[i], recty1);
                cr->line_to (rectx1 + (rectx2-rectx1) * vert_ratios[i], recty2);
                cr->move_to (rectx1, recty1 + (recty2-recty1) * horiz_ratios[i]);
                cr->line_to (rectx2, recty1 + (recty2-recty1) * horiz_ratios[i]);
                cr->stroke ();
                ds.resize (0);
                cr->set_dash (ds, 0);
            }           
        }
        else {
            int corners_from[4][2];
            int corners_to[4][2];
            int mindim = MIN (rectx2-rectx1+1, recty2-recty1+1);
            corners_from[0][0] = rectx1;
            corners_from[0][1] = recty1;
            corners_to[0][0]   = rectx1 + mindim;
            corners_to[0][1]   = recty1 + mindim;
            corners_from[1][0] = rectx1;
            corners_from[1][1] = recty2;
            corners_to[1][0]   = rectx1 + mindim;
            corners_to[1][1]   = recty2 - mindim;
            corners_from[2][0] = rectx2;
            corners_from[2][1] = recty1;
            corners_to[2][0]   = rectx2 - mindim;
            corners_to[2][1]   = recty1 + mindim;
            corners_from[3][0] = rectx2;
            corners_from[3][1] = recty2;
            corners_to[3][0]   = rectx2 - mindim;
            corners_to[3][1]   = recty2 - mindim;
            for (int i=0; i<4; i++) {
                cr->set_source_rgb (1.0, 1.0, 1.0);
                cr->move_to (corners_from[i][0], corners_from[i][1]);
                cr->line_to (corners_to[i][0], corners_to[i][1]);
                cr->stroke ();
                cr->set_source_rgb (0.0, 0.0, 0.0);
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_dash (ds, 0);
                cr->move_to (corners_from[i][0], corners_from[i][1]);
                cr->line_to (corners_to[i][0], corners_to[i][1]);
                cr->stroke ();
                ds.resize (0);
                cr->set_dash (ds, 0);
            }
        }
    }
    cr->reset_clip ();
}
