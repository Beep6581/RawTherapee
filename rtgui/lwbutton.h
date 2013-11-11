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
#ifndef _LWBUTTON_
#define _LWBUTTON_

#include <gtkmm.h>

class LWButton;
class LWButtonListener {

    public:
        virtual ~LWButtonListener () {}
        virtual void buttonPressed (LWButton* button, int actionCode, void* actionData) {}
        virtual void redrawNeeded  (LWButton* button) {}
};

class LWButton {
    
    public:
        enum Alignment {Left, Right, Top, Bottom, Center};
        enum State { Normal, Over, Pressed_In, Pressed_Out};
        
    private:
        int xpos, ypos, w, h;
        Alignment halign, valign;
        Cairo::RefPtr<Cairo::ImageSurface> icon;
        double bgr, bgg, bgb;
        double fgr, fgg, fgb;
        State state;
        LWButtonListener* listener;
        int actionCode;
        void* actionData;
        Glib::ustring toolTip;
               
    public:
        LWButton (Cairo::RefPtr<Cairo::ImageSurface> i, int aCode, void* aData, Alignment ha=Left, Alignment va=Center, Glib::ustring tooltip="");
    
        void    getSize             (int& minw, int& minh);
        void    getAlignment        (Alignment& ha, Alignment& va);
        void    setPosition         (int x, int y);
        void    getPosition         (int& x, int& y);
        bool    inside              (int x, int y);
        void    setIcon             (Cairo::RefPtr<Cairo::ImageSurface> i);
        Cairo::RefPtr<Cairo::ImageSurface> getIcon ();
        void    setColors           (const Gdk::Color& bg, const Gdk::Color& fg);
        void    setToolTip          (const Glib::ustring& tooltip);    
    
        bool    motionNotify        (int x, int y);
        bool    pressNotify         (int x, int y);
        bool    releaseNotify       (int x, int y);

        Glib::ustring getToolTip (int x, int y);
        
        void    setButtonListener   (LWButtonListener* bl) { listener = bl; }

        void    redraw              (Cairo::RefPtr<Cairo::Context> context);
};

#endif
