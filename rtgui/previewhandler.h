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
#ifndef _PREVIEWHANDLER_
#define _PREVIEWHANDLER_

#include <list>

#include <gtkmm.h>

#include "threadutils.h"
#include "guiutils.h"

#include "../rtengine/rtengine.h"

class PreviewListener
{

public:
    virtual ~PreviewListener () {}
    virtual void previewImageChanged () {}
};

class PreviewHandler;
struct PreviewHandlerIdleHelper {
    PreviewHandler* phandler;
    bool destroyed;
    int pending;
};

class PreviewHandler : public rtengine::PreviewImageListener
{
private:
    friend int setImageUI   (void* data);
    friend int delImageUI   (void* data);
    friend int imageReadyUI (void* data);

    IdleRegister idle_register;

protected:
    rtengine::IImage8* image;
    rtengine::procparams::CropParams cropParams;
    double previewScale;
    PreviewHandlerIdleHelper* pih;
    std::list<PreviewListener*> listeners;
    MyMutex previewImgMutex;
    Glib::RefPtr<Gdk::Pixbuf> previewImg;

public:

    PreviewHandler ();
    virtual ~PreviewHandler ();

    void addPreviewImageListener (PreviewListener* l)
    {
        listeners.push_back (l);
    }

    // previewimagelistener
    void setImage   (rtengine::IImage8* img, double scale, rtengine::procparams::CropParams cp);
    void delImage   (rtengine::IImage8* img);
    void imageReady (rtengine::procparams::CropParams cp);

    // this function is called when a new preview image arrives from rtengine
    void previewImageChanged ();

    // with this function it is possible to ask for a rough approximation of a (possibly zoomed) crop of the image
    Glib::RefPtr<Gdk::Pixbuf>           getRoughImage (int x, int y, int w, int h, double zoom);
    Glib::RefPtr<Gdk::Pixbuf>           getRoughImage (int desiredW, int desiredH, double& zoom);
    rtengine::procparams::CropParams    getCropParams ()
    {
        return cropParams;
    }
};

#endif
