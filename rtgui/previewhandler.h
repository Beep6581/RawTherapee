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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <list>
#include <memory>

#include <gtkmm.h>

#include "guiutils.h"
#include "threadutils.h"

#include "rtengine/noncopyable.h"
#include "rtengine/rtengine.h"

class PreviewListener
{
public:
    virtual ~PreviewListener() = default;
    virtual void previewImageChanged() = 0;
};

class PreviewHandler;

struct PreviewHandlerIdleHelper {
    PreviewHandler* phandler;
    bool destroyed;
    int pending;
};

class PreviewHandler final : public rtengine::PreviewImageListener, public rtengine::NonCopyable
{
private:
    friend int setImageUI   (void* data);
    friend int delImageUI   (void* data);
    friend int imageReadyUI (void* data);

    IdleRegister idle_register;

protected:
    rtengine::IImage8* image;
    const std::unique_ptr<rtengine::procparams::CropParams> cropParams;
    double previewScale;
    PreviewHandlerIdleHelper* pih;
    std::list<PreviewListener*> listeners;
    MyMutex previewImgMutex;
    Glib::RefPtr<Gdk::Pixbuf> previewImg;

public:

    PreviewHandler ();
    ~PreviewHandler () override;

    void addPreviewImageListener (PreviewListener* l)
    {
        listeners.push_back (l);
    }

    // previewimagelistener
    void setImage(rtengine::IImage8* img, double scale, const rtengine::procparams::CropParams& cp) override;
    void delImage(rtengine::IImage8* img) override;
    void imageReady(const rtengine::procparams::CropParams& cp) override;

    // this function is called when a new preview image arrives from rtengine
    void previewImageChanged ();

    // with this function it is possible to ask for a rough approximation of a (possibly zoomed) crop of the image
    Glib::RefPtr<Gdk::Pixbuf>           getRoughImage (int x, int y, int w, int h, double zoom);
    Glib::RefPtr<Gdk::Pixbuf>           getRoughImage (int desiredW, int desiredH, double& zoom);
    rtengine::procparams::CropParams    getCropParams ();
};
