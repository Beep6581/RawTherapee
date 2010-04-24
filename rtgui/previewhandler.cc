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
#include <previewhandler.h>
#include <gtkmm.h>
#include <rtengine.h>

using namespace rtengine;
using namespace rtengine::procparams;

PreviewHandler::PreviewHandler () : image(NULL) {
    
    pih = new PreviewHandlerIdleHelper;
    pih->phandler = this;
    pih->destroyed = false;
    pih->pending = 0;
};

PreviewHandler::~PreviewHandler () {

    if (pih->pending)
        pih->destroyed = true;
    else
        delete pih;
};

//----------------previewimagelistener functions--------------------

struct iaimgpar {
    IImage8* image;
    PreviewHandlerIdleHelper* pih;
    double scale;
    CropParams cp;
};

int iasetimage (void* data) {

    gdk_threads_enter ();

    iaimgpar* iap = (iaimgpar*)data;
    PreviewHandlerIdleHelper* pih = iap->pih;

    if (pih->destroyed) {
        if (pih->pending == 1)
            delete pih;
        else    
            pih->pending--;
        delete iap;
        gdk_threads_leave ();
        return 0;
    }

    if (pih->phandler->image) {
        IImage8* temp = pih->phandler->image;
        temp->getMutex().lock ();
        pih->phandler->image = iap->image;
        temp->getMutex().unlock ();
    }
    else    
        pih->phandler->image = iap->image;
    pih->phandler->cropParams = iap->cp;
    pih->phandler->previewScale = iap->scale;
    pih->pending--;
    delete iap;

    gdk_threads_leave ();
    
    return 0;
}

void PreviewHandler::setImage (rtengine::IImage8* i, double scale, rtengine::procparams::CropParams cp) { 

    pih->pending++;

    iaimgpar* iap = new iaimgpar;
    iap->image      = i;
    iap->pih        = pih;
    iap->scale      = scale;
    iap->cp         = cp;

    g_idle_add (iasetimage, iap);
}

int iadelimage (void* data) {

    gdk_threads_enter ();

    iaimgpar* iap = (iaimgpar*)data;
    PreviewHandlerIdleHelper* pih = iap->pih;

    if (pih->destroyed) {
        if (pih->pending == 1)
            delete pih;
        else    
            pih->pending--;
        delete iap;
        gdk_threads_leave ();
        return 0;
    }
    
    if (pih->phandler->image) {
        IImage8* temp = pih->phandler->image;
        temp->getMutex().lock ();
        pih->phandler->image = NULL;
        temp->getMutex().unlock ();
    }    
    iap->image->free ();
    pih->phandler->previewImgMutex.lock ();
    pih->phandler->previewImg.clear ();
    pih->phandler->previewImgMutex.unlock ();

    pih->pending--;
    delete iap;

    gdk_threads_leave ();

    return 0;
}

void PreviewHandler::delImage (IImage8* i) { 

    pih->pending++;

    iaimgpar* iap = new iaimgpar;
    iap->image    = i;
    iap->pih = pih;

    g_idle_add (iadelimage, iap);
}

int imready (void* data) {

    gdk_threads_enter ();
    
    iaimgpar* iap = (iaimgpar*)data;
    PreviewHandlerIdleHelper* pih = iap->pih;
    
    if (pih->destroyed) {
        if (pih->pending == 1)
            delete pih;
        else    
            pih->pending--;
        delete iap;
        gdk_threads_leave ();
        return 0;
    }

    pih->phandler->previewImgMutex.lock ();
    pih->phandler->previewImg = Gdk::Pixbuf::create_from_data (pih->phandler->image->getData(), Gdk::COLORSPACE_RGB, false, 8, pih->phandler->image->getWidth(), pih->phandler->image->getHeight(), 3*pih->phandler->image->getWidth());    
    pih->phandler->previewImgMutex.unlock ();
    pih->phandler->cropParams = iap->cp;
    pih->phandler->previewImageChanged ();
    pih->pending--;
    delete iap;

    gdk_threads_leave ();

    return 0;
}

void PreviewHandler::imageReady (CropParams cp) {

    pih->pending++;
    iaimgpar* iap = new iaimgpar;
    iap->pih      = pih;
    iap->cp       = cp;
	g_idle_add (imready, iap);
}

Glib::RefPtr<Gdk::Pixbuf> PreviewHandler::getRoughImage (int x, int y, int w, int h, double zoom) {

    Glib::RefPtr<Gdk::Pixbuf> resPixbuf;
    previewImgMutex.lock ();
    if (previewImg) {
        double totalZoom = zoom*previewScale;
        if (w>previewImg->get_width()*totalZoom)
            w = image->getWidth()*totalZoom;
        if (h>previewImg->get_height()*totalZoom)
            h = image->getHeight()*totalZoom;           
        int ix = x*zoom;
        int iy = y*zoom;
        if (ix<0)
            ix = 0;
        if (iy<0)
            iy = 0;
        if ((ix+w)/totalZoom>previewImg->get_width())
            ix = previewImg->get_width()*totalZoom - w;
        if ((iy+h)/totalZoom>previewImg->get_height())
            iy = previewImg->get_height()*totalZoom - h;

        resPixbuf = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, w, h);
        previewImg->scale (resPixbuf, 0, 0, w, h, -ix, -iy, totalZoom, totalZoom, Gdk::INTERP_NEAREST);
    }
    previewImgMutex.unlock ();
    return resPixbuf;
}

Glib::RefPtr<Gdk::Pixbuf> PreviewHandler::getRoughImage (int desiredW, int desiredH, double& zoom_) {

    Glib::RefPtr<Gdk::Pixbuf> resPixbuf;
    previewImgMutex.lock ();
    if (previewImg) {
        double zoom1 = (double)desiredW / previewImg->get_width();
        double zoom2 = (double)desiredH / previewImg->get_height();
        double zoom = zoom1<zoom2 ? zoom1 : zoom2;
        
        resPixbuf = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, image->getWidth()*zoom, image->getHeight()*zoom);
        previewImg->scale (resPixbuf, 0, 0, previewImg->get_width()*zoom, previewImg->get_height()*zoom, 0, 0, zoom, zoom, Gdk::INTERP_NEAREST);
        zoom_ = zoom / previewScale;
    }
    previewImgMutex.unlock ();
    return resPixbuf;
}

void PreviewHandler::previewImageChanged () {
    
    for (std::list<PreviewListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
        (*i)->previewImageChanged ();
}
