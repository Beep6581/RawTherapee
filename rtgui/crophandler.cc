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
#include <crophandler.h>

using namespace rtengine;

CropHandler::CropHandler () 
    : crop(NULL), listener(NULL), cropimg(NULL), ipc(NULL),
    cx(0), cy(0), cw(0), ch(0), 
    cropX(0), cropY(0), cropW(0), cropH(0),
    zoom(1000), enabled(false) {
        
    chi = new CropHandlerIdleHelper;
    chi->destroyed = false;
    chi->pending = 0;
    chi->cropHandler = this;
}

CropHandler::~CropHandler () {

    if (ipc)
        ipc->delSizeListener (this);

    setEnabled (false);
    if (crop)
        crop->destroy ();
    cimg.lock ();
    if (chi->pending)
        chi->destroyed = true;
    else
        delete chi;
    cimg.unlock ();
}
        
void CropHandler::newImage (StagedImageProcessor* ipc_) {

    ipc = ipc_;
    cx = 0;
    cy = 0;
    
    if (!ipc)
	    return;
	
	crop = ipc->createCrop ();
    ipc->setSizeListener (this);
    crop->setListener (enabled ? this : NULL);
    initial = true;
}

void CropHandler::sizeChanged (int x, int y, int ow, int oh) {  // the ipc notifies it to keep track size changes like rotation

    compDim ();

// this should be put into an idle source!!!
/*    if (listener)
        listener->cropWindowChanged ();
    */
}

double CropHandler::getFitZoom () {
    
    if (ipc) {
        double z1 = (double) wh / ipc->getFullHeight ();
        double z2 = (double) ww / ipc->getFullWidth ();
        return z1<z2 ? z1 : z2;
    }
    else
        return 1.0;
}

void CropHandler::setZoom (int z, int centerx, int centery) {

    int x = cx + cw / 2;
    int y = cy + ch / 2;

    if (centerx>=0)
            x = centerx;
    if (centery>=0)
            y = centery;
            
    zoom = z;
    if (zoom>=1000) {
        cw = ww * 1000 / zoom;
        ch = wh * 1000 / zoom;
    }
    else {
        cw = ww * zoom;
        ch = wh * zoom;
    }
    cx = x - cw / 2;
    cy = y - ch / 2;

    compDim ();
    if (enabled)
        update ();
}


void CropHandler::setWSize (int w, int h) {

    ww = w;
    wh = h;
    if (zoom>=1000) {
        cw = ww * 1000 / zoom;
        ch = wh * 1000 / zoom;
    }
    else {
        cw = ww * zoom;
        ch = wh * zoom;
    }

    compDim ();
    if (enabled)
        update ();   
}

void CropHandler::getWSize (int& w, int &h) {

    w = ww;
    h = wh;
}

void CropHandler::setPosition (int x, int y, bool update_) {

    cx = x;
    cy = y;

    compDim ();
    if (enabled && update_)
        update ();   
}

void CropHandler::getPosition (int& x, int& y) {

    x = cropX;
    y = cropY;
}


int createpixbufs (void* data) {

    gdk_threads_enter ();

    CropHandlerIdleHelper* chi = (CropHandlerIdleHelper*) data;
    if (chi->destroyed) {
        if (chi->pending == 1)
            delete chi;
        else    
            chi->pending--;
        gdk_threads_leave ();
        return 0;
    }
   
    CropHandler* ch = chi->cropHandler;

    ch->cimg.lock ();
    ch->cropPixbuf.clear ();

    if (!ch->enabled) {
        delete [] ch->cropimg;
        ch->cropimg = NULL;
        ch->cimg.unlock ();
        gdk_threads_leave ();
        return 0;
    }    

    if (ch->cropimg) {
        if (ch->cix==ch->cropX && ch->ciy==ch->cropY && ch->ciw==ch->cropW && ch->cih==ch->cropH && ch->cis==(ch->zoom>=1000?1:ch->zoom)) {
            // calculate final image size
            int czoom = ch->zoom<1000 ? 1000 : ch->zoom;
            int imw = ch->cropimg_width * czoom / 1000;
            int imh = ch->cropimg_height * czoom / 1000;
            if (imw>ch->ww)
                imw = ch->ww;
            if (imh>ch->wh)
                imh = ch->wh;

            Glib::RefPtr<Gdk::Pixbuf> tmpPixbuf = Gdk::Pixbuf::create_from_data (ch->cropimg, Gdk::COLORSPACE_RGB, false, 8, ch->cropimg_width, ch->cropimg_height, 3*ch->cropimg_width);
            ch->cropPixbuf = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, imw, imh);
            tmpPixbuf->scale (ch->cropPixbuf, 0, 0, imw, imh, 0, 0, czoom/1000.0, czoom/1000.0, Gdk::INTERP_NEAREST);
            tmpPixbuf.clear ();
        }
        delete [] ch->cropimg;
        ch->cropimg = NULL;
    }
    ch->cimg.unlock ();
    if (ch->listener) {
        ch->listener->cropImageUpdated ();
        if (ch->initial) {
            ch->listener->initialImageArrived ();
            ch->initial = false;
        }
    }
  
    chi->pending--;
    
    gdk_threads_leave ();
    return 0;
}

void CropHandler::setDetailedCrop (IImage8* im, rtengine::procparams::CropParams cp, int ax, int ay, int aw, int ah, int askip) {

   if (!enabled)
        return;

    cimg.lock ();

    cropParams = cp;

    cropPixbuf.clear ();
    if (cropimg)
        delete [] cropimg;
    cropimg = NULL;
    
    if (ax==cropX && ay==cropY && aw==cropW && ah==cropH && askip==(zoom>=1000?1:zoom)) {
        cropimg_width = im->getWidth ();
        cropimg_height = im->getHeight ();
        cropimg = new unsigned char [3*cropimg_width*cropimg_height];
        memcpy (cropimg, im->getData(), 3*cropimg_width*cropimg_height);
        cix = ax;
        ciy = ay;
        ciw = aw;
        cih = ah;
        cis = askip;
        chi->pending++;
        g_idle_add (createpixbufs, chi);
    }
    cimg.unlock ();
 }

bool CropHandler::getWindow (int& cwx, int& cwy, int& cww, int& cwh, int& cskip) { 
    
    cwx = cropX;
    cwy = cropY;
    cww = cropW;
    cwh = cropH;

    // hack: if called before first size allocation the size will be 0
    if (cww<10)
        cww = 10;
    if (cwh<32)
        cwh = 32;

    cskip = zoom>=1000 ? 1 : zoom;

    return true; 
}

void CropHandler::update () {

    if (crop) {
//        crop->setWindow (cropX, cropY, cropW, cropH, zoom>=1000 ? 1 : zoom); --> we use the "getWindow" hook instead of setting the size before
        crop->setListener (this);
        cropPixbuf.clear ();
        Glib::Thread::create(sigc::mem_fun(*crop, &DetailedCrop::fullUpdate), 0, false, true, Glib::THREAD_PRIORITY_NORMAL);    
    }
}

void CropHandler::setEnabled (bool e) {

    enabled = e;
    if (!enabled) {
        if (crop)
            crop->setListener (NULL);
        cimg.lock ();
        delete [] cropimg;
        cropimg = NULL;
        cropPixbuf.clear ();
        cimg.unlock ();
    }
    else 
        update ();
}

bool CropHandler::getEnabled () {

    return enabled;
}

void CropHandler::getSize (int& w, int& h) {
    
    w = cropW;
    h = cropH;
}

void CropHandler::compDim () {

    cropX = cx;
    cropY = cy;
    cropW = cw;
    cropH = ch;

    cutRectToImgBounds (cropX, cropY, cropW, cropH);
}

void CropHandler::cutRectToImgBounds (int& x, int& y, int& w, int& h) {

    if (ipc) {
        if (w > ipc->getFullWidth())
            w = ipc->getFullWidth();
        if (h > ipc->getFullHeight())
            h = ipc->getFullHeight();
        if (x < 0) 
            x = 0;
        if (y < 0) 
            y = 0;
        if (x + w >= ipc->getFullWidth()) 
            x = ipc->getFullWidth() - w;
        if (y + h >= ipc->getFullHeight()) 
            y = ipc->getFullHeight() - h;
    }
}
