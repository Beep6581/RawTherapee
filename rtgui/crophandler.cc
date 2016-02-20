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
#include "crophandler.h"
#undef THREAD_PRIORITY_NORMAL

#include <cstring>
#include "guiutils.h"
#include "cropwindow.h"
#include "../rtengine/dcrop.h"
#include "../rtengine/refreshmap.h"
#include "../rtengine/rt_math.h"

using namespace rtengine;

CropHandler::CropHandler ()
    : zoom(10), ww(0), wh(0), imx(-1), imy(-1), imw(0), imh(0), cax(-1), cay(-1),
      cx(0), cy(0), cw(0), ch(0), cropX(0), cropY(0), cropW(0), cropH(0), enabled(false),
      cropimg(NULL), cropimgtrue(NULL), cropimg_width(0), cropimg_height(0),
      initial(false), isLowUpdatePriority(false), ipc(NULL), crop(NULL),
      displayHandler(NULL)
{

    chi = new CropHandlerIdleHelper;
    chi->destroyed = false;
    chi->pending = 0;
    chi->cropHandler = this;
}

CropHandler::~CropHandler ()
{

    if (ipc) {
        ipc->delSizeListener (this);
    }

    setEnabled (false);

    if (crop) {
        //crop->destroy ();
        delete crop; // will do the same than destroy, plus delete the object
        crop = NULL;
    }

    cimg.lock ();

    if (chi->pending) {
        chi->destroyed = true;
    } else {
        delete chi;
    }

    cimg.unlock ();
}

void CropHandler::setEditSubscriber (EditSubscriber* newSubscriber)
{
    (static_cast<rtengine::Crop *>(crop))->setEditSubscriber(newSubscriber);
}

void CropHandler::newImage (StagedImageProcessor* ipc_, bool isDetailWindow)
{

    ipc = ipc_;
    cx = 0;
    cy = 0;

    if (!ipc) {
        return;
    }

    EditDataProvider *editDataProvider = NULL;
    CropWindow *cropWin = displayHandler ? static_cast<CropWindow*>(displayHandler) : NULL;

    if (cropWin) {
        editDataProvider = cropWin->getImageArea();
    }

    crop = ipc->createCrop (editDataProvider, isDetailWindow);
    ipc->setSizeListener (this);
    crop->setListener (enabled ? this : NULL);
    initial = true;
}

void CropHandler::sizeChanged (int x, int y, int ow, int oh)    // the ipc notifies it to keep track size changes like rotation
{

    compDim ();

// this should be put into an idle source!!!
    /*    if (listener)
            listener->cropWindowChanged ();
        */
}

double CropHandler::getFitCropZoom ()
{
    double z1 = (double) wh / cropParams.h;
    double z2 = (double) ww / cropParams.w;
    return z1 < z2 ? z1 : z2;
}

double CropHandler::getFitZoom ()
{

    if (ipc) {
        double z1 = (double) wh / ipc->getFullHeight ();
        double z2 = (double) ww / ipc->getFullWidth ();
        return z1 < z2 ? z1 : z2;
    } else {
        return 1.0;
    }
}

void CropHandler::setZoom (int z, int centerx, int centery)
{
    assert (ipc);

    float oldScale = zoom >= 1000 ? float(zoom / 1000) : 1.f / float(zoom);
    float newScale = z >= 1000 ? float(z / 1000) : 1.f / float(z);

    int oldcax = cax;
    int oldcay = cay;

    if (centerx == -1) {
        cax = ipc->getFullWidth () / 2;
    } else {
        float distToAnchor = float(cax - centerx);
        distToAnchor = distToAnchor / newScale * oldScale;
        cax = centerx + int(distToAnchor);
    }

    if (centery == -1) {
        cay = ipc->getFullHeight () / 2;
    } else {
        float distToAnchor = float(cay - centery);
        distToAnchor = distToAnchor / newScale * oldScale;
        cay = centery + int(distToAnchor);
    }

    // maybe demosaic etc. if we cross the border to >100%
    bool needsFullRefresh = (z >= 1000 && zoom < 1000);

    zoom = z;

    if (zoom >= 1000) {
        cw = ww * 1000 / zoom;
        ch = wh * 1000 / zoom;
    } else {
        cw = ww * zoom;
        ch = wh * zoom;
    }

    cx = cax - cw / 2;
    cy = cay - ch / 2;


    int oldCropX = cropX;
    int oldCropY = cropY;
    int oldCropW = cropW;
    int oldCropH = cropH;

    compDim ();

    if (enabled && (oldcax != cax || oldcay != cay || oldCropX != cropX || oldCropY != cropY || oldCropW != cropW || oldCropH != cropH)) {
        if (needsFullRefresh) {
            ipc->startProcessing(M_HIGHQUAL);
        } else {
            update ();
        }
    }
}


void CropHandler::setWSize (int w, int h)
{

    ww = w;
    wh = h;

    if (zoom >= 1000) {
        cw = ww * 1000 / zoom;
        ch = wh * 1000 / zoom;
    } else {
        cw = ww * zoom;
        ch = wh * zoom;
    }

    compDim ();

    if (enabled) {
        update ();
    }
}

void CropHandler::getWSize (int& w, int &h)
{

    w = ww;
    h = wh;
}

void CropHandler::getAnchorPosition (int& x, int& y)
{
    x = cax;
    y = cay;
}

void CropHandler::setAnchorPosition (int x, int y, bool update_)
{
    cax = x;
    cay = y;

    compDim ();

    if (enabled && update_) {
        update ();
    }
}

void CropHandler::moveAnchor (int deltaX, int deltaY, bool update_)
{
    cax += deltaX;
    cay += deltaY;

    compDim ();

    if (enabled && update_) {
        update ();
    }
}

void CropHandler::centerAnchor (bool update_)
{
    assert (ipc);

    // Computes the crop's size and position given the anchor's position and display size

    cax = ipc->getFullWidth() / 2;
    cay = ipc->getFullHeight() / 2;

    compDim ();

    if (enabled && update_) {
        update ();
    }
}

void CropHandler::getPosition (int& x, int& y)
{

    x = cropX;
    y = cropY;
}


int createpixbufs (void* data)
{

    CropHandlerIdleHelper* chi = static_cast<CropHandlerIdleHelper*>(data);

    if (chi->destroyed) {
        if (chi->pending == 1) {
            delete chi;
        } else {
            chi->pending--;
        }

        return 0;
    }

    CropHandler* ch = chi->cropHandler;

    ch->cimg.lock ();
    ch->cropPixbuf.clear ();

    if (!ch->enabled) {
        delete [] ch->cropimg;
        ch->cropimg = NULL;
        delete [] ch->cropimgtrue;
        ch->cropimgtrue = NULL;
        ch->cimg.unlock ();
        return 0;
    }

    if (ch->cropimg) {
        if (ch->cix == ch->cropX && ch->ciy == ch->cropY && ch->ciw == ch->cropW && ch->cih == ch->cropH && ch->cis == (ch->zoom >= 1000 ? 1 : ch->zoom)) {
            // calculate final image size
            int czoom = ch->zoom < 1000 ? 1000 : ch->zoom;
            int imw = ch->cropimg_width * czoom / 1000;
            int imh = ch->cropimg_height * czoom / 1000;

            if (imw > ch->ww) {
                imw = ch->ww;
            }

            if (imh > ch->wh) {
                imh = ch->wh;
            }

            Glib::RefPtr<Gdk::Pixbuf> tmpPixbuf = Gdk::Pixbuf::create_from_data (ch->cropimg, Gdk::COLORSPACE_RGB, false, 8, ch->cropimg_width, 2 * ch->cropimg_height, 3 * ch->cropimg_width);
            ch->cropPixbuf = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, imw, imh);
            tmpPixbuf->scale (ch->cropPixbuf, 0, 0, imw, imh, 0, 0, czoom / 1000.0, czoom / 1000.0, Gdk::INTERP_NEAREST);
            tmpPixbuf.clear ();

            Glib::RefPtr<Gdk::Pixbuf> tmpPixbuftrue = Gdk::Pixbuf::create_from_data (ch->cropimgtrue, Gdk::COLORSPACE_RGB, false, 8, ch->cropimg_width, 2 * ch->cropimg_height, 3 * ch->cropimg_width);
            ch->cropPixbuftrue = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, imw, imh);
            tmpPixbuftrue->scale (ch->cropPixbuftrue, 0, 0, imw, imh, 0, 0, czoom / 1000.0, czoom / 1000.0, Gdk::INTERP_NEAREST);
            tmpPixbuftrue.clear ();
        }

        delete [] ch->cropimg;
        ch->cropimg = NULL;
        delete [] ch->cropimgtrue;
        ch->cropimgtrue = NULL;
    }

    ch->cimg.unlock ();

    if (ch->displayHandler) {
        ch->displayHandler->cropImageUpdated ();

        if (ch->initial) {
            ch->displayHandler->initialImageArrived ();
            ch->initial = false;
        }
    }

    chi->pending--;

    return 0;
}

void CropHandler::setDetailedCrop (IImage8* im, IImage8* imtrue, rtengine::procparams::ColorManagementParams cmp,
                                   rtengine::procparams::CropParams cp, int ax, int ay, int aw, int ah, int askip)
{

    if (!enabled) {
        return;
    }

    cimg.lock ();

    cropParams = cp;
    colorParams = cmp;

    cropPixbuf.clear ();

    if (cropimg) {
        delete [] cropimg;
    }

    cropimg = NULL;

    if (cropimgtrue) {
        delete [] cropimgtrue;
    }

    cropimgtrue = NULL;

    if (ax == cropX && ay == cropY && aw == cropW && ah == cropH && askip == (zoom >= 1000 ? 1 : zoom)) {
        cropimg_width = im->getWidth ();
        cropimg_height = im->getHeight ();
        cropimg = new unsigned char [3 * cropimg_width * cropimg_height];
        cropimgtrue = new unsigned char [3 * cropimg_width * cropimg_height];
        memcpy (cropimg, im->getData(), 3 * cropimg_width * cropimg_height);
        memcpy (cropimgtrue, imtrue->getData(), 3 * cropimg_width * cropimg_height);
        cix = ax;
        ciy = ay;
        ciw = aw;
        cih = ah;
        cis = askip;
        chi->pending++;
        add_idle (createpixbufs, chi);
    }

    cimg.unlock ();
}

bool CropHandler::getWindow (int& cwx, int& cwy, int& cww, int& cwh, int& cskip)
{

    cwx = cropX;
    cwy = cropY;
    cww = cropW;
    cwh = cropH;

    // hack: if called before first size allocation the size will be 0
    if (cww < 10) {
        cww = 10;
    }

    if (cwh < 32) {
        cwh = 32;
    }

    cskip = zoom >= 1000 ? 1 : zoom;

    return true;
}

void CropHandler::update ()
{

    if (crop && enabled) {
//        crop->setWindow (cropX, cropY, cropW, cropH, zoom>=1000 ? 1 : zoom); --> we use the "getWindow" hook instead of setting the size before
        crop->setListener (this);
        cropPixbuf.clear ();

        // To save threads, try to mark "needUpdate" without a thread first
        if (crop->tryUpdate()) {
            if (isLowUpdatePriority) {
                Glib::Thread::create(sigc::mem_fun(*crop, &DetailedCrop::fullUpdate), 0, false, true, Glib::THREAD_PRIORITY_LOW);
            } else {
                Glib::Thread::create(sigc::mem_fun(*crop, &DetailedCrop::fullUpdate), false );
            }
        }
    }
}

void CropHandler::setEnabled (bool e)
{

    enabled = e;

    if (!enabled) {
        if (crop) {
            crop->setListener (NULL);
        }

        cimg.lock ();
        delete [] cropimg;
        cropimg = NULL;
        delete [] cropimgtrue;
        cropimgtrue = NULL;
        cropPixbuf.clear ();
        cimg.unlock ();
    } else {
        update ();
    }
}

bool CropHandler::getEnabled ()
{

    return enabled;
}

void CropHandler::getSize (int& w, int& h)
{

    w = cropW;
    h = cropH;
}

void CropHandler::getFullImageSize (int& w, int& h)
{
    if (ipc) {
        w = ipc->getFullWidth ();
        h = ipc->getFullHeight ();
    } else {
        w = h = 0;
    }
}

void CropHandler::compDim ()
{
    assert (ipc && displayHandler);

    // Computes the crop's size and position given the anchor's position and display size

    int fullW = ipc->getFullWidth();
    int fullH = ipc->getFullHeight();
    int imgX = -1, imgY = -1;
    //int scaledFullW, scaledFullH;
    int scaledCAX, scaledCAY;
    int wwImgSpace;
    int whImgSpace;

    cax = rtengine::LIM(cax, 0, fullW-1);
    cay = rtengine::LIM(cay, 0, fullH-1);

    if (zoom >= 1000) {
        wwImgSpace = int(float(ww) / float(zoom/1000) + 0.5f);
        whImgSpace = int(float(wh) / float(zoom/1000) + 0.5f);
        //scaledFullW = fullW * (zoom/1000);
        //scaledFullH = fullH * (zoom/1000);
        scaledCAX = cax * (zoom/1000);
        scaledCAY = cay * (zoom/1000);
    } else {
        wwImgSpace = int(float(ww) * float(zoom) + 0.5f);
        whImgSpace = int(float(wh) * float(zoom) + 0.5f);
        //scaledFullW = fullW / zoom;
        //scaledFullH = fullH / zoom;
        scaledCAX = cax / zoom;
        scaledCAY = cay / zoom;
    }

    imgX = ww / 2 - scaledCAX;
    if (imgX < 0) {
        imgX = 0;
    }
    imgY = wh / 2 - scaledCAY;
    if (imgY < 0) {
        imgY = 0;
    }

    cropX = cax - (wwImgSpace/2);
    cropY = cay - (whImgSpace/2);
    cropW = wwImgSpace;
    cropH = whImgSpace;

    if (cropX + cropW > fullW) {
        cropW = fullW - cropX;
    }

    if (cropY + cropH > fullH) {
        cropH = fullH - cropY;
    }

    if (cropX < 0) {
        cropW += cropX;
        cropX = 0;
    }

    if (cropY < 0) {
        cropH += cropY;
        cropY = 0;
    }

    // Should be good already, but this will correct eventual rounding error

    if (cropW > fullW) {
        cropW = fullW;
    }

    if (cropH > fullH) {
        cropH = fullH;
    }

    displayHandler->setDisplayPosition(imgX, imgY);
}
