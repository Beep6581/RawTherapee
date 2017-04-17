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
#include "imagearea.h"
#include "../rtengine/dcrop.h"
#include "../rtengine/refreshmap.h"
#include "../rtengine/rt_math.h"

using namespace rtengine;

CropHandler::CropHandler ()
    : zoom(100), ww(0), wh(0), imx(-1), imy(-1), imw(0), imh(0), cax(-1), cay(-1),
      cx(0), cy(0), cw(0), ch(0), cropX(0), cropY(0), cropW(0), cropH(0), enabled(false),
      cropimg(nullptr), cropimgtrue(nullptr), cropimg_width(0), cropimg_height(0),
      cix(0), ciy(0), ciw(0), cih(0), cis(1),
      initial(false), isLowUpdatePriority(false), ipc(nullptr), crop(nullptr),
      displayHandler(nullptr)
{

    idle_helper = new IdleHelper;
    idle_helper->destroyed = false;
    idle_helper->pending = 0;
    idle_helper->cropHandler = this;
}

CropHandler::~CropHandler ()
{
    idle_register.destroy();

    if (ipc) {
        ipc->delSizeListener (this);
    }

    setEnabled (false);

    if (crop) {
        //crop->destroy ();
        delete crop; // will do the same than destroy, plus delete the object
        crop = nullptr;
    }

    cimg.lock ();

    if (idle_helper->pending) {
        idle_helper->destroyed = true;
    } else {
        delete idle_helper;
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

    EditDataProvider *editDataProvider = nullptr;
    CropWindow *cropWin = displayHandler ? static_cast<CropWindow*>(displayHandler) : nullptr;

    if (cropWin) {
        editDataProvider = cropWin->getImageArea();
    }

    crop = ipc->createCrop (editDataProvider, isDetailWindow);
    ipc->setSizeListener (this);
    crop->setListener (enabled ? this : nullptr);
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

bool CropHandler::isFullDisplay ()
{
    int w, h;
    getFullImageSize(w, h);
    if (!w) {
        return false;
    }
    return cropW == w && cropH == h;
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

    int oldZoom = zoom;
    float oldScale = zoom >= 1000 ? float(zoom / 1000) : 10.f / float(zoom);
    float newScale = z >= 1000 ? float(z / 1000) : 10.f / float(z);

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
        cw = ww * (zoom / 10);
        ch = wh * (zoom / 10);
    }

    cx = cax - cw / 2;
    cy = cay - ch / 2;


    int oldCropX = cropX;
    int oldCropY = cropY;
    int oldCropW = cropW;
    int oldCropH = cropH;

    compDim ();

    if (enabled && (oldZoom != zoom || oldcax != cax || oldcay != cay || oldCropX != cropX || oldCropY != cropY || oldCropW != cropW || oldCropH != cropH)) {
        if (needsFullRefresh) {
            cropPixbuf.clear ();
            ipc->startProcessing(M_HIGHQUAL);
        } else {
            update ();
        }
    }
}

float CropHandler::getZoomFactor ()
{
    if (zoom >= 1000) {
        return zoom / 1000;
    } else {
        return 10.f / (float)zoom;
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
        cw = ww * (zoom / 10);
        ch = wh * (zoom / 10);
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

    cropimg = nullptr;

    if (cropimgtrue) {
        delete [] cropimgtrue;
    }

    cropimgtrue = nullptr;

    if (ax == cropX && ay == cropY && aw == cropW && ah == cropH && askip == (zoom >= 1000 ? 1 : zoom / 10)) {
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
        idle_helper->pending++;

        const auto func = [](gpointer data) -> gboolean {
            IdleHelper* const idle_helper = static_cast<IdleHelper*>(data);

            if (idle_helper->destroyed) {
                if (idle_helper->pending == 1) {
                    delete idle_helper;
                } else {
                    idle_helper->pending--;
                }

                return FALSE;
            }

            CropHandler* ch = idle_helper->cropHandler;

            ch->cimg.lock ();
            ch->cropPixbuf.clear ();

            if (!ch->enabled) {
                delete [] ch->cropimg;
                ch->cropimg = nullptr;
                delete [] ch->cropimgtrue;
                ch->cropimgtrue = nullptr;
                ch->cimg.unlock ();
                return FALSE;
            }

            if (ch->cropimg) {
                if (ch->cix == ch->cropX && ch->ciy == ch->cropY && ch->ciw == ch->cropW && ch->cih == ch->cropH && ch->cis == (ch->zoom >= 1000 ? 1 : ch->zoom / 10)) {
                    // calculate final image size
                    float czoom = ch->zoom >= 1000 ?
                        ch->zoom / 1000.f :
                        float((ch->zoom/10) * 10) / float(ch->zoom);
                    int imw = ch->cropimg_width * czoom;
                    int imh = ch->cropimg_height * czoom;

                    if (imw > ch->ww) {
                        imw = ch->ww;
                    }

                    if (imh > ch->wh) {
                        imh = ch->wh;
                    }

                    Glib::RefPtr<Gdk::Pixbuf> tmpPixbuf = Gdk::Pixbuf::create_from_data (ch->cropimg, Gdk::COLORSPACE_RGB, false, 8, ch->cropimg_width, ch->cropimg_height, 3 * ch->cropimg_width);
                    ch->cropPixbuf = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, imw, imh);
                    tmpPixbuf->scale (ch->cropPixbuf, 0, 0, imw, imh, 0, 0, czoom, czoom, Gdk::INTERP_TILES);
                    tmpPixbuf.clear ();

                    Glib::RefPtr<Gdk::Pixbuf> tmpPixbuftrue = Gdk::Pixbuf::create_from_data (ch->cropimgtrue, Gdk::COLORSPACE_RGB, false, 8, ch->cropimg_width, ch->cropimg_height, 3 * ch->cropimg_width);
                    ch->cropPixbuftrue = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, imw, imh);
                    tmpPixbuftrue->scale (ch->cropPixbuftrue, 0, 0, imw, imh, 0, 0, czoom, czoom, Gdk::INTERP_TILES);
                    tmpPixbuftrue.clear ();
                }

                delete [] ch->cropimg;
                ch->cropimg = nullptr;
                delete [] ch->cropimgtrue;
                ch->cropimgtrue = nullptr;
            }

            ch->cimg.unlock ();

            if (ch->displayHandler) {
                ch->displayHandler->cropImageUpdated ();

                if (ch->initial) {
                    ch->displayHandler->initialImageArrived ();
                    ch->initial = false;
                }
            }

            idle_helper->pending--;

            return FALSE;
        };

        idle_register.add(func, idle_helper);
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

    cskip = zoom >= 1000 ? 1 : zoom/10;

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
            crop->setListener (nullptr);
        }

        cimg.lock ();
        delete [] cropimg;
        cropimg = nullptr;
        delete [] cropimgtrue;
        cropimgtrue = nullptr;
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

void CropHandler::colorPick (const rtengine::Coord &pickerPos, float &r, float &g, float &b, float &rpreview, float &gpreview, float &bpreview, LockableColorPicker::Size size)
{

    if (!cropPixbuf || !cropPixbuftrue) {
        r = g = b = 0.f;
        rpreview = gpreview = bpreview = 0.f;
        return;
    }

    int xSize = (int)size;
    int ySize = (int)size;
    int pixbufW = cropPixbuftrue->get_width();
    int pixbufH = cropPixbuftrue->get_height();
    rtengine::Coord topLeftPos(pickerPos.x - xSize/2, pickerPos.y - ySize/2);

    if (topLeftPos.x > pixbufW || topLeftPos.y > pixbufH || topLeftPos.x + xSize < 0 || topLeftPos.y + ySize < 0) {
        return;
    }

    // Store the position of the center of the picker
    int radius = (int)size / 2;

    // X/Width clip
    if (topLeftPos.x < 0) {
        xSize += topLeftPos.x;
        topLeftPos.x = 0;
    }
    if (topLeftPos.x + xSize > pixbufW) {
        xSize = pixbufW - topLeftPos.x;
    }
    // Y/Height clip
    if (topLeftPos.y < 0) {
        ySize += topLeftPos.y;
        topLeftPos.y = 0;
    }
    if (topLeftPos.y + ySize > pixbufH) {
        ySize = pixbufH - topLeftPos.y;
    }

    // Accumulating the data
    std::uint32_t r2=0, g2=0, b2=0;
    std::uint32_t count = 0;
    const guint8* data = cropPixbuftrue->get_pixels();
    for (int j = topLeftPos.y ; j < topLeftPos.y + ySize ; ++j) {
        const guint8* data2 = data + cropPixbuftrue->get_rowstride()*j;
        for (int i = topLeftPos.x ; i < topLeftPos.x + xSize ; ++i) {
            const guint8* data3 = data2 + i*3;
            rtengine::Coord currPos(i, j);
            rtengine::Coord delta = pickerPos - currPos;
            rtengine::PolarCoord p(delta);
            if (p.radius <= radius) {
                r2 += *data3;
                g2 += *(data3+1);
                b2 += *(data3+2);
                ++count;
            }
        }
    }

    // Averaging
    r = (float)r2 / (float)count / 255.f;
    g = (float)g2 / (float)count / 255.f;
    b = (float)b2 / (float)count / 255.f;

    // Accumulating the data
    r2=0, g2=0, b2=0;
    count = 0;
    data = cropPixbuf->get_pixels();
    for (int j = topLeftPos.y ; j < topLeftPos.y + ySize ; ++j) {
        const guint8* data2 = data + cropPixbuf->get_rowstride()*j;
        for (int i = topLeftPos.x ; i < topLeftPos.x + xSize ; ++i) {
            const guint8* data3 = data2 + i*3;
            rtengine::Coord currPos(i, j);
            rtengine::Coord delta = pickerPos - currPos;
            rtengine::PolarCoord p(delta);
            if (p.radius <= radius) {
                r2 += *data3;
                g2 += *(data3+1);
                b2 += *(data3+2);
                ++count;
            }
        }
    }
    // Averaging
    rpreview = (float)r2 / (float)count / 255.f;
    gpreview = (float)g2 / (float)count / 255.f;
    bpreview = (float)b2 / (float)count / 255.f;
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
        wwImgSpace = int(float(ww) * (float(zoom)/10.f) + 0.5f);
        whImgSpace = int(float(wh) * (float(zoom)/10.f) + 0.5f);
        //scaledFullW = fullW / zoom;
        //scaledFullH = fullH / zoom;
        scaledCAX = int(float(cax) / (float(zoom)/10.f));
        scaledCAY = int(float(cay) / (float(zoom)/10.f));
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
