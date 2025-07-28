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
#include "previewhandler.h"
#include <gtkmm.h>
#include "rtengine/rtengine.h"
#include "rtengine/procparams.h"
#include "rtscalable.h"

using namespace rtengine;
using namespace rtengine::procparams;

PreviewHandler::PreviewHandler () :
    image(nullptr),
    cropParams(new procparams::CropParams),
    previewScale(1.)
{

    pih = new PreviewHandlerIdleHelper;
    pih->phandler = this;
    pih->destroyed = false;
    pih->pending = 0;
}

PreviewHandler::~PreviewHandler ()
{
    idle_register.destroy();

    if (pih->pending) {
        pih->destroyed = true;
    } else {
        delete pih;
    }
}

//----------------previewimagelistener functions--------------------

void PreviewHandler::setImage(rtengine::IImage8* i, double scale, const rtengine::procparams::CropParams& cp)
{
    pih->pending++;

    idle_register.add(
        [this, i, scale, cp]() -> bool
        {
            if (pih->destroyed) {
                if (pih->pending == 1) {
                    delete pih;
                } else {
                    --pih->pending;
                }

                return false;
            }

            if (pih->phandler->image) {
                IImage8* const oldImg = pih->phandler->image;

                oldImg->getMutex().lock();
                pih->phandler->image = i;
                oldImg->getMutex().unlock();
            } else {
                pih->phandler->image = i;
            }

            *pih->phandler->cropParams = cp;
            pih->phandler->previewScale = scale;
            --pih->pending;

            return false;
        }
    );
}


void PreviewHandler::delImage(IImage8* i)
{
    pih->pending++;

    idle_register.add(
        [this, i]() -> bool
        {
            if (pih->destroyed) {
                if (pih->pending == 1) {
                    delete pih;
                } else {
                    --pih->pending;
                }

                return false;
            }

            if (pih->phandler->image) {
                IImage8* oldImg = pih->phandler->image;
                oldImg->getMutex().lock();
                pih->phandler->image = nullptr;
                oldImg->getMutex().unlock();
            }

            delete i;
            pih->phandler->previewImgMutex.lock();
            pih->phandler->previewImg.clear();
            pih->phandler->previewImgMutex.unlock();

            --pih->pending;

            return false;
        }
    );
}

void PreviewHandler::imageReady(const rtengine::procparams::CropParams& cp)
{
    pih->pending++;

    idle_register.add(
        [this, cp]() -> bool
        {
            if (pih->destroyed) {
                if (pih->pending == 1) {
                    delete pih;
                } else {
                    --pih->pending;
                }

                return false;
            }

            pih->phandler->previewImgMutex.lock();
            pih->phandler->previewImg = Gdk::Pixbuf::create_from_data(pih->phandler->image->getData(), Gdk::COLORSPACE_RGB, false, 8, pih->phandler->image->getWidth(), pih->phandler->image->getHeight(), 3 * pih->phandler->image->getWidth());
            pih->phandler->previewImgMutex.unlock ();

            *pih->phandler->cropParams = cp;
            pih->phandler->previewImageChanged ();
            --pih->pending;

            return false;
        }
    );
}

Glib::RefPtr<Gdk::Pixbuf> PreviewHandler::getRoughImage (
    ImageCoord pos, hidpi::ScaledDeviceSize desiredSize, double zoom)
{
    MyMutex::MyLock lock(previewImgMutex);

    Glib::RefPtr<Gdk::Pixbuf> resPixbuf;

    if (previewImg) {
        double totalZoom = zoom * previewScale;

        int w = desiredSize.width;
        int h = desiredSize.height;

        if (w > previewImg->get_width()*totalZoom) {
            w = image->getWidth() * totalZoom;
        }

        if (h > previewImg->get_height()*totalZoom) {
            h = image->getHeight() * totalZoom;
        }

        pos.x *= zoom;
        pos.y *= zoom;

        w = rtengine::LIM<int>(w, 0, int(previewImg->get_width() * totalZoom) - pos.x);
        h = rtengine::LIM<int>(h, 0, int(previewImg->get_height() * totalZoom) - pos.y);

        resPixbuf = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, w, h);
        previewImg->scale (resPixbuf, 0, 0, w, h, -pos.x, -pos.y, totalZoom, totalZoom, Gdk::INTERP_NEAREST);
    }

    return resPixbuf;
}

hidpi::DevicePixbuf PreviewHandler::getRoughImage (hidpi::LogicalSize desiredSize,
                                                   int deviceScale, double& outLogicalZoom)
{
    MyMutex::MyLock lock(previewImgMutex);

    hidpi::DevicePixbuf result;

    if (previewImg) {
        double zoom1 = (double)max(desiredSize.width, 20) / previewImg->get_width(); // too small values lead to extremely increased processing time in scale function, Issue 2783
        double zoom2 = (double)max(desiredSize.height, 20) / previewImg->get_height();
        double zoom = zoom1 < zoom2 ? zoom1 : zoom2;

        outLogicalZoom = zoom / previewScale;
        zoom = zoom * deviceScale;

        auto pixbuf = Gdk::Pixbuf::create (Gdk::COLORSPACE_RGB, false, 8, image->getWidth() * zoom, image->getHeight() * zoom);
        previewImg->scale (pixbuf, 0, 0, previewImg->get_width()*zoom, previewImg->get_height()*zoom, 0, 0, zoom, zoom, Gdk::INTERP_BILINEAR);

        result = hidpi::DevicePixbuf(pixbuf, deviceScale);
    }

    return result;
}

void PreviewHandler::previewImageChanged ()
{

    for (std::list<PreviewListener*>::iterator i = listeners.begin(); i != listeners.end(); ++i) {
        (*i)->previewImageChanged ();
    }
}

rtengine::procparams::CropParams PreviewHandler::getCropParams()
{
    return *cropParams;
}
