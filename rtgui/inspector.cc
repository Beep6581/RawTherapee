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
#include "inspector.h"
#include "guiutils.h"
#include <gtkmm.h>
#include "cursormanager.h"
#include "guiutils.h"
#include "options.h"
#include "../rtengine/safegtk.h"
#include "../rtengine/previewimage.h"

extern Options options;

InspectorBuffer::InspectorBuffer(const Glib::ustring &imagePath) : currTransform(0), fromRaw(false)
{
    if (!imagePath.empty() && safe_file_test(imagePath, Glib::FILE_TEST_EXISTS) && !safe_file_test(imagePath, Glib::FILE_TEST_IS_DIR)) {
        imgPath = imagePath;

        // generate thumbnail image
        Glib::ustring ext = getExtension (imagePath);

        if (ext == "") {
            imgPath.clear();
            return;
        }

        rtengine::PreviewImage pi(imagePath, ext, rtengine::PreviewImage::PIM_EmbeddedOrRaw);
        Cairo::RefPtr<Cairo::ImageSurface> imageSurface = pi.getImage();

        if (imageSurface) {
            imgBuffer.setSurface(imageSurface);
            fromRaw = true;
        } else {
            imgPath.clear();
        }
    }
}

/*
InspectorBuffer::~InspectorBuffer() {
}
*/

int InspectorBuffer::infoFromImage (const Glib::ustring& fname)
{

    rtengine::ImageMetaData* idata = rtengine::ImageMetaData::fromFile (fname, NULL);

    if (!idata) {
        return 0;
    }

    int deg = 0;

    if (idata->hasExif()) {
        if      (idata->getOrientation() == "Rotate 90 CW" ) {
            deg = 90;
        } else if (idata->getOrientation() == "Rotate 180"   ) {
            deg = 180;
        } else if (idata->getOrientation() == "Rotate 270 CW") {
            deg = 270;
        }
    }

    delete idata;
    return deg;
}

Inspector::Inspector () : currImage(NULL), zoom(0.0), active(false) {}

Inspector::~Inspector()
{
    deleteBuffers();
}

bool Inspector::on_expose_event (GdkEventExpose* event)
{

    Glib::RefPtr<Gdk::Window> win = get_window();

    if (!win) {
        return false;
    }

    if (!active) {
        active = true;
    }


    // cleanup the region


    if (currImage && currImage->imgBuffer.surfaceCreated()) {
        // this will eventually create/update the off-screen pixmap

        // compute the displayed area
        rtengine::Coord availableSize;
        rtengine::Coord topLeft;
        rtengine::Coord displayedSize;
        rtengine::Coord dest(0, 0);
        win->get_size(availableSize.x, availableSize.y);
        int imW = currImage->imgBuffer.getWidth();
        int imH = currImage->imgBuffer.getHeight();

        if (imW < availableSize.x) {
            // center the image in the available space along X
            topLeft.x = 0;
            displayedSize.x = availableSize.x;
            dest.x = (availableSize.x - imW) / 2;
        } else {
            // partial image display
            // double clamp
            topLeft.x = center.x + availableSize.x / 2;
            topLeft.x = rtengine::min<int>(topLeft.x, imW);
            topLeft.x -= availableSize.x;
            topLeft.x = rtengine::max<int>(topLeft.x, 0);
        }

        if (imH < availableSize.y) {
            // center the image in the available space along Y
            topLeft.y = 0;
            displayedSize.y = availableSize.y;
            dest.y = (availableSize.y - imH) / 2;
        } else {
            // partial image display
            // double clamp
            topLeft.y = center.y + availableSize.y / 2;
            topLeft.y = rtengine::min<int>(topLeft.y, imH);
            topLeft.y -= availableSize.y;
            topLeft.y = rtengine::max<int>(topLeft.y, 0);
        }

        //printf("center: %d, %d   (img: %d, %d)  (availableSize: %d, %d)  (topLeft: %d, %d)\n", center.x, center.y, imW, imH, availableSize.x, availableSize.y, topLeft.x, topLeft.y);

        // define the destination area
        currImage->imgBuffer.setDrawRectangle(win, dest.x, dest.y, rtengine::min<int>(availableSize.x - dest.x, imW), rtengine::min<int>(availableSize.y - dest.y, imH), false);
        currImage->imgBuffer.setSrcOffset(topLeft.x, topLeft.y);

        if (!currImage->imgBuffer.surfaceCreated()) {
            return false;
        }

        // Draw!

        Gdk::Color c;
        Cairo::RefPtr<Cairo::Context> cr = win->create_cairo_context();
        Glib::RefPtr<Gtk::Style> style = get_style();

        // draw the background
        c = style->get_bg (Gtk::STATE_NORMAL);
        cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
        cr->set_line_width (0);
        cr->rectangle (0, 0, availableSize.x, availableSize.y);
        cr->fill ();

        currImage->imgBuffer.copySurface(win);

        // draw the frame
        c = style->get_fg (Gtk::STATE_NORMAL);
        cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
        cr->set_line_width (1);
        cr->rectangle (0.5, 0.5, availableSize.x - 1, availableSize.y - 1);
        cr->stroke ();
    }

    return true;
}

void Inspector::mouseMove (rtengine::Coord2D pos, int transform)
{
    if (!active) {
        return;
    }

    if (currImage) {
        center.set(int(rtengine::LIM01(pos.x)*double(currImage->imgBuffer.getWidth())), int(rtengine::LIM01(pos.y)*double(currImage->imgBuffer.getHeight())));
    } else {
        center.set(0, 0);
    }

    queue_draw();
}

void Inspector::switchImage (const Glib::ustring &fullPath)
{
    if (!active) {
        return;
    }

    // we first check the size of the list, it may have been changed in Preference
    if (images.size() > size_t(options.maxInspectorBuffers)) {
        // deleting the last entries
        for (size_t i = images.size() - 1; i > size_t(options.maxInspectorBuffers - 1); --i) {
            delete images.at(i);
            images.at(i) = NULL;
        }

        // resizing down
        images.resize(options.maxInspectorBuffers);
    }

    if (fullPath.empty()) {
        currImage = NULL;
        queue_draw();
        return;
    } else {
        bool found = false;

        for (size_t i = 0; i < images.size(); ++i) {
            if (images.at(i) != NULL && images.at(i)->imgPath == fullPath) {
                currImage = images.at(i);

                // rolling the list 1 step to the beginning
                for (size_t j = i; j < images.size() - 1; ++j) {
                    images.at(j) = images.at(j + 1);
                }

                images.at(images.size() - 1) = currImage; // move the last used image to the tail
                found = true;
                break;
            }
        }

        if (!found) {
            if (images.size() == size_t(options.maxInspectorBuffers)) {
                // The list is full, delete the first entry
                delete images.at(0);
                images.erase(images.begin());
            }

            // Loading a new image
            InspectorBuffer *iBuffer = new InspectorBuffer(fullPath);

            // and add it to the tail
            if (!iBuffer->imgPath.empty()) {
                images.push_back(iBuffer);
                currImage = images.at(images.size() - 1);
            } else {
                currImage = NULL;
            }
        }
    }
}

void Inspector::deleteBuffers ()
{
    for (size_t i = 0; i < images.size(); ++i) {
        if (images.at(i) != NULL) {
            delete images.at(i);
            images.at(i) = NULL;
        }
    }

    images.resize(0);
    currImage = NULL;
}

void Inspector::flushBuffers ()
{
    if (!active) {
        return;
    }

    deleteBuffers();
}

void Inspector::setActive(bool state)
{
    if (!state) {
        flushBuffers();
    }

    active = state;
}

