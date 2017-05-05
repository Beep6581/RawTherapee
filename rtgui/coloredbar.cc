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

#include "coloredbar.h"
#include "../rtengine/utils.h"

ColoredBar::ColoredBar (eRTOrientation orient)
{
    orientation = orient;
    dirty = true;
    this->x = this->y = this->w = this->h = 0;
}

bool ColoredBar::setDrawRectangle(int newX, int newY, int newW, int newH, bool updateBackBufferSize)
{
    return BackBuffer::setDrawRectangle(Cairo::FORMAT_ARGB32, newX, newY, newW, newH, updateBackBufferSize);
}

/*
 * Redraw the bar to a Cairo::ImageSurface
 */
void ColoredBar::expose(Gtk::DrawingArea &drawingArea, Cairo::RefPtr<Cairo::ImageSurface> destSurface)
{
    // look out if the Surface has to be redrawn
    if (!surfaceCreated() || !destSurface) {
        return;
    }

    updateBackBuffer(drawingArea);
    Gdk::Rectangle rect(x, y, w, h);
    copySurface(destSurface, &rect);
}

/*
 * Redraw the bar to a Gdk::Window
 */
void ColoredBar::expose(Gtk::DrawingArea &drawingArea, Glib::RefPtr<Gdk::Window> destWindow)
{
    // look out if the Surface has to be redrawn
    if (!surfaceCreated() || !destWindow) {
        return;
    }

    updateBackBuffer(drawingArea);
    Gdk::Rectangle rect(x, y, w, h);
    copySurface(destWindow, &rect);
}

void ColoredBar::expose(Gtk::DrawingArea &drawingArea, const Cairo::RefPtr< Cairo::Context> &cr)
{
    // look out if the Surface has to be redrawn
    if (!surfaceCreated()) {
        return;
    }

    updateBackBuffer(drawingArea);
    Gdk::Rectangle rect(x, y, w, h);
    copySurface(cr, &rect);
}

/*
 * Redraw the bar to a Gdk::Window
 */
void ColoredBar::expose(Gtk::DrawingArea &drawingArea, BackBuffer *backBuffer)
{
    // look out if the Surface has to be redrawn
    if (!surfaceCreated() || !backBuffer) {
        return;
    }

    updateBackBuffer(drawingArea);
    Gdk::Rectangle rect(x, y, w, h);
    copySurface(backBuffer, &rect);
}

void ColoredBar::updateBackBuffer(Gtk::DrawingArea &drawingArea)
{
    if (isDirty()) {
        Cairo::RefPtr<Cairo::Context> cr = getContext();

        // the bar has to be drawn to the Surface first
        if (!bgGradient.empty()) {
            // a gradient has been set, we use it
            cr->set_line_width(0.);

            // gradient background
            Cairo::RefPtr< Cairo::LinearGradient > bggradient;

            switch (orientation) {
            case (RTO_Left2Right):
                bggradient = Cairo::LinearGradient::create (0., 0., double(w), 0.);
                break;

            case (RTO_Right2Left):
                bggradient = Cairo::LinearGradient::create (double(w), 0., 0., 0.);
                break;

            case (RTO_Bottom2Top):
                bggradient = Cairo::LinearGradient::create (0., double(h), 0., 0.);
                break;

            case (RTO_Top2Bottom):
            default:
                bggradient = Cairo::LinearGradient::create (0., 0., 0., double(h));
                break;
            }

            for (std::vector<GradientMilestone>::iterator i = bgGradient.begin(); i != bgGradient.end(); ++i) {
                bggradient->add_color_stop_rgb (i->position, i->r, i->g, i->b);
            }

            cr->set_source (bggradient);
            cr->rectangle(0, 0, w, h);
            cr->fill();
        } else {
            // ask the ColorProvider to provide colors :) for each pixels
            if (colorProvider) {
                surface->flush();
                
                unsigned char *surfaceData = surface->get_data();

                cr->set_antialias(Cairo::ANTIALIAS_NONE);
                cr->set_line_width(1.);

                switch (orientation) {
                case (RTO_Left2Right):
                    for (int py = 0; py < h; ++py) {
                        for (int px = 0; px < w; ++px) {
                        	unsigned char *pixel = surfaceData + (py * w + px) * 4;
                            double x_ = double(      px);
                            //double y_ = double((h-1)-py);  unused
                            double x01 = x_        / double(w - 1);
                            double y01 = double(py) / double(h - 1);
                            colorProvider->colorForValue (x01, y01, CCET_BACKGROUND, colorCallerId, this);

                            rtengine::poke01_d(pixel, ccRed, ccGreen, ccBlue);
                        }
                    }

                    break;

                case (RTO_Right2Left):
                    for (int py = 0; py < h; ++py) {
                        for (int px = 0; px < w; ++px) {
                        	unsigned char *pixel = surfaceData + (py * w + px) * 4;
                            //double x_ = double((w-1)-px);  unused
                            //double y_ = double((h-1)-py);  unused
                            double x01 = double(px) / double(w - 1);
                            double y01 = double(py) / double(h - 1);
                            colorProvider->colorForValue (x01, y01, CCET_BACKGROUND, colorCallerId, this);

                            rtengine::poke01_d(pixel, ccRed, ccGreen, ccBlue);
                        }
                    }

                    break;

                case (RTO_Bottom2Top):
                    for (int py = 0; py < h; ++py) {
                        for (int px = 0; px < w; ++px) {
                        	unsigned char *pixel = surfaceData + (py * w + px) * 4;
                            //double x_ = double((w-1)-px);  unused
                            //double y_ = double((h-1)-py);  unused
                            double x01 = double(px) / double(w - 1);
                            double y01 = double(py) / double(h - 1);
                            colorProvider->colorForValue (y01, x01, CCET_BACKGROUND, colorCallerId, this);

                            rtengine::poke01_d(pixel, ccRed, ccGreen, ccBlue);
                        }
                    }

                    break;

                case (RTO_Top2Bottom):
                default:
                    for (int py = 0; py < h; ++py) {
                        for (int px = 0; px < w; ++px) {
                        	unsigned char *pixel = surfaceData + (py * w + px) * 4;
                            double x_ = double(      px);
                            double y_ = double(      py);
                            double x01 = x_ / double(w - 1);
                            double y01 = y_ / double(h - 1);
                            colorProvider->colorForValue (y01, x01, CCET_BACKGROUND, colorCallerId, this);

                            rtengine::poke01_d(pixel, ccRed, ccGreen, ccBlue);
                        }
                    }

                    break;
                }

                surface->mark_dirty();
            }
        }

        // has it been updated or not, we assume that the Surface has been correctly set (we don't handle allocation error)
        setDirty(false);
    }
}

void ColoredBar::setBgGradient (const std::vector<GradientMilestone> &milestones)
{
    bgGradient = milestones;
    setDirty(true);
}

void ColoredBar::clearBgGradient ()
{
    bgGradient.clear();
    setDirty(true);
}

bool ColoredBar::canGetColors()
{
    return colorProvider != nullptr || bgGradient.size() > 0;
}
