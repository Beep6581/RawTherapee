/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Jean-Christophe FRISCH <natureh.510@gmail.com>
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

#include "rtscalable.h"
#include "../rtengine/icons.h"

double RTScalable::dpi = 0.;
int RTScalable::scale = 0;
Gtk::TextDirection RTScalable::direction = Gtk::TextDirection::TEXT_DIR_NONE;

void RTScalable::setDPInScale (const double newDPI, const int newScale)
{
    if (scale != newScale || (scale == 1 && dpi != newDPI)) {
        // reload all images
        scale = newScale;
        // HOMBRE: On windows, if scale = 2, the dpi is non significant, i.e. should be considered = 192 ; don't know for linux/macos
        dpi = scale == 2 ? 192 : newDPI;
        printf("RTScalable::setDPInScale  /  New scale = %d & new DPI = %.3f (%.3f asked) -> Reloading all RTScalable\n", scale, dpi, newDPI);
    }
}

double RTScalable::getDPI ()
{
    return dpi;
}

int RTScalable::getScale ()
{
    return scale;
}

Gtk::TextDirection RTScalable::getDirection()
{
    return direction;
}

void RTScalable::init(Gtk::Window *window)
{
    printf("RTScalable::init\n");
    setDPInScale(window->get_screen()->get_resolution(), window->get_scale_factor());
}

void RTScalable::resizeImage(Cairo::RefPtr<Cairo::ImageSurface> &surf, double factor)
{
    int newWidth = int((double)surf->get_width() * factor);
    int newHeight = int((double)surf->get_height() * factor);
    printf("Resizing from %dx%d to %dx%d (factor: %.5f)\n", surf->get_width(), surf->get_height(), newWidth, newHeight, factor);
    Cairo::RefPtr<Cairo::ImageSurface> surf2 = Cairo::ImageSurface::create(surf->get_format(), newWidth, newHeight);
    Cairo::RefPtr<Cairo::Context> c = Cairo::Context::create(surf2);
    c->scale (factor, factor);
    c->set_source(surf, 0., 0.);
    //Cairo::RefPtr<Cairo::SurfacePattern> p = Cairo::SurfacePattern::create(surf);
    //p->set_filter (Cairo::FILTER_BILINEAR);
    c->paint ();
    surf2->flush();

    surf = surf2;
}
