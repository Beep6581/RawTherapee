/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2018 Jean-Christophe FRISCH <natureh.510@gmail.com>
 *  Copyright (c) 2022 Pierre CABRERA <pierre.cab@gmail.com>
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

#include "rtimage.h"

#include "rtsurface.h"

std::map<std::pair<Glib::ustring, Gtk::IconSize>, std::shared_ptr<RTSurface>> RTImageCache::cache;

std::shared_ptr<RTSurface> RTImageCache::getCachedSurface(const Glib::ustring &icon_name, const Gtk::IconSize icon_size)
{
    // Look for an existing cached icon
    const auto key = std::pair<Glib::ustring, Gtk::IconSize>(icon_name, icon_size);
    const auto item = cache.find(key);

    if (item != cache.end()) { // A cached icon exists
        return item->second;
    } else { // Create the icon
        auto surface = std::shared_ptr<RTSurface>(new RTSurface(icon_name, icon_size));

        // Add the surface to the cache if the icon exist
        if (surface) {
            cache.insert({key, surface});
        }

        return surface;
    }
}

void RTImageCache::updateCache()
{
    // Iterate over cache to updated RTSurface
    for (auto const& item : cache) {
        item.second->updateSurface();
    }
}

RTImage::RTImage () {}

RTImage::RTImage (const Glib::ustring& iconName, const Gtk::IconSize iconSize) :
    sigc::trackable(),
    Glib::ObjectBase(),
    Gtk::Image(),
    size(iconSize),
    icon_name(iconName),
    g_icon(Glib::RefPtr<const Gio::Icon>())
{
    // Set surface from icon cache
    surface = RTImageCache::getCachedSurface(this->icon_name, this->size);

    // Add it to the RTImage if surface exists
    if (surface) {
        set(surface->get());
    }
}

RTImage::RTImage (const Glib::RefPtr<const Gio::Icon>& gIcon, const Gtk::IconSize iconSize) :
    Gtk::Image(),
    size(iconSize),
    icon_name(""),
    g_icon(gIcon)
{
    // Configure RTImage based on g_icon
    set(this->g_icon, this->size);
}

void RTImage::set_from_icon_name(const Glib::ustring& iconName)
{
    set_from_icon_name(iconName, this->size);
}

void RTImage::set_from_icon_name(const Glib::ustring& iconName, const Gtk::IconSize iconSize)
{
    this->icon_name = iconName;
    this->size = iconSize;

    // Set surface from icon cache
    surface = RTImageCache::getCachedSurface(this->icon_name, this->size);

    // Add it to the RTImage if previously chosen
    if (surface) {
        set(surface->get());
    }

    // Unset Gio::Icon if previously chosen
    if (this->g_icon) {
        g_icon = Glib::RefPtr<const Gio::Icon>();
    }
}

void RTImage::set_from_gicon(const Glib::RefPtr<const Gio::Icon>& gIcon)
{
    set_from_gicon(gIcon, this->size);
}

void RTImage::set_from_gicon(const Glib::RefPtr<const Gio::Icon>& gIcon, const Gtk::IconSize iconSize)
{
    this->g_icon = gIcon;
    this->size = iconSize;

    // Set image from Gio::Icon
    set(this->g_icon, this->size);

    // Unset surface if previously chosen
    this->icon_name = "";

    if (surface) {
        surface = std::shared_ptr<RTSurface>();
    }
}

int RTImage::get_width()
{
    if (surface) {
        return surface->getWidth();
    } else if (g_icon) {
        Gtk::Image::get_width();
    }

    return -1;
}

int RTImage::get_height()
{
    if (surface) {
        return surface->getHeight();
    } else if (g_icon) {
        Gtk::Image::get_height();
    }

    return -1;
}
