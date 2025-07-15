/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2024 Daniel Gao <daniel.gao.work@gmail.com>
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

#include "rtscalable.h"
#include "hidpi.h"

#include <cairomm/pattern.h>
#include <cairomm/surface.h>
#include <gdkmm/pixbuf.h>

#include <utility>

namespace hidpi {

DeviceCoord LogicalCoord::scaleToDevice(int device_scale) const {
    DeviceCoord device = {};
    device.x = x * device_scale;
    device.y = y * device_scale;
    return device;
}

LogicalCoord LogicalCoord::operator+(const LogicalCoord& other) const {
    return LogicalCoord(x + other.x, y + other.y);
}

LogicalCoord LogicalCoord::operator-(const LogicalCoord& other) const {
    return LogicalCoord(x - other.x, y - other.y);
}

DeviceCoord DeviceCoord::operator+(const DeviceSize& other) const {
    return DeviceCoord(x + other.width, y + other.height);
}

DeviceCoord DeviceCoord::operator-(const DeviceSize& other) const {
    return DeviceCoord(x - other.width, y - other.height);
}

LogicalSize LogicalSize::forWidget(const Gtk::Widget* widget) {
    LogicalSize result = {};
    result.width = widget->get_width();
    result.height = widget->get_height();
    return result;
}

ScaledDeviceSize LogicalSize::scaleToDevice(int device_scale) const {
    ScaledDeviceSize result = {};
    result.width = width * device_scale;
    result.height = height * device_scale;
    result.device_scale = device_scale;
    return result;
}

ScaledDeviceSize ScaledDeviceSize::forWidget(const Gtk::Widget* widget) {
    int scale = RTScalable::getScaleForWidget(widget);

    ScaledDeviceSize result = {};
    result.width = widget->get_width() * scale;
    result.height = widget->get_height() * scale;
    result.device_scale = scale;
    return result;
}

DevicePixbuf::DevicePixbuf() : m_device_scale(1) {}

DevicePixbuf::DevicePixbuf(const Glib::RefPtr<Gdk::Pixbuf>& ptr, int device_scale)
        : m_pixbuf(ptr), m_device_scale(device_scale) {}

DevicePixbuf::DevicePixbuf(const DevicePixbuf& other)
        : m_pixbuf(other.m_pixbuf), m_device_scale(other.m_device_scale) {}

DevicePixbuf& DevicePixbuf::operator=(const DevicePixbuf& other) {
    m_pixbuf = other.m_pixbuf;
    m_device_scale = other.m_device_scale;
    return *this;
}

DevicePixbuf::DevicePixbuf(DevicePixbuf&& other) : m_device_scale(1) { swap(*this, other); }

DevicePixbuf& DevicePixbuf::operator=(DevicePixbuf&& other) {
    swap(*this, other);
    return *this;
}

ScaledDeviceSize DevicePixbuf::size() const {
    ScaledDeviceSize size = {};
    if (m_pixbuf) {
        size.width = m_pixbuf->get_width();
        size.height = m_pixbuf->get_height();
    } else {
        size.width = 0;
        size.height = 0;
    }
    size.device_scale = m_device_scale;
    return size;
}

void swap(DevicePixbuf& lhs, DevicePixbuf& rhs) {
    using std::swap;
    swap(lhs.m_pixbuf, rhs.m_pixbuf);
    swap(lhs.m_device_scale, rhs.m_device_scale);
}

Cairo::RefPtr<Cairo::SurfacePattern>
getSourceForSurface(const Cairo::RefPtr<Cairo::Context>& context) {
    Cairo::RefPtr<Cairo::SurfacePattern> result;

    Cairo::RefPtr<Cairo::Pattern> src = context->get_source();
    if (!src) return result;

    auto type = src->get_type();
    if (type != Cairo::PATTERN_TYPE_SURFACE) return result;

    result = Cairo::RefPtr<Cairo::SurfacePattern>::cast_static(src);
    return result;
}

void getDeviceScale(const Cairo::RefPtr<Cairo::Surface>& surface,
                    double& x_scale, double& y_scale) {
    cairo_surface_t* cobj = surface->cobj();
    cairo_surface_get_device_scale(cobj, &x_scale, &y_scale);
}

void setDeviceScale(const Cairo::RefPtr<Cairo::Surface>& surface,
                    int x_scale, int y_scale) {
    cairo_surface_t* cobj = surface->cobj();
    cairo_surface_set_device_scale(cobj, x_scale, y_scale);
}

void getDeviceScale(const Cairo::RefPtr<Cairo::ImageSurface>& surface,
                    double& x_scale, double& y_scale) {
    cairo_surface_t* cobj = surface->cobj();
    cairo_surface_get_device_scale(cobj, &x_scale, &y_scale);
}

void setDeviceScale(const Cairo::RefPtr<Cairo::ImageSurface>& surface,
                    int x_scale, int y_scale) {
    cairo_surface_t* cobj = surface->cobj();
    cairo_surface_set_device_scale(cobj, x_scale, y_scale);
}

}  // namespace hidpi
