/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2026 Daniel Gao <daniel.gao.work@gmail.com>
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

#include "model.h"

#include "cursormanager.h"

#include "rtengine/math/pointvec.h"

using namespace rt;
using namespace rt::canvas;

Session::Session()
    : m_modifiers(GdkModifierType(0)),
      m_cursor_shape(CSArrow),
      m_pan_zoom_flags(PanZoomFlags::ALL)
{
    regenerateTransforms();
}

void Session::setCameraPos(WorldPoint pos)
{
    if (m_camera.pos == pos) return;

    m_camera.pos = pos;
    regenerateTransforms();
    m_events.signal_camera_update.emit();
    queueDraw();
}

void Session::setCameraZoom(double zoom)
{
    if (m_camera.zoom == zoom) return;

    m_camera.zoom = zoom;
    regenerateTransforms();
    m_events.signal_camera_update.emit();
    queueDraw();
}

void Session::setCameraPosZoom(WorldPoint pos, double zoom)
{
    if (m_camera.pos == pos && m_camera.zoom == zoom) return;

    m_camera.pos = pos;
    m_camera.zoom = zoom;
    regenerateTransforms();
    m_events.signal_camera_update.emit();
    queueDraw();
}

void Session::setCameraSize(WidgetSize size)
{
    if (m_camera.size == size) return;

    m_camera.size = size;
    regenerateTransforms();
    m_events.signal_camera_update.emit();
    queueDraw();
}

void Session::setDeviceScale(int device_scale)
{
    if (m_camera.device_scale == device_scale) return;

    m_camera.device_scale = device_scale;
    regenerateTransforms();
    m_events.signal_camera_update.emit();
    queueDraw();
}

void Session::setCameraBounds(const geom::IntBBox& content, CameraBounds bounds)
{
    if (m_bound_mode == bounds) return;
    m_bound_mode = bounds;
    setCamera(content, m_camera);
}

void Session::setCamera(const geom::IntBBox& content, const CameraState& new_state)
{
    auto adjust_to_image = [&]() {
        const double content_width = content.width();
        const double content_height = content.height();
        if (content_width <= 0 || content_height <= 0) {
            return;
        }

        // Prevent image edge from crossing center of camera
        const double min_x = content.min().x;
        const double max_x = content.max().x;
        if (new_state.pos.x.value() < min_x) {
            m_camera.pos.x = WorldScalar(min_x);
        } else if (new_state.pos.x.value() > max_x) {
            m_camera.pos.x = WorldScalar(max_x);
        }

        const double min_y = content.min().y;
        const double max_y = content.max().y;
        if (new_state.pos.y.value() < min_y) {
            m_camera.pos.y = WorldScalar(min_y);
        } else if (new_state.pos.y.value() > max_y) {
            m_camera.pos.y = WorldScalar(max_y);
        }
    };

    auto adjust_to_fill = [&]() {
        const double content_width = content.width();
        const double content_height = content.height();
        if (content_width <= 0 || content_height <= 0) {
            return;
        }

        // Adjust zoom to at least fill the screen
        const double width = new_state.size.width.value() * new_state.device_scale;
        const double height = new_state.size.height.value() * new_state.device_scale;
        const double min_zoom = std::max(width / content_width, height / content_height);
        if (new_state.zoom < min_zoom) {
            m_camera.zoom = min_zoom;
        }

        // Prevent image edge from crossing widget edges
        const double margin_x = width / m_camera.zoom / 2.0;
        const double min_x = content.min().x + margin_x;
        const double max_x = content.max().x - margin_x;
        if (new_state.pos.x.value() < min_x) {
            m_camera.pos.x = WorldScalar(min_x);
        } else if (new_state.pos.x.value() > max_x) {
            m_camera.pos.x = WorldScalar(max_x);
        }

        const double margin_y = height / m_camera.zoom / 2.0;
        const double min_y = content.min().y + margin_y;
        const double max_y = content.max().y - margin_y;
        if (new_state.pos.y.value() < min_y) {
            m_camera.pos.y = WorldScalar(min_y);
        } else if (new_state.pos.y.value() > max_y) {
            m_camera.pos.y = WorldScalar(max_y);
        }
    };

    m_camera = new_state;

    switch (m_bound_mode) {
        case CameraBounds::IMAGE:
        {
            adjust_to_image();
            break;
        }
        case CameraBounds::FILL:
        {
            adjust_to_fill();
            break;
        }
        case CameraBounds::NONE:
        default:
            break;
    }

    regenerateTransforms();
    m_events.signal_camera_update.emit();
    queueDraw();
}

void Session::changeCursorShape(rt::optional<CursorShape> shape)
{
    if (shape) {
        m_cursor_shape = *shape;
    }
    m_events.signal_change_cursor.emit(shape);
}

void Session::zoom11()
{
    setCameraZoom(1.0);
}

void Session::zoomFit(WorldPoint top_left, WorldSize img_size)
{
    WorldPoint center = top_left + (img_size.asVec()) / 2.0;

    double device_scale = m_camera.device_scale;
    double bounds_x = m_camera.size.width.value();
    double bounds_y = m_camera.size.height.value();

    bounds_x = std::max(bounds_x - 20, bounds_x * 0.95);
    bounds_y = std::max(bounds_y - 20, bounds_y * 0.95);

    // Convert to world space without zoom
    bounds_x *= device_scale;
    bounds_y *= device_scale;

    double zoom_x = bounds_x / img_size.width.value();
    double zoom_y = bounds_y / img_size.height.value();
    double zoom = std::min(zoom_x, zoom_y);

    setCameraPosZoom(center, zoom);
}

void Session::onWindowFocusLost(CanvasModel* model)
{
    m_modifiers = GdkModifierType(0);
}

void Session::regenerateTransforms()
{
    m_world_to_widget = SpaceTransform<WorldSpace, WidgetSpace>::build(m_camera);
    m_widget_to_world = SpaceTransform<WidgetSpace, WorldSpace>::build(m_camera);
}

bool ImageModel::isInsideImage(WorldPoint pos) const
{
    geom::BBox bbox(geom::Point(0, 0),
                    geom::Point(m_img_size.width.value(), m_img_size.height.value()));
    return bbox.contains(static_cast<geom::Point>(pos));
}

bool CanvasModel::isCursorInsideImage() const
{
    WorldPoint pos = m_session.widgetToWorldTransform()(m_session.cursorPos());
    return m_image_model.isInsideImage(pos);
}

void CanvasModel::setCameraPos(WorldPoint pos)
{
    CameraState camera = m_session.camera();
    if (camera.pos == pos) return;

    camera.pos = pos;
    m_session.setCamera(buildImageBBox(), camera);
}

void CanvasModel::setCameraZoom(double zoom)
{
    CameraState camera = m_session.camera();
    if (camera.zoom == zoom) return;

    camera.zoom = zoom;
    m_session.setCamera(buildImageBBox(), camera);
}

void CanvasModel::setCameraPosZoom(WorldPoint pos, double zoom)
{
    CameraState camera = m_session.camera();
    if (camera.pos == pos && camera.zoom == zoom) return;

    camera.pos = pos;
    camera.zoom = zoom;
    m_session.setCamera(buildImageBBox(), camera);
}

void CanvasModel::setCameraSize(WidgetSize size)
{
    CameraState camera = m_session.camera();
    if (camera.size == size) return;

    camera.size = size;
    m_session.setCamera(buildImageBBox(), camera);
}

void CanvasModel::setCameraBounds(Session::CameraBounds bounds)
{
    if (m_session.cameraBounds() == bounds) return;
    m_session.setCameraBounds(buildImageBBox(), bounds);
}

geom::IntBBox CanvasModel::buildImageBBox() const
{
    return geom::IntBBox(
        geom::IntPoint(),
        static_cast<geom::IntPoint>(m_image_model.fullSize().asPoint()));
}
