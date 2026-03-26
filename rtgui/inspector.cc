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

#include "inspector.h"

#include "canvas/canvas.h"
#include "canvas/model.h"
#include "canvas/render.h"
#include "guiutils.h"
#include "multilangmgr.h"
#include "options.h"
#include "pathutils.h"
#include "rtscalable.h"

#include "rtengine/previewimage.h"
#include "rtengine/rt_math.h"
#include "rtengine/rtapp.h"
#include "rtengine/util/cpp.h"

using namespace rt::canvas;

constexpr bool NO_PADDING = false;

struct InspectorBuffer
{
    Glib::ustring filepath;
    Cairo::RefPtr<Cairo::ImageSurface> surface;

    InspectorBuffer(const Glib::ustring& image_path,
                    const Cairo::RefPtr<Cairo::ImageSurface>& image_surface)
        : filepath(image_path), surface(image_surface)
    {
    }
};

Inspector::Inspector()
    : m_canvas_model(rt::make_unique<CanvasModel>()),
      m_renderer(rt::make_unique<InspectorRenderer>()),
      m_curr_image(nullptr),
      m_is_active(false),
      m_is_pinned(false),
      m_fit_to_screen(false),
      m_is_initialized(false),
      m_is_window_fullscreen(false),
      m_is_window_showing(false)
{
    set_name("Inspector");

    m_canvas = rt::make_managed<Canvas>(m_canvas_model.get());
    m_canvas->setRenderer(m_renderer.get());
    onPreferencesChanged();  // Configure pan zoom based on options
    m_canvas->signal_widget_size_update.connect(
        sigc::mem_fun(*this, &Inspector::onCanvasSizeChanged));
    pack_start(*m_canvas, true, true);

    m_canvas_model->session().setCameraBounds(Session::CameraBounds::FILL_OR_FIT);

    const auto& options = App::get().options();
    if (options.inspectorWindow) {
        m_window = rt::make_unique<Gtk::Window>();

        m_window->set_name("InspectorWindow");
        m_window->set_title("RawTherapee " + M("INSPECTOR_WINDOW_TITLE"));
        m_window->set_visible(false);

        m_window->signal_hide().connect(sigc::mem_fun(*this, &Inspector::onWindowHide));
        m_window->signal_window_state_event().connect(
            sigc::mem_fun(*this, &Inspector::onWindowStateEvent));
        m_window->signal_focus_out_event().connect(
            sigc::mem_fun(*this, &Inspector::onWindowFocusOut));

        m_window->add_events(Gdk::BUTTON_PRESS_MASK | Gdk::KEY_PRESS_MASK
                             | Gdk::FOCUS_CHANGE_MASK);

        m_click_controller = rt::gtk4::GestureClick::create(m_window.get());
        m_key_controller = rt::gtk4::EventControllerKey::create(m_window.get());

        m_click_controller->signal_pressed().connect(
            sigc::mem_fun(*this, &Inspector::onButtonPressed));
        m_key_controller->set_propagation_phase(GTK_PHASE_CAPTURE);
        m_key_controller->signal_key_pressed().connect(
            sigc::mem_fun(*this, &Inspector::onKeyPressed));
        m_key_controller->signal_key_released().connect(
            sigc::mem_fun(*this, &Inspector::onKeyReleased));

        m_window->add(*this);
        m_window->set_size_request(500, 500);
        m_window->fullscreen();

        m_canvas->enablePanZoom(true);

        m_is_initialized = false;  // Delay init to avoid flickering on some systems
        m_is_active = true;  // Always track inspected thumbnails
    } else {
        m_renderer->setDrawFrame(true);
    }

    App::get().signal_preferences_changed().connect(
        sigc::mem_fun(*this, &Inspector::onPreferencesChanged));
}

Inspector::~Inspector() = default;

void Inspector::showWindow(bool pinned, bool scaled)
{
    if (!m_window || m_is_window_showing) return;

    if (!m_is_initialized) {
        m_window->show_all();
        m_is_initialized = true;
    }

    // The window must be set to visible before calling switchImage() otherwise
    // an internal check ignores the update...
    m_window->set_visible(true);
    m_is_window_showing = true;

    m_is_pinned = pinned;
    m_fit_to_screen = scaled;

    // Update content when becoming visible
    switchImage(m_next_image_path);
    mouseMove(m_next_image_pos);
}

void Inspector::onButtonPressed(int n_press, double x, double y)
{
    if (!m_window) return;

    if (!m_is_pinned) {
        // Pin window with mouse click
        m_is_pinned = true;
    }
}

bool Inspector::onKeyPressed(guint keyval, guint keycode, GdkModifierType state)
{
    if (!m_window) return false;

    switch (keyval) {
        case GDK_KEY_z:
        case GDK_KEY_F:
            if (m_is_pinned) {
                m_fit_to_screen = false;
                if (m_zoomed_pos) {
                    m_canvas_model->setCameraPosZoom(*m_zoomed_pos, 1.0);
                } else {
                    m_canvas_model->setCameraZoom(1.0);
                }
                m_canvas_model->session().queueDraw();
            }
            return true;
        case GDK_KEY_f:
            if (m_is_pinned) {
                m_fit_to_screen = true;
                m_zoomed_pos = m_canvas_model->session().camera().pos;
                m_canvas_model->zoomFit(NO_PADDING);
                m_canvas_model->session().queueDraw();
            }
            return true;
        case GDK_KEY_F11:
            // Toggle fullscreen
            if (m_is_window_fullscreen) {
                m_window->unfullscreen();
            } else {
                m_window->fullscreen();
            }
            m_is_window_fullscreen = !m_is_window_fullscreen;
            return true;
        case GDK_KEY_Escape:
            // Hide window
            if (m_is_pinned) {
                m_zoomed_pos = m_canvas_model->session().camera().pos;
            }
            m_is_pinned = false;
            m_window->set_visible(false);
            clearImage();
            return true;
    }

    if (m_is_pinned) {
        return m_canvas->onKeyPressed(keyval, keycode, state);
    }

    return false;
}

void Inspector::onKeyReleased(guint keyval, guint keycode, GdkModifierType state)
{
    if (!m_window) return;

    if (!m_is_pinned) {
        switch (keyval) {
            case GDK_KEY_f:
            case GDK_KEY_F:
                m_window->set_visible(false);
                clearImage();
            default:
                break;
        }
    } else {
        m_canvas->onKeyReleased(keyval, keycode, state);
    }
}

bool Inspector::onWindowStateEvent(GdkEventWindowState* event)
{
    if (!m_window->get_window() || m_window->get_window()->gobj() != event->window) {
        return false;
    }

    m_is_window_fullscreen = event->new_window_state & GDK_WINDOW_STATE_FULLSCREEN;

    return true;
}

bool Inspector::onWindowFocusOut(GdkEventFocus* event)
{
    m_canvas_model->session().onWindowFocusLost(m_canvas_model.get());
    return false;
}

void Inspector::onCanvasSizeChanged()
{
    m_canvas_model->zoomFit(NO_PADDING);
}

void Inspector::onPreferencesChanged()
{
    const auto& options = App::get().options();

    m_canvas->setScrollMode(options.zoomOnScroll ? ScrollMode::ZOOM : ScrollMode::PAN);
    m_canvas->setSmoothScrollDirection(
        options.reverseScrollDir
        ? ScrollDirection::REVERSE : ScrollDirection::NATURAL);
    m_canvas->setSmoothScrollSensitivity(
        static_cast<double>(options.smoothScrollSensitivity)
        / Options::SMOOTH_SCROLL_SENSITIVITY_FACTOR);
    m_canvas->setSmoothScrollPanSensitivity(
        static_cast<double>(options.smoothScrollPanSensitivity)
        / Options::SMOOTH_SCROLL_PAN_SENSITIVITY_FACTOR);
}

void Inspector::mouseMove(rtengine::Coord2D pos)
{
    if (!m_is_active) return;

    m_next_image_pos = pos;

    // Skip actual update of content when not visible
    if (m_window && !m_window->get_visible()) return;

    if (m_fit_to_screen) return;

    if (m_curr_image) {
        double x = static_cast<double>(m_curr_image->surface->get_width());
        double y = static_cast<double>(m_curr_image->surface->get_height());

        x *= rtengine::LIM01(pos.x);
        y *= rtengine::LIM01(pos.y);

        m_canvas_model->setCameraPos(WorldPoint{WorldScalar(x), WorldScalar(y)});
        m_canvas_model->session().queueDraw();
    }
}

void Inspector::switchImage(const Glib::ustring& full_path)
{
    if (!m_is_active) return;

    if (m_delay_connection.connected()) {
        m_delay_connection.disconnect();
    }

    m_next_image_path = full_path;

    // Skip actual update of content when not visible
    if (m_window && !m_window->get_visible()) return;

    const auto& options = App::get().options();
    if (!options.inspectorDelay) {
        doSwitchImage();
    } else {
        m_delay_connection = Glib::signal_timeout().connect(
            sigc::mem_fun(*this, &Inspector::doSwitchImage), options.inspectorDelay);
    }
}


bool Inspector::doSwitchImage()
{
    // Update buffer cache based on preferences
    const size_t max_cache_size = []() {
        int val = App::get().options().maxInspectorBuffers;
        return val > 0 ? static_cast<size_t>(val) : 1;
    }();

    if (m_images.size() > max_cache_size) {
        const size_t num_dropped = m_images.size() - max_cache_size;

        // Drop oldest entries
        for (size_t i = 0; i < num_dropped; i++) {
            if (m_curr_image == m_images.at(i).get()) {
                // Shouldn't be possible but guard for memory safety
                changeCurrImage(nullptr);
            }
            m_images.at(i) = std::move(m_images.at(num_dropped + i));
        }
        m_images.resize(max_cache_size);
    }

    if (m_next_image_path.empty()) return true;

    for (size_t i = 0; i < m_images.size(); ++i) {
        if (!m_images[i] || m_images[i]->filepath != m_next_image_path) continue;

        // Rotate towards front by 1
        std::unique_ptr<InspectorBuffer> tmp = std::move(m_images[i]);
        for (size_t j = i; j < m_images.size() - 1; ++j) {
            m_images.at(j) = std::move(m_images.at(j + 1));
        }

        // Move the last used image to the tail
        changeCurrImage(tmp.get());
        m_images.back() = std::move(tmp);
        return true;
    }

    // Loading a new image
    Glib::ustring ext = getExtension(m_next_image_path);
    if (ext.empty()) {
        changeCurrImage(nullptr);
        m_next_image_path.clear();
        return true;
    }

    rtengine::PreviewImage pi(m_next_image_path, ext,
                              rtengine::PreviewImage::PIM_EmbeddedOrRaw);
    Cairo::RefPtr<Cairo::ImageSurface> surface = pi.getImage();
    if (!surface) {
        changeCurrImage(nullptr);
        m_next_image_path.clear();
        return true;
    }

    // Add loaded image to tail
    auto buffer = rt::make_unique<InspectorBuffer>(m_next_image_path, surface);
    changeCurrImage(buffer.get());
    if (m_images.size() == max_cache_size) {
        m_images.erase(m_images.begin());  // Delete the oldest entry
    }
    m_images.emplace_back(std::move(buffer));
    m_next_image_path.clear();
    return true;
}

void Inspector::changeCurrImage(InspectorBuffer* buffer)
{
    if (m_curr_image == buffer) {
        if (!m_canvas_model->image().imageSurface()) {
            showImage();
        }
        return;
    }

    m_curr_image = buffer;
    m_zoomed_pos = rt::nullopt;

    if (!buffer || !buffer->surface) {
        clearImage();
        return;
    }

    showImage();
}

void Inspector::showImage()
{
    if (!m_curr_image || !m_curr_image->surface) return;

    IntWorldSize img_size;
    img_size.width = IntWorldScalar(m_curr_image->surface->get_width());
    img_size.height = IntWorldScalar(m_curr_image->surface->get_height());
    m_canvas_model->image().setImageSurface(m_curr_image->surface, img_size);

    if (m_fit_to_screen) {
        m_canvas_model->zoomFit(NO_PADDING);
    } else {
        double x = static_cast<double>(m_curr_image->surface->get_width());
        double y = static_cast<double>(m_curr_image->surface->get_height());

        x *= rtengine::LIM01(m_next_image_pos.x);
        y *= rtengine::LIM01(m_next_image_pos.y);

        WorldPoint new_pos{WorldScalar(x), WorldScalar(y)};
        m_canvas_model->setCameraPosZoom(new_pos, 1.0);
    }

    m_canvas_model->session().queueDraw();
}

void Inspector::clearImage()
{
    m_canvas_model->image().setImageSurface(
        Cairo::RefPtr<Cairo::ImageSurface>{}, IntWorldSize{});
    m_canvas_model->setCameraPosZoom(WorldPoint{}, 1.0);
    m_canvas_model->session().queueDraw();
}

void Inspector::flushBuffers()
{
    changeCurrImage(nullptr);
    m_images.clear();
}

void Inspector::setActive(bool state)
{
    if (!state) {
        flushBuffers();
    }

    if (!m_window) {
        m_is_active = state;
    }
}
