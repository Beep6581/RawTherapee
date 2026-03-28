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

#include "gtk4.h"

#include "rtengine/util/optional.h"

#include <fmt/format.h>
#include <gtkmm/widget.h>

#ifdef GDK_WINDOWING_QUARTZ
    #include <gdk/gdkquartz.h>
#endif
#ifdef GDK_WINDOWING_WAYLAND
    #include <gdk/gdkwayland.h>
#endif
#ifdef GDK_WINDOWING_WIN32
    #include <gdk/gdkwin32.h>
#endif
#ifdef GDK_WINDOWING_X11
    #include <gdk/gdkx.h>
#endif

#include <cmath>
#include <stdexcept>

namespace {

bool assumeSurfaceScrollDevice(const GdkEvent* event)
{
    GdkDevice* device = gdk_event_get_source_device(event);
    GdkInputSource input_source = gdk_device_get_source(device);

    switch (input_source) {
        case GDK_SOURCE_PEN:
        case GDK_SOURCE_TOUCHSCREEN:
        case GDK_SOURCE_TOUCHPAD:
        case GDK_SOURCE_TRACKPOINT:
            return true;
        default:
            return false;
    }
}

}  // namespace

namespace rt {
namespace gtk4 {

Glib::RefPtr<GestureClick> GestureClick::create(Gtk::Widget* widget)
{
    return Glib::RefPtr<GestureClick>(new GestureClick(widget));
}

GestureClick::GestureClick(Gtk::Widget* widget)
    : Glib::ObjectBase("RtGestureClick"),
      m_controller(gtk_gesture_multi_press_new(GTK_WIDGET(widget->gobj())))
{
    g_signal_connect(GTK_GESTURE_MULTI_PRESS(m_controller), "pressed",
                     G_CALLBACK(&GestureClick::hook_pressed), this);
    g_signal_connect(GTK_GESTURE_MULTI_PRESS(m_controller), "released",
                     G_CALLBACK(&GestureClick::hook_released), this);
    g_signal_connect(GTK_GESTURE_MULTI_PRESS(m_controller), "stopped",
                     G_CALLBACK(&GestureClick::hook_stopped), this);
}

GestureClick::~GestureClick()
{
    g_object_unref(m_controller);
    m_controller = nullptr;
}

GtkGestureMultiPress* GestureClick::gobj()
{
    return GTK_GESTURE_MULTI_PRESS(m_controller);
}

GtkPropagationPhase GestureClick::get_propagation_phase() const
{
    return gtk_event_controller_get_propagation_phase(
        GTK_EVENT_CONTROLLER(m_controller));
}

void GestureClick::set_propagation_phase(GtkPropagationPhase phase)
{
    gtk_event_controller_set_propagation_phase(
        GTK_EVENT_CONTROLLER(m_controller), phase);
}

unsigned int GestureClick::get_current_button() const
{
    return gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(m_controller));
}

unsigned int GestureClick::get_button() const
{
    return gtk_gesture_single_get_button(GTK_GESTURE_SINGLE(m_controller));
}

void GestureClick::set_button(unsigned int button)
{
    return gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(m_controller), button);
}

void GestureClick::hook_pressed(GtkGestureMultiPress*, gint n_press, gdouble x,
                                gdouble y, gpointer user_data)
{
    static_cast<GestureClick*>(user_data)->m_pressed.emit(n_press, x, y);
}

void GestureClick::hook_released(GtkGestureMultiPress*, gint n_press, gdouble x,
                                 gdouble y, gpointer user_data)
{
    static_cast<GestureClick*>(user_data)->m_released.emit(n_press, x, y);
}

void GestureClick::hook_stopped(GtkGestureMultiPress*, gpointer user_data)
{
    static_cast<GestureClick*>(user_data)->m_stopped.emit();
}

Glib::RefPtr<EventControllerKey> EventControllerKey::create(Gtk::Widget* widget)
{
    return Glib::RefPtr<EventControllerKey>(new EventControllerKey(widget));
}

EventControllerKey::EventControllerKey(Gtk::Widget* widget)
    : Glib::ObjectBase("RtEventControllerKey"),
      m_controller(gtk_event_controller_key_new(GTK_WIDGET(widget->gobj())))
{
    g_signal_connect(GTK_EVENT_CONTROLLER_KEY(m_controller), "key-pressed",
                     G_CALLBACK(&EventControllerKey::hook_pressed), this);
    g_signal_connect(GTK_EVENT_CONTROLLER_KEY(m_controller), "key-released",
                     G_CALLBACK(&EventControllerKey::hook_released), this);
    g_signal_connect(GTK_EVENT_CONTROLLER_KEY(m_controller), "modifiers",
                     G_CALLBACK(&EventControllerKey::hook_modifiers), this);
}

EventControllerKey::~EventControllerKey()
{
    g_object_unref(m_controller);
    m_controller = nullptr;
}

GtkEventControllerKey* EventControllerKey::gobj()
{
    return GTK_EVENT_CONTROLLER_KEY(m_controller);
}

GtkPropagationPhase EventControllerKey::get_propagation_phase() const
{
    return gtk_event_controller_get_propagation_phase(
        GTK_EVENT_CONTROLLER(m_controller));
}

void EventControllerKey::set_propagation_phase(GtkPropagationPhase phase)
{
    gtk_event_controller_set_propagation_phase(
        GTK_EVENT_CONTROLLER(m_controller), phase);
}

bool EventControllerKey::hook_pressed(
    GtkEventControllerKey*, guint keyval, guint keycode, GdkModifierType state,
    gpointer user_data)
{
    return static_cast<EventControllerKey*>(user_data)
        ->m_pressed.emit(keyval, keycode, state);
}

void EventControllerKey::hook_released(
    GtkEventControllerKey*, guint keyval, guint keycode, GdkModifierType state,
    gpointer user_data)
{
    static_cast<EventControllerKey*>(user_data)->m_released.emit(keyval, keycode, state);
}

bool EventControllerKey::hook_modifiers(GtkEventControllerKey*, GdkModifierType state,
                                        gpointer user_data)
{
    return static_cast<EventControllerKey*>(user_data)->m_modifiers.emit(state);
}

Glib::RefPtr<EventControllerMotion> EventControllerMotion::create(Gtk::Widget* widget)
{
    return Glib::RefPtr<EventControllerMotion>(new EventControllerMotion(widget));
}

EventControllerMotion::EventControllerMotion(Gtk::Widget* widget)
    : Glib::ObjectBase("RtEventControllerMotion"),
      m_controller(gtk_event_controller_motion_new(GTK_WIDGET(widget->gobj())))
{
    g_signal_connect(GTK_EVENT_CONTROLLER_MOTION(m_controller), "enter",
                     G_CALLBACK(&EventControllerMotion::hook_enter), this);
    g_signal_connect(GTK_EVENT_CONTROLLER_MOTION(m_controller), "motion",
                     G_CALLBACK(&EventControllerMotion::hook_motion), this);
    g_signal_connect(GTK_EVENT_CONTROLLER_MOTION(m_controller), "leave",
                     G_CALLBACK(&EventControllerMotion::hook_leave), this);
}

EventControllerMotion::~EventControllerMotion()
{
    g_object_unref(m_controller);
    m_controller = nullptr;
}

GtkEventControllerMotion* EventControllerMotion::gobj()
{
    return GTK_EVENT_CONTROLLER_MOTION(m_controller);
}

GtkPropagationPhase EventControllerMotion::get_propagation_phase() const
{
    return gtk_event_controller_get_propagation_phase(
        GTK_EVENT_CONTROLLER(m_controller));
}

void EventControllerMotion::set_propagation_phase(GtkPropagationPhase phase)
{
    gtk_event_controller_set_propagation_phase(
        GTK_EVENT_CONTROLLER(m_controller), phase);
}

void EventControllerMotion::hook_enter(GtkEventControllerMotion*, gdouble x, gdouble y,
                                       gpointer user_data)
{
    static_cast<EventControllerMotion*>(user_data)->m_enter.emit(x, y);
}

void EventControllerMotion::hook_motion(GtkEventControllerMotion*, gdouble x, gdouble y,
                                        gpointer user_data)
{
    static_cast<EventControllerMotion*>(user_data)->m_motion.emit(x, y);
}

void EventControllerMotion::hook_leave(GtkEventControllerMotion*, gpointer user_data)
{
    static_cast<EventControllerMotion*>(user_data)->m_leave.emit();
}

HeuristicEventControllerScroll::HeuristicEventControllerScroll(Gtk::Widget* widget,
                                                               Flags flags)
    : Glib::ObjectBase("RtHeuristicEventControllerScroll"),
      m_widget(widget),
      m_flags(flags),
      m_scroll_unit(ScrollUnit::WHEEL),
      m_dx_accum(0),
      m_dy_accum(0),
      m_is_active(false)
{
}

void HeuristicEventControllerScroll::beginSmoothScrollIfNeeded()
{
    if (m_is_active) return;

    m_signal_begin.emit();
    m_is_active = true;
}

void HeuristicEventControllerScroll::endScroll()
{
    if (!m_is_active) return;

    m_signal_end.emit();
    m_is_active = false;
    m_dx_accum = 0;
    m_dy_accum = 0;
}

void HeuristicEventControllerScroll::populateDeltasByDirection(
    GdkScrollDirection dir, double& dx, double& dy)
{
    switch (dir) {
        case GDK_SCROLL_UP:
            dx = 0;
            dy = -1;
            break;
        case GDK_SCROLL_DOWN:
            dx = 0;
            dy = 1;
            break;
        case GDK_SCROLL_LEFT:
            dx = -1;
            dy = 0;
            break;
        case GDK_SCROLL_RIGHT:
            dx = 1;
            dy = 0;
            break;
        default:
            dx = 0;
            dy = 0;
            break;
    }

    filterDeltas(dx, dy);
}

void HeuristicEventControllerScroll::filterDeltas(double& dx, double& dy)
{
    if (rt::none(m_flags & Flags::VERTICAL)) dy = 0;
    if (rt::none(m_flags & Flags::HORIZONTAL)) dx = 0;
}

void HeuristicEventControllerScroll::makeDiscrete(double& dx, double& dy)
{
    m_dx_accum += dx;
    m_dy_accum += dy;
    dx = 0;
    dy = 0;

    if (m_scroll_unit == rt::gtk4::ScrollUnit::SURFACE) {
        // From GTK 4's implementation of EventControllerScroll
        constexpr double SURFACE_UNIT_DISCRETE_MAPPING = 10;

        dx = std::floor(m_dx_accum / SURFACE_UNIT_DISCRETE_MAPPING);
        dy = std::floor(m_dy_accum / SURFACE_UNIT_DISCRETE_MAPPING);

        m_dx_accum -= dx * SURFACE_UNIT_DISCRETE_MAPPING;
        m_dy_accum -= dy * SURFACE_UNIT_DISCRETE_MAPPING;

        m_scroll_unit = rt::gtk4::ScrollUnit::WHEEL;
    } else {
        if (std::abs(m_dx_accum) >= 1) {
            double steps = std::floor(m_dx_accum);
            m_dx_accum -= steps;
            dx = steps;
        }
        if (std::abs(m_dy_accum) >= 1) {
            double steps = std::floor(m_dy_accum);
            m_dy_accum -= steps;
            dy = steps;
        }
    }
}

bool HeuristicEventControllerScroll::onEvent(GdkEvent* event)
{
    if (gdk_event_get_event_type(event) != GDK_SCROLL) return false;
    if (rt::none(m_flags & Flags::BOTH_AXES)) return false;

    GdkDisplay* display = gdk_display_get_default();

    // clang-format off
#ifdef GDK_WINDOWING_QUARTZ
    if (GDK_IS_QUARTZ_DISPLAY(display)) {
        return handleEventMacOS(event);
    }
    else
#endif
#ifdef GDK_WINDOWING_WAYLAND
    if (GDK_IS_WAYLAND_DISPLAY(display)) {
        return handleEventWayland(event);
    }
    else
#endif
#ifdef GDK_WINDOWING_WIN32
    if (GDK_IS_WIN32_DISPLAY(display)) {
        return handleEventWin32(event);
    }
    else
#endif
#ifdef GDK_WINDOWING_X11
    if (GDK_IS_X11_DISPLAY(display)) {
        return handleEventX11(event);
    }
    else
#endif
    {
        return handleEventFallback(event);
    }
    // clang-format on
}

// Reference: https://gitlab.gnome.org/GNOME/gtk/-/blob/gtk-3-24/gdk/quartz/gdkevents-quartz.c
//
// Implementation based on gtk-3-24 gdkevents-quartz.c fill_scroll_event()
bool HeuristicEventControllerScroll::handleEventMacOS(const GdkEvent* event)
{
    auto scroll_event = reinterpret_cast<const GdkEventScroll*>(event);

    double dx = 0;
    double dy = 0;

    if (scroll_event->direction == GDK_SCROLL_SMOOTH) {
        m_scroll_unit = ScrollUnit::SURFACE;

        beginSmoothScrollIfNeeded();

        dx = scroll_event->delta_x;
        dy = scroll_event->delta_y;
        filterDeltas(dx, dy);
    } else {
        m_scroll_unit = ScrollUnit::WHEEL;

        populateDeltasByDirection(scroll_event->direction, dx, dy);
        // In MacOS, the deltas are absolute values for discrete mouse wheel
        // detents. These values represent the "speed" of a scroll. If a user
        // scrolls quickly, MacOS interprets that as the user wanting to scroll
        // more in a single detent.
        //
        // TODO: Option to disable acceleration
        dx *= scroll_event->delta_x;
        dy *= scroll_event->delta_y;
    }

    bool handled = false;

    if (dx != 0 || dy != 0) {
        handled = m_signal_scroll.emit(dx, dy);
    }

    if (m_is_active && scroll_event->is_stop) {
        endScroll();
        // Allow scroll stop to propagate up widget hierarchy
        handled = false;
    }

    return handled;
}

bool HeuristicEventControllerScroll::handleEventWayland(const GdkEvent* event)
{
    // libinput's default mouse detent is 15 "units" which translates to a
    // delta of 1.5 or -1.5 in GTK 3. This value does not change with display
    // scaling (e.g. fractional scaling of 125% or device scale of 2).
    constexpr double LIBINPUT_DETENT = 1.5;

    return handleEventWithInferredSingleMouseDetents(event, LIBINPUT_DETENT);
}

bool HeuristicEventControllerScroll::handleEventWin32(const GdkEvent* event)
{
    return false;
}

bool HeuristicEventControllerScroll::handleEventX11(const GdkEvent* event)
{
    // XInput2's mouse detent is 10 "units" which translates to a delta of 1
    // or -1 in GTK 3. This value does not change with display scaling (e.g.
    // fractional scaling of 125% or device scale of 2).
    constexpr double XINPUT2_DETENT = 1;

    return handleEventWithInferredSingleMouseDetents(event, XINPUT2_DETENT);
}

bool HeuristicEventControllerScroll::handleEventWithInferredSingleMouseDetents(
    const GdkEvent* event, double detent_delta)
{
    auto scroll_event = reinterpret_cast<const GdkEventScroll*>(event);

    double dx = scroll_event->delta_x;
    double dy = scroll_event->delta_y;

    // In GTK 4, the scroll unit is provided by the GTK event controller.
    // This information is not available in GTK 3. However, we can guess based
    // on heuristics if a smooth scroll event represents a mouse detent.
    auto inferred_mouse_detent = [&]() -> rt::optional<GdkScrollDirection>
    {
        if (scroll_event->direction != GDK_SCROLL_SMOOTH) return rt::nullopt;
        // If we are already in a smooth scroll, don't suddenly switch to mouse
        // detent scrolling.
        if (m_is_active) return rt::nullopt;
        if (scroll_event->is_stop) return rt::nullopt;

        if (assumeSurfaceScrollDevice(event)) return rt::nullopt;

        // A mouse detent cannot be in multiple directions at the same time
        if ((dx != 0) && (dy != 0)) return rt::nullopt;

        if ((dx == 0) && rt::any(m_flags & Flags::VERTICAL)) {
            if (dy == detent_delta) {
                return GDK_SCROLL_DOWN;
            } else if (dy == -detent_delta) {
                return GDK_SCROLL_UP;
            }
        } else if ((dy == 0) && rt::any(m_flags & Flags::HORIZONTAL)) {
            if (dx == detent_delta) {
                return GDK_SCROLL_RIGHT;
            } else if (dx == -detent_delta) {
                return GDK_SCROLL_LEFT;
            }
        }

        return rt::nullopt;
    }();

    if (inferred_mouse_detent) {
        m_scroll_unit = ScrollUnit::WHEEL;

        // Emit a begin signal but don't set m_is_active (i.e. don't start a
        // smooth scroll)
        m_signal_begin.emit();

        populateDeltasByDirection(*inferred_mouse_detent, dx, dy);
    } else if (scroll_event->direction == GDK_SCROLL_SMOOTH) {
        m_scroll_unit = ScrollUnit::SURFACE;

        beginSmoothScrollIfNeeded();

        filterDeltas(dx, dy);
    } else {
        m_scroll_unit = ScrollUnit::WHEEL;

        populateDeltasByDirection(scroll_event->direction, dx, dy);
    }

    bool handled = false;

    if (dx != 0 || dy != 0) {
        handled = m_signal_scroll.emit(dx, dy);
    }

    if (m_is_active && scroll_event->is_stop) {
        endScroll();
        // Allow scroll stop to propagate up widget hierarchy
        handled = false;
    }

    return handled;
}

// Force scroll events to be discrete
bool HeuristicEventControllerScroll::handleEventFallback(GdkEvent* event)
{
    auto scroll_event = reinterpret_cast<GdkEventScroll*>(event);

    m_scroll_unit = ScrollUnit::WHEEL;

    double dx = 0;
    double dy = 0;

    if (scroll_event->direction == GDK_SCROLL_SMOOTH) {
        beginSmoothScrollIfNeeded();

        dx = scroll_event->delta_x;
        dy = scroll_event->delta_y;
        filterDeltas(dx, dy);
        makeDiscrete(dx, dy);
    } else {
        populateDeltasByDirection(scroll_event->direction, dx, dy);
    }

    bool handled = false;

    if (dx != 0 || dy != 0) {
        handled = m_signal_scroll.emit(dx, dy);
        if (!handled) {
            // The discretized value was not handled so propagate the coalesced
            // scroll event deltas.
            scroll_event->delta_x = dx;
            scroll_event->delta_y = dy;
        }
    } else if (scroll_event->direction == GDK_SCROLL_SMOOTH) {
        // Capture all coalesced smooth scroll events
        handled = m_is_active;
    }

    if (m_is_active && scroll_event->is_stop) {
        endScroll();
        // Allow scroll stop to propagate up widget hierarchy
        handled = false;
    }

    return handled;
}

}  // namespace gtk4
}  // namespace rt
