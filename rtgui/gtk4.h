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

#pragma once

#include "rtengine/util/enum.h"

#include <glibmm/object.h>
#include <glibmm/refptr.h>
#include <glibmm/signalproxy.h>
#include <gtk/gtk.h>
#include <sigc++/sigc++.h>

namespace Gtk {
class Widget;
}

namespace rt {
namespace gtk4 {

// TODO(GTK4): Replace with proper Gtk::Gesture*/EventController* gtkmm wrappers.
//
// WARNING: You must store the RefPtr of the event controllers as adding the
//          controller to the widget does not increase the ref count of the
//          controller in GTK 3.

class GestureClick : public Glib::Object
{
public:
    static Glib::RefPtr<GestureClick> create(Gtk::Widget* widget);

    ~GestureClick();

    GtkGestureMultiPress* gobj();

    GtkPropagationPhase get_propagation_phase() const;
    void set_propagation_phase(GtkPropagationPhase phase);

    unsigned int get_current_button() const;
    unsigned int get_button() const;
    void set_button(unsigned int button = 0);

    sigc::signal<void(int, double, double)>& signal_pressed() { return m_pressed; }
    sigc::signal<void(int, double, double)>& signal_released() { return m_released; }
    sigc::signal<void()>& signal_stopped() { return m_stopped; }
private:
    static void hook_pressed(GtkGestureMultiPress*, gint n_press, gdouble x, gdouble y,
                             gpointer user_data);
    static void hook_released(GtkGestureMultiPress*, gint n_press, gdouble x, gdouble y,
                              gpointer user_data);
    static void hook_stopped(GtkGestureMultiPress*, gpointer user_data);

    GestureClick(Gtk::Widget* widget);

    sigc::signal<void(int, double, double)> m_pressed;
    sigc::signal<void(int, double, double)> m_released;
    sigc::signal<void()> m_stopped;
    GtkGesture* m_controller;
};

class EventControllerKey : public Glib::Object
{
public:
    static Glib::RefPtr<EventControllerKey> create(Gtk::Widget* widget);

    ~EventControllerKey();

    GtkEventControllerKey* gobj();

    GtkPropagationPhase get_propagation_phase() const;
    void set_propagation_phase(GtkPropagationPhase phase);

    sigc::signal<bool(guint, guint, GdkModifierType)>&
    signal_key_pressed() { return m_pressed; }
    sigc::signal<void(guint, guint, GdkModifierType)>&
    signal_key_released() { return m_released; }
    sigc::signal<bool(GdkModifierType)>& signal_modifiers() { return m_modifiers; }

private:
    static bool hook_pressed(GtkEventControllerKey*, guint keyval, guint keycode,
                             GdkModifierType state, gpointer user_data);
    static void hook_released(GtkEventControllerKey*, guint keyval, guint keycode,
                              GdkModifierType state, gpointer user_data);
    static bool hook_modifiers(GtkEventControllerKey*, GdkModifierType state,
                               gpointer user_data);

    EventControllerKey(Gtk::Widget* widget);

    sigc::signal<bool(guint, guint, GdkModifierType)> m_pressed;
    sigc::signal<void(guint, guint, GdkModifierType)> m_released;
    sigc::signal<bool(GdkModifierType)> m_modifiers;
    GtkEventController* m_controller;
};

class EventControllerMotion : public Glib::Object
{
public:
    static Glib::RefPtr<EventControllerMotion> create(Gtk::Widget* widget);

    ~EventControllerMotion();

    GtkEventControllerMotion* gobj();

    GtkPropagationPhase get_propagation_phase() const;
    void set_propagation_phase(GtkPropagationPhase phase);

    sigc::signal<void(double, double)>& signal_enter() { return m_enter; }
    sigc::signal<void(double, double)>& signal_motion() { return m_motion; }
    sigc::signal<void()>& signal_leave() { return m_leave; }

private:
    static void hook_enter(GtkEventControllerMotion*, gdouble x, gdouble y,
                           gpointer user_data);
    static void hook_motion(GtkEventControllerMotion*, gdouble x, gdouble y,
                            gpointer user_data);
    static void hook_leave(GtkEventControllerMotion*, gpointer user_data);

    EventControllerMotion(Gtk::Widget* widget);

    sigc::signal<void(double, double)> m_enter;
    sigc::signal<void(double, double)> m_motion;
    sigc::signal<void()> m_leave;
    GtkEventController* m_controller;
};

enum class ScrollUnit { WHEEL, SURFACE };

// When porting to GTK 4, replace uses with Gtk::EventControllerScroll
class HeuristicEventControllerScroll : public Glib::Object
{
public:
    enum class Flags {
        NONE = 0,
        VERTICAL = (1 << 0),
        HORIZONTAL = (1 << 1),
        BOTH_AXES = VERTICAL | HORIZONTAL,
    };

    HeuristicEventControllerScroll(Gtk::Widget* widget, Flags flags);

    // Widgets using this backport should call this event handler in their
    // generic event handling.
    bool onEvent(GdkEvent* event);

    ScrollUnit get_scroll_unit() const { return m_scroll_unit; }
    Flags get_flags() const { return m_flags; }
    void set_flags(Flags flags) { m_flags = flags; }

    sigc::signal<void()>& signal_scroll_begin() { return m_signal_begin; }
    sigc::signal<bool(double, double)>& signal_scroll() { return m_signal_scroll; }
    sigc::signal<void()>& signal_scroll_end() { return m_signal_end; }

private:
    void beginSmoothScrollIfNeeded();
    void endScroll();

    void populateDeltasByDirection(GdkScrollDirection dir, double& dx, double& dy);
    void filterDeltas(double& dx, double& dy);
    void makeDiscrete(double& dx, double& dy);

    bool handleEventMacOS(const GdkEvent* event);
    bool handleEventWayland(const GdkEvent* event);
    bool handleEventWin32(const GdkEvent* event);
    bool handleEventX11(const GdkEvent* event);

    bool handleEventWithInferredSingleMouseDetents(const GdkEvent* event,
                                                   double detent_delta);
    bool handleEventFallback(GdkEvent* event);

    sigc::signal<bool(double, double)> m_signal_scroll;
    sigc::signal<void()> m_signal_begin;
    sigc::signal<void()> m_signal_end;

    Gtk::Widget* m_widget;
    Flags m_flags;
    ScrollUnit m_scroll_unit;
    double m_dx_accum;
    double m_dy_accum;
    bool m_is_active;
};

}  // namespace gtk4
}  // namespace rt

template <>
struct rt::EnumAsBitflags<rt::gtk4::HeuristicEventControllerScroll::Flags>
    : std::true_type {};
