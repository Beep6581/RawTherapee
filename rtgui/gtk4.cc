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

#include <gtkmm/widget.h>

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

Glib::RefPtr<EventControllerScroll> EventControllerScroll::create(
    Gtk::Widget* widget, GtkEventControllerScrollFlags flags)
{
    return Glib::RefPtr<EventControllerScroll>(new EventControllerScroll(widget, flags));
}

EventControllerScroll::EventControllerScroll(Gtk::Widget* widget,
                                             GtkEventControllerScrollFlags flags)
    : Glib::ObjectBase("RtEventControllerScroll"),
      m_controller(gtk_event_controller_scroll_new(GTK_WIDGET(widget->gobj()), flags))
{
    g_signal_connect(GTK_EVENT_CONTROLLER_SCROLL(m_controller), "scroll-begin",
                     G_CALLBACK(&EventControllerScroll::hook_begin), this);
    g_signal_connect(GTK_EVENT_CONTROLLER_SCROLL(m_controller), "scroll",
                     G_CALLBACK(&EventControllerScroll::hook_scroll), this);
    g_signal_connect(GTK_EVENT_CONTROLLER_SCROLL(m_controller), "scroll-end",
                     G_CALLBACK(&EventControllerScroll::hook_end), this);
}

EventControllerScroll::~EventControllerScroll()
{
    g_object_unref(m_controller);
    m_controller = nullptr;
}

GtkEventControllerScroll* EventControllerScroll::gobj()
{
    return GTK_EVENT_CONTROLLER_SCROLL(m_controller);
}

GtkPropagationPhase EventControllerScroll::get_propagation_phase() const
{
    return gtk_event_controller_get_propagation_phase(
        GTK_EVENT_CONTROLLER(m_controller));
}

void EventControllerScroll::set_propagation_phase(GtkPropagationPhase phase)
{
    gtk_event_controller_set_propagation_phase(
        GTK_EVENT_CONTROLLER(m_controller), phase);
}

GtkEventControllerScrollFlags EventControllerScroll::get_flags() const
{
    return gtk_event_controller_scroll_get_flags(
        GTK_EVENT_CONTROLLER_SCROLL(m_controller));
}

void EventControllerScroll::set_flags(GtkEventControllerScrollFlags flags)
{
    return gtk_event_controller_scroll_set_flags(
        GTK_EVENT_CONTROLLER_SCROLL(m_controller), flags);
}

void EventControllerScroll::hook_scroll(GtkEventControllerScroll*, gdouble dx,
                                        gdouble dy, gpointer user_data)
{
    static_cast<EventControllerScroll*>(user_data)->m_scroll.emit(dx, dy);
}

void EventControllerScroll::hook_begin(GtkEventControllerScroll*, gpointer user_data)
{
    static_cast<EventControllerScroll*>(user_data)->m_begin.emit();
}

void EventControllerScroll::hook_end(GtkEventControllerScroll*, gpointer user_data)
{
    static_cast<EventControllerScroll*>(user_data)->m_end.emit();
}

}  // namespace gtk4
}  // namespace rt
