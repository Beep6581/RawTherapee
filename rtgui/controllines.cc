/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2020 Lawrence Lee <billee@ucdavis.edu>
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
#include <memory>

#include "controllines.h"
#include "editcallbacks.h"
#include "editwidgets.h"
#include "rtsurface.h"

#include "../rtengine/perspectivecorrection.h"

using namespace rtengine;

::ControlLine::~ControlLine() = default;

ControlLineManager::ControlLineManager():
    EditSubscriber(ET_OBJECTS),
    canvas_area(new Rectangle()),
    cursor(CSHandOpen),
    draw_mode(false),
    drawing_line(false),
    edited(false),
    prev_obj(-1),
    selected_object(-1)
{
    canvas_area->filled = true;
    canvas_area->topLeft = Coord(0, 0);
    mouseOverGeometry.push_back(canvas_area.get());

    line_icon_h = Cairo::RefPtr<RTSurface>(new RTSurface(
            "bidirectional-arrow-horizontal-hicontrast.png"));
    line_icon_v = Cairo::RefPtr<RTSurface>(new RTSurface(
            "bidirectional-arrow-vertical-hicontrast.png"));
    line_icon_h_prelight = Cairo::RefPtr<RTSurface>(new RTSurface(
                               "bidirectional-arrow-horizontal-prelight.png"));
    line_icon_v_prelight = Cairo::RefPtr<RTSurface>(new RTSurface(
                               "bidirectional-arrow-vertical-prelight.png"));
}

ControlLineManager::~ControlLineManager() = default;

void ControlLineManager::setActive(bool active)
{
    EditDataProvider* provider = getEditProvider();

    if (!provider || (this == provider->getCurrSubscriber()) == active) {
        return;
    }

    if (active) {
        subscribe();

        int ih, iw;
        provider->getImageSize(iw, ih);
        canvas_area->bottomRight = Coord(iw, ih);
    } else {
        unsubscribe();
    }
}

void ControlLineManager::setDrawMode(bool draw)
{
    draw_mode = draw;
}

size_t ControlLineManager::size(void) const
{
    return control_lines.size();
}

bool ControlLineManager::button1Pressed(int modifierKey)
{
    EditDataProvider* dataProvider = getEditProvider();

    if (!dataProvider) {
        return false;
    }

    drag_delta = Coord(0, 0);

    const int object = dataProvider->getObject();

    if (object > 0) { // A control line.
        if (object % ::ControlLine::OBJ_COUNT == 2) { // Icon.
            action = Action::PICKING;
        } else {
            selected_object = object;
            action = Action::DRAGGING;
        }
    } else if (draw_mode && (modifierKey & GDK_CONTROL_MASK)) { // Add new line.
        addLine(dataProvider->posImage, dataProvider->posImage);
        drawing_line = true;
        selected_object = mouseOverGeometry.size() - 1; // Select endpoint.
        action = Action::DRAGGING;
    }

    return true;
}

bool ControlLineManager::button1Released(void)
{
    action = Action::NONE;

    if (selected_object > 0) {
        mouseOverGeometry[selected_object]->state = Geometry::NORMAL;
    }

    edited = true;
    callbacks->lineChanged();
    drawing_line = false;
    selected_object = -1;
    return false;
}

bool ControlLineManager::button3Pressed(int modifierKey)
{
    EditDataProvider* provider = getEditProvider();

    action = Action::NONE;

    if (!provider || provider->getObject() < 1) {
        return false;
    }

    action = Action::PICKING;
    return false;
}

bool ControlLineManager::pick1(bool picked)
{
    action = Action::NONE;

    if (!picked) {
        return false;
    }

    EditDataProvider* provider = getEditProvider();

    if (!provider || provider->getObject() % ::ControlLine::OBJ_COUNT != 2) {
        return false;
    }

    // Change line type.
    int object_id = provider->getObject();
    ::ControlLine& line =
        *control_lines[(object_id - 1) / ::ControlLine::OBJ_COUNT];

    if (line.type == rtengine::ControlLine::HORIZONTAL) {
        line.icon = line.icon_v;
        line.type = rtengine::ControlLine::VERTICAL;
    } else if (line.type == rtengine::ControlLine::VERTICAL) {
        line.icon = line.icon_h;
        line.type = rtengine::ControlLine::HORIZONTAL;
    }

    visibleGeometry[object_id - 1] = line.icon.get();

    edited = true;
    callbacks->lineChanged();

    return true;
}

bool ControlLineManager::pick3(bool picked)
{
    action = Action::NONE;

    if (!picked) {
        return false;
    }

    EditDataProvider* provider = getEditProvider();

    if (!provider) {
        return false;
    }

    removeLine((provider->getObject() - 1) / ::ControlLine::OBJ_COUNT);
    prev_obj = -1;
    selected_object = -1;
    return false;
}

bool ControlLineManager::drag1(int modifierKey)
{
    EditDataProvider* provider = getEditProvider();

    if (!provider || selected_object < 1) {
        return false;
    }

    ::ControlLine& control_line =
        *control_lines[(selected_object - 1) / ::ControlLine::OBJ_COUNT];
    // 0 == end, 1 == line, 2 == icon, 3 == begin
    int component = selected_object % ::ControlLine::OBJ_COUNT;
    Coord mouse = provider->posImage + provider->deltaImage;
    Coord delta = provider->deltaImage - drag_delta;
    int ih, iw;
    provider->getImageSize(iw, ih);

    switch (component) {
        case (0): // end
            control_line.end->center = mouse;
            control_line.end->center.clip(iw, ih);
            control_line.line->end = control_line.end->center;
            control_line.end->state = Geometry::DRAGGED;
            break;

        case (1): { // line
            // Constrain delta so the end stays above the image.
            Coord new_delta = control_line.end->center + delta;
            new_delta.clip(iw, ih);
            new_delta -= control_line.end->center;
            // Constrain delta so the beginning stays above the image.
            new_delta += control_line.begin->center;
            new_delta.clip(iw, ih);
            new_delta -= control_line.begin->center;
            // Move all objects in the control line.
            control_line.end->center += new_delta;
            control_line.begin->center += new_delta;
            control_line.line->end = control_line.end->center;
            control_line.line->begin = control_line.begin->center;
            drag_delta += new_delta;
            control_line.line->state = Geometry::DRAGGED;
            break;
        }

        case (3): // begin
            control_line.begin->center = mouse;
            control_line.begin->center.clip(iw, ih);
            control_line.line->begin = control_line.begin->center;
            control_line.begin->state = Geometry::DRAGGED;
            break;
    }

    control_line.icon_h->position.x = (control_line.begin->center.x +
                                       control_line.end->center.x) / 2;
    control_line.icon_h->position.y = (control_line.begin->center.y +
                                       control_line.end->center.y) / 2;
    control_line.icon_v->position.x = control_line.icon_h->position.x;
    control_line.icon_v->position.y = control_line.icon_h->position.y;

    if (drawing_line) {
        autoSetLineType(selected_object);
    }

    return false;
}

bool ControlLineManager::getEdited(void) const
{
    return edited;
}

CursorShape ControlLineManager::getCursor(int objectID) const
{
    return cursor;
}

bool ControlLineManager::mouseOver(int modifierKey)
{
    EditDataProvider* provider = getEditProvider();

    if (!provider) {
        return false;
    }

    int cur_obj = provider->getObject();

    if (cur_obj == 0) { // Canvas
        if (draw_mode && modifierKey & GDK_CONTROL_MASK) {
            cursor = CSCrosshair;
        } else {
            cursor = CSHandOpen;
        }
    } else if (cur_obj < 0) { // Nothing
        cursor = CSArrow;
    } else if (cur_obj % ::ControlLine::OBJ_COUNT == 2) { // Icon
        visibleGeometry[cur_obj - 1]->state = Geometry::PRELIGHT;
        cursor = CSArrow;
    } else { // Object
        visibleGeometry[cur_obj - 1]->state = Geometry::PRELIGHT;
        cursor = CSMove2D;
    }

    if (prev_obj != cur_obj && prev_obj > 0) {
        visibleGeometry[prev_obj - 1]->state = Geometry::NORMAL;
    }

    prev_obj = cur_obj;

    return true;
}

void ControlLineManager::switchOffEditMode(void)
{
    if (callbacks) {
        callbacks->switchOffEditMode();
    }
}

void ControlLineManager::setEdited(bool edited)
{
    this->edited = edited;
}

void ControlLineManager::setEditProvider(EditDataProvider* provider)
{
    EditSubscriber::setEditProvider(provider);
}

void ControlLineManager::setLines(const std::vector<rtengine::ControlLine>&
                                  lines)
{
    removeAll();

    for (auto&& line : lines) {
        Coord start(line.x1, line.y1);
        Coord end(line.x2, line.y2);
        addLine(start, end, line.type);
    }
}

void ControlLineManager::addLine(Coord begin, Coord end,
                                 rtengine::ControlLine::Type type)
{
    constexpr int line_width = 2;
    constexpr int handle_radius = 6;
    std::unique_ptr<Line> line;
    std::shared_ptr<OPIcon> icon_h, icon_v;
    std::unique_ptr<Circle> begin_c, end_c;

    line = std::unique_ptr<Line>(new Line());
    line->datum = Geometry::IMAGE;
    line->innerLineWidth = line_width;
    line->begin = begin;
    line->end = end;

    const Cairo::RefPtr<RTSurface> null_surface =
        Cairo::RefPtr<RTSurface>(nullptr);

    icon_h = std::make_shared<OPIcon>(line_icon_h, null_surface,
                                      line_icon_h_prelight,
                                      null_surface, null_surface,
                                      Geometry::DP_CENTERCENTER);
    icon_h->position = Coord((begin.x + end.x) / 2, (begin.y + end.y) / 2);

    icon_v = std::make_shared<OPIcon>(line_icon_v, null_surface,
                                      line_icon_v_prelight,
                                      null_surface, null_surface,
                                      Geometry::DP_CENTERCENTER);
    icon_v->position = Coord((begin.x + end.x) / 2, (begin.y + end.y) / 2);

    begin_c = std::unique_ptr<Circle>(new Circle());
    begin_c->datum = Geometry::IMAGE;
    begin_c->filled = true;
    begin_c->radius = handle_radius;
    begin_c->center = begin;

    end_c = std::unique_ptr<Circle>(new Circle());
    end_c->datum = Geometry::IMAGE;
    end_c->filled = true;
    end_c->radius = handle_radius;
    end_c->center = end;

    std::unique_ptr<::ControlLine> control_line(new ::ControlLine());
    control_line->begin = std::move(begin_c);
    control_line->end = std::move(end_c);
    control_line->icon_h = icon_h;
    control_line->icon_v = icon_v;

    if (type == rtengine::ControlLine::HORIZONTAL) {
        control_line->icon = icon_h;
    } else {
        control_line->icon = icon_v;
    }

    control_line->line = std::move(line);
    control_line->type = type;

    EditSubscriber::visibleGeometry.push_back(control_line->line.get());
    EditSubscriber::visibleGeometry.push_back(control_line->icon.get());
    EditSubscriber::visibleGeometry.push_back(control_line->begin.get());
    EditSubscriber::visibleGeometry.push_back(control_line->end.get());

    EditSubscriber::mouseOverGeometry.push_back(control_line->line.get());
    EditSubscriber::mouseOverGeometry.push_back(control_line->icon.get());
    EditSubscriber::mouseOverGeometry.push_back(control_line->begin.get());
    EditSubscriber::mouseOverGeometry.push_back(control_line->end.get());

    control_lines.push_back(std::move(control_line));
}

void ControlLineManager::autoSetLineType(int object_id)
{
    int line_id = (object_id - 1) / ::ControlLine::OBJ_COUNT;
    ::ControlLine& line = *control_lines[line_id];

    int dx = line.begin->center.x - line.end->center.x;
    int dy = line.begin->center.y - line.end->center.y;

    if (dx < 0) {
        dx = -dx;
    }

    if (dy < 0) {
        dy = -dy;
    }

    rtengine::ControlLine::Type type;
    std::shared_ptr<OPIcon> icon;

    if (dx > dy) { // More horizontal than vertical.
        type = rtengine::ControlLine::HORIZONTAL;
        icon = line.icon_h;
    } else {
        type = rtengine::ControlLine::VERTICAL;
        icon = line.icon_v;
    }

    if (type != line.type) { // Need to update line type.
        line.type = type;
        line.icon = icon;
        visibleGeometry[line_id * ::ControlLine::OBJ_COUNT + 1] =
            line.icon.get();
    }
}

void ControlLineManager::removeAll(void)
{
    visibleGeometry.clear();
    mouseOverGeometry.erase(mouseOverGeometry.begin() + 1,
                            mouseOverGeometry.end());
    control_lines.clear();
    prev_obj = -1;
    selected_object = -1;
    edited = true;
    callbacks->lineChanged();
}

void ControlLineManager::removeLine(size_t line_id)
{
    if (line_id >= control_lines.size()) {
        return;
    }

    visibleGeometry.erase(
        visibleGeometry.begin() + ::ControlLine::OBJ_COUNT * line_id,
        visibleGeometry.begin() + ::ControlLine::OBJ_COUNT * line_id
        + ::ControlLine::OBJ_COUNT
    );
    mouseOverGeometry.erase(
        mouseOverGeometry.begin() + ::ControlLine::OBJ_COUNT * line_id + 1,
        mouseOverGeometry.begin() + ::ControlLine::OBJ_COUNT * line_id
        + ::ControlLine::OBJ_COUNT + 1
    );
    control_lines.erase(control_lines.begin() + line_id);

    edited = true;
    callbacks->lineChanged();
}

void ControlLineManager::toControlLines(std::vector<rtengine::ControlLine>&
                                        converted) const
{
    converted.clear();
    converted.resize(control_lines.size());

    for (unsigned int i = 0; i < control_lines.size(); i++) {
        converted[i].x1 = control_lines[i]->begin->center.x;
        converted[i].y1 = control_lines[i]->begin->center.y;
        converted[i].x2 = control_lines[i]->end->center.x;
        converted[i].y2 = control_lines[i]->end->center.y;
        converted[i].type = control_lines[i]->type;
    }
}
