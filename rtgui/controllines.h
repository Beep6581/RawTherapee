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
#pragma once

#include <memory>

#include "editcallbacks.h"
#include "../rtengine/perspectivecorrection.h"

class Circle;
class Line;
class OPIcon;
class Rectangle;
class RTSurface;

struct ControlLine {
    static constexpr int OBJ_COUNT = 4;
    std::unique_ptr<Line> line;
    std::shared_ptr<OPIcon> icon;
    std::shared_ptr<OPIcon> icon_h, icon_v;
    std::unique_ptr<Circle> begin, end;
    rtengine::ControlLine::Type type;

    ~ControlLine();
};

class ControlLineManager: EditSubscriber
{

protected:
    /** Hidden object for capturing mouse events. */
    std::unique_ptr<Rectangle> canvas_area;
    rtengine::Coord drag_delta;
    std::vector<std::unique_ptr<ControlLine>> control_lines;
    CursorShape cursor;
    bool draw_mode;
    bool drawing_line;
    bool edited;
    Cairo::RefPtr<RTSurface> line_icon_h, line_icon_v;
    Cairo::RefPtr<RTSurface> line_icon_h_prelight, line_icon_v_prelight;
    int prev_obj;
    int selected_object;

    void addLine(rtengine::Coord begin, rtengine::Coord end,
                 rtengine::ControlLine::Type type = rtengine::ControlLine::VERTICAL);
    /**
     * Set the line type of the line containing the object according to the
     * line's angle.
     *
     * If the line is within 45 degrees of a perfectly vertical
     * line, inclusive, the line type is set to vertical. Otherwise, horizontal.
     */
    void autoSetLineType(int object_id);
    void removeLine(size_t line_id);

public:
    class Callbacks
    {
    public:
        virtual ~Callbacks() {};
        /** Called when a line changed (added, removed, moved, etc.). */
        virtual void lineChanged(void) {};
        /** Called when the EditSubscriber's switchOffEditMode is called. */
        virtual void switchOffEditMode(void) {};
    };

    /** Callbacks to invoke. */
    std::shared_ptr<Callbacks> callbacks;

    ControlLineManager();
    ~ControlLineManager();

    bool getEdited(void) const;
    void removeAll(void);
    /** Sets whether or not the lines are visible and interact-able. */
    void setActive(bool active);
    /** Set whether or not lines can be drawn and deleted. */
    void setDrawMode(bool draw);
    void setEdited(bool edited);
    void setEditProvider(EditDataProvider* provider);
    void setLines(const std::vector<rtengine::ControlLine>& lines);
    /** Returns the number of lines. */
    size_t size(void) const;
    /**
     * Allocates a new array and populates it with copies of the control lines.
     */
    void toControlLines(std::vector<rtengine::ControlLine>& converted) const;

    // EditSubscriber overrides
    bool button1Pressed(int modifierKey) override;
    bool button1Released(void) override;
    bool button3Pressed(int modifierKey) override;
    bool pick1(bool picked) override;
    bool pick3(bool picked) override;
    bool drag1(int modifierKey) override;
    CursorShape getCursor(int objectID) const override;
    bool mouseOver(int modifierKey) override;
    void switchOffEditMode(void) override;
};
