/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Romei <aldrop8@gmail.com>
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

#include <array>

#include <gtkmm.h>

#include "adjuster.h"
#include "editcallbacks.h"
#include "guiutils.h"
#include "toolpanel.h"

class FilmNegProvider
{
public:
    virtual ~FilmNegProvider() = default;

    virtual bool getFilmNegativeExponents(rtengine::Coord spotA, rtengine::Coord spotB, std::array<float, 3>& newExps) = 0;
    virtual bool getRawSpotValues(rtengine::Coord spot, int spotSize, std::array<float, 3>& rawValues) = 0;
};

class FilmNegative final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public EditSubscriber,
    public rtengine::FilmNegListener
{
public:
    FilmNegative();
    ~FilmNegative() override;

    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode(bool batchMode) override;

    void adjusterChanged(Adjuster* a, double newval) override;
    void enabledChanged() override;

    void filmBaseValuesChanged(std::array<float, 3> rgb) override;

    void setFilmNegProvider(FilmNegProvider* provider);

    void setEditProvider(EditDataProvider* provider) override;

    // EditSubscriber interface
    CursorShape getCursor(int objectID) const override;
    bool mouseOver(int modifierKey) override;
    bool button1Pressed(int modifierKey) override;
    bool button1Released() override;
    void switchOffEditMode() override;

private:
    void editToggled();
    void baseSpotToggled();

    const rtengine::ProcEvent evFilmNegativeExponents;
    const rtengine::ProcEvent evFilmNegativeEnabled;
    const rtengine::ProcEvent evFilmBaseValues;

    std::vector<rtengine::Coord> refSpotCoords;

    std::array<float, 3> filmBaseValues;

    FilmNegProvider* fnp;

    Adjuster* const greenExp;
    Adjuster* const redRatio;
    Adjuster* const blueRatio;

    Gtk::Grid* const spotgrid;
    Gtk::ToggleButton* const spotbutton;

    Gtk::Label* const filmBaseLabel;
    Gtk::Label* const filmBaseValuesLabel;
    Gtk::ToggleButton* const filmBaseSpotButton;

};
