/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2021 Lawrence Lee
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

#include <gtkmm/box.h>

class Options;

/**
 * Widget for configuring the location of tools in the tool panel tabs.
 */
class ToolLocationPreference : public Gtk::Box
{
private:
    struct Impl;
    std::unique_ptr<Impl> impl;

public:
    /**
     * Constructs a tool location preference widget that gets values from and
     * updates the provided options object.
     */
    explicit ToolLocationPreference(Options &options);
    /**
     * Updates the options object associated with this object with the current
     * favorites preferences.
     */
    void updateOptions();
};
