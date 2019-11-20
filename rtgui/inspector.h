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
#pragma once

#include <gtkmm.h>

#include "guiutils.h"

#include "../rtengine/coord.h"
#include "../rtengine/coord2d.h"

class InspectorBuffer
{
//private:
//    int infoFromImage (const Glib::ustring& fname);

public:
    BackBuffer imgBuffer;
    Glib::ustring imgPath;
    int currTransform;  // coarse rotation from RT, not from shot orientation
    bool fromRaw;

    explicit InspectorBuffer(const Glib::ustring &imgagePath);
    //~InspectorBuffer();
};

class Inspector : public Gtk::DrawingArea
{

private:
    rtengine::Coord center;
    std::vector<InspectorBuffer*> images;
    InspectorBuffer* currImage;
    double zoom;
    bool active;

    sigc::connection delayconn;
    Glib::ustring next_image_path;

    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr) override;
    void deleteBuffers();

    bool doSwitchImage();

public:
    Inspector();
    ~Inspector() override;

    /** @brief Mouse movement to a new position
     * @param pos Location of the mouse, in percentage (i.e. [0;1] range) relative to the full size image ; -1,-1 == out of the image
     * @param transform H/V flip and coarse rotation transformation
     */
    void mouseMove (rtengine::Coord2D pos, int transform);

    /** @brief A new image is being flown over
     * @param fullPath Full path of the image that is being hovered inspect, or an empty string if out of any image.
     */
    void switchImage (const Glib::ustring &fullPath);

    /** @brief Set the new coarse rotation transformation
     * @param transform A semi-bitfield coarse transformation using #defines from iimage.h
     */
    void setTransformation (int transform);

    /** @brief Use this method to flush all image buffer whenever the Inspector panel is hidden
     */
    void flushBuffers ();

    /** @brief Set the inspector on/off
     * @param state true if to activate the Inspector, false to disable it and flush the buffers
     */
    void setActive(bool state);

    /** @brief Get the on/off state
     */
    bool isActive() const
    {
        return active;
    };

    Gtk::SizeRequestMode get_request_mode_vfunc () const override;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const override;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const override;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const override;
    void get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const override;

};
