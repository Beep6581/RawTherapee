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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _COLOREDBAR_
#define _COLOREDBAR_

#include "colorprovider.h"
#include "guiutils.h"

/*
 * Parent class for all colored bar type; a ColorProvider has to be set
 * thanks to "setColorProvider" to be able to display colors inside the bar
 */
class ColoredBar : public BackBuffer, public ColorCaller
{

private:
    void draw();

protected:
    eRTOrientation orientation;
    std::vector<GradientMilestone> bgGradient;

public:
    explicit ColoredBar (eRTOrientation orient);

    void expose(Glib::RefPtr<Gdk::Window> destWindow);
    void expose(Cairo::RefPtr<Cairo::ImageSurface> destSurface);
    void expose(BackBuffer *backBuffer);

    bool canGetColors();

    // Method for convenience; if no Gradient provided, the ColoredBar will ask colors on a per pixel basis
    void setBgGradient (const std::vector<GradientMilestone> &milestones);
    // by clearing the gradient, the ColorProvider will have to provide colors on a per pixel basis if a ColorProvider
    // has been set, through ColorProvider::colorForValue on next ColoredBar::expose
    void clearBgGradient ();
};

#endif
