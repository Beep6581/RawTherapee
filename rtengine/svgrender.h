/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2025 Daniel Gao <daniel.gao.work@gmail.com>
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

#include <cairomm/surface.h>

#include <stdexcept>
#include <string>

namespace rtengine {

class SvgRenderException : public std::runtime_error {
public:
    explicit SvgRenderException(const std::string& message)
        : std::runtime_error(message) {}
};

Cairo::RefPtr<Cairo::ImageSurface> renderSvg(const std::string& data, int logical_width,
                                             int logical_height, int device_scale);

} // namespace rtengine
