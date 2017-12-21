/* -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio <alberto.griggio@gmail.com>
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
#pragma once

#include <string>
#include <unordered_map>
#include "../rtengine/refreshmap.h"


class ProcEventMapper {
public:
    static ProcEventMapper *getInstance();
    rtengine::ProcEvent newEvent(int action, const std::string &history_msg="");
    const std::string &getHistoryMsg(rtengine::ProcEvent event) const;

private:
    ProcEventMapper();
    
    std::unordered_map<int, std::string> history_msgs_;
};
