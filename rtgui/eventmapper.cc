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

#include "eventmapper.h"


ProcEventMapper::ProcEventMapper()
{
    for (int event = 0; event < rtengine::NUMOFEVENTS; ++event) {
        history_msgs_[event] = "HISTORY_MSG_" + std::to_string(event + 1);
    }
}


ProcEventMapper *ProcEventMapper::getInstance()
{
    static ProcEventMapper instance;
    return &instance;
}


rtengine::ProcEvent ProcEventMapper::newEvent(int action, const std::string &history_msg)
{
    rtengine::ProcEvent event = rtengine::RefreshMapper::getInstance()->newEvent();
    rtengine::RefreshMapper::getInstance()->mapEvent(event, action);    

    if (history_msg.empty()) {
        history_msgs_[event] = "HISTORY_MSG_" + std::to_string(event + 1);
    } else {
        history_msgs_[event] = history_msg;
    }
    
    return event;
}


const std::string &ProcEventMapper::getHistoryMsg(rtengine::ProcEvent event) const
{
    static std::string empty;
    auto it = history_msgs_.find(event);
    if (it == history_msgs_.end()) {
        return empty;
    } else {
        return it->second;
    }
}
