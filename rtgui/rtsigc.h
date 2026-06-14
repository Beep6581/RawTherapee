/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2024 Daniel Gao <daniel.gao.work@gmail.com>
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

#include <sigc++/connection.h>

class RtScopedConnection
{
public:
    RtScopedConnection() = default;
    RtScopedConnection(sigc::connection&& conn) noexcept : m_connection(std::move(conn)) {}

    RtScopedConnection(const RtScopedConnection&) = delete;
    RtScopedConnection(RtScopedConnection&& other)
    {
        m_connection.disconnect();
        m_connection = other.m_connection;
        other.m_connection = {};
    }

    ~RtScopedConnection() { m_connection.disconnect(); }

    RtScopedConnection& operator=(const RtScopedConnection&) = delete;
    RtScopedConnection& operator=(RtScopedConnection&& other)
    {
        m_connection.disconnect();
        m_connection = other.m_connection;
        other.m_connection = {};
        return *this;
    }

private:
    sigc::connection m_connection;
};
