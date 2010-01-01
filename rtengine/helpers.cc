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
#include <helpers.h>
#include <string.h>

using namespace rtengine;
using namespace rtengine::procparams;

std::string*   create_stdstring (char* txt) { return new std::string (txt); }
Glib::ustring* create_gustring (char* txt)  { return new Glib::ustring (txt); }
void set_stdstring (std::string& s, char* txt)   { s = txt; }
void set_gustring  (Glib::ustring& s, char* txt) { s = txt; }
char* get_stdstring (std::string& s)   { return strdup (s.c_str()); }
char* get_gustring  (Glib::ustring& s) { return strdup (s.c_str()); }

void clear_double_vector  (std::vector<double>& d) { d.clear(); }
void append_double_vector (std::vector<double>& d, double dd) { d.push_back(dd); }
void resize_double_vector (std::vector<double>& d, int size) { d.resize (size); }
double get_double_vector (std::vector<double>& d, int i) { return d[i]; }
void set_double_vector (std::vector<double>& d, int i, double x) { d[i] = x; }
int size_double_vector (std::vector<double>& d) { return d.size(); }

void clear_gustring_vector  (std::vector<Glib::ustring>& d) { d.clear(); }
void append_gustring_vector (std::vector<Glib::ustring>& d, const Glib::ustring& dd) { d.push_back(dd); }
void resize_gustring_vector (std::vector<Glib::ustring>& d, int size) { d.resize (size); }
Glib::ustring& get_gustring_vector (std::vector<Glib::ustring>& d, int i) { return d[i]; }
void set_gustring_vector (std::vector<Glib::ustring>& d, int i, const Glib::ustring& x) { d[i] = x; }
int size_gustring_vector (std::vector<Glib::ustring>& d) { return d.size(); }

void clear_exifpair_vector  (std::vector<ExifPair>& d) { d.clear(); }
void append_exifpair_vector (std::vector<ExifPair>& d, const ExifPair& dd) { d.push_back(dd); }
void resize_exifpair_vector (std::vector<ExifPair>& d, int size) { d.resize (size); }
ExifPair& get_exifpair_vector (std::vector<ExifPair>& d, int i) { return d[i]; }
void set_exifpair_vector (std::vector<ExifPair>& d, int i, const ExifPair& x) { d[i] = x; }
int size_exifpair_vector (std::vector<ExifPair>& d) { return d.size(); }

void clear_iptcpair_vector  (std::vector<IPTCPair>& d) { d.clear(); }
void append_iptcpair_vector (std::vector<IPTCPair>& d, const IPTCPair& dd) { d.push_back(dd); }
void resize_iptcpair_vector (std::vector<IPTCPair>& d, int size) { d.resize (size); }
IPTCPair& get_iptcpair_vector (std::vector<IPTCPair>& d, int i) { return d[i]; }
void set_iptcpair_vector (std::vector<IPTCPair>& d, int i, const IPTCPair& x) { d[i] = x; }
int size_iptcpair_vector (std::vector<IPTCPair>& d) { return d.size(); }

int sizeof_stdstring () { return sizeof(std::string); }
int sizeof_gustring ()  { return sizeof(Glib::ustring); }
int sizeof_double_vector ()  { return sizeof(std::vector<double>); }
int sizeof_gustring_vector ()  { return sizeof(std::vector<Glib::ustring>); }
int sizeof_exifpair_vector ()  { return sizeof(std::vector<ExifPair>); }
int sizeof_iptcpair_vector ()  { return sizeof(std::vector<IPTCPair>); }


