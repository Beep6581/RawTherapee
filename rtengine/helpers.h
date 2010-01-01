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
#include <string>
#include <procparams.h>
#include <vector>

using namespace rtengine;
using namespace rtengine::procparams;

std::string*   create_stdstring (char* txt);
Glib::ustring* create_gustring (char* txt);
void set_stdstring (std::string& s, char* txt);
void set_gustring  (Glib::ustring& s, char* txt);
char* get_stdstring (std::string& s);
char* get_gustring  (Glib::ustring& s);

void clear_double_vector  (std::vector<double>& d);
void append_double_vector (std::vector<double>& d, double dd);
void resize_double_vector (std::vector<double>& d, int size);
double get_double_vector (std::vector<double>& d, int i);
void set_double_vector (std::vector<double>& d, int i, double x);
int size_double_vector (std::vector<double>& d);

void clear_gustring_vector  (std::vector<Glib::ustring>& d);
void append_gustring_vector (std::vector<Glib::ustring>& d, const Glib::ustring& dd);
void resize_gustring_vector (std::vector<Glib::ustring>& d, int size);
Glib::ustring& get_gustring_vector (std::vector<Glib::ustring>& d, int i);
void set_gustring_vector (std::vector<Glib::ustring>& d, int i, const Glib::ustring& x);
int size_gustring_vector (std::vector<Glib::ustring>& d);

void clear_exifpair_vector  (std::vector<ExifPair>& d);
void append_exifpair_vector (std::vector<ExifPair>& d, const ExifPair& dd);
void resize_exifpair_vector (std::vector<ExifPair>& d, int size);
ExifPair& get_exifpair_vector (std::vector<ExifPair>& d, int i);
void set_exifpair_vector (std::vector<ExifPair>& d, int i, const ExifPair& x);
int size_exifpair_vector (std::vector<ExifPair>& d);

void clear_iptcpair_vector  (std::vector<IPTCPair>& d);
void append_iptcpair_vector (std::vector<IPTCPair>& d, const IPTCPair& dd);
void resize_iptcpair_vector (std::vector<IPTCPair>& d, int size);
IPTCPair& get_iptcpair_vector (std::vector<IPTCPair>& d, int i);
void set_iptcpair_vector (std::vector<IPTCPair>& d, int i, const IPTCPair& x);
int size_iptcpair_vector (std::vector<IPTCPair>& d);

int sizeof_stdstring ();
int sizeof_gustring ();
int sizeof_double_vector ();
int sizeof_gustring_vector ();
int sizeof_exifpair_vector ();
int sizeof_iptcpair_vector ();

