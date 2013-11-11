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
#ifndef _IPTCPAIRS_
#define _IPTCPAIRS_


struct IptcPair {
    IptcTag tag;
    size_t size;
    Glib::ustring field;
};

const IptcPair strTags[] = {
	{IPTC_TAG_CAPTION, 2000, "Caption"},
	{IPTC_TAG_WRITER_EDITOR, 32, "CaptionWriter"},
	{IPTC_TAG_HEADLINE, 256, "Headline"},
	{IPTC_TAG_SPECIAL_INSTRUCTIONS, 256, "Instructions"},
	{IPTC_TAG_CATEGORY, 3, "Category"},
	{IPTC_TAG_BYLINE, 32, "Author"},
	{IPTC_TAG_BYLINE_TITLE, 32, "AuthorsPosition"},
	{IPTC_TAG_CREDIT, 32, "Credit"},
	{IPTC_TAG_SOURCE, 32, "Source"},
	{IPTC_TAG_COPYRIGHT_NOTICE, 128, "Copyright"},
	{IPTC_TAG_CITY, 32, "City"},
	{IPTC_TAG_STATE, 32, "Province"},
	{IPTC_TAG_COUNTRY_NAME, 64, "Country"},
	{IPTC_TAG_OBJECT_NAME, 64, "Title"},
	{IPTC_TAG_ORIG_TRANS_REF, 32, "TransReference"},
	{IPTC_TAG_DATE_CREATED, 8, "DateCreated"}
};

#endif

