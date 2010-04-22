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
#include <myfile.h>
#include <cstdarg>
#include <glibmm.h>
#ifdef RAWZOR_SUPPORT
#include <rwz_sdk.h>
#endif

IMFILE* fopen (const char* fname) {

	FILE* f = fopen (fname, "rb");
    if (!f)
        return NULL;
    IMFILE* mf = new IMFILE;
	fseek (f, 0, SEEK_END);
	mf->size = ftell (f);
	mf->data = new char [mf->size];
	fseek (f, 0, SEEK_SET);
	fread (mf->data, 1, mf->size, f);
	fclose (f);
	mf->pos = 0;
	mf->eof = false;
#ifdef RAWZOR_SUPPORT
    // RAWZOR support begin
    bool rawzor = false;
    Glib::ustring bname = Glib::path_get_basename(fname);
    int lastdot = bname.find_last_of ('.');
    if (lastdot!=bname.npos)
        rawzor = bname.substr (lastdot).casefold() == Glib::ustring(".rwz").casefold();

    if (rawzor) {
        int realSize = 0;
        if (!m_rwz_check (mf->data, mf->size, &realSize)) {
            char* realData = new char [realSize];
            m_rwz_decompress (mf->data, mf->size, realData, realSize);
            delete [] mf->data;
            mf->data = realData;
            mf->size = realSize;
        }
    }
    // RAWZOR support end
#endif

	return mf;
}

IMFILE* gfopen (const char* fname) {

	FILE* f = g_fopen (fname, "rb");
    if (!f)
        return NULL;
    IMFILE* mf = new IMFILE;
	fseek (f, 0, SEEK_END);
	mf->size = ftell (f);
	mf->data = new char [mf->size];
	fseek (f, 0, SEEK_SET);
	fread (mf->data, 1, mf->size, f);
	fclose (f);
	mf->pos = 0;
	mf->eof = false;

#ifdef RAWZOR_SUPPORT
    // RAWZOR support begin
    bool rawzor = false;
    Glib::ustring bname = Glib::path_get_basename(fname);
    int lastdot = bname.find_last_of ('.');
    if (lastdot!=bname.npos)
        rawzor = bname.substr (lastdot).casefold() == Glib::ustring(".rwz").casefold();

    if (rawzor) {
        int realSize = 0;
        if (!m_rwz_check (mf->data, mf->size, &realSize)) {
            char* realData = new char [realSize];
            m_rwz_decompress (mf->data, mf->size, realData, realSize);
            delete [] mf->data;
            mf->data = realData;
            mf->size = realSize;
        }
    }
    // RAWZOR support end
#endif
	return mf;
}

IMFILE* fopen (unsigned* buf, int size) {

	IMFILE* mf = new IMFILE;
	mf->size = size;
	mf->data = new char [mf->size];
	memcpy (mf->data, buf, size);
	mf->pos = 0;
	mf->eof = false;
	return mf;
}

void fclose (IMFILE* f) {

	delete [] f->data;
	delete f;
}

int fscanf (IMFILE* f, const char* s ...) {

	va_list ap;
	return sscanf (f->data, s, ap);
}

char* fgets (char* s, int n, IMFILE* f) {

	if (f->pos>=f->size) {
		f->eof = true;
		return NULL;
	}
	int i = 0;
	do s[i++] = f->data[f->pos++];		
	while (i<n && f->pos<f->size);
	return s;
}
