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
#ifndef _MYFILE_
#define _MYFILE_

#include <glib/gstdio.h>
#include <cstdio>
#include <cstring>
struct IMFILE {
	int fd;
	int pos;
	int size;
	char* data;
	bool eof;
};

IMFILE* fopen (const char* fname);
IMFILE* gfopen (const char* fname);IMFILE* fopen (unsigned* buf, int size);
void fclose (IMFILE* f);
inline int ftell (IMFILE* f) {

	return f->pos;
}

inline int feof (IMFILE* f) {

	return f->eof;
}

inline void fseek (IMFILE* f, int p, int how) {
	int fpos = f->pos;

	if (how==SEEK_SET)
		f->pos = p;
	else if (how==SEEK_CUR)
		f->pos += p;
	else if (how==SEEK_END)
		f->pos = f->size-p;

	if (f->pos < 0  || f->pos> f->size)
		f->pos = fpos;
}

inline int fgetc (IMFILE* f) {

	if (f->pos<f->size)
		return (unsigned char)f->data[f->pos++];
	f->eof = true;
	return EOF;
}

inline int getc (IMFILE* f) {

	if (f->pos<f->size)
		return (unsigned char)f->data[f->pos++];
	f->eof = true;
	return EOF;
}

inline int fread (void* dst, int es, int count, IMFILE* f) {

	int s = es*count;
	int avail = f->size - f->pos;
	if (s<=avail) {
		memcpy (dst, f->data+f->pos, s);
		f->pos += s;
		return count;
	}
	else {
		memcpy (dst, f->data+f->pos, avail);
		f->pos += avail;
		f->eof = true;
		return avail/es;
	}
}

inline unsigned char* fdata(int offset, IMFILE* f) {
	return (unsigned char*)f->data + offset;
}

int fscanf (IMFILE* f, const char* s ...);
char* fgets (char* s, int n, IMFILE* f);
#endif

