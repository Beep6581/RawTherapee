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
#include "rtengine.h"
struct IMFILE {
    int fd;
    ssize_t pos;
    ssize_t size;
    char* data;
    bool eof;
    rtengine::ProgressListener *plistener;
    double progress_range;
    ssize_t progress_next;
    ssize_t progress_current;
};

/*
  Functions for progress bar updates
  Note: progress bar is not intended to be exact, eg if you read same data over and over again progress
  will potentially reach 100% before you're finished.
 */
void imfile_set_plistener(IMFILE *f, rtengine::ProgressListener *plistener, double progress_range);
void imfile_update_progress(IMFILE *f);

IMFILE* fopen (const char* fname);
IMFILE* gfopen (const char* fname);
IMFILE* fopen (unsigned* buf, int size);
void fclose (IMFILE* f);
inline int ftell (IMFILE* f)
{

    return f->pos;
}

inline int feof (IMFILE* f)
{

    return f->eof;
}

inline void fseek (IMFILE* f, int p, int how)
{
    int fpos = f->pos;

    if (how == SEEK_SET) {
        f->pos = p;
    } else if (how == SEEK_CUR) {
        f->pos += p;
    } else if (how == SEEK_END) {
        if(p <= 0 && -p <= f->size) {
            f->pos = f->size + p;
        }
        return;
    }

    if (f->pos < 0  || f->pos > f->size) {
        f->pos = fpos;
    }
}

inline int fgetc (IMFILE* f)
{

    if (LIKELY(f->pos < f->size)) {
        if (f->plistener && ++f->progress_current >= f->progress_next) {
            imfile_update_progress(f);
        }

        return (unsigned char)f->data[f->pos++];
    }

    f->eof = true;
    return EOF;
}

inline int getc (IMFILE* f)
{

    return fgetc(f);
}

inline int fread (void* dst, int es, int count, IMFILE* f)
{

    int s = es * count;
    int avail = f->size - f->pos;

    if (s <= avail) {
        memcpy (dst, f->data + f->pos, s);
        f->pos += s;

        if (f->plistener) {
            f->progress_current += s;

            if (f->progress_current >= f->progress_next) {
                imfile_update_progress(f);
            }
        }

        return count;
    } else {
        memcpy (dst, f->data + f->pos, avail);
        f->pos += avail;
        f->eof = true;
        return avail / es;
    }
}

inline unsigned char* fdata(int offset, IMFILE* f)
{
    return (unsigned char*)f->data + offset;
}

int fscanf (IMFILE* f, const char* s ...);
char* fgets (char* s, int n, IMFILE* f);

#endif

