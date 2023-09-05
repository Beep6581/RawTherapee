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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "myfile.h"
#include <cstdarg>
#include "rtengine.h"
// get mmap() sorted out
#ifdef MYFILE_MMAP

#ifdef _WIN32

#include <fcntl.h>
#include <windows.h>

// dummy values
#define MAP_PRIVATE 1
#define PROT_READ 1
#define MAP_FAILED (void *)-1

#ifdef __GNUC__ // silence warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

void* mmap(void *start, size_t length, int prot, int flags, int fd, off_t offset)
{
    HANDLE handle = CreateFileMapping((HANDLE)_get_osfhandle(fd), NULL, PAGE_WRITECOPY, 0, 0, NULL);

    if (handle != NULL) {
        start = MapViewOfFile(handle, FILE_MAP_COPY, 0, offset, length);
        CloseHandle(handle);
        return start;
    }

    return MAP_FAILED;
}

int munmap(void *start, size_t length)
{
    UnmapViewOfFile(start);
    return 0;
}
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#else // _WIN32

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

#endif // _WIN32
#endif // MYFILE_MMAP

#ifdef MYFILE_MMAP

rtengine::IMFILE* rtengine::fopen (const char* fname)
{
    int fd;

#ifdef _WIN32

    fd = -1;
    // First convert UTF8 to UTF16, then use Windows function to open the file and convert back to file descriptor.
    std::unique_ptr<wchar_t, GFreeFunc> wfname (reinterpret_cast<wchar_t*>(g_utf8_to_utf16 (fname, -1, NULL, NULL, NULL)), g_free);

    HANDLE hFile = CreateFileW (wfname.get (), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile != INVALID_HANDLE_VALUE) {
        fd = _open_osfhandle((intptr_t)hFile, 0);
    }

#else

    fd = ::g_open (fname, O_RDONLY);

#endif

    if ( fd < 0 ) {
        return nullptr;
    }

    struct stat stat_buffer;

    if ( fstat(fd, &stat_buffer) < 0 ) {
        printf("no stat\n");
        close (fd);
        return nullptr;
    }

    void* data = mmap(nullptr, stat_buffer.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

    if ( data == MAP_FAILED ) {
        printf("no mmap %s\n", fname);
        close(fd);
        return nullptr;
    }

    IMFILE* mf = new IMFILE;

    memset(mf, 0, sizeof(*mf));
    mf->fd = fd;
    mf->pos = 0;
    mf->size = stat_buffer.st_size;
    mf->data = (char*)data;
    mf->eof = false;

    return mf;
}

rtengine::IMFILE* rtengine::gfopen (const char* fname)
{
    return fopen(fname);
}
#else

rtengine::IMFILE* rtengine::fopen (const char* fname)
{

    FILE* f = g_fopen (fname, "rb");

    if (!f) {
        return NULL;
    }

    IMFILE* mf = new IMFILE;
    memset(mf, 0, sizeof(*mf));
    fseek (f, 0, SEEK_END);
    mf->size = ftell (f);
    mf->data = new char [mf->size];
    fseek (f, 0, SEEK_SET);
    fread (mf->data, 1, mf->size, f);
    fclose (f);
    mf->pos = 0;
    mf->eof = false;

    return mf;
}

rtengine::IMFILE* rtengine::gfopen (const char* fname)
{

    FILE* f = g_fopen (fname, "rb");

    if (!f) {
        return NULL;
    }

    IMFILE* mf = new IMFILE;
    memset(mf, 0, sizeof(*mf));
    fseek (f, 0, SEEK_END);
    mf->size = ftell (f);
    mf->data = new char [mf->size];
    fseek (f, 0, SEEK_SET);
    fread (mf->data, 1, mf->size, f);
    fclose (f);
    mf->pos = 0;
    mf->eof = false;

    return mf;
}
#endif //MYFILE_MMAP

rtengine::IMFILE* rtengine::fopen (unsigned* buf, int size)
{

    IMFILE* mf = new IMFILE;
    memset(mf, 0, sizeof(*mf));
    mf->fd = -1;
    mf->size = size;
    mf->data = new char [mf->size];
    memcpy ((void*)mf->data, buf, size);
    mf->pos = 0;
    mf->eof = false;
    return mf;
}

void rtengine::fclose (IMFILE* f)
{
#ifdef MYFILE_MMAP

    if ( f->fd == -1 ) {
        delete [] f->data;
    } else {
        munmap((void*)f->data, f->size);
        close(f->fd);
    }

#else
    delete [] f->data;
#endif
    delete f;
}

int rtengine::fscanf (IMFILE* f, const char* s ...)
{
    // fscanf not easily wrapped since we have no terminating \0 at end
    // of file data and vsscanf() won't tell us how many characters that
    // were parsed. However, only dcraw.cc code use it and only for "%f" and
    // "%d", so we make a dummy fscanf here just to support dcraw case.
    char buf[51], *endptr = nullptr;
    int copy_sz = f->size - f->pos;

    if (copy_sz >= static_cast<int>(sizeof(buf))) {
        copy_sz = sizeof(buf) - 1;
    }

    memcpy(buf, &f->data[f->pos], copy_sz);
    buf[copy_sz] = '\0';
    va_list ap;
    va_start (ap, s);

    if (strcmp(s, "%d") == 0) {
        int i = strtol(buf, &endptr, 10);

        if (endptr == buf) {
            va_end (ap);
            return 0;
        }

        int *pi = va_arg(ap, int*);
        *pi = i;
    } else if (strcmp(s, "%f") == 0) {
        float f = strtof(buf, &endptr);

        if (endptr == buf) {
            va_end (ap);
            return 0;
        }

        float *pf = va_arg(ap, float*);
        *pf = f;
    }

    va_end (ap);
    f->pos += endptr - buf;
    return 1;
}


char* rtengine::fgets (char* s, int n, IMFILE* f)
{

    if (f->pos >= f->size) {
        f->eof = true;
        return nullptr;
    }

    int i = 0;

    do {
        s[i++] = f->data[f->pos++];
    } while (i < n && f->pos < f->size);

    return s;
}

void rtengine::imfile_set_plistener(IMFILE *f, rtengine::ProgressListener *plistener, double progress_range)
{
    f->plistener = plistener;
    f->progress_range = progress_range;
    f->progress_next = f->size / 10 + 1;
    f->progress_current = 0;
}

void rtengine::imfile_update_progress(IMFILE *f)
{
    if (!f->plistener || f->progress_current < f->progress_next) {
        return;
    }

    do {
        f->progress_next += f->size / 10 + 1;
    } while (f->progress_next < f->progress_current);

    double p = (double)f->progress_current / f->size;

    if (p > 1.0) {
        /* this can happen if same bytes are read over and over again. Progress bar is not intended
           to be exact, just give some progress indication for normal raw file access patterns */
        p = 1.0;
    }

    f->plistener->setProgress(p * f->progress_range);
}
