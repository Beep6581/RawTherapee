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
#include <safegtk.h>
#ifdef RAWZOR_SUPPORT
#include <rwz_sdk.h>
#endif

// get mmap() sorted out
#ifdef MYFILE_MMAP

#ifdef WIN32

#include <fcntl.h>
#include <windows.h>

// dummy values
#define MAP_PRIVATE 1
#define PROT_READ 1

void* mmap(void *start, size_t length, int prot, int flags, int fd, off_t offset)
{
	HANDLE handle = CreateFileMapping((HANDLE)_get_osfhandle(fd), NULL, PAGE_WRITECOPY, 0, 0, NULL);

	if (handle != NULL) {
		start = MapViewOfFile(handle, FILE_MAP_COPY, 0, offset, length);
		CloseHandle(handle);
	}

	return start;
}

int munmap(void *start, size_t length)
{
	UnmapViewOfFile(start);
	return 0;
}

#else // WIN32

#include <fcntl.h>
#include <sys/mman.h>

#endif // WIN32
#endif // MYFILE_MMAP

#ifdef MYFILE_MMAP

IMFILE* fopen (const char* fname)
{
	int fd = safe_open_ReadOnly(fname);
	if ( fd < 0 )
		return 0;

	struct stat stat_buffer;
	if ( fstat(fd,&stat_buffer) < 0 )
	{
		printf("no stat\n");
		close(fd);
		return 0;
	}

	void* data = mmap(0,stat_buffer.st_size,PROT_READ,MAP_PRIVATE,fd,0);
	if ( data == 0 )
	{
		printf("no mmap\n");
		close(fd);
		return 0;
	}

	IMFILE* mf = new IMFILE;

	mf->fd = fd;
	mf->pos = 0;
	mf->size = stat_buffer.st_size;
	mf->data = (char*)data;
	mf->eof = false;

	return mf;
}

IMFILE* gfopen (const char* fname)
{
	return fopen(fname);
}
#else

IMFILE* fopen (const char* fname) {

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
#endif //MYFILE_MMAP

IMFILE* fopen (unsigned* buf, int size) {

	IMFILE* mf = new IMFILE;
	mf->fd = -1;
	mf->size = size;
	mf->data = new char [mf->size];
	memcpy ((void*)mf->data, buf, size);
	mf->pos = 0;
	mf->eof = false;
	return mf;
}

void fclose (IMFILE* f) {
#ifdef MYFILE_MMAP
	if ( f->fd == -1 )
	{
		delete [] f->data;
	}
	else
	{
		munmap((void*)f->data,f->size);
		close(f->fd);
	}
#else
	delete [] f->data;
#endif
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
