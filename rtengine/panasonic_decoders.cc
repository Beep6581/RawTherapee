/*
 *  This file is part of RawTherapee.
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

#include <iostream>
#include "dcraw.h"

// Code adapted from libraw
/* -*- C++ -*-
 * Copyright 2019 LibRaw LLC (info@libraw.org)
 *
 LibRaw is free software; you can redistribute it and/or modify
 it under the terms of the one of two licenses as you choose:

1. GNU LESSER GENERAL PUBLIC LICENSE version 2.1
   (See file LICENSE.LGPL provided in LibRaw distribution archive for details).

2. COMMON DEVELOPMENT AND DISTRIBUTION LICENSE (CDDL) Version 1.0
   (See file LICENSE.CDDL provided in LibRaw distribution archive for details).

*/

unsigned DCraw::pana_bits_t::operator() (int nbits, unsigned *bytes)
{
    int byte;

    if (!nbits && !bytes) {
        return vbits=0;
    }
    if (!vbits) {
        fread (buf+load_flags, 1, 0x4000-load_flags, ifp);
        fread (buf, 1, load_flags, ifp);
    }
    if (encoding == 5) {
        for (byte = 0; byte < 16; byte++)
        {
          bytes[byte] = buf[vbits++];
          vbits &= 0x3FFF;
        }
        return 0;
    } else {
        vbits = (vbits - nbits) & 0x1ffff;
        byte = vbits >> 3 ^ 0x3ff0;
        return (buf[byte] | buf[byte+1] << 8) >> (vbits & 7) & ~(-1 << nbits);
    }
}

class pana_cs6_page_decoder
{
    unsigned int pixelbuffer[14], lastoffset, maxoffset;
    unsigned char current, *buffer;
public:
    pana_cs6_page_decoder(unsigned char *_buffer, unsigned int bsize)
      : pixelbuffer{}, lastoffset(0), maxoffset(bsize), current(0), buffer(_buffer)
    {
    }
    void read_page(); // will throw IO error if not enough space in buffer
    unsigned int nextpixel()
    {
        return current < 14 ? pixelbuffer[current++] : 0;
    }
};

#define wbuffer(i) ((unsigned short)buffer[lastoffset + 15 - i])

void pana_cs6_page_decoder::read_page()
{
    if (!buffer || (maxoffset - lastoffset < 16))
        ;
    pixelbuffer[0] = (wbuffer(0) << 6) | (wbuffer(1) >> 2); // 14 bit
    pixelbuffer[1] = (((wbuffer(1) & 0x3) << 12) | (wbuffer(2) << 4) | (wbuffer(3) >> 4)) & 0x3fff;
    pixelbuffer[2] = (wbuffer(3) >> 2) & 0x3;
    pixelbuffer[3] = ((wbuffer(3) & 0x3) << 8) | wbuffer(4);
    pixelbuffer[4] = (wbuffer(5) << 2) | (wbuffer(6) >> 6);
    pixelbuffer[5] = ((wbuffer(6) & 0x3f) << 4) | (wbuffer(7) >> 4);
    pixelbuffer[6] = (wbuffer(7) >> 2) & 0x3;
    pixelbuffer[7] = ((wbuffer(7) & 0x3) << 8) | wbuffer(8);
    pixelbuffer[8] = ((wbuffer(9) << 2) & 0x3fc) | (wbuffer(10) >> 6);
    pixelbuffer[9] = ((wbuffer(10) << 4) | (wbuffer(11) >> 4)) & 0x3ff;
    pixelbuffer[10] = (wbuffer(11) >> 2) & 0x3;
    pixelbuffer[11] = ((wbuffer(11) & 0x3) << 8) | wbuffer(12);
    pixelbuffer[12] = (((wbuffer(13) << 2) & 0x3fc) | wbuffer(14) >> 6) & 0x3ff;
    pixelbuffer[13] = ((wbuffer(14) << 4) | (wbuffer(15) >> 4)) & 0x3ff;
    current = 0;
    lastoffset += 16;
}
#undef wbuffer

void DCraw::panasonic_load_raw()
{
    int enc_blck_size = RT_pana_info.bpp == 12 ? 10 : 9;
    if (RT_pana_info.encoding == 5) {
        pana_bits_t pana_bits(ifp, load_flags, RT_pana_info.encoding);
        pana_bits(0, 0);
        unsigned bytes[16] = {};
        for (int row = 0; row < raw_height; ++row) {
            ushort* raw_block_data = raw_image + row * raw_width;

            for (int col = 0; col < raw_width; col += enc_blck_size) {
                pana_bits(0, bytes);

                if (RT_pana_info.bpp == 12) {
                    raw_block_data[col] = ((bytes[1] & 0xF) << 8) + bytes[0];
                    raw_block_data[col + 1] = 16 * bytes[2] + (bytes[1] >> 4);
                    raw_block_data[col + 2] = ((bytes[4] & 0xF) << 8) + bytes[3];
                    raw_block_data[col + 3] = 16 * bytes[5] + (bytes[4] >> 4);
                    raw_block_data[col + 4] = ((bytes[7] & 0xF) << 8) + bytes[6];
                    raw_block_data[col + 5] = 16 * bytes[8] + (bytes[7] >> 4);
                    raw_block_data[col + 6] = ((bytes[10] & 0xF) << 8) + bytes[9];
                    raw_block_data[col + 7] = 16 * bytes[11] + (bytes[10] >> 4);
                    raw_block_data[col + 8] = ((bytes[13] & 0xF) << 8) + bytes[12];
                    raw_block_data[col + 9] = 16 * bytes[14] + (bytes[13] >> 4);
                }
                else if (RT_pana_info.bpp == 14) {
                    raw_block_data[col] = bytes[0] + ((bytes[1] & 0x3F) << 8);
                    raw_block_data[col + 1] = (bytes[1] >> 6) + 4 * (bytes[2]) + ((bytes[3] & 0xF) << 10);
                    raw_block_data[col + 2] = (bytes[3] >> 4) + 16 * (bytes[4]) + ((bytes[5] & 3) << 12);
                    raw_block_data[col + 3] = ((bytes[5] & 0xFC) >> 2) + (bytes[6] << 6);
                    raw_block_data[col + 4] = bytes[7] + ((bytes[8] & 0x3F) << 8);
                    raw_block_data[col + 5] = (bytes[8] >> 6) + 4 * bytes[9] + ((bytes[10] & 0xF) << 10);
                    raw_block_data[col + 6] = (bytes[10] >> 4) + 16 * bytes[11] + ((bytes[12] & 3) << 12);
                    raw_block_data[col + 7] = ((bytes[12] & 0xFC) >> 2) + (bytes[13] << 6);
                    raw_block_data[col + 8] = bytes[14] + ((bytes[15] & 0x3F) << 8);
                }
            }
        }
    } else if (RT_pana_info.encoding == 6) {
        panasonicC6_load_raw();
    } else if (RT_pana_info.encoding == 7) {
        panasonicC7_load_raw();
    } else {
        pana_bits_t pana_bits(ifp, load_flags, RT_pana_info.encoding);
        pana_bits(0, 0);
        int sh = 0, pred[2], nonz[2];
        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < raw_width; ++col) {
                int i;
                if ((i = col % 14) == 0) {
                    pred[0] = pred[1] = nonz[0] = nonz[1] = 0;
                }
                if (i % 3 == 2) {
                    sh = 4 >> (3 - pana_bits(2));
                }
                if (nonz[i & 1]) {
                    int j;
                    if ((j = pana_bits(8))) {
                        if ((pred[i & 1] -= 0x80 << sh) < 0 || sh == 4) {
                            pred[i & 1] &= ~(-1 << sh);
                        }
                        pred[i & 1] += j << sh;
                    }
                } else if ((nonz[i & 1] = pana_bits(8)) || i > 11) {
                    pred[i & 1] = nonz[i & 1] << 4 | pana_bits(4);
                }
                if ((raw_image[(row)*raw_width+(col)] = pred[col & 1]) > 4098 && col < width) {
                    derror();
                }
            }
        }
    }
}

void DCraw::panasonicC6_load_raw()
{
    constexpr int rowstep = 16;
    const int blocksperrow = raw_width / 11;
    const int rowbytes = blocksperrow * 16;
    unsigned char *iobuf = (unsigned char *)malloc(rowbytes * rowstep);
    merror(iobuf, "panasonicC6_load_raw()");

    for (int row = 0; row < raw_height - rowstep + 1; row += rowstep) {
        const int rowstoread = MIN(rowstep, raw_height - row);
        fread(iobuf, rowbytes, rowstoread, ifp);
        pana_cs6_page_decoder page(iobuf, rowbytes * rowstoread);
        for (int crow = 0, col = 0; crow < rowstoread; ++crow, col = 0) {
            unsigned short *rowptr = &raw_image[(row + crow) * raw_width];
            for (int rblock = 0; rblock < blocksperrow; rblock++) {
                page.read_page();
                unsigned oddeven[2] = {0, 0}, nonzero[2] = {0, 0};
                unsigned pmul = 0, pixel_base = 0;
                for (int pix = 0; pix < 11; ++pix) {
                    if (pix % 3 == 2) {
                        unsigned base = page.nextpixel();
                        if (base > 3) {
                            derror();
                        }
                        if (base == 3) {
                            base = 4;
                        }
                        pixel_base = 0x200 << base;
                        pmul = 1 << base;
                    }
                    unsigned epixel = page.nextpixel();
                    if (oddeven[pix % 2]) {
                        epixel *= pmul;
                        if (pixel_base < 0x2000 && nonzero[pix % 2] > pixel_base) {
                            epixel += nonzero[pix % 2] - pixel_base;
                        }
                        nonzero[pix % 2] = epixel;
                    } else {
                        oddeven[pix % 2] = epixel;
                        if (epixel) {
                            nonzero[pix % 2] = epixel;
                        } else {
                            epixel = nonzero[pix % 2];
                        }
                    }
                    const unsigned spix = epixel - 0xf;
                    if (spix <= 0xffff) {
                        rowptr[col++] = spix & 0xffff;
                    } else {
                        epixel = (((signed int)(epixel + 0x7ffffff1)) >> 0x1f);
                        rowptr[col++] = epixel & 0x3fff;
                    }
                }
            }
        }
    }
    free(iobuf);
    tiff_bps = RT_pana_info.bpp;
}

void DCraw::panasonicC7_load_raw()
{
    constexpr int rowstep = 16;
    const int pixperblock = RT_pana_info.bpp == 14 ? 9 : 10;
    const int rowbytes = raw_width / pixperblock * 16;

    unsigned char *iobuf = (unsigned char *)malloc(rowbytes * rowstep);
    merror(iobuf, "panasonicC7_load_raw()");
    for (int row = 0; row < raw_height - rowstep + 1; row += rowstep) {
        const int rowstoread = MIN(rowstep, raw_height - row);
        fread (iobuf, rowbytes, rowstoread, ifp);
        unsigned char *bytes = iobuf;
        for (int crow = 0; crow < rowstoread; crow++) {
            ushort *rowptr = &raw_image[(row + crow) * raw_width];
            for (int col = 0; col < raw_width - pixperblock + 1; col += pixperblock, bytes += 16) {
                if (RT_pana_info.bpp == 14) {
                    rowptr[col] = bytes[0] + ((bytes[1] & 0x3F) << 8);
                    rowptr[col + 1] = (bytes[1] >> 6) + 4 * (bytes[2]) + ((bytes[3] & 0xF) << 10);
                    rowptr[col + 2] = (bytes[3] >> 4) + 16 * (bytes[4]) + ((bytes[5] & 3) << 12);
                    rowptr[col + 3] = ((bytes[5] & 0xFC) >> 2) + (bytes[6] << 6);
                    rowptr[col + 4] = bytes[7] + ((bytes[8] & 0x3F) << 8);
                    rowptr[col + 5] = (bytes[8] >> 6) + 4 * bytes[9] + ((bytes[10] & 0xF) << 10);
                    rowptr[col + 6] = (bytes[10] >> 4) + 16 * bytes[11] + ((bytes[12] & 3) << 12);
                    rowptr[col + 7] = ((bytes[12] & 0xFC) >> 2) + (bytes[13] << 6);
                    rowptr[col + 8] = bytes[14] + ((bytes[15] & 0x3F) << 8);
                } else if (RT_pana_info.bpp == 12) { // have not seen in the wild yet
                    rowptr[col] = ((bytes[1] & 0xF) << 8) + bytes[0];
                    rowptr[col + 1] = 16 * bytes[2] + (bytes[1] >> 4);
                    rowptr[col + 2] = ((bytes[4] & 0xF) << 8) + bytes[3];
                    rowptr[col + 3] = 16 * bytes[5] + (bytes[4] >> 4);
                    rowptr[col + 4] = ((bytes[7] & 0xF) << 8) + bytes[6];
                    rowptr[col + 5] = 16 * bytes[8] + (bytes[7] >> 4);
                    rowptr[col + 6] = ((bytes[10] & 0xF) << 8) + bytes[9];
                    rowptr[col + 7] = 16 * bytes[11] + (bytes[10] >> 4);
                    rowptr[col + 8] = ((bytes[13] & 0xF) << 8) + bytes[12];
                    rowptr[col + 9] = 16 * bytes[14] + (bytes[13] >> 4);
                }
            }
        }
    }
    free(iobuf);
    tiff_bps = RT_pana_info.bpp;
}
