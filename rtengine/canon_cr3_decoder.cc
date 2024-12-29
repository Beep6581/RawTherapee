/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 RawTherapee development team
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

// Code adapted from ART
// https://bitbucket.org/agriggio/art/
/*
 *
 *  This file is part of ART.
 *
 *  ART is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ART.  If not, see <https://www.gnu.org/licenses/>.
*/

// Code adapted from libraw
// https://github.com/LibRaw/LibRaw/
/*
 * File: libraw_crxdec.cpp
 * Copyright (C) 2018-2019 Alexey Danilchenko
 * Copyright (C) 2019 Alex Tutubalin, LibRaw LLC
 *
   Canon CR3 file decoder

LibRaw is free software; you can redistribute it and/or modify
it under the terms of the one of two licenses as you choose:

1. GNU LESSER GENERAL PUBLIC LICENSE version 2.1
   (See file LICENSE.LGPL provided in LibRaw distribution archive for details).

2. COMMON DEVELOPMENT AND DISTRIBUTION LICENSE (CDDL) Version 1.0
   (See file LICENSE.CDDL provided in LibRaw distribution archive for details).

 */

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <vector>

#include "dcraw.h"

#include "rt_math.h"

void DCraw::parse_canon_cr3()
{
    unsigned long long szAtomList = ifp->size;
    short nesting = -1;
    short nTrack = -1;
    short TrackType;
    char AtomNameStack[128];
    strncpy(make, "Canon", sizeof(make));

    const int err = parseCR3(0, szAtomList, nesting, AtomNameStack, nTrack, TrackType);

    if ((err == 0 || err == -14) && nTrack >= 0) { // no error, or too deep nesting
        selectCRXTrack(nTrack);
    }
}

void DCraw::selectCRXTrack(unsigned short maxTrack)
{
    std::int64_t bitcounts[CanonCR3Data::CRXTRACKS_MAXCOUNT] = {};
    std::int64_t maxbitcount = 0;
    std::uint32_t maxjpegbytes = 0;

    memset(bitcounts, 0, sizeof(bitcounts));

    for (unsigned int i = 0; i <= maxTrack && i < RT_canon_CR3_data.CRXTRACKS_MAXCOUNT; ++i) {
        CanonCR3Data::crx_data_header_t* const d = &RT_canon_CR3_data.crx_header[i];

        if (d->MediaType == 1) { // RAW
            bitcounts[i] = std::int64_t(d->nBits) * std::int64_t(d->f_width) * std::int64_t(d->f_height);
            if (bitcounts[i] > maxbitcount) {
                maxbitcount = bitcounts[i];
            }
        } else if (d->MediaType == 2) { // JPEG
            if (d->MediaSize > maxjpegbytes) {
                maxjpegbytes = d->MediaSize;
                thumb_offset = d->MediaOffset;
                thumb_length = d->MediaSize;
            }
        }
    }

    if (maxbitcount < 8) {
        return;
    }

    bool has_framei = false;
    unsigned int framei = 0;
    unsigned int framecnt = 0;

    for (unsigned int i = 0; i <= maxTrack && i < RT_canon_CR3_data.CRXTRACKS_MAXCOUNT; ++i) {
        if (bitcounts[i] == maxbitcount) {
            if (framecnt <= shot_select) {
                has_framei = true;
                framei = i;
            }
            framecnt++;
        }
    }

    is_raw = framecnt;

    if (has_framei) {
        CanonCR3Data::crx_data_header_t* const d = &RT_canon_CR3_data.crx_header[framei];
        data_offset = d->MediaOffset;
        // data_size = d->MediaSize;
        raw_width = d->f_width;
        raw_height = d->f_height;
        load_raw = &DCraw::crxLoadRaw;

        switch (d->cfaLayout) {
            case 0: {
                filters = 0x94949494;
                break;
            }
            case 1: {
                filters = 0x61616161;
                break;
            }
            case 2: {
                filters = 0x49494949;
                break;
            }
            case 3: {
                filters = 0x16161616;
                break;
            }
        }

        RT_canon_CR3_data.crx_track_selected = framei;

        int tiff_idx = -1;
        std::int64_t tpixels = 0;

        for (unsigned int i = 0; i < tiff_nifds; ++i) {
            if (std::int64_t(tiff_ifd[i].height) * std::int64_t(tiff_ifd[i].height) > tpixels) {
                tpixels = std::int64_t(tiff_ifd[i].height) * std::int64_t(tiff_ifd[i].height);
                tiff_idx = i;
            }
        }

        if (tiff_idx >= 0) {
            flip = tiff_ifd[tiff_idx].flip;
        }
    }
}

int DCraw::parseCR3(
    unsigned long long oAtomList,
    unsigned long long szAtomList,
    short& nesting,
    char* AtomNameStack,
    short& nTrack,
    short& TrackType
)
{
    /*
    Atom starts with 4 bytes for Atom size and 4 bytes containing Atom name
    Atom size includes the length of the header and the size of all "contained"
    Atoms if Atom size == 1, Atom has the extended size stored in 8 bytes located
    after the Atom name if Atom size == 0, it is the last top-level Atom extending
    to the end of the file Atom name is often a 4 symbol mnemonic, but can be a
    4-byte integer
    */
    const char UIID_Canon[17] =
        "\x85\xc0\xb6\x87\x82\x0f\x11\xe0\x81\x11\xf4\xce\x46\x2b\x6a\x48";
//  const char UIID_Preview[17] =
//      "\xea\xf4\x2b\x5e\x1c\x98\x4b\x88\xb9\xfb\xb7\xdc\x40\x6e\x4d\x16";

    /*
    AtomType = 0 - unknown: "unk."
    AtomType = 1 - container atom: "cont"
    AtomType = 2 - leaf atom: "leaf"
    AtomType = 3 - can be container, can be leaf: "both"
    */
//  const char sAtomeType[4][5] = {"unk.", "cont", "leaf", "both"};
    short AtomType;
    static const struct {
        char AtomName[5];
        short AtomType;
    } AtomNamesList[] = {
        {"dinf", 1},
        {"edts", 1},
        {"fiin", 1},
        {"ipro", 1},
        {"iprp", 1},
        {"mdia", 1},
        {"meco", 1},
        {"mere", 1},
        {"mfra", 1},
        {"minf", 1},
        {"moof", 1},
        {"moov", 1},
        {"mvex", 1},
        {"paen", 1},
        {"schi", 1},
        {"sinf", 1},
        {"skip", 1},
        {"stbl", 1},
        {"stsd", 1},
        {"strk", 1},
        {"tapt", 1},
        {"traf", 1},
        {"trak", 1},

        {"cdsc", 2},
        {"colr", 2},
        {"dimg", 2},
        // {"dref", 2},
        {"free", 2},
        {"frma", 2},
        {"ftyp", 2},
        {"hdlr", 2},
        {"hvcC", 2},
        {"iinf", 2},
        {"iloc", 2},
        {"infe", 2},
        {"ipco", 2},
        {"ipma", 2},
        {"iref", 2},
        {"irot", 2},
        {"ispe", 2},
        {"meta", 2},
        {"mvhd", 2},
        {"pitm", 2},
        {"pixi", 2},
        {"schm", 2},
        {"thmb", 2},
        {"tkhd", 2},
        {"url ", 2},
        {"urn ", 2},

        {"CCTP", 1},
        {"CRAW", 1},

        {"JPEG", 2},
        {"CDI1", 2},
        {"CMP1", 2},

        {"CNCV", 2},
        {"CCDT", 2},
        {"CTBO", 2},
        {"CMT1", 2},
        {"CMT2", 2},
        {"CMT3", 2},
        {"CMT4", 2},
        {"THMB", 2},
        {"co64", 2},
        {"mdat", 2},
        {"mdhd", 2},
        {"nmhd", 2},
        {"stsc", 2},
        {"stsz", 2},
        {"stts", 2},
        {"vmhd", 2},

        {"dref", 3},
        {"uuid", 3},
    };

    const char sHandlerType[5][5] = {
        "unk.",
        "soun",
        "vide",
        "hint",
        "meta"
    };

    int err = 0;

    unsigned short tL; // Atom length represented in 4 or 8 bytes
    char nmAtom[5]; // Atom name
    unsigned long long oAtom;
    unsigned long long szAtom; // Atom offset and Atom size
    unsigned long long oAtomContent;
    unsigned long long szAtomContent; // offset and size of Atom content
    unsigned long long lHdr;

    char UIID[16];
    uchar CMP1[85];
    char HandlerType[5];
    char MediaFormatID[5];
//  unsigned int ImageWidth, ImageHeight;
    unsigned long relpos_inDir;
    unsigned long relpos_inBox;
    unsigned int szItem;
    unsigned int Tag;
    unsigned int lTag;
    unsigned short tItem;

    nmAtom[0] = MediaFormatID[0] = nmAtom[4] = MediaFormatID[4] = '\0';
    strncpy(HandlerType, sHandlerType[0], sizeof(HandlerType));
//  ImageWidth = ImageHeight = 0U;
    oAtom = oAtomList;
    ++nesting;

    if (nesting > 31) {
        return -14; // too deep nesting
    }

    short s_order = order;

    const auto is_bad_header =
        [this]() -> bool
        {
            return
                (
                    order != 0x4D4D
                    && order != 0x4949
                )
                || get2() != 0x002A
                || get4() != 0x00000008;
        };

    while ((oAtom + 8) <= (oAtomList + szAtomList)) {
        lHdr = 0U;
        err = 0;
        order = 0x4D4D;
        fseek(ifp, oAtom, SEEK_SET);
        szAtom = get4();
        for (unsigned int c = 0; c < 4; ++c) {
            nmAtom[c] = AtomNameStack[nesting * 4 + c] = fgetc(ifp);
        }
        AtomNameStack[(nesting + 1) * 4] = '\0';
        tL = 4;
        AtomType = 0;

        for (const auto& atom : AtomNamesList) {
            if (!strcmp(nmAtom, atom.AtomName)) {
                AtomType = atom.AtomType;
                break;
            }
        }

        if (!AtomType) {
            err = 1;
        }

        if (szAtom == 0) {
            if (nesting != 0) {
                err = -2;
                goto fin;
            }

            szAtom = szAtomList - oAtom;
            oAtomContent = oAtom + 8;
            szAtomContent = szAtom - 8;
        } else if (szAtom == 1) {
            if ((oAtom + 16) > (oAtomList + szAtomList)) {
                err = -3;
                goto fin;
            }

            tL = 8;
            szAtom = (static_cast<unsigned long long>(get4()) << 32) | get4();
            oAtomContent = oAtom + 16;
            szAtomContent = szAtom - 16;
        } else {
            oAtomContent = oAtom + 8;
            szAtomContent = szAtom - 8;
        }

        if (!strcmp(nmAtom, "trak")) {
            nTrack++;
            TrackType = 0;

            if (nTrack >= RT_canon_CR3_data.CRXTRACKS_MAXCOUNT) {
                break;
            }
        }

        if (!strcmp(AtomNameStack, "moovuuid")) {
            lHdr = 16;
            fread(UIID, 1, lHdr, ifp);

            if (!strncmp(UIID, UIID_Canon, lHdr)) {
                AtomType = 1;
            } else {
                fseek(ifp, -lHdr, SEEK_CUR);
            }
        } else if (!strcmp(AtomNameStack, "moovuuidCCTP")) {
            lHdr = 12;
        } else if (!strcmp(AtomNameStack, "moovuuidCMT1")) {
            const short q_order = order;
            order = get2();

            if (tL != 4 || is_bad_header()) {
                err = -4;
                goto fin;
            }

            parse_tiff_ifd(oAtomContent);
            order = q_order;
        } else if (!strcmp(AtomNameStack, "moovuuidCMT2")) {
            const short q_order = order;
            order = get2();

            if (tL != 4 || is_bad_header()) {
                err = -5;
                goto fin;
            }

            parse_exif(oAtomContent);
            order = q_order;
        } else if (!strcmp(AtomNameStack, "moovuuidCMT3")) {
            const short q_order = order;
            order = get2();

            if (tL != 4 || is_bad_header()) {
                err = -6;
                goto fin;
            }

            fseek(ifp, -12L, SEEK_CUR);
            parse_makernote(oAtomContent, 0);
            order = q_order;
        } else if (!strcmp(AtomNameStack, "moovuuidCMT4")) {
            const short q_order = order;
            order = get2();

            if (tL != 4 || is_bad_header()) {
                err = -6;
                goto fin;
            }

            const long off = ftell(ifp);
            parse_gps(oAtomContent);
            fseek(ifp, off, SEEK_SET);
//      parse_gps_libraw(oAtomContent);
            order = q_order;
        } else if (!strcmp(AtomNameStack, "moovtrakmdiahdlr")) {
            fseek(ifp, 8, SEEK_CUR);
            for (unsigned int c = 0; c < 4; ++c) {
                HandlerType[c] = fgetc(ifp);
            }

            for (unsigned int c = 1; c < sizeof(sHandlerType) / sizeof(*sHandlerType); ++c) {
                if (!strcmp(HandlerType, sHandlerType[c])) {
                    TrackType = c;
                    break;
                }
            }
        } else if (!strcmp(AtomNameStack, "moovtrakmdiaminfstblstsd")) {
            if (szAtomContent >= 16) {
                fseek(ifp, 12, SEEK_CUR);
                lHdr = 8;
            } else {
                err = -7;
                goto fin;
            }

            for (unsigned int c = 0; c < 4; ++c) {
                MediaFormatID[c] = fgetc(ifp);
            }

            if (TrackType == 2 && !strcmp(MediaFormatID, "CRAW")) {
                if (szAtomContent >= 44) {
                    fseek(ifp, 24, SEEK_CUR);
                } else {
                    err = -8;
                    goto fin;
                }
            } else {
                AtomType = 2; // only continue for CRAW
                lHdr = 0;
            }

            /*ImageWidth = */ get2();
            /*ImageHeight = */ get2();
        } else if (!strcmp(AtomNameStack, "moovtrakmdiaminfstblstsdCRAW")) {
            lHdr = 82;
        } else if (!strcmp(AtomNameStack, "moovtrakmdiaminfstblstsdCRAWCMP1")) {
            int read_size = szAtomContent > 85 ? 85 : szAtomContent;
            if (szAtomContent >= 40) {
                fread(CMP1, 1, read_size, ifp);
            } else {
                err = -7;
                goto fin;
            }

            if (crxParseImageHeader(CMP1, nTrack, read_size)) {
                RT_canon_CR3_data.crx_header[nTrack].MediaType = 1;
            }
        } else if (!strcmp(AtomNameStack, "moovtrakmdiaminfstblstsdCRAWJPEG")) {
            RT_canon_CR3_data.crx_header[nTrack].MediaType = 2;
        } else if (!strcmp(AtomNameStack, "moovtrakmdiaminfstblstsz")) {
            if (szAtomContent == 12) {
                fseek(ifp, 4, SEEK_CUR);
            } else if (szAtomContent == 16) {
                fseek(ifp, 12, SEEK_CUR);
            } else {
                err = -9;
                goto fin;
            }

            RT_canon_CR3_data.crx_header[nTrack].MediaSize = get4();
        } else if (!strcmp(AtomNameStack, "moovtrakmdiaminfstblco64")) {
            if (szAtomContent == 16) {
                fseek(ifp, 8, SEEK_CUR);
            } else {
                err = -10;
                goto fin;
            }

            RT_canon_CR3_data.crx_header[nTrack].MediaOffset = (static_cast<unsigned long long>(get4()) << 32) | get4();
        }

        if (
            nTrack >= 0 && nTrack < RT_canon_CR3_data.CRXTRACKS_MAXCOUNT &&
            RT_canon_CR3_data.crx_header[nTrack].MediaSize
            && RT_canon_CR3_data.crx_header[nTrack].MediaOffset
            && oAtom + szAtom >= oAtomList + szAtomList
            && !strncmp(AtomNameStack, "moovtrakmdiaminfstbl", 20)
        ) {
            if (TrackType == 4 && !strcmp(MediaFormatID, "CTMD")) {
                order = 0x4949;
                relpos_inDir = 0;

                while (relpos_inDir + 6 < RT_canon_CR3_data.crx_header[nTrack].MediaSize) {
                    fseek(ifp, RT_canon_CR3_data.crx_header[nTrack].MediaOffset + relpos_inDir, SEEK_SET);
                    szItem = get4();
                    tItem = get2();

                    if ((relpos_inDir + szItem) > RT_canon_CR3_data.crx_header[nTrack].MediaSize) {
                        err = -11;
                        goto fin;
                    }

                    if (
                        tItem == 7
                        || tItem == 8
                        || tItem == 9
                    ) {
                        relpos_inBox = relpos_inDir + 12;

                        while (relpos_inBox + 8 < relpos_inDir + szItem) {
                            fseek(ifp, RT_canon_CR3_data.crx_header[nTrack].MediaOffset + relpos_inBox, SEEK_SET);
                            lTag = get4();
                            Tag = get4();

                            if (lTag < 8) {
                                err = -12;
                                goto fin;
                            } else if (relpos_inBox + lTag > relpos_inDir + szItem) {
                                err = -11;
                                goto fin;
                            }

                            if (
                                Tag == 0x927C
                                && (
                                    tItem == 7
                                    || tItem == 8
                                )
                            ) {
                                fseek(ifp, RT_canon_CR3_data.crx_header[nTrack].MediaOffset + relpos_inBox + 8, SEEK_SET);
                                const short q_order = order;
                                order = get2();

                                if (is_bad_header()) {
                                    err = -13;
                                    goto fin;
                                }

                                fseek(ifp, -8, SEEK_CUR);
                                RT_canon_CR3_data.CR3_CTMDtag = 1;
                                parse_makernote(RT_canon_CR3_data.crx_header[nTrack].MediaOffset + relpos_inBox + 8, 0);
                                RT_canon_CR3_data.CR3_CTMDtag = 0;
                                order = q_order;
                            }

                            relpos_inBox += lTag;
                        }
                    }
                    if (!szItem) {
                        goto fin;
                    }
                    relpos_inDir += szItem;
                }

                order = 0x4D4D;
            }
        }

        if (AtomType == 1) {
            err = parseCR3(oAtomContent + lHdr, szAtomContent - lHdr, nesting, AtomNameStack, nTrack, TrackType);

            if (err) {
                goto fin;
            }
        }

        oAtom += szAtom;
    }

fin:
    --nesting;

    if (nesting >= 0) {
        AtomNameStack[nesting * 4] = '\0';
    }

    order = s_order;
    return err;
}

// -----------------------------------------------------------------------------

namespace
{

unsigned int sgetn(int n, unsigned char* s)
{
    unsigned int result = 0;

    while (n-- > 0) {
        result = (result << 8) | (*s++);
    }

    return result;
}

// this should be divisible by 4
constexpr std::uint64_t CRX_BUF_SIZE = 0x10000;

#if !defined (_WIN32) || (defined (__GNUC__) && !defined (__INTRINSIC_SPECIAL__BitScanReverse))
/* __INTRINSIC_SPECIAL__BitScanReverse found in MinGW32-W64 v7.30 headers, may be there is a better solution? */
inline void _BitScanReverse(std::uint32_t* Index, unsigned long Mask)
{
    *Index = sizeof(unsigned long) * 8 - 1 - __builtin_clzl(Mask);
}
std::uint32_t _byteswap_ulong(std::uint32_t x)
{
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return x;
#else
    return __builtin_bswap32(x);
#endif
}
#endif

struct LibRaw_abstract_datastream {
    rtengine::IMFILE* ifp;

    void lock()
    {
    }
    void unlock()
    {
    }
    void seek(long p, int how)
    {
        fseek(ifp, p, how);
    }
    int read(void* dst, int es, int count)
    {
        return fread(dst, es, count, ifp);
    }
};

struct CrxBitstream {
    std::uint8_t mdatBuf[CRX_BUF_SIZE];
    std::uint64_t mdatSize;
    std::uint64_t curBufOffset;
    std::uint32_t curPos;
    std::uint32_t curBufSize;
    std::uint32_t bitData;
    std::int32_t bitsLeft;
    LibRaw_abstract_datastream* input;
};

struct CrxBandParam {
    CrxBitstream bitStream;
    std::int16_t subbandWidth;
    std::int16_t subbandHeight;
    std::int32_t roundedBitsMask;
    std::int32_t roundedBits;
    std::int16_t curLine;
    std::int32_t* lineBuf0;
    std::int32_t* lineBuf1;
    std::int32_t* lineBuf2;
    std::int32_t sParam;
    std::int32_t kParam;
    std::int32_t* paramData;
    std::int32_t* nonProgrData;
    bool supportsPartial;
};

struct CrxWaveletTransform {
    std::int32_t* subband0Buf;
    std::int32_t* subband1Buf;
    std::int32_t* subband2Buf;
    std::int32_t* subband3Buf;
    std::int32_t* lineBuf[8];
    std::int16_t curLine;
    std::int16_t curH;
    std::int8_t fltTapH;
    std::int16_t height;
    std::int16_t width;
};

struct CrxSubband {
    CrxBandParam* bandParam;
    std::uint64_t mdatOffset;
    std::uint8_t* bandBuf;
    std::uint16_t width;
    std::uint16_t height;
    std::int32_t qParam;
    std::int32_t kParam;
    std::int32_t qStepBase;
    std::uint32_t qStepMult;
    bool supportsPartial;
    std::int32_t bandSize;
    std::uint64_t dataSize;
    std::int64_t dataOffset;
    short rowStartAddOn;
    short rowEndAddOn;
    short colStartAddOn;
    short colEndAddOn;
    short levelShift;
};

struct CrxPlaneComp {
    std::uint8_t* compBuf;
    CrxSubband* subBands;
    CrxWaveletTransform* waveletTransform;
    std::int8_t compNumber;
    std::int64_t dataOffset;
    std::int32_t compSize;
    bool supportsPartial;
    std::int32_t roundedBitsMask;
    std::int8_t tileFlag;
};

struct CrxQStep {
    std::uint32_t *qStepTbl;
    int width;
    int height;
};

struct CrxTile {
    CrxPlaneComp* comps;
    std::int8_t tileFlag;
    std::int8_t tileNumber;
    std::int64_t dataOffset;
    std::int32_t tileSize;
    std::uint16_t width;
    std::uint16_t height;
    bool hasQPData;
    CrxQStep *qStep;
    std::uint32_t mdatQPDataSize;
    std::uint16_t mdatExtraSize;
};

struct CrxImage {
    std::uint8_t nPlanes;
    std::uint16_t planeWidth;
    std::uint16_t planeHeight;
    std::uint8_t samplePrecision;
    std::uint8_t medianBits;
    std::uint8_t subbandCount;
    std::uint8_t levels;
    std::uint8_t nBits;
    std::uint8_t encType;
    std::uint8_t tileCols;
    std::uint8_t tileRows;
    CrxTile* tiles;
    std::uint64_t mdatOffset;
    std::uint64_t mdatSize;
    std::int16_t* outBufs[4]; // one per plane
    std::int16_t* planeBuf;
    LibRaw_abstract_datastream* input;
};

enum TileFlags {
    E_HAS_TILES_ON_THE_RIGHT = 1,
    E_HAS_TILES_ON_THE_LEFT = 2,
    E_HAS_TILES_ON_THE_BOTTOM = 4,
    E_HAS_TILES_ON_THE_TOP = 8
};

const std::int32_t exCoefNumTbl[144] = {
1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0,
0, 0, 1, 2, 2, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 2, 2,
1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2, 1, 1, 1,
1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

constexpr std::int32_t q_step_tbl[6] = {0x28, 0x2D, 0x33, 0x39, 0x40, 0x48};

const std::uint32_t JS[32] = {
    0x0001, 0x0001, 0x0001, 0x0001, 0x0002, 0x0002, 0x0002, 0x0002,
    0x0004, 0x0004, 0x0004, 0x0004, 0x0008, 0x0008, 0x0008, 0x0008,
    0x0010, 0x0010, 0x0020, 0x0020, 0x0040, 0x0040, 0x0080, 0x0080,
    0x0100, 0x0200, 0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000
};

const std::uint32_t J[32] = {
    0x0, 0x0, 0x0, 0x0, 0x1, 0x1, 0x1, 0x1, 0x2, 0x2, 0x2,
    0x2, 0x3, 0x3, 0x3, 0x3, 0x4, 0x4, 0x5, 0x5, 0x6, 0x6,
    0x7, 0x7, 0x8, 0x9, 0xA, 0xB, 0xC, 0xD, 0xE, 0xF
};

inline void crxFillBuffer(CrxBitstream* bitStrm)
{
    if (bitStrm->curPos >= bitStrm->curBufSize && bitStrm->mdatSize) {
        bitStrm->curPos = 0;
        bitStrm->curBufOffset += bitStrm->curBufSize;
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
#ifndef _OPENMP
            bitStrm->input->lock();
#endif
            bitStrm->input->seek(bitStrm->curBufOffset, SEEK_SET);
            bitStrm->curBufSize = bitStrm->input->read(bitStrm->mdatBuf, 1, std::min(bitStrm->mdatSize, CRX_BUF_SIZE));
#ifndef _OPENMP
            bitStrm->input->unlock();
#endif

            if (bitStrm->curBufSize < 1) {  // nothing read
                throw std::runtime_error("Unexpected end of file in CRX bitstream");
            }

            bitStrm->mdatSize -= bitStrm->curBufSize;
        }
    }
}

inline int crxBitstreamGetZeros(CrxBitstream* bitStrm)
{
    std::uint32_t nonZeroBit = 0;
    std::uint64_t nextData = 0;
    std::int32_t result = 0;

    if (bitStrm->bitData) {
        _BitScanReverse(&nonZeroBit, bitStrm->bitData);
        result = 31 - nonZeroBit;
        bitStrm->bitData <<= 32 - nonZeroBit;
        bitStrm->bitsLeft -= 32 - nonZeroBit;
    } else {
        std::uint32_t bitsLeft = bitStrm->bitsLeft;

        while (true) {
            while (bitStrm->curPos + 4 <= bitStrm->curBufSize) {
                nextData = _byteswap_ulong(*reinterpret_cast<std::uint32_t*>(bitStrm->mdatBuf + bitStrm->curPos));
                bitStrm->curPos += 4;
                crxFillBuffer(bitStrm);

                if (nextData) {
                    _BitScanReverse(&nonZeroBit, static_cast<std::uint32_t>(nextData));
                    result = bitsLeft + 31 - nonZeroBit;
                    bitStrm->bitData = nextData << (32 - nonZeroBit);
                    bitStrm->bitsLeft = nonZeroBit;
                    return result;
                }

                bitsLeft += 32;
            }

            if (bitStrm->curBufSize < bitStrm->curPos + 1) {
                break; // error
            }

            nextData = bitStrm->mdatBuf[bitStrm->curPos++];
            crxFillBuffer(bitStrm);

            if (nextData) {
                break;
            }

            bitsLeft += 8;
        }

        _BitScanReverse(&nonZeroBit, static_cast<std::uint32_t>(nextData));
        result = static_cast<std::uint32_t>(bitsLeft + 7 - nonZeroBit);
        bitStrm->bitData = nextData << (32 - nonZeroBit);
        bitStrm->bitsLeft = nonZeroBit;
    }

    return result;
}

inline std::uint32_t crxBitstreamGetBits(CrxBitstream* bitStrm, int bits)
{
    int bitsLeft = bitStrm->bitsLeft;
    std::uint32_t bitData = bitStrm->bitData;
    std::uint32_t nextWord;
    std::uint8_t nextByte;
    std::uint32_t result;

    if (bitsLeft < bits) {
        // get them from stream
        if (bitStrm->curPos + 4 <= bitStrm->curBufSize) {
            nextWord = _byteswap_ulong(*reinterpret_cast<std::uint32_t*>(bitStrm->mdatBuf + bitStrm->curPos));
            bitStrm->curPos += 4;
            crxFillBuffer(bitStrm);
            bitStrm->bitsLeft = 32 - (bits - bitsLeft);
            result = ((nextWord >> bitsLeft) | bitData) >> (32 - bits);
            bitStrm->bitData = nextWord << (bits - bitsLeft);
            return result;
        }

        // less than a word left - read byte at a time
        do {
            if (bitStrm->curPos >= bitStrm->curBufSize) {
                break; // error
            }

            bitsLeft += 8;
            nextByte = bitStrm->mdatBuf[bitStrm->curPos++];
            crxFillBuffer(bitStrm);
            bitData |= nextByte << (32 - bitsLeft);
        } while (bitsLeft < bits);
    }

    result = bitData >> (32 - bits); // 32-bits
    bitStrm->bitData = bitData << bits;
    bitStrm->bitsLeft = bitsLeft - bits;

    return result;
}

inline std::int32_t crxPrediction(std::int32_t left, std::int32_t top, std::int32_t deltaH, std::int32_t deltaV) 
{
    std::int32_t symb[4] = {left + deltaH, left + deltaH, left, top};
    return symb[(((deltaV < 0) ^ (deltaH < 0)) << 1) + ((left < top) ^ (deltaH < 0))];
}

inline std::int32_t crxPredictKParameter(std::int32_t prevK, std::int32_t bitCode, std::int32_t maxVal = 0)
{
    const std::int32_t newKParam =
        prevK
        - (bitCode < (1 << prevK >> 1))
        + ((bitCode >> prevK) > 2) + ((bitCode >> prevK) > 5);

    return
        !maxVal || newKParam < maxVal
            ? newKParam
            : maxVal;
}

inline void crxDecodeSymbolL1(CrxBandParam* param, bool doMedianPrediction, bool notEOL = false)
{
    if (doMedianPrediction) {
        const std::int32_t delta = param->lineBuf0[1] - param->lineBuf0[0];
        const std::int32_t symb[4] = {
            delta + param->lineBuf1[0],
            delta + param->lineBuf1[0],
            param->lineBuf1[0],
            param->lineBuf0[1]
        };

        param->lineBuf1[1] = symb[
            (((param->lineBuf0[0] < param->lineBuf1[0]) ^ (delta < 0)) << 1)
            + ((param->lineBuf1[0] < param->lineBuf0[1]) ^ (delta < 0))
        ];
    } else {
        param->lineBuf1[1] = param->lineBuf0[1];
    }

    // get next error symbol
    std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

    if (bitCode >= 41) {
        bitCode = crxBitstreamGetBits(&param->bitStream, 21);
    } else if (param->kParam) {
        bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
    }

    // add converted (+/-) error code to predicted value
    param->lineBuf1[1] += -(bitCode & 1) ^ (bitCode >> 1);

    // for not end of the line - use one symbol ahead to estimate next K
    if (notEOL) {
        const std::int32_t nextDelta = (param->lineBuf0[2] - param->lineBuf0[1]) << 1;
        bitCode = (bitCode + std::abs(nextDelta)) >> 1;
        ++param->lineBuf0;
    }

    // update K parameter
    param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);

    ++param->lineBuf1;
}

bool crxDecodeLine(CrxBandParam* param)
{
    int length = param->subbandWidth;

    param->lineBuf1[0] = param->lineBuf0[1];

    for (; length > 1; --length) {
        if (param->lineBuf1[0] != param->lineBuf0[1] || param->lineBuf1[0] != param->lineBuf0[2]) {
            crxDecodeSymbolL1(param, true, true);
        } else {
            if (crxBitstreamGetBits(&param->bitStream, 1)) {
                int nSyms = 1;

                while (crxBitstreamGetBits(&param->bitStream, 1)) {
                    nSyms += JS[param->sParam];

                    if (nSyms > length) {
                        nSyms = length;
                        break;
                    }

                    if (param->sParam < 31) {
                        ++param->sParam;
                    }

                    if (nSyms == length) {
                        break;
                    }
                }

                if (nSyms < length) {
                    if (J[param->sParam]) {
                        nSyms += crxBitstreamGetBits(&param->bitStream, J[param->sParam]);
                    }

                    if (param->sParam > 0) {
                        --param->sParam;
                    }

                    if (nSyms > length) {
                        return false;
                    }
                }

                length -= nSyms;

                // copy symbol nSyms times
                param->lineBuf0 += nSyms;

                // copy symbol nSyms times
                while (nSyms-- > 0) {
                    param->lineBuf1[1] = param->lineBuf1[0];
                    ++param->lineBuf1;
                }
            }

            if (length > 0) {
                crxDecodeSymbolL1(param, false, length > 1);
            }
        }
    }

    if (length == 1) {
        crxDecodeSymbolL1(param, true, false);
    }

    param->lineBuf1[1] = param->lineBuf1[0] + 1;

    return true;
}

inline void crxDecodeSymbolL1Rounded(CrxBandParam* param, bool doSym = true, bool doCode = true)
{
    std::int32_t sym = param->lineBuf0[1];

    if (doSym) {
        // calculate the next symbol gradient
        const std::int32_t deltaH = param->lineBuf0[1] - param->lineBuf0[0];
        const std::int32_t symb[4] = {
            deltaH + param->lineBuf1[0],
            deltaH + param->lineBuf1[0],
            param->lineBuf1[0],
            param->lineBuf0[1]
        };
        sym = symb[
            (((param->lineBuf0[0] < param->lineBuf1[0]) ^ (deltaH < 0)) << 1)
            + ((param->lineBuf1[0] < param->lineBuf0[1]) ^ (deltaH < 0))
        ];
    }

    std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

    if (bitCode >= 41) {
        bitCode = crxBitstreamGetBits(&param->bitStream, 21);
    } else if (param->kParam) {
        bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
    }

    std::int32_t code = -(bitCode & 1) ^ (bitCode >> 1);
    param->lineBuf1[1] = param->roundedBitsMask * 2 * code + (code >> 31) + sym;

    if (doCode) {
        if (param->lineBuf0[2] > param->lineBuf0[1]) {
            code = (param->lineBuf0[2] - param->lineBuf0[1] + param->roundedBitsMask - 1) >> param->roundedBits;
        } else {
            code = -((param->lineBuf0[1] - param->lineBuf0[2] + param->roundedBitsMask) >> param->roundedBits);
        }

        param->kParam = crxPredictKParameter(param->kParam, (bitCode + 2 * std::abs(code)) >> 1, 15);
    } else {
        param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
    }

    ++param->lineBuf1;
}

bool crxDecodeLineRounded(CrxBandParam* param)
{
    bool valueReached = false;

    param->lineBuf0[0] = param->lineBuf0[1];
    param->lineBuf1[0] = param->lineBuf0[1];
    std::int32_t length = param->subbandWidth;

    for (; length > 1; --length) {
        if (std::abs(param->lineBuf0[2] - param->lineBuf0[1]) > param->roundedBitsMask) {
            crxDecodeSymbolL1Rounded(param);
            ++param->lineBuf0;
            valueReached = true;
        } else if (valueReached || std::abs(param->lineBuf0[0] - param->lineBuf1[0]) > param->roundedBitsMask) {
            crxDecodeSymbolL1Rounded(param);
            ++param->lineBuf0;
            valueReached = false;
        } else {
            int nSyms = 0;

            if (crxBitstreamGetBits(&param->bitStream, 1)) {
                nSyms = 1;

                while (crxBitstreamGetBits(&param->bitStream, 1)) {
                    nSyms += JS[param->sParam];

                    if (nSyms > length) {
                        nSyms = length;
                        break;
                    }

                    if (param->sParam < 31) {
                        ++param->sParam;
                    }

                    if (nSyms == length) {
                        break;
                    }
                }

                if (nSyms < length) {
                    if (J[param->sParam]) {
                        nSyms += crxBitstreamGetBits(&param->bitStream, J[param->sParam]);
                    }

                    if (param->sParam > 0) {
                        --param->sParam;
                    }
                }

                if (nSyms > length) {
                    return false;
                }
            }

            length -= nSyms;

            // copy symbol nSyms times
            param->lineBuf0 += nSyms;

            // copy symbol nSyms times
            while (nSyms-- > 0) {
                param->lineBuf1[1] = param->lineBuf1[0];
                ++param->lineBuf1;
            }

            if (length > 1) {
                crxDecodeSymbolL1Rounded(param, false);
                ++param->lineBuf0;
                valueReached = std::abs(param->lineBuf0[1] - param->lineBuf0[0]) > param->roundedBitsMask;
            } else if (length == 1) {
                crxDecodeSymbolL1Rounded(param, false, false);
            }
        }
    }

    if (length == 1) {
        crxDecodeSymbolL1Rounded(param, true, false);
    }

    param->lineBuf1[1] = param->lineBuf1[0] + 1;

    return true;
}

bool crxDecodeLineNoRefPrevLine(CrxBandParam* param)
{
    std::int32_t i = 0;

    for (; i < param->subbandWidth - 1; ++i) {
        if (param->lineBuf0[i + 2] || param->lineBuf0[i + 1] || param->lineBuf1[i]) {
            std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

            if (bitCode >= 41) {
                bitCode = crxBitstreamGetBits(&param->bitStream, 21);
            } else if (param->kParam) {
                bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
            }

            param->lineBuf1[i + 1] = -(bitCode & 1) ^ (bitCode >> 1);
            param->kParam = crxPredictKParameter(param->kParam, bitCode);

            if (param->lineBuf2[i + 1] - param->kParam <= 1) {
                if (param->kParam >= 15) {
                    param->kParam = 15;
                }
            } else {
                ++param->kParam;
            }
        } else {
            int nSyms = 0;

            if (crxBitstreamGetBits(&param->bitStream, 1)) {
                nSyms = 1;

                if (i != param->subbandWidth - 1) {
                    while (crxBitstreamGetBits(&param->bitStream, 1)) {
                        nSyms += JS[param->sParam];

                        if (i + nSyms > param->subbandWidth) {
                            nSyms = param->subbandWidth - i;
                            break;
                        }

                        if (param->sParam < 31) {
                            ++param->sParam;
                        }

                        if (i + nSyms == param->subbandWidth) {
                            break;
                        }
                    }

                    if (i + nSyms < param->subbandWidth) {
                        if (J[param->sParam]) {
                            nSyms += crxBitstreamGetBits(&param->bitStream, J[param->sParam]);
                        }

                        if (param->sParam > 0) {
                            --param->sParam;
                        }
                    }

                    if (i + nSyms > param->subbandWidth) {
                        return false;
                    }
                }
            } else if (i > param->subbandWidth) {
                return false;
            }

            if (nSyms > 0) {
                memset(param->lineBuf1 + i + 1, 0, nSyms * sizeof(std::int32_t));
                memset(param->lineBuf2 + i, 0, nSyms * sizeof(std::int32_t));
                i += nSyms;
            }

            if (i >= param->subbandWidth - 1) {
                if (i == param->subbandWidth - 1) {
                    std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

                    if (bitCode >= 41) {
                        bitCode = crxBitstreamGetBits(&param->bitStream, 21);
                    } else if (param->kParam) {
                        bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
                    }

                    param->lineBuf1[i + 1] = -((bitCode + 1) & 1) ^ ((bitCode + 1) >> 1);
                    param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
                    param->lineBuf2[i] = param->kParam;
                }

                continue;
            } else {
                std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

                if (bitCode >= 41) {
                    bitCode = crxBitstreamGetBits(&param->bitStream, 21);
                } else if (param->kParam) {
                    bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
                }

                param->lineBuf1[i + 1] = -((bitCode + 1) & 1) ^ ((bitCode + 1) >> 1);
                param->kParam = crxPredictKParameter(param->kParam, bitCode);

                if (param->lineBuf2[i + 1] - param->kParam <= 1) {
                    if (param->kParam >= 15) {
                        param->kParam = 15;
                    }
                } else {
                    ++param->kParam;
                }
            }
        }

        param->lineBuf2[i] = param->kParam;
    }

    if (i == param->subbandWidth - 1) {
        std::int32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

        if (bitCode >= 41) {
            bitCode = crxBitstreamGetBits(&param->bitStream, 21);
        } else if (param->kParam) {
            bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
        }

        param->lineBuf1[i + 1] = -(bitCode & 1) ^ (bitCode >> 1);
        param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
        param->lineBuf2[i] = param->kParam;
    }

    return true;
}

bool crxDecodeTopLine(CrxBandParam* param)
{
    param->lineBuf1[0] = 0;

    std::int32_t length = param->subbandWidth;

    // read the line from bitstream
    for (; length > 1; --length) {
        if (param->lineBuf1[0]) {
            param->lineBuf1[1] = param->lineBuf1[0];
        } else {
            if (crxBitstreamGetBits(&param->bitStream, 1)) {
                int nSyms = 1;

                while (crxBitstreamGetBits(&param->bitStream, 1)) {
                    nSyms += JS[param->sParam];

                    if (nSyms > length) {
                        nSyms = length;
                        break;
                    }

                    if (param->sParam < 31) {
                        ++param->sParam;
                    }

                    if (nSyms == length) {
                        break;
                    }
                }

                if (nSyms < length) {
                    if (J[param->sParam]) {
                        nSyms += crxBitstreamGetBits(&param->bitStream, J[param->sParam]);
                    }

                    if (param->sParam > 0) {
                        --param->sParam;
                    }

                    if (nSyms > length) {
                        return false;
                    }
                }

                length -= nSyms;

                // copy symbol nSyms times
                while (nSyms-- > 0) {
                    param->lineBuf1[1] = param->lineBuf1[0];
                    ++param->lineBuf1;
                }

                if (length <= 0) {
                    break;
                }
            }

            param->lineBuf1[1] = 0;
        }

        std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

        if (bitCode >= 41) {
            bitCode = crxBitstreamGetBits(&param->bitStream, 21);
        } else if (param->kParam) {
            bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
        }

        param->lineBuf1[1] += -(bitCode & 1) ^ (bitCode >> 1);
        param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
        ++param->lineBuf1;
    }

    if (length == 1) {
        param->lineBuf1[1] = param->lineBuf1[0];
        std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

        if (bitCode >= 41) {
            bitCode = crxBitstreamGetBits(&param->bitStream, 21);
        } else if (param->kParam) {
            bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
        }

        param->lineBuf1[1] += -(bitCode & 1) ^ (bitCode >> 1);
        param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
        ++param->lineBuf1;
    }

    param->lineBuf1[1] = param->lineBuf1[0] + 1;

    return true;
}

bool crxDecodeTopLineRounded(CrxBandParam* param)
{
    param->lineBuf1[0] = 0;

    std::int32_t length = param->subbandWidth;

    // read the line from bitstream
    for (; length > 1; --length) {
        if (std::abs(param->lineBuf1[0]) > param->roundedBitsMask) {
            param->lineBuf1[1] = param->lineBuf1[0];
        } else {
            int nSyms = 0;

            if (crxBitstreamGetBits(&param->bitStream, 1)) {
                nSyms = 1;

                while (crxBitstreamGetBits(&param->bitStream, 1)) {
                    nSyms += JS[param->sParam];

                    if (nSyms > length) {
                        nSyms = length;
                        break;
                    }

                    if (param->sParam < 31) {
                        ++param->sParam;
                    }

                    if (nSyms == length) {
                        break;
                    }
                }

                if (nSyms < length) {
                    if (J[param->sParam]) {
                        nSyms += crxBitstreamGetBits(&param->bitStream, J[param->sParam]);
                    }

                    if (param->sParam > 0) {
                        --param->sParam;
                    }

                    if (nSyms > length) {
                        return false;
                    }
                }
            }

            length -= nSyms;

            // copy symbol nSyms times
            while (nSyms-- > 0) {
                param->lineBuf1[1] = param->lineBuf1[0];
                ++param->lineBuf1;
            }

            if (length <= 0) {
                break;
            }

            param->lineBuf1[1] = 0;
        }

        std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

        if (bitCode >= 41) {
            bitCode = crxBitstreamGetBits(&param->bitStream, 21);
        } else if (param->kParam) {
            bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
        }

        const std::int32_t sVal = -(bitCode & 1) ^ (bitCode >> 1);
        param->lineBuf1[1] += param->roundedBitsMask * 2 * sVal + (sVal >> 31);
        param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
        ++param->lineBuf1;
    }

    if (length == 1) {
        std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

        if (bitCode >= 41) {
            bitCode = crxBitstreamGetBits(&param->bitStream, 21);
        } else if (param->kParam) {
            bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
        }

        const std::int32_t sVal = -(bitCode & 1) ^ (bitCode >> 1);
        param->lineBuf1[1] += param->roundedBitsMask * 2 * sVal + (sVal >> 31);
        param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
        ++param->lineBuf1;
    }

    param->lineBuf1[1] = param->lineBuf1[0] + 1;

    return true;
}

bool crxDecodeTopLineNoRefPrevLine(CrxBandParam* param)
{
    param->lineBuf0[0] = 0;
    param->lineBuf1[0] = 0;
    std::int32_t length = param->subbandWidth;

    for (; length > 1; --length) {
        if (param->lineBuf1[0]) {
            std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

            if (bitCode >= 41) {
                bitCode = crxBitstreamGetBits(&param->bitStream, 21);
            } else if (param->kParam) {
                bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
            }

            param->lineBuf1[1] = -(bitCode & 1) ^ (bitCode >> 1);
            param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
        } else {
            int nSyms = 0;

            if (crxBitstreamGetBits(&param->bitStream, 1)) {
                nSyms = 1;

                while (crxBitstreamGetBits(&param->bitStream, 1)) {
                    nSyms += JS[param->sParam];

                    if (nSyms > length) {
                        nSyms = length;
                        break;
                    }

                    if (param->sParam < 31) {
                        ++param->sParam;
                    }

                    if (nSyms == length) {
                        break;
                    }
                }

                if (nSyms < length) {
                    if (J[param->sParam]) {
                        nSyms += crxBitstreamGetBits(&param->bitStream, J[param->sParam]);
                    }

                    if (param->sParam > 0) {
                        --param->sParam;
                    }

                    if (nSyms > length) {
                        return false;
                    }
                }
            }

            length -= nSyms;

            // copy symbol nSyms times
            while (nSyms-- > 0) {
                param->lineBuf2[0] = 0;
                param->lineBuf1[1] = 0;
                ++param->lineBuf1;
                ++param->lineBuf2;
            }

            if (length <= 0) {
                break;
            }

            std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

            if (bitCode >= 41) {
                bitCode = crxBitstreamGetBits(&param->bitStream, 21);
            } else if (param->kParam) {
                bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
            }

            param->lineBuf1[1] = -((bitCode + 1) & 1) ^ ((bitCode + 1) >> 1);
            param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
        }

        param->lineBuf2[0] = param->kParam;
        ++param->lineBuf2;
        ++param->lineBuf1;
    }

    if (length == 1) {
        std::uint32_t bitCode = crxBitstreamGetZeros(&param->bitStream);

        if (bitCode >= 41) {
            bitCode = crxBitstreamGetBits(&param->bitStream, 21);
        } else if (param->kParam) {
            bitCode = crxBitstreamGetBits(&param->bitStream, param->kParam) | (bitCode << param->kParam);
        }

        param->lineBuf1[1] = -(bitCode & 1) ^ (bitCode >> 1);
        param->kParam = crxPredictKParameter(param->kParam, bitCode, 15);
        param->lineBuf2[0] = param->kParam;
        ++param->lineBuf1;
    }

    param->lineBuf1[1] = 0;

    return true;
}

bool crxDecodeLine(CrxBandParam* param, std::uint8_t* bandBuf)
{
    if (!param || !bandBuf) {
        return false;
    }
    if (param->curLine >= param->subbandHeight) {
        return false;
    }
    if (param->curLine == 0) {
        const std::int32_t lineLength = param->subbandWidth + 2;

        param->sParam = 0;
        param->kParam = 0;

        if (param->supportsPartial) {
            if (param->roundedBitsMask <= 0) {
                param->lineBuf0 = reinterpret_cast<std::int32_t*>(param->paramData);
                param->lineBuf1 = param->lineBuf0 + lineLength;
                const std::int32_t* const lineBuf = param->lineBuf1 + 1;

                if (!crxDecodeTopLine(param)) {
                    return false;
                }

                memcpy(bandBuf, lineBuf, param->subbandWidth * sizeof(std::int32_t));
                ++param->curLine;
            } else {
                param->roundedBits = 1;

                if (param->roundedBitsMask & ~1) {
                    while (param->roundedBitsMask >> param->roundedBits) {
                        ++param->roundedBits;
                    }
                }

                param->lineBuf0 = reinterpret_cast<std::int32_t*>(param->paramData);
                param->lineBuf1 = param->lineBuf0 + lineLength;
                const std::int32_t* const lineBuf = param->lineBuf1 + 1;

                if (!crxDecodeTopLineRounded(param)) {
                    return false;
                }

                memcpy(bandBuf, lineBuf, param->subbandWidth * sizeof(std::int32_t));
                ++param->curLine;
            }
        } else {
            param->lineBuf2 = reinterpret_cast<std::int32_t*>(param->nonProgrData);
            param->lineBuf0 = reinterpret_cast<std::int32_t*>(param->paramData);
            param->lineBuf1 = param->lineBuf0 + lineLength;
            const std::int32_t* const lineBuf = param->lineBuf1 + 1;

            if (!crxDecodeTopLineNoRefPrevLine(param)) {
                return false;
            }

            memcpy(bandBuf, lineBuf, param->subbandWidth * sizeof(std::int32_t));
            ++param->curLine;
        }
    } else if (!param->supportsPartial) {
        const std::int32_t lineLength = param->subbandWidth + 2;
        param->lineBuf2 = reinterpret_cast<std::int32_t*>(param->nonProgrData);

        if (param->curLine & 1) {
            param->lineBuf1 = reinterpret_cast<std::int32_t*>(param->paramData);
            param->lineBuf0 = param->lineBuf1 + lineLength;
        } else {
            param->lineBuf0 = reinterpret_cast<std::int32_t*>(param->paramData);
            param->lineBuf1 = param->lineBuf0 + lineLength;
        }

        const std::int32_t* const lineBuf = param->lineBuf1 + 1;

        if (!crxDecodeLineNoRefPrevLine(param)) {
            return false;
        }

        memcpy(bandBuf, lineBuf, param->subbandWidth * sizeof(std::int32_t));
        ++param->curLine;
    } else if (param->roundedBitsMask <= 0) {
        const std::int32_t lineLength = param->subbandWidth + 2;

        if (param->curLine & 1) {
            param->lineBuf1 = reinterpret_cast<std::int32_t*>(param->paramData);
            param->lineBuf0 = param->lineBuf1 + lineLength;
        } else {
            param->lineBuf0 = reinterpret_cast<std::int32_t*>(param->paramData);
            param->lineBuf1 = param->lineBuf0 + lineLength;
        }

        const std::int32_t* const lineBuf = param->lineBuf1 + 1;

        if (!crxDecodeLine(param)) {
            return false;
        }

        memcpy(bandBuf, lineBuf, param->subbandWidth * sizeof(std::int32_t));
        ++param->curLine;
    } else {
        const std::int32_t lineLength = param->subbandWidth + 2;

        if (param->curLine & 1) {
            param->lineBuf1 = reinterpret_cast<std::int32_t*>(param->paramData);
            param->lineBuf0 = param->lineBuf1 + lineLength;
        } else {
            param->lineBuf0 = reinterpret_cast<std::int32_t*>(param->paramData);
            param->lineBuf1 = param->lineBuf0 + lineLength;
        }

        const std::int32_t* const lineBuf = param->lineBuf1 + 1;

        if (!crxDecodeLineRounded(param)) {
            return false;
        }

        memcpy(bandBuf, lineBuf, param->subbandWidth * sizeof(std::int32_t));
        ++param->curLine;
    }

    return true;
}

inline int getSubbandRow(CrxSubband *band, int row)
{
    return row < band->rowStartAddOn
            ? 0
            : (row < band->height - band->rowEndAddOn ? row - band->rowEndAddOn
                        : band->height - band->rowEndAddOn - band->rowStartAddOn - 1);
}

bool crxDecodeLineWithIQuantization(CrxSubband* subband, CrxQStep *qStep)
{
    if (!subband->dataSize) {
        memset(subband->bandBuf, 0, subband->bandSize);
        return true;
    }

    if (subband->supportsPartial && !qStep) {
        std::uint32_t bitCode = crxBitstreamGetZeros(&subband->bandParam->bitStream);

        if (bitCode >= 23) {
            bitCode = crxBitstreamGetBits(&subband->bandParam->bitStream, 8);
        } else if (subband->kParam) {
            bitCode = crxBitstreamGetBits(&subband->bandParam->bitStream, subband->kParam) | (bitCode << subband->kParam);
        }

        subband->qParam += -(bitCode & 1) ^ (bitCode >> 1); // converting encoded to signed integer
        subband->kParam = crxPredictKParameter(subband->kParam, bitCode);

        if (subband->kParam > 7) {
            return false;
        }
    }

    if (!crxDecodeLine(subband->bandParam, subband->bandBuf)) {
        return false;
    }

    if (subband->width <= 0) {
        return true;
    }

    // update subband buffers
    std::int32_t* const bandBuf = reinterpret_cast<std::int32_t*>(subband->bandBuf);
    if (qStep) {
        // new version
        std::uint32_t *qStepTblPtr = &qStep->qStepTbl[qStep->width * getSubbandRow(subband, subband->bandParam->curLine - 1)];
        for (int i = 0; i < subband->colStartAddOn; ++i) {
            int32_t quantVal = subband->qStepBase + ((qStepTblPtr[0] * subband->qStepMult) >> 3);
            bandBuf[i] *= rtengine::LIM(quantVal, 1, 0x168000);
        }

        for (int i = subband->colStartAddOn; i < subband->width - subband->colEndAddOn; ++i) {
            int32_t quantVal = subband->qStepBase + ((qStepTblPtr[(i - subband->colStartAddOn) >> subband->levelShift] * subband->qStepMult) >> 3);
            bandBuf[i] *= rtengine::LIM(quantVal, 1, 0x168000);
        }
        int lastIdx = (subband->width - subband->colEndAddOn - subband->colStartAddOn - 1) >> subband->levelShift;
        for (int i = subband->width - subband->colEndAddOn; i < subband->width; ++i)
        {
            int32_t quantVal = subband->qStepBase + ((qStepTblPtr[lastIdx] * subband->qStepMult) >> 3);
            bandBuf[i] *= rtengine::LIM(quantVal, 1, 0x168000);
        }
    } else {
        // prev. version
        std::int32_t qScale = q_step_tbl[subband->qParam % 6] >> (6 - subband->qParam / 6);
        if (subband->qParam / 6 >= 6) {
            qScale = q_step_tbl[subband->qParam % 6] * (1 << (subband->qParam / 6 + 26));
        }

        if (qScale != 1) {
            for (std::int32_t i = 0; i < subband->width; ++i) {
                bandBuf[i] *= qScale;
            }
        }
    }

    return true;
}

void crxHorizontal53(
    std::int32_t* lineBufLA,
    std::int32_t* lineBufLB,
    CrxWaveletTransform* wavelet,
    std::uint32_t tileFlag
)
{
    std::int32_t* band0Buf = wavelet->subband0Buf;
    std::int32_t* band1Buf = wavelet->subband1Buf;
    std::int32_t* band2Buf = wavelet->subband2Buf;
    std::int32_t* band3Buf = wavelet->subband3Buf;

    if (wavelet->width <= 1) {
        lineBufLA[0] = band0Buf[0];
        lineBufLB[0] = band2Buf[0];
    } else {
        if (tileFlag & E_HAS_TILES_ON_THE_LEFT) {
            lineBufLA[0] = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
            lineBufLB[0] = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
            ++band1Buf;
            ++band3Buf;
        } else {
            lineBufLA[0] = band0Buf[0] - ((band1Buf[0] + 1) >> 1);
            lineBufLB[0] = band2Buf[0] - ((band3Buf[0] + 1) >> 1);
        }

        ++band0Buf;
        ++band2Buf;

        for (int i = 0; i < wavelet->width - 3; i += 2) {
            std::int32_t delta = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
            lineBufLA[1] = band1Buf[0] + ((delta + lineBufLA[0]) >> 1);
            lineBufLA[2] = delta;

            delta = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
            lineBufLB[1] = band3Buf[0] + ((delta + lineBufLB[0]) >> 1);
            lineBufLB[2] = delta;

            ++band0Buf;
            ++band1Buf;
            ++band2Buf;
            ++band3Buf;
            lineBufLA += 2;
            lineBufLB += 2;
        }

        if (tileFlag & E_HAS_TILES_ON_THE_RIGHT) {
            const std::int32_t deltaA = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
            lineBufLA[1] = band1Buf[0] + ((deltaA + lineBufLA[0]) >> 1);

            const std::int32_t deltaB = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
            lineBufLB[1] = band3Buf[0] + ((deltaB + lineBufLB[0]) >> 1);

            if (wavelet->width & 1) {
                lineBufLA[2] = deltaA;
                lineBufLB[2] = deltaB;
            }
        } else if (wavelet->width & 1) {
            lineBufLA[1] = band1Buf[0] + ((lineBufLA[0] + band0Buf[0] - ((band1Buf[0] + 1) >> 1)) >> 1);
            lineBufLA[2] = band0Buf[0] - ((band1Buf[0] + 1) >> 1);

            lineBufLB[1] = band3Buf[0] + ((lineBufLB[0] + band2Buf[0] - ((band3Buf[0] + 1) >> 1)) >> 1);
            lineBufLB[2] = band2Buf[0] - ((band3Buf[0] + 1) >> 1);
        } else {
            lineBufLA[1] = lineBufLA[0] + band1Buf[0];
            lineBufLB[1] = lineBufLB[0] + band3Buf[0];
        }
    }
}

std::int32_t* crxIdwt53FilterGetLine(CrxPlaneComp* comp, std::int32_t level)
{
    std::int32_t* const result = comp->waveletTransform[level].lineBuf[
        (comp->waveletTransform[level].fltTapH - comp->waveletTransform[level].curH + 5) % 5 + 3
    ];
    --comp->waveletTransform[level].curH;
    return result;
}

bool crxIdwt53FilterDecode(CrxPlaneComp* comp, std::int32_t level, CrxQStep *qStep)
{
    if (comp->waveletTransform[level].curH) {
        return true;
    }

    CrxSubband* const sband = comp->subBands + 3 * level;
    CrxQStep* qStepLevel = qStep ? qStep + level : 0;

    if (comp->waveletTransform[level].height - 3 <= comp->waveletTransform[level].curLine && !(comp->tileFlag & E_HAS_TILES_ON_THE_BOTTOM)) {
        if (comp->waveletTransform[level].height & 1) {
            if (level) {
                if (!crxIdwt53FilterDecode(comp, level - 1, qStep)) {
                    return false;
                }
            } else if (!crxDecodeLineWithIQuantization(sband, qStepLevel)) {
                return false;
            }

            if (!crxDecodeLineWithIQuantization(sband + 1, qStepLevel)) {
                return false;
            }
        }
    } else {
        if (level) {
            if (!crxIdwt53FilterDecode(comp, level - 1, qStep)) {
                return false;
            }
        } else if (!crxDecodeLineWithIQuantization(sband, qStepLevel)) { // LL band
            return false;
        }

        if (
            !crxDecodeLineWithIQuantization(sband + 1, qStepLevel) // HL band
            || !crxDecodeLineWithIQuantization(sband + 2, qStepLevel) // LH band
            || !crxDecodeLineWithIQuantization(sband + 3, qStepLevel) // HH band
        ) {
            return false;
        }
    }

    return true;
}

bool crxIdwt53FilterTransform(CrxPlaneComp* comp, std::uint32_t level)
{
    CrxWaveletTransform* const wavelet = comp->waveletTransform + level;

    if (wavelet->curH) {
        return true;
    }

    if (wavelet->curLine >= wavelet->height - 3) {
        if (!(comp->tileFlag & E_HAS_TILES_ON_THE_BOTTOM)) {
            if (wavelet->height & 1) {
                if (level) {
                    if (!wavelet[-1].curH) {
                        if (!crxIdwt53FilterTransform(comp, level - 1)) {
                            return false;
                        }
                    }

                    wavelet->subband0Buf = crxIdwt53FilterGetLine(comp, level - 1);
                }

                const std::int32_t* band0Buf = wavelet->subband0Buf;
                const std::int32_t* band1Buf = wavelet->subband1Buf;
                const std::int32_t* const lineBufH0 = wavelet->lineBuf[wavelet->fltTapH + 3];
                std::int32_t* const lineBufH1 = wavelet->lineBuf[(wavelet->fltTapH + 1) % 5 + 3];
                std::int32_t* const lineBufH2 = wavelet->lineBuf[(wavelet->fltTapH + 2) % 5 + 3];

                std::int32_t* lineBufL0 = wavelet->lineBuf[0];
                std::int32_t* lineBufL1 = wavelet->lineBuf[1];
                wavelet->lineBuf[1] = wavelet->lineBuf[2];
                wavelet->lineBuf[2] = lineBufL1;

                // process L bands
                if (wavelet->width <= 1) {
                    lineBufL0[0] = band0Buf[0];
                } else {
                    if (comp->tileFlag & E_HAS_TILES_ON_THE_LEFT) {
                        lineBufL0[0] = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                        ++band1Buf;
                    } else {
                        lineBufL0[0] = band0Buf[0] - ((band1Buf[0] + 1) >> 1);
                    }

                    ++band0Buf;

                    for (int i = 0; i < wavelet->width - 3; i += 2) {
                        const std::int32_t delta = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                        lineBufL0[1] = band1Buf[0] + ((lineBufL0[0] + delta) >> 1);
                        lineBufL0[2] = delta;
                        ++band0Buf;
                        ++band1Buf;
                        lineBufL0 += 2;
                    }

                    if (comp->tileFlag & E_HAS_TILES_ON_THE_RIGHT) {
                        const std::int32_t delta = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                        lineBufL0[1] = band1Buf[0] + ((lineBufL0[0] + delta) >> 1);

                        if (wavelet->width & 1) {
                            lineBufL0[2] = delta;
                        }
                    } else if (wavelet->width & 1) {
                        const std::int32_t delta = band0Buf[0] - ((band1Buf[0] + 1) >> 1);
                        lineBufL0[1] = band1Buf[0] + ((lineBufL0[0] + delta) >> 1);
                        lineBufL0[2] = delta;
                    } else {
                        lineBufL0[1] = band1Buf[0] + lineBufL0[0];
                    }
                }

                // process H bands
                lineBufL0 = wavelet->lineBuf[0];
                lineBufL1 = wavelet->lineBuf[1];

                for (std::int32_t i = 0; i < wavelet->width; ++i) {
                    const std::int32_t delta = lineBufL0[i] - ((lineBufL1[i] + 1) >> 1);
                    lineBufH1[i] = lineBufL1[i] + ((delta + lineBufH0[i]) >> 1);
                    lineBufH2[i] = delta;
                }

                wavelet->curH += 3;
                wavelet->curLine += 3;
                wavelet->fltTapH = (wavelet->fltTapH + 3) % 5;
            } else {
                std::int32_t* const lineBufL2 = wavelet->lineBuf[2];
                const std::int32_t* const lineBufH0 = wavelet->lineBuf[wavelet->fltTapH + 3];
                std::int32_t* const lineBufH1 = wavelet->lineBuf[(wavelet->fltTapH + 1) % 5 + 3];
                wavelet->lineBuf[1] = lineBufL2;
                wavelet->lineBuf[2] = wavelet->lineBuf[1];

                for (std::int32_t i = 0; i < wavelet->width; ++i) {
                    lineBufH1[i] = lineBufH0[i] + lineBufL2[i];
                }

                wavelet->curH += 2;
                wavelet->curLine += 2;
                wavelet->fltTapH = (wavelet->fltTapH + 2) % 5;
            }
        }
    } else {
        if (level) {
            if (!wavelet[-1].curH && !crxIdwt53FilterTransform(comp, level - 1)) {
                return false;
            }

            wavelet->subband0Buf = crxIdwt53FilterGetLine(comp, level - 1);
        }

        const std::int32_t* band0Buf = wavelet->subband0Buf;
        const std::int32_t* band1Buf = wavelet->subband1Buf;
        const std::int32_t* band2Buf = wavelet->subband2Buf;
        const std::int32_t* band3Buf = wavelet->subband3Buf;

        std::int32_t* lineBufL0 = wavelet->lineBuf[0];
        std::int32_t* lineBufL1 = wavelet->lineBuf[1];
        const std::int32_t* const lineBufH0 = wavelet->lineBuf[wavelet->fltTapH + 3];
        std::int32_t* const lineBufH1 = wavelet->lineBuf[(wavelet->fltTapH + 1) % 5 + 3];
        std::int32_t* const lineBufH2 = wavelet->lineBuf[(wavelet->fltTapH + 2) % 5 + 3];

        wavelet->lineBuf[1] = wavelet->lineBuf[2];
        wavelet->lineBuf[2] = lineBufL1;

        // process L bands
        if (wavelet->width <= 1) {
            lineBufL0[0] = band0Buf[0];
            lineBufL1[0] = band2Buf[0];
        } else {
            if (comp->tileFlag & E_HAS_TILES_ON_THE_LEFT) {
                lineBufL0[0] = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                lineBufL1[0] = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
                ++band1Buf;
                ++band3Buf;
            } else {
                lineBufL0[0] = band0Buf[0] - ((band1Buf[0] + 1) >> 1);
                lineBufL1[0] = band2Buf[0] - ((band3Buf[0] + 1) >> 1);
            }

            ++band0Buf;
            ++band2Buf;

            for (int i = 0; i < wavelet->width - 3; i += 2) {
                std::int32_t delta = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                lineBufL0[1] = band1Buf[0] + ((delta + lineBufL0[0]) >> 1);
                lineBufL0[2] = delta;

                delta = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
                lineBufL1[1] = band3Buf[0] + ((delta + lineBufL1[0]) >> 1);
                lineBufL1[2] = delta;

                ++band0Buf;
                ++band1Buf;
                ++band2Buf;
                ++band3Buf;
                lineBufL0 += 2;
                lineBufL1 += 2;
            }

            if (comp->tileFlag & E_HAS_TILES_ON_THE_RIGHT) {
                const std::int32_t deltaA = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                lineBufL0[1] = band1Buf[0] + ((deltaA + lineBufL0[0]) >> 1);

                const std::int32_t deltaB = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
                lineBufL1[1] = band3Buf[0] + ((deltaB + lineBufL1[0]) >> 1);

                if (wavelet->width & 1) {
                    lineBufL0[2] = deltaA;
                    lineBufL1[2] = deltaB;
                }
            } else if (wavelet->width & 1) {
                std::int32_t delta = band0Buf[0] - ((band1Buf[0] + 1) >> 1);
                lineBufL0[1] = band1Buf[0] + ((delta + lineBufL0[0]) >> 1);
                lineBufL0[2] = delta;

                delta = band2Buf[0] - ((band3Buf[0] + 1) >> 1);
                lineBufL1[1] = band3Buf[0] + ((delta + lineBufL1[0]) >> 1);
                lineBufL1[2] = delta;
            } else {
                lineBufL0[1] = lineBufL0[0] + band1Buf[0];
                lineBufL1[1] = lineBufL1[0] + band3Buf[0];
            }
        }

        // process H bands
        lineBufL0 = wavelet->lineBuf[0];
        lineBufL1 = wavelet->lineBuf[1];
        const std::int32_t* lineBufL2 = wavelet->lineBuf[2];

        for (std::int32_t i = 0; i < wavelet->width; ++i) {
            const std::int32_t delta = lineBufL0[i] - ((lineBufL2[i] + lineBufL1[i] + 2) >> 2);
            lineBufH1[i] = lineBufL1[i] + ((delta + lineBufH0[i]) >> 1);
            lineBufH2[i] = delta;
        }

        if (wavelet->curLine >= wavelet->height - 3 && (wavelet->height & 1)) {
            wavelet->curH += 3;
            wavelet->curLine += 3;
            wavelet->fltTapH = (wavelet->fltTapH + 3) % 5;
        } else {
            wavelet->curH += 2;
            wavelet->curLine += 2;
            wavelet->fltTapH = (wavelet->fltTapH + 2) % 5;
        }
    }

    return true;
}

bool crxIdwt53FilterInitialize(CrxPlaneComp* comp, std::int32_t prevLevel, CrxQStep *qStep)
{
    if (prevLevel == 0) {
        return true;
    }

    for (int curLevel = 0, curBand = 0; curLevel < prevLevel; curLevel++, curBand += 3) {
        CrxQStep* qStepLevel = qStep ? qStep + curLevel : 0;
        CrxWaveletTransform* const wavelet = comp->waveletTransform + curLevel;

        if (curLevel) {
            wavelet[0].subband0Buf = crxIdwt53FilterGetLine(comp, curLevel - 1);
        } else if (!crxDecodeLineWithIQuantization(comp->subBands + curBand, qStepLevel)) {
            return false;
        }

        std::int32_t* lineBufH0 = wavelet->lineBuf[wavelet->fltTapH + 3];

        if (wavelet->height > 1) {
            if (
                !crxDecodeLineWithIQuantization(comp->subBands + curBand + 1, qStepLevel)
                || !crxDecodeLineWithIQuantization(comp->subBands + curBand + 2, qStepLevel)
                || !crxDecodeLineWithIQuantization(comp->subBands + curBand + 3, qStepLevel)
            ) {
                return false;
            }

            std::int32_t* const lineBufL0 = wavelet->lineBuf[0];
            const std::int32_t* const lineBufL1 = wavelet->lineBuf[1];
            std::int32_t* lineBufL2 = wavelet->lineBuf[2];

            if (comp->tileFlag & E_HAS_TILES_ON_THE_TOP) {
                crxHorizontal53(lineBufL0, wavelet->lineBuf[1], wavelet, comp->tileFlag);

                if (!crxDecodeLineWithIQuantization(comp->subBands + curBand + 3, qStepLevel)|| !crxDecodeLineWithIQuantization(comp->subBands + curBand + 2, qStepLevel)) {
                    return false;
                }

                const std::int32_t* band2Buf = wavelet->subband2Buf;
                const std::int32_t* band3Buf = wavelet->subband3Buf;

                // process L band
                if (wavelet->width <= 1) {
                    lineBufL2[0] = band2Buf[0];
                } else {
                    if (comp->tileFlag & E_HAS_TILES_ON_THE_LEFT) {
                        lineBufL2[0] = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
                        ++band3Buf;
                    } else {
                        lineBufL2[0] = band2Buf[0] - ((band3Buf[0] + 1) >> 1);
                    }

                    ++band2Buf;

                    for (int i = 0; i < wavelet->width - 3; i += 2) {
                        const std::int32_t delta = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
                        lineBufL2[1] = band3Buf[0] + ((lineBufL2[0] + delta) >> 1);
                        lineBufL2[2] = delta;

                        ++band2Buf;
                        ++band3Buf;
                        lineBufL2 += 2;
                    }

                    if (comp->tileFlag & E_HAS_TILES_ON_THE_RIGHT) {
                        const std::int32_t delta = band2Buf[0] - ((band3Buf[0] + band3Buf[1] + 2) >> 2);
                        lineBufL2[1] = band3Buf[0] + ((lineBufL2[0] + delta) >> 1);

                        if (wavelet->width & 1) {
                            lineBufL2[2] = delta;
                        }
                    } else if (wavelet->width & 1) {
                        const std::int32_t delta = band2Buf[0] - ((band3Buf[0] + 1) >> 1);

                        lineBufL2[1] = band3Buf[0] + ((lineBufL2[0] + delta) >> 1);
                        lineBufL2[2] = delta;
                    } else {
                        lineBufL2[1] = band3Buf[0] + lineBufL2[0];
                    }
                }

                // process H band
                for (std::int32_t i = 0; i < wavelet->width; ++i) {
                    lineBufH0[i] = lineBufL0[i] - ((lineBufL1[i] + lineBufL2[i] + 2) >> 2);
                }
            } else {
                crxHorizontal53(lineBufL0, wavelet->lineBuf[2], wavelet, comp->tileFlag);

                for (int i = 0; i < wavelet->width; ++i) {
                    lineBufH0[i] = lineBufL0[i] - ((lineBufL2[i] + 1) >> 1);
                }
            }

            if (!crxIdwt53FilterDecode(comp, curLevel, qStep) || !crxIdwt53FilterTransform(comp, curLevel)) {
                return false;
            }
        } else {
            if (!crxDecodeLineWithIQuantization(comp->subBands + curBand + 1, qStepLevel)) {
                return false;
            }

            const std::int32_t* band0Buf = wavelet->subband0Buf;
            const std::int32_t* band1Buf = wavelet->subband1Buf;

            // process H band
            if (wavelet->width <= 1) {
                lineBufH0[0] = band0Buf[0];
            } else {
                if (comp->tileFlag & E_HAS_TILES_ON_THE_LEFT) {
                    lineBufH0[0] = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                    ++band1Buf;
                } else {
                    lineBufH0[0] = band0Buf[0] - ((band1Buf[0] + 1) >> 1);
                }

                ++band0Buf;

                for (int i = 0; i < wavelet->width - 3; i += 2) {
                    std::int32_t delta = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                    lineBufH0[1] = band1Buf[0] + ((lineBufH0[0] + delta) >> 1);
                    lineBufH0[2] = delta;

                    ++band0Buf;
                    ++band1Buf;
                    lineBufH0 += 2;
                }

                if (comp->tileFlag & E_HAS_TILES_ON_THE_RIGHT) {
                    const std::int32_t delta = band0Buf[0] - ((band1Buf[0] + band1Buf[1] + 2) >> 2);
                    lineBufH0[1] = band1Buf[0] + ((lineBufH0[0] + delta) >> 1);
                    lineBufH0[2] = delta;
                } else if (wavelet->width & 1) {
                    const std::int32_t delta = band0Buf[0] - ((band1Buf[0] + 1) >> 1);
                    lineBufH0[1] = band1Buf[0] + ((lineBufH0[0] + delta) >> 1);
                    lineBufH0[2] = delta;
                } else {
                    lineBufH0[1] = band1Buf[0] + lineBufH0[0];
                }
            }

            ++wavelet->curLine;
            ++wavelet->curH;
            wavelet->fltTapH = (wavelet->fltTapH + 1) % 5;
        }
    }

    return true;
}

void crxFreeSubbandData(CrxImage* image, CrxPlaneComp* comp)
{
    if (comp->compBuf) {
        free(comp->compBuf);
        comp->compBuf = nullptr;
    }

    if (!comp->subBands) {
        return;
    }

    for (std::int32_t i = 0; i < image->subbandCount; ++i) {
        if (comp->subBands[i].bandParam) {
            free(comp->subBands[i].bandParam);
            comp->subBands[i].bandParam = nullptr;
        }

        comp->subBands[i].bandBuf = nullptr;
        comp->subBands[i].bandSize = 0;
    }
}

void crxConvertPlaneLine(
    CrxImage* img,
    int imageRow,
    int imageCol = 0,
    int plane = 0,
    const std::int32_t* lineData = nullptr,
    int lineLength = 0
)
{
    if (lineData) {
        std::uint64_t rawOffset = 4 * img->planeWidth * imageRow + 2 * imageCol;

        if (img->encType == 1) {
            const std::int32_t maxVal = 1 << (img->nBits - 1);
            const std::int32_t minVal = -maxVal;

            for (int i = 0; i < lineLength; ++i) {
                img->outBufs[plane][rawOffset + 2 * i] = rtengine::LIM(lineData[i], minVal, maxVal - 1);
            }
        } else if (img->encType == 3) {
            // copy to intermediate planeBuf
            rawOffset = plane * img->planeWidth * img->planeHeight + img->planeWidth * imageRow + imageCol;

            for (int i = 0; i < lineLength; ++i) {
                img->planeBuf[rawOffset + i] = lineData[i];
            }
        } else if (img->nPlanes == 4) {
            const std::int32_t median = 1 << (img->nBits - 1);
            const std::int32_t maxVal = (1 << img->nBits) - 1;

            for (int i = 0; i < lineLength; ++i) {
                img->outBufs[plane][rawOffset + 2 * i] = rtengine::LIM(median + lineData[i], 0, maxVal);
            }
        } else if (img->nPlanes == 1) {
            const std::int32_t maxVal = (1 << img->nBits) - 1;
            const std::int32_t median = 1 << (img->nBits - 1);

            rawOffset = img->planeWidth * imageRow + imageCol;

            for (int i = 0; i < lineLength; ++i) {
                img->outBufs[0][rawOffset + i] = rtengine::LIM(median + lineData[i], 0, maxVal);
            }
        }
    } else if (img->encType == 3 && img->planeBuf) {
        const std::int32_t planeSize = img->planeWidth * img->planeHeight;
        const std::int16_t* const plane0 = img->planeBuf + imageRow * img->planeWidth;
        const std::int16_t* const plane1 = plane0 + planeSize;
        const std::int16_t* const plane2 = plane1 + planeSize;
        const std::int16_t* const plane3 = plane2 + planeSize;

        const std::int32_t median = 1 << (img->medianBits - 1) << 10;
        const std::int32_t maxVal = (1 << img->medianBits) - 1;
        const std::uint32_t rawLineOffset = 4 * img->planeWidth * imageRow;

        // for this stage - all except imageRow is ignored
        for (int i = 0; i < img->planeWidth; ++i) {
            std::int32_t gr = median + (plane0[i] << 10) - 168 * plane1[i] - 585 * plane3[i];

            if (gr < 0) {
                gr = -(((std::abs(gr) + 512) >> 9) & ~1);
            } else {
                gr = ((std::abs(gr) + 512) >> 9) & ~1;
            }

            // Essentially R = round(median + P0 + 1.474*P3)
            std::int32_t val = (median + (plane0[i] << 10) + 1510 * plane3[i] + 512) >> 10;
            img->outBufs[0][rawLineOffset + 2 * i] = rtengine::LIM(val, 0, maxVal);
            // Essentially G1 = round(median + P0 + P2 - 0.164*P1 - 0.571*P3)
            val = (plane2[i] + gr + 1) >> 1;
            img->outBufs[1][rawLineOffset + 2 * i] = rtengine::LIM(val, 0, maxVal);
            // Essentially G1 = round(median + P0 - P2 - 0.164*P1 - 0.571*P3)
            val = (gr - plane2[i] + 1) >> 1;
            img->outBufs[2][rawLineOffset + 2 * i] = rtengine::LIM(val, 0, maxVal);
            // Essentially B = round(median + P0 + 1.881*P1)
            val = (median + (plane0[i] << 10) + 1927 * plane1[i] + 512) >> 10;
            img->outBufs[3][rawLineOffset + 2 * i] = rtengine::LIM(val, 0, maxVal);
        }
    }
}

bool crxParamInit(
    CrxBandParam** param,
    std::uint64_t subbandMdatOffset,
    std::uint64_t subbandDataSize,
    std::uint32_t subbandWidth,
    std::uint32_t subbandHeight,
    bool supportsPartial,
    std::uint32_t roundedBitsMask,
    LibRaw_abstract_datastream* input
)
{
    const std::int32_t progrDataSize =
        supportsPartial
            ? 0
            : sizeof(std::int32_t) * subbandWidth;
    const std::int32_t paramLength = 2 * subbandWidth + 4;

    std::uint8_t* paramBuf = static_cast<std::uint8_t*>(calloc(1, sizeof(CrxBandParam) + sizeof(std::int32_t) * paramLength + progrDataSize));

    if (!paramBuf) {
        return false;
    }

    *param = reinterpret_cast<CrxBandParam*>(paramBuf);

    paramBuf += sizeof(CrxBandParam);

    (*param)->paramData = reinterpret_cast<std::int32_t*>(paramBuf);
    (*param)->nonProgrData =
        progrDataSize
            ? (*param)->paramData + paramLength
            : nullptr;
    (*param)->subbandWidth = subbandWidth;
    (*param)->subbandHeight = subbandHeight;
    (*param)->roundedBits = 0;
    (*param)->curLine = 0;
    (*param)->roundedBitsMask = roundedBitsMask;
    (*param)->supportsPartial = supportsPartial;
    (*param)->bitStream.bitData = 0;
    (*param)->bitStream.bitsLeft = 0;
    (*param)->bitStream.mdatSize = subbandDataSize;
    (*param)->bitStream.curPos = 0;
    (*param)->bitStream.curBufSize = 0;
    (*param)->bitStream.curBufOffset = subbandMdatOffset;
    (*param)->bitStream.input = input;

    crxFillBuffer(&(*param)->bitStream);

    return true;
}

bool crxSetupSubbandData(
    CrxImage* img,
    CrxPlaneComp* planeComp,
    const CrxTile* tile,
    std::uint32_t mdatOffset
)
{
    long compDataSize = 0;
    long waveletDataOffset = 0;
    long compCoeffDataOffset = 0;
    const std::int32_t toSubbands = 3 * img->levels + 1;

    CrxSubband* const subbands = planeComp->subBands;

    // calculate sizes
    for (std::int32_t subbandNum = 0; subbandNum < toSubbands; ++subbandNum) {
        subbands[subbandNum].bandSize = subbands[subbandNum].width * sizeof(std::int32_t); // 4 bytes
        compDataSize += subbands[subbandNum].bandSize;
    }

    if (img->levels) {
        const std::int32_t encLevels =
            img->levels
                ? img->levels
                : 1;
        waveletDataOffset = (compDataSize + 7) & ~7;
        compDataSize = (sizeof(CrxWaveletTransform) * encLevels + waveletDataOffset + 7) & ~7;
        compCoeffDataOffset = compDataSize;

        // calc wavelet line buffer sizes (always at one level up from current)
        for (int level = 0; level < img->levels; ++level) {
            if (level < img->levels - 1) {
                compDataSize += 8 * sizeof(std::int32_t) * planeComp->subBands[3 * (level + 1) + 2].width;
            } else {
                compDataSize += 8 * sizeof(std::int32_t) * tile->width;
            }
        }
    }

    // buffer allocation
    planeComp->compBuf = static_cast<std::uint8_t*>(malloc(compDataSize));

    if (!planeComp->compBuf) {
        return false;
    }

    // subbands buffer and sizes initialisation
    const std::uint64_t subbandMdatOffset = img->mdatOffset + mdatOffset;
    std::uint8_t* subbandBuf = planeComp->compBuf;

    for (std::int32_t subbandNum = 0; subbandNum < toSubbands; ++subbandNum) {
        subbands[subbandNum].bandBuf = subbandBuf;
        subbandBuf += subbands[subbandNum].bandSize;
        subbands[subbandNum].mdatOffset = subbandMdatOffset + subbands[subbandNum].dataOffset;
    }

    // wavelet data initialisation
    if (img->levels) {
        CrxWaveletTransform* const waveletTransforms = reinterpret_cast<CrxWaveletTransform*>(planeComp->compBuf + waveletDataOffset);
        std::int32_t* paramData = reinterpret_cast<std::int32_t*>(planeComp->compBuf + compCoeffDataOffset);

        planeComp->waveletTransform = waveletTransforms;
        waveletTransforms[0].subband0Buf = reinterpret_cast<std::int32_t*>(subbands->bandBuf);

        for (int level = 0; level < img->levels; ++level) {
            const std::int32_t band = 3 * level + 1;

            std::int32_t transformWidth = 0;

            if (level >= img->levels - 1) {
                waveletTransforms[level].height = tile->height;
                transformWidth = tile->width;
            } else {
                waveletTransforms[level].height = subbands[band + 3].height;
                transformWidth = subbands[band + 4].width;
            }

            waveletTransforms[level].width = transformWidth;
            waveletTransforms[level].lineBuf[0] = paramData;
            waveletTransforms[level].lineBuf[1] = waveletTransforms[level].lineBuf[0] + transformWidth;
            waveletTransforms[level].lineBuf[2] = waveletTransforms[level].lineBuf[1] + transformWidth;
            waveletTransforms[level].lineBuf[3] = waveletTransforms[level].lineBuf[2] + transformWidth;
            waveletTransforms[level].lineBuf[4] = waveletTransforms[level].lineBuf[3] + transformWidth;
            waveletTransforms[level].lineBuf[5] = waveletTransforms[level].lineBuf[4] + transformWidth;
            waveletTransforms[level].lineBuf[6] = waveletTransforms[level].lineBuf[5] + transformWidth;
            waveletTransforms[level].lineBuf[7] = waveletTransforms[level].lineBuf[6] + transformWidth;
            waveletTransforms[level].curLine = 0;
            waveletTransforms[level].curH = 0;
            waveletTransforms[level].fltTapH = 0;
            waveletTransforms[level].subband1Buf = reinterpret_cast<std::int32_t*>(subbands[band].bandBuf);
            waveletTransforms[level].subband2Buf = reinterpret_cast<std::int32_t*>(subbands[band + 1].bandBuf);
            waveletTransforms[level].subband3Buf = reinterpret_cast<std::int32_t*>(subbands[band + 2].bandBuf);

            paramData = waveletTransforms[level].lineBuf[7] + transformWidth;
        }
    }

    // decoding params and bitstream initialisation
    for (std::int32_t subbandNum = 0; subbandNum < toSubbands; ++subbandNum) {
        if (subbands[subbandNum].dataSize) {
            bool supportsPartial = false;
            std::uint32_t roundedBitsMask = 0;

            if (planeComp->supportsPartial && subbandNum == 0) {
                roundedBitsMask = planeComp->roundedBitsMask;
                supportsPartial = true;
            }

            if (
                !crxParamInit(
                    &subbands[subbandNum].bandParam,
                    subbands[subbandNum].mdatOffset,
                    subbands[subbandNum].dataSize,
                    subbands[subbandNum].width,
                    subbands[subbandNum].height,
                    supportsPartial,
                    roundedBitsMask,
                    img->input
                )
            ) {
                return false;
            }
        }
    }

    return true;
}

} // namespace

bool DCraw::crxDecodePlane(void* p, std::uint32_t planeNumber)
{
    CrxImage* const img = static_cast<CrxImage*>(p);
    int imageRow = 0;

    for (int tRow = 0; tRow < img->tileRows; ++tRow) {
        int imageCol = 0;

        for (int tCol = 0; tCol < img->tileCols; ++tCol) {
            const CrxTile* const tile = img->tiles + tRow * img->tileRows + tCol;
            CrxPlaneComp* const planeComp = tile->comps + planeNumber;
            const std::uint64_t tileMdatOffset = tile->dataOffset + tile->mdatQPDataSize + tile->mdatExtraSize + planeComp->dataOffset;

            // decode single tile
            if (!crxSetupSubbandData(img, planeComp, tile, tileMdatOffset)) {
                return false;
            }

            if (img->levels) {
                if (!crxIdwt53FilterInitialize(planeComp, img->levels, tile->qStep)) {
                    return false;
                }

                for (int i = 0; i < tile->height; ++i) {
                    if (!crxIdwt53FilterDecode(planeComp, img->levels - 1, tile->qStep) || !crxIdwt53FilterTransform(planeComp, img->levels - 1)) {
                        return false;
                    }

                    const std::int32_t* const lineData = crxIdwt53FilterGetLine(planeComp, img->levels - 1);
                    crxConvertPlaneLine(img, imageRow + i, imageCol, planeNumber, lineData, tile->width);
                }
            } else {
                // we have the only subband in this case
                if (!planeComp->subBands->dataSize) {
                    memset(planeComp->subBands->bandBuf, 0, planeComp->subBands->bandSize);
                    return true;
                }

                for (int i = 0; i < tile->height; ++i) {
                    if (!crxDecodeLine(planeComp->subBands->bandParam, planeComp->subBands->bandBuf)) {
                        return false;
                    }

                    const std::int32_t* const lineData = reinterpret_cast<std::int32_t*>(planeComp->subBands->bandBuf);
                    crxConvertPlaneLine(img, imageRow + i, imageCol, planeNumber, lineData, tile->width);
                }
            }

            imageCol += tile->width;
        }

        imageRow += img->tiles[tRow * img->tileCols].height;
    }
    return true;
}

namespace
{

using crx_data_header_t = DCraw::CanonCR3Data::crx_data_header_t;

std::uint32_t crxReadQP(CrxBitstream *bitStrm, std::int32_t kParam)
{
    std::uint32_t qp = crxBitstreamGetZeros(bitStrm);
    if (qp >= 23)
        qp = crxBitstreamGetBits(bitStrm, 8);
    else if (kParam)
        qp = crxBitstreamGetBits(bitStrm, kParam) | (qp << kParam);
    return qp;
}

void crxDecodeGolombTop(CrxBitstream *bitStrm, std::int32_t width, std::int32_t *lineBuf, std::int32_t *kParam)
{
    lineBuf[0] = 0;
    while (width-- > 0)
    {
        lineBuf[1] = lineBuf[0];
        std::uint32_t qp = crxReadQP(bitStrm, *kParam);
        lineBuf[1] += -(qp & 1) ^ (qp >> 1);
        *kParam = crxPredictKParameter(*kParam, qp, 7);
        ++lineBuf;
    }
    lineBuf[1] = lineBuf[0] + 1;
}

void crxDecodeGolombNormal(CrxBitstream *bitStrm, std::int32_t width, std::int32_t *lineBuf0, std::int32_t *lineBuf1, std::int32_t *kParam)
{
    lineBuf1[0] = lineBuf0[1];
    std::int32_t deltaH = lineBuf0[1] - lineBuf0[0];
    while (width-- > 0)
    {
        lineBuf1[1] = crxPrediction(lineBuf1[0], lineBuf0[1], deltaH, lineBuf0[0] - lineBuf1[0]);
        std::uint32_t qp = crxReadQP(bitStrm, *kParam);
        lineBuf1[1] += -(qp & 1) ^ (qp >> 1);
        if (width) {
            deltaH = lineBuf0[2] - lineBuf0[1];
            *kParam = crxPredictKParameter(*kParam, (qp + 2 * std::abs(deltaH)) >> 1, 7);
            ++lineBuf0;
        } else {
            *kParam = crxPredictKParameter(*kParam, qp, 7);
        }
        ++lineBuf1;
    }
    lineBuf1[1] = lineBuf1[0] + 1;
}

bool crxMakeQStep(CrxImage *img, CrxTile *tile, std::int32_t *qpTable, std::uint32_t totalQP)
{
    if (img->levels > 3 || img->levels < 1) {
        return false;
    }
    int qpWidth = (tile->width >> 3) + ((tile->width & 7) != 0);
    int qpHeight = (tile->height >> 1) + (tile->height & 1);
    int qpHeight4 = (tile->height >> 2) + ((tile->height & 3) != 0);
    int qpHeight8 = (tile->height >> 3) + ((tile->height & 7) != 0);
    std::size_t totalHeight = qpHeight;
    if (img->levels > 1) {
        totalHeight += qpHeight4;
    }
    if (img->levels > 2) {
        totalHeight += qpHeight8;
    }

    tile->qStep = static_cast<CrxQStep*>(
        malloc(totalHeight * qpWidth * sizeof(std::uint32_t) + img->levels * sizeof(CrxQStep))
    );

    if (!tile->qStep) {
        return false;
    }
    std::uint32_t *qStepTbl = (std::uint32_t *)(tile->qStep + img->levels);
    CrxQStep *qStep = tile->qStep;
    switch (img->levels) {
    case 3:
        qStep->qStepTbl = qStepTbl;
        qStep->width = qpWidth;
        qStep->height = qpHeight8;
        for (int qpRow = 0; qpRow < qpHeight8; ++qpRow) {
            int row0Idx = qpWidth * std::min(4 * qpRow, qpHeight - 1);
            int row1Idx = qpWidth * std::min(4 * qpRow + 1, qpHeight - 1);
            int row2Idx = qpWidth * std::min(4 * qpRow + 2, qpHeight - 1);
            int row3Idx = qpWidth * std::min(4 * qpRow + 3, qpHeight - 1);

            for (int qpCol = 0; qpCol < qpWidth; ++qpCol, ++qStepTbl) {
                std::int32_t quantVal = qpTable[row0Idx++] + qpTable[row1Idx++] + qpTable[row2Idx++] + qpTable[row3Idx++];
                // not sure about this nonsense - why is it not just avg like with 2 levels?
                quantVal = ((quantVal < 0) * 3 + quantVal) >> 2;
                if (quantVal / 6 >= 6)
                    *qStepTbl = q_step_tbl[quantVal % 6] * (1 << (quantVal / 6 + 26));
                else
                    *qStepTbl = q_step_tbl[quantVal % 6] >> (6 - quantVal / 6);
            }
        }
        // continue to the next level - we always decode all levels
        ++qStep;
        // fall through
    case 2:
        qStep->qStepTbl = qStepTbl;
        qStep->width = qpWidth;
        qStep->height = qpHeight4;
        for (int qpRow = 0; qpRow < qpHeight4; ++qpRow) {
            int row0Idx = qpWidth * std::min(2 * qpRow, qpHeight - 1);
            int row1Idx = qpWidth * std::min(2 * qpRow + 1, qpHeight - 1);

            for (int qpCol = 0; qpCol < qpWidth; ++qpCol, ++qStepTbl) {
                std::int32_t quantVal = (qpTable[row0Idx++] + qpTable[row1Idx++]) / 2;
                if (quantVal / 6 >= 6)
                    *qStepTbl = q_step_tbl[quantVal % 6] * (1 << (quantVal / 6 + 26));
                else
                    *qStepTbl = q_step_tbl[quantVal % 6] >> (6 - quantVal / 6);
            }
        }
        // continue to the next level - we always decode all levels
        ++qStep;
        // fall through
    case 1:
        qStep->qStepTbl = qStepTbl;
        qStep->width = qpWidth;
        qStep->height = qpHeight;
        for (int qpRow = 0; qpRow < qpHeight; ++qpRow) {
            for (int qpCol = 0; qpCol < qpWidth; ++qpCol, ++qStepTbl, ++qpTable) {
                if (*qpTable / 6 >= 6)
                    *qStepTbl = q_step_tbl[*qpTable % 6] * (1 << (*qpTable / 6 + 26));
                else
                    *qStepTbl = q_step_tbl[*qpTable % 6] >> (6 - *qpTable / 6);
            }
        }
        break;
    }
    return true;
}

inline void crxSetupSubbandIdx(crx_data_header_t *hdr, CrxImage *img, CrxSubband *band, int level,
                                      short colStartIdx, short bandWidthExCoef, short rowStartIdx,
                                      short bandHeightExCoef)
{
  if (hdr->version == 0x200)
  {
    band->rowStartAddOn = rowStartIdx;
    band->rowEndAddOn = bandHeightExCoef;
    band->colStartAddOn = colStartIdx;
    band->colEndAddOn = bandWidthExCoef;
    band->levelShift = 3 - level;
  }
  else
  {
    band->rowStartAddOn = 0;
    band->rowEndAddOn = 0;
    band->colStartAddOn = 0;
    band->colEndAddOn = 0;
    band->levelShift = 0;
  }
}

bool crxReadSubbandHeaders( // Combined with crxProcessSubbands function
    crx_data_header_t* hdr,
    CrxImage* img,
    CrxTile* tile,
    CrxPlaneComp* comp,
    std::uint8_t** subbandMdatPtr,
    std::int32_t* mdatSize
)
{
    CrxSubband* band = comp->subBands + img->subbandCount - 1; // set to last band
    std::uint32_t bandHeight = tile->height;
    std::uint32_t bandWidth = tile->width;
    std::int32_t bandWidthExCoef = 0;
    std::int32_t bandHeightExCoef = 0;

    if (img->levels) {
        // Build up subband sequences to crxDecode to a level in a header

        // Coefficient structure is a bit unclear and convoluted:
        //   3 levels max - 8 groups (for tile width rounded to 8 bytes)
        //                  of 3 band per level 4 sets of coefficients for each
        const std::int32_t* rowExCoef = exCoefNumTbl + 0x30 * (img->levels - 1) + 6 * (tile->width & 7);
        const std::int32_t* colExCoef = exCoefNumTbl + 0x30 * (img->levels - 1) + 6 * (tile->height & 7);

        for (int level = 0; level < img->levels; ++level) {
            const std::int32_t widthOddPixel = bandWidth & 1;
            const std::int32_t heightOddPixel = bandHeight & 1;
            bandWidth = (widthOddPixel + bandWidth) >> 1;
            bandHeight = (heightOddPixel + bandHeight) >> 1;

            std::int32_t bandWidthExCoef0 = 0;
            std::int32_t bandWidthExCoef1 = 0;
            std::int32_t bandHeightExCoef0 = 0;
            std::int32_t bandHeightExCoef1 = 0;
            std::int32_t colStartIdx = 0;
            std::int32_t rowStartIdx = 0;

            if (tile->tileFlag & E_HAS_TILES_ON_THE_RIGHT) {
                bandWidthExCoef0 = rowExCoef[2 * level];
                bandWidthExCoef1 = rowExCoef[2 * level + 1];
            }

            if (tile->tileFlag & E_HAS_TILES_ON_THE_LEFT) {
                ++bandWidthExCoef0;
                colStartIdx = 1;
            }

            if (tile->tileFlag & E_HAS_TILES_ON_THE_BOTTOM) {
                bandHeightExCoef0 = colExCoef[2 * level];
                bandHeightExCoef1 = colExCoef[2 * level + 1];
            }

            if (tile->tileFlag & E_HAS_TILES_ON_THE_TOP) {
                ++bandHeightExCoef0;
                rowStartIdx = 1;
            }

            band[0].width = bandWidth + bandWidthExCoef0 - widthOddPixel;
            band[0].height = bandHeight + bandHeightExCoef0 - heightOddPixel;
            crxSetupSubbandIdx(hdr, img, band, level + 1, colStartIdx, bandWidthExCoef0 - colStartIdx, rowStartIdx, bandHeightExCoef0 - rowStartIdx);

            band[-1].width = bandWidth + bandWidthExCoef1;
            band[-1].height = bandHeight + bandHeightExCoef0 - heightOddPixel;
            crxSetupSubbandIdx(hdr, img, band - 1, level + 1, 0, bandWidthExCoef1, rowStartIdx, bandHeightExCoef0 - rowStartIdx);

            band[-2].width = bandWidth + bandWidthExCoef0 - widthOddPixel;
            band[-2].height = bandHeight + bandHeightExCoef1;
            crxSetupSubbandIdx(hdr, img, band - 2, level + 1, colStartIdx, bandWidthExCoef0 - colStartIdx, 0, bandHeightExCoef1);

            band -= 3;
        }

        bandWidthExCoef = 0;
        bandHeightExCoef = 0;

        if (tile->tileFlag & E_HAS_TILES_ON_THE_RIGHT) {
            bandWidthExCoef = rowExCoef[2 * img->levels - 1];
        }

        if (tile->tileFlag & E_HAS_TILES_ON_THE_BOTTOM) {
            bandHeightExCoef = colExCoef[2 * img->levels - 1];
        }
    }

    band->width = bandWidthExCoef + bandWidth;
    band->height = bandHeightExCoef + bandHeight;

    if (img->levels) {
        crxSetupSubbandIdx(hdr, img, band, img->levels, 0, bandWidthExCoef, 0, bandHeightExCoef);
    }
    // End of crxProcessSubbands

    // Begin of crxReadSubbandHeaders
    if (!img->subbandCount) {
        return true;
    }

    std::int32_t subbandOffset = 0;
    band = comp->subBands;

    for (unsigned int curSubband = 0; curSubband < img->subbandCount; curSubband++, band++) {
        if (*mdatSize < 4) {
            return false;
        }

        int hdrSign = sgetn(2, *subbandMdatPtr);
        int hdrSize = sgetn(2, *subbandMdatPtr + 2);

        if (*mdatSize < hdrSize + 4) {
            return false;
        }
        if ((hdrSign != 0xFF03 || hdrSize != 8) && (hdrSign != 0xFF13 || hdrSize != 16)) {
            return false;
        }

        const std::uint32_t subbandSize = sgetn(4, *subbandMdatPtr + 4);
        if (curSubband != ((*subbandMdatPtr)[8] & 0xF0) >> 4) {
            band->dataSize = subbandSize;
            return false;
        }

        band->dataOffset = subbandOffset;
        band->kParam = 0;
        band->bandParam = 0;
        band->bandBuf = 0;
        band->bandSize = 0;

        if (hdrSign == 0xFF03) {
            // old header
            std::uint32_t bitData = sgetn(4, *subbandMdatPtr + 8);
            band->dataSize = subbandSize - (bitData & 0x7FFFF);
            band->supportsPartial = bitData & 0x8000000;
            band->qParam = (bitData >> 19) & 0xFF;
            band->qStepBase = 0;
            band->qStepMult = 0;
        } else {
            // new header
            if (sgetn(2, *subbandMdatPtr + 8) & 0xFFF) {
                // partial and qParam are not supported
                return false;
            }
            if (sgetn(2, *subbandMdatPtr + 18)) {
                // new header terminated by 2 zero bytes
                return false;
            }
            band->supportsPartial = false;
            band->qParam = 0;
            band->dataSize = subbandSize - sgetn(2, *subbandMdatPtr + 16);
            band->qStepBase = sgetn(4, *subbandMdatPtr + 12);
            band->qStepMult = sgetn(2, *subbandMdatPtr + 10);
        }

        subbandOffset += subbandSize;

        *subbandMdatPtr += hdrSize + 4;
        *mdatSize -= hdrSize + 4;
    }

    return true;
}

bool crxReadImageHeaders(
    crx_data_header_t* hdr,
    CrxImage* img,
    std::uint8_t* mdatPtr,
    std::uint32_t mdatHdrSize
)
{
    const unsigned int nTiles = img->tileRows * img->tileCols;

    if (!nTiles) {
        return false;
    }

    if (!img->tiles) {
        img->tiles = static_cast<CrxTile*>(
            malloc(
                sizeof(CrxTile) * nTiles
                + sizeof(CrxPlaneComp) * nTiles * img->nPlanes
                + sizeof(CrxSubband) * nTiles * img->nPlanes * img->subbandCount
            )
        );

        if (!img->tiles) {
            return false;
        }

        // memory areas in allocated chunk
        CrxTile* tile = img->tiles;
        CrxPlaneComp* const comps = reinterpret_cast<CrxPlaneComp*>(tile + nTiles);
        CrxSubband* const bands = reinterpret_cast<CrxSubband*>(comps + img->nPlanes * nTiles);

        for (unsigned int curTile = 0; curTile < nTiles; curTile++, tile++) {
            tile->tileFlag = 0; // tile neighbouring flags
            tile->tileNumber = curTile;
            tile->tileSize = 0;
            tile->comps = comps + curTile * img->nPlanes;

            if ((curTile + 1) % img->tileCols) {
                // not the last tile in a tile row
                tile->width = hdr->tileWidth;

                if (img->tileCols > 1) {
                    tile->tileFlag = E_HAS_TILES_ON_THE_RIGHT;

                    if (curTile % img->tileCols) {
                        // not the first tile in tile row
                        tile->tileFlag |= E_HAS_TILES_ON_THE_LEFT;
                    }
                }
            } else {
                // last tile in a tile row
                tile->width = img->planeWidth - hdr->tileWidth * (img->tileCols - 1);

                if (img->tileCols > 1) {
                    tile->tileFlag = E_HAS_TILES_ON_THE_LEFT;
                }
            }

            if (curTile < nTiles - img->tileCols) {
                // in first tile row
                tile->height = hdr->tileHeight;

                if (img->tileRows > 1) {
                    tile->tileFlag |= E_HAS_TILES_ON_THE_BOTTOM;

                    if (curTile >= img->tileCols) {
                        tile->tileFlag |= E_HAS_TILES_ON_THE_TOP;
                    }
                }
            } else {
                // non first tile row
                tile->height = img->planeHeight - hdr->tileHeight * (img->tileRows - 1);

                if (img->tileRows > 1) {
                    tile->tileFlag |= E_HAS_TILES_ON_THE_TOP;
                }
            }

            if (img->nPlanes) {
                CrxPlaneComp* comp = tile->comps;
                CrxSubband* band = bands + curTile * img->nPlanes * img->subbandCount;

                for (int curComp = 0; curComp < img->nPlanes; curComp++, comp++) {
                    comp->compNumber = curComp;
                    comp->supportsPartial = true;
                    comp->tileFlag = tile->tileFlag;
                    comp->subBands = band;
                    comp->compBuf = nullptr;
                    comp->waveletTransform = nullptr;

                    if (img->subbandCount) {
                        for (int curBand = 0; curBand < img->subbandCount; curBand++, band++) {
                            band->supportsPartial = false;
                            band->qParam = 4;
                            band->bandParam = nullptr;
                            band->dataSize = 0;
                        }
                    }
                }
            }
        }
    }

    std::uint32_t tileOffset = 0;
    std::int32_t dataSize = mdatHdrSize;
    std::uint8_t* dataPtr = mdatPtr;
    CrxTile* tile = img->tiles;

    for (unsigned int curTile = 0; curTile < nTiles; curTile++, tile++) {
        if (dataSize < 4) {
            return false;
        }

        int hdrSign = sgetn(2, dataPtr);
        int hdrSize = sgetn(2, dataPtr + 2);
        if ((hdrSign != 0xFF01 || hdrSize != 8) && (hdrSign != 0xFF11 || (hdrSize != 8 && hdrSize != 16))) {
            return false;
        }
        if (dataSize < hdrSize + 4) {
            return false;
        }
        int tailSign = sgetn(2, dataPtr + 10);
        if ((hdrSize == 8 && tailSign) || (hdrSize == 16 && tailSign != 0x4000)) {
            return false;
        }
        if (sgetn(2, dataPtr + 8) != static_cast<unsigned>(curTile)) {
            return false;
        }

        dataSize -= hdrSize + 4;

        tile->tileSize = sgetn(4, dataPtr + 4);
        tile->dataOffset = tileOffset;
        tile->qStep = 0;

        if (hdrSize == 16) {
          // extended header data - terminated by 0 bytes
          if (sgetn(2, dataPtr + 18) != 0) {
              return false;
          }
          tile->hasQPData = true;
          tile->mdatQPDataSize = sgetn(4, dataPtr + 12);
          tile->mdatExtraSize = sgetn(2, dataPtr + 16);
        } else {
          tile->hasQPData = false;
          tile->mdatQPDataSize = 0;
          tile->mdatExtraSize = 0;
        }

        dataPtr += hdrSize + 4;
        tileOffset += tile->tileSize;

        std::uint32_t compOffset = 0;
        CrxPlaneComp* comp = tile->comps;

        for (int compNum = 0; compNum < img->nPlanes; compNum++, comp++) {
            if (dataSize < 0xC) {
                return false;
            }

            hdrSign = sgetn(2, dataPtr);
            hdrSize = sgetn(2, dataPtr + 2);
            if ((hdrSign != 0xFF02 && hdrSign != 0xFF12) || hdrSize != 8) {
                return false;
            }
            if (compNum != dataPtr[8] >> 4) {
                return false;
            }
            if (sgetn(3, dataPtr + 9) != 0) {
                return false;
            }

            comp->compSize = sgetn(4, dataPtr + 4);

            std::int32_t compHdrRoundedBits = (dataPtr[8] >> 1) & 3;
            comp->supportsPartial = (dataPtr[8] & 8) != 0;

            comp->dataOffset = compOffset;
            comp->tileFlag = tile->tileFlag;

            compOffset += comp->compSize;
            dataSize -= 0xC;
            dataPtr += 0xC;

            comp->roundedBitsMask = 0;

            if (compHdrRoundedBits) {
                if (img->levels || !comp->supportsPartial) {
                    return false;
                }

                comp->roundedBitsMask = 1 << (compHdrRoundedBits - 1);
            }

            if (!crxReadSubbandHeaders(hdr, img, tile, comp, &dataPtr, &dataSize)) {
                return false;
            }
        }
    }

    if (hdr->version != 0x200) {
        return true;
    }

    tile = img->tiles;
    for (unsigned int curTile = 0; curTile < nTiles; ++curTile, ++tile) {
        if (tile->hasQPData) {
            CrxBitstream bitStrm;
            bitStrm.bitData = 0;
            bitStrm.bitsLeft = 0;
            bitStrm.curPos = 0;
            bitStrm.curBufSize = 0;
            bitStrm.mdatSize = tile->mdatQPDataSize;
            bitStrm.curBufOffset = img->mdatOffset + tile->dataOffset;
            bitStrm.input = img->input;

            crxFillBuffer(&bitStrm);

            unsigned int qpWidth = (tile->width >> 3) + ((tile->width & 7) != 0);
            unsigned int qpHeight = (tile->height >> 1) + (tile->height & 1);
            unsigned long totalQP = static_cast<std::size_t>(qpHeight) * qpWidth;

            try {
                std::vector<std::int32_t> qpTable(totalQP + 2 * (qpWidth + 2));
                std::int32_t *qpCurElem = qpTable.data();
                // 2 lines padded with extra pixels at the start and at the end
                std::int32_t *qpLineBuf = qpTable.data() + totalQP;
                std::int32_t kParam = 0;
                for (unsigned qpRow = 0; qpRow < qpHeight; ++qpRow) {
                    std::int32_t *qpLine0 = qpRow & 1 ? qpLineBuf + qpWidth + 2 : qpLineBuf;
                    std::int32_t *qpLine1 = qpRow & 1 ? qpLineBuf : qpLineBuf + qpWidth + 2;

                    if (qpRow) {
                        crxDecodeGolombNormal(&bitStrm, qpWidth, qpLine0, qpLine1, &kParam);
                    } else {
                        crxDecodeGolombTop(&bitStrm, qpWidth, qpLine1, &kParam);
                    }

                    for (unsigned qpCol = 0; qpCol < qpWidth; ++qpCol) {
                        *qpCurElem++ = qpLine1[qpCol + 1] + 4;
                    }
                }

                // now we read QP data - build tile QStep
                if (!crxMakeQStep(img, tile, qpTable.data(), totalQP)) {
                    return false;
                }
            } catch (...) {
                return false;
            }
        }
    }

    return true;
}

bool crxSetupImageData(
    crx_data_header_t* hdr,
    CrxImage* img,
    std::int16_t* outBuf,
    std::uint64_t mdatOffset,
    std::uint32_t mdatSize,
    std::uint8_t* mdatHdrPtr,
    std::int32_t mdatHdrSize
)
{
    constexpr bool IncrBitTable[16] = {
        false, false, false, false, false, false, false, false, false, true, true, false, false, false, true, false
    };

    img->planeWidth = hdr->f_width;
    img->planeHeight = hdr->f_height;

    if (
        hdr->tileWidth < 0x16
        || hdr->tileHeight < 0x16
        || img->planeWidth > 0x7FFF
        || img->planeHeight > 0x7FFF
    ) {
        return false;
    }

    img->tileCols = (img->planeWidth + hdr->tileWidth - 1) / hdr->tileWidth;
    img->tileRows = (img->planeHeight + hdr->tileHeight - 1) / hdr->tileHeight;

    if (
        img->planeWidth - hdr->tileWidth * (img->tileCols - 1) < 0x16 ||
        img->planeHeight - hdr->tileHeight * (img->tileRows - 1) < 0x16
    ) {
        return false;
    }

    img->tiles = nullptr;
    img->levels = hdr->imageLevels;
    img->subbandCount = 3 * img->levels + 1; // 3 bands per level + one last LL
    img->nPlanes = hdr->nPlanes;
    img->nBits = hdr->nBits;
    img->encType = hdr->encType;
    img->samplePrecision = hdr->nBits + IncrBitTable[4 * hdr->encType + 2] + 1;
    img->mdatOffset = mdatOffset + hdr->mdatHdrSize;
    img->mdatSize = mdatSize;
    img->planeBuf = nullptr;
    img->outBufs[0] = img->outBufs[1] = img->outBufs[2] = img->outBufs[3] = nullptr;
    img->medianBits = hdr->medianBits;

    // The encoding type 3 needs all 4 planes to be decoded to generate row of
    // RGGB values. It seems to be using some other colour space for raw encoding
    // It is a massive buffer so ideallly it will need a different approach:
    // decode planes line by line and convert single line then without
    // intermediate plane buffer. At the moment though it's too many changes so
    // left as is.
    if (img->encType == 3 && img->nPlanes == 4 && img->nBits > 8) {
        img->planeBuf = static_cast<std::int16_t*>(
            malloc(static_cast<std::size_t>(img->planeHeight) * img->planeWidth * img->nPlanes * ((img->samplePrecision + 7) >> 3))
        );

        if (!img->planeBuf) {
            return false;
        }
    }

    const std::int32_t rowSize = 2 * img->planeWidth;

    if (img->nPlanes == 1) {
        img->outBufs[0] = outBuf;
    } else {
        switch (hdr->cfaLayout) {
            case 0: {
                // R G
                // G B
                img->outBufs[0] = outBuf;
                img->outBufs[1] = outBuf + 1;
                img->outBufs[2] = outBuf + rowSize;
                img->outBufs[3] = img->outBufs[2] + 1;
                break;
            }

            case 1: {
                // G R
                // B G
                img->outBufs[1] = outBuf;
                img->outBufs[0] = outBuf + 1;
                img->outBufs[3] = outBuf + rowSize;
                img->outBufs[2] = img->outBufs[3] + 1;
                break;
            }

            case 2: {
                // G B
                // R G
                img->outBufs[2] = outBuf;
                img->outBufs[3] = outBuf + 1;
                img->outBufs[0] = outBuf + rowSize;
                img->outBufs[1] = img->outBufs[0] + 1;
                break;
            }

            case 3: {
                // B G
                // G R
                img->outBufs[3] = outBuf;
                img->outBufs[2] = outBuf + 1;
                img->outBufs[1] = outBuf + rowSize;
                img->outBufs[0] = img->outBufs[1] + 1;
                break;
            }
        }
    }

    // read header
    return crxReadImageHeaders(hdr, img, mdatHdrPtr, mdatHdrSize);
}

void crxFreeImageData(CrxImage* img)
{
    if (img->tiles) {
        CrxTile* const tile = img->tiles;
        const int nTiles = img->tileRows * img->tileCols;

        for (std::int32_t curTile = 0; curTile < nTiles; ++curTile) {
            if (tile[curTile].comps) {
                for (std::int32_t curPlane = 0; curPlane < img->nPlanes; ++curPlane) {
                    crxFreeSubbandData(img, tile[curTile].comps + curPlane);
                }
            }
            if (tile[curTile].qStep) {
                free(tile[curTile].qStep);
            }
        }

        free(img->tiles);
        img->tiles = nullptr;
    }

    if (img->planeBuf) {
        free(img->planeBuf);
        img->planeBuf = nullptr;
    }
}

}   // namespace

void DCraw::crxLoadDecodeLoop(void* img, int nPlanes)
{
#ifdef _OPENMP
    bool results[4]; // nPlanes is always <= 4
    #pragma omp parallel for

    for (std::int32_t plane = 0; plane < nPlanes; ++plane) {
        results[plane] = crxDecodePlane(img, plane);
    }

    for (std::int32_t plane = 0; plane < nPlanes; ++plane) {
        if (!results[plane]) {
            derror();
        }
    }

#else

    for (std::int32_t plane = 0; plane < nPlanes; ++plane) {
        if (!crxDecodePlane(img, plane)) {
            derror();
        }
    }

#endif
}

void DCraw::crxConvertPlaneLineDf(void* p, int imageRow)
{
    crxConvertPlaneLine(static_cast<CrxImage*>(p), imageRow);
}

void DCraw::crxLoadFinalizeLoopE3(void* p, int planeHeight)
{
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 0; i < planeHeight; ++i) {
        crxConvertPlaneLineDf(p, i);
    }
}

void DCraw::crxLoadRaw()
{
    CrxImage img;

    if (RT_canon_CR3_data.crx_track_selected < 0 || 
        RT_canon_CR3_data.crx_track_selected >= RT_canon_CR3_data.CRXTRACKS_MAXCOUNT) {
        derror();
    }

    crx_data_header_t hdr = RT_canon_CR3_data.crx_header[RT_canon_CR3_data.crx_track_selected];

    LibRaw_abstract_datastream input = { ifp };
    img.input = &input; // libraw_internal_data.internal_data.input;

    // update sizes for the planes
    if (hdr.nPlanes == 4) {
        hdr.f_width >>= 1;
        hdr.f_height >>= 1;
        hdr.tileWidth >>= 1;
        hdr.tileHeight >>= 1;
    }

    //  /*imgdata.color.*/maximum = (1 << hdr.nBits) - 1;
    tiff_bps = hdr.nBits;

    std::uint8_t* const hdrBuf = static_cast<std::uint8_t*>(malloc(hdr.mdatHdrSize * 2));

    // read image header
#ifdef _OPENMP
    #pragma omp critical
#endif
    {
#ifndef _OPENMP
        /*libraw_internal_data.internal_data.input->*/ input.lock();
#endif
        /*libraw_internal_data.internal_data.input->*/ input.seek(data_offset, SEEK_SET);
        /*libraw_internal_data.internal_data.input->*/ input.read(hdrBuf, 1, hdr.mdatHdrSize);
#ifndef _OPENMP
        /*libraw_internal_data.internal_data.input->*/ input.unlock();
#endif
    }

    // parse and setup the image data
    if (!crxSetupImageData(&hdr, &img, reinterpret_cast<std::int16_t*>(raw_image), hdr.MediaOffset /*data_offset*/, hdr.MediaSize /*RT_canon_CR3_data.data_size*/, hdrBuf, hdr.mdatHdrSize*2)) {
        derror();
    }

    free(hdrBuf);

    crxLoadDecodeLoop(&img, hdr.nPlanes);

    if (img.encType == 3) {
        crxLoadFinalizeLoopE3(&img, img.planeHeight);
    }

    crxFreeImageData(&img);
}

bool DCraw::crxParseImageHeader(uchar* cmp1TagData, int nTrack, int size)
{
    if (nTrack < 0 || nTrack >= RT_canon_CR3_data.CRXTRACKS_MAXCOUNT) {
        return false;
    }

    if (!cmp1TagData) {
        return false;
    }

    crx_data_header_t* const hdr = &RT_canon_CR3_data.crx_header[nTrack];

    hdr->version = sgetn(2, cmp1TagData + 4);
    hdr->f_width = sgetn(4, cmp1TagData + 8);
    hdr->f_height = sgetn(4, cmp1TagData + 12);
    hdr->tileWidth = sgetn(4, cmp1TagData + 16);
    hdr->tileHeight = sgetn(4, cmp1TagData + 20);
    hdr->nBits = cmp1TagData[24];
    hdr->nPlanes = cmp1TagData[25] >> 4;
    hdr->cfaLayout = cmp1TagData[25] & 0xF;
    hdr->encType = cmp1TagData[26] >> 4;
    hdr->imageLevels = cmp1TagData[26] & 0xF;
    hdr->hasTileCols = cmp1TagData[27] >> 7;
    hdr->hasTileRows = (cmp1TagData[27] >> 6) & 1;
    hdr->mdatHdrSize = sgetn(4, cmp1TagData + 28);
    int extHeader = cmp1TagData[32] >> 7;
    int useMedianBits = 0;
    hdr->medianBits = hdr->nBits;

    if (extHeader && size >= 56 && hdr->nPlanes == 4)
        useMedianBits = cmp1TagData[56] >> 6 & 1;

    if (useMedianBits && size >= 84)
        hdr->medianBits = cmp1TagData[84];

    // validation
    if ((hdr->version != 0x100 && hdr->version != 0x200) || !hdr->mdatHdrSize) {
        return false;
    }

    if (hdr->encType == 1) {
        if (hdr->nBits > 15) {
            return false;
        }
    } else {
        if (hdr->encType && hdr->encType != 3) {
            return false;
        }

        if (hdr->nBits > 14) {
            return false;
        }
    }

    if (hdr->nPlanes == 1) {
        if (hdr->cfaLayout || hdr->encType || hdr->nBits != 8) {
            return false;
        }
    } else if (
        hdr->nPlanes != 4
        || hdr->f_width & 1
        || hdr->f_height & 1
        || hdr->tileWidth & 1
        || hdr->tileHeight & 1
        || hdr->cfaLayout > 3
        || hdr->nBits == 8
    ) {
        return false;
    }

    if (hdr->tileWidth > hdr->f_width || hdr->tileHeight > hdr->f_height) {
        return false;
    }

    if (hdr->imageLevels > 3 || hdr->hasTileCols > 1 || hdr->hasTileRows > 1) {
        return false;
    }

    return true;
}
