/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2010 Fabio Suprani <ffsup2@yahoo.it>
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

#pragma once

#include <cstdint>
#include <iostream>

#include "myfile.h"
#include <csetjmp>
#include "settings.h"

class DCraw
{
public:
	typedef unsigned short ushort;
	typedef unsigned char uchar;
	typedef unsigned short (*dcrawImage_t)[4];
#ifdef _WIN32
	typedef __int64 INT64;
	typedef unsigned __int64 UINT64;
#else
	typedef long long INT64;
	typedef unsigned long long UINT64;
#endif


	DCraw()
    :exif_base(-1)
    ,ciff_base(-1)
    ,ciff_len(0)
    ,ifp(nullptr),ofp(nullptr)
    ,order(0x4949)
    ,ifname(nullptr)
    ,meta_data(nullptr)
    ,shot_select(0)
    ,multi_out(0)
    ,row_padding(0)
	,float_raw_image(nullptr)
    ,image(nullptr)
    ,bright(1.)
    ,half_size(0),four_color_rgb(0),document_mode(0),highlight(0)
    ,verbose(0)
    ,use_auto_wb(0),use_camera_wb(0),use_camera_matrix(1)
    ,output_color(1),output_bps(8),output_tiff(0),med_passes(0),no_auto_bright(0)
    ,RT_whitelevel_from_constant(ThreeValBool::X)
    ,RT_blacklevel_from_constant(ThreeValBool::X)
    ,RT_matrix_from_constant(ThreeValBool::X)
    ,RT_baseline_exposure(0)
	,getbithuff(this,ifp,zero_after_ff)
	,nikbithuff(ifp)
    {
        shrink=0;
        memset(&hbd, 0, sizeof(hbd));
        aber[0]=aber[1]=aber[2]=aber[3]=1;
        gamm[0]=0.45;gamm[1]=4.5;gamm[2]=gamm[3]=gamm[4]=gamm[5]=0;
        user_mul[0]=user_mul[1]=user_mul[2]=user_mul[3]=0;
        greybox[0]=greybox[1]=0; greybox[2]=greybox[3]= UINT_MAX;
        RT_canon_CR3_data.CR3_CTMDtag = 0;
    }

protected:
    enum class CropMode : std::uint_fast16_t {                                   // RT
        NA = 0,                                                                  // RT
        FullFrameOnGfx = 1,                                                      // RT
        SportsFinderMode = 2,                                                    // RT
        ElectronicShutter1_25xCrop = 4                                           // RT
    };                                                                           // RT
    // stores the cropdata read from the file                                       RT
    struct CropData {                                                            // RT
        std::uint_fast16_t  width,                                               // RT
                            height,                                              // RT
                            top_margin,                                          // RT
                            left_margin;                                         // RT
        CropMode crop_mode = CropMode::NA;                                       // RT
    } read_crop;                                                                 // RT
    int exif_base, ciff_base, ciff_len;
    rtengine::IMFILE *ifp;
    FILE *ofp;
    short order;
    const char *ifname;
    char *meta_data;
    int xtrans[6][6],xtrans_abs[6][6];
    char cdesc[5], desc[512], make[64], model[64], model2[64], model3[64], artist[64];
    float flash_used, canon_ev, iso_speed, shutter, aperture, focal_len;
    time_t timestamp;
    unsigned shot_order, kodak_cbpp, filters, exif_cfa, unique_id;
    off_t    strip_offset, data_offset;
    off_t    thumb_offset, meta_offset, profile_offset;
    unsigned thumb_length, meta_length, profile_length;
    unsigned thumb_misc, *oprof, fuji_layout, shot_select, multi_out;
    unsigned tiff_nifds, tiff_samples, tiff_bps, tiff_compress;
    unsigned black, cblack[4102], maximum, mix_green, raw_color, zero_is_bad;
    unsigned zero_after_ff, is_raw, dng_version, is_foveon, data_error;
    unsigned tile_width, tile_length, gpsdata[32], load_flags, row_padding;

    struct fuji_q_table
    {
      int8_t *q_table; /* quantization table */
      int    raw_bits;
      int    total_values;
      int    max_grad;    // sdp val
      int    q_grad_mult; // quant_gradient multiplier
      int    q_base;
    };

    struct fuji_compressed_params
    {
        struct      fuji_q_table qt[4];
        void        *buf;
        int         max_bits;
        int         min_value;
        int         max_value;
        ushort      line_width;
    };

    struct int_pair {
        int value1;
        int value2;
    };

    enum _xt_lines
    {
        _R0=0,_R1,_R2,_R3,_R4,
        _G0,_G1,_G2,_G3,_G4,_G5,_G6,_G7,
        _B0,_B1,_B2,_B3,_B4,
        _ltotal
    };

    // tables of gradients for single sample level
    struct fuji_grads
    {
      int_pair grads[41];
      int_pair lossy_grads[3][5];
    };

    struct fuji_compressed_block {
        int         cur_bit;         // current bit being read (from left to right)
        int         cur_pos;         // current position in a buffer
        INT64       cur_buf_offset;  // offset of this buffer in a file
        unsigned	max_read_size;	 // Amount of data to be read
        int         cur_buf_size;    // buffer size
        uchar       *cur_buf;        // currently read block
        int         fillbytes;          // Counter to add extra byte for block size N*16
        rtengine::IMFILE      *input;
        fuji_grads even[3]; // tables of even gradients
        fuji_grads odd[3];  // tables of odd gradients
        ushort		*linealloc;
        ushort      *linebuf[_ltotal];
    };

    /**
     * Metadata for merged pixel-shift image.
     */
    struct MergedPixelshift
    {
        bool is_merged_pixelshift = false;
        unsigned sub_frame_shot_select;
    };

    struct SonyMeta
    {
        /// SR2SubIFD black levels tag 0x7310 exists.
        bool sr2subifd_black2 = false;
    };

    int fuji_total_lines, fuji_total_blocks, fuji_block_width, fuji_bits, fuji_raw_type, fuji_lossless;

    ushort raw_height, raw_width, height, width, top_margin, left_margin;
    ushort shrink, iheight, iwidth, fuji_width, thumb_width, thumb_height;
    unsigned raw_size;
    ushort *raw_image;
    float * float_raw_image;
    ushort white[8][8], curve[0x10000], cr2_slice[3], sraw_mul[4];
    int mask[8][4], flip, tiff_flip, colors;
    double pixel_aspect;
    double aber[4];
    double gamm[6];
    dcrawImage_t image;
    float bright, threshold, user_mul[4];

    int half_size, four_color_rgb, document_mode, highlight;
    int verbose, use_auto_wb, use_camera_wb, use_camera_matrix;
    int output_color, output_bps, output_tiff, med_passes;
    int no_auto_bright;
    unsigned greybox[4] ;
    enum class ThreeValBool { X = -1, F, T };
    ThreeValBool RT_whitelevel_from_constant;
    ThreeValBool RT_blacklevel_from_constant;
    ThreeValBool RT_matrix_from_constant;
    std::string RT_software;
    double RT_baseline_exposure;
    struct MergedPixelshift merged_pixelshift;
    struct SonyMeta sony_meta;

public:
    struct CanonCR3Data {
        // contents of tag CMP1 for relevant track in CR3 file
        struct crx_data_header_t {
            int32_t version;
            int32_t f_width;
            int32_t f_height;
            int32_t tileWidth;
            int32_t tileHeight;
            int32_t nBits;
            int32_t nPlanes;
            int32_t cfaLayout;
            int32_t encType;
            int32_t imageLevels;
            int32_t hasTileCols;
            int32_t hasTileRows;
            int32_t mdatHdrSize;
            int32_t medianBits;
            // Not from header, but from datastream
            uint32_t MediaSize;
            int64_t MediaOffset;
            uint32_t MediaType; /* 1 -> /C/RAW, 2-> JPEG */
        };
        static constexpr int CRXTRACKS_MAXCOUNT = 16;
        crx_data_header_t crx_header[CRXTRACKS_MAXCOUNT];
        int crx_track_selected;
        short CR3_CTMDtag;
    };

    struct PanasonicRW2Info {
        struct v8_tags_t
        {
            uint32_t tag39[6];
            uint16_t tag3A[6];
            uint16_t tag3B;
            uint16_t initial[4];
            uint16_t tag40a[17], tag40b[17], tag41[17];
            uint16_t stripe_count; // 0x42
            uint16_t tag43;
            int64_t  stripe_offsets[5]; //0x44
            uint16_t stripe_left[5]; // 0x45
            uint32_t stripe_compressed_size[5]; //0x46
            uint16_t stripe_width[5]; //0x47
            uint16_t stripe_height[5];
        };

        ushort bpp;
        ushort encoding;
        v8_tags_t v8tags;
        PanasonicRW2Info(): bpp(0), encoding(0), v8tags{} {}
    };

    bool isBayer() const
    {
        return (filters != 0 && filters != 9);
    }

    struct CanonLevelsData {
        unsigned cblack[4];
        unsigned white;
        bool black_ok;
        bool white_ok;
        CanonLevelsData(): cblack{0}, white{0}, black_ok(false), white_ok(false) {}
    };

protected:
    CanonCR3Data RT_canon_CR3_data;

    CanonLevelsData RT_canon_levels_data;

    PanasonicRW2Info RT_pana_info;

    float cam_mul[4], pre_mul[4], cmatrix[3][4], rgb_cam[3][4];

    void (DCraw::*write_thumb)();
    void (DCraw::*load_raw)(), (DCraw::*thumb_load_raw)();
    jmp_buf failure;

    unsigned huff[1024]; // static inside foveon_decoder

    struct decode {
      struct decode *branch[2];
      int leaf;
    } first_decode[2048], *second_decode, *free_decode;

    struct tiff_ifd {
      int new_sub_file_type, width, height, bps, comp, phint, offset, flip, samples, bytes;
      int tile_width, tile_length, sample_format, predictor;
      float shutter;
    } tiff_ifd[10];

    struct ph1 {
	  int format, key_off, tag_21a;
	  int black, split_col, black_col, split_row, black_row;
	  float tag_210;
    } ph1;
    struct hbd {
        off_t levels, unknown1, flatfield;
    } hbd;

    struct jhead {
      int algo, bits, high, wide, clrs, sraw, psv, restart, vpred[6];
      ushort quant[64], idct[64], *huff[20], *free[20], *row;
    };

    struct tiff_tag {
      ushort tag, type;
      int count;
      union { char c[4]; short s[2]; int i; } val;
    };

    struct tiff_hdr {
      ushort order, magic;
      int ifd;
      ushort pad, ntag;
      struct tiff_tag tag[23];
      int nextifd;
      ushort pad2, nexif;
      struct tiff_tag exif[4];
      ushort pad3, ngps;
      struct tiff_tag gpst[10];
      short bps[4];
      int rat[10];
      unsigned gps[26];
      char desc[512], make[64], model[64], soft[32], date[20], artist[64];
    };
protected:

int fcol (int row, int col);
void merror (void *ptr, const char *where);
void derror();
inline void derror(bool condition) {if(UNLIKELY(condition)) ++data_error;}
ushort sget2 (uchar *s);
ushort get2();
unsigned sget4 (uchar *s);
unsigned get4();
unsigned getint (int type);
float int_to_float (int i);
double getreal (int type);
void read_shorts (ushort *pixel, int count);
void cubic_spline(const int *x_, const int *y_, const int len);
void canon_600_fixed_wb (int temp);
int canon_600_color (int ratio[2], int mar);
void canon_600_auto_wb();
void canon_600_coeff();
void canon_600_load_raw();
void canon_600_correct();
int canon_s2is();
void redcine_load_raw();
void parse_redcine();

// getbithuff(int nbits, ushort *huff);
class getbithuff_t
{
public:
   getbithuff_t(DCraw *p,rtengine::IMFILE *&i, unsigned &z):parent(p),bitbuf(0),vbits(0),reset(0),ifp(i),zero_after_ff(z){}
   unsigned operator()(int nbits, ushort *huff);

private:
   void derror(){
	   parent->derror();
   }
   DCraw *parent;
   unsigned bitbuf;
   int vbits, reset;
   rtengine::IMFILE *&ifp;
   unsigned &zero_after_ff;
};
getbithuff_t getbithuff;

class nikbithuff_t
{
public:
   explicit nikbithuff_t(rtengine::IMFILE *&i):bitbuf(0),errors(0),vbits(0),ifp(i){}
   void operator()() {bitbuf = vbits = 0;};
   unsigned operator()(int nbits, ushort *huff);
   unsigned errorCount() { return errors; }
private:
   inline bool derror(bool condition){
       if (UNLIKELY(condition)) {
           ++errors;
       }
	   return condition;
   }
   unsigned bitbuf, errors;
   int vbits;
   rtengine::IMFILE *&ifp;
};
nikbithuff_t nikbithuff;

ushort * make_decoder_ref (const uchar **source);
ushort * make_decoder (const uchar *source);
void crw_init_tables (unsigned table, ushort *huff[2]);
int canon_has_lowbits();
void canon_load_raw();
int ljpeg_start (struct jhead *jh, int info_only);
void ljpeg_end (struct jhead *jh);
int ljpeg_diff (ushort *huff);
ushort * ljpeg_row (int jrow, struct jhead *jh);
void lossless_jpeg_load_raw();
void ljpeg_idct (struct jhead *jh);


void canon_sraw_load_raw();
void adobe_copy_pixel (unsigned row, unsigned col, ushort **rp);
void lossless_dng_load_raw();
void packed_dng_load_raw();
void deflate_dng_load_raw();
void init_fuji_main_qtable(fuji_compressed_params *params, uchar q_base);
void init_fuji_main_grads(const fuji_compressed_params *params, fuji_compressed_block *info);
void init_fuji_compr(struct fuji_compressed_params* info);
void fuji_fill_buffer(struct fuji_compressed_block *info);
void init_fuji_block(struct fuji_compressed_block* info, const struct fuji_compressed_params *params, INT64 raw_offset, unsigned dsize);
void copy_line_to_xtrans(struct fuji_compressed_block* info, int cur_line, int cur_block, int cur_block_width);
void copy_line_to_bayer(struct fuji_compressed_block* info, int cur_line, int cur_block, int cur_block_width);
void fuji_zerobits(struct fuji_compressed_block* info, int *count);
void fuji_read_code(struct fuji_compressed_block* info, int *data, int bits_to_read);
int fuji_decode_sample_even(struct fuji_compressed_block* info, const struct fuji_compressed_params* params, ushort* line_buf, int pos, struct fuji_grads* grad_params);
int fuji_decode_sample_odd(struct fuji_compressed_block* info, const struct fuji_compressed_params* params, ushort* line_buf, int pos, struct fuji_grads* grad_params);
void fuji_decode_interpolation_even(int line_width, ushort* line_buf, int pos);
void fuji_extend_generic(ushort *linebuf[_ltotal], int line_width, int start, int end);
void fuji_extend_red(ushort *linebuf[_ltotal], int line_width);
void fuji_extend_green(ushort *linebuf[_ltotal], int line_width);
void fuji_extend_blue(ushort *linebuf[_ltotal], int line_width);
void xtrans_decode_block(struct fuji_compressed_block* info, const struct fuji_compressed_params *params);
void fuji_bayer_decode_block(struct fuji_compressed_block* info, const struct fuji_compressed_params *params);
void fuji_decode_strip (fuji_compressed_params* params, int cur_block, INT64 raw_offset, unsigned dsize, uchar *q_bases);
void fuji_compressed_load_raw();
void fuji_decode_loop(fuji_compressed_params* common_info, int count, INT64* raw_block_offsets, unsigned *block_sizes, uchar *q_bases);
void parse_fuji_compressed_header();
void fuji_14bit_load_raw();
void pana8_decode_loop(void *data);
bool pana8_decode_strip(void* data, int stream);
void pentax_load_raw();
void nikon_load_raw();
int nikon_is_compressed();
int nikon_e995();
int nikon_e2100();
void nikon_3700();
int minolta_z2();
void ppm_thumb();
void ppm16_thumb();
void layer_thumb();
void rollei_thumb();
void rollei_load_raw();
int raw (unsigned row, unsigned col);
void phase_one_flat_field (int is_float, int nc);
void phase_one_correct();
void phase_one_load_raw();
void nikon_yuv_load_raw();
void kodak_c330_load_raw();
void kodak_c603_load_raw();
void samsung3_load_raw();
void parse_qt (int end);

// ph1_bithuff(int nbits, ushort *huff);
class ph1_bithuff_t {
public:
   ph1_bithuff_t(DCraw *p, rtengine::IMFILE *i, short &o):order(o),ifp(i),bitbuf(0),vbits(0){}
   unsigned operator()(int nbits, ushort *huff);
   unsigned operator()(int nbits);
   unsigned operator()();
    ushort get2() {
        uchar str[2] = { 0xff,0xff };
        fread (str, 1, 2, ifp);
        if (order == 0x4949) { /* "II" means little-endian */
            return str[0] | str[1] << 8;
        } else { /* "MM" means big-endian */
            return str[0] << 8 | str[1];
        }
    }
private:
    inline unsigned get4() {
        unsigned val = 0xffffff;
        uchar* str = (uchar*)&val;
        fread (str, 1, 4, ifp);
        if (order == 0x4949) {
#if __BYTE_ORDER__==__ORDER_LITTLE_ENDIAN__
            return val;
#else
            return str[0] | str[1] << 8 | str[2] << 16 | str[3] << 24;
#endif
        } else {
#if __BYTE_ORDER__==__ORDER_LITTLE_ENDIAN__
            return str[0] << 24 | str[1] << 16 | str[2] << 8 | str[3];
#else
            return val;
#endif
        }
   }

   short &order;
   rtengine::IMFILE* const ifp;
   UINT64 bitbuf;
   int vbits;
};

void phase_one_load_raw_c();
void hasselblad_correct();
void parse_hasselblad_gain();
void hasselblad_load_raw();
void leaf_hdr_load_raw();
void unpacked_load_raw();
void unpacked_load_raw_FujiDBP();
void sinar_4shot_load_raw();
void imacon_full_load_raw();
void packed_load_raw();
void nokia_load_raw();

class pana_bits_t{
public:
   pana_bits_t(rtengine::IMFILE *i, unsigned &u, unsigned enc):
    ifp(i), load_flags(u), vbits(0), encoding(enc) {}
   unsigned operator()(int nbits, unsigned *bytes=nullptr);
private:
   rtengine::IMFILE *ifp;
   unsigned &load_flags;
   uchar buf[0x4000];
   int vbits;
   unsigned encoding;
};

void panasonicC6_load_raw();
void panasonicC7_load_raw();
void panasonicC8_load_raw();

void canon_rmf_load_raw();
void panasonic_load_raw();
void olympus_load_raw();
void minolta_rd175_load_raw();
void quicktake_100_load_raw();
void kodak_radc_load_raw();
void samsung_load_raw();
void samsung2_load_raw();
void kodak_jpeg_load_raw();
void lossy_dng_load_raw();
void kodak_dc120_load_raw();
void eight_bit_load_raw();
void kodak_yrgb_load_raw();
void kodak_262_load_raw();
int  kodak_65000_decode (short *out, int bsize);
void kodak_65000_load_raw();
void kodak_ycbcr_load_raw();
void kodak_rgb_load_raw();
void kodak_thumb_load_raw();

// sony_decrypt(unsigned *data, int len, int start, int key);
class sony_decrypt_t{
public:
   explicit sony_decrypt_t() : p(0) {}
   void operator()(unsigned *data, int len, int start, int key);
private:
   unsigned pad[128], p;
};
sony_decrypt_t sony_decrypt;

void sony_load_raw();
void sony_arw_load_raw();
void sony_arw2_load_raw();
void sony_arq_load_raw(); // RT
void smal_decode_segment (unsigned seg[2][2], int holes);
void smal_v6_load_raw();

int median4 (int *p);
void fill_holes (int holes);
void smal_v9_load_raw();

void foveon_decoder (unsigned size, unsigned code);
void foveon_thumb();
void foveon_sd_load_raw();
void foveon_load_camf();
void foveon_load_raw();
void foveon_huff (ushort *huff);
void foveon_dp_load_raw();
const char * foveon_camf_param (const char *block, const char *param);
void *foveon_camf_matrix (unsigned dim[3], const char *name);
int foveon_fixed (void *ptr, int size, const char *name);
float foveon_avg (short *pix, int range[2], float cfilt);
short *foveon_make_curve (double max, double mul, double filt);
void foveon_make_curves(short **curvep, float dq[3], float div[3], float filt);
int foveon_apply_curve (short *curve, int i);
void foveon_interpolate();

void gamma_curve (double pwr, double ts, int mode, int imax);
void pseudoinverse (double (*in)[3], double (*out)[3], int size);
void cam_xyz_coeff (float rgb_cam[3][4], double cam_xyz[4][3]);
void pre_interpolate();
void crop_masked_pixels();

void tiff_get (unsigned base,	unsigned *tag, unsigned *type, unsigned *len, unsigned *save);
void parse_thumb_note (int base, unsigned toff, unsigned tlen);
int  parse_tiff_ifd (int base);
void parse_makernote (int base, int uptag);
void get_timestamp (int reversed);
void parse_exif (int base);
void parse_gps (int base);
void romm_coeff (float romm_cam[3][3]);
void parse_mos (int offset);
void linear_table (unsigned len);
void parse_kodak_ifd (int base);

void parse_minolta (int base);
int  parse_tiff (int base);
void apply_tiff();
void parse_external_jpeg();
void ciff_block_1030();
void parse_ciff (int offset, int length, int depth);
void parse_rollei();
void parse_sinar_ia();
void parse_phase_one (int base);
void parse_fuji (int offset);
int  parse_jpeg (int offset);
void parse_riff();
void parse_smal (int offset, int fsize);
void parse_cine();
char *foveon_gets (int offset, char *str, int len);
void parse_foveon();

void adobe_coeff (const char *make, const char *model);
void simple_coeff (int index);
short guess_byte_order (int words);
float find_green (int bps, int bite, int off0, int off1);
void identify();
void jpeg_thumb() {}  // not needed
bool dcraw_coeff_overrides(const char make[], const char model[], int iso_speed, short trans[12], int *black_level, int *white_level);
void shiftXtransMatrix( const int offsy, const int offsx) {
    char xtransTemp[6][6];
    for(int row = 0;row < 6; row++) {
        for(int col = 0;col < 6; col++) {
            xtransTemp[row][col] = xtrans[(row+offsy)%6][(col+offsx)%6];
        }
    }
    for(int row = 0;row < 6; row++) {
        for(int col = 0;col < 6; col++) {
            xtrans[row][col] = xtransTemp[row][col];
        }
    }
}

void nikon_14bit_load_raw(); // ported from LibRaw

//-----------------------------------------------------------------------------
// Canon CR3 support ported from LibRaw
//-----------------------------------------------------------------------------
void parse_canon_cr3();
void selectCRXTrack(unsigned short maxTrack);
int parseCR3(unsigned long long oAtomList,
             unsigned long long szAtomList, short &nesting,
             char *AtomNameStack, short &nTrack, short &TrackType);
bool crxDecodePlane(void *p, uint32_t planeNumber);
void crxLoadDecodeLoop(void *img, int nPlanes);
void crxConvertPlaneLineDf(void *p, int imageRow);
void crxLoadFinalizeLoopE3(void *p, int planeHeight);
void crxLoadRaw();
bool crxParseImageHeader(uchar *cmp1TagData, int nTrack, int size);
//-----------------------------------------------------------------------------

};
