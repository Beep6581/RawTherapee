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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DCRAW_H
#define DCRAW_H

#include "myfile.h"
#include <csetjmp>


class DCraw
{
public:
	typedef unsigned short ushort;
	typedef unsigned char uchar;
	typedef unsigned short (*dcrawImage_t)[4];
#ifdef WIN32
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
    ,ifp(NULL),ofp(NULL)
    ,order(0x4949)
    ,ifname(NULL)
    ,meta_data(NULL)
    ,shot_select(0),multi_out(0)
    ,image(NULL)
    ,bright(1.),threshold(0.)
    ,half_size(0),four_color_rgb(0),document_mode(0),highlight(0)
    ,verbose(0)
    ,RT_whitelevel_from_constant(0)
    ,RT_blacklevel_from_constant(0)
    ,RT_matrix_from_constant(0)
    ,use_auto_wb(0),use_camera_wb(0),use_camera_matrix(1)
    ,output_color(1),output_bps(8),output_tiff(0),med_passes(0),no_auto_bright(0)
	,getbithuff(this,ifp,zero_after_ff)
	,ph1_bithuff(this,ifp,order)	
	,pana_bits(ifp,load_flags)
    {
        aber[0]=aber[1]=aber[2]=aber[3]=1;
        gamm[0]=0.45;gamm[1]=4.5;gamm[2]=gamm[3]=gamm[4]=gamm[5]=0;
        user_mul[0]=user_mul[1]=user_mul[2]=user_mul[3]=0;
        greybox[0]=greybox[1]=0; greybox[2]=greybox[3]= UINT_MAX;
    }

    //int main (int argc, const char **argv);
protected:
    int exif_base, ciff_base, ciff_len;
    IMFILE *ifp;
    FILE *ofp;
    short order;
    const char *ifname;
    char *meta_data, xtrans[6][6];
    char cdesc[5], desc[512], make[64], model[64], model2[64], artist[64];
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
    unsigned tile_width, tile_length, gpsdata[32], load_flags;
    ushort raw_height, raw_width, height, width, top_margin, left_margin;
    ushort shrink, iheight, iwidth, fuji_width, thumb_width, thumb_height;
    ushort *raw_image;
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
    int RT_matrix_from_constant;
    int RT_blacklevel_from_constant;
    int RT_whitelevel_from_constant;

    float cam_mul[4], pre_mul[4], cmatrix[3][4], rgb_cam[3][4];
    
    int histogram[4][0x2000];
    void (DCraw::*write_thumb)(), (DCraw::*write_fun)();
    void (DCraw::*load_raw)(), (DCraw::*thumb_load_raw)();
    jmp_buf failure;

    unsigned huff[1024]; // static inside foveon_decoder

    struct decode {
      struct decode *branch[2];
      int leaf;
    } first_decode[2048], *second_decode, *free_decode;

    struct tiff_ifd {
      int width, height, bps, comp, phint, offset, flip, samples, bytes;
      int tile_width, tile_length;
    } tiff_ifd[10];

    struct ph1 {
      int format, key_off, black, black_off, split_col, tag_21a;
      float tag_210;
    } ph1;

    struct jhead {
      int bits, high, wide, clrs, sraw, psv, restart, vpred[6];
      ushort *huff[6], *free[4], *row;
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
ushort sget2 (uchar *s);
ushort get2();
unsigned sget4 (uchar *s);
unsigned get4();
unsigned getint (int type);
float int_to_float (int i);
double getreal (int type);
void read_shorts (ushort *pixel, int count);
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
   getbithuff_t(DCraw *p,IMFILE *&i, unsigned &z):parent(p),bitbuf(0),vbits(0),reset(0),ifp(i),zero_after_ff(z){}
   unsigned operator()(int nbits, ushort *huff);

private:
   void derror(){
	   parent->derror();
   }
   DCraw *parent;
   unsigned bitbuf;
   int vbits, reset;
   IMFILE *&ifp;
   unsigned &zero_after_ff;
};
getbithuff_t getbithuff;

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

void canon_sraw_load_raw();
void adobe_copy_pixel (unsigned row, unsigned col, ushort **rp);
void lossless_dng_load_raw();
void packed_dng_load_raw();
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

// ph1_bithuff(int nbits, ushort *huff);
class ph1_bithuff_t {
public:
   ph1_bithuff_t(DCraw *p,IMFILE *&i,short &o):parent(p),order(o),ifp(i),bitbuf(0),vbits(0){}
   unsigned operator()(int nbits, ushort *huff);
private:
   unsigned get4(){
	 return parent->get4();
   }
   DCraw *parent;
   short &order;
   IMFILE *&ifp;
   UINT64 bitbuf;
   int vbits;
};
ph1_bithuff_t ph1_bithuff;

void phase_one_load_raw_c();
void hasselblad_load_raw();
void leaf_hdr_load_raw();
void unpacked_load_raw();
void sinar_4shot_load_raw();
void imacon_full_load_raw();
void packed_load_raw();
void nokia_load_raw();

// pana_bits(int nbits);
class pana_bits_t{
public:
   pana_bits_t(IMFILE *&i,unsigned &u):ifp(i),load_flags(u){}
   unsigned operator()(int nbits);
private:
   IMFILE *&ifp;
   unsigned &load_flags;
   uchar buf[0x4000];
   int vbits;
};
pana_bits_t pana_bits;

void canon_rmf_load_raw();
void panasonic_load_raw();
void olympus_load_raw();
void minolta_rd175_load_raw();
void quicktake_100_load_raw();
void kodak_radc_load_raw();
void samsung_load_raw();

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
   void operator()(unsigned *data, int len, int start, int key);
private:
   unsigned pad[128], p;
};
sony_decrypt_t sony_decrypt;

void sony_load_raw();
void sony_arw_load_raw();
void sony_arw2_load_raw();
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

void xtrans_interpolate (int passes);
void cielab (ushort rgb[3], short lab[3]);

void remove_zeroes();
void bad_pixels (const char *cfname);
void subtract (const char *fname);
void gamma_curve (double pwr, double ts, int mode, int imax);
void pseudoinverse (double (*in)[3], double (*out)[3], int size);
void cam_xyz_coeff (float rgb_cam[3][4], double cam_xyz[4][3]);
void hat_transform (float *temp, float *base, int st, int size, int sc);
void wavelet_denoise();
void scale_colors();
void pre_interpolate();
void border_interpolate (int border);
void median_filter();
void blend_highlights();
void recover_highlights();
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
void apply_profile (const char *input, const char *output);
void jpeg_thumb() {}  // not needed
bool dcraw_coeff_overrides(const char make[], const char model[], int iso_speed, short trans[12], int *black_level, int *white_level);

};


#endif //DCRAW_H
