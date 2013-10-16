/*
 *  This file is part of RawTherapee.
 *
 *  Created on: 20/nov/2010
 */

#include "rawimage.h"
#include "settings.h"
#include "colortemp.h"
#include "camconst.h"
#include "utils.h"
#ifdef WIN32
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif
#include "safegtk.h"

namespace rtengine{

extern const Settings* settings;

RawImage::RawImage(  const Glib::ustring name )
:data(NULL)
,prefilters(0)
,filename(name)
,profile_data(NULL)
,allocation(NULL)
{
	memset(maximum_c4, 0, sizeof(maximum_c4));
}

RawImage::~RawImage()
{
   if(ifp)
	   fclose(ifp);
   if( image )
		free(image);
   if(allocation){ delete [] allocation; allocation=NULL;}
   if(data){ delete [] data; data=NULL;}
   if(profile_data){ delete [] profile_data; profile_data=NULL;}
}

/* Similar to dcraw scale_colors for coeff. calculation, but without actual pixels scaling.
 * need pixels in data[][] available
 */
void RawImage::get_colorsCoeff( float *pre_mul_, float *scale_mul_, float *cblack_)

{
	unsigned  row, col, x, y, c, sum[8];
	unsigned  W = this->get_width();
	unsigned  H = this->get_height();
	int val;
	double dsum[8], dmin, dmax;

	for (int c = 0; c < 4; c++){
		cblack_[c] = (float) this->get_cblack(c);
		pre_mul_[c] = this->get_pre_mul(c);
	}
	if ( this->get_cam_mul(0) == -1 ) {
		memset(dsum, 0, sizeof dsum);
		for (row = 0; row < H; row += 8)
			for (col = 0; col < W ; col += 8) {
				memset(sum, 0, sizeof sum);
				for (y = row; y < row + 8 && y < H; y++)
					for (x = col; x < col + 8 && x < W; x++)
						for (int c = 0; c < 3; c++) {
							if (this->isBayer()) {
								c = FC(y, x);
								val = data[y][x];
							} else
								val = data[y][3*x+c];
							if (val > this->get_white(c) - 25)
								goto skip_block;
							if ((val -= cblack_[c]) < 0)
								val = 0;
							sum[c] += val;
							sum[c + 4]++;
							if ( this->isBayer())
								break;
						}
				for (c = 0; c < 8; c++)
					dsum[c] += sum[c];
skip_block: ;
			}
		for (int c = 0; c < 4; c++)
			if (dsum[c])
				pre_mul_[c] = dsum[c + 4] / dsum[c];
	}else{
		memset(sum, 0, sizeof sum);
		for (row = 0; row < 8; row++)
			for (col = 0; col < 8; col++) {
				int c = FC(row, col);
				if ((val = white[row][col] - cblack_[c]) > 0)
					sum[c] += val;
				sum[c + 4]++;
			}
		if (sum[0] && sum[1] && sum[2] && sum[3])
			for (int c = 0; c < 4; c++)
				pre_mul_[c] = (float) sum[c + 4] / sum[c];
		else if (this->get_cam_mul(0) && this->get_cam_mul(2)){
			pre_mul_[0] = this->get_cam_mul(0);
			pre_mul_[1] = this->get_cam_mul(1);
			pre_mul_[2] = this->get_cam_mul(2);
			pre_mul_[3] = this->get_cam_mul(3);
		}else
			fprintf(stderr, "Cannot use camera white balance.\n");
	}
	if (pre_mul_[3] == 0)
		pre_mul_[3] = this->get_colors() < 4 ? pre_mul_[1] : 1;
	for (dmin = DBL_MAX, dmax = c = 0; c < 4; c++) {
		if (dmin > pre_mul_[c])
			dmin = pre_mul_[c];
		if (dmax < pre_mul_[c])
			dmax = pre_mul_[c];
	}

	for (c = 0; c < 4; c++) {
		int sat = this->get_white(c) - this->get_cblack(c);
		scale_mul_[c] = (pre_mul_[c] /= dmax) * 65535.0 / sat;
		}
}

int RawImage::loadRaw (bool loadData, bool closeFile)
{
  ifname = filename.c_str();
  image = NULL;
  verbose = settings->verbose;
  oprof = NULL;

  ifp = gfopen (ifname);  // Maps to either file map or direct fopen
  if (!ifp) return 3;

  thumb_length = 0;
  thumb_offset = 0;
  thumb_load_raw = 0;
  use_camera_wb = 0;
  highlight = 1;
  half_size = 0;
  raw_image = 0;

  //***************** Read ALL raw file info
  identify ();
  if (!is_raw) {
    fclose(ifp);
    ifp=NULL;
    return 2;
  }

  if (flip==5)
     this->rotate_deg = 270;
  else if (flip==3)
     this->rotate_deg = 180;
  else if (flip==6)
     this->rotate_deg = 90;
  else if (flip % 90 == 0 && flip < 360)
     this->rotate_deg = flip;
  else
     this->rotate_deg = 0;

  if( loadData ){

	  use_camera_wb = 1;
	  shrink = 0;
	  if (settings->verbose) printf ("Loading %s %s image from %s...\n", make, model, filename.c_str());
	  iheight = height;
	  iwidth  = width;

      if (filters || colors == 1) {
        raw_image = (ushort *) calloc ((raw_height+7)*raw_width, 2);
        merror (raw_image, "main()");
      }

	  // dcraw needs this global variable to hold pixel data
	  image = (dcrawImage_t)calloc (height*width*sizeof *image + meta_length, 1);
	  meta_data = (char *) (image + height*width);
	  if(!image)
		  return 200;

	  if (setjmp (failure)) {
          if (image) { free (image); image=NULL; }
          if (raw_image) { free(raw_image); raw_image=NULL; }
		  fclose(ifp); ifp=NULL;
		  return 100;
	  }

	  // Load raw pixels data
	  fseek (ifp, data_offset, SEEK_SET);
	  (this->*load_raw)();

      if (raw_image) {
        crop_masked_pixels();
        free (raw_image);
        raw_image=NULL;
      }

	  // Load embedded profile
	  if (profile_length) {
	    profile_data = new char[profile_length];
		fseek ( ifp, profile_offset, SEEK_SET);
		fread ( profile_data, 1, profile_length, ifp);
	  }

      /*
        Setting the black level, there are three sources:
        dcraw single value 'black' or multi-value 'cblack', can be calculated or come
        from a hard-coded table or come from a stored value in the raw file, and
        finally there's 'black_c4' which are table values provided by RT camera constants.
        Any of these may or may not be set.

        We reduce these sources to one four channel black level, and do this by picking
        the highest found.
      */
      int black_c4[4] = { -1, -1, -1, -1 };

      CameraConstantsStore* ccs = CameraConstantsStore::getInstance();
      CameraConst *cc = ccs->get(make, model);
      if (cc) {
           for (int i = 0; i < 4; i++) {
                 if (RT_blacklevel_from_constant) {
                     black_c4[i] = cc->get_BlackLevel(i, iso_speed);
                 }
                 // load 4 channel white level here, will be used if available
                 if (RT_whitelevel_from_constant) {
                 maximum_c4[i] = cc->get_WhiteLevel(i, iso_speed);
             }
         }
      }
      if (black_c4[0] == -1) {
          // RT constants not set, bring in the DCRAW single channel black constant
          for (int c=0; c < 4; c++) black_c4[c] = black;
      }
      for (int c=0; c < 4; c++) {
          if (cblack[c] < black_c4[c]) {
              cblack[c] = black_c4[c];
          }
      }
  }

  if ( closeFile ) {
      fclose(ifp); ifp=NULL;
  }

  return 0;
}

unsigned short** RawImage::compress_image()
{
	if( !image )
		return NULL;
	if (filters) {
		if (!allocation) {
			allocation = new unsigned short[height * width];
			data = new unsigned short*[height];
			for (int i = 0; i < height; i++)
				data[i] = allocation + i * width;
		}
	} else {
		if (!allocation) {
			allocation = new unsigned short[3 * height * width];
			data = new unsigned short*[height];
			for (int i = 0; i < height; i++)
				data[i] = allocation + 3 * i * width;
		}
	}

    // copy pixel raw data: the compressed format earns space
	if (filters != 0) {
        #pragma omp parallel for
		for (int row = 0; row < height; row++)
			for (int col = 0; col < width; col++)
				this->data[row][col] = image[row * width + col][FC(row, col)];
	} else {
        #pragma omp parallel for
		for (int row = 0; row < height; row++)
			for (int col = 0; col < width; col++) {
				this->data[row][3 * col + 0] = image[row * width + col][0];
				this->data[row][3 * col + 1] = image[row * width + col][1];
				this->data[row][3 * col + 2] = image[row * width + col][2];
			}
	}
    free(image); // we don't need this anymore
    image=NULL;
    return data;
}

bool 
RawImage::is_supportedThumb() const 
{
    return ( (thumb_width * thumb_height) > 0 &&
    		   ( write_thumb == &rtengine::RawImage::jpeg_thumb ||
                 write_thumb == &rtengine::RawImage::ppm_thumb ||
                 thumb_load_raw == &rtengine::RawImage::kodak_thumb_load_raw ));
}

bool 
RawImage::get_thumbSwap() const
{
    return ((order == 0x4949) == (ntohs(0x1234) == 0x1234)) ? true : false; 
}

} //namespace rtengine

bool
DCraw::dcraw_coeff_overrides(const char make[], const char model[], const int iso_speed, short trans[12], int *black_level, int *white_level)
{
    static const struct {
        const char *prefix;
        int black_level, white_level; // set to -1 for no change
        short trans[12]; // set first value to 0 for no change
    } table[] = {

        { "Canon EOS 5D Mark III", 0, 0x3a98,  /* RT */
          { 6722,-635,-963,-4287,12460,2028,-908,2162,5668 } },
        { "Canon EOS 5D", 0, 0xe6c, /* RT */
          { 6319,-793,-614,-5809,13342,2738,-1132,1559,7971 } },
        { "Canon EOS 6D", 0, 0x3c82,
          { 7034,-804,-1014,-4420,12564,2058,-851,1994,5758 } },
        { "Canon EOS 7D", 0, 0x3510, /* RT - Colin Walker */
          { 5962,-171,-732,-4189,12307,2099,-911,1981,6304 } },
        { "Canon EOS 20D", 0, 0xfff,  /* RT */
          { 7590,-1646,-673,-4697,12411,2568,-627,1118,7295 } },
        { "Canon EOS 40D", 0, 0x3f60,  /* RT */
          { 6070,-531,-883,-5763,13647,2315,-1533,2582,6801 } },
        { "Canon EOS 60D", 0, 0x2ff7, /* RT - Colin Walker */
          { 5678,-179,-718,-4389,12381,2243,-869,1819,6380 } },
        { "Canon EOS 450D", 0, 0x390d, /* RT */
          { 6246,-1272,-523,-5075,12357,3075,-1035,1825,7333 } },
        { "Canon EOS 550D", 0, 0x3dd7, /* RT - Lebedev*/
          { 6519,-772,-703,-4994,12737,2519,-1387,2492,6175 } },
        { "Canon EOS-1D Mark III", 0, 0x3bb0,  /* RT */
          { 7406,-1592,-646,-4856,12457,2698,-432,726,7921 } },
        { "Canon PowerShot G10", 0, 0,  /* RT */
          { 12535,-5030,-796,-2711,10134,3006,-413,1605,5264 } },
        { "Canon PowerShot G12", 0, 0,  /* RT */
          { 12222,-4097,-1380,-2876,11016,2130,-888,1630,4434 } },


        { "Fujifilm X100", 0, 0,  /* RT - Colin Walker */
          { 10841,-3288,-807,-4652,12552,2344,-642,1355,7206 } },


        { "Nikon D200", 0, 0xfbc, /* RT */
          { 8498,-2633,-295,-5423,12869,2860,-777,1077,8124 } },
        { "Nikon D3000", 0, 0, /* RT */
          { 9211,-2521,-104,-6487,14280,2394,-754,1122,8033 } },
        { "Nikon D3100", 0, 0, /* RT */
          { 7729,-2212,-481,-5709,13148,2858,-1295,1908,8936 } },
        { "Nikon D3S", 0, 0, /* RT */
          { 8792,-2663,-344,-5221,12764,2752,-1491,2165,8121 } },
        { "Nikon D5200", 0, 0, // color matrix copied from D5200 DNG D65 matrix
          { 8322,-3112,-1047,-6367,14342,2179,-988,1638,6394 } },
        { "Nikon D7000", 0, 0, /* RT - Tanveer(tsk1979) */
          { 7530,-1942,-255,-4318,11390,3362,-926,1694,7649 } },
        { "Nikon D7100", 0, 0, // color matrix and WP copied from D7100 DNG D65 matrix
          { 8322,-3112,-1047,-6367,14342,2179,-988,1638,6394 } },
        { "Nikon D700", 0, 0, /* RT */
          { 8364,-2503,-352,-6307,14026,2492,-1134,1512,8156 } },
        { "Nikon COOLPIX A", 0, 0x3e00, // color matrix and WP copied from "COOLPIX A" DNG D65 matrix
          { 8198,-2239,-724,-4871,12389,2798,-1043,205,7181 } },


        { "Olympus E-30", 0, 0xfbc,  /* RT - Colin Walker */
          { 8510,-2355,-693,-4819,12520,2578,-1029,2067,7752 } },
        { "Olympus E-5", 0, 0xeec, /* RT - Colin Walker */
          { 9732,-2629,-999,-4899,12931,2173,-1243,2353,7457 } },
        { "Olympus E-P1", 0, 0xffd, /* RT - Colin Walker */
          { 8834,-2344,-804,-4691,12503,2448,-978,1919,7603 } },
        { "Olympus E-P2", 0, 0xffd, /* RT - Colin Walker */
          { 7758,-1619,-800,-5002,12886,2349,-985,1964,8305 } },
        { "Olympus E-P3", 0, 0, /* RT - Colin Walker */
          { 7041,-1794,-336,-3790,11192,2984,-1364,2625,6217 } },
        { "Olympus E-PL1s", 0, 0, /* RT - Colin Walker */
          { 9010,-2271,-838,-4792,12753,2263,-1059,2058,7589 } },
        { "Olympus E-PL1", 0, 0, /* RT - Colin Walker */
          { 9010,-2271,-838,-4792,12753,2263,-1059,2058,7589 } },
        { "Olympus E-PL2", 0, 0, /* RT - Colin Walker */
          { 11975,-3351,-1184,-4500,12639,2061,-1230,2353,7009 } },
        { "Olympus E-PL3", 0, 0, /* RT - Colin Walker */
          { 7041,-1794,-336,-3790,11192,2984,-1364,2625,6217 } },
        { "Olympus XZ-1", 0, 0, /* RT - Colin Walker */
          { 8665,-2247,-762,-2424,10372,2382,-1011,2286,5189 } },


        { "Panasonic DMC-FZ150", 143, 0xfff,  /* RT */
          { 10435,-3208,-72,-2293,10506,2067,-486,1725,4682 } },
        { "Panasonic DMC-G10", 15, 0xf3c, /* RT - Colin Walker */
          { 8310,-1811,-960,-4941,12990,2151,-1378,2468,6860 } },
        { "Panasonic DMC-G1", 15, 0xf94,    /* RT - Colin Walker*/
          { 7477,-1615,-651,-5016,12769,2506,-1380,2475,7240 } },
        { "Panasonic DMC-G2", 15, 0xf3c, /* RT - Colin Walker */
          { 8310,-1811,-960,-4941,12990,2151,-1378,2468,6860 } },
        { "Panasonic DMC-G3", 143, 0xfff,   /* RT - Colin Walker */
          { 6051,-1406,-671,-4015,11505,2868,-1654,2667,6219 } },
        { "Panasonic DMC-G5", 143, 0xfff,   /* RT */
          { 7122,-2092,-419,-4643,11769,3283,-1363,2413,5944 } },
        { "Panasonic DMC-GF1", 15, 0xf92, /* RT - Colin Walker */
          { 7863,-2080,-668,-4623,12331,2578,-1020,2066,7266 } },
        { "Panasonic DMC-GF2", 143, 0xfff, /* RT - Colin Walker */
          { 7694,-1791,-745,-4917,12818,2332,-1221,2322,7197 } },
        { "Panasonic DMC-GF3", 143, 0xfff, /* RT - Colin Walker */
          { 8074,-1846,-861,-5026,12999,2239,-1320,2375,7422 } },
        { "Panasonic DMC-GH1", 15, 0xf92,  /* RT - Colin Walker */
          { 6360,-1557,-375,-4201,11504,3086,-1378,2518,5843 } },
        { "Panasonic DMC-GH2", 15, 0xf95, /* RT - Colin Walker */
//        { 6855,-1765,-456,-4223,11600,2996,-1450,2602,5761 } }, disabled
          { 7780,-2410,-806,-3913,11724,2484,-1018,2390,5298 } }, // dcraw original


        { "Pentax K200D", 0, 0,  /* RT */
          { 10962,-4428,-542,-5486,13023,2748,-569,842,8390 } },


        { "Leica Camera AG M9 Digital Camera", 0, 0,  /* RT */
          { 7181,-1706,-55,-3557,11409,2450,-621,2072,7533 } },


        { "Sony DSLR-A700", 126, 0, /* RT */
          { 6509,-1333,-137,-6171,13621,2824,-1490,2226,6952 } },
        { "Sony DSLR-A900", 128, 0, /* RT */
          { 5715,-1433,-410,-5603,12937,2989,-644,1247,8372 } },
        { "SONY NEX-3", 128, 0, /* RT - Colin Walker */
          { 5145,-741,-123,-4915,12310,2945,-794,1489,6906 } },
        { "SONY NEX-5", 128, 0, /* RT - Colin Walker */
          { 5154,-716,-115,-5065,12506,2882,-988,1715,6800 } },
        { "Sony NEX-5N", 128, 0, /* RT - Colin Walker */
          { 5130,-1055,-269,-4473,11797,3050,-701,1310,7121 } },
        { "Sony NEX-5R", 128, 0,
          { 6129,-1545,-418,-4930,12490,2743,-977,1693,6615 } },
        { "SONY NEX-C3", 128, 0,  /* RT - Colin Walker */
          { 5130,-1055,-269,-4473,11797,3050,-701,1310,7121 } },
        { "Sony SLT-A77", 128, 0,  /* RT - Colin Walker */
          { 5126,-830,-261,-4788,12196,2934,-948,1602,7068 } },
    };

    *black_level = -1;
    *white_level = -1;
    memset(trans, 0, sizeof(trans));

    // indicate that DCRAW wants these from constants (rather than having loaded these from RAW file
    // note: this is simplified so far, in some cases dcraw calls this when it has say the black level
    // from file, but then we will not provide any black level in the tables. This case is mainly just
    // to avoid loading table values if we have loaded a DNG conversion of a raw file (which already
    // have constants stored in the file).
    RT_whitelevel_from_constant = 1;
    RT_blacklevel_from_constant = 1;
    RT_matrix_from_constant = 1;

    { // test if we have any information in the camera constants store, if so we take that.
        rtengine::CameraConstantsStore* ccs = rtengine::CameraConstantsStore::getInstance();
        rtengine::CameraConst *cc = ccs->get(make, model);
        if (cc) {
            *black_level = cc->get_BlackLevel(0, iso_speed);
            *white_level = cc->get_WhiteLevel(0, iso_speed);
            if (cc->has_dcrawMatrix()) {
                const short *mx = cc->get_dcrawMatrix();
                for (int j = 0; j < 12; j++) {
                    trans[j] = mx[j];
                }
            }
            return true;
        }
    }

    char name[strlen(make) + strlen(model) + 32];
    sprintf(name, "%s %s", make, model);
    for (int i = 0; i < sizeof table / sizeof(table[0]); i++) {
        if (strcasecmp(name, table[i].prefix) == 0) {
            *black_level = table[i].black_level;
            *white_level = table[i].white_level;
            for (int j = 0; j < 12; j++) {
                trans[j] = table[i].trans[j];
            }
            return true;
        }
    }
    return false;
}
