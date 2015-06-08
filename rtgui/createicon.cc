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
#include <png.h>

void png_read_data(png_structp png_ptr, png_bytep data, png_size_t length) {
   png_size_t check;

   /* fread() returns 0 on error, so it is OK to store this in a png_size_t
    * instead of an int, which is what fread() actually returns.
    */
   check = (png_size_t)fread(data, (png_size_t)1, length, (FILE *)png_ptr->io_ptr);

   if (check != length)
   {
      png_error(png_ptr, "Read Error");
   }
}

void png_write_data(png_structp png_ptr, png_bytep data, png_size_t length) {
   png_uint_32 check;

   check = fwrite(data, 1, length, (FILE *)(png_ptr->io_ptr));
   if (check != length)
   {
      png_error(png_ptr, "Write Error");
   }
}

void png_flush(png_structp png_ptr) {
   FILE *io_ptr;
   io_ptr = (FILE *)CVT_PTR((png_ptr->io_ptr));
   if (io_ptr != NULL)
      fflush(io_ptr);
}


unsigned char* loadPNG  (char* fname, int& w, int& h) {

    FILE *file = fopen (fname,"rb");

	unsigned char header[8];
	fread (header, 1, 8, file);

	png_structp png = png_create_read_struct (PNG_LIBPNG_VER_STRING, 0, 0, 0);
	png_infop info = png_create_info_struct (png);
	png_infop end_info = png_create_info_struct (png);
	if (setjmp (png_jmpbuf(png))) {
		png_destroy_read_struct (&png, &info, &end_info);
		fclose (file);
		return NULL;
    }
	//set up png read
    png_set_read_fn (png, file, png_read_data);
	png_set_sig_bytes (png,8);

	png_read_info(png,info);

	unsigned long width,height;
	int bit_depth,color_type,interlace_type,compression_type,filter_method;
	png_get_IHDR(png,info,&width,&height,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);

	//converting to 32bpp format
	if (color_type==PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
	
	if (color_type==PNG_COLOR_TYPE_GRAY || color_type==PNG_COLOR_TYPE_GRAY_ALPHA)
          png_set_gray_to_rgb(png);

	if (png_get_valid(png,info,PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);

	if (interlace_type!=PNG_INTERLACE_NONE) {
    	  png_destroy_read_struct (&png, &info, &end_info);
		  fclose (file);
          return NULL;
    }

    if (color_type & PNG_COLOR_MASK_ALPHA)
        png_set_strip_alpha(png);

	//setting gamma
	double gamma;
	if (png_get_gAMA(png,info,&gamma))
		png_set_gamma(png, 2.0, gamma);
	else
		png_set_gamma(png,2.0, 0.45455);

    int bps = 8;

	//updating png info struct
	png_read_update_info(png,info);
	png_get_IHDR(png,info,&width,&height,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);

        if (color_type & PNG_COLOR_MASK_ALPHA)
          png_set_strip_alpha(png);

	png_read_update_info(png,info);
	png_get_IHDR(png,info,&width,&height,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);


    unsigned char* data = new unsigned char[width*height*3];
    int rowlen = width*3;
    unsigned char *row = new unsigned char [rowlen];

	for (unsigned int i=0;i<height;i++) {

  	    png_read_row (png, (png_byte*)row, NULL);
        memcpy (data+3*i*width, row, 3*width);
    }

	png_read_end (png, 0);
	png_destroy_read_struct (&png, &info, &end_info);
	
	delete [] row;
	fclose(file);
    
    w = width;
    h = height;
    return data;
}

void savePNG  (char* fname, unsigned char* data1, unsigned char* data2, int w, int h) {

    FILE* file = fopen(fname,"wb");

	png_structp png = png_create_write_struct (PNG_LIBPNG_VER_STRING,0,0,0);
	png_infop info = png_create_info_struct(png);

	if (setjmp(png_jmpbuf(png))) {
		png_destroy_write_struct (&png,&info);
		fclose(file);
		return;
    }

	png_set_write_fn (png, file, png_write_data, png_flush);	
	png_set_compression_level(png,6);

    int width = w;
    int height = h;
    int bps = 8;

	png_set_IHDR(png, info, width, height, bps, PNG_COLOR_TYPE_RGB_ALPHA,
		PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_BASE);


    int rowlen = width*4;
    unsigned char *row = new unsigned char [rowlen];

	png_write_info(png,info);
	for (unsigned int i=0;i<height;i++) {
        for (int j=0; j<width; j++) {
            int ofs = 3*(width*i+j);
            unsigned char alpha = data2[ofs] - data1[ofs];
            if (i==8 && j==8)
                printf ("alpha=%d pix=%d\n",(int)alpha, (int)(data1[ofs+0] / (1.0 - alpha/255.0)));
            if (alpha<255) {
                row[4*j+0] = data1[ofs+0] / (1.0 - alpha/255.0);
                row[4*j+1] = data1[ofs+1] / (1.0 - alpha/255.0);
                row[4*j+2] = data1[ofs+2] / (1.0 - alpha/255.0);
            }
            else {
            
            }
            row[4*j+3] = 255-alpha;
        }
        png_write_row (png, (png_byte*)row);
    }

	png_write_end(png,info);
	png_destroy_write_struct(&png,&info);

    delete [] row;
	fclose (file);
}

int main (int argc, char* argv[]) {

        int w, h;
        unsigned char* data1 = loadPNG (argv[1], w, h);
        unsigned char* data2 = loadPNG (argv[2], w, h);
        savePNG (argv[3], data1, data2, w, h);
        
}
