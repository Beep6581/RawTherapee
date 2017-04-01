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

#include "previewimage.h"
#include "iimage.h"
#include "utils.h"
#include "iimage.h"
#include "rtthumbnail.h"
#include "rawimagesource.h"

using namespace rtengine;
using namespace procparams;

PreviewImage::PreviewImage (const Glib::ustring &fname, const Glib::ustring &ext, const PreviewImageMode mode)
{
    rtengine::Thumbnail* tpp = nullptr;

    if (mode == PIM_EmbeddedPreviewOnly || mode == PIM_EmbeddedOrRaw) {

        const unsigned char *data = nullptr;

        int width = -1, height = -1;

        if (ext.lowercase() == "jpg" || ext.lowercase() == "jpeg") {
            // int deg = infoFromImage (fname);
            tpp = rtengine::Thumbnail::loadFromImage (fname, width, height, 1, 1., true);

            if (tpp) {
                data = tpp->getImage8Data();
            }
        } else if (ext.lowercase() == "png") {
            tpp = rtengine::Thumbnail::loadFromImage (fname, width, height, 1, 1., true);

            if (tpp) {
                data = tpp->getImage8Data();
            }
        } else if (ext.lowercase() == "tif" || ext.lowercase() == "tiff") {
            // int deg = infoFromImage (fname);
            tpp = rtengine::Thumbnail::loadFromImage (fname, width, height, 1, 1., true);

            if (tpp) {
                data = tpp->getImage8Data();
            }
        } else {
            rtengine::RawMetaDataLocation ri;
            tpp = rtengine::Thumbnail::loadQuickFromRaw (fname, ri, width, height, 1, true, true);

            if (tpp) {
                data = tpp->getImage8Data();
            }
        }

        if (tpp) {
            if (data) {
                int w, h;
                double scale = 1.;

                tpp->getDimensions(w, h, scale);

                previewImage = Cairo::ImageSurface::create(Cairo::FORMAT_RGB24, w, h);
                previewImage->flush();

                #pragma omp parallel
                {
                    const unsigned char *src;
                    unsigned char *dst;
                    #pragma omp for schedule(static,10)

                    for (unsigned int i = 0; i < (unsigned int)(h); ++i) {
                        src = data + i * w * 3;
                        dst = previewImage->get_data() + i * w * 4;

                        for (unsigned int j = 0; j < (unsigned int)(w); ++j) {
                            unsigned char r = *(src++);
                            unsigned char g = *(src++);
                            unsigned char b = *(src++);

                            poke255_uc(dst, r, g, b);
                        }
                    }
                }
                previewImage->mark_dirty();
            }
        }
    }

    if ((mode == PIM_EmbeddedOrRaw && !tpp) || mode == PIM_ForceRaw) {
        RawImageSource rawImage;
        int error = rawImage.load(fname, true);

        if (!error) {
            rtengine::Image8 *output = nullptr;
            const unsigned char *data = nullptr;
            int fw, fh;

            procparams::ProcParams params;
            ColorTemp wb = rawImage.getWB ();
            rawImage.getFullSize (fw, fh, TR_NONE);
            PreviewProps pp (0, 0, fw, fh, 1);
            params.icm.input = Glib::ustring("(embedded)");
            params.raw.bayersensor.method = RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::fast];
            params.raw.deadPixelFilter = false;
            params.raw.ca_autocorrect = false;
            params.raw.xtranssensor.method = RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::fast];
            rawImage.preprocess(params.raw, params.lensProf, params.coarse);
            rawImage.demosaic(params.raw);
            Imagefloat* image = new rtengine::Imagefloat (fw, fh);
            rawImage.getImage (wb, TR_NONE, image, pp, params.toneCurve, params.icm, params.raw);
            output = new Image8(fw, fh);
            rawImage.convertColorSpace(image, params.icm, wb);
            #pragma omp parallel for schedule(dynamic, 10)
            for (int i = 0; i < fh; ++i)
                for (int j = 0; j < fw; ++j) {
                    image->r(i, j) = Color::gamma2curve[image->r(i, j)];
                    image->g(i, j) = Color::gamma2curve[image->g(i, j)];
                    image->b(i, j) = Color::gamma2curve[image->b(i, j)];
                }


            image->resizeImgTo<Image8>(fw, fh, TI_Nearest, output);
            data = output->getData();


            if (data) {
                int w, h;
             //   double scale = 1.;
                w = output->getWidth();
                h = output->getHeight();
                previewImage = Cairo::ImageSurface::create(Cairo::FORMAT_RGB24, w, h);
                previewImage->flush();

                #pragma omp parallel
                {
                    const unsigned char *src;
                    unsigned char *dst;
                    #pragma omp for schedule(static,10)

                    for (unsigned int i = 0; i < (unsigned int)(h); i++) {
                        src = data + i * w * 3;
                        dst = previewImage->get_data() + i * w * 4;

                        for (unsigned int j = 0; j < (unsigned int)(w); j++) {
                            unsigned char r = *(src++);
                            unsigned char g = *(src++);
                            unsigned char b = *(src++);

                            poke255_uc(dst, r, g, b);
                        }
                    }
                }

                if (output) {
                    delete output;
                }

                previewImage->mark_dirty();
            }
        }
    }

    if (tpp) {
        delete tpp;
    }

}

Cairo::RefPtr<Cairo::ImageSurface> PreviewImage::getImage()
{
    return previewImage;
}
