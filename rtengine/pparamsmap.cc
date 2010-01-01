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
#include <rtengine.h>

int main () {

    printf ("Member offsets in 'ProcParams':\n");

    rtengine::procparams::ProcParams p;

    int d = (int) &p;

    printf ("%d\n", (int)&p.toneCurve.curve - d);
    printf ("%d\n", (int)&p.toneCurve.brightness - d);
    printf ("%d\n", (int)&p.toneCurve.black - d);
    printf ("%d\n", (int)&p.toneCurve.contrast - d);
    printf ("%d\n", (int)&p.toneCurve.shcompr - d);
    printf ("%d\n", (int)&p.toneCurve.hlcompr - d);
    printf ("%d\n", (int)&p.toneCurve.autoexp - d);
    printf ("%d\n", (int)&p.toneCurve.clip - d);
    printf ("%d\n", (int)&p.toneCurve.expcomp - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.lumaCurve.curve - d);
    printf ("%d\n", (int)&p.lumaCurve.brightness - d);
    printf ("%d\n", (int)&p.lumaCurve.black - d);
    printf ("%d\n", (int)&p.lumaCurve.contrast - d);
    printf ("%d\n", (int)&p.lumaCurve.shcompr - d);
    printf ("%d\n", (int)&p.lumaCurve.hlcompr - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.sharpening.enabled - d);
    printf ("%d\n", (int)&p.sharpening.radius - d);
    printf ("%d\n", (int)&p.sharpening.amount - d);
    printf ("%d\n", (int)&p.sharpening.threshold - d);
    printf ("%d\n", (int)&p.sharpening.edgesonly - d);
    printf ("%d\n", (int)&p.sharpening.edges_radius - d);
    printf ("%d\n", (int)&p.sharpening.edges_tolerance - d);
    printf ("%d\n", (int)&p.sharpening.halocontrol - d);
    printf ("%d\n", (int)&p.sharpening.halocontrol_amount - d);
    printf ("%d\n", (int)&p.sharpening.method - d);
    printf ("%d\n", (int)&p.sharpening.deconvamount - d);
    printf ("%d\n", (int)&p.sharpening.deconvradius - d);
    printf ("%d\n", (int)&p.sharpening.deconviter - d);
    printf ("%d\n", (int)&p.sharpening.deconvdamping - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.colorBoost.a - d);
    printf ("%d\n", (int)&p.colorBoost.b - d);
    printf ("%d\n", (int)&p.colorBoost.avoidclip - d);
    printf ("%d\n", (int)&p.colorBoost.enable_saturationlimiter - d);
    printf ("%d\n", (int)&p.colorBoost.saturationlimit - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.wb.method - d);
    printf ("%d\n", (int)&p.wb.temperature - d);
    printf ("%d\n", (int)&p.wb.green - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.colorShift.a - d);
    printf ("%d\n", (int)&p.colorShift.b - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.lumaDenoise.enabled - d);
    printf ("%d\n", (int)&p.lumaDenoise.radius - d);
    printf ("%d\n", (int)&p.lumaDenoise.edgetolerance - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.colorDenoise.enabled - d);
    printf ("%d\n", (int)&p.colorDenoise.radius - d);
    printf ("%d\n", (int)&p.colorDenoise.edgetolerance - d);
    printf ("%d\n", (int)&p.colorDenoise.edgesensitive - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.sh.enabled - d);
    printf ("%d\n", (int)&p.sh.highlights - d);
    printf ("%d\n", (int)&p.sh.htonalwidth - d);
    printf ("%d\n", (int)&p.sh.shadows - d);
    printf ("%d\n", (int)&p.sh.stonalwidth - d);
    printf ("%d\n", (int)&p.sh.localcontrast - d);
    printf ("%d\n", (int)&p.sh.radius - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.crop.enabled - d);
    printf ("%d\n", (int)&p.crop.x - d);
    printf ("%d\n", (int)&p.crop.y - d);
    printf ("%d\n", (int)&p.crop.w - d);
    printf ("%d\n", (int)&p.crop.h - d);
    printf ("%d\n", (int)&p.crop.fixratio - d);
    printf ("%d\n", (int)&p.crop.ratio - d);
    printf ("%d\n", (int)&p.crop.orientation - d);
    printf ("%d\n", (int)&p.crop.guide - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.coarse.rotate - d);
    printf ("%d\n", (int)&p.coarse.hflip - d);
    printf ("%d\n", (int)&p.coarse.vflip - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.rotate.degree - d);
    printf ("%d\n", (int)&p.rotate.fill - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.distortion.amount - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.cacorrection.red - d);
    printf ("%d\n", (int)&p.cacorrection.blue - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.vignetting.amount - d);
    printf ("%d\n", (int)&p.vignetting.radius - d);
    printf ("---------\n");
    printf ("%d\n",  (int)&p.chmixer.red[0] - d);
    printf ("%d\n",  (int)&p.chmixer.green[0] - d);
    printf ("%d\n",  (int)&p.chmixer.blue[0] - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.hlrecovery.enabled - d);
    printf ("%d\n", (int)&p.hlrecovery.method - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.resize.scale - d);
    printf ("%d\n", (int)&p.resize.method - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.icm.input - d);
    printf ("%d\n", (int)&p.icm.gammaOnInput - d);
    printf ("%d\n", (int)&p.icm.working - d);
    printf ("%d\n", (int)&p.icm.output - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.exif - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.iptc - d);
    printf ("---------\n");
    printf ("%d\n", (int)&p.version - d);
    printf ("---------\n");
    printf ("sizeof 'ProcParams' = %d\n", sizeof(rtengine::procparams::ProcParams));

    printf ("Member offsets in 'Settings':\n");
    rtengine::Settings s;

    d = (int) &s;

    printf ("%d\n", (int)&s.dualThreadEnabled - d);
    printf ("%d\n", (int)&s.demosaicMethod - d);
    printf ("%d\n", (int)&s.colorCorrectionSteps - d);
    printf ("%d\n", (int)&s.iccDirectory - d);
    printf ("%d\n", (int)&s.colorimetricIntent - d);
    printf ("%d\n", (int)&s.monitorProfile - d);
    printf ("---------\n");
    printf ("sizeof 'Settings' = %d\n", sizeof(rtengine::Settings));

    printf ("Member offsets in 'RawMetaDataLocation':\n");
    rtengine::RawMetaDataLocation r;

    d = (int) &r;

    printf ("%d\n", (int)&r.exifBase - d);
    printf ("%d\n", (int)&r.ciffBase - d);
    printf ("%d\n", (int)&r.ciffLength - d);
    printf ("---------\n");
    printf ("sizeof 'RawMetaDataLocation' = %d\n", sizeof(rtengine::RawMetaDataLocation));


}

