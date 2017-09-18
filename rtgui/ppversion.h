#ifndef _PPVERSION_
#define _PPVERSION_

// This number has to be incremented whenever the PP3 file format is modified or the behaviour of a tool changes
#define PPVERSION 327
#define PPVERSION_AEXP 301 //value of PPVERSION when auto exposure algorithm was modified

/*
  Log of version changes
   327  2017-09-15
        [Profiles Lens Correction] Added Lensfun
   326  2015-07-26
        [Exposure] Added 'Perceptual' tone curve mode
   325  2015-07-23
        [Exposure] [RGB Curves] [B&W] Normalized RGB pipeline curve gammas to sRGB (before it was a mix between sRGB and 1.0 and depended on file format)
   323  2015-10-05
        [Exposure] Added 'Luminance' tone curve mode
   322  2015-01-31
        [Wavelet] new tool using wavelet levels
   321  2014-08-17
        [Film Simulation] new  tool using HALDCLUT files
   320  2014-07-02  (yes, same version number... this is an error due to a wrong version number set in comment of previous change)
        New [RAW Bayer] and [RAW X-Trans] sections, with some parameters transfered from [RAW] to [RAW Bayer]
   320  2014-03-29
        [ColorToning] new tool for color toning
   319  2014-02-11
        Hue skin for Contrast by detail levels
   318  2014-02-10
        Vignetting Correction bug makes hard transitions for positive Amount values, Issue 2241
   317  2014-01-19
        changes to behaviour of LC curve, Issue 2209
   315  2013-12-12
        add LH et HH curve to lab mode
   313  2013-11-19
        add CL curve to lab mode
   312  2013-11-08
        added numerous changes to [channel mixer]
   311  2013-11-07
        [Gradient] new tool (gradient/graduated filter
        [PCVignette] new tool (vignette filter)
   310  2013-09-16
        Defringing /Threshold - changed calculation, issue 1801
   307  2013-03-16
        [Perspective] Horizontal and Vertical changed from int to double
        added  [Directional Pyramid Denoising] Method, Redchro, Bluechro
        added [RGB Curves] LumaMode
 */

#endif
