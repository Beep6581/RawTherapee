#ifndef _PPVERSION_
#define _PPVERSION_

// This number have to be incremented whenever the PP3 file format is modified
#define PPVERSION 312
#define PPVERSION_AEXP 301 //value of PPVERSION when auto exposure algorithm was modified

/*
  Log of version changes
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
