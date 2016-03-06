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
#include "rt_math.h"

#define SORT3(a1,a2,a3,b1,b2,b3) \
  {  \
  b2 = min(a1,a2);\
  b1 = min(b2,a3);\
  b3 = max(a1,a2);\
  b2 = max(b2, min(b3,a3));\
  b3 = max(b3,a3);\
 }


#define NETWORKSORT4OF6(s0,s1,s2,s3,s4,s5,d0,d1,d2,d3,d4,d5,temp) \
{\
d1 = min(s1,s2);\
d2 = max(s1,s2);\
d0 = min(s0,d2);\
d2 = max(s0,d2);\
temp = min(d0,d1);\
d1 = max(d0,d1);\
d0 = temp;\
d4 = min(s4,s5);\
d5 = max(s4,s5);\
temp = min(s3,d5);\
d5 = max(s3,d5);\
d3 = temp;\
temp = min(d3,d4);\
d4 = max(d3,d4);\
d3 = temp;\
d3 = max(d0,d3);\
temp = min(d1,d4);\
d4 = max(d1,d4);\
d1 = temp;\
d2 = min(d2,d5);\
temp = min(d2,d4);\
d4 = max(d2,d4);\
d2 = temp;\
temp = min(d1,d3);\
d3 = max(d1,d3);\
d1 = temp;\
temp = min(d2,d3);\
d3 = max(d2,d3);\
d2 = temp;\
}

#define MEDIAN7(s0,s1,s2,s3,s4,s5,s6,t0,t1,t2,t3,t4,t5,t6,median) \
{\
t0 = min(s0,s5);\
t5 = max(s0,s5);\
t3 = max(t0,s3);\
t0 = min(t0,s3);\
t1 = min(s1,s6);\
t6 = max(s1,s6);\
t2 = min(s2,s4);\
t4 = max(s2,s4);\
t1 = max(t0,t1);\
median = min(t3,t5);\
t5 = max(t3,t5);\
t3 = median;\
median = min(t2,t6);\
t6 = max(t2,t6);\
t2 = median;\
t3 = max(t2,t3);\
t3 = min(t3,t6);\
t4 = min(t4,t5);\
median = min(t1,t4);\
t4 = max(t1,t4);\
t1 = median;\
t3 = max(t1,t3);\
median = min(t3,t4);\
}
