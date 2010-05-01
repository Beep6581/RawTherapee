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
#define MINMAX3(a,b,c,min,max) \
{ \
if ((a)<(b)) { \
  if ((b)<(c)) { \
    (min) = (a); \
    (max) = (c); \
  } \
  else { \
    (max) = (b); \
    if ((a)<(c)) \
      (min) = (a); \
    else \
      (min) = (c); \
  } \
} else { \
  if ((b)>(c)) { \
    (min) = (c); \
    (max) = (a); \
  } \
  else { \
    (min) = (b); \
    if ((a)>(c)) \
      (max) = (a); \
    else \
      (max) = (c); \
  } \
} \
}

#define MIN3(a,b,c,min) \
{ \
if ((a)<(b)) { \
  if ((a)<(c))  \
    (min) = (a); \
  else  \
    (min) = (c); \
} else { \
  if ((b)>(c))  \
    (min) = (c); \
  else \
    (min) = (b); \
} \
} 

#define MAX3(a,b,c,min) \
{ \
if ((a)>(b)) { \
  if ((a)>(c))  \
    (max) = (a); \
  else  \
    (max) = (c); \
} else { \
  if ((b)<(c))  \
    (max) = (c); \
  else \
    (max) = (b); \
} \
} 
