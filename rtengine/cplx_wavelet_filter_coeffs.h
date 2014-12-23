/*
 *  This file is part of RawTherapee.
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
 *
 *  2012 Emil Martinec <ejmartin@uchicago.edu>
 */


namespace rtengine {

const int Daub4_len=6;
const int Daub4_offset=2;
const float Daub4_anal[2][6] ALIGNED16 = {//analysis filter
	{0, 0, 0.34150635, 0.59150635, 0.15849365, -0.091506351}, 
	{-0.091506351, -0.15849365, 0.59150635, -0.34150635, 0, 0}};

};

