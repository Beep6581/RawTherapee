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
 *  2014 Jacques Desmis <jdesmis@gmail.com>
 */


namespace rtengine {

const int Daub4_offset=2;

const float Daub4_anal0[2][4] ALIGNED16 = {//analysis filter 2  Hall
		{0.f, 0.f, 0.5f, 0.5f}, 
		{-0.5f, 0.5f,  0.f, 0.f}};

const float Daub4_anal[2][6] ALIGNED16 = {//Daub4
		{0.f, 0.f, 0.34150635f, 0.59150635f, 0.15849365f, -0.091506351f}, 
		{-0.091506351f, -0.15849365f, 0.59150635f, -0.34150635f, 0.f, 0.f}};

const float Daub4_anal8[2][8] ALIGNED16 = {//Daub6
		{0.f, 0.f, 0.235233605f, 0.57055846f, 0.3251825f, -0.09546721f, -0.060416105f, 0.02490875f}, 
		{-0.02490875f,  -0.060416105f, 0.09546721f, 0.3251825f, -0.57055846f , 0.235233605f, 0.f, 0.f}};
		
const float Daub4_anal12[2][12] ALIGNED16 = {//Daub10
		{0.f, 0.f, 0.11320949f, 0.42697177f, 0.51216347f, 0.09788348f, -0.171328355f, -0.022800565f, 0.054851325f, -0.0044134f, -0.008895935f, 0.002358714f}, 
		{-0.002358714f,  -0.008895935f, 0.0044134f, 0.054851325f, 0.022800565f , -0.171328355f, -0.09788348f, 0.51216347f, -0.42697177f, 0.11320949f, 0.f, 0.f}};
		
const float Daub4_anal16[2][16] ALIGNED16 = {//Daub 14	
		{0.f, 0.f, 0.055049715f, 0.28039564f, 0.515574245f, 0.33218624f, -0.10175691f, -0.158417505f, 0.05042335f, 0.057001725f, -0.026891225f, -0.01171997f, 0.008874895f, 0.0003037575f, -0.0012739524f, 0.0002501134f}, 
		{-0.0002501134f,  -0.0012739524f, -0.0003037575f, 0.008874895f, 0.01171997f , -0.026891225f, -0.057001725f, 0.05042335f, 0.158417505f, -0.10175691f, -0.33218624f, 0.515574245f, -0.28039564f, 0.055049715f, 0.f, 0.f}};

// if necessary ?? we can add D20 !!		
};

