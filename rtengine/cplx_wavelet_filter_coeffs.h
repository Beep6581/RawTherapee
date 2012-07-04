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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//#include "improcfun.h"

#ifdef _OPENMP
#include <omp.h>
#endif


namespace rtengine {


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	
/*const int AntonB_len = 12;//length of filter
const int AntonB_offset = 6;//offset
	
const float AntonB_anal[2][2][12] = {//analysis filter
		{{0, -0.08838834764832, 0.08838834764832, 0.69587998903400, 0.69587998903400, 
			0.08838834764832, -0.08838834764832, 0.01122679215254, 0.01122679215254, 0},
		{0, 0, 0, 0.04563588155712, -0.02877176311425, -0.29563588155712 , 
			0.55754352622850, -0.29563588155713, -0.02877176311425, 0.04563588155712, 0, 0}},
		{{0 , 0 , 0.02674875741081, -0.01686411844287, -0.07822326652899, 0.26686411844288, 
			0.60294901823636, 0.26686411844287, -0.07822326652899, -0.01686411844287, 0.02674875741081, 0},
		{0 , 0 , 0, 0 , 0.04563588155712, -0.02877176311425, 
			-0.29563588155712 , 0.55754352622850, -0.29563588155713, -0.02877176311425, 0.04563588155712 , 0}} };
	
const float AntonB_synth[2][2][12] = {//synthesis filter
		{{0 , 0 , 0, -0.04563588155712, -0.02877176311425, 0.29563588155712, 
			0.55754352622850, 0.29563588155713, -0.02877176311425, -0.04563588155712, 0, 0},
		{0, 0.02674875741081, 0.01686411844287, -0.07822326652899, -0.26686411844288 , 0.60294901823636, 
			-0.26686411844287, -0.07822326652899, 0.01686411844287, 0.02674875741081, 0, 0}},
		{{0 , 0, -0.04563588155712, -0.02877176311425, 0.29563588155712 , 0.55754352622850 , 
			0.29563588155713, -0.02877176311425, -0.04563588155712, 0, 0 , 0},
		{0.02674875741081 , 0.01686411844287, -0.07822326652899, -0.26686411844288 , 0.60294901823636, -0.26686411844287, 
			-0.07822326652899, 0.01686411844287 , 0.02674875741081 , 0 , 0, 0}} };*/

	
/*	for (int i=0; i<4; i++)
		for (int n=0; n<12; n++) {
			AntonB.synth[i][n] *= 2;
		}*/

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	


const int FSFarras_len=10;//length of filter
const int FSFarras_offset=4;//offset

const float FSFarras_anal[2][2][10] = {//analysis filter
	{{0, -0.08838834764832, 0.08838834764832, 0.69587998903400, 0.69587998903400, 0.08838834764832, -0.08838834764832, 0.01122679215254 , 0.01122679215254, 0},
		{ 0, -0.01122679215254, 0.01122679215254, 0.08838834764832, 0.08838834764832, -0.69587998903400, 0.69587998903400, -0.08838834764832, -0.08838834764832, 0}},
	{{0.01122679215254, 0.01122679215254, -0.08838834764832, 0.08838834764832, 0.69587998903400, 0.69587998903400, 0.08838834764832, -0.08838834764832, 0, 0},
		{0, 0, -0.08838834764832, -0.08838834764832, 0.69587998903400, -0.69587998903400, 0.08838834764832, 0.08838834764832, 0.01122679215254, -0.01122679215254}} };

//synthesis filter is the reverse (see cplx_wavelet_dec.h)


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

/*
 % Kingsbury Q-filters for the dual-tree complex DWT
 %
 % af{i},i=1,2-analysis filters for tree i
 % sf{i},i=1,2-synthesis filters for tree i
 % note:af{2} is the reverse of af{1}
 % ordering is {af[1],af[2],sf[1],sf[2]}
 % REFERENCE:% N.G.Kingsbury,"A dual-tree complex wavelet
 %    transform with improved orthogonality and symmetry
 %    properties",Proceedings of the IEEE Int.Conf.on
 % Image Proc.(ICIP),2000  */

const int Kingsbury_len=10;//length of filter
const int Kingsbury_offset=4;//offset

const float Kingsbury_anal[2][2][10] =	{//analysis filter
	{{0.03516384000000, 0, -0.08832942000000, 0.23389032000000, 0.76027237000000, 0.58751830000000, 0, -0.11430184000000 , 0, 0},
		{ 0, 0, -0.11430184000000, 0, 0.58751830000000, -0.76027237000000, 0.23389032000000, 0.08832942000000, 0, -0.03516384000000}},
	{{0, 0, -0.11430184000000, 0, 0.58751830000000, 0.76027237000000, 0.23389032000000, -0.08832942000000, 0, 0.03516384000000},
		{-0.03516384000000, 0, 0.08832942000000, 0.23389032000000, -0.76027237000000, 0.58751830000000, 0, -0.11430184000000, 0, 0}} };

//synthesis filter is the reverse (see cplx_wavelet_dec.h)

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

const int Haar_len=2;//length of filter
const int Haar_offset=1;//offset

const float Haar_anal[2][2] =	{{0.5,0.5}, {0.5,-0.5}};//analysis filter

//synthesis filter is the reverse (see cplx_wavelet_dec.h)

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

const int LeGall_len=6;
const int LeGall_offset=2;
const float LeGall_anal[2][6] = {{0, 0.25, 0.5, 0.25, 0, 0}, {0, -0.125, -0.25, 0.75, -0.25, -0.125}};
const float LeGall_synth[2][6] = {{-0.125, 0.25, 0.75, 0.25, -0.125, 0}, {0, 0, -0.25, 0.5, -0.25, 0}};
	
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

const int Daub4_len=6;
const int Daub4_offset=2;
const float Daub4_anal[2][6] = {//analysis filter
	{0, 0, 0.34150635, 0.59150635, 0.15849365, -0.091506351}, 
	{-0.091506351, -0.15849365, 0.59150635, -0.34150635, 0, 0}};

};

