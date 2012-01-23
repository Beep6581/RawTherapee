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
 *  2010 Ilya Popov <ilia_popov@rambler.ru>
 *  2012 Emil Martinec <ejmartin@uchicago.edu>
 */

#ifndef CPLX_WAVELET_DEC_H_INCLUDED
#define CPLX_WAVELET_DEC_H_INCLUDED

#include <cstddef>
#include <math.h>


#include "cplx_wavelet_level.h"
#include "cplx_wavelet_filter_coeffs.h"

namespace rtengine {


// %%%%%%%%%%%%%%%%%%%%%%%%%%%

template <typename A, typename B>
void copy_out(A ** a, B * b, size_t datalen)
{
	for (size_t j=0; j<datalen; j++) {
		b[j] = static_cast<B> (0.25f*(a[0][j]+a[1][j]+a[2][j]+a[3][j]));
	}
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%


class cplx_wavelet_decomposition
{
public:

    typedef float internal_type;

private:

    static const int maxlevels = 8;//should be greater than any conceivable order of decimation
    
    int lvltot;
    size_t m_w, m_h;//dimensions
    size_t m_w1, m_h1;
	
	int first_lev_len, first_lev_offset;
	//multi_array2D<float,2> first_lev_anal;
	//multi_array2D<float,2> first_lev_synth;
	float *first_lev_anal;
	float *first_lev_synth;
	
	int wavfilt_len, wavfilt_offset;
	//multi_array2D<float,2> wavfilt_anal;
	//multi_array2D<float,2> wavfilt_synth;
	float *wavfilt_anal;
	float *wavfilt_synth;

    cplx_wavelet_level<internal_type> * dual_tree_coeffs[maxlevels][4];//m_c in old code
    
public:

    template<typename E>
    cplx_wavelet_decomposition(E * src, int width, int height, int maxlvl);
    
    ~cplx_wavelet_decomposition();
    
    template<typename E>
    void reconstruct(E * dst);
	
		
};

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
/*template<typename E>
cplx_wavelet_decomposition::cplx_wavelet_decomposition(E * src, int width, int height, int maxlvl)
: lvltot(0), m_w(w), m_h(h), m_w1(0), m_h1(0)
{
    m_w1 = w;
    m_h1 = h;
    
    m_c[0] = new cplx_wavelet_level<internal_type>(src, m_w1, m_h1, FSFarras);
    lvltot = 1;
    
    while(lvltot < maxlevels)
    {
        m_c[level] = new cplx_wavelet_level<internal_type>(m_c[lvltot-1]->data[0], m_c[lvltot-1]->width(), 
														   m_c[lvltot-1]->height(), Kingsbury);
        lvltot ++;
    }
}*/


/*template<typename E, typename L>
void cplx_wavelet_decomposition::reconstruct(E * dst)
{
    noop<internal_type> n;

    for(int level = lvltot - 1; level > 0; level--)
    {
        int alpha = 1024 + 10 * c[level];
        m_c[level]->reconstruct(m_c[level-1]->lowfreq(), alpha, n);
    }
    
    int alpha = 1024 + 10 * c[0];
    m_c[0]->reconstruct(dst, alpha, l);
}*/
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename E>
	cplx_wavelet_decomposition::cplx_wavelet_decomposition(E * src, int width, int height, int maxlvl)
	: lvltot(0), m_w(width), m_h(height), m_w1(0), m_h1(0)
	{
		m_w1 = width;
		m_h1 = height;
		
		//initialize wavelet filters
		
		first_lev_len = FSFarras_len;
		first_lev_offset = FSFarras_offset;
		//multi_array2D<float,2> first_lev_anal(2,first_lev_len);
		//multi_array2D<float,2> first_lev_synth(2,first_lev_len);
		float *first_level_anal = new float[4*first_lev_len];
		float *first_level_synth = new float[4*first_lev_len];

		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				for (int i=0; i<first_lev_len; i++) {
					//first_lev_anal[n][m][i]  = FSFarras_anal[n][m][i];
					//first_lev_synth[n][m][i] = FSFarras_anal[n][m][first_lev_len-1-i];
					first_lev_anal[first_lev_len*(2*n+m)+i]  = FSFarras_anal[n][m][i];
					first_lev_synth[first_lev_len*(2*n+m)+i] = FSFarras_anal[n][m][first_lev_len-1-i];
				}
			}
		}
		
		wavfilt_len = Kingsbury_len;
		wavfilt_offset = Kingsbury_offset;
		//multi_array2D<float,2> wavfilt_anal(2,Kingsbury_len);
		//multi_array2D<float,2> wavfilt_synth(2,Kingsbury_len);
		float *wavfilt_anal = new float[4*wavfilt_len];
		float *wavfilt_synth = new float[4*wavfilt_len];
		
		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				for (int i=0; i<wavfilt_len; i++) {
					//wavfilt_anal[n][m][i]  = Kingsbury_anal[n][m][i];
					//wavfilt_synth[n][m][i] = Kingsbury_anal[n][m][wavfilt_len-1-i];
					wavfilt_anal[wavfilt_len*(2*n+m)+i]  = Kingsbury_anal[n][m][i];
					wavfilt_synth[wavfilt_len*(2*n+m)+i] = Kingsbury_anal[n][m][first_lev_len-1-i];
				}
			}
		}
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
		// Initialize wavelet coeffs
		
		/*CplxWavelet AntonB = {
			12,//length of filter
			6,//offset
			
			{//analysis filter
				{{0, -0.08838834764832, 0.08838834764832, 0.69587998903400, 0.69587998903400, 
					0.08838834764832, -0.08838834764832, 0.01122679215254, 0.01122679215254, 0},
					{0, 0, 0, 0.04563588155712, -0.02877176311425, -0.29563588155712 , 
						0.55754352622850, -0.29563588155713, -0.02877176311425, 0.04563588155712, 0, 0}},
				{{0 , 0 , 0.02674875741081, -0.01686411844287, -0.07822326652899, 0.26686411844288, 
					0.60294901823636, 0.26686411844287, -0.07822326652899, -0.01686411844287, 0.02674875741081, 0},
					{0 , 0 , 0, 0 , 0.04563588155712, -0.02877176311425, 
						-0.29563588155712 , 0.55754352622850, -0.29563588155713, -0.02877176311425, 0.04563588155712 , 0}} },
			
			{//synthesis filter
				{{0 , 0 , 0, -0.04563588155712, -0.02877176311425, 0.29563588155712, 
					0.55754352622850, 0.29563588155713, -0.02877176311425, -0.04563588155712, 0, 0},
					{0, 0.02674875741081, 0.01686411844287, -0.07822326652899, -0.26686411844288 , 0.60294901823636, 
						-0.26686411844287, -0.07822326652899, 0.01686411844287, 0.02674875741081, 0, 0}},
				{{0 , 0, -0.04563588155712, -0.02877176311425, 0.29563588155712 , 0.55754352622850 , 
					0.29563588155713, -0.02877176311425, -0.04563588155712, 0, 0 , 0},
					{0.02674875741081 , 0.01686411844287, -0.07822326652899, -0.26686411844288 , 0.60294901823636, -0.26686411844287, 
						-0.07822326652899, 0.01686411844287 , 0.02674875741081 , 0 , 0, 0}} }
		};*/
		
		/*for (int i=0; i<4; i++)
			for (int n=0; n<12; n++) {
				AntonB.synth[i][n] *= 2;
			}*/
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
		// Initialize wavelet coeffs
		
		//CplxWavelet FSFarras = {
		/*int FSFarras_len = 10;//length of filter
		int FSFarras_offset = 5;//offset
			
			float FSFarras_anal[2][2][10] = {//analysis filter
				{{0, -0.08838834764832, 0.08838834764832, 0.69587998903400, 0.69587998903400, 0.08838834764832, -0.08838834764832, 0.01122679215254 , 0.01122679215254, 0},
					{ 0, -0.01122679215254, 0.01122679215254, 0.08838834764832, 0.08838834764832, -0.69587998903400, 0.69587998903400, -0.08838834764832, -0.08838834764832, 0}},
				{{0.01122679215254, 0.01122679215254, -0.08838834764832, 0.08838834764832, 0.69587998903400, 0.69587998903400, 0.08838834764832, -0.08838834764832, 0, 0},
					{0, 0, -0.08838834764832, -0.08838834764832, 0.69587998903400, -0.69587998903400, 0.08838834764832, 0.08838834764832, 0.01122679215254, -0.01122679215254}} };
		float FSFarras_synth[2][2][10];*/
		//};
		
		/*for (int i=0; i<4; i++)
			for (int n=0; n<10; n++) {
				FSFarras_synth[i][n] = FSFarras_anal[i][9-n];
			}*/
		
		//sf = Reverse[af, 3];
		
		
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
		// Initialize wavelet coeffs

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
		
		//CplxWavelet Kingsbury {
		/*int Kingsbury_len = 10;//length of filter
		int Kingsbury_offset = 5;//offset
			
			float Kingsbury_anal[2][2][10] = {//analysis filter
				{{0.03516384000000, 0, -0.08832942000000, 0.23389032000000, 0.76027237000000, 0.58751830000000, 0, -0.11430184000000 , 0, 0},
					{ 0, 0, -0.11430184000000, 0, 0.58751830000000, -0.76027237000000, 0.23389032000000, 0.08832942000000, 0, -0.03516384000000}},
				{{0, 0, -0.11430184000000, 0, 0.58751830000000, 0.76027237000000, 0.23389032000000, -0.08832942000000, 0, 0.03516384000000},
					{-0.03516384000000, 0, 0.08832942000000, 0.23389032000000, -0.76027237000000, 0.58751830000000, 0, -0.11430184000000, 0, 0}} };
		float Kingsbury_synth[2][2][10];*/

		//};
		
		/*for (int i=0; i<4; i++)
			for (int n=0; n<10; n++) {
				Kingsbury_synth[i][n] = Kingsbury_anal[i][9-n];
			}*/
		
		//sf = Reverse[af, 3];
		
		
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
		
		
		// data structure is dual_tree_coeffs[scale][2*n+m=2*(Re/Im)+dir][channel={lo,hi1,hi2,hi3}][pixel_array]
		
		
		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				
				//dual_tree_coeffs[0][2*n+m] = new cplx_wavelet_level<internal_type>(src, first_lev_anal[n], first_lev_anal[m], first_lev_len, first_lev_offset);
				dual_tree_coeffs[0][2*n+m] = new cplx_wavelet_level<internal_type>(src, m_w, m_h, first_lev_anal+first_lev_len*2*n, \
																				   first_lev_anal+first_lev_len*2*m, first_lev_len, first_lev_offset);
				lvltot=1;
				while(lvltot < maxlevels) {
					//dual_tree_coeffs[lvltot][2*n+m] = new cplx_wavelet_level<internal_type>(dual_tree_coeffs[lvltot-1][2*n+m]->lopass()/*lopass*/, \
																							wavfilt_anal[n], wavfilt_anal[m], wavfilt_len, wavfilt_offset);
					dual_tree_coeffs[lvltot][2*n+m] = new cplx_wavelet_level<internal_type>(dual_tree_coeffs[lvltot-1][2*n+m]->lopass()/*lopass*/, \
																							dual_tree_coeffs[lvltot-1][2*n+m]->width(), \
																							dual_tree_coeffs[lvltot-1][2*n+m]->height(), \
																							wavfilt_anal+wavfilt_len*2*n, wavfilt_anal+wavfilt_len*2*m, wavfilt_len, wavfilt_offset);
					lvltot++;
				}
			}
		}
		
		
		//rotate detail coefficients
		float root2 = sqrt(2);
		for (int lvl=0; lvl<lvltot; lvl++) {
			int Wlvl = dual_tree_coeffs[lvl][0]->width();
			int Hlvl = dual_tree_coeffs[lvl][0]->height();
			for (int i=0; i<Wlvl*Hlvl; i++) {//pixel
				for (int m=1; m<4; m++) {//detail coefficients only
					float wavtmp = (dual_tree_coeffs[lvl][0]->wavcoeffs[m][i] + dual_tree_coeffs[lvl][3]->wavcoeffs[m][i])/root2;
					dual_tree_coeffs[lvl][3]->wavcoeffs[m][i] = (dual_tree_coeffs[lvl][0]->wavcoeffs[m][i] - dual_tree_coeffs[lvl][3]->wavcoeffs[m][i])/root2;
					dual_tree_coeffs[lvl][0]->wavcoeffs[m][i] = wavtmp;
					
					wavtmp = (dual_tree_coeffs[lvl][1]->wavcoeffs[m][i] + dual_tree_coeffs[lvl][2]->wavcoeffs[m][i])/root2;
					dual_tree_coeffs[lvl][2]->wavcoeffs[m][i] = (dual_tree_coeffs[lvl][1]->wavcoeffs[m][i] - dual_tree_coeffs[lvl][2]->wavcoeffs[m][i])/root2;
					dual_tree_coeffs[lvl][1]->wavcoeffs[m][i] = wavtmp;
				}
			}
		}
		
	}
	
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	
	
	/* function y=reconstruct(w,J,Fsf,sf) */
	
	template<typename E>
	void cplx_wavelet_decomposition::reconstruct(E * dst) { 
		
		// data structure is wavcoeffs[scale][2*n+m=2*(Re/Im)+dir][channel={lo,hi1,hi2,hi3}][pixel_array]
		
		//rotate detail coefficients
		float root2 = sqrt(2);
		for (int lvl=0; lvl<lvltot; lvl++) {
			int Wlvl = dual_tree_coeffs[lvl][0]->width();
			int Hlvl = dual_tree_coeffs[lvl][0]->height();
			for (int i=0; i<Wlvl*Hlvl; i++) {//pixel
				for (int m=1; m<4; m++) {//detail coefficients only
					float wavtmp = (dual_tree_coeffs[lvl][0]->wavcoeffs[m][i] + dual_tree_coeffs[lvl][3]->wavcoeffs[m][i])/root2;
					dual_tree_coeffs[lvl][3]->wavcoeffs[m][i] = (dual_tree_coeffs[lvl][0]->wavcoeffs[m][i] - dual_tree_coeffs[lvl][3]->wavcoeffs[m][i])/root2;
					dual_tree_coeffs[lvl][0]->wavcoeffs[m][i] = wavtmp;
					
					wavtmp = (dual_tree_coeffs[lvl][1]->wavcoeffs[m][i] + dual_tree_coeffs[lvl][2]->wavcoeffs[m][i])/root2;
					dual_tree_coeffs[lvl][2]->wavcoeffs[m][i] = (dual_tree_coeffs[lvl][1]->wavcoeffs[m][i] - dual_tree_coeffs[lvl][2]->wavcoeffs[m][i])/root2;
					dual_tree_coeffs[lvl][1]->wavcoeffs[m][i] = wavtmp;
				}
			}
		}
		
		//y = ConstantArray[0, {vsizetmp, hsizetmp}];
		array2D<internal_type> tmp(4,m_w*m_h);
		
		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				for (int lvl=lvltot-1; lvl>0; lvl--) {
					//m_c[level]->reconstruct(m_c[level-1]->lowfreq(), alpha, n);
					//dual_tree_coeffs[lvl][2*n+m]->reconstruct_level(dual_tree_coeffs[lvl-1][2*n+m]->wavcoeffs[0], wavfilt_synth[n], wavfilt_synth[m], wavfilt_len, wavfilt_offset);
					dual_tree_coeffs[lvl][2*n+m]->reconstruct_level(dual_tree_coeffs[lvl-1][2*n+m]->wavcoeffs[0], wavfilt_synth+wavfilt_len*2*n, \
																	wavfilt_synth+wavfilt_len*2*m, wavfilt_len, wavfilt_offset);
				}
				//dual_tree_coeffs[0][2*n+m]->reconstruct_level(tmp[2*n+m], first_lev_synth[n], first_lev_synth[m], first_lev_len, first_lev_offset);
				dual_tree_coeffs[0][2*n+m]->reconstruct_level(tmp[2*n+m], first_lev_synth+wavfilt_len*2*n, first_lev_synth+wavfilt_len*2*m, first_lev_len, first_lev_offset);
			}
		}
		
		copy_out(tmp,dst,m_w*m_h);
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


};

#endif
