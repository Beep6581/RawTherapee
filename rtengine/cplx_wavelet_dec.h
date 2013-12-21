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
	
template <typename A>
void copy_out(A * a, A * b, size_t datalen)
{// for standard wavelet decomposition
	memcpy(b, a, datalen*sizeof(A));
}

template <typename A, typename B>
void copy_out(A ** a, B * b, size_t datalen)
{// for complex wavelet decomposition
	for (size_t j=0; j<datalen; j++) {
		b[j] = static_cast<B> (0.25*(a[0][j]+a[1][j]+a[2][j]+a[3][j]));
	}
}
	
template <typename A, typename B>
void copy_out(A * a, B * b, size_t datalen)
{// for standard wavelet decomposition
	for (size_t j=0; j<datalen; j++) {
		b[j] = static_cast<B> (a[j]);
	}
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%


class cplx_wavelet_decomposition
{
public:

    typedef float internal_type;

private:

  //  static const int maxlevels = 8;//should be greater than any conceivable order of decimation
    static const int maxlevels = 9;//should be greater than any conceivable order of decimation
    
    int lvltot, subsamp;
    size_t m_w, m_h;//dimensions
	
	int first_lev_len, first_lev_offset;
	float *first_lev_anal;
	float *first_lev_synth;
	
	int wavfilt_len, wavfilt_offset;
	float *wavfilt_anal;
	float *wavfilt_synth;
	
	int testfilt_len, testfilt_offset;
	float *testfilt_anal;
	float *testfilt_synth;

    wavelet_level<internal_type> * dual_tree[maxlevels][4];
    
public:

    template<typename E>
    cplx_wavelet_decomposition(E * src, int width, int height, int maxlvl, int subsampling);
    
    ~cplx_wavelet_decomposition();
	
	internal_type ** level_coeffs(int level, int branch) const
	{
		return dual_tree[level][branch]->subbands();
	}
	
	int level_W(int level, int branch) const
	{
		return dual_tree[level][branch]->width();
	}
    
	int level_H(int level, int branch) const
	{
		return dual_tree[level][branch]->height();
	}
	
	int level_pad(int level, int branch) const
	{
		return dual_tree[level][branch]->padding();
	}
	
	int maxlevel() const
	{
		return lvltot;
	}
	
    template<typename E>
    void reconstruct(E * dst);
	
};
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename E>
	cplx_wavelet_decomposition::cplx_wavelet_decomposition(E * src, int width, int height, int maxlvl, int subsampling)
	: lvltot(0), subsamp(subsampling), m_w(width), m_h(height)
	{
		
		//initialize wavelet filters
		
		first_lev_len = FSFarras_len;
		first_lev_offset = FSFarras_offset;
		first_lev_anal = new float[4*first_lev_len];
		first_lev_synth = new float[4*first_lev_len];

		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				for (int i=0; i<first_lev_len; i++) {
					first_lev_anal[first_lev_len*(2*n+m)+i]  = FSFarras_anal[n][m][i]/sqrt(2);
					first_lev_synth[first_lev_len*(2*n+m)+i] = FSFarras_anal[n][m][first_lev_len-1-i]/sqrt(2);
				}
			}
		}
		
		/*first_lev_len = AntonB_len;
		first_lev_offset = AntonB_offset;
		first_lev_anal = new float[4*first_lev_len];
		first_lev_synth = new float[4*first_lev_len];
		
		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				for (int i=0; i<first_lev_len; i++) {
					first_lev_anal[first_lev_len*(2*n+m)+i]  = AntonB_anal[n][m][i];
					first_lev_synth[first_lev_len*(2*n+m)+i] = 2*AntonB_synth[n][m][i];
				}
			}
		}*/
		
		wavfilt_len = Kingsbury_len;
		wavfilt_offset = Kingsbury_offset;
		wavfilt_anal = new float[4*wavfilt_len];
		wavfilt_synth = new float[4*wavfilt_len];
		
		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				for (int i=0; i<wavfilt_len; i++) {
					wavfilt_anal[wavfilt_len*(2*n+m)+i]  = Kingsbury_anal[n][m][i]/sqrt(2);
					wavfilt_synth[wavfilt_len*(2*n+m)+i] = Kingsbury_anal[n][m][first_lev_len-1-i]/sqrt(2);
				}
			}
		}
		
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

		
		// after coefficient rotation, data structure is:
		// dual_tree[scale][2*n+m=2*(Re/Im)+dir][channel={lo,hi1,hi2,hi3}][pixel_array]
		
		//srand((unsigned)time(0));
		//for (int i=0; i<m_w*m_h; i++ )
		//	src[i] = (float)rand()/(float)RAND_MAX;

		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				lvltot=0;
				float padding = 0;//1<<(maxlvl-1);
				dual_tree[0][2*n+m] = new wavelet_level<internal_type>(src, lvltot, subsamp, padding, m_w, m_h, first_lev_anal+first_lev_len*2*n, \
																				   first_lev_anal+first_lev_len*2*m, first_lev_len, first_lev_offset);
				while(lvltot < maxlvl) {
					lvltot++;
					dual_tree[lvltot][2*n+m] = new wavelet_level<internal_type>(dual_tree[lvltot-1][2*n+m]->lopass()/*lopass*/, lvltot, subsamp, 0/*no padding*/, \
																							dual_tree[lvltot-1][2*n+m]->width(), \
																							dual_tree[lvltot-1][2*n+m]->height(), \
																							wavfilt_anal+wavfilt_len*2*n, wavfilt_anal+wavfilt_len*2*m, \
																							wavfilt_len, wavfilt_offset);
				}
			}
		}
		
		
		//rotate detail coefficients
		float coeffave[5][4][3];

		float root2 = sqrt(2);
		for (int lvl=0; lvl<lvltot; lvl++) {
			int Wlvl = dual_tree[lvl][0]->width();
			int Hlvl = dual_tree[lvl][0]->height();
			for (int n=0; n<4; n++) 
				for (int m=1; m<4; m++)
					coeffave[lvl][n][m-1]=0;

			for (int m=1; m<4; m++) {//detail coefficients only
				for (int i=0; i<Wlvl*Hlvl; i++) {//pixel
					
					float wavtmp = (dual_tree[lvl][0]->wavcoeffs[m][i] + dual_tree[lvl][3]->wavcoeffs[m][i])/root2;
					dual_tree[lvl][3]->wavcoeffs[m][i] = (dual_tree[lvl][0]->wavcoeffs[m][i] - dual_tree[lvl][3]->wavcoeffs[m][i])/root2;
					dual_tree[lvl][0]->wavcoeffs[m][i] = wavtmp;
					
					wavtmp = (dual_tree[lvl][1]->wavcoeffs[m][i] + dual_tree[lvl][2]->wavcoeffs[m][i])/root2;
					dual_tree[lvl][2]->wavcoeffs[m][i] = (dual_tree[lvl][1]->wavcoeffs[m][i] - dual_tree[lvl][2]->wavcoeffs[m][i])/root2;
					dual_tree[lvl][1]->wavcoeffs[m][i] = wavtmp;
					
					for (int n=0; n<4; n++) coeffave[lvl][n][m-1] += fabs(dual_tree[lvl][n]->wavcoeffs[m][i]);
				}
			}
			for (int n=0; n<4; n++) 
				for (int i=0; i<3; i++) 
					coeffave[lvl][n][i] /= Wlvl*Hlvl;
		}

	}
	
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	
		
	template<typename E>
	void cplx_wavelet_decomposition::reconstruct(E * dst) { 
		
		// data structure is wavcoeffs[scale][2*n+m=2*(Re/Im)+dir][channel={lo,hi1,hi2,hi3}][pixel_array]
		
		//rotate detail coefficients
		float root2 = sqrt(2);
		for (int lvl=0; lvl<lvltot; lvl++) {
			int Wlvl = dual_tree[lvl][0]->width();
			int Hlvl = dual_tree[lvl][0]->height();
			for (int i=0; i<Wlvl*Hlvl; i++) {//pixel
				for (int m=1; m<4; m++) {//detail coefficients only
					float wavtmp = (dual_tree[lvl][0]->wavcoeffs[m][i] + dual_tree[lvl][3]->wavcoeffs[m][i])/root2;
					dual_tree[lvl][3]->wavcoeffs[m][i] = (dual_tree[lvl][0]->wavcoeffs[m][i] - dual_tree[lvl][3]->wavcoeffs[m][i])/root2;
					dual_tree[lvl][0]->wavcoeffs[m][i] = wavtmp;
					
					wavtmp = (dual_tree[lvl][1]->wavcoeffs[m][i] + dual_tree[lvl][2]->wavcoeffs[m][i])/root2;
					dual_tree[lvl][2]->wavcoeffs[m][i] = (dual_tree[lvl][1]->wavcoeffs[m][i] - dual_tree[lvl][2]->wavcoeffs[m][i])/root2;
					dual_tree[lvl][1]->wavcoeffs[m][i] = wavtmp;
				}
			}
		}
		
		internal_type ** tmp = new internal_type *[4];
		for (int i=0; i<4; i++) {
			tmp[i] = new internal_type[m_w*m_h];
		}
		
		for (int n=0; n<2; n++) {
			for (int m=0; m<2; m++) {
				int skip=1<<(lvltot-1);
				for (int lvl=lvltot-1; lvl>0; lvl--) {
					dual_tree[lvl][2*n+m]->reconstruct_level(dual_tree[lvl-1][2*n+m]->wavcoeffs[0], wavfilt_synth+wavfilt_len*2*n, \
																	wavfilt_synth+wavfilt_len*2*m, wavfilt_len, wavfilt_offset);
					skip /=2;
				}
				dual_tree[0][2*n+m]->reconstruct_level(tmp[2*n+m], first_lev_synth+first_lev_len*2*n, 
															  first_lev_synth+first_lev_len*2*m, first_lev_len, first_lev_offset);
			}
		}
				
		copy_out(tmp,dst,m_w*m_h);
		
		for (int i=0; i<4; i++) {
			delete[] tmp[i];
		}
		delete[] tmp;
		
		
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	
	class wavelet_decomposition
	{
	public:
		
		typedef float internal_type;
		
	private:
		
		static const int maxlevels = 9;//should be greater than any conceivable order of decimation
		
		int lvltot, subsamp;
		size_t m_w, m_h;//dimensions
		
		int wavfilt_len, wavfilt_offset;
		float *wavfilt_anal;
		float *wavfilt_synth;
		
		int testfilt_len, testfilt_offset;
		float *testfilt_anal;
		float *testfilt_synth;
		
		wavelet_level<internal_type> * wavelet_decomp[maxlevels];
		
	public:
		
		template<typename E>
		wavelet_decomposition(E * src, int width, int height, int maxlvl, int subsampling);
		
		~wavelet_decomposition();
		
		internal_type ** level_coeffs(int level) const
		{
			return wavelet_decomp[level]->subbands();
		}
		
		int level_W(int level) const
		{
			return wavelet_decomp[level]->width();
		}
		
		int level_H(int level) const
		{
			return wavelet_decomp[level]->height();
		}
		
		int level_pad(int level) const
		{
			return wavelet_decomp[level]->padding();
		}
		
		int level_stride(int level) const
		{
			return wavelet_decomp[level]->stride();
		}
		
		int maxlevel() const
		{
			return lvltot;
		}
		
		int subsample() const
		{
			return subsamp;
		}
		
		template<typename E>
		void reconstruct(E * dst);
		
	};
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename E>
	wavelet_decomposition::wavelet_decomposition(E * src, int width, int height, int maxlvl, int subsampling)
	: lvltot(0), subsamp(subsampling), m_w(width), m_h(height)
	{
		
		//initialize wavelet filters
		
		wavfilt_len = Daub4_len;
		wavfilt_offset = Daub4_offset;
		wavfilt_anal = new float[2*wavfilt_len];
		wavfilt_synth = new float[2*wavfilt_len];
		
		for (int n=0; n<2; n++) {
			for (int i=0; i<wavfilt_len; i++) {
				wavfilt_anal[wavfilt_len*(n)+i]  = Daub4_anal[n][i];
				wavfilt_synth[wavfilt_len*(n)+i] = Daub4_anal[n][wavfilt_len-1-i];
				//n=0 lopass, n=1 hipass
			}
		}
		
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
		
		
		// after coefficient rotation, data structure is:
		// wavelet_decomp[scale][channel={lo,hi1,hi2,hi3}][pixel_array]
		
		//srand((unsigned)time(0));
		//for (int i=0; i<m_w*m_h; i++ )
		//	src[i] = (float)rand()/(float)RAND_MAX;
		
		int padding = 0;//pow(2, maxlvl);//must be a multiple of two
		lvltot=0;
		wavelet_decomp[lvltot] = new wavelet_level<internal_type>(src, lvltot/*level*/, subsamp, padding/*padding*/, m_w, m_h, \
																  wavfilt_anal, wavfilt_anal, wavfilt_len, wavfilt_offset);
		while(lvltot < maxlvl) {
			lvltot++;
			wavelet_decomp[lvltot] = new wavelet_level<internal_type>(wavelet_decomp[lvltot-1]->lopass()/*lopass*/, lvltot/*level*/, subsamp, 0/*no padding*/, \
																	  wavelet_decomp[lvltot-1]->width(), wavelet_decomp[lvltot-1]->height(), \
																	  wavfilt_anal, wavfilt_anal, wavfilt_len, wavfilt_offset);
		}
		
		
	}
	
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	
	
	template<typename E>
	void wavelet_decomposition::reconstruct(E * dst) { 
		
		// data structure is wavcoeffs[scale][channel={lo,hi1,hi2,hi3}][pixel_array]
		
		//int skip=1<<(lvltot-1);
		for (int lvl=lvltot-1; lvl>0; lvl--) {
			wavelet_decomp[lvl]->reconstruct_level(wavelet_decomp[lvl-1]->wavcoeffs[0], wavfilt_synth, wavfilt_synth, wavfilt_len, wavfilt_offset);
			//skip /=2;
		}
		
		internal_type * tmp = new internal_type[m_w*m_h];

		wavelet_decomp[0]->reconstruct_level(tmp, wavfilt_synth, wavfilt_synth, wavfilt_len, wavfilt_offset);
		
		copy_out(tmp,dst,m_w*m_h);
		
		delete[] tmp;
		
		
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	

};

#endif
