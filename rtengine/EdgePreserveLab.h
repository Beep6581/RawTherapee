#pragma once
/*
The EdgePreserveLab files contain standard C++ (standard except the first line) code for creating and, to a
limited extent (create your own uses!), messing with multi scale edge preserving decompositions of a 32 bit single channel
image. As a byproduct it contains a lot of linear algebra which can be useful for optimization problems that
you want to solve in rectangles on rectangular grids.

Anyway. Basically, this is an implementation of what's presented in the following papers:
	Edge-Preserving Decompositions for Multi-Scale Tone and Detail Manipulation
	An Iterative Solution Method for Linear Systems of Which the Coefficient Matrix is a Symetric M-Matrix
	Color correction for tone mapping
	Wikipedia, the free encyclopedia

First one is most of what matters, next two are details, last everything else. I did a few things differently, especially:
	Reformulated the minimization with finite elements instead of finite differences. This results in better conditioning,
	slightly better accuracy (less artifacts), the possibility of a better picked edge stopping function, but more memory consumption.

	A single rotationally invariant edge stopping function is used instead of two non-invariant ones.

	Incomplete Cholseky factorization instead of Szeliski's LAHBF. Slower, but not subject to any patents.

	For tone mapping, original images are decomposed instead of their logarithms, and just one decomposition is made;
	I find that this way works plenty good (theirs isn't better or worse... just different) and is simpler.

Written by ben_pcc in Portland, Oregon, USA; code modified by Emil Martinec.

 // EdgePreserveLab.h and EdgePreserveLab.cpp are free software: 
 // you can redistribute it and/or modify it
 //	under the terms of the GNU General Public License as published by
 //	the Free Software Foundation, either version 3 of the License, or
 //	(at your option) any later version.
 //
 //	This program is distributed in the hope that it will be useful,
 //	but WITHOUT ANY WARRANTY; without even the implied warranty of
 //	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 //	GNU General Public License for more details.
 //
 //	You should have received a copy of the GNU General Public License
 //	along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/



#include <cmath>
#include <stdio.h>
#include <string.h>
#include "EdgePreservingDecomposition.h"

class EdgePreserveLab{
public:
	EdgePreserveLab(unsigned int width, unsigned int height);
	~EdgePreserveLab();

	//Create an edge preserving blur of Source. Will create and return, or fill into Blur if not NULL. In place not ok.
	//If UseBlurForEdgeStop is true, supplied not NULL Blur is used to calculate the edge stopping function instead of Source.
	float *CreateBlur(float *Source, float LScale, float abScale, float EdgeStoppingLuma, float EdgeStoppingChroma, unsigned int Iterates, float *Blur = NULL, bool UseBlurForEdgeStop = false);

	//Iterates CreateBlur such that the smoothness term approaches a specific norm via iteratively reweighted least squares. In place not ok.
	float *CreateIteratedBlur(float *Source, float LScale, float abScale, float EdgeStoppingLuma, float EdgeStoppingChroma, unsigned int Iterates, unsigned int Reweightings, float *Blur = NULL);

	/*Lowers global contrast while preserving or boosting local contrast. Can fill into Compressed. The smaller Compression
	the more compression is applied, with Compression = 1 giving no effect and above 1 the opposite effect. You can totally
	use Compression = 1 and play with DetailBoost for some really sweet unsharp masking. If working on luma/grey, consider giving it a logarithm.	
	In place calculation to save memory (Source == Compressed) is totally ok. Reweightings > 0 invokes CreateIteratedBlur instead of CreateBlur. */
	float *CompressDynamicRange(float *Source, float LScale = 1.0f, float abScale = 5.0f, float EdgeStoppingLuma = 1.0f, float EdgeStoppingChroma = 1.0f, 
								float CompressionExponent = 0.8f, float DetailBoost = 0.1f, unsigned int Iterates = 20, 
								unsigned int Reweightings = 0, float *Compressed = NULL);

private:
	MultiDiagonalSymmetricMatrix *A;	//The equations are simple enough to not mandate a matrix class, but fast solution NEEDS a complicated preconditioner.
	unsigned int w, h, n;

	//Convenient access to the data in A.
	float *a0, *a_1, *a_w, *a_w_1, *a_w1;
};
