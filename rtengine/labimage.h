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
#ifndef _LABIMAGE_H_
#define _LABIMAGE_H_

namespace rtengine
{

class LabImage
{
private:
    bool fromImage;
    void allocLab(int w, int h)
    {
        L = new float*[H];
        a = new float*[H];
        b = new float*[H];

        data = new float [W * H * 3];
        float * index = data;

        for (int i = 0; i < H; i++) {
            L[i] = index + i * W;
        }

        index += W * H;

        for (int i = 0; i < H; i++) {
            a[i] = index + i * W;
        }

        index += W * H;

        for (int i = 0; i < H; i++) {
            b[i] = index + i * W;
        }
    };
public:
    int W, H;
    float * data;
    float** L;
    float** a;
    float** b;

    LabImage (int w, int h);
    ~LabImage ();

    //Copies image data in Img into this instance.
    void CopyFrom(LabImage *Img);
    void getPipetteData (float &L, float &a, float &b, int posX, int posY, int squareSize);
    void deleteLab( )
    {
        if (!fromImage) {
            delete [] L;
            delete [] a;
            delete [] b;
            delete [] data;
        }
    }
    void reallocLab( )
    {
        allocLab(W, H);
    };

};

}
#endif
