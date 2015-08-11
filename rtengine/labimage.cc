#include "labimage.h"
#include <memory.h>
namespace rtengine
{

LabImage::LabImage (int w, int h) : fromImage(false), W(w), H(h)
{
    allocLab(w, h);
}

LabImage::~LabImage ()
{
    deleteLab();
}

void LabImage::CopyFrom(LabImage *Img)
{
    memcpy(data, Img->data, W * H * 3 * sizeof(float));
}

void LabImage::getPipetteData (float &v1, float &v2, float &v3, int posX, int posY, int squareSize)
{
    float accumulator_L = 0.f;
    float accumulator_a = 0.f;
    float accumulator_b = 0.f;
    unsigned long int n = 0;
    int halfSquare = squareSize / 2;

    for (int iy = posY - halfSquare; iy < posY - halfSquare + squareSize; ++iy) {
        for (int ix = posX - halfSquare; ix < posX - halfSquare + squareSize; ++ix) {
            if (ix >= 0 && iy >= 0 && ix < W && iy < H) {
                accumulator_L += L[iy][ix];
                accumulator_a += a[iy][ix];
                accumulator_b += b[iy][ix];
                ++n;
            }
        }
    }

    v1 = n ? accumulator_L / float(n) : 0.f;
    v2 = n ? accumulator_a / float(n) : 0.f;
    v3 = n ? accumulator_b / float(n) : 0.f;
}

}
