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

#include "improcfun.h"
#include "alpha.h"
#include "procparams.h"
#include "imagesource.h"
#include <iostream>

namespace
{

// "ceil" rounding
template<typename T>
constexpr T skips(T a, T b)
{
    return a / b + static_cast<bool>(a % b);
}

}

namespace rtengine
{

/* Code taken from Gimp 2.8.10 and converted for RawTherapee by Jean-Christophe FRISCH (aka Hombre) on 02.19.2014
 *
 * ORIGINAL NOTES
 *
 *     The method used here is similar to the lighting invariant correction
 *     method but slightly different: we do not divide the RGB components,
 *     but substract them I2 = I0 - I1, where I0 is the sample image to be
 *     corrected, I1 is the reference pattern. Then we solve DeltaI=0
 *     (Laplace) with I2 Dirichlet conditions at the borders of the
 *     mask. The solver is a unoptimized red/black checker Gauss-Siedel
 *     with an over-relaxation factor of 1.8. It can benefit from a
 *     multi-grid evaluation of an initial solution before the main
 *     iteration loop.
 *
 *     I reduced the convergence criteria to 0.1% (0.001) as we are
 *     dealing here with RGB integer components, more is overkill.
 *
 *     Jean-Yves Couleaud cjyves@free.fr
 */

/* Original Algorithm Design:
 *
 * T. Georgiev, "Photoshop Healing Brush: a Tool for Seamless Cloning
 * http://www.tgeorgiev.net/Photoshop_Healing.pdf
 */

#if 1

class SpotBox {

public:
    enum class Type {
        SOURCE,
        TARGET,
        FINAL
    };

    struct Rectangle {
    public:
        int x1;
        int y1;
        int x2;
        int y2;

        Rectangle() : x1(0), y1(0), x2(0), y2(0) {}
        Rectangle(const Rectangle &other) : x1(other.x1), y1(other.y1), x2(other.x2), y2(other.y2) {}
        Rectangle(int X1, int Y1, int X2, int Y2) : x1(X1), y1(Y1), x2(X2), y2(Y2) {}

        bool intersects(const Rectangle &other) const {
            return (other.x1 <= x2 && other.x2 >= x1)
                && (other.y1 <= y2 && other.y2 >= y1);
        }

        bool getIntersection(const Rectangle &other, std::unique_ptr<Rectangle> &intersection) const {
            if (intersects(other)) {
                if (!intersection) {
                    intersection.reset(new Rectangle());
                }
                intersection->x1 = rtengine::max(x1, other.x1);
                intersection->x2 = rtengine::min(x2, other.x2);
                intersection->y1 = rtengine::max(y1, other.y1);
                intersection->y2 = rtengine::min(y2, other.y2);

                if (intersection->x1 > intersection->x2 || intersection->y1 > intersection->y2) {
                    intersection.release();
                    return false;
                }
                return true;
            }
            if (intersection) {
                intersection.release();
            }
            return false;
        }

        Rectangle& operator+=(const Coord &v) {
            x1 += v.x;
            y1 += v.y;
            x2 += v.x;
            y2 += v.y;
            return *this;
        }

        Rectangle& operator-=(const Coord &v) {
            x1 -= v.x;
            y1 -= v.y;
            x2 -= v.x;
            y2 -= v.y;
            return *this;
        }

        Rectangle& operator/=(const int &v) {
            if (v == 1) {
                return *this;
            }

            /*
            float fv = float(v);
            x1 = int(float(x1) / fv + 0.5f);
            y1 = int(float(y1) / fv + 0.5f);
            x2 = int(float(x2) / fv + 0.5f);
            y2 = int(float(y2) / fv + 0.5f);
            */

            // Aletrnate rounding possibility
            int w = x2 - x1 + 1;
            int h = x2 - x1 + 1;
            w = w / v + (w % v > 0);
            h = h / v + (h % v > 0);
            x1 /= v;
            y1 /= v;
            x2 = x1 + w - 1;
            y2 = y1 + h - 1;

            return *this;
        }
};

private:
    Type type;
    Imagefloat* image;

public:
    // top/left and bottom/right coordinates of the spot in image space (at some point divided by scale factor)
    Rectangle spotArea;
    // top/left and bottom/right coordinates of the spot in scaled image space (on borders, imgArea won't cover spotArea)
    Rectangle imgArea;
    // top/left and bottom/right coordinates of useful part of the image in scaled image space (rounding error workaround)
    Rectangle intersectionArea;
    float radius;
    float featherRadius;

    SpotBox (int tl_x, int tl_y, int br_x, int br_y, int radius, int feather_radius, Imagefloat* image, Type type) :
       type(type),
       image(image),
       spotArea(tl_x, tl_y, br_x, br_y),
       imgArea(spotArea),
       intersectionArea(),
       radius(radius),
       featherRadius(feather_radius)
    {}

    SpotBox (int tl_x, int tl_y, int radius, int feather_radius, Imagefloat* image, Type type) :
       type(type),
       image(image),
       spotArea(tl_x, tl_y, image ? tl_x + image->getWidth() - 1 : 0, image ? tl_y + image->getHeight() - 1 : 0),
       imgArea(spotArea),
       intersectionArea(),
       radius(radius),
       featherRadius(feather_radius)
    {}

    SpotBox (SpotEntry &spot, Type type) :
        type(type),
        image(nullptr),
        intersectionArea(),
        radius(spot.radius),
        featherRadius(int(spot.getFeatherRadius() + 0.5f))  // rounding to int before resizing
    {
        spotArea.x1 = int ((type == Type::SOURCE ? spot.sourcePos.x : spot.targetPos.x) - featherRadius);
        spotArea.x2 = int ((type == Type::SOURCE ? spot.sourcePos.x : spot.targetPos.x) + featherRadius);
        spotArea.y1 = int ((type == Type::SOURCE ? spot.sourcePos.y : spot.targetPos.y) - featherRadius);
        spotArea.y2 = int ((type == Type::SOURCE ? spot.sourcePos.y : spot.targetPos.y) + featherRadius);
        imgArea = spotArea;
    }

    ~SpotBox() {
        if (image && type != Type::FINAL) {
            delete image;
        }
    }

    SpotBox& operator /=(const int& v) {
        if (v == 1) {
            return *this;
        }
        spotArea /= v;
        imgArea /= v;
        radius = radius / float(v);
        featherRadius = getWidth() / 2.f;
        // intersectionArea doesn't need resize, because it's set after resizing
        return *this;
    }

    int getWidth() {
        return spotArea.x2 - spotArea.x1 + 1;
    }

    int getHeight() {
        return spotArea.y2 - spotArea.y1 + 1;
    }

    int getImageWidth() {
        return imgArea.x2 - imgArea.x1 + 1;
    }

    int getImageHeight() {
        return imgArea.y2 - imgArea.y1 + 1;
    }

    int getIntersectionWidth() {
        return intersectionArea.x2 - intersectionArea.x1 + 1;
    }

    int getIntersectionHeight() {
        return intersectionArea.y2 - intersectionArea.y1 + 1;
    }

    bool checkImageSize() {
        if (!image || getImageWidth() != image->getWidth() || getImageHeight() != image->getHeight()) {
            return false;
        }
        return true;
    }

    void tuneImageSize() {
        if (!image) {
            return;
        }
        if (getImageWidth() > image->getWidth()) {
            imgArea.x2 = imgArea.x1 + image->getWidth() - 1;
        }
        if (getImageHeight() > image->getHeight()) {
            imgArea.y2 = imgArea.y1 + image->getHeight() - 1;
        }
    }

    Imagefloat *getImage() {  // TODO: this should send back a const value, but getImage don't want it to be const...
        return image;
    }

    void allocImage() {
        int newW = imgArea.x2 - imgArea.x1 + 1;
        int newH = imgArea.y2 - imgArea.y1 + 1;

        if (image && type != Type::FINAL && (image->getWidth() != newW || image->getHeight() != newH)) {
            delete image;
            image = nullptr;
        }
        if (image == nullptr) {
            image = new Imagefloat(newW, newH);
        }
    }

    bool spotIntersects(const SpotBox &other) const {
        return spotArea.intersects(other.spotArea);
    }

    bool getSpotIntersection(const SpotBox &other, std::unique_ptr<Rectangle> &intersection) const {
        return spotArea.getIntersection(other.spotArea, intersection);
    }

    bool imageIntersects(const SpotBox &other, bool atDestLocation=false) const {
        if (atDestLocation) {
            Coord v(other.spotArea.x1 - spotArea.x1, other.spotArea.y1 - spotArea.y1);
            Rectangle imgArea2(imgArea.x1, imgArea.y1, imgArea.x2, imgArea.y2);
            imgArea2 += v;
            return imgArea2.intersects(other.imgArea);
        }
        return imgArea.intersects(other.imgArea);
    }

    bool mutuallyClipImageArea(SpotBox &other) {
        Coord v(other.spotArea.x1 - spotArea.x1, other.spotArea.y1 - spotArea.y1);
        Rectangle imgArea2 = imgArea;
        imgArea2 += v;
        std::unique_ptr<Rectangle> intersection;
        if (!imgArea2.getIntersection(other.imgArea, intersection)) {
            return false;
        }
        other.intersectionArea = *intersection;
        Coord v2(-v.x, -v.y);
        *intersection -= v;
        intersectionArea = *intersection;
        return true;
    }

    bool setIntersectionWith(const SpotBox &other) {
        if (!spotIntersects(other)) {
            return false;
        }
        imgArea.x1 = rtengine::max(spotArea.x1, other.spotArea.x1);
        imgArea.x2 = rtengine::min(spotArea.x2, other.spotArea.x2);
        imgArea.y1 = rtengine::max(spotArea.y1, other.spotArea.y1);
        imgArea.y2 = rtengine::min(spotArea.y2, other.spotArea.y2);
        if (imgArea.x1 > imgArea.x2 || imgArea.y1 > imgArea.y2) {
            return false;
        }
        return true;
    }

    bool processIntersectionWith(SpotBox &destBox) {
        Imagefloat *dstImg = destBox.image;

        if (image == nullptr || dstImg == nullptr) {
            std::cerr << "One of the source or destination SpotBox image is missing !" << std::endl;
            return false;
        }

        /*
        printf("#%llu: spotBox(X1:%d, Y1:%d, X2:%d, Y2:%d, W:%d, H:%d)  img(W:%d, H:%d)\n", i,
                spotArea.x1, spotArea.y1, spotArea.x2, spotArea.y2,
                getWidth(), getHeight(),
                destBox.getWidth(), destBox.getHeight()
                );
        */

#ifndef NDEBUG
        printf("[processIntersectionWith] srcSpotBox @ %p(%d, %d) - img @ %p(%d, %d)  /  dstSpotBox @ %p(%d, %d) - img @ %p(%d, %d)\n",
                this, getWidth(), getHeight(), image, destBox.getWidth(), destBox.getHeight(),
                &destBox, destBox.getWidth(), destBox.getHeight(), destBox.image, dstImg->getWidth(), dstImg->getHeight());
        printf("                   [spot] srcSpotBox (%d, %d) -> (%d, %d)  /  dstSpotBox (%d, %d) -> (%d, %d)\n",
                spotArea.x1, spotArea.y1, spotArea.x2, spotArea.y2,
                destBox.spotArea.x1, destBox.spotArea.y1, destBox.spotArea.x2, destBox.spotArea.y2);
        printf("                   [img ] srcSpotBox (%d, %d) -> (%d, %d)  /  dstSpotBox (%d, %d) -> (%d, %d)\n",
                imgArea.x1, imgArea.y1, imgArea.x2, imgArea.y2,
                destBox.imgArea.x1, destBox.imgArea.y1, destBox.imgArea.x2, destBox.imgArea.y2);
        printf("                   [intr] srcSpotBox (%d, %d) -> (%d, %d)  /  dstSpotBox (%d, %d) -> (%d, %d)\n",
                intersectionArea.x1, intersectionArea.y1, intersectionArea.x2, intersectionArea.y2,
                destBox.intersectionArea.x1, destBox.intersectionArea.y1, destBox.intersectionArea.x2, destBox.intersectionArea.y2);
        printf("                   radius = %.3f\n", radius);
#endif

        int srcImgY = intersectionArea.y1 - imgArea.y1;
        int dstImgY = destBox.intersectionArea.y1 - destBox.imgArea.y1;
        for (int y = intersectionArea.y1; y <= intersectionArea.y2; ++y) {
            float  dy = float(y - spotArea.y1) - featherRadius;

            int srcImgX = intersectionArea.x1 - imgArea.x1;
            int dstImgX = destBox.intersectionArea.x1 - destBox.imgArea.x1;
            for (int x = intersectionArea.x1; x <= intersectionArea.x2; ++x) {
                float dx = float(x - spotArea.x1) - featherRadius;
                float r = sqrt(dx * dx + dy * dy);

                if (r >= featherRadius) {
                    ++srcImgX;
                    ++dstImgX;
                    continue;
                }
                if (r <= radius) {
                    dstImg->r(dstImgY, dstImgX) = image->r(srcImgY, srcImgX);
                    dstImg->g(dstImgY, dstImgX) = image->g(srcImgY, srcImgX);
                    dstImg->b(dstImgY, dstImgX) = image->b(srcImgY, srcImgX);
                } else {
                    float opacity = (featherRadius - r) / (featherRadius - radius);
                    dstImg->r(dstImgY, dstImgX) = (image->r(srcImgY, srcImgX) - dstImg->r(dstImgY, dstImgX)) * opacity + dstImg->r(dstImgY,dstImgX);
                    dstImg->g(dstImgY, dstImgX) = (image->g(srcImgY, srcImgX) - dstImg->g(dstImgY, dstImgX)) * opacity + dstImg->g(dstImgY,dstImgX);
                    dstImg->b(dstImgY, dstImgX) = (image->b(srcImgY, srcImgX) - dstImg->b(dstImgY, dstImgX)) * opacity + dstImg->b(dstImgY,dstImgX);
                }
                ++srcImgX;
                ++dstImgX;
            }
            ++srcImgY;
            ++dstImgY;
        }

        return true;
    }

    // Copy the intersecting part
    bool copyImgTo(SpotBox &destBox) {
        Imagefloat *destImg = destBox.image;

        if (image == nullptr || destImg == nullptr) {
            std::cerr << "One of the source or destination SpotBox image is missing !" << std::endl;
            return false;
        }

        /*
        printf("#%llu: spotBox(X1:%d, Y1:%d, X2:%d, Y2:%d, W:%d, H:%d)  img(W:%d, H:%d)\n", i,
                spotArea.x1, spotArea.y1, spotArea.x2, spotArea.y2,
                getWidth(), getHeight(),
                destBox.getWidth(), destBox.getHeight()
                );
        */

        std::unique_ptr<Rectangle> intersection;

        if (!intersectionArea.getIntersection(destBox.intersectionArea, intersection)) {
            return false;
        }

        Imagefloat *srcImg = image;
        Imagefloat *dstImg = destBox.image;

#ifndef NDEBUG
        printf("[copyImgTo] srcSpotBox @ %p(%d, %d) - img @ %p(%d, %d)  /  dstSpotBox @ %p(%d, %d) - img @ %p(%d, %d)\n",
                this, getWidth(), getHeight(), image, srcImg->getWidth(), srcImg->getHeight(),
                &destBox, destBox.getWidth(), destBox.getHeight(), destBox.image, dstImg->getWidth(), dstImg->getHeight());
        printf("            [spot] srcSpotBox (%d, %d) -> (%d, %d)  /  dstSpotBox (%d, %d) -> (%d, %d)\n",
                spotArea.x1, spotArea.y1, spotArea.x2, spotArea.y2,
                destBox.spotArea.x1, destBox.spotArea.y1, destBox.spotArea.x2, destBox.spotArea.y2);
        printf("            [img ] srcSpotBox (%d, %d) -> (%d, %d)  /  dstSpotBox (%d, %d) -> (%d, %d)\n",
                imgArea.x1, imgArea.y1, imgArea.x2, imgArea.y2,
                destBox.imgArea.x1, destBox.imgArea.y1, destBox.imgArea.x2, destBox.imgArea.y2);
        printf("            [intr] intersection (%d, %d) -> (%d, %d)\n",
                intersection->x1, intersection->y1, intersection->x2, intersection->y2);
#endif

        int srcImgY = intersection->y1 - imgArea.y1;
        int dstImgY = intersection->y1 - destBox.imgArea.y1;
        for (int y = intersection->y1; y <= intersection->y2; ++y) {
            int srcImgX = intersection->x1 - imgArea.x1;
            int dstImgX = intersection->x1 - destBox.imgArea.x1;

            for (int x = intersection->x1; x <= intersection->x2; ++x) {
                dstImg->r(dstImgY, dstImgX) = srcImg->r(srcImgY, srcImgX);
                dstImg->g(dstImgY, dstImgX) = srcImg->g(srcImgY, srcImgX);
                dstImg->b(dstImgY, dstImgX) = srcImg->b(srcImgY, srcImgX);
                ++srcImgX;
                ++dstImgX;
            }
            ++srcImgY;
            ++dstImgY;
        }

        return true;
    }
};

void ImProcFunctions::removeSpots (Imagefloat* img, ImageSource* imgsrc, const std::vector<SpotEntry> &entries, const PreviewProps &pp, const ColorTemp &currWB, int tr)
{
    //Get the clipped image areas (src & dst) from the source image

    std::vector< std::shared_ptr<SpotBox> > srcSpotBoxs;
    std::vector< std::shared_ptr<SpotBox> > dstSpotBoxs;
    int fullImgWidth = 0;
    int fullImgHeight = 0;
    imgsrc->getFullSize(fullImgWidth, fullImgHeight, tr);
    SpotBox fullImageBox(0, 0, fullImgWidth - 1, fullImgHeight - 1, 0, 0, nullptr, SpotBox::Type::FINAL);
    SpotBox cropBox(pp.getX(), pp.getY(),
                    pp.getX() + pp.getWidth() - 1, pp.getY() + pp.getHeight() - 1,
                    0, 0, img, SpotBox::Type::FINAL);

    std::set<int> visibleSpots;   // list of dest spots intersecting the preview's crop
    int i = 0;

    for (auto entry : params->spot.entries) {
        std::shared_ptr<SpotBox> srcSpotBox(new SpotBox(entry,  SpotBox::Type::SOURCE));
        std::shared_ptr<SpotBox> dstSpotBox(new SpotBox(entry,  SpotBox::Type::TARGET));
        if (   !srcSpotBox->setIntersectionWith(fullImageBox)
            || !dstSpotBox->setIntersectionWith(fullImageBox)
            || !srcSpotBox->imageIntersects(*dstSpotBox, true))
        {
            continue;
            ++i;
        }

        // If spot intersect the preview image, add it to the visible spots
        if (dstSpotBox->spotIntersects(cropBox)) {
#ifndef NDEBUG
            printf("Spot visible insere: %d\n", i);
#endif
            visibleSpots.insert(i);
        }
        ++i;

        // Source area
        PreviewProps spp(srcSpotBox->imgArea.x1, srcSpotBox->imgArea.y1,
                         srcSpotBox->getImageWidth(), srcSpotBox->getImageHeight(), pp.getSkip());
        int w = 0;
        int h = 0;
        imgsrc->getSize(spp, w, h);
        *srcSpotBox /= pp.getSkip();
        srcSpotBox->allocImage();
        Imagefloat *_image_ = srcSpotBox->getImage();
        for (int y = 0; y < (int)_image_->getHeight(); ++y) {
            for (int x = 0; x < (int)_image_->getWidth(); ++x) {
                _image_->r(y, x) = 60000.f;
                _image_->g(y, x) = 500.f;
                _image_->b(y, x) = 500.f;
            }
        }

#ifndef NDEBUG
        printf("~~~~~ Image Size Before : %d, %d\n", _image_->getWidth(), _image_->getHeight());
#endif
        imgsrc->getImage(currWB, tr, srcSpotBox->getImage(), spp, params->toneCurve, params->raw);
#ifndef NDEBUG
        printf("~~~~~ Image Size After  : %d, %d\n", _image_->getWidth(), _image_->getHeight());
#endif
        assert(srcSpotBox->checkImageSize());
        //printf("  /  src size: %d,%d", currSrcSpot->getWidth(), currSrcSpot->getHeight());


        // Destination area
        spp.set(dstSpotBox->imgArea.x1, dstSpotBox->imgArea.y1, dstSpotBox->getImageWidth(),
                dstSpotBox->getImageHeight(), pp.getSkip());
        *dstSpotBox /= pp.getSkip();
        dstSpotBox->allocImage();
        _image_ = dstSpotBox->getImage();
        for (int y = 0; y < (int)_image_->getHeight(); ++y) {
            for (int x = 0; x < (int)_image_->getWidth(); ++x) {
                _image_->r(y, x) = 500.f;
                _image_->g(y, x) = 500.f;
                _image_->b(y, x) = 60000.f;
            }
        }
#ifndef NDEBUG
        printf("~~~~~ Image Size Before : %d, %d\n", _image_->getWidth(), _image_->getHeight());
#endif
        imgsrc->getImage(currWB, tr, dstSpotBox->getImage(), spp, params->toneCurve, params->raw);
#ifndef NDEBUG
        printf("~~~~~ Image Size After  : %d, %d\n", _image_->getWidth(), _image_->getHeight());
#endif
        assert(dstSpotBox->checkImageSize());
        //printf("  /  src size: %d,%d", currDstSpot->getWidth(), currDstSpot->getHeight());

        // Update the intersectionArea between src and dest
        if (srcSpotBox->mutuallyClipImageArea(*dstSpotBox)) {
            srcSpotBoxs.push_back(srcSpotBox);
            dstSpotBoxs.push_back(dstSpotBox);

#ifndef NDEBUG
            printf("Taille de l'image pour spot #%d:  src(%d, %d)  /  dst(%d, %d)\n", i-1,
                    srcSpotBox->getImageWidth(), srcSpotBox->getImageHeight(),
                    dstSpotBox->getImageWidth(), dstSpotBox->getImageHeight()
                    );
#endif
        }

    }

    // Construct list of upstream dependancies

    std::set<int> requiredSpots = visibleSpots;  // starting point, visible spots are necessarilly required spots
#ifndef NDEBUG
    printf("Required spots: ");
#endif
    for (auto i = requiredSpots.rbegin(); i != requiredSpots.rend(); i++) {
        int spotNbr = *i;
#ifndef NDEBUG
        printf("%d ", spotNbr);
#endif
        requiredSpots.insert(spotNbr);
        if (spotNbr > 0) {
            for (int j = spotNbr - 1; j >= 0; --j) {
                if ((srcSpotBoxs.at(spotNbr))->imageIntersects(*dstSpotBoxs.at(j))) {
                    requiredSpots.insert(spotNbr);
#ifndef NDEBUG
                    printf("%d ", spotNbr);
#endif
                }
            }
        }
    }
#ifndef NDEBUG
    printf("\n");
#endif

    // Process spots and copy them downstream

    for (auto i = requiredSpots.begin(); i != requiredSpots.end(); i++) {
        // Process
#ifndef NDEBUG
        printf("========>>> Processing spot #%d\n", *i);
#endif
        srcSpotBoxs.at(*i)->processIntersectionWith(*dstSpotBoxs.at(*i));

        // Propagate
        std::set<int> positiveSpots;  // For DEBUG purpose only !
        auto j = i;
        ++j;
        while (j != requiredSpots.end()) {
            bool intersectionFound = false;
            int i_ = *i;
            int j_ = *j;
#ifndef NDEBUG
            printf("Propagating to SRC #%d ?\n", j_);
#endif
            intersectionFound |= dstSpotBoxs.at(i_)->copyImgTo(*srcSpotBoxs.at(j_));
#ifndef NDEBUG
            printf("Propagating to DSTC #%d ?\n", j_);
#endif
            intersectionFound |= dstSpotBoxs.at(i_)->copyImgTo(*dstSpotBoxs.at(j_));
            if (intersectionFound) {
                positiveSpots.insert(j_);
            }
            ++j;
        }
#ifndef NDEBUG
        printf("||| DEST spot #%d propagated to : ", *i);
        for (auto g : positiveSpots) {
            printf("%d ", g);
        }
        printf("\n");
#endif
    }

    // Copy the dest spot to the preview image
#ifndef NDEBUG
    printf("cropBox (%d, %d) -> (%d, %d)",
            cropBox.imgArea.x1, cropBox.imgArea.y1,
            cropBox.imgArea.x2, cropBox.imgArea.y2);
#endif
    cropBox /= pp.getSkip();
#ifndef NDEBUG
    printf("  ->  [/skip] (%d, %d) -> (%d, %d)",
            cropBox.imgArea.x1, cropBox.imgArea.y1,
            cropBox.imgArea.x2, cropBox.imgArea.y2);
#endif
    cropBox.tuneImageSize();
    cropBox.intersectionArea = cropBox.imgArea;
#ifndef NDEBUG
    printf("  ->  [/tuned] (%d, %d) -> (%d, %d)\n",
            cropBox.imgArea.x1, cropBox.imgArea.y1,
            cropBox.imgArea.x2, cropBox.imgArea.y2);
#endif

    int f = 0;
    for (auto i : visibleSpots) {
        f += dstSpotBoxs.at(i)->copyImgTo(cropBox) ? 1 : 0;
    }
#ifndef NDEBUG
    printf("Nombre de copie finale : %d\n", f);
#endif

    /*
    printf("#--: spotBox(X1:%d, Y1:%d, X2:%d, Y2:%d, W:%d, H:%d)  img(W:%d, H:%d)\n\n",
            cropBox.spotArea.x1, cropBox.spotArea.y1, cropBox.spotArea.x2, cropBox.spotArea.y2,
            cropBox.getWidth(), cropBox.getHeight(),
            cropBox.image->getWidth(), cropBox.image->getHeight()
            );
    */
}

#endif





#if 0
void ImProcFunctions::removeSpots (Imagefloat* img, ImageSource* imgsrc, const std::vector<SpotEntry> &entries, const PreviewProps &pp, const ColorTemp &currWB, int tr)
{

    Alpha mask;


    for (const auto entry : entries) {
        float featherRadius = entry.getFeatherRadius();
        int scaledFeatherRadius = featherRadius / pp.getSkip ();

        SpotBox srcBox(entry, SpotBox::Type::SOURCE);
        srcBox.translate(-pp.getX(), -pp.getY());
        srcBox /= float (pp.getSkip());

        SpotBox dstBox(entry, SpotBox::Type::TARGET);
        dstBox.translate(-pp.getX(), -pp.getY());
        dstBox /= float (pp.getSkip());

        //printf("  -> X: %04d > %04d\n  -> Y: %04d > %04d\n", dst_XMin, dst_XMax, dst_YMin, dst_YMax);

        // scaled spot is too small, we do not preview it
        if (scaledFeatherRadius < 2 && pp.getSkip() != 1) {
#ifndef NDEBUG

            if (options.rtSettings.verbose) {
                printf ("Skipping spot located at %d x %d, too small for the preview zoom rate\n", entry.sourcePos.x, entry.sourcePos.y);
            }

#endif
            continue;
        }

        // skipping entries totally transparent
        if (entry.opacity == 0.) {
#ifndef NDEBUG

            if (options.rtSettings.verbose) {
                printf ("Skipping spot located at %d x %d: opacity=%.3f\n", entry.sourcePos.x, entry.sourcePos.y, entry.opacity);
            }

            continue;
#endif
        }

        // skipping entries where the source circle isn't inside the image bounds, even partially
        if (src_XMin < 0 || src_XMax >= img->getWidth() || src_YMin < 0 || src_YMax >= img->getHeight()) {
#ifndef NDEBUG

            if (options.rtSettings.verbose) {
                printf ("Skipping spot located at %d x %d, from the data at %d x %d, radius=%d, feather=%.3f, opacity=%.3f: source out of bounds\n", entry.sourcePos.x, entry.sourcePos.y, entry.targetPos.x, entry.targetPos.y, entry.radius, entry.feather, entry.opacity);
                printf ("%d < 0 || %d >= %d || %d < 0 || %d >= %d\n",
                        src_XMin, src_XMax, img->getWidth(), src_YMin, src_YMax, img->getHeight());
            }

#endif
            continue;
        }

        // skipping entries where the dest circle is completely outside the image bounds
        /*
        if (dst_XMin >= img->getWidth() || dst_XMax <= 0 || dst_YMin >= img->getHeight() || dst_YMax <= 0) {
#ifndef NDEBUG

            if (options.rtSettings.verbose) {
                printf ("Skipping spot located at %d x %d, from the data at %d x %d, radius=%d, feather=%.3f, opacity=%.3f: source out of bounds\n", entry.sourcePos.x, entry.sourcePos.y, entry.targetPos.x, entry.targetPos.y, entry.radius, entry.feather, entry.opacity);
                printf ("%d >= %d || %d <= 0 || %d >= %d || %d <= 0\n",
                        dst_XMin, img->getWidth(), dst_XMax, dst_YMin, img->getHeight(), dst_YMax);
            }

#endif
            continue;
        }
        */

        // ----------------- Core function -----------------

#if 0
        int scaledPPX = pp.getX() / pp.skip;
        int scaledPPY = pp.getY() / pp.skip;
        int scaledPPW = pp.getWidth() / pp.skip + (pp.getWidth() % pp.getSkip() > 0);
        int scaledPPH = pp.getHeight() / pp.skip + (pp.getHeight() % pp.skip > 0);

        int sizeX = dst_XMax - dst_XMin + 1;
        int sizeY = dst_YMax - dst_YMin + 1;

        Imagefloat matrix (sizeX, sizeY);
        Imagefloat solution (sizeX, sizeY);

        // allocate the mask and draw it
        mask.setSize (sizeX, sizeY);
        {
            Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create (mask.getSurface());

            // clear the bitmap
            cr->set_source_rgba (0., 0., 0., 0.);
            cr->rectangle (0., 0., sizeX, sizeY);
            cr->set_line_width (0.);
            cr->fill();

            // draw the mask
            cr->set_antialias (Cairo::ANTIALIAS_GRAY);
            cr->set_line_width (featherRadius);
            double gradientCenterX = double (sizeX) / 2.;
            double gradientCenterY = double (sizeY) / 2.;
            {
                Cairo::RefPtr<Cairo::RadialGradient> radialGradient = Cairo::RadialGradient::create (
                            gradientCenterX, gradientCenterY, radius,
                            gradientCenterX, gradientCenterY, featherRadius
                        );
                radialGradient->add_color_stop_rgb (0., 0., 0., 1.);
                radialGradient->add_color_stop_rgb (1., 0., 0., 0.);
                cr->set_source_rgba (0., 0., 0., 1.);
                cr->mask (radialGradient);
                cr->rectangle (0., 0., sizeX, sizeY);
                cr->fill();
            }
        }

        // copy the src part to a temporary buffer to avoid possible self modified source
        Imagefloat *srcBuff = img->copySubRegion (srcX, srcY, sizeX, sizeY);


        // subtract pattern to image and store the result as a double in matrix
        for (int i = 0, i2 = dst_YMin;     i2 < sizeY - 1;     ++i, ++i2) {
            for (int j = 0, j2 = dst_XMin;     i2 < sizeX - 1;     ++j, ++j2) {
                matrix.r (i, j) = img->r (i2, j2) - srcBuff->r (i, j);
                matrix.g (i, j) = img->g (i2, j2) - srcBuff->g (i, j);
                matrix.b (i, j) = img->b (i2, j2) - srcBuff->b (i, j);
            }
        }


        // FIXME: is a faster implementation needed?
#define EPSILON   0.001
#define MAX_ITER  500

        // repeat until convergence or max iterations
        for (int n = 0; n < MAX_ITER; ++n) {

            printf ("<<< n=#%d\n", n);
            // ----------------------------------------------------------------

            /* Perform one iteration of the Laplace solver for matrix.  Store the
             * result in solution and get the square of the cumulative error
             * of the solution.
             */
            int i, j;
            double tmp, diff;
            double sqr_err_r = 0.0;
            double sqr_err_g = 0.0;
            double sqr_err_b = 0.0;
            const double w = 1.80 * 0.25; /* Over-relaxation = 1.8 */

            // we use a red/black checker model of the discretization grid

            // do reds
            for (i = 0; i < matrix.getHeight(); ++i) {
                for (j = i % 2; j < matrix.getWidth(); j += 2) {
                    printf ("/%d,%d", j, i);

                    if ((0 == mask (i, j)) || (i == 0) || (i == (matrix.getHeight() - 1)) || (j == 0) || (j == (matrix.getWidth() - 1))) {
                        // do nothing at the boundary or outside mask
                        solution.r (i, j) = matrix.r (i, j);
                        solution.g (i, j) = matrix.g (i, j);
                        solution.b (i, j) = matrix.b (i, j);
                    } else {
                        // Use Gauss Siedel to get the correction factor then over-relax it
                        tmp = solution.r (i, j);
                        solution.r (i, j) = (matrix.r (i, j) + w *
                                             (
                                                 matrix.r (i,  j - 1) +  // west
                                                 matrix.r (i,  j + 1) +  // east
                                                 matrix.r (i - 1,  j) +  // north
                                                 matrix.r (i + 1,  j) - 4.0 * matrix.r (i, j) // south
                                             )
                                            );

                        diff = solution.r (i, j) - tmp;
                        sqr_err_r += diff * diff;


                        tmp = solution.g (i, j);
                        solution.g (i, j) = (matrix.g (i, j) + w *
                                             (
                                                 matrix.g (i,  j - 1) +  // west
                                                 matrix.g (i,  j + 1) +  // east
                                                 matrix.g (i - 1,  j) +  // north
                                                 matrix.g (i + 1,  j) - 4.0 * matrix.g (i, j) // south
                                             )
                                            );

                        diff = solution.g (i, j) - tmp;
                        sqr_err_g += diff * diff;



                        tmp = solution.b (i, j);
                        solution.b (i, j) = (matrix.b (i, j) + w *
                                             (
                                                 matrix.b (i,  j - 1) +  // west
                                                 matrix.b (i,  j + 1) +  // east
                                                 matrix.b (i - 1,  j) +  // north
                                                 matrix.b (i + 1,  j) - 4.0 * matrix.b (i, j) // south
                                             )
                                            );

                        diff = solution.b (i, j) - tmp;
                        sqr_err_b += diff * diff;

                    }
                }
            }


            /* Do blacks
            *
            * As we've done the reds earlier, we can use them right now to
            * accelerate the convergence. So we have "solution" in the solver
            * instead of "matrix" above
            */
            for (i = 0; i < matrix.getHeight(); i++) {
                for (j = (i % 2) ? 0 : 1; j < matrix.getWidth(); j += 2) {
                    printf (":%d,%d", j, i);

                    if ((0 == mask (i, j)) || (i == 0) || (i == (matrix.getHeight() - 1)) || (j == 0) || (j == (matrix.getWidth() - 1))) {
                        // do nothing at the boundary or outside mask
                        solution.r (i, j) = matrix.r (i, j);
                        solution.g (i, j) = matrix.g (i, j);
                        solution.b (i, j) = matrix.b (i, j);
                    } else {
                        // Use Gauss Siedel to get the correction factor then over-relax it
                        tmp = solution.r (i, j);
                        solution.r (i, j) = (matrix.r (i, j) + w *
                                             (
                                                 matrix.r (i,  j - 1) +  // west
                                                 matrix.r (i,  j + 1) +  // east
                                                 matrix.r (i - 1,  j) +  // north
                                                 matrix.r (i + 1,  j) - 4.0 * matrix.r (i, j) // south
                                             )
                                            );

                        diff = solution.r (i, j) - tmp;
                        sqr_err_r += diff * diff;



                        tmp = solution.g (i, j);
                        solution.g (i, j) = (matrix.g (i, j) + w *
                                             (
                                                 matrix.g (i,  j - 1) +  // west
                                                 matrix.g (i,  j + 1) +  // east
                                                 matrix.g (i - 1,  j) +  // north
                                                 matrix.g (i + 1,  j) - 4.0 * matrix.g (i, j) // south
                                             )
                                            );

                        diff = solution.g (i, j) - tmp;
                        sqr_err_g += diff * diff;



                        tmp = solution.b (i, j);
                        solution.b (i, j) = (matrix.b (i, j) + w *
                                             (
                                                 matrix.b (i,  j - 1) +  // west
                                                 matrix.b (i,  j + 1) +  // east
                                                 matrix.b (i - 1,  j) +  // north
                                                 matrix.b (i + 1,  j) - 4.0 * matrix.b (i, j) // south
                                             )
                                            );

                        diff = solution.b (i, j) - tmp;
                        sqr_err_b += diff * diff;
                    }
                }
            }

            // ----------------------------------------------------------------

            // copy solution to matrix
            solution.copyData (&matrix);

            if (sqr_err_r < EPSILON && sqr_err_g < EPSILON && sqr_err_b < EPSILON) {
                break;
            }

            printf ("\n>>> n=#%d\n", n);
        }

        printf ("\n");
#endif




        // add solution to original image and store in tempPR
        for     (int i = 0, i2 = dst_YMin; i2 < dst_YMax - 1; ++i, ++i2) {
            if (i2 < 0 || i2 >= img->getHeight()) {
                continue;
            }

            for (int j = 0, j2 = dst_XMin; j2 < dst_XMax - 1; ++j, ++j2) {
                if (j2 < 0 || j2 >= img->getWidth()) {
                    continue;
                }

                //float c2 = float (mask (i, j)) / 255.f;
                //float c1 = 1.f - c2;
                //resultPR->r(i,j) = (unsigned char) CLAMP0255 ( ROUND( double(first->r(i,j)) + double(secondPR->r(i,j)) ) );


                img->r (i2, j2) = 65535.0f; //img->r(i2,j2)*c1 + srcBuff->r(i,j)*c2;
                img->g (i2, j2) = 0.0f; //img->g(i2,j2)*c1 + srcBuff->g(i,j)*c2;
                img->b (i2, j2) = 0.0f; //img->b(i2,j2)*c1 + srcBuff->b(i,j)*c2;
                /*
                img->r(i2,j2) = img->r(i2,j2)*c1 + (solution.r(i,j) + srcBuff->r(i,j))*c2;
                img->g(i2,j2) = img->g(i2,j2)*c1 + (solution.g(i,j) + srcBuff->g(i,j))*c2;
                img->b(i2,j2) = img->b(i2,j2)*c1 + (solution.b(i,j) + srcBuff->b(i,j))*c2;
                */
            }
        }

    }
}
#endif

}

