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
#ifndef _CURVELISTENER_
#define _CURVELISTENER_

class CurveEditor;

class CurveListener {

    private:
        bool multi;
    public:
        CurveListener() : multi(false) {}
        virtual ~CurveListener() {}
        virtual void curveChanged () {}
        virtual void curveChanged (CurveEditor* ce) {}
        void setMulti(bool value) { multi = value; }
        bool isMulti() { return multi; }

        /** @brief Ask the reset curve for a given curve type
         * @param ce CurveEditor that we want to reset
         * @param curve Actual curve for the return value. The actual curve type (given by the first value of the vector)
         *              should be kept the same. Change the curve type if REALLY necessary! */
        virtual bool getResetCurve(CurveEditor *ce, std::vector<double> &curve) { return false; };

        /** @brief Blend pipette values from its different channels into a single value
        If the buffer has more than one channel and one channel, this method will blend them together.
        @param chan1 first channel's value
        @param chan2 second channel's value
        @param chan3 third channel's value
        @return the blended value */
        virtual float blendPipetteValues(float chan1, float chan2, float chan3) {
            float retVal = 0.f;
            int n = 0;
            if (chan1 != -1.f) {
                retVal += chan1;
                ++n;
            }
            if (chan2 != -1.f) {
                retVal += chan2;
                ++n;
            }
            if (chan3 != -1.f) {
                retVal += chan3;
                ++n;
            }
            if (n>1)
                retVal /= n;
            else if (!n)
                retVal = -1.f;
            return retVal;
        }
};

#endif
