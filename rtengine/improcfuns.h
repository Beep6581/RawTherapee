/*
 * improcfuns.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef IMPROCFUNS_H_
#define IMPROCFUNS_H_

namespace rtengine {

class ImProcFunctions {

    public:
        static void calcAutoExp (unsigned int* histogram, int histcompr, double expcomp, double clip, double& br, int& bl);

};

}

#endif /* IMPROCFUNS_H_ */
