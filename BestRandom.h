//
// Created by f4b3r on 06/04/20.
//

#ifndef BESTRANDOM_H
#define BESTRANDOM_H

#include "mersenneTwister.h"


class BestRandom {

private:

    std::vector<int> values;

public:

    BestRandom();

    /**
     *
     * @param min = lower bound of set
     * @param max = upper bound of set
     */

    BestRandom(int min, int max);

    /**
     *
     * @return a integer value randomly extracted from the remaining instances of set.
     */

    int random();

};


#endif //BESTRANDOM_H
