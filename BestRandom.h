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

    BestRandom(int min, int max);

    int random();

};


#endif //BESTRANDOM_H
