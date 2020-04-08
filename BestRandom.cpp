//
// Created by f4b3r on 06/04/20.
//

#include "BestRandom.h"
#include <iostream>

BestRandom::BestRandom() {}

BestRandom::BestRandom(int min, int max) {
    for (int i = min; i < max; ++i) {
        this->values.push_back(i);
    }
}

int BestRandom::random() {

    MersenneTwister rand;
    int position = rand.getRandomInteger() % this->values.size();
    int random_value = this->values[position];
    this->values[position] = this->values.back();
    this->values.resize(this->values.size()-1);
    return random_value;

}