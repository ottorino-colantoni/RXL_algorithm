//
// Created by f4b3r on 31/03/20.
//

#include "SamPG.h"

SamPG::SamPG(){};

SamPG::SamPG(int k){
    this->num_samples = k;
}

void SamPG::createForest(NetworKit::Graph* graph){
    MersenneTwister random;
    Dijkstra* dijkstra;
    for(int i = 0; i<this->num_samples; i++){
        this->samplesForest.push_back(new Tree(random.getRandomInteger() % graph->numberOfNodes(), graph->numberOfNodes()));
        dijkstra->runDijkstra(samplesForest.back(), graph);
    }
}

int SamPG::maxDescNode(){
    std::vector<int> counters(this->samplesForest[0]->getDescVect().size(),0);
    for (int i = 0; i < this->num_samples ; ++i) {
        for(int j = 0; j < this->samplesForest[0]->getDescVect().size(); j++){
            counters[j] += this->samplesForest[i]->getDescVect()[j];
        }
    }
    return std::max_element(counters.begin(), counters.end());
}

