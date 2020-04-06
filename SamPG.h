//
// Created by f4b3r on 31/03/20.
//

#include "mersenneTwister.h"
#include "Dijkstra.h"
#include <networkit/graph/Graph.hpp>
#include "Labeling.h"
#include "BestRandom.h"

#ifndef SAMPG_H
#define SAMPG_H


class SamPG {

private:

    BestRandom* random_roots;
    int totalNodes;


public:


    int num_samples;
    std::vector<Tree*> samplesForest;

    SamPG();

    SamPG(int k);

    int getTotalNodes();

    void createForest(NetworKit::Graph* graph);

    int maxDescNode();
	
    void encreaseForest(int Samples,NetworKit::Graph* graph,Labeling* labeling);

    void updateForest(int node);

};


#endif //SAMPG_H
