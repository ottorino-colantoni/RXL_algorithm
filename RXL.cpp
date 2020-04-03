//
// Created by f4b3r on 03/04/20.
//


#include "Labeling.h"
#include "Auxiliary.h"
#include "CustomDataTypes.h"
#include "Labeling_Tools.h"
#include "Dijkstra.h"
#include "SamPG.h"
#include "Tree.h"
#include "networkit/graph/Graph.hpp"



int main(){

    int k = 50;
    NetworKit::Graph* graph;
    Auxiliary::read("example.hist", false, &graph);

    SamPG* spg = new SamPG(k);
    spg->createForest(graph);

    for(int i = 0; i<graph->numberOfNodes(); i++){

        int max = spg->maxDescNode();



    }

}