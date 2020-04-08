//
// Created by f4b3r on 03/04/20.
//


#include "Labeling.h"
#include "Auxiliary.h"
#include "CustomDataTypes.h"
#include "Labeling_Tools.h"
#include "Dijkstra.h"
#include "SamPG.h"
#include "networkit/graph/Graph.hpp"
#include <stdio.h>
#include "mytimer.h"
#include <omp.h>



int main(){
    int k = 50;
	int c = 16;
    NetworKit::Graph *graph;
    Auxiliary::read("graph.hist", false, &graph);
    SamPG *spg = new SamPG(k,c);
    spg->createForest(graph);
	

  

    std::pair <std::vector<custom_node>, std::vector<custom_node>> keeper;
    keeper.second.resize(graph->upperNodeIdBound());
    Labeling *labeling = new Labeling(graph->isDirected());
    Labeling_Tools *lt = new Labeling_Tools(graph, labeling, keeper);
    int max;
	int maxprec;

    for (int i = 0; i < graph->numberOfNodes()-k; i++) {
        max = spg->maxDescNode();
        std::cout << "nodo con più discendenti " << max << "\n";
        lt->add_node_to_keeper(max, i);
        lt->weighted_build_RXL();
        spg->updateForest(max);

        if(i % 1 == 0) {
            spg->encreaseForest(1, graph, labeling);
        }
        std::cout<<i<<"-";
    }

    std::cout << "numero di label  :" << labeling->getNumberOfLabelEntries() << "\n";
    std::cout << "distanza nodo 3-0 :" << labeling->query(3, 1) << "\n";

}

