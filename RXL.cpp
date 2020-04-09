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
    int k = 16;
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

    for (int i = 0; i < graph->numberOfNodes(); i++) {
        max = spg->maxDescNode();
        if(max == 54){
            int ciao = 0;
        }
        std::cout << "nodo con più discendenti " << max << "\n";
        lt->add_node_to_keeper(max, i);
        lt->weighted_build_RXL();
        spg->updateForest(max);
        if(spg->getTotalNodes()<10*k*graph->numberOfNodes()) {
            std::cout<<"SUUUUIIIIIIIIIIII"<<" numero nodi: "<<spg->getTotalNodes()<<"\n";
            spg->encreaseForest(1, graph, labeling);
            for(int i = 0; i < c; i++){
                if(spg->counters[i][54]>0){
                    bool PD = true;
                }
            }
        }
        std::cout<<i<<"-";
    }

    std::cout << "numero di label  :" << labeling->getNumberOfLabelEntries() << "\n";
    std::cout << "distanza nodo 3-0 :" << labeling->query(3, 1) << "\n";

}

