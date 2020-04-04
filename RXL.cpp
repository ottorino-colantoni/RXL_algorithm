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




int main(){

    int k = 15;
    NetworKit::Graph* graph;
    Auxiliary::read("graph1.hist", false, &graph);

    SamPG* spg = new SamPG(k);
    spg->createForest(graph);
    std::pair <std::vector<custom_node>, std::vector<custom_node>> keeper;
    keeper.second.resize(graph->upperNodeIdBound());
    Labeling *labeling = new Labeling(graph->isDirected());
    Labeling_Tools *lt = new Labeling_Tools(graph, labeling, keeper);
    int max;


    for (int i = 0; i < graph->numberOfNodes(); ++i) {

        max = spg->maxDescNode();
        std::cout << "nodo con più discendenti " << max << "\n";

        lt->add_node_to_keeper(max, i);

        lt->weighted_build_RXL();

        spg->updateForest(max);

    }


    std::cout<<"numero di label dopo prima iterazione :"<<labeling->getNumberOfLabelEntries()<< "\n";

    labeling->printInLabels();
        

	
	




    

}