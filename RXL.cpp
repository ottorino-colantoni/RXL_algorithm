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

    int k = 10;
    NetworKit::Graph* graph;
    Auxiliary::read("graph1.hist", false, &graph);

    SamPG* spg = new SamPG(k);
    spg->createForest(graph);

    int max = spg->maxDescNode();

    std::cout << "nodo con più discendenti " << max <<"\n"; 

    custom_node new_node = max;
	
    std::pair<std::vector<custom_node>,std::vector<custom_node>> keeper;

    keeper.first.push_back(new_node);
	
    keeper.second.resize(graph->upperNodeIdBound());

    keeper.second[max] = 0;

    Labeling* labeling = new Labeling(graph->isDirected());

    Labeling_Tools* lt = new Labeling_Tools(graph,labeling,keeper);
	

    std::cout<<"numero di label dopo prima iterazione :"<<labeling->getNumberOfLabelEntries()<< "\n";

    labeling->printInLabels();
        

	
	




    

}