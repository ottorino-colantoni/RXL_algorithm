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

    int k = 5;
    NetworKit::Graph* graph;
    Auxiliary::read("graph1.hist", false, &graph);

    SamPG* spg = new SamPG(k);
    spg->createForest(graph);
    std::pair <std::vector<custom_node>, std::vector<custom_node>> keeper;
    keeper.second.resize(graph->upperNodeIdBound());
    Labeling *labeling = new Labeling(graph->isDirected());
    Labeling_Tools *lt = new Labeling_Tools(graph, labeling, keeper);
    int max;


    for (int i = 0; i < graph->numberOfNodes(); i++) {

        max = spg->maxDescNode();
        std::cout << "nodo con più discendenti " << max << "\n";
        lt->add_node_to_keeper(max, i);

        lt->weighted_build_RXL();

        spg->updateForest(max);

		spg->encreaseForest(5,graph,labeling);

    }


    std::cout<<"numero di label  :"<<labeling->getNumberOfLabelEntries()<< "\n";
 	std::cout<<"distanza nodo 3-0 :"<<labeling->query(3,0)<< "\n";
    labeling->printInLabels();

	

    /*
    NetworKit::Graph* graph;
    Auxiliary::read("graph1.hist", false, &graph);
    Dijkstra d;
    Tree* t = new Tree(0, graph->numberOfNodes());
    d.runDijkstra(t, graph);
    t->printTree(t->getRoot());
    std::vector<int> c = t->getDescVect();
    for (int i = 0; i < c.size() ; ++i) {
        std::cout<<"nodo: "<<i<<" disc: "<<c[i]<<"\n";
    }

    treeNode* n = t->DFS(3);
    t->deleteSubTree(n);

    t->printTree(t->getRoot());
    c = t->getDescVect();
    for (int i = 0; i < c.size() ; ++i) {
        std::cout<<"nodo: "<<i<<" disc: "<<c[i]<<"\n";
    }
    */



    

}