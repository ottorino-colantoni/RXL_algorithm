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
    NetworKit::Graph *graph;
    Auxiliary::read("graph.hist", false, &graph);
    SamPG *spg = new SamPG(k);
    spg->createForest(graph);
    std::pair <std::vector<custom_node>, std::vector<custom_node>> keeper;
    keeper.second.resize(graph->upperNodeIdBound());
    Labeling *labeling = new Labeling(graph->isDirected());
    Labeling_Tools *lt = new Labeling_Tools(graph, labeling, keeper);
    int max;

    for (int i = 0; i < graph->numberOfNodes(); i++) {
        mytimer timercounter;

        timercounter.restart();
        max = spg->maxDescNode();
        /*if (max == maxprec){
            break;
        }*/
        std::cout<<"tempo max: "<<timercounter.elapsed()<<"\n";
        //std::cout << "nodo con più discendenti " << max << "\n";
        std::cout<<"iterazione"<<i<<"\n";
        lt->add_node_to_keeper(max, i);
        timercounter.restart();
        lt->weighted_build_RXL();
        std::cout<<"tempo label: "<<timercounter.elapsed()<<"\n";
        //labeling->printInLabels();
        timercounter.restart();
        spg->updateForest(max);
        std::cout<<"tempo update: "<<timercounter.elapsed()<<"\n";
        std::cout<<spg->getTotalNodes()<<"\n";
        if(i % 2 == 0) {
            timercounter.restart();
            spg->encreaseForest(1, graph, labeling);
            std::cout<<"tempo encrease: "<<timercounter.elapsed()<<"\n";
        }
    }


    std::cout << "numero di label  :" << labeling->getNumberOfLabelEntries() << "\n";
    std::cout << "distanza nodo 3-0 :" << labeling->query(1, 0) << "\n";
    //labeling->printInLabels();


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