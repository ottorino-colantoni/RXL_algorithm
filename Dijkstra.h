//
// Created by f4b3r on 31/03/20.
//

#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include "Tree.h"
#include <networkit/graph/Graph.hpp>

struct heap_data {
    int node;
    int prio;
    heap_data(int n,int w)
    {
        node = n;
        prio = w;
    }

    bool operator<(heap_data const & rhs) const
    {
        //MIN HEAP WITHOUT CHANGIN SIGN
        return this->prio > rhs.prio;
    }
};

class Dijkstra {

public:

    Dijkstra();

    void runDijkstra(Tree* treeDijkstra, NetworKit::Graph* graph);

};

#endif //DIJKSTRA_H


