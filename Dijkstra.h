#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include "Tree.h"
#include <networkit/graph/Graph.hpp>
#include "Labeling.h"


struct heap_dijkstra {
    treeNode* node;
    int prio;
    heap_dijkstra(treeNode* n,int w)
    {
        node = n;
        prio = w;
    }

    bool operator<(heap_dijkstra const & rhs) const
    {
        //MIN HEAP WITHOUT CHANGIN SIGN
        return this->prio > rhs.prio;
    }
};


class Dijkstra {

public:

    Dijkstra();

    /**
     *
     * @param treeDijkstra = tree's pointer to be built
     * @param graph  = graph's pointer
     * @param pruned = if true the function will realize the Pruned Dijkstra, if false otherwise
     * @param index = pointer to the set of labels used to prune the vertices
     *
     * This function realizes the Dijkstra's algorithm and, moreover, it creates on the treeDikstra the resulting shortest
     * paths tree. Moreover, if pruned is setted to true and a set of labels is passed, it realizes the pruned Dijkstra's
     * algorithm to build a pruned tree.
     */

    void runDijkstra(Tree* treeDijkstra, NetworKit::Graph* graph, bool pruned = false, Labeling* index = nullptr);

};

#endif //DIJKSTRA_H


