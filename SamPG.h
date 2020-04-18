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

    std::vector<Tree*> samplesForest; // pointers to trees used as samples
    std::vector<bool> already_drawn; // vector of boolean values, if a node has been already drawn the corresponding
                                    // value is setted to true, otherwise is setted to false
    std::vector<std::vector<int>> counters; // matrix of counters
    BestRandom* random_roots; // set of integer values used to randomly choose nodes, not already used, as tree root
    int num_samples; // number of samples used
    int num_counters; // number of counters used
    int totalNodes; // total number of nodes in the forest
    bool EndGeneration; // boolean value used to understand if no more new trees are needed
    int last_node_checked; // integer value that keeps the ID of last node used to generate labels
                          // (it's used to boost the maxDescNode function)
    NetworKit::Graph* graph; // Graph


public:



    SamPG();

    /**
     *
     * @param k = number of samples to create at the beginning
     * @param c = number of counters to use
     * @param graph
     */

    SamPG(int k, int c, NetworKit::Graph* graph);

    /**
     *
     * @return the number of total nodes in the set of samples
     *
     */
    int getTotalNodes();

    /**
     * This function creates the forest of tree samples using the parameters setted in the constructor.
     */

    void createForest();

    /**
     * @return the ID of node with max number of descendants in the set of samples.
     */

    int maxDescNode();

    /**
     *
     * @param n_samples = number of new tree to introduce in the forest
     * @param labeling = set of labels used in the Pruned Dijkstra's algorithm
     *
     * This function creates n_samples of new trees using the Pruned Dijkstra's algorithm and so pruning the vertex
     * already covered with previous labels.
     */
	
    void encreaseForest(int n_samples,Labeling* labeling);

    /**
     *
     * @param node = ID of node to eliminate from all the samples of forest
     *
     * This function removes from all the samples the subtrees rooted on the node corresponding to the passed ID.
     */

    void updateForest(int node);

    /**
     *
     * @return true if the random trees generation is ended (so it has reached the max number of trees setted), false
     * otherwise.
     */
    bool isEnded();

    /**
     *
     * @return number of samples in the forest.
     */
    int getNumSamples();

};


#endif //SAMPG_H
