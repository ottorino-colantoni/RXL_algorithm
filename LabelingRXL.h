//
// Created by f4b3r on 09/04/20.
//

#ifndef LABELINGRXL_H
#define LABELINGRXL_H



struct heap_Labeling {
    int node;
    int prio;
    heap_Labeling(int n,int w)
    {
        node = n;
        prio = w;
    }

    bool operator<(heap_Labeling const & rhs) const
    {
        //MIN HEAP WITHOUT CHANGIN SIGN
        return this->prio > rhs.prio;
    }
};

class LabelingRXL {

private:

    std::vector<int> order;
    std::vector<int> reverse_order;
    std::vector<std::vector<std::pair<int,int>>> Labels;

public:


    LabelingRXL();

    std::vector<std::vector<std::pair<int,int>>> getLabels();

    void createLabels(int node, NetworKit::Graph* graph);

    int query(int node1, int node2);

    int getNumberOfLabels();

    void printLabels();

    int node_to_index(int node);

    int index_to node(int index);



};


#endif //LABELINGRXL_H
