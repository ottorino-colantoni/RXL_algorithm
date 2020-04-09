//
// Created by f4b3r on 09/04/20.
//

#include "LabelingRXL.h"
#define INF 0x3F3F3F


LabelingRXL::LabelingRXL() {}

std::vector<std::vector<std::pair<int,int>>> LabelingRXL::getLabels() {
    return this->Labels;
}

int LabelingRXL::getNumberOfLabels() {
    int sum = 0;
    for (int i = 0; i < this->Labels.size() ; i++) {
        sum += this->Labels[i].size();
    }
    return sum;
}

int LabelingRXL::node_to_index(int node){
    return this->reverse_order[node];
}

int LabelingRXL::index_to_node(int index){
    return this->order[index];
}

void LabelingRXL::printLabels() {

    for (int i = 0; i < this->Labels.size() ; i++) {
        for (int j = 0; j < this->Labels[i].size() ; j++) {
            std::cout<<"Label nodo: "<<i<<" "<<this->Labels[i][j].first<<" Distanza: "<<this->Labels[i][j].second<<"\n";
        }
    }
}

void LabelingRXL::createLabels(int node, NetworKit::Graph *graph) {
    this->reverse_order[node] = this->order.size();
    this->order.push_back(node);
    boost::heap::fibonacci_heap<heap_Labeling>* pq = new boost::heap::fibonacci_heap<heap_Labeling>();
    boost::heap::fibonacci_heap<heap_Labeling>::handle_type* handles =
            new boost::heap::fibonacci_heap<heap_Labeling>::handle_type[graph->upperNodeIdBound()];
    int source = node;
    std::vector<int> distances;
    distances.resize(graph->numberOfNodes(), INF);
    distances[source] = 0;
    handles[source] = pq->push(heap_Labeling(source, 0));

    while (!pq->empty()){
        int current_node = pq->top().node;
        int current_distance = pq->top().prio;
        pq->pop();
        if(query(source, current_node)<=current_distance){
            continue;
        }
        std::pair<int,int> new_label;
        new_label.first = source;
        new_label.second = current_distance;
        this->Labels[current_node].push_back(new_label);
        graph->forNeighborsOf(current->node, [&](int v) {
            int wuv = graph -> weight(current_node,v)
            if(distances[v] == INF){
                distances[v] = distances[current] + wuv;
                handles[v] = pq->push(heap_Labeling(v, distances[v]));
            }
            else if(current_distance + wuv < distances[v]){
                distances[v] = current_distance + wuv;
                (*handles[v]).prio = distances[v];
                pq->decrease(handles[v]);
            }
        }
    }
}

int LabelingRXL::query(int node1, int node2) {

}