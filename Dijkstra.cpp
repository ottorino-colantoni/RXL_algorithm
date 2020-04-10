#include "boost/heap/fibonacci_heap.hpp"
#include <vector>
#include "Dijkstra.h"
#define INF 0x3F3F3F



Dijkstra::Dijkstra() {};


void Dijkstra::runDijkstra(Tree* treeDijkstra, NetworKit::Graph* graph, bool pruned , Labeling* index ) {

    boost::heap::fibonacci_heap<heap_dijkstra>* pq = new boost::heap::fibonacci_heap<heap_dijkstra>();
    boost::heap::fibonacci_heap<heap_dijkstra>::handle_type* handles =
            new boost::heap::fibonacci_heap<heap_dijkstra>::handle_type[graph->upperNodeIdBound()];
    std::vector<int> distances;
    int size= graph->numberOfNodes();
    distances.resize(size,INF);
    int source = treeDijkstra->getRoot()->ID;
    distances[source] = 0;
    handles[source] = pq->push(heap_dijkstra(treeDijkstra->getRoot(),0));
    int wuv;
    while(pq->empty() == false){
        treeNode* current = pq->top().node;
        int distance = pq->top().prio;
        pq->pop();
        graph->forNeighborsOf(current->ID, [&](int v) {
            wuv = graph -> weight(current->ID,v);
            if(distances[v] == INF) {
                if(!pruned || (pruned && index->query(source,v)>(distance+wuv))) {
                    distances[v] = distance + wuv;
                    treeNode *new_child = new treeNode();
                    new_child->ID = v;
                    handles[v] = pq->push(heap_dijkstra(new_child, distances[v]));
                    treeDijkstra->addNode(new_child, current);
                }
            }
            else if(wuv+distance< distances[v]) {
                distances[v] = distance + wuv;
                (*handles[v]).prio = distances[v];
                pq->decrease(handles[v]);
                treeDijkstra->updateFather(current,(*handles[v]).node);
            }
        });
    }
}