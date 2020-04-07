#include <boost/program_options.hpp>
#include "boost/heap/fibonacci_heap.hpp"
#include <omp.h>
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
    int source = treeDijkstra->getRoot()->node;
    distances[source] = 0;
    //creazione array
    //inserimento radice nell'array
    //inserimento radice nel vettore handles
    handles[source] = pq->push(heap_dijkstra(treeDijkstra->getRoot(),0));
    // variabile per peso arco
    int wuv;
    //inizio ciclo per dijkstra
    while(pq->empty() == false){
        // estraggo il nodo  v con chiave minima
        treeNode* current = pq->top().node;
        //estraggo la distanza dalla sorgente s al nodo u
        int distance = pq->top().prio;
        //rimozione elemento dallo heap
        pq->pop();
        if(pruned && current->node!=source){
            if(index->query(source, current->node) <= distance){
                //std::cout<<"QUERY: "<<index->query(source, current)<<"\n";
                //std::cout<<"Distanza: "<<distance<<"\n";
                continue;
            }
        }
        graph->forNeighborsOf(current->node, [&](int v) {
            wuv = graph -> weight(current->node,v);
            if(distances[v] == INF) {
                distances[v] = distance + wuv;
                // il nodo è stato trovato per la prima volta quindi faccio addNode poi lo aggiungo alla tabella hash
                treeNode* new_child = new treeNode();
                new_child->node=v;
                handles[v] = pq->push(heap_dijkstra(new_child, distances[v]));
                treeDijkstra->addNode(new_child,current);
            }
            else if(wuv+distance< distances[v]) {
                distances[v] = distance + wuv;
                (*handles[v]).prio = distances[v];
                pq->decrease(handles[v]);
                //aggiorno il padre di v e rimuovo il figlio v dalla lista del padre precedemente salvato;
                treeDijkstra->updateFather(current,(*handles[v]).node);
            }
        });
    }
}