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
    treeDijkstra->direct_acc[source]= treeDijkstra->getRoot();
    //inserimento radice nel vettore handles
    handles[source] = pq->emplace(heap_dijkstra(source,0));
    // variabile per peso arco
    int wuv;
    //inizio ciclo per dijkstra
    while(pq->empty() == false){
        // estraggo il nodo  v con chiave minima
        int current = pq->top().node;
        //estraggo la distanza dalla sorgente s al nodo u
        int distance = pq->top().prio;
        //rimozione elemento dallo heap
        pq->pop();
        if(pruned && current!=source){
            if(index->query(source, current) <= distance){
                //std::cout<<"QUERY: "<<index->query(source, current)<<"\n";
                //std::cout<<"Distanza: "<<distance<<"\n";
                continue;
            }
        }
        graph->forNeighborsOf(current, [&](int v) {
            wuv = graph -> weight(current,v);
            if(distances[v] == INF) {
                distances[v] = distance + wuv;
                handles[v] = pq->push(heap_dijkstra(v, distances[v]));
                // il nodo è stato trovato per la prima volta quindi faccio addNode poi lo aggiungo alla tabella hash
                treeNode* new_child= new treeNode();
                new_child->node=v;
                treeDijkstra->addNode(new_child,treeDijkstra->direct_acc[current]);
                treeDijkstra->direct_acc[v]=new_child;
            }
            else if(wuv+distance< distances[v]) {
                distances[v] = distance + wuv;
                pq->erase(handles[v]);
                handles[v] = pq->push(heap_dijkstra(v, distances[v]));
                //aggiorno il padre di v e rimuovo il figlio v dalla lista del padre precedemente salvato;
                treeDijkstra->updateFather(treeDijkstra->direct_acc[current],treeDijkstra->direct_acc[v]);
            }
        });
    }
    treeDijkstra->computeDescendants(treeDijkstra->getRoot());
}