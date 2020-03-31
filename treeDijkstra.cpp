#include "Auxiliary.h"
#include <boost/program_options.hpp>
#include "Tree.h"
#include "boost/heap/fibonacci_heap.hpp"
#include <omp.h>
#include "mersenneTwister.h"
#include <string>
#include <iostream>
#include <vector>
#define INF 0x3F3F3F
#define PATH_FILE "log.csv"

int numDecreaseKey = 0;

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


int main(int argc, char** argv) {

    //declare supported options
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");

    desc.add_options()
            ("graph_location,g", po::value<std::string>(), "Input Graph File Location (in case of -a 1)");


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify( vm);

    int ne = -1;
    std::string graph_location;

    if (vm.empty()){
        std::cout << desc << "\n";
        throw std::runtime_error("Empty options");
    }

    if (vm.count("graph_location"))
        graph_location = vm["graph_location"].as<std::string>();


    if(graph_location == ""){
        std::cout << desc << "\n";
        throw std::runtime_error("wrong graph_location");
    }

    MersenneTwister random;
    NetworKit::Graph* graph;
    Auxiliary::read(graph_location,0,&graph); //creazione del grafo leggendo il documento undirected.hist
    boost::heap::fibonacci_heap<heap_data>* pq = new boost::heap::fibonacci_heap<heap_data>();
    boost::heap::fibonacci_heap<heap_data>::handle_type* handles =
            new boost::heap::fibonacci_heap<heap_data>::handle_type[graph->upperNodeIdBound()];
    std::vector<int> distances;
    int source = random.getRandomInteger() % graph->numberOfNodes();
    distances.resize(graph->numberOfNodes(),INF);
    distances[source] = 0;
    //creazione albero
    Tree* treeDijkstra= new Tree(source);
    //creazione array
    int size= graph->numberOfNodes();
    std::vector<treeNode*> support;
    support.resize(size);
    //inserimento radice nell'array
    support[source]= treeDijkstra->getRoot();
    //inserimento radice nel vettore handles
    handles[source] = pq->emplace(heap_data(source,0));
    // variabile per peso arco
    int wuv;
    //inizio ciclo per dijkstra
    while(pq->empty() == false){

        std::cout << "dimensione array "<< support.size()<< "\n";
        // estraggo il nodo  v con chiave minima
        int current = pq->top().node;
        std::cout << "nodo estratto"<< current <<"\n";
        //estraggo la distanza dalla sorgente s al nodo u
        int distance = pq->top().prio;
        std::cout << distance <<"\n";
        //rimozione elemento dallo heap
        pq->pop();
        graph->forNeighborsOf(current, [&](int v) {
            wuv = graph -> weight(current,v);
            if(distances[v] == INF) {
                distances[v] = distance + wuv;
                handles[v] = pq->push(heap_data(v, distances[v]));
                std::cout << "nodo estratto"<< support[current]->node <<"\n";
                std::cout << "vicino di current  "<< v << "\n";
                std::cout << "distanza del nuovo nodo  "<< distances[v] << "\n";
                // il nodo è stato trovato per la prima volta quindi faccio addNode poi lo aggiungo alla tabella hash
                std::cout << "questo è il valore di current  "<< current << "\n";
                std::cout << "nodo estratto dallo heap nella stessa iterazione  "<< support[current]->node << "\n";
                treeNode* new_child= new treeNode();
                new_child->node=v;
                treeDijkstra->addNode(new_child,support[current]);
                support[v]=new_child;
                std::cout << "nodo appena creato  "<< support[v]->node << "\n";
                std::cout << "padre nodo appena creato  "<< support[v]->father->node << "\n";
                std::cout << "valore nodo current dopo addnode  "<< support[current]->node << "\n";
            }
            else if(wuv+distance< distances[v]) {
                distances[v] = distance + wuv;
                std::cout << "questa è la nuova distanza" << distances[v] << "\n";
                pq->erase(handles[v]);
                handles[v] = pq->push(heap_data(v, distances[v]));
                std::cout << "nuova distanza nello heap"<< (*handles[v]).prio << "\n";
                //aggiorno il padre di v e rimuovo il figlio v dalla lista del padre precedemente salvato;
                std::cout << support[v]->father->node<<"\n";
                treeDijkstra->updateFather(support[current],support[v]);
                std::cout << support[v]->father->node<<"\n";
            }
        });
    }
    std::cout<<"questa è la struttura dell'albero"<< "\n";
    treeDijkstra->printTree(treeDijkstra->getRoot());

    return EXIT_SUCCESS;

}