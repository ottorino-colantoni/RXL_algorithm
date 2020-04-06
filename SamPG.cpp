//
// Created by f4b3r on 31/03/20.
//

#include "SamPG.h"
#include "mytimer.h"


SamPG::SamPG(){};

SamPG::SamPG(int k){
    this->num_samples = k;
}

void SamPG::createForest(NetworKit::Graph* graph){
    this->random_roots= new BestRandom(0, graph->numberOfNodes());
    Dijkstra* dijkstra;
    for(int i = 0; i<this->num_samples; i++){
        Tree* new_tree = new Tree(this->random_roots->random(), graph->numberOfNodes());
        this->samplesForest.push_back(new_tree);
        dijkstra->runDijkstra(samplesForest.back(), graph);
        this->totalNodes += new_tree->desc_vect[new_tree->getRoot()->node]+1;
    }
}

int SamPG::maxDescNode(){
    int counter=0;
    int roundNode;
    int max = 0;
    for (int i = 0; i < this->samplesForest[0]->desc_vect.size() ; i++) {
        for(int j = 0; j < this->num_samples; j++){
            counter += this->samplesForest[j]->desc_vect[i];
        }
        if(counter>=max){
            max = counter;
            roundNode = i;
        }
        counter = 0;
    }
    return roundNode;

}

// Funzione per incrementare il numero di alberi campionati.  | Serve una versione di dijkstra modificata |

void SamPG::encreaseForest(int Samples, NetworKit::Graph* graph, Labeling* index) {

      Dijkstra *dijkstra;
      for (int i = 0; i < Samples; i++) {
          Tree* new_tree = new Tree(this->random_roots->random(), graph->numberOfNodes());
          this->samplesForest.push_back(new_tree);
          dijkstra->runDijkstra(samplesForest.back(), graph, true, index);
          this->totalNodes += new_tree->desc_vect[new_tree->getRoot()->node]+1;
         // std::cout<<"Albero generato: \n";
          //this->samplesForest.back()->printTree(this->samplesForest.back()->getRoot());
      }
      this->num_samples += Samples;
  }


void SamPG::updateForest(int node){

    treeNode* maxNode;
        for (int i = 0; i < num_samples; i++) {
            maxNode = this->samplesForest[i]->direct_acc[node];
            //std::cout<<"PRIMA"<<"\n";
            //samplesForest[i]->printTree(samplesForest[i]->getRoot());
            this->totalNodes -= samplesForest[i]->desc_vect[node];

            samplesForest[i]->deleteSubTree(maxNode);
            //std::cout<<"DOPO"<<"\n";
            //samplesForest[i]->printTree(samplesForest[i]->getRoot());
        }

}

int SamPG::getTotalNodes() {
    return this->totalNodes;
}



























