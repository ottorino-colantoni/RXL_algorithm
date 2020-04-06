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
    MersenneTwister random;
    Dijkstra* dijkstra;
    for(int i = 0; i<this->num_samples; i++){
        this->samplesForest.push_back(new Tree(random.getRandomInteger() % graph->numberOfNodes(), graph->numberOfNodes()));
        dijkstra->runDijkstra(samplesForest.back(), graph);
    }
}

int SamPG::maxDescNode(){
    std::vector<int> counters(this->samplesForest[0]->desc_vect.size(),0);
    int roundNode;
    int max = 0;
    mytimer t;
    t.restart();
    for (int i = 0; i < this->samplesForest[0]->desc_vect.size() ; i++) {
        for(int j = 0; j < this->samplesForest.size(); j++){
            counters[i] += this->samplesForest[j]->desc_vect[i];
        }
        if(counters[i]>=max){
            max = counters[i];
            roundNode = i;
        }
    }
    std::cout<<"tempo double for max: "<<t.elapsed();
    return roundNode;

}

// Funzione per incrementare il numero di alberi campionati.  | Serve una versione di dijkstra modificata |

void SamPG::encreaseForest(int Samples, NetworKit::Graph* graph, Labeling* index) {


      MersenneTwister random;
      Dijkstra *dijkstra;
      for (int i = 0; i < Samples; i++) {
          this->samplesForest.push_back(
                  new Tree(random.getRandomInteger() % graph->numberOfNodes(), graph->numberOfNodes()));
          dijkstra->runDijkstra(samplesForest.back(), graph, true, index);
      }
      this->num_samples += Samples;
  }


void SamPG::updateForest(int node){

    treeNode* maxNode;
        for (int i = 0; i < num_samples; i++) {
            maxNode = this->samplesForest[i]->direct_acc[node];
            //std::cout<<"PRIMA"<<"\n";
            //samplesForest[i]->printTree(samplesForest[i]->getRoot());
            samplesForest[i]->deleteSubTree(maxNode);
            //std::cout<<"DOPO"<<"\n";
            //samplesForest[i]->printTree(samplesForest[i]->getRoot());
        }

}



























