//
// Created by f4b3r on 31/03/20.
//

#include "SamPG.h"

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
    std::vector<int> counters(this->samplesForest[0]->getDescVect().size(),0);
    for (int i = 0; i < this->num_samples ; ++i) {
        for(int j = 0; j < this->samplesForest[0]->getDescVect().size(); j++){
            counters[j] += this->samplesForest[i]->getDescVect()[j];
        }
    }

    int roundNode = 0;
    for(int y=1;y<counters.size();t++){
        if(counters[y-1]<counters[y]){
            roundNode=y;
        }
    }

	return roundNode;

}

// Funzione per incrementare il numero di alberi campionati.  | Serve una versione di dijkstra modificata |

void SamPG::encreaseForest(int Samples, NetworKit::Graph* graph, Labeling* index) {


      MersenneTwister random;
      Dijkstra *dijkstra;
      for (int i = 0; i < this->Samples; i++) {
          this->samplesForest.push_back(
                  new Tree(random.getRandomInteger() % graph->numberOfNodes(), graph->numberOfNodes()));
          dijkstra->runDijkstra(samplesForest.back(), graph, true, index);
      }
      this->num_samples += Samples;
  }


void updateForest(int node){

	for(int i=0; i<num_samples; i++){
		treeNode* maxNode;
		maxNode=this->sampleForest[i]->DFS(node);
		sampleForest[i]->deleteSubTree(maxNode);

	}
	
}



























