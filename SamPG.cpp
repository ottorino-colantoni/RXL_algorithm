//
// Created by f4b3r on 31/03/20.
//

#include "SamPG.h"
#include "mytimer.h"


SamPG::SamPG(){};

SamPG::SamPG(int k, int c){
    this->num_samples = k;
	this->num_counters= c;
	this->counters.resize(c);
	this->zero=false;
	this->totalNodes = 0;
}

void SamPG::createForest(NetworKit::Graph* graph){
    this->random_roots= new BestRandom(0, graph->numberOfNodes());
	for(int i=0; i<this->num_counters;i++){
		this->counters[i].resize(graph->numberOfNodes(),0);
	}
    Dijkstra* dijkstra;
    for(int i = 0; i<this->num_samples; i++){
        Tree* new_tree = new Tree(this->random_roots->random(), graph->numberOfNodes());
        this->samplesForest.push_back(new_tree);
        dijkstra->runDijkstra(samplesForest.back(), graph);
		new_tree->computeDescendants(new_tree->getRoot(),this->counters, i, this->num_counters);
		this->totalNodes += new_tree->getRoot()->num_of_descendants+1;
    }
}

int SamPG::maxDescNode(){
  
   int useless_value1=0;
   int max=0;
   int temporarymax=0;
   int roundNode=0;
   int roundNode_for_zero=0;
   int max_for_zero=0;

   for(int i = 0; i< counters[0].size(); i++){

		for(int j=0; j <this->num_counters; j++){

		    temporarymax += this->counters[j][i];
		    if(this->counters[j][i]>useless_value1){useless_value1= counters[j][i];}

		}

		if(temporarymax> max_for_zero){
	        max_for_zero=temporarymax;
	        roundNode_for_zero=i;
	    }

	    temporarymax-=useless_value1;
	    temporarymax = (temporarymax/(this->num_counters-1));
	    if(temporarymax> max){
	        max=temporarymax;
	        roundNode=i;
	    }

	temporarymax=0;
	useless_value1=0;
	}


  if(roundNode == 0 ){this->zero=true;}
  if(zero && roundNode == 0){return roundNode_for_zero;}
  else{return roundNode;}

}

// Funzione per incrementare il numero di alberi campionati.  | Serve una versione di dijkstra modificata |

void SamPG::encreaseForest(int Samples, NetworKit::Graph* graph, Labeling* index) {

      Dijkstra *dijkstra;
      for (int i = 0; i < Samples; i++) {
          Tree* new_tree = new Tree(this->random_roots->random(), graph->numberOfNodes());
          this->samplesForest.push_back(new_tree);
          dijkstra->runDijkstra(this->samplesForest.back(), graph, true, index);
		  int j = this->num_samples + i;
		  new_tree->computeDescendants(new_tree->getRoot(),counters, j, this->num_counters);
		  this->totalNodes += new_tree->getRoot()->num_of_descendants + 1;
      }
      this->num_samples += Samples;
  }


void SamPG::updateForest(int node){

    treeNode* maxNode;
    for (int i = 0; i < samplesForest.size(); i++) {
        maxNode = this->samplesForest[i]->DFS(node);
        samplesForest[i]->deleteSubTree(maxNode, this->counters, i, this->num_counters);
    }
}

int SamPG::getTotalNodes() {
    return this->totalNodes;
}


























