//
// Created by f4b3r on 31/03/20.
//

#include "SamPG.h"
#include "mytimer.h"


SamPG::SamPG(){};

SamPG::SamPG(int k, int c){
    this->num_samples = k;
	this->num_counters= c;
	counters.resize(c);
	
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
       // this->totalNodes += new_tree->desc_vect[new_tree->getRoot()->node]+1;
    }
}

int SamPG::maxDescNode(){
  
   int useless_value1=0;
   int max=0;
   int temporarymax=0;
   int roundNode=0;

   for(int i = 0; i< counters[0].size(); i++){

		for(int j=0; j <this->num_counters; j++){

		    temporarymax += counters[j][i];
		    if(counters[j][i]>useless_value1){useless_value1= counters[j][i];}

		}
	    temporarymax-=useless_value1;
	    temporarymax = (temporarymax/(this->num_counters-1));
	    std::cout<<"tempMax: "<<temporarymax<<"\n";
	    if(temporarymax> max){
	        max=temporarymax;
	        roundNode=i;
	    }

	temporarymax=0;
	useless_value1=0;

	}

   assert(roundNode>=0 && roundNode<counters[0].size());
   std::cout<<roundNode<<"\n";

  return roundNode;


}

// Funzione per incrementare il numero di alberi campionati.  | Serve una versione di dijkstra modificata |

void SamPG::encreaseForest(int Samples, NetworKit::Graph* graph, Labeling* index) {

      Dijkstra *dijkstra;
      for (int i = 0; i < Samples; i++) {
          Tree* new_tree = new Tree(this->random_roots->random(), graph->numberOfNodes());
          this->samplesForest.push_back(new_tree);
          dijkstra->runDijkstra(samplesForest.back(), graph, true, index);
          //this->totalNodes += new_tree->desc_vect[new_tree->getRoot()->node]+1;
		  int j=this->num_samples + i;
		  new_tree->computeDescendants(new_tree->getRoot(),counters, j, this->num_counters);
          //std::cout<<"Albero generato: \n";
          //this->samplesForest.back()->printTree(this->samplesForest.back()->getRoot());
      }
      this->num_samples += Samples;
  }


void SamPG::updateForest(int node){

    treeNode* maxNode;
    int i = 0;
    while(i<num_samples){
        if(samplesForest[i]->getRoot()->num_of_descendants != 0) {
            std::cout<<"numero di discendenti radice "<<i<<": "<<samplesForest[i]->getRoot()->num_of_descendants<<"\n";
            maxNode = this->samplesForest[i]->DFS(node);
            samplesForest[i]->deleteSubTree(maxNode, this->counters, i, this->num_counters);
            i++;
        }
        else{
            std::cout<<"VOGLIO CANCELLARE UN ALBERO!!\n";
            samplesForest[i] = samplesForest.back();
            this->num_samples -=1;
            samplesForest.resize(this->num_samples);
        }
    }
    std::cout<<"numero alberi nella foresta: "<<num_samples<<"\n";
}

int SamPG::getTotalNodes() {
    return this->totalNodes;
}



























