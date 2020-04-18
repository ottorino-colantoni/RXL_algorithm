//
// Created by f4b3r on 31/03/20.
//

#include "SamPG.h"
#include "mytimer.h"
#include <omp.h>

SamPG::SamPG(){};

SamPG::SamPG(int k, int c, NetworKit::Graph* g){
    this->graph = g;
    this->num_samples = k;
	this->num_counters= c;
	this->counters.resize(c);
	this->totalNodes = 0;
	this->EndGeneration = false;
	this->last_node_checked = 0;
    this->already_drawn.resize(this->graph->numberOfNodes(), false);
    this->random_roots= new BestRandom(0, this->graph->numberOfNodes());
    for(int i=0; i<this->num_counters;i++){
        this->counters[i].resize(graph->numberOfNodes(),0);
    }

}

void SamPG::createForest(){
    Dijkstra* dijkstra;
    for(int i = 0; i<this->num_samples; i++){
        Tree* new_tree = new Tree(this->random_roots->random());
        this->samplesForest.push_back(new_tree);
        dijkstra->runDijkstra(this->samplesForest.back(), this->graph);
		new_tree->computeDescendants(new_tree->getRoot(), this->counters[i % this->num_counters]);
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

   if(this->EndGeneration == false) {
       for (int i = 0; i < counters[0].size(); i++) {

           for (int j = 0; j < this->num_counters; j++) {

               temporarymax += this->counters[j][i];
               if (this->counters[j][i] > useless_value1) { useless_value1 = counters[j][i]; }

           }

           if (temporarymax > max_for_zero) {
               max_for_zero = temporarymax;
               roundNode_for_zero = i;
           }
           if(num_counters>1) {
               temporarymax -= useless_value1;
               temporarymax = (temporarymax / (this->num_counters - 1));
           }
           if (temporarymax > max) {
               max = temporarymax;
               roundNode = i;
           }

           temporarymax = 0;
           useless_value1 = 0;
       }
   }

   if(this->already_drawn[roundNode] == false){
       this->already_drawn[roundNode] = true;
   }
   else{
       if(roundNode == 0 && roundNode_for_zero != 0){
           roundNode = roundNode_for_zero;
           this->already_drawn[roundNode] = true;
       }
       else{
           this->EndGeneration = true;
           for(int i = this->last_node_checked; i< already_drawn.size(); i++){
               if(already_drawn[i] == false){
                   this->last_node_checked = i;
                   roundNode = i;
                   this->already_drawn[i] = true;
                   break;
               }

           }
       }
   }
    return  roundNode;
}

// Funzione per incrementare il numero di alberi campionati.  | Serve una versione di dijkstra modificata |

void SamPG::encreaseForest(int Samples, Labeling* index) {

      Dijkstra *dijkstra;
    for (int i = 0; i < Samples; i++) {
          Tree* new_tree = new Tree(this->random_roots->random());
          this->samplesForest.push_back(new_tree);
          dijkstra->runDijkstra(this->samplesForest.back(), this->graph, true, index);
		  int j = this->num_samples + i;
		  new_tree->computeDescendants(new_tree->getRoot(),this->counters[j % this->num_counters]);
		  this->totalNodes += new_tree->getRoot()->num_of_descendants + 1;
      }
      this->num_samples += Samples;
  }


void SamPG::updateForest(int node){

    treeNode* maxNode;
    for (int i = 0; i < samplesForest.size(); i++) {
        maxNode = this->samplesForest[i]->BFS(node);
        samplesForest[i]->deleteSubTree(maxNode, this->counters[i % this->num_counters]);
    }
}

int SamPG::getTotalNodes() {
    return this->totalNodes;
}

bool SamPG::isEnded(){
    return this->EndGeneration;
}

int SamPG::getNumSamples(){
    return this->num_samples;
};



























