//
// Created by f4b3r on 28/03/20.
//

#include "Tree.h"
#include <iostream>
#include <queue>
//#include </usr/lib/include/boost/lockfree/queque.hpp>




Tree::Tree(){}



//funziona
Tree::Tree(int source, int size){
    treeNode r;
    r.node = source;
    r.father = NULL;
    this->treeRoot = r;
}



//funziona
treeNode* Tree::getRoot(){
    return &treeRoot;
}

//funziona
void Tree::addNode(treeNode* node, treeNode* father){

    node->father = father;
    father->children.push_back(node);
}

void Tree::encreaseDescendants(treeNode* node){
    node->num_of_descendants += 1;
    if(node->father != NULL){ // se il nodo non è la radice, cioè non ha un genitore
        encreaseDescendants(node->father);
    }
}

void Tree::decreaseDescendants(treeNode* node ,int numdesc, std::vector<std::vector<int>> &counters, int j, int c){
    node->num_of_descendants -= numdesc;
	counters[j % c][node->node] -= numdesc;
    if(node->father != NULL){
        decreaseDescendants(node->father,numdesc, counters, j, c);
    }
}




//funziona
void Tree::computeDescendants(treeNode* node, std::vector<std::vector<int>> &counters, int j,int c){
    if(!node->children.empty()){
        for(int i = 0; i<node->children.size(); i++){
            computeDescendants(node->children[i], counters,j,c);
        }
    }
    int descendants = node->children.size();
    for(int i = 0; i < node->children.size(); i++){
        descendants += node->children[i]->num_of_descendants;
    }
	counters[j % c][node->node] += descendants;
    node->num_of_descendants = descendants;
}


//funziona
void Tree::removeChild(treeNode* child){

    if(child->father != NULL) {
        int size = child->father->children.size();
        for (int i = 0; i < size; i++) {
            if (child->node == child->father->children[i]->node){
                child->father->children.erase(child->father->children.begin()+i);
                child->father = NULL;
                break;
            }
        }
    }

}


//funziona
void Tree::updateFather(treeNode* father, treeNode* child){


    removeChild(child);
    child->father = father;
    father->children.push_back(child);
}




// funziona
void Tree::printTree(treeNode* root){

    std::cout << root->node << "-";
    for(int i=0; i<root->children.size(); i++){
        std::cout<<root->children[i]->node;
    }

    std::cout<<"\n";

    for(int i=0;i< root->children.size();i++){
        printTree(root->children[i]);
    }



}


treeNode* Tree::DFS(int node){

	bool trovato=false;
	treeNode* selectedNode;
	std::queue<treeNode*> queue;
	queue.push(this->getRoot());
	
	while(! queue.empty() & ! trovato){
		selectedNode=queue.front();
		queue.pop();

		if(selectedNode->node == node){
			trovato = true;;
		}
		else{
			for(int i=0;i<selectedNode->children.size();i++){

				queue.push(selectedNode->children[i]);
		
			}

		}
	
	}

	//out of cicle
	return selectedNode;

}

// Questa funzione elimina tutto l'albero radicato nel nodo root compresa la radice
void Tree::removeTree(treeNode* root, std::vector<std::vector<int>> &counters,int j,int c){

    if(!root->children.empty()){
        for(int i = 0; i<root->children.size(); i++){
            removeTree(root->children[i], counters, j, c);
        }
    }
	//std::cout<<"discendnti di : "<<root->node<<"da rimuovere"<<root->num_of_descendants<<"\n";
	counters[j % c][root->node] -= root->num_of_descendants;
    root->num_of_descendants = 0;
    root = NULL;


}

// Questa funzione rimuove tutti i sottoalberi radicati in un nodo
void Tree::deleteSubTree(treeNode* root, std::vector<std::vector<int>> &counters,int j,int c){
    if(root != NULL) {
        for (int i = 0; i < root->children.size(); i++) {
            removeTree(root->children[i], counters,j,c);

        }
        decreaseDescendants(root, root->num_of_descendants, counters,j,c);
        root->children.resize(0);
    }
}


