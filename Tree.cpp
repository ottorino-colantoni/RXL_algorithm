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
    this->desc_vect.resize(size);
}



//funziona
treeNode* Tree::getRoot(){
    return &treeRoot;
}

std::vector<int> Tree::getDescVect(){
    return this->desc_vect;
}



//funziona
void Tree::addNode(treeNode* node, treeNode* father){

    node->father = father;
    father->children.push_back(node);
    encreaseDescendants(node->father);
}

void Tree::encreaseDescendants(treeNode* node){
    node->num_of_descendants += 1;
    this->desc_vect[node->node] = node->num_of_descendants;
    if(node->father != NULL){ // se il nodo non è la radice, cioè non ha un genitore
        encreaseDescendants(node->father);
    }
}

void Tree::decreaseDescendants(treeNode* node ,int numdesc){
    node->num_of_descendants -= numdesc;
    this->desc_vect[node->node] = node->num_of_descendants;
    if(node->father != NULL){
        decreaseDescendants(node->father,numdesc);
    }
}




//funziona
void Tree::computeDescendants(treeNode* node){
    if(!node->children.empty()){
        for(int i = 0; i<node->children.size(); i++){
            computeDescendants(node->children[i]);
        }
    }
    int descendants = node->children.size();
    for(int i = 0; i < node->children.size(); i++){
        descendants += node->children[i]->num_of_descendants;
    }
    node->num_of_descendants = descendants;
}


//funziona
void Tree::removeChild(treeNode* child){

    if(child->father != NULL) {
        int size = child->father->children.size();
        for (int i = 0; i < size; i++) {
            if (child->node == child->father->children[i]->node){
                child->father->children.erase(child->father->children.begin()+i);
                decreaseDescendants(child->father, 1);
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
    encreaseDescendants(child->father);
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
void Tree::removeTree(treeNode* root){

    if(!root->children.empty()){
        for(int i = 0; i<root->children.size(); i++){
            removeTree(root->children[i]);
        }
    }

    this->desc_vect[root->node] = 0;
    delete root;


}

// Questa funzione rimuove tutti i sottoalberi radicati in un nodo
void Tree::deleteSubTree(treeNode* root){

    for(int i=0;i<root->children.size();i++)
    {
        removeTree(root->children[i]);

    }
    decreaseDescendants(root, root->num_of_descendants);
    root->num_of_descendants = 0;
    root->children.resize(0);
}


