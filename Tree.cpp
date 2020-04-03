//
// Created by f4b3r on 28/03/20.
//

#include "Tree.h"
#include <iostream>
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
    encreaseDescendants(node);
}

void Tree::encreaseDescendants(treeNode* node){
    if(node->father != NULL){ // se il nodo non è la radice, cioè non ha un genitore
        node->father->num_of_descendants += 1;
        this->desc_vect[node->father->node] = node->father->num_of_descendants;
        encreaseDescendants(node->father);
    }
}

void Tree::decreaseDescendants(treeNode* node){
    if(node->father != NULL){
        node->father->num_of_descendants -= 1;
        this->desc_vect[node->father->node] = node->father->num_of_descendants;
        decreaseDescendants(node->father);
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
                decreaseDescendants(child);
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
    encreaseDescendants(child);
}




// funziona
void Tree::printTree(treeNode* root){

    std::cout << root->node << "\n";

    for(int i=0;i< root->children.size();i++){
        printTree(root->children[i]);
    }



}


treeNode* Tree::DFS(int node){

	bool trovato= false;
	treeNode* selectedNode;
	/*boost::lockfree::queque<treeNode*> queque;
	queque.push(this->getRoot());
	
	while(! queque.empty() & !trovato){
		queque.pop(selectedNode);

		if(selectedNode->node = node){
			break;
		}
		else{
			for(int i=0;i<selectedNode->children.size();i++){

				queque.push(selectedNode->children[i]);
		
			}

		}
	
	}*/

	//out of cicle
	return selectedNode;


}

// va richiamata sui figli non sul nodo altrimenti lo elimina
void Tree::removeSubTree(treeNode* root){

 	if(!root->children.empty()){
        for(int i = 0; i<root->children.size(); i++){
            removeSubTree(root->children[i]);
        }
        }

	std::cout << "cancello il nodo" << root->node <<"\n";
	delete root;
		
}

// va richiamata sul nodo radice del sottoalbero
void Tree::deleteSubTree(treeNode* root){

	for(int i=0;i<root->children.size();i++)
        {
	removeSubTree(root->children[i]);
		
	}
	root->children.resize(0);
}


/*int main(){

    Tree* t = new Tree(5,5);
    treeNode* root = t->getRoot();
    std::vector<treeNode*> nodes;
    for(int i = 0; i<10; i++){
        treeNode* n= new treeNode();
        n->node = i;
        nodes.push_back(n);
        t->addNode(nodes[i], root);
    }
		
	treeNode* n= new treeNode();
        n->node = 50;
        t->addNode(n, nodes[2]);
	std::cout<< "prima \n";
        t->printTree(t->getRoot());
	
	t->deleteSubTree(t->getRoot());
	
	
	std::cout<< "dopo \n";
    t->printTree(t->getRoot());



}*/


