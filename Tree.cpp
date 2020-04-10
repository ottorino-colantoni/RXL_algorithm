#include "Tree.h"
#include <iostream>
#include <queue>




Tree::Tree(){}


Tree::Tree(int source){
    treeNode r;
    r.ID = source;
    r.father = NULL;
    this->treeRoot = r;
}


treeNode* Tree::getRoot(){
    return &treeRoot;
}


void Tree::addNode(treeNode* node, treeNode* father){

    node->father = father;
    father->children.push_back(node);
}


void Tree::decreaseDescendants(treeNode* node ,int numdesc, std::vector<int> &counter){
    node->num_of_descendants -= numdesc;
    counter[node->ID] -= numdesc;
    if(node->father != NULL){
        decreaseDescendants(node->father,numdesc, counter);
    }
}


void Tree::computeDescendants(treeNode* node, std::vector<int> &counter){
    if(!node->children.empty()){
        for(int i = 0; i<node->children.size(); i++){
            computeDescendants(node->children[i], counter);
        }
    }
    int descendants = node->children.size();
    for(int i = 0; i < node->children.size(); i++){
        descendants += node->children[i]->num_of_descendants;
    }
    counter[node->ID] += descendants;
    node->num_of_descendants = descendants;
}


void Tree::removeChild(treeNode* child){

    if(child->father != NULL) {
        int size = child->father->children.size();
        for (int i = 0; i < size; i++) {
            if (child->ID == child->father->children[i]->ID){
                child->father->children.erase(child->father->children.begin()+i);
                child->father = NULL;
                break;
            }
        }
    }

}


void Tree::updateFather(treeNode* father, treeNode* child){


    removeChild(child);
    child->father = father;
    father->children.push_back(child);
}


void Tree::printTree(treeNode* root){

    std::cout << root->ID << "-";
    for(int i=0; i<root->children.size(); i++){
        std::cout<<root->children[i]->ID;
    }

    std::cout<<"\n";

    for(int i=0;i< root->children.size();i++){
        printTree(root->children[i]);
    }
}


treeNode* Tree::BFS(int node){

	bool trovato=false;
	treeNode* selectedNode;
	std::queue<treeNode*> queue;
	queue.push(this->getRoot());
	
	while(! queue.empty() & ! trovato){
		selectedNode=queue.front();
		queue.pop();

		if(selectedNode->ID == node){
			trovato = true;;
		}
		else{
			for(int i=0;i<selectedNode->children.size();i++){

				queue.push(selectedNode->children[i]);
		
			}

		}
	
	}
	return selectedNode;

}


void Tree::removeTree(treeNode* root, std::vector<int> &counter){

    if(!root->children.empty()){
        for(int i = 0; i<root->children.size(); i++){
            removeTree(root->children[i], counter);
        }
    }
    counter[root->ID] -= root->num_of_descendants;
    root->num_of_descendants = 0;
    root = NULL;


}


void Tree::deleteSubTree(treeNode* root, std::vector<int> &counter){
    if(root != NULL) {
        for (int i = 0; i < root->children.size(); i++) {
            removeTree(root->children[i], counter);
        }
        decreaseDescendants(root, root->num_of_descendants, counter);
        root->children.resize(0);
    }
}


