//
// Created by f4b3r on 28/03/20.
//

#include "Tree.h"
#include <iostream>




Tree::Tree(){}



//funziona
Tree::Tree(int source){
    treeNode r;
    r.node=source;
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

    std::cout << root->node << "\n";

    for(int i=0;i< root->children.size();i++){
        printTree(root->children[i]);
    }



}


/*int main(){

    Tree* t = new Tree(5);
    treeNode* root = t->getRoot();
    std::vector<treeNode> nodes;
    for(int i = 0; i<10; i++){
        treeNode n;
        n.node = i;
        nodes.push_back(n);
        t->addNode(&nodes[i], root);
    }
    for(int i = 0; i<root->children.size(); i++){
        std::cout<<root->children[i]->node<<"\n";
    }
    t->removeChild(&nodes[5]);
    for(int i = 0; i<root->children.size(); i++){
        std::cout<<root->children[i]->node<<"\n";
    }
}

*/
