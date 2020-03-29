//
// Created by f4b3r on 28/03/20.
//

#include <Tree.h>

Tree::Tree(){}

Tree::Tree(int root){
    treeNode* r;
    r.node = root;
}


treeNode* getRoot(){
    return this->r
}

void addDescendant(treeNode* father, treeNode child){
    father.children.push_back(child);
}

treeNode* addNode(int node, treeNode* father){

    treeNode new_node;
    new_node.node = node;
    new_node.father = father;
    this->addDescendant(father, new_node);
    return &new_node;
}

void computeDescendants(treeNode node){
    if(!node.children.empty()){
        for(int i = 0; i<node.children.size(); i++){
            computeDescendants(node.children[i]);
        }
    }
    int descendants = node.children.size()
    for(int i = 0; i < node.children.size(); i++){
        descendants += node.children.num_of_descendants
    }
    node.children.num_of_descendants = descendants
}

void removeChild(treeNode* father, treeNode* child){

    for(int i=0; i<father.children.size(); i++){
        if(father.children[i].node == child.node){
            father.children.erase(i);
        }
    }

}

void updateFather(treeNode* father, treeNode* child){

    this -> removeChild(father, child);
    child.father = father;
}