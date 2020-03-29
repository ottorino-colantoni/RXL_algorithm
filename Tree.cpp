//
// Created by f4b3r on 28/03/20.
//

#include <Tree.h>

Tree::Tree(){}

Tree::Tree(int root){
    treeNode s;
    s.node = root;
}


void addNode(int node, treeNode* father){

    treeNode new_node;
    new_node.node = node;
    new_node.father = father;
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