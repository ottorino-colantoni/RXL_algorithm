//
// Created by f4b3r on 28/03/20.
//


#include <vector>

struct treeNode{
    int node;
    treeNode *father;
    std::vector<treeNode*> children;
    int num_of_descendants;
};


#ifndef TREE_H
#define TREE_H


class Tree{

private:

    treeNode treeRoot;


public:

    Tree();

    Tree(int root);

    treeNode* getRoot();

    void addNode(treeNode* node, treeNode* father);

    //void addDescendant(treeNode* father, treeNode* child);

    void computeDescendants(treeNode* node);

    //void removeSubtree(int node);

    void removeChild(treeNode* child);

    void updateFather(treeNode* father, treeNode* child);

    void printTree(treeNode* root);
};

#endif //EXAMPLE_TREE_H