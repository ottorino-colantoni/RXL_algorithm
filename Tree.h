//
// Created by f4b3r on 28/03/20.
//


#include <vector>
#include <iostream>

struct treeNode{
    int node;
    treeNode *father;
    std::vector<treeNode*> children;
    int num_of_descendants = 0;
};


#ifndef TREE_H
#define TREE_H


class Tree{

private:

    treeNode treeRoot;

public:


    std::vector<treeNode*> direct_acc;

    Tree();

    Tree(int root, int size);

    treeNode* getRoot();

    void encreaseDescendants(treeNode* node);

    void decreaseDescendants(treeNode* node, int num_of_descendants, std::vector<std::vector<int>> &counters,int k,int c);

    void addNode(treeNode* node, treeNode* father);

    //void addDescendant(treeNode* father, treeNode* child);

    void computeDescendants(treeNode* node , std::vector<std::vector<int>> &counters, int k, int c);

    //void removeSubtree(int node);

    void removeChild(treeNode* child);

    void updateFather(treeNode* father, treeNode* child);

    void printTree(treeNode* root);

    treeNode* DFS(int node);

    void removeTree(treeNode* root , std::vector<std::vector<int>> &counters, int k, int c);

    void deleteSubTree(treeNode* root , std::vector<std::vector<int>> &counters, int k, int c);
};

#endif //EXAMPLE_TREE_H