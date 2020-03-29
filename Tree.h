//
// Created by f4b3r on 28/03/20.
//

#ifndef EXAMPLE_TREE_H
#define EXAMPLE_TREE_H

#endif //EXAMPLE_TREE_H

struct treeNode{
    int node;
    treeNode *father;
    Vector<treeNode*> children;
    int num_of_descendants;
};


class Tree{

private:

    treeNode root;


public:

    Tree();

    Tree(int root);

    void addNode(int node, treeNode* father);

    void computeDescendants();

    void removeSubtree(int node);

    void removeChild(treeNode* father, treeNode* child);

    void updateFather(treeNode* father, treeNode* child);
};
