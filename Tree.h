#include <vector>
#include <iostream>



struct treeNode{
    int ID; // node's ID
    treeNode *father; // pointer to father
    std::vector<treeNode*> children; // vector of pointers to children
    int num_of_descendants = 0; // number of node's descendants
};


#ifndef TREE_H
#define TREE_H


class Tree{


private:

    treeNode treeRoot;

public:

    Tree();

    /**
     * @param root = root node ID
     */

    Tree(int root);

    /**
     * @return pointer tree's root
     */

    treeNode* getRoot();

    /**
     *
     * @param node = pointer to node from which start to decrease descendants
     * @param numdesc = number of descendants to be subtracted
     * @param counter = vector of counters
     *
     * This function decreases the number of descendants from all the ancestors starting from the node passed to
     * the root node. If a counter is passed the function updates it according to the decrease operation.
     *
     */


    void decreaseDescendants(treeNode* node, int numdesc, std::vector<int> &counter);

    /**
     *
     * @param node  = pointer to node to add in the tree
     * @param father = pointer to the father of the new child
     *
     * This function set as child of the father the node passed.
     */

    void addNode(treeNode* node, treeNode* father);

    /**
     *
     * @param node = node from which start to calculate the tree's descendants
     * @param counter = vector of counters
     *
     * This function calculates the descendants of all tree's nodes, starting from the node passed. If a vector is passed
     * the function updates it according to the values of descendans computed.
     */

    void computeDescendants(treeNode* node , std::vector<int> &counter);

    /**
     *
     * @param child = node to remove from the tree
     *
     * This function removes the node passed from the tree.
     */

    void removeChild(treeNode* child);

    /**
     *
     * @param father = new father's pointer to assign
     * @param child = node to change the father
     *
     * This function updates the father of the child node passed.
     */

    void updateFather(treeNode* father, treeNode* child);

    /**
     *
     * @param root = pointer to node from which start the print
     *
     * This function prints the tree starting from the root passed.
     */
    void printTree(treeNode* root);

    /**
     *
     * @param node = node ID to search
     * @return = a pointer to the node searched
     *
     * This function implements the Breadth First Search algorithm.
     */

    treeNode* BFS(int node);

    /**
     *
     * @param root = root node of the tree to remove
     * @param counter = vector of counters
     *
     * This function removes the tree rooted on the node passed. If a counter is passed it is updated according to the
     * changes of number of descendants.
     */

    void removeTree(treeNode* root , std::vector<int> &counter);

    /**
     *
     * @param root = root node of the subtrees to remove
     * @param counter = vector of counters
     *
     * This function removes all the subtrees rooted on children nodes of the passed node.
     */
    void deleteSubTree(treeNode* root , std::vector<int> &counter);
};

#endif //TREE_H