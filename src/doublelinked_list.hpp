#ifndef _DL_LIST_H
#define _DL_LIST_H

#include "util.hpp"

struct DoubleLinkedNode
{
    int key;
    int prev, next;
    DoubleLinkedNode(int _k, int _p, int _n): key(_k), prev(_p), next(_n) {};
};

class DoubleLinkedList
{
public:
    DoubleLinkedList(int _capacity, int _range);
    ~DoubleLinkedList();
    void add(int key); // add a node to the list's tail, in O(1).
    void del(int key); // delete a node by key, in O(1).
    int get_head(); // return the key of the first node.
    int get_tail(); // return the key of the last node.
    int pop_head(); // delete the first node.
    int pop_tail(); // delete the last node.
    void print(); // scan the list from head to tail.
    void print_inverse(); // scan the list from tail to head.

private:
    DoubleLinkedNode *nodes = NULL;
    int *key2pos = NULL;
    int node_idx;
    int head, tail;
    int capacity, range;
};

#endif