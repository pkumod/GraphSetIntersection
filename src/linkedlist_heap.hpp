#ifndef _L_HEAP_H
#define _L_HEAP_H

#include "util.hpp"

struct LHNode
{
    int state, val;
    int prev, next;

    LHNode() {};
    LHNode(int _s, int _v, int _p, int _n): state(_s), val(_v), prev(_p), next(_n) {};
};

struct LHHeader
{
    int first, last;
    int label;
    
    LHHeader() {};
    LHHeader(int _f, int _l, int _lb): first(_f), last(_l), label(_lb) {};
};

struct LHUpdateInfo
{
    int key;
    int old_val;

    LHUpdateInfo() {};
    LHUpdateInfo(int _k, int _v): key(_k), old_val(_v) {};
};

class LinkedListHeap
{
public:
    LinkedListHeap(int _capacity);
    ~LinkedListHeap();

    void inc(int key);
    int top();
    int pop();
    void del(int key);
    void reset();
    void adjust();
    bool in_heap(int key);
    bool is_top_zero();
    int get_top_val();
    int get_size();

    void print();
    bool check();
    bool check_size();
private:
    const int RESET_LABEL_MASK = 0x7fffffff;
    const int UPDATE_LABEL_MASK = 0x80000000;

    int capacity = 0, size = 0, reset_label = 0;
    int head = -1, tail = -1;

    LHNode *linkedlist;
    LHHeader *headers;

    LHUpdateInfo *up_list;
    int up_list_cnt = 0;

    void up(int key, int old_val);
};

#endif