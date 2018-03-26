#ifndef _M_HEAP_H
#define _M_HEAP_H

#include "util.hpp"

struct MHNode
{
    int key;
    int label;
    double val;

    MHNode() {};
    MHNode(int _k, int _l, double _v): key(_k), label(_l), val(_v) {};
};

struct UpdateInfo
{
    int pos; // pos == 0 means the node is not in heap.
    int state;
    double utd_val; // up-to-date values of each key.

    UpdateInfo() {};
    UpdateInfo(int _p, int _s, double _v): pos(_p), state(_s), utd_val(_v) {};
};

class ModifiedHeap
{
public:
    ModifiedHeap(int max_size);
    ~ModifiedHeap();

    void inc(int key, double delta = 1.0);
    MHNode top();
    int pop();
    void del(int key);
    void reset(); // reset nodes' values in heap.
    void adjust();
    int get_size();
    bool in_heap(int key);
    
private:
    int size = 0, reset_label = 0;
    MHNode *heap;
    UpdateInfo *up_info;
    int *up_list, up_list_cnt;

    bool is_zero(int pos);
    void up(int pos);
    void down(int pos);
    void swap_mhnode(int pos_x, int pos_y);

};

#endif