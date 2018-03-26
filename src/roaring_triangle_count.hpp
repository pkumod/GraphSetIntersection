#ifndef _RBM_TRIANGLE_COUNT_H
#define _RBM_TRIANGLE_COUNT_H

#include "util.hpp"
#include "./roaring/roaring.hh"

class RoaringTriangleCount
{
public:
    int v_num;
    long long e_num;   
   
    RoaringTriangleCount();
    void build(const EdgeVector& _e_v);
    int triangle_count();
    
private:
    EdgeVector edge_vec;
    std::vector<Roaring> graph;
};


#endif