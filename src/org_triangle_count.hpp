#ifndef _ORG_TRIANGLE_COUNT_H
#define _ORG_TRIANGLE_COUNT_H

#include "util.hpp"
#include "set_operation.hpp"

class OrgTriangleCount
{
public:
    int v_num;
    long long e_num;   
   
    OrgTriangleCount();
    ~OrgTriangleCount();
    void build(const EdgeVector& _e_v);
    int triangle_count();
    
private:
    EdgeVector edge_vec;
    std::vector<UVertex> graph;
    int *pool_edges = NULL;

    // int intersection_count_8x(int* set_a, int size_a, int* set_b, int size_b);
    // int intersection_count_8x2(int* set_a, int size_a, int* set_b, int size_b);
};

// extern double intersect_time;
// extern unsigned long long intersect_cnt;
#endif