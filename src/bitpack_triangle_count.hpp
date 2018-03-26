#ifndef _BITPACK_TRIANGLE_COUNT_H
#define _BITPACK_TRIANGLE_COUNT_H

#include "util.hpp"
#include "set_operation.hpp"

class BPTriangleCount
{
public:
    int v_num;
    long long e_num;   

    BPTriangleCount();
    ~BPTriangleCount();
    void build(const EdgeVector& _e_v);
    int triangle_count();
    int triangle_count_mt(int thread_num); // multi-threading using Intel TBB.
    // void calc_triangle_on_edge(long long e_idx_l, long long e_idx_r, int& res);

// private:
    EdgeVector edge_vec;
    std::vector<UVertex> graph;
    int *pool_base = NULL;
    PackState *pool_state = NULL;
    
    // int bp_intersection_count_2x(int* bases_a, PackState* states_a, int size_a,
    //     int* bases_b, PackState* states_b, int size_b);
    // int bp_intersection_count_8x(int* bases_a, PackState* states_a, int size_a,
    //     int* bases_b, PackState* states_b, int size_b); 
    // int bp_intersection_count_8x2(int* bases_a, PackState* states_a, int size_a,
    //     int* bases_b, PackState* states_b, int size_b); 
};

#endif