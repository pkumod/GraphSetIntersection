#ifndef _BITPACK_MAXIMAL_CLIQUE_H
#define _BITPACK_MAXIMAL_CLIQUE_H

#include "util.hpp"
#include "set_operation.hpp"

class BPMaximalClique
{
public:
    int v_num, p_num;
    long long e_num;

    BPMaximalClique();
    ~BPMaximalClique();
    void build(const EdgeVector& _e_v);
    std::vector<int> degeneracy_order();

    int maximal_clique_degen();

    void save_answers(const char* file_path);

private:
    EdgeVector edge_vec;
    std::vector<UVertex> graph;
    std::vector<int> org_deg;
    int *pool_base = NULL;
    PackState *pool_state = NULL;

    int *sets_base = NULL;
    PackState *sets_state = NULL;

    int *pool_mc = NULL, pool_mc_idx = 0, mc_num = 0;
   
    int max_pool_sets_idx = 0, maximum_clique_size = 0;
    int intersect_call_time = 0, big_intersect_call_time = 0;

    void Tomita(std::vector<int>& R, UVertex P, UVertex X);
};

#endif