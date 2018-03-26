#ifndef _ORG_MAXIMAL_CLIQUE_H
#define _ORG_MAXIMAL_CLIQUE_H

#include "util.hpp"
#include "set_operation.hpp"

class OrgMaximalClique
{
public:
    int v_num;
    long long e_num;

    OrgMaximalClique();
    ~OrgMaximalClique();

    void build(const EdgeVector& _e_v);
    std::vector<int> degeneracy_order();
    int maximal_clique_bk();
    int maximal_clique_pivot();
    int maximal_clique_degen();

    void save_answers(const char* file_path);

    // double intersect_time = 0.0;
    // unsigned long long intersect_cnt = 0;
    // struct timeval time_start;
    // struct timeval time_end;
private:
    EdgeVector edge_vec;
    std::vector<UVertex> graph;
    int *pool_edges = NULL;
    int *pool_sets = NULL;
    int *pool_mc = NULL, pool_mc_idx = 0, mc_num = 0;

    int *temp_set = NULL;
    int max_pool_sets_idx = 0, maximum_clique_size = 0;
    int intersect_call_time = 0, big_intersect_call_time = 0;

    void BronKerbosch(std::vector<int>& R, UVertex P, UVertex X, int mem_idx);
    void Tomita(std::vector<int>& R, UVertex P, UVertex X);
};


#endif