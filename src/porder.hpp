#ifndef _PORDER_H
#define _PORDER_H

#include "util.hpp"
#include "metis.h"

struct NbrNode
{
    int vid;
    bool type; // 0 for out-neighbors, 1 for in-neighbors.
    NbrNode(int _v, int _t): vid(_v), type(_t) {};
};

struct NodeWithGain
{
    int u;
    double gain;
    NodeWithGain(): u(0), gain(0.0) {};
    NodeWithGain(int _u, double _g): u(_u), gain(_u) {};
    bool operator < (const NodeWithGain& p) const {return gain > p.gain;};
};

class POrder
{
public:
    int v_num; 
    long long e_num;
    int p_num; // pack num: ceil(v_num/PACK_WIDTH).
    std::vector<int> org2newid; // org_id: i --> new_id: org2newid[i].
    void load_org_graph(EdgeVector _e_v);
    EdgeVector hybrid_bfsdeg();
    EdgeVector greedy_naive();
    EdgeVector greedy_mheap();
    EdgeVector mloggapa_order(); // MLOGGAPA order using Recursive Graph Bisection algorithm
    EdgeVector metis_order(); // METIS partitioning algorithm
    EdgeVector slashburn_order(); // SlashBurn algorithm
    EdgeVector bfsr_order(); // Recursive BFS order algorithm
    EdgeVector dfs_order(); // DFS order algorithm
    void set_alpha_by_deg();
    void set_alpha(double *_a_out, double *_a_in);
    
    int leaf_node_count();
    std::vector<int> select_bignode(double deg_ratio);
    double comp_ratio();

    POrder();
    ~POrder();

private:
    EdgeVector edge_vec;
    std::vector<DVertex> graph;
    std::vector<int> outedge, inedge;
    std::vector<NbrNode> nbr;

    int *new_id;
    double *alpha_out, *alpha_in;

    void build(); // build graph's adjacent lists by edge_vec/new_id.
    void sort_nbr();
    void deg_order();
    void rcm_order();
    void bfs_order();
    void deg_desc_order();

    void dfs(int u, int* nid, int& cur_idx);
    // for Recursive Graph Bisection Order
    void graph_bisection(NodeWithGain *nodes, int *set_label, double *gain,
        int *q_l_nodes, int *q_r_nodes, int tot_num, int level, int& cur_label);
    void graph_bisection2(NodeWithGain *nodes, int *set_label, double *gain,
        int tot_num, int& cur_label);
    void bfsr_bisection(int *nodes, int tot_num, int *que, int *visited, int& vis_label);
};

#endif