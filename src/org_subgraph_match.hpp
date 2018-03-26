#ifndef _ORG_SUBGRAPH_MATCH_H
#define _ORG_SUBGRAPH_MATCH_H

#include "util.hpp"
#include "set_operation.hpp"

struct LabelSubgraph
{
    int v_num, e_num;
    EdgeVector edge_vec;
    std::vector<int> vertex2label;

    LabelSubgraph(int _v, int _e): v_num(_v), e_num(_e) {};
    void print() const;
};


class OrgSubGraphMatch
{
public:
    int v_num, l_num;
    long long e_num;

    OrgSubGraphMatch();
    ~OrgSubGraphMatch();

    void build(const EdgeVector& _e_v, const std::vector<int> _v_l);
    std::vector<std::vector<int>> subgraph_matching(const LabelSubgraph& q);

private:
    EdgeVector edge_vec;
    std::vector<UVertex> graph, labels;
    std::vector<int> vertex2label;    
    int *pool_edges = NULL, *pool_vertices = NULL;
    int *temp_join_res[2] = {NULL, NULL}, temp_join_size[2];

    void gen_join_order(int u, std::vector<std::vector<int>>& q_edge_list,
            std::vector<int>& visited, std::vector<int>& join_order);

};



#endif