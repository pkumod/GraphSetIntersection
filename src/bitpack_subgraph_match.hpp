#ifndef _BITPACK_SUBGRAPH_MATCH_H
#define _BITPACK_SUBGRAPH_MATCH_H

#include "util.hpp"
#include "set_operation.hpp"
#include "org_subgraph_match.hpp"

class BPSubGraphMatch
{
public:
    int v_num, l_num;
    long long e_num;

    BPSubGraphMatch();
    ~BPSubGraphMatch();

    void build(const EdgeVector& _e_v, const std::vector<int> _v_l);
    std::vector<std::vector<int>> subgraph_matching(const LabelSubgraph& q);

private:
    EdgeVector edge_vec;
    std::vector<UVertex> graph, labels;
    std::vector<int> vertex2label;
    int *pool_base = NULL, *label_base = NULL;
    PackState *pool_state = NULL, *label_state = NULL;
};

#endif