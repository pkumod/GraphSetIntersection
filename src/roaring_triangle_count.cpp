#include "roaring_triangle_count.hpp"

RoaringTriangleCount::RoaringTriangleCount()
{
    v_num = 0;
    e_num = 0;
}

void RoaringTriangleCount::build(const EdgeVector& _e_v)
{
    EdgeVector rev_edge_vec;
    rev_edge_vec.reserve(_e_v.size() / 2);
    edge_vec.reserve(_e_v.size() / 2);
    for (const auto& e : _e_v) {
        if (e.first < e.second) {
            v_num = std::max(v_num, e.second);
            edge_vec.push_back(e);
        } else if (e.first > e.second) {
            rev_edge_vec.push_back(e);
        }
    }
    v_num++;
    
    
    std::sort(edge_vec.begin(), edge_vec.end(), edge_idpair_cmp);
    edge_vec.erase(std::unique(edge_vec.begin(), edge_vec.end()), edge_vec.end());
    std::sort(rev_edge_vec.begin(), rev_edge_vec.end(), edge_idpair_cmp);
    rev_edge_vec.erase(std::unique(rev_edge_vec.begin(), rev_edge_vec.end()), rev_edge_vec.end());
    
    e_num = (long long)edge_vec.size();
    graph.resize(v_num);

    for (auto& e : rev_edge_vec) graph[e.first].add(e.second);

    printf("v_num=%d e_num=%lld\n", v_num, e_num);
}

int RoaringTriangleCount::triangle_count()
{
    int res = 0;
    for (const auto& e : edge_vec) {
        const Roaring& Nu = graph[e.first];
        const Roaring& Nv = graph[e.second];
        res += Nu.and_cardinality(Nv);
    }   

    return res;
}
