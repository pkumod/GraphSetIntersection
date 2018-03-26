#ifndef _RBM_MAXIMAL_CLIQUE_H
#define _RBM_MAXIMAL_CLIQUE_H

#include "util.hpp"
#include "./roaring/roaring.hh"

class RoaringMaximalClique
{
public:
    int v_num;
    long long e_num;

    RoaringMaximalClique();

    void build(const EdgeVector& _e_v);
    std::vector<int> degeneracy_order();
    int maximal_clique_degen();

    void save_answers(const char* file_path);
private:
    EdgeVector edge_vec;
    std::vector<Roaring> graph;

    std::vector<Roaring> maximal_cliques;
    int mc_num = 0, maximum_clique_size = 0;

    void Tomita(Roaring& R, Roaring& P, Roaring& X);
};

#endif