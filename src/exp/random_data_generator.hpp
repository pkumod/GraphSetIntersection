#ifndef _RAND_DATA_GEN_H
#define _RAND_DATA_GEN_H

#include "../util.hpp"

void gen_id_list(int len, double skew_ratio, double selectivity, double density,
    int*& set_a, int*& set_b, int& size_a, int& size_b);

#endif