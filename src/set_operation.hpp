#ifndef _SET_OPERATION_H
#define _SET_OPERATION_H

#include "util.hpp"

int intersect(int *set_a, int size_a, int *set_b, int size_b, int *set_c);
int intersect_count(int *set_a, int size_a, int *set_b, int size_b);

int intersect_scalar2x(int *set_a, int size_a, int *set_b, int size_b, int *set_c);
int intersect_scalar2x_count(int* set_a, int size_a, int* set_b, int size_b);

int intersect_simd4x(int *set_a, int size_a, int *set_b, int size_b, int *set_c);
int intersect_simd4x_count(int* set_a, int size_a, int* set_b, int size_b);

int intersect_filter_simd4x(int *set_a, int size_a, int *set_b, int size_b, int *set_c);
int intersect_filter_simd4x_count(int* set_a, int size_a, int* set_b, int size_b);

int bp_intersect(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c);
int bp_intersect_count(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b);

int bp_intersect_scalar2x(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c);
int bp_intersect_scalar2x_count(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b);

int bp_intersect_simd4x(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c);
int bp_intersect_simd4x_count(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b);

int bp_intersect_filter_simd4x(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c);
int bp_intersect_filter_simd4x_count(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b);

int merge(int *set_a, int size_a, int *set_b, int size_b, int *set_c);
// merge one element in-place.
int bp_merge_one(int* bases_a, PackState* states_a, int size_a,
            int v_base, PackState v_bit);

int bp_subtract_visited(int* bases_a, PackState* states_a, int size_a,
            PackState* visited, int* bases_c, PackState* states_c);
int bp_subtract_unvisited(int* bases_a, PackState* states_a, int size_a,
            PackState* visited, int* bases_c, PackState* states_c);

int bp_subtract_visited_simd4x(int* bases_a, PackState* states_a, int size_a,
            PackState* visited, int* bases_c, PackState* states_c);
int bp_subtract_unvisited_simd4x(int* bases_a, PackState* states_a, int size_a,
            PackState* visited, int* bases_c, PackState* states_c);

int subtract(int *set_a, int size_a, int *set_b, int size_b, int *set_c);
int bp_subtract(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c);

extern unsigned long long inter_cnt, no_match_cnt, byte_check_cnt[4], cmp_cnt, multimatch_cnt, skew_cnt, low_select_cnt;
#endif

