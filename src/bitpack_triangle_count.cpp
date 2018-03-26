#include "bitpack_triangle_count.hpp"


BPTriangleCount::BPTriangleCount()
{
    v_num = 0;
    e_num = 0;
    align_malloc((void**)&pool_base, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
    align_malloc((void**)&pool_state, 32, sizeof(PackState) * PACK_NODE_POOL_SIZE);
}

BPTriangleCount::~BPTriangleCount()
{
    free(pool_state);
    free(pool_base);
}

void BPTriangleCount::build(const EdgeVector& _e_v)
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
     
    int cur_packnode_idx = -1;
    int prev_u = -1;
    for (auto& e : rev_edge_vec) {
        PackBase v_base = (e.second >> PACK_SHIFT);
        PackState v_bit = ((PackState)1 << (e.second & PACK_MASK));
        if (e.first != prev_u) {
            prev_u = e.first;            
            graph[e.first].start = ++cur_packnode_idx;
            graph[e.first].deg++;
            pool_base[cur_packnode_idx] = v_base;
            pool_state[cur_packnode_idx] = v_bit;
        } else {            
            if (pool_base[cur_packnode_idx] == v_base) {
                pool_state[cur_packnode_idx] |= v_bit;
            } else {
                graph[e.first].deg++;
                pool_base[++cur_packnode_idx] = v_base;
                pool_state[cur_packnode_idx] = v_bit;
            }
        }
    }
    cur_packnode_idx++;

    double comp_ratio = (double)cur_packnode_idx / e_num;
    printf("comp_ratio=%d/%lld, %.4f\n", cur_packnode_idx, e_num, comp_ratio);
}

int BPTriangleCount::triangle_count()
{
    int res = 0;
    // int total_len = 0, long_len = 0;
    for (const auto& e : edge_vec) {
        const UVertex& u = graph[e.first];
        const UVertex& v = graph[e.second];

        // gettimeofday(&time_start, NULL);
#if SIMD_STATE == 2
        res += bp_intersect_scalar2x_count(pool_base + u.start, pool_state + u.start, u.deg,
                pool_base + v.start, pool_state + v.start, v.deg);        
#elif SIMD_STATE == 4
#if SIMD_MODE == 0
        res += bp_intersect_simd4x_count(pool_base + u.start, pool_state + u.start, u.deg,
                pool_base + v.start, pool_state + v.start, v.deg);
#else
        res += bp_intersect_filter_simd4x_count(pool_base + u.start, pool_state + u.start, u.deg,
                pool_base + v.start, pool_state + v.start, v.deg);
#endif
#else
        res += bp_intersect_count(pool_base + u.start, pool_state + u.start, u.deg,
                pool_base + v.start, pool_state + v.start, v.deg);        
#endif

        // gettimeofday(&time_end, NULL);
        // total_len += u.deg + v.deg;
        // if (u.deg < 16 && v.deg < 16) long_len += u.deg + v.deg;
    }
    
    // printf("long_ratio: %.3f%%(%d/%d)\n", (double)long_len/total_len * 100,
    //         long_len, total_len);
    return res;
}

// struct BPTCQuery
// {
//     int *base_u, *base_v;
//     PackState* state_u, *state_v;
//     int size_u, size_v;
//     BPTCQuery() {};
//     BPTCQuery(int *_b_u, int *_b_v, PackState *_s_u, PackState *_s_v,
//         int _l_u, int _l_v): base_u(_b_u), base_v(_b_v),
//         state_u(_s_u), state_v(_s_v), size_u(_l_u), size_v(_l_v) {};
// };

// class ConcurrentTCProcessor
// {
//     BPTriangleCount *bptc;
//     int *res;
    
// public:
//     void operator() (const blocked_range<size_t>& r) const
//     {
//         for (size_t i=r.begin();i!=r.end();++i)
//         {
//             long long e_idx_l = bptc->e_num / (r.end()) * i;
//             long long e_idx_r = bptc->e_num / (r.end()) * (i + 1) - 1;
//             if (i + 1 == r.end()) e_idx_r = bptc->e_num - 1;
//             // printf("i=%d pid=%d\n", i, pid);

//             bptc->calc_triangle_on_edge(e_idx_l, e_idx_r, res[i]);            

//         }
//     }
//     ConcurrentTCProcessor(BPTriangleCount *_bptc, int *_res)
//     {
//         bptc = _bptc;
//         res = _res;
//     }

// };

int con_thread_num = 1;
int *con_res = NULL;
BPTriangleCount *bptc;

void* con_calc_triangle(void *con)
{
    long id;
    id = (unsigned long long)(con);
    long long e_idx_l = bptc->edge_vec.size() / con_thread_num * id;
    long long e_idx_r = bptc->edge_vec.size() / con_thread_num * (id + 1) - 1;
    if (id + 1 == con_thread_num) e_idx_r = bptc->edge_vec.size() - 1;

    for (long long i = e_idx_l; i < e_idx_r; ++i) {
        const auto& e = bptc->edge_vec[i];
        const UVertex& u = bptc->graph[e.first];
        const UVertex& v = bptc->graph[e.second];

#if SIMD_STATE == 2
        con_res[id] += bp_intersect_scalar2x_count(bptc->pool_base + u.start, bptc->pool_state + u.start, u.deg,
                bptc->pool_base + v.start, bptc->pool_state + v.start, v.deg);        
#elif SIMD_STATE == 4
#if SIMD_MODE == 0
        con_res[id] += bp_intersect_simd4x_count(bptc->pool_base + u.start, bptc->pool_state + u.start, u.deg,
                bptc->pool_base + v.start, bptc->pool_state + v.start, v.deg);
#else
        con_res[id] += bp_intersect_filter_simd4x_count(bptc->pool_base + u.start, bptc->pool_state + u.start, u.deg,
                bptc->pool_base + v.start, bptc->pool_state + v.start, v.deg);
#endif
#else
        con_res[id] += bp_intersect_count(bptc->pool_base + u.start, bptc->pool_state + u.start, u.deg,
                bptc->pool_base + v.start, bptc->pool_state + v.start, v.deg);        
#endif
    }   
    pthread_exit(NULL);
}


int BPTriangleCount::triangle_count_mt(int thread_num)
{
    con_thread_num = thread_num; 
    align_malloc((void**)&con_res, 32, sizeof(int) * thread_num);
    printf("thread_num=%d\n", thread_num);
    memset(con_res, 0, sizeof(int) * thread_num);
    bptc = this;

    pthread_t *pt = (pthread_t *)malloc(con_thread_num * sizeof(pthread_t));
    struct timeval time_start;
    struct timeval time_end;
    gettimeofday(&time_start, NULL);
    for (long a = 0; a < con_thread_num; a++)
        pthread_create(&pt[a], NULL, con_calc_triangle, (void*)a);
    for (long a = 0; a < con_thread_num; a++)
        pthread_join(pt[a], NULL);
    gettimeofday(&time_end, NULL);
    double tc_mt_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("tc_mt_time=%.3fms\n", tc_mt_time);
    free(pt);
    // parallel_for(blocked_range<size_t>(0, thread_num),
    //         ConcurrentTCProcessor(this, res));

    int tot_res = 0;
    for (int i = 0; i < thread_num; ++i) tot_res += con_res[i];
    return tot_res;
}



// void BPTriangleCount::calc_triangle_on_edge(long long e_idx_l, long long e_idx_r, int& res)
// {
//     for (long long i = e_idx_l; i < e_idx_r; ++i) {
//         const auto& e = edge_vec[i];
//         const UVertex& u = graph[e.first];
//         const UVertex& v = graph[e.second];

// #if SIMD_STATE == 2
//         res += bp_intersect_scalar2x_count(pool_base + u.start, pool_state + u.start, u.deg,
//                 pool_base + v.start, pool_state + v.start, v.deg);        
// #elif SIMD_STATE == 4
// #if SIMD_MODE == 0
//         res += bp_intersect_simd4x_count(pool_base + u.start, pool_state + u.start, u.deg,
//                 pool_base + v.start, pool_state + v.start, v.deg);
// #else
//         res += bp_intersect_filter_simd4x_count(pool_base + u.start, pool_state + u.start, u.deg,
//                 pool_base + v.start, pool_state + v.start, v.deg);
// #endif
// #else
//         res += bp_intersect_count(pool_base + u.start, pool_state + u.start, u.deg,
//                 pool_base + v.start, pool_state + v.start, v.deg);        
// #endif
//     }
// }

// inline int BPTriangleCount::bp_intersection_count_2x(int* bases_a, PackState* states_a, int size_a,
//         int* bases_b, PackState* states_b, int size_b)
// {
//     int i = 0, j = 0, res = 0;
//     int ds_a = size_a - (size_a & 1);
//     int ds_b = size_b - (size_b & 1);

//     while (i < ds_a && j < ds_b) {
//         __m64 base_a = *(__m64*)(bases_a + i);
//         __m64 base_b = *(__m64*)(bases_b + j);
//         __m64 state_a = *(__m64*)(states_a + i);
//         __m64 state_b = *(__m64*)(states_b + j);

//         int a_max = bases_a[i + 1];
//         int b_max = bases_b[j + 1];
//         if (a_max == b_max) {
//             i += 2;
//             j += 2;
//         } else if (a_max < b_max) {
//             i += 2;
//         } else {
//             j += 2;
//         }

//         __m64 and_all = _mm_setzero_si64();

//         // shift0:        
//         __m64 cmp_mask0 = _mm_cmpeq_pi32(base_a, base_b);
//         __m64 and_0 = _mm_and_si64(state_a, state_b);
//         and_all = _mm_or_si64(and_all, _mm_and_si64(and_0, cmp_mask0));

//         // shift1:
//         __m64 base_sf1 = _mm_shuffle_pi16(base_b, 0x004e);
//         __m64 state_sf1 = _mm_shuffle_pi16(state_b, 0x004e);
//         __m64 cmp_mask1 = _mm_cmpeq_pi32(base_a, base_sf1);
//         __m64 and_1 = _mm_and_si64(state_a, state_sf1);
//         and_all = _mm_or_si64(and_all, _mm_and_si64(and_1, cmp_mask1));

//         res += __builtin_popcountl((long long)and_all);
//     }

//     while (i < size_a && j < size_b) {
//         if (bases_a[i] == bases_b[j]) {
//             res += __builtin_popcount(states_a[i] & states_b[j]);
//             i++; j++;
//         } else if (bases_a[i] < bases_b[j]) {
//             i++;
//         } else {
//             j++;
//         }
//     }

//     return res;
// }

// inline int BPTriangleCount::bp_intersection_count_8x(int* bases_a, PackState* states_a, int size_a,
//         int* bases_b, PackState* states_b, int size_b)
// {
//     int i = 0, j = 0, res = 0;
//     int os_a = size_a - (size_a & 7);
//     int os_b = size_b - (size_b & 7);
//     uint64_t bits[4] __attribute__((aligned(32)));

//     while (i < os_a && j < os_b) {
//         __m256i base_a = _mm256_lddqu_si256((__m256i*)(bases_a + i));
//         __m256i base_b = _mm256_lddqu_si256((__m256i*)(bases_b + j));
//         __m256i state_a = _mm256_lddqu_si256((__m256i*)(states_a + i));
//         __m256i state_b = _mm256_lddqu_si256((__m256i*)(states_b + j));

//         int a_max = bases_a[i + 7];
//         int b_max = bases_b[j + 7];
//         i += (a_max <= b_max) * 8;
//         j += (b_max <= a_max) * 8;

//         __m256i and_all = _mm256_setzero_si256();

//         // shift0:
//         __m256i cmp_mask0 = _mm256_cmpeq_epi32(base_a, base_b);
//         __m256i and_0 = _mm256_and_si256(state_a, state_b);
//         and_all = _mm256_or_si256(and_all, _mm256_and_si256(and_0, cmp_mask0));

//         // shift1:
//         __m256 base_sf1 = _mm256_permute_ps((__m256)base_b, cyclic_shift1);
//         __m256 state_sf1 = _mm256_permute_ps((__m256)state_b, cyclic_shift1);
//         __m256i cmp_mask1 = _mm256_cmpeq_epi32(base_a, (__m256i)base_sf1);
//         __m256i and_1 = _mm256_and_si256(state_a, (__m256i)state_sf1);
//         and_all = _mm256_or_si256(and_all, _mm256_and_si256(and_1, cmp_mask1));
//         // shift2:
//         __m256 base_sf2 = _mm256_permute_ps((__m256)base_b, cyclic_shift2);
//         __m256 state_sf2 = _mm256_permute_ps((__m256)state_b, cyclic_shift2);
//         __m256i cmp_mask2 = _mm256_cmpeq_epi32(base_a, (__m256i)base_sf2);
//         __m256i and_2 = _mm256_and_si256(state_a, (__m256i)state_sf2);
//         and_all = _mm256_or_si256(and_all, _mm256_and_si256(and_2, cmp_mask2));

//         // shift3:
//         __m256 base_sf3 = _mm256_permute_ps((__m256)base_b, cyclic_shift3);
//         __m256 state_sf3 = _mm256_permute_ps((__m256)state_b, cyclic_shift3);
//         __m256i cmp_mask3 = _mm256_cmpeq_epi32(base_a, (__m256i)base_sf3);
//         __m256i and_3 = _mm256_and_si256(state_a, (__m256i)state_sf3);
//         and_all = _mm256_or_si256(and_all, _mm256_and_si256(and_3, cmp_mask3));

//         // shift4:
//         __m256i base_sf4 = _mm256_permute2f128_si256(base_b, base_b, 1);
//         __m256i state_sf4 = _mm256_permute2f128_si256(state_b, state_b, 1);
//         __m256i cmp_mask4 = _mm256_cmpeq_epi32(base_a, base_sf4);
//         __m256i and_4 = _mm256_and_si256(state_a, state_sf4);
//         and_all = _mm256_or_si256(and_all, _mm256_and_si256(and_4, cmp_mask4));

//         // shift5:
//         __m256 base_sf5 = _mm256_permute_ps((__m256)base_sf4, cyclic_shift1);
//         __m256 state_sf5 = _mm256_permute_ps((__m256)state_sf4, cyclic_shift1);
//         __m256i cmp_mask5 = _mm256_cmpeq_epi32(base_a, (__m256i)base_sf5);
//         __m256i and_5 = _mm256_and_si256(state_a, (__m256i)state_sf5);
//         and_all = _mm256_or_si256(and_all, _mm256_and_si256(and_5, cmp_mask5));

//         // shift6:
//         __m256 base_sf6 = _mm256_permute_ps((__m256)base_sf4, cyclic_shift2);
//         __m256 state_sf6 = _mm256_permute_ps((__m256)state_sf4, cyclic_shift2);
//         __m256i cmp_mask6 = _mm256_cmpeq_epi32(base_a, (__m256i)base_sf6);
//         __m256i and_6 = _mm256_and_si256(state_a, (__m256i)state_sf6);
//         and_all = _mm256_or_si256(and_all, _mm256_and_si256(and_6, cmp_mask6));

//         // shift7:
//         __m256 base_sf7 = _mm256_permute_ps((__m256)base_sf4, cyclic_shift3);
//         __m256 state_sf7 = _mm256_permute_ps((__m256)state_sf4, cyclic_shift3);
//         __m256i cmp_mask7 = _mm256_cmpeq_epi32(base_a, (__m256i)base_sf7);
//         __m256i and_7 = _mm256_and_si256(state_a, (__m256i)state_sf7);
//         and_all = _mm256_or_si256(and_all, _mm256_and_si256(and_7, cmp_mask7));

//         // popcnt:
//         _mm256_store_si256((__m256i*)bits, and_all);
//         res += __builtin_popcountl(bits[0]);
//         res += __builtin_popcountl(bits[1]);
//         res += __builtin_popcountl(bits[2]);
//         res += __builtin_popcountl(bits[3]);
//         // res += popcnt((const void*) bits, 32);
//     }

//     // res += bp_intersection_count_4x(bases_a + i, states_a + i, size_a - i,
//     //         bases_b + j, states_b + j, size_b - j);

//     while (i < size_a && j < size_b) {
//         if (bases_a[i] == bases_b[j]) {
//             res += __builtin_popcount(states_a[i] & states_b[j]);
//             i++; j++;
//         } else if (bases_a[i] < bases_b[j]) {
//             i++;
//         } else {
//             j++;
//         }
//     }

//     return res;
// }

// inline int BPTriangleCount::bp_intersection_count_8x2(int* bases_a, PackState* states_a, int size_a,
//         int* bases_b, PackState* states_b, int size_b)
// {
//     int i = 0, j = 0, res = 0;
//     int os_a = size_a - (size_a & 7);
//     int os_b = size_b - (size_b & 7);
//     uint64_t bits[8] __attribute__((aligned(32)));

//     while (i < os_a && j < os_b) {
//         __m256i base_a_0 = _mm256_lddqu_si256((__m256i*)(bases_a + i));
//         __m256i base_b_0 = _mm256_lddqu_si256((__m256i*)(bases_b + j));
//         __m256i state_a_0 = _mm256_lddqu_si256((__m256i*)(states_a + i));
//         __m256i state_b_0 = _mm256_lddqu_si256((__m256i*)(states_b + j));

//         int a_max = bases_a[i + 7];
//         int b_max = bases_b[j + 7];
//         i += (a_max <= b_max) * 8;
//         j += (b_max <= a_max) * 8;

//         // shift0:
//         __m256i mask0 = _mm256_cmpeq_epi32(base_a_0, base_b_0);

//         // shift1:
//         __m256 base_b_1 = _mm256_permute_ps((__m256)base_b_0, cyclic_shift1);
//         __m256 state_b_1 = _mm256_permute_ps((__m256)state_b_0, cyclic_shift1);
//         __m256i mask1 = _mm256_cmpeq_epi32(base_a_0, (__m256i)base_b_1);

//         // shift2:
//         __m256 base_b_2 = _mm256_permute_ps((__m256)base_b_0, cyclic_shift2);
//         __m256 state_b_2 = _mm256_permute_ps((__m256)state_b_0, cyclic_shift2);
//         __m256i mask2 = _mm256_cmpeq_epi32(base_a_0, (__m256i)base_b_2);

//         // shift3:
//         __m256 base_b_3 = _mm256_permute_ps((__m256)base_b_0, cyclic_shift3);
//         __m256 state_b_3 = _mm256_permute_ps((__m256)state_b_0, cyclic_shift3);
//         __m256i mask3 = _mm256_cmpeq_epi32(base_a_0, (__m256i)base_b_3);

//         __m256i and_1 = _mm256_or_si256(
//             _mm256_or_si256(
//                 _mm256_and_si256(mask0, _mm256_and_si256(state_a_0, state_b_0)),
//                 _mm256_and_si256(mask1, _mm256_and_si256(state_a_0, (__m256i)state_b_1))
//             ),
//             _mm256_or_si256(
//                 _mm256_and_si256(mask2, _mm256_and_si256(state_a_0, (__m256i)state_b_2)),
//                 _mm256_and_si256(mask3, _mm256_and_si256(state_a_0, (__m256i)state_b_3))
//             )
//         );

//         // _mm256_store_si256((__m256i*)bits, and_1);        
//         // res += __builtin_popcountl(bits[0]);
//         // res += __builtin_popcountl(bits[1]);
//         // res += __builtin_popcountl(bits[2]);
//         // res += __builtin_popcountl(bits[3]);
//         _mm256_store_si256((__m256i*)bits, and_1);

//         // shift4:
//         __m256i base_a_1 = _mm256_permute2f128_si256(base_a_0, base_a_0, 1);
//         __m256i state_a_1 = _mm256_permute2f128_si256(state_a_0, state_a_0, 1);
//         __m256i mask4 = _mm256_cmpeq_epi32(base_a_1, base_b_0);

//         // shift5:
//         __m256i mask5 = _mm256_cmpeq_epi32(base_a_1, (__m256i)base_b_1);

//         // shift6:
//         __m256i mask6 = _mm256_cmpeq_epi32(base_a_1, (__m256i)base_b_2);

//         // shift7:
//         __m256i mask7 = _mm256_cmpeq_epi32(base_a_1, (__m256i)base_b_3);

//         __m256i and_2 = _mm256_or_si256(
//             _mm256_or_si256(
//                 _mm256_and_si256(mask4, _mm256_and_si256(state_a_1, state_b_0)),
//                 _mm256_and_si256(mask5, _mm256_and_si256(state_a_1, (__m256i)state_b_1))
//             ),
//             _mm256_or_si256(
//                 _mm256_and_si256(mask6, _mm256_and_si256(state_a_1, (__m256i)state_b_2)),
//                 _mm256_and_si256(mask7, _mm256_and_si256(state_a_1, (__m256i)state_b_3))
//             )
//         );
//         // popcnt:
//         // _mm256_store_si256((__m256i*)bits, and_2);
//         // res += __builtin_popcountl(bits[0]);
//         // res += __builtin_popcountl(bits[1]);
//         // res += __builtin_popcountl(bits[2]);
//         // res += __builtin_popcountl(bits[3]);

//         _mm256_store_si256((__m256i*)(bits + 4), and_2);
//         res += popcnt((const void*) bits, 64);
//     }

//     while (i < size_a && j < size_b) {
//         if (bases_a[i] == bases_b[j]) {
//             res += __builtin_popcount(states_a[i] & states_b[j]);
//             i++; j++;
//         } else if (bases_a[i] < bases_b[j]) {
//             i++;
//         } else {
//             j++;
//         }
//     }

//     return res;
// }