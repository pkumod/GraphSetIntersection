#include "org_triangle_count.hpp"

OrgTriangleCount::OrgTriangleCount()
{
    v_num = 0;
    e_num = 0;
    align_malloc((void**)&pool_edges, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
}

OrgTriangleCount::~OrgTriangleCount()
{
    free(pool_edges);
}

void OrgTriangleCount::build(const EdgeVector& _e_v)
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

    int cur_node_idx = 0;
    int prev_u = -1;
    for (auto& e : rev_edge_vec) {
        if (e.first != prev_u) {
            prev_u = e.first;
            // while (cur_node_idx % 4 != 0) cur_node_idx++;
            graph[e.first].start = cur_node_idx;
        }
        graph[e.first].deg++;
        pool_edges[cur_node_idx++] = e.second;       
    }

    printf("v_num=%d e_num=%lld\n", v_num, e_num);
}

// extern double intersect_time = 0.0;
// extern unsigned long long intersect_cnt = 0;
int OrgTriangleCount::triangle_count()
{
    int res = 0;    

    // struct timeval time_start;
    // struct timeval time_end;
    // int count = 0;
    for (const auto& e : edge_vec) {
        const UVertex& u = graph[e.first];
        const UVertex& v = graph[e.second];

        // intersect_cnt++;
        // gettimeofday(&time_start, NULL);
#if SIMD_STATE == 2
        res += intersect_scalar2x_count(pool_edges + u.start, u.deg,
                pool_edges + v.start, v.deg);
#elif SIMD_STATE == 4
        res += intersect_simd4x_count(pool_edges + u.start, u.deg,
                pool_edges + v.start, v.deg);
#else
        res += intersect_count(pool_edges + u.start, u.deg,
                pool_edges + v.start, v.deg);
#endif
        // if (++count % 10000 == 0) printf("count=%d\n", count);
    // gettimeofday(&time_end, NULL);
    // intersect_time += (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    } 
    // printf("intersect_cnt=%llu intersect_time=%.3fms\n", intersect_cnt, intersect_time);

    return res;
}

// inline int OrgTriangleCount::intersection_count_8x(int* set_a, int size_a, int* set_b, int size_b)
// {
//     int i = 0, j = 0, res = 0;
//     int qs_a = size_a - (size_a & 7);
//     int qs_b = size_b - (size_b & 7);

//     while (i < qs_a && j < qs_b) {
//         __m256i v_a = _mm256_lddqu_si256((__m256i*)(set_a + i));
//         __m256i v_b = _mm256_lddqu_si256((__m256i*)(set_b + j));

//         int a_max = set_a[i + 7];
//         int b_max = set_b[j + 7];
//         i += (a_max <= b_max) * 8;
//         j += (b_max <= a_max) * 8;

//         constexpr int cyclic_shift1 = _MM_SHUFFLE(0,3,2,1);
//         constexpr int cyclic_shift2 = _MM_SHUFFLE(2,1,0,3);
//         constexpr int cyclic_shift3 = _MM_SHUFFLE(1,0,3,2);

//         __m256i cmp_mask0 = _mm256_cmpeq_epi32(v_a, v_b); 
//         __m256  rot1 = _mm256_permute_ps((__m256)v_b, cyclic_shift1);
//         __m256i cmp_mask1 = _mm256_cmpeq_epi32(v_a, (__m256i)rot1);
//         __m256  rot2 = _mm256_permute_ps((__m256)v_b, cyclic_shift2);
//         __m256i cmp_mask2 = _mm256_cmpeq_epi32(v_a, (__m256i)rot2);
//         __m256  rot3 = _mm256_permute_ps((__m256)v_b, cyclic_shift3);
//         __m256i cmp_mask3 = _mm256_cmpeq_epi32(v_a, (__m256i)rot3);

//         __m256 rot4 = _mm256_permute2f128_ps((__m256)v_b, (__m256)v_b, 1);
//         __m256i cmp_mask4 = _mm256_cmpeq_epi32(v_a, (__m256i)rot4);
//         __m256 rot5 = _mm256_permute_ps(rot4, cyclic_shift1);
//         __m256i cmp_mask5 = _mm256_cmpeq_epi32(v_a, (__m256i)rot5);
//         __m256 rot6 = _mm256_permute_ps(rot4, cyclic_shift2);
//         __m256i cmp_mask6 = _mm256_cmpeq_epi32(v_a, (__m256i)rot6);
//         __m256 rot7 = _mm256_permute_ps(rot4, cyclic_shift3);
//         __m256i cmp_mask7 = _mm256_cmpeq_epi32(v_a, (__m256i)rot7);

//         __m256i cmp_mask = _mm256_or_si256(
//             _mm256_or_si256(
//                 _mm256_or_si256(cmp_mask0, cmp_mask1),
//                 _mm256_or_si256(cmp_mask2, cmp_mask3)
//             ),
//             _mm256_or_si256(
//                 _mm256_or_si256(cmp_mask4, cmp_mask5),
//                 _mm256_or_si256(cmp_mask6, cmp_mask7)
//             )
//         );
//         int32_t mask = _mm256_movemask_ps((__m256)cmp_mask);
//         res += _mm_popcnt_u32(mask);
//     }

//     while (i < size_a && j < size_b) {
//         if (set_a[i] == set_b[j]) {
//             res++;
//             i++; j++;
//         } else if (set_a[i] < set_b[j]) {
//             i++;
//         } else {
//             j++;
//         }
//     }
          
//     return res;
// }

// inline int OrgTriangleCount::intersection_count_8x2(int* set_a, int size_a, int* set_b, int size_b)
// {
//     int i = 0, j = 0, res = 0;
//     int os_a = size_a - (size_a & 7);
//     int os_b = size_b - (size_b & 7);

//     while (i < os_a && j < os_b) {
//         __m256i v_a_0 = _mm256_lddqu_si256((__m256i*)(set_a + i));
//         __m256i v_b_0 = _mm256_lddqu_si256((__m256i*)(set_b + j));

//         int a_max = set_a[i + 7];
//         int b_max = set_b[j + 7];
//         i += (a_max <= b_max) * 8;
//         j += (b_max <= a_max) * 8;

//         constexpr int cyclic_shift1 = _MM_SHUFFLE(0,3,2,1);
//         constexpr int cyclic_shift2 = _MM_SHUFFLE(2,1,0,3);
//         constexpr int cyclic_shift3 = _MM_SHUFFLE(1,0,3,2);

//         __m256i cmp_mask0 = _mm256_cmpeq_epi32(v_a_0, v_b_0); 
//         __m256  v_b_1 = _mm256_permute_ps((__m256)v_b_0, cyclic_shift1);
//         __m256i cmp_mask1 = _mm256_cmpeq_epi32(v_a_0, (__m256i)v_b_1);
//         __m256  v_b_2 = _mm256_permute_ps((__m256)v_b_0, cyclic_shift2);
//         __m256i cmp_mask2 = _mm256_cmpeq_epi32(v_a_0, (__m256i)v_b_2);
//         __m256  v_b_3 = _mm256_permute_ps((__m256)v_b_0, cyclic_shift3);
//         __m256i cmp_mask3 = _mm256_cmpeq_epi32(v_a_0, (__m256i)v_b_3);

//         __m256i v_a_1 = _mm256_permute2f128_si256(v_a_0, v_a_0, 1);
//         __m256i cmp_mask4 = _mm256_cmpeq_epi32(v_a_1, (__m256i)v_b_0);
//         __m256i cmp_mask5 = _mm256_cmpeq_epi32(v_a_1, (__m256i)v_b_1);
//         __m256i cmp_mask6 = _mm256_cmpeq_epi32(v_a_1, (__m256i)v_b_2);
//         __m256i cmp_mask7 = _mm256_cmpeq_epi32(v_a_1, (__m256i)v_b_3);

//         __m256i cmp_mask_l = _mm256_or_si256(
//                 _mm256_or_si256(cmp_mask0, cmp_mask1),
//                 _mm256_or_si256(cmp_mask2, cmp_mask3)
//             );
//         __m256i cmp_mask_h = _mm256_or_si256(
//                 _mm256_or_si256(cmp_mask4, cmp_mask5),
//                 _mm256_or_si256(cmp_mask6, cmp_mask7)
//             );
 
//         int32_t mask_l = _mm256_movemask_ps((__m256)cmp_mask_l);
//         int32_t mask_h = _mm256_movemask_ps((__m256)cmp_mask_h);
//         res += _mm_popcnt_u32(mask_l);
//         res += _mm_popcnt_u32(mask_h);
//     }

//     while (i < size_a && j < size_b) {
//         if (set_a[i] == set_b[j]) {
//             res++;
//             i++; j++;
//         } else if (set_a[i] < set_b[j]) {
//             i++;
//         } else {
//             j++;
//         }
//     }
          
//     return res;
// }