#include "bitpack_maximal_clique.hpp"

BPMaximalClique::BPMaximalClique()
{
    v_num = 0;
    e_num = 0;
    align_malloc((void**)&pool_base, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
    align_malloc((void**)&pool_state, 32, sizeof(PackState) * PACK_NODE_POOL_SIZE);
    align_malloc((void**)&pool_mc, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
}

BPMaximalClique::~BPMaximalClique()
{
    free(pool_state);
    free(pool_base);
    free(pool_mc);
}

void BPMaximalClique::build(const EdgeVector& _e_v)
{
    edge_vec.reserve(_e_v.size());
    for (auto& e : _e_v) if (e.first != e.second) edge_vec.push_back(e);

    std::sort(edge_vec.begin(), edge_vec.end(), edge_idpair_cmp);
    edge_vec.erase(std::unique(edge_vec.begin(), edge_vec.end()), edge_vec.end());

    for (auto& e : edge_vec) {
        v_num = std::max(v_num, e.first);
        v_num = std::max(v_num, e.second);
    }
    v_num++;
    e_num = (long long)edge_vec.size();
    p_num = v_num / PACK_WIDTH;
    if (v_num % PACK_WIDTH != 0) p_num++;

    graph.resize(v_num);
    org_deg.reserve(v_num);
    for (int i = 0; i < v_num; ++i) org_deg[i] = 0;

    int cur_packnode_idx = -1;
    int prev_u = -1;
    for (auto& e : edge_vec) {
        int v_base = (e.second >> PACK_SHIFT);
        PackState v_bit = ((PackState)1 << (e.second & PACK_MASK));
        org_deg[e.first]++;
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

int BPMaximalClique::maximal_clique_degen()
{
    align_malloc((void**)&sets_base, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
    align_malloc((void**)&sets_state, 32, sizeof(PackState) * PACK_NODE_POOL_SIZE);
    max_pool_sets_idx = 0; maximum_clique_size = 0;
    pool_mc_idx = 0; mc_num = 0;
    intersect_call_time = 0; big_intersect_call_time = 0;
    
    PackState *visited_state = new PackState[p_num];
    memset(visited_state, 0, sizeof(PackState) * p_num);    
    std::vector<int> dorder = degeneracy_order();
    
    std::vector<int> R;
    R.reserve(2048);
    R.push_back(-1);
    for (auto v : dorder) {
        R[0] = v;
        UVertex P(0, 0);
#if SIMD_STATE == 4
        P.deg = bp_subtract_visited_simd4x(pool_base + graph[v].start,
            pool_state + graph[v].start, graph[v].deg,
            visited_state, sets_base + P.start, sets_state + P.start);
#else
        P.deg = bp_subtract_visited(pool_base + graph[v].start,
            pool_state + graph[v].start, graph[v].deg,
            visited_state, sets_base + P.start, sets_state + P.start);
#endif
        UVertex X(P.deg, 0);
#if SIMD_STATE == 4
        X.deg = bp_subtract_unvisited_simd4x(pool_base + graph[v].start,
            pool_state + graph[v].start, graph[v].deg,
            visited_state, sets_base + X.start, sets_state + X.start);
#else
        X.deg = bp_subtract_unvisited(pool_base + graph[v].start,
            pool_state + graph[v].start, graph[v].deg,
            visited_state, sets_base + X.start, sets_state + X.start);
#endif
        // printf("P.deg=%d X.deg=%d\n", P.deg, X.deg);
        Tomita(R, P, X);
        int v_base = (v >> PACK_SHIFT);
        PackState v_bit = ((PackState)1 << (v & PACK_MASK));
        visited_state[v_base] |= v_bit;
        // printf("v=%d done.\n", v);
    }
    R.pop_back();

    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    // printf("intersect_call_ratio=%.3f%%(%d/%d)\n",
    //     big_intersect_call_time * 100.0 / intersect_call_time,
    //     big_intersect_call_time, intersect_call_time);
       
    free(sets_base);
    free(sets_state);
    delete []visited_state;

    return mc_num;
}

void BPMaximalClique::Tomita(std::vector<int>& R, UVertex P, UVertex X)
{
    if (P.deg == 0) {
        if (X.deg == 0) {
            // std::sort(c.begin(), c.end());
            memcpy(pool_mc + pool_mc_idx, R.data(), R.size() * sizeof(int));
            pool_mc_idx += R.size();
            pool_mc[pool_mc_idx++] = -1; // split different maximal cliques by -1. 
            mc_num++;
            max_pool_sets_idx = std::max(max_pool_sets_idx, X.start + X.deg);
            maximum_clique_size = std::max(maximum_clique_size, (int)R.size());
        }
        return;
    }

    // heuristic 2: choose a pivot u from P or X to maximize |N(u)|.
    // int u = -1, max_deg = -1;
    // for (int i = 0; i < P.deg; ++i) {
    //     int cu_high = (sets_base[P.start + i] << PACK_SHIFT);
    //     int cu_state = sets_state[P.start + i];
    //     while (cu_state) {
    //         int cu = (cu_high | __builtin_ctz(cu_state));
    //         cu_state &= (cu_state - 1);
    //         // if (max_deg < graph[cu].deg) {
    //         //     max_deg = graph[cu].deg;
    //         //     u = cu;
    //         // }
    //         if (max_deg < org_deg[cu]) {
    //             max_deg = org_deg[cu];
    //             u = cu;
    //         }
    //     }
    // }
    // for (int i = 0; i < X.deg; ++i) {
    //     int cu_high = (sets_base[X.start + i] << PACK_SHIFT);
    //     int cu_state = sets_state[X.start + i];
    //     while (cu_state) {
    //         int cu = (cu_high | __builtin_ctz(cu_state));
    //         cu_state &= (cu_state - 1);
    //         // if (max_deg < graph[cu].deg) {
    //         //     max_deg = graph[cu].deg;
    //         //     u = cu;
    //         // }
    //         if (max_deg < org_deg[cu]) {
    //             max_deg = org_deg[cu];
    //             u = cu;
    //         }
    //     }        
    // }
    // heuristic 3: choose the first vertex from P as the pivot.
    int u;
    if (X.deg > 0) u = (sets_base[X.start] << PACK_SHIFT) | __builtin_ctz(sets_state[X.start]);
    else u = (sets_base[P.start] << PACK_SHIFT) | __builtin_ctz(sets_state[P.start]);

    int N_u_idx = graph[u].start, N_u_end = graph[u].start + graph[u].deg;

    R.push_back(-1); // enlarge R by a placeholder(-1).
    for (int i = 0; i < P.deg; ++i) {
        int v_base = sets_base[P.start + i];
        int v_state = sets_state[P.start + i];
        int v_high = (v_base << PACK_SHIFT);
        while (N_u_idx != N_u_end && pool_base[N_u_idx] < v_base) N_u_idx++;
        if (N_u_idx != N_u_end && pool_base[N_u_idx] == v_base) {
            v_state &= (~pool_state[N_u_idx]);
            N_u_idx++;
        }
                
        while (v_state) {
            int v = (v_high | __builtin_ctz(v_state));
            v_state &= (v_state - 1);

            R.back() = v;
            int newPstart = X.start + X.deg;

#if SIMD_STATE == 2
            int newPdeg = bp_intersect_scalar2x(sets_base + P.start, sets_state + P.start, P.deg,
                    pool_base + graph[v].start, pool_state + graph[v].start, graph[v].deg,
                    sets_base + newPstart, sets_state + newPstart);
#elif SIMD_STATE == 4
#if SIMD_MODE == 0
            int newPdeg = bp_intersect_simd4x(sets_base + P.start, sets_state + P.start, P.deg,
                    pool_base + graph[v].start, pool_state + graph[v].start, graph[v].deg,
                    sets_base + newPstart, sets_state + newPstart);
#else
            int newPdeg = bp_intersect_filter_simd4x(sets_base + P.start, sets_state + P.start, P.deg,
                    pool_base + graph[v].start, pool_state + graph[v].start, graph[v].deg,
                    sets_base + newPstart, sets_state + newPstart);
#endif            
#else
            int newPdeg = bp_intersect(sets_base + P.start, sets_state + P.start, P.deg,
                    pool_base + graph[v].start, pool_state + graph[v].start, graph[v].deg,
                    sets_base + newPstart, sets_state + newPstart);
#endif
            UVertex newP(newPstart, newPdeg);
            int newXstart = newPstart + newPdeg;
#if SIMD_STATE == 2
            int newXdeg = bp_intersect_scalar2x(sets_base + X.start, sets_state + X.start, X.deg,
                    pool_base + graph[v].start, pool_state + graph[v].start, graph[v].deg,
                    sets_base + newXstart, sets_state + newXstart);
#elif SIMD_STATE == 4
#if SIMD_MODE == 0
            int newXdeg = bp_intersect_simd4x(sets_base + X.start, sets_state + X.start, X.deg,
                    pool_base + graph[v].start, pool_state + graph[v].start, graph[v].deg,
                    sets_base + newXstart, sets_state + newXstart);
#else
            int newXdeg = bp_intersect_filter_simd4x(sets_base + X.start, sets_state + X.start, X.deg,
                    pool_base + graph[v].start, pool_state + graph[v].start, graph[v].deg,
                    sets_base + newXstart, sets_state + newXstart);
#endif           
#else
            int newXdeg = bp_intersect(sets_base + X.start, sets_state + X.start, X.deg,
                    pool_base + graph[v].start, pool_state + graph[v].start, graph[v].deg,
                    sets_base + newXstart, sets_state + newXstart);
#endif
            UVertex newX(newXstart, newXdeg);

            Tomita(R, newP, newX);

            // move v from P to X.
            PackState v_bit = ((PackState)1 << (v & PACK_MASK));
            sets_state[P.start + i] &= (~v_bit);
            X.deg = bp_merge_one(sets_base + X.start, sets_state + X.start, X.deg,
                    v_base, v_bit);
        }
    }
    R.pop_back();
}

std::vector<int> BPMaximalClique::degeneracy_order()
{
    std::vector<int> deg(v_num);
    int md = 0;
    for (int i = 0; i < v_num; ++i) {
        deg[i] = 0;
        for (int j = 0; j < graph[i].deg; ++j)
            deg[i] += _mm_popcnt_u32(pool_state[graph[i].start + j]);
        md = std::max(md, deg[i]);
    }

    std::vector<int> bin(md + 1);
    for (int i = 0; i <= md; ++i) bin[i] = 0;
    for (int i = 0; i < v_num; ++i) bin[deg[i]]++;

    int start = 0;
    for (int i = 0; i <= md; ++i) {
        int num = bin[i];
        bin[i] = start;
        start += num;
    }

    std::vector<int> vert(v_num), pos(v_num);
    for (int i = 0; i < v_num; ++i) {
        pos[i] = bin[deg[i]];
        vert[pos[i]] = i;
        bin[deg[i]]++;
    }
    for (int i = md; i > 0; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    std::vector<int> degen_order;
    degen_order.reserve(v_num);
    int degeneracy = 0;
    for (int i = 0; i < v_num; ++i) {
        int v = vert[i];
        degen_order.push_back(v);
        degeneracy = std::max(degeneracy, deg[v]);
        for (int j = 0; j < graph[v].deg; ++j) {
            int u_high = (pool_base[graph[v].start + j] << PACK_SHIFT);
            PackState u_state = pool_state[graph[v].start + j];
            while (u_state) {
                int u = (u_high | __builtin_ctz(u_state));
                u_state &= (u_state - 1);
                if (deg[u] > deg[v]) {
                    int du = deg[u], pu = pos[u];
                    int pw = bin[du], w = vert[pw];
                    if (u != w) {
                        pos[u] = pw; vert[pu] = w;
                        pos[w] = pu; vert[pw] = u;
                    }
                    bin[du]++;
                    deg[u]--;
                }
            }
        }
    }
    
    printf("degeneracy=%d\n", degeneracy);

    return degen_order;
}

void BPMaximalClique::save_answers(const char* file_path)
{
    FILE *fp = fopen(file_path, "w");
    if (fp == NULL) {
        std::cout << "fail to create " << file_path << std::endl;
        quit();
    }

    for (int i = 0; i < pool_mc_idx; ++i)
        if (pool_mc[i] == -1) fprintf(fp, "\n");
        else fprintf(fp, "%d ", pool_mc[i]);

    fclose(fp);
}
