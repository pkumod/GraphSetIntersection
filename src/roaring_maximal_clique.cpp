#include "roaring_maximal_clique.hpp"

RoaringMaximalClique::RoaringMaximalClique()
{
    v_num = 0;
    e_num = 0;
}

void RoaringMaximalClique::build(const EdgeVector& _e_v)
{
    edge_vec.reserve(_e_v.size());
    for (auto& e : _e_v) if (e.first != e.second) edge_vec.push_back(e);
    // edge_vec = _e_v;

    std::sort(edge_vec.begin(), edge_vec.end(), edge_idpair_cmp);
    edge_vec.erase(std::unique(edge_vec.begin(), edge_vec.end()), edge_vec.end());

    for (auto& e : edge_vec) {
        v_num = std::max(v_num, e.first);
        v_num = std::max(v_num, e.second);
    }
    v_num++;
    e_num = (long long)edge_vec.size();

    graph.resize(v_num);

    for (auto& e : edge_vec) graph[e.first].add(e.second);

    printf("v_num=%d e_num=%lld\n", v_num, e_num);
}

int RoaringMaximalClique::maximal_clique_degen()
{
    maximum_clique_size = 0;
    mc_num = 0;

    bool *visited = new bool[v_num];
    memset(visited, 0, v_num * sizeof(bool));
    std::vector<int> dorder = degeneracy_order();

    Roaring R;
    for (auto v : dorder) {
        Roaring P, X;
        R.add(v);
        for(auto iter = graph[v].begin(); iter != graph[v].end(); iter++) {
            int u = *iter;
            if (visited[u] == 0)
                P.add(u);
            else
                X.add(u);
        }

        Tomita(R, P, X);
        visited[v] = 1;
        R.remove(v);
    }

    printf("maximum_clique_size=%d\n", maximum_clique_size);   

    delete []visited;

    return mc_num;
}

void RoaringMaximalClique::Tomita(Roaring& R, Roaring& P, Roaring& X)
{
    if (P.isEmpty()) { // report R as a maximal clique.
        if (X.isEmpty()) {
            maximal_cliques.push_back(R);
            mc_num++;
            maximum_clique_size = std::max(maximum_clique_size, (int)R.cardinality());
            // printf("R=%s\n", R.toString().c_str());
            // printf("mc_num=%d\n", mc_num);
        }
        return;
    }

    // heuristic 3: choose the first vertex from P as the pivot.
    int u;
    if (X.isEmpty()) u = P.minimum();
    else u = X.minimum();    

    Roaring tmpP = P - graph[u];
    // only need to enumerate verices in P - N(u).
    for(auto iter = tmpP.begin(); iter != tmpP.end(); iter++) {
        int v = *iter;
        // if (graph[u].contains(v)) continue; // skip vertices in N(u).
        R.add(v);
        Roaring newP = (P & graph[v]);
        Roaring newX = (X & graph[v]);
        Tomita(R, newP, newX);

        // move v from P to X.
        P.remove(v);
        X.add(v);
        R.remove(v);        
    }
}

std::vector<int> RoaringMaximalClique::degeneracy_order()
{
    std::vector<int> deg(v_num);
    int md = 0;
    for (int i = 0; i < v_num; ++i) {
        deg[i] = graph[i].cardinality();
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
        for(auto iter = graph[v].begin(); iter != graph[v].end(); iter++) {
            int u = *iter;
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
    
    printf("degeneracy=%d\n", degeneracy);

    return degen_order;
}

void RoaringMaximalClique::save_answers(const char* file_path)
{
    FILE *fp = fopen(file_path, "w");
    if (fp == NULL) {
        std::cout << "fail to create " << file_path << std::endl;
        quit();
    }

    for (const auto& mc : maximal_cliques) {
        fprintf(fp, "%s\n", mc.toString().c_str());
    }

    fclose(fp);
}