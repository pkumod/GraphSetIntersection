#include "util.hpp"
#include "org_subgraph_match.hpp"

using namespace std;

struct timeval time_start;
struct timeval time_end;

string graph_edges_file = "../data/test_subgraph_match/youtube_undirected_GRO.txt";
string graph_labels_file = "../data/test_subgraph_match/youtube_label_GRO.txt";
string sm_queries_file = "../data/test_subgraph_match/youtube_undirected_GRO.sm.queries";

OrgSubGraphMatch sm;

vector<LabelSubgraph> load_queries(const std::string path)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == NULL) {
        std::cout << "fail to open " << path << std::endl;
        quit();
    }

    vector<LabelSubgraph> queries;
    int v_num, e_num;
    while (fscanf(fp, "%d%d", &v_num, &e_num) != EOF) {
        LabelSubgraph q(v_num, e_num);
        q.vertex2label.reserve(v_num);
        q.edge_vec.reserve(e_num);
        int l, u, v;
        for (int i = 0; i < v_num; ++i) {
            fscanf(fp, "%d", &l);
            q.vertex2label.push_back(l);
        }
        for (int i = 0; i < e_num; ++i) {
            fscanf(fp, "%d%d", &u, &v);
            if (u > v) swap(u, v);
            q.edge_vec.push_back(Edge(u, v));
        }
        queries.push_back(q);            
    }

    return queries;
}

vector<int> load_labels(const std::string path)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == NULL) {
        std::cout << "fail to open " << path << std::endl;
        quit();
    }

    vector<int> labels;
    int u, l;
    while (fscanf(fp, "%d%d", &u, &l) != EOF) labels.push_back(l);

    return labels;
}

void save_answers(const char* file_path, const vector<vector<vector<int>>>& answers)
{
    FILE *fp = fopen(file_path, "w");
    if (fp == NULL) {
        std::cout << "fail to create " << file_path << std::endl;
        quit();
    }

    int ans_cnt = 0;
    for (const auto& ans : answers) {
        fprintf(fp, "ans%d=%lu\n", ans_cnt++, ans.size());
        // for (const auto& rec : ans) {
        //     for (const auto& u : rec)
        //         fprintf(fp, "%d ", u);
        //     fprintf(fp, "\n");
        // }        
    }
}

int main(int argc, char* argv[])
{
    if (argc > 3) {
        graph_edges_file = string(argv[1]);
        graph_labels_file = string(argv[2]);
        sm_queries_file = string(argv[3]);
    }

    auto edge_vec = load_graph(graph_edges_file);
    auto labels = load_labels(graph_labels_file);
    auto queries = load_queries(sm_queries_file);
    printf("load_graph done.\n");

    sm.build(edge_vec, labels);
    printf("build done.\n");
    
    vector<vector<vector<int>>> answers;
    gettimeofday(&time_start, NULL);
    for (const auto& q : queries) {
        auto ans = sm.subgraph_matching(q);
        answers.push_back(ans);
    }
    gettimeofday(&time_end, NULL);
    double query_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    if (query_time > 10000.0) printf("query_time=%.3fs\n", query_time / 1000.0);
    else printf("query_time=%.3fms\n", query_time);

    if (argc > 4) save_answers(argv[4], answers);

    return 0;
}