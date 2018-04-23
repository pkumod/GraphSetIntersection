#include <immintrin.h>
#include <sys/time.h>
#include <cstdio>
#include <string>
#include <algorithm>
#include "util.hpp"
#include "org_maximal_clique.hpp"
#include "bitpack_maximal_clique.hpp"
#include "roaring_maximal_clique.hpp"
using namespace std;

struct timeval time_start;
struct timeval time_end;

string graph_file_path = "../data/youtube_cont_GRO.txt";

// OrgMaximalClique mc;
// RoaringMaximalClique mc;
BPMaximalClique mc;

int main(int argc, char* argv[])
{
    if (argc > 1) graph_file_path = std::string(argv[1]);
    
    auto edge_vec = load_graph(graph_file_path);

    printf("load_graph done.\n");

    gettimeofday(&time_start, NULL);
    mc.build(edge_vec);

    gettimeofday(&time_end, NULL);
    double build_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("build_time=%.3fms\n", build_time);  
    // mc.degeneracy_order();
    gettimeofday(&time_start, NULL);
    // int mc_num = mc.maximal_clique_bk();
    // int mc_num  = mc.maximal_clique_pivot();
    int mc_num = mc.maximal_clique_degen();
    gettimeofday(&time_end, NULL);
    double list_mc_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    // if (list_mc_time > 10000.0) printf("list_mc_time=%.3fs\n", list_mc_time / 1000.0);
    // else printf("list_mc_time=%.3fms\n", list_mc_time);

    printf("MC on %s:\n", graph_file_path.c_str());
    printf("mc_num=%d time=%.3fms\n", mc_num, list_mc_time);
    printf("cmp_cnt=%llu\n", cmp_cnt);

    // printf("intersect_cnt=%llu intersect_time=%.3fms\n", mc.intersect_cnt, mc.intersect_time);
    // double intersect_ratio = mc.intersect_time / list_mc_time;
    // printf("intersect_ratio=%.3f\n", intersect_ratio);
    
    if (argc > 2) mc.save_answers(argv[2]);

    // printf("All done.\n");

    // printf("cmp_cnt=%d, no_match_cnt=%d byte_check_cnt=", cmp_cnt, no_match_cnt);
    // for (int i = 0; i < 4; ++i) printf("%d ", byte_check_cnt[i]);
    // printf("\n");
    
    // double multimatch_proportion = (double) multimatch_cnt / cmp_cnt * 100.0;
    // double skew_proportion = (double) skew_cnt / inter_cnt * 100.0;
    // double low_select_proportion = (double) low_select_cnt / (inter_cnt - skew_cnt) * 100.0;
    // double nomatch_proportion = (double) no_match_cnt / cmp_cnt * 100.0;
    // printf("multi_pro= %.2f\nskew_pro= %.2f\nlow_select_pro=%.2f\n",
    //         multimatch_proportion, skew_proportion, low_select_proportion);
    // printf("multi_pro= %.2f\nnomatch_pro=%.2f\n", multimatch_proportion, nomatch_proportion);

    // uint8_t *elems = (uint8_t *)shuffle_mask;
    // for (int i = 0; i < 16; ++i) {
    //     // printf("%d: ", i);
    //     for (int j = 0; j < 16; ++j)
    //         printf("%d, ", elems[i * 16 + j]);
    //     printf("\n");
    // }
    
    return 0;
}