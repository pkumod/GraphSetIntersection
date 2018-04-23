#include "util.hpp"
#include "bitpack_triangle_count.hpp"
#include "org_triangle_count.hpp"
#include "roaring_triangle_count.hpp"
// #include "tbb/task_scheduler_init.h"
using namespace std;


struct timeval time_start;
struct timeval time_end;

string graph_file_path = "../data/youtube_cont_GRO.txt";

BPTriangleCount tc;
// RoaringTriangleCount tc;
// OrgTriangleCount tc;

int main(int argc, char* argv[])
{
    if (argc > 1) graph_file_path = std::string(argv[1]);

    auto edge_vec = load_graph(graph_file_path);
    printf("load_graph done.\n");

    gettimeofday(&time_start, NULL);
    tc.build(edge_vec);
    gettimeofday(&time_end, NULL);
    double build_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("build_time=%.3fms\n", build_time);
    
    gettimeofday(&time_start, NULL);
    int triangle_num = 0;
    // triangle_num = tc.triangle_count();
    // for (int i = 0; i < 10; ++i) 
        triangle_num = tc.triangle_count();
        // if (thread_num == 0) triangle_num = tc.triangle_count();
        // else {
        //     // tbb::task_scheduler_init init(thread_num);
        //     triangle_num = tc.triangle_count_mt(thread_num);
        // }
    gettimeofday(&time_end, NULL);
    double triangle_count_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    // triangle_count_time /= 10;
    printf("TC on %s:\n", graph_file_path.c_str());
    printf("triangle_num=%d time=%.3fms\n", triangle_num, triangle_count_time);
    printf("cmp_cnt=%llu\n", cmp_cnt);
    // double intersect_ratio = intersect_time / triangle_count_time;
    // printf("intersect_ratio=%.3f\n", intersect_ratio);
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

    return 0;
}