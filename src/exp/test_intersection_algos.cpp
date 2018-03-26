#include "../util.hpp"
#include "random_data_generator.hpp"
#include "../intersection_algos.hpp"
#include "../roaring/roaring.hh"
using namespace std;

struct timeval time_start;
struct timeval time_end;

typedef int (* InterUINTAlgoFuncPtr)(int *, int,
                int *, int, int *);

typedef int (* InterBSRAlgoFuncPtr)(int *, int *, int,
                int *, int *, int, int *, int *);

struct InterAlgoFunc
{
    InterUINTAlgoFuncPtr uint_func_ptr;
    InterBSRAlgoFuncPtr bsr_func_ptr;
    int encoding; // 0: INT list; 1: BSR; 2:RoaringBitmap.
    std::string algo_name;    
    InterAlgoFunc(): uint_func_ptr(NULL), bsr_func_ptr(NULL),
            encoding(0), algo_name("null") {};
};

InterAlgoFunc get_intersection_algo_func(int algo_id)
{
    InterAlgoFunc intersection_algo_func;
       
    switch (algo_id) {
        case 0:
            intersection_algo_func.uint_func_ptr = &intersect_scalarmerge_uint;
            intersection_algo_func.algo_name ="ScalarMerge";
            break;
        case 1:
            intersection_algo_func.bsr_func_ptr = &intersect_scalarmerge_bsr;
            intersection_algo_func.algo_name ="ScalarMerge+BSR";
            intersection_algo_func.encoding = 1;
            break;
        case 10:
            intersection_algo_func.uint_func_ptr = &intersect_shuffle_uint_b4;
            intersection_algo_func.algo_name ="Shuffling";
            break;
        case 11:
            intersection_algo_func.bsr_func_ptr = &intersect_shuffle_bsr_b4;
            intersection_algo_func.algo_name ="Shuffling+BSR";
            intersection_algo_func.encoding = 1;
            break;
        case 20:
            intersection_algo_func.uint_func_ptr = &intersect_bmiss_uint_b4;
            intersection_algo_func.algo_name ="BMiss";
            break;
        case 21:
            intersection_algo_func.uint_func_ptr = &intersect_bmiss_uint_sttni_b8;
            intersection_algo_func.algo_name ="BMiss+STTNI";
            break;
        case 30:
            intersection_algo_func.uint_func_ptr = &intersect_hierainter_uint_sttni;
            intersection_algo_func.algo_name ="HieraInter";
            break;            
        case 40:
            intersection_algo_func.uint_func_ptr = &intersect_qfilter_uint_b4_v2;
            intersection_algo_func.algo_name ="QFilter";
            break;
        case 41:
            intersection_algo_func.bsr_func_ptr = &intersect_qfilter_bsr_b4_v2;
            intersection_algo_func.algo_name ="QFilter+BSR";
            intersection_algo_func.encoding = 1;
            break;
        case 50:
            intersection_algo_func.uint_func_ptr = &intersect_scalargalloping_uint;
            intersection_algo_func.algo_name ="ScalarGalloping";
            break;
        case 51:
            intersection_algo_func.bsr_func_ptr = &intersect_scalargalloping_bsr;
            intersection_algo_func.algo_name ="ScalarGalloping+BSR";
            intersection_algo_func.encoding = 1;
            break;
        case 60:
            intersection_algo_func.uint_func_ptr = &intersect_simdgalloping_uint;
            intersection_algo_func.algo_name ="SIMDGalloping";
            break;
        case 61:
            intersection_algo_func.bsr_func_ptr = &intersect_simdgalloping_bsr;
            intersection_algo_func.algo_name ="SIMDGalloping+BSR";
            intersection_algo_func.encoding = 1;
            break;
        case 70:
            intersection_algo_func.bsr_func_ptr = NULL;
            intersection_algo_func.algo_name ="Roaring";
            intersection_algo_func.encoding = 2;
            break;
        default:
            intersection_algo_func.uint_func_ptr = &intersect_scalarmerge_uint;
            intersection_algo_func.algo_name ="ScalarMerge";
            break;
    }
    return intersection_algo_func;  
}


void check_result(int* set_a, int size_a, int* set_b, int size_b,
        int* set_c_y, int size_c_y)
{
    int *set_c_x = NULL, size_c_x = 0;
    align_malloc((void**)&set_c_x, 32, sizeof(int) * std::min(size_a, size_b));
    size_c_x = intersect_scalarmerge_uint(set_a, size_a, set_b, size_b, set_c_x);
    bool same = (size_c_x == size_c_y);
    for (int i = 0; i < std::min(size_c_x, size_c_y); ++i)
        if (set_c_x[i] != set_c_y[i]) {
            same = false;
            break;
        }

    if (same) printf("same result!\n");
    else {
        printf("size_c_x=%d, size_c_y=%d\n", size_c_x, size_c_y);
        printf("set_c_x:\n");
        for (int i = 0; i < size_c_x; ++i)
            printf("%d, ", set_c_x[i]);
        printf("\n************\n");
        printf("set_c_y:\n");
        for (int i = 0; i < size_c_y; ++i)
            printf("%d, ", set_c_y[i]);
        printf("\n");        
    }
}

int main(int argc, char* argv[])
{ 
    // default parameters:
    int len = 4000000;
    double skew_ratio = 1.0;
    double selectivity = 0.1;
    double density = 0.01;

    bool check_res = false;
    int n = 3;
    int algo_id = 0;

    /*
        00: ScalarMerge;
        01: ScalarMerge +BSR;  
        10: Shuffling; (by Kastov)
        11: Shuffling+BSR;
        20: BMiss (block size = 4); (by Inoue)
        21: BMiss+STTNI (block size = 8); (by Inoue)
        30: HieraInter (with STTNI, by Schlegel)
        40: QFilter (by Han, our algorithm)
        41: QFilter+BSR
        50: ScalarGalloping;
        51: ScalarGalloping+BSR; 
        60: SIMDGalloping;
        61: SIMDGalloping+BSR;
    */

    int i;
    if ((i = arg_pos((char *)"-len", argc, argv)) > 0)
        len = atoi(argv[i + 1]);
    if ((i = arg_pos((char *)"-skew", argc, argv)) > 0)
        skew_ratio = atof(argv[i + 1]);
    if ((i = arg_pos((char *)"-select", argc, argv)) > 0)
        selectivity = atof(argv[i + 1]);    
    if ((i = arg_pos((char *)"-dense", argc, argv)) > 0)
        density = atof(argv[i + 1]);

    if ((i = arg_pos((char *)"-check", argc, argv)) > 0)
        check_res = (atoi(argv[i + 1]) == 1) ? true : false;
    if ((i = arg_pos((char *)"-n", argc, argv)) > 0)
        n = atoi(argv[i + 1]); 
    if ((i = arg_pos((char *)"-algo", argc, argv)) > 0)
        algo_id = atoi(argv[i + 1]);

    int *set_a = NULL, *set_b = NULL, size_a = 0, size_b = 0;    
    gen_id_list(len, skew_ratio, selectivity, density,
            set_a, set_b, size_a, size_b);
    int *set_c = NULL, size_c = 0;
    align_malloc((void**)&set_c, 32, sizeof(int) * std::min(size_a, size_b));

    InterAlgoFunc inter_algo = get_intersection_algo_func(algo_id);
    printf("%s\n", inter_algo.algo_name.c_str());

    if (inter_algo.encoding == 1) {
        int *bases_a, *states_a, *bases_b, *states_b, *bases_c, *states_c;
        int card_a = 0, card_b = 0, card_c = 0;
        align_malloc((void**)&bases_a, 32, sizeof(int) * size_a);
        align_malloc((void**)&states_a, 32, sizeof(int) * size_a);
        align_malloc((void**)&bases_b, 32, sizeof(int) * size_b);
        align_malloc((void**)&states_b, 32, sizeof(int) * size_b);
        align_malloc((void**)&bases_c, 32, sizeof(int) * std::min(size_a, size_b));
        align_malloc((void**)&states_c, 32, sizeof(int) * std::min(size_a, size_b));
        card_a = offline_uint_trans_bsr(set_a, size_a, bases_a, states_a);
        card_b = offline_uint_trans_bsr(set_b, size_b, bases_b, states_b);

        gettimeofday(&time_start, NULL);
        for (int t = 0; t < n; ++t) {
            card_c = (*(inter_algo.bsr_func_ptr))(bases_a, states_a, card_a,
                    bases_b, states_b, card_b, bases_c, states_c);
        }
        gettimeofday(&time_end, NULL);
        printf("card_a=%d card_b=%d card_c=%d\n", card_a, card_b, card_c);
        size_c = offline_bsr_trans_uint(bases_c, states_c, card_c, set_c);
    } else if (inter_algo.encoding == 2) {
        Roaring roar_a(size_a, (const uint32_t*)set_a);
        Roaring roar_b(size_b, (const uint32_t*)set_b);
        // roar_a.runOptimize(); roar_a.shrinkToFit();
        // roar_b.runOptimize(); roar_b.shrinkToFit();
        Roaring roar_c;
        gettimeofday(&time_start, NULL);
        for (int t = 0; t < n; ++t) {
            roar_c = roar_a & roar_b;
        }
        gettimeofday(&time_end, NULL);
        size_c = roar_c.cardinality();
        roar_c.toUint32Array((uint32_t*)set_c);
    } else {
        gettimeofday(&time_start, NULL);
        for (int t = 0; t < n; ++t) {
            size_c = (*(inter_algo.uint_func_ptr))(set_a, size_a,
                    set_b, size_b, set_c);           
        }
        gettimeofday(&time_end, NULL);        
    }
    
    printf("size_c=%d\n", size_c);
    double run_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 +
            (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    run_time /= n;
    double ele_per_usec = (size_a + size_b) / (run_time * 1000.0);
    printf("ele_per_usec=%.2f ", ele_per_usec);
    if (run_time > 10000.0) printf("run_time=%.3fs\n", run_time / 1000.0);
    else if (run_time < 0.1) printf("run_time=%.3fus\n", run_time * 1000.0);
    else printf("run_time=%.3fms\n", run_time);

    if (check_res) {
      check_result(set_a, size_a, set_b, size_b, set_c, size_c); 
    }
    
    return 0;
}