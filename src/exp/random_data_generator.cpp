#include "random_data_generator.hpp"
#include <unordered_set>

void gen_id_list(int len, double skew_ratio, double selectivity, double density,
    int*& set_a, int*& set_b, int& size_a, int& size_b)
{
    size_a = len;
    size_b = (int)(len / skew_ratio); 
    int size_c = (int) (len * selectivity); 
    int range = (int) (std::min((uint64_t)( (size_a + size_b - size_c) / density),
                (uint64_t)RAND_MAX));     

    align_malloc((void**)&set_a, 32, sizeof(int) * size_a);
    align_malloc((void**)&set_b, 32, sizeof(int) * size_b);

    std::unordered_set<int> ele_set;
    // srand(time(NULL));

    int x;
    for (int i = 0; i < size_a; ++i) {        
        do {
            x = rand() % range;
        } while (ele_set.find(x) != ele_set.end());
        ele_set.insert(x);
        set_a[i] = x;
    }

    for (int i = 0; i < size_b; ++i) {
        if (i < size_c) {
            set_b[i] = set_a[i];
        } else {
            do {
                x = rand() % range;
            } while (ele_set.find(x) != ele_set.end());
            ele_set.insert(x);
            set_b[i] = x;
        }
    }

    std::sort(set_a, set_a + size_a);
    std::sort(set_b, set_b + size_b);

    printf("size_a=%d, size_b=%d, size_c=%d, range=%d\n", size_a, size_b, size_c, range);
    printf("gen_id_list() done.\n");
}