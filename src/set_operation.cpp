#include "set_operation.hpp"

// #define likely(x)       __builtin_expect(!!(x), 1)
// #define unlikely(x)     __builtin_expect(!!(x), 0)

constexpr int cyclic_shift1 = _MM_SHUFFLE(0,3,2,1); //rotating right
constexpr int cyclic_shift2 = _MM_SHUFFLE(2,1,0,3); //rotating left
constexpr int cyclic_shift3 = _MM_SHUFFLE(1,0,3,2); //between
static const __m128i all_zero_si128 = _mm_setzero_si128();
static const __m128i all_one_si128 = _mm_set_epi32(0xffffffff, 0xffffffff,
        0xffffffff, 0xffffffff);
static const uint8_t shuffle_pi8_array[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 
    8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 
    12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 
    8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
};
// static const ShuffleDict shuffle_mask_arr;
static const __m128i *shuffle_mask = (__m128i*)(shuffle_pi8_array);

int * prepare_byte_check_mask_dict()
{
    int * mask = new int[65536];

    auto trans_c_s = [](const int c) -> int {
        switch (c) {
            case 0: return -1; // no match
            case 1: return 0;
            case 2: return 1;
            case 4: return 2;
            case 8: return 3;
            default: return 4; // multiple matches.
        }
    };

    for (int x = 0; x < 65536; ++x) {        
        int c0 = (x & 0xf), c1 = ((x >> 4) & 0xf);
        int c2 = ((x >> 8) & 0xf), c3 = ((x >> 12) & 0xf);
        int s0 = trans_c_s(c0), s1= trans_c_s(c1);
        int s2 = trans_c_s(c2), s3 = trans_c_s(c3);
        
        bool is_multiple_match = (s0 == 4) || (s1 == 4) ||
                (s2 == 4) || (s3 == 4);
        if (is_multiple_match) {
            mask[x] = -1;
            continue;
        }
        bool is_no_match = (s0 == -1) && (s1 == -1) &&
                (s2 == -1) && (s3 == -1);
        if (is_no_match) {
            mask[x] = -2;
            continue;
        }
        if (s0 == -1) s0 = 0; if (s1 == -1) s1 = 1;
        if (s2 == -1) s2 = 2; if (s3 == -1) s3 = 3;
        mask[x] = (s0) | (s1 << 2) | (s2 << 4) | (s3 << 6);        
    }

    return mask;
}
static const int *byte_check_mask_dict = prepare_byte_check_mask_dict();

uint8_t * prepare_match_shuffle_dict()
{
    uint8_t * dict = new uint8_t[4096];

    for (int x = 0; x < 256; ++x) {
        for (int i = 0; i < 4; ++i) {
            uint8_t c = (x >> (i << 1)) & 3; // c = 0, 1, 2, 3
            int pos = x * 16 + i * 4;
            for (uint8_t j = 0; j < 4; ++j)
                dict[pos + j] = c * 4 + j;
        }
    }

    return dict;
}
static const __m128i *match_shuffle_dict = (__m128i*)prepare_match_shuffle_dict();

static const uint8_t byte_check_group_a_pi8[64] = {
    0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 8, 8, 12, 12, 12, 12,
    1, 1, 1, 1, 5, 5, 5, 5, 9, 9, 9, 9, 13, 13, 13, 13,
    2, 2, 2, 2, 6, 6, 6, 6, 10, 10, 10, 10, 14, 14, 14, 14,
    3, 3, 3, 3, 7, 7, 7, 7, 11, 11, 11, 11, 15, 15, 15, 15,
};
static const uint8_t byte_check_group_b_pi8[64] = {
    0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12,
    1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13,
    2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 14,
    3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15, 3, 7, 11, 15,
};
static const __m128i *byte_check_group_a_order = (__m128i*)(byte_check_group_a_pi8);
static const __m128i *byte_check_group_b_order = (__m128i*)(byte_check_group_b_pi8);
// static const __m128i byte_check_group_perm_a = _mm_set_epi32(0x0c0c0c0c, 0x08080808,
//         0x04040404, 0x00000000);
// static const __m128i byte_check_group_perm_b = _mm_set_epi32(0x0c080400, 0x0c080400,
//         0x0c080400, 0x0c080400);
// static const __m128i all_one_pi8 = _mm_set_epi32(0x11111111, 0x11111111, 0x11111111, 0x11111111);

int intersect(int *set_a, int size_a, int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;   
}

int intersect_count(int *set_a, int size_a, int *set_b, int size_b)
{
    int i = 0, j = 0, res = 0;
    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            res++;
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return res;
}

int intersect_scalar2x(int *set_a, int size_a, int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int ds_a = size_a - (size_a & 1);
    int ds_b = size_b - (size_b & 1);

    if (ds_a > 0 && ds_b > 0) {
        while (true) {
            int Adat0 = set_a[i], Adat1 = set_a[i + 1];
            int Bdat0 = set_b[j], Bdat1 = set_b[j + 1];

            if (Adat0 == Bdat0) {
                set_c[size_c++] = Adat0;
            } else if (Adat0 == Bdat1) {
                set_c[size_c++] = Adat0;
                goto advanceB;
            } else if (Adat1 == Bdat0) {
                set_c[size_c++] = Adat1;
                goto advanceA;
            }

            if (Adat1 == Bdat1) {
                set_c[size_c++] = Adat1;
                goto advanceAB;
            } else if (Adat1 > Bdat1)  {
                goto advanceB;
            } else {
                goto advanceA;
            }

    advanceA:
            i += 2;
            if (i >= ds_a) break;
            else continue;
    advanceB:
            j += 2;
            if (j >= ds_b) break;
            else continue;
    advanceAB:
            i += 2; j += 2;
            if (i >= ds_a || j >= ds_b) break;
            else continue;
        }        
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int intersect_scalar2x_count(int* set_a, int size_a, int* set_b, int size_b)
{
    int i = 0, j = 0, res = 0;
    int ds_a = size_a - (size_a & 1);
    int ds_b = size_b - (size_b & 1);

    if (ds_a > 0 && ds_b > 0) {
        while (true) {
            int Adat0 = set_a[i], Adat1 = set_a[i + 1];
            int Bdat0 = set_b[j], Bdat1 = set_b[j + 1];

            if (Adat0 == Bdat0) {
                res++;
            } else if (Adat0 == Bdat1) {
                res++;
                goto advanceB;
            } else if (Adat1 == Bdat0) {
                res++;
                goto advanceA;
            }

            if (Adat1 == Bdat1) {
                res++;
                goto advanceAB;
            } else if (Adat1 > Bdat1)  {
                goto advanceB;
            } else {
                goto advanceA;
            }

    advanceA:
            i += 2;
            if (i >= ds_a) break;
            else continue;
    advanceB:
            j += 2;
            if (j >= ds_b) break;
            else continue;
    advanceAB:
            i += 2; j += 2;
            if (i >= ds_a || j >= ds_b) break;
            else continue;
        }        
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            res++;
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }
    
    return res;
}

int intersect_simd4x(int *set_a, int size_a, int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));

        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        }

        __m128i cmp_mask0 = _mm_cmpeq_epi32(v_a, v_b); // pairwise comparison
        __m128i rot1 = _mm_shuffle_epi32(v_b, cyclic_shift1);   // shuffling
        __m128i cmp_mask1 = _mm_cmpeq_epi32(v_a, rot1);
        __m128i rot2 = _mm_shuffle_epi32(v_b, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(v_a, rot2);
        __m128i rot3 = _mm_shuffle_epi32(v_b, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(v_a, rot3);
        __m128i cmp_mask = _mm_or_si128(
                _mm_or_si128(cmp_mask0, cmp_mask1),
                _mm_or_si128(cmp_mask2, cmp_mask3));

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        __m128i p = _mm_shuffle_epi8(v_a, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(set_c + size_c), p);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;    
}

int intersect_simd4x_count(int* set_a, int size_a, int* set_b, int size_b)
{
    int i = 0, j = 0, res = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));

        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        }

        __m128i cmp_mask0 = _mm_cmpeq_epi32(v_a, v_b); // pairwise comparison
        __m128i rot1 = _mm_shuffle_epi32(v_b, cyclic_shift1);   // shuffling
        __m128i cmp_mask1 = _mm_cmpeq_epi32(v_a, rot1);
        __m128i rot2 = _mm_shuffle_epi32(v_b, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(v_a, rot2);
        __m128i rot3 = _mm_shuffle_epi32(v_b, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(v_a, rot3);
        __m128i cmp_mask = _mm_or_si128(
                _mm_or_si128(cmp_mask0, cmp_mask1),
                _mm_or_si128(cmp_mask2, cmp_mask3));

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        res += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            res++;
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return res;
}

int intersect_filter_simd4x(int *set_a, int size_a,
            int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    // for (int x = 0; x < 65536; x += 16)
    //     _mm_prefetch((char*) (byte_check_mask_dict + x), _MM_HINT_T2);
    // for (int x = 0; x < 256; x += 4)
    //     _mm_prefetch((char*) (match_shuffle_dict + x), _MM_HINT_T0);
    // _mm_prefetch((char*) (byte_check_group_a_order), _MM_HINT_T0);
    // _mm_prefetch((char*) (byte_check_group_b_order), _MM_HINT_T0);
    
    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));
        
        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        }
      
        __m128i byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[0]);
        __m128i byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[0]);
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];
        if (__builtin_expect(ms_order == -1, 0)) {
            byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[1]);
            byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[1]);
            byte_check_mask = _mm_and_si128(byte_check_mask,
                    _mm_cmpeq_epi8(byte_group_a, byte_group_b));
            bc_mask = _mm_movemask_epi8(byte_check_mask);
            ms_order = byte_check_mask_dict[bc_mask];
            
            if (__builtin_expect(ms_order == -1, 0)) {
                byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[2]);
                byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[2]);
                byte_check_mask = _mm_and_si128(byte_check_mask,
                        _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                bc_mask = _mm_movemask_epi8(byte_check_mask);
                ms_order = byte_check_mask_dict[bc_mask];
                
                if (__builtin_expect(ms_order == -1, 0)) {
                    byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[3]);
                    byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[3]);
                    byte_check_mask = _mm_and_si128(byte_check_mask,
                            _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                    bc_mask = _mm_movemask_epi8(byte_check_mask);
                    ms_order = byte_check_mask_dict[bc_mask];                    
                }
            }
        }
        if (ms_order == -2) continue; // "no match" in this two block.

        __m128i sf_v_b = _mm_shuffle_epi8(v_b, match_shuffle_dict[ms_order]);
        __m128i cmp_mask = _mm_cmpeq_epi32(v_a, sf_v_b);

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        __m128i p = _mm_shuffle_epi8(v_a, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(set_c + size_c), p);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int intersect_filter_simd4x_count(int* set_a, int size_a,
            int* set_b, int size_b)
{
    int i = 0, j = 0, res = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));
        
        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
        }

        __m128i byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[0]);
        __m128i byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[0]);
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];
        if (__builtin_expect(ms_order == -1, 0)) {
            byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[1]);
            byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[1]);
            byte_check_mask = _mm_and_si128(byte_check_mask,
                    _mm_cmpeq_epi8(byte_group_a, byte_group_b));
            bc_mask = _mm_movemask_epi8(byte_check_mask);
            ms_order = byte_check_mask_dict[bc_mask];
            
            if (__builtin_expect(ms_order == -1, 0)) {
                byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[2]);
                byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[2]);
                byte_check_mask = _mm_and_si128(byte_check_mask,
                        _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                bc_mask = _mm_movemask_epi8(byte_check_mask);
                ms_order = byte_check_mask_dict[bc_mask];
                
                if (__builtin_expect(ms_order == -1, 0)) {
                    byte_group_a = _mm_shuffle_epi8(v_a, byte_check_group_a_order[3]);
                    byte_group_b = _mm_shuffle_epi8(v_b, byte_check_group_b_order[3]);
                    byte_check_mask = _mm_and_si128(byte_check_mask,
                            _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                    bc_mask = _mm_movemask_epi8(byte_check_mask);
                    ms_order = byte_check_mask_dict[bc_mask];                    
                }
            }
        }
        if (ms_order == -2) continue; // "no match" in this two block.

        __m128i sf_v_b = _mm_shuffle_epi8(v_b, match_shuffle_dict[ms_order]);
        __m128i cmp_mask = _mm_cmpeq_epi32(v_a, sf_v_b);

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        res += _mm_popcnt_u32(mask);
    }
    
    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            res++;
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return res;
}


int bp_intersect(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c)
{
    int i = 0, j = 0, size_c = 0;
    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            // bases_c[size_c] = bases_a[i];
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    return size_c;    
}

int bp_intersect_count(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b)
{
    int i = 0, j = 0, res = 0;
    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            res += _mm_popcnt_u32(states_a[i] & states_b[j]);
            i++; j++;
        } else if (bases_a[i] < bases_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return res;    
}

int bp_intersect_scalar2x(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c)
{
    int i = 0, j = 0, size_c = 0;

    int ds_a = size_a - (size_a & 1);
    int ds_b = size_b - (size_b & 1);

    if (ds_a > 0 && ds_b > 0) {
        while (true) {
            int Adat0 = bases_a[i], Adat1 = bases_a[i + 1];
            int Bdat0 = bases_b[j], Bdat1 = bases_b[j + 1];

            if (Adat0 == Bdat0) {
                states_c[size_c] = states_a[i] & states_b[j];
                if (states_c[size_c] != 0) bases_c[size_c++] = Adat0;
            } else if (Adat0 == Bdat1) {
                states_c[size_c] = states_a[i] & states_b[j + 1];
                if (states_c[size_c] != 0) bases_c[size_c++] = Adat0;
                goto advanceB;            
            } else if (Adat1 == Bdat0) {
                states_c[size_c] = states_a[i + 1] & states_b[j];
                if (states_c[size_c] != 0) bases_c[size_c++] = Adat1;          
                goto advanceA;
            }

            if (Adat1 == Bdat1) {
                states_c[size_c] = states_a[i + 1] & states_b[j + 1];
                if (states_c[size_c] != 0) bases_c[size_c++] = Adat1;
                goto advanceAB;
            } else if (Adat1 > Bdat1) {
                goto advanceB;
            }  else {
                goto advanceA;
            }

    advanceA:
            i += 2;            
            if (i >= ds_a) break;
            else continue;
    advanceB:
            j += 2;            
            if (j >= ds_b) break;
            else continue;
    advanceAB:
            i += 2; j += 2;
            if (i >= ds_a || j >= ds_b) break;
            else continue;
        }       
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int bp_intersect_scalar2x_count(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b)
{
    int i = 0, j = 0, res = 0;

    int ds_a = size_a - (size_a & 1);
    int ds_b = size_b - (size_b & 1);

    if (ds_a > 0 && ds_b > 0) {
        while (true) {
            int Adat0 = bases_a[i], Adat1 = bases_a[i + 1];
            int Bdat0 = bases_b[j], Bdat1 = bases_b[j + 1];

            if (Adat0 == Bdat0) {
                res += _mm_popcnt_u32(states_a[i] & states_b[j]);
            } else if (Adat0 == Bdat1) {
                res += _mm_popcnt_u32(states_a[i] & states_b[j + 1]);
                goto advanceB;            
            } else if (Adat1 == Bdat0) {
                res += _mm_popcnt_u32(states_a[i + 1] & states_b[j]);           
                goto advanceA;
            }

            if (Adat1 == Bdat1) {
                res += _mm_popcnt_u32(states_a[i + 1] & states_b[j + 1]);
                goto advanceAB;
            } else if (Adat1 > Bdat1) {
                goto advanceB;
            }  else {
                goto advanceA;
            }

    advanceA:
            i += 2;            
            if (i >= ds_a) break;
            else continue;
    advanceB:
            j += 2;            
            if (j >= ds_b) break;
            else continue;
    advanceAB:
            i += 2; j += 2;
            if (i >= ds_a || j >= ds_b) break;
            else continue;
        }       
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            res += _mm_popcnt_u32(states_a[i] & states_b[j]);
            i++; j++;
        } else if (bases_a[i] < bases_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return res;
}

int bp_intersect_simd4x(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);
    
    while (i < qs_a && j < qs_b) {
        __m128i base_a = _mm_lddqu_si128((__m128i*)(bases_a + i));
        __m128i base_b = _mm_lddqu_si128((__m128i*)(bases_b + j));
        __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
        __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));
        
        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
            // if (a_max < bases_b[j]) continue;
        } else {
            j += 4;
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
            // if (b_max < bases_a[i]) continue;
        }        

        // shift0:
        __m128i cmp_mask0 = _mm_cmpeq_epi32(base_a, base_b);
        __m128i state_c0 = _mm_and_si128(
                _mm_and_si128(state_a, state_b), cmp_mask0);

        // shift1:
        __m128i base_sf1 = _mm_shuffle_epi32(base_b, cyclic_shift1);
        __m128i state_sf1 = _mm_shuffle_epi32(state_b, cyclic_shift1);
        __m128i cmp_mask1 = _mm_cmpeq_epi32(base_a, base_sf1);
        __m128i state_c1 = _mm_and_si128(
                _mm_and_si128(state_a, state_sf1), cmp_mask1);

        // shift2:
        __m128i base_sf2 = _mm_shuffle_epi32(base_b, cyclic_shift2);
        __m128i state_sf2 = _mm_shuffle_epi32(state_b, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(base_a, base_sf2);
        __m128i state_c2 = _mm_and_si128(
                _mm_and_si128(state_a, state_sf2), cmp_mask2);

        // shift3:
        __m128i base_sf3 = _mm_shuffle_epi32(base_b, cyclic_shift3);
        __m128i state_sf3 = _mm_shuffle_epi32(state_b, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(base_a, base_sf3);
        __m128i state_c3 = _mm_and_si128(
                _mm_and_si128(state_a, state_sf3), cmp_mask3);

        __m128i state_all = _mm_or_si128(
                _mm_or_si128(state_c0, state_c1),
                _mm_or_si128(state_c2, state_c3)
                );
        __m128i cmp_mask = _mm_or_si128(
                _mm_or_si128(cmp_mask0, cmp_mask1),
                _mm_or_si128(cmp_mask2, cmp_mask3)
                );
        __m128i state_mask = _mm_cmpeq_epi32(state_all, all_zero_si128);
        int mask = (_mm_movemask_ps((__m128)cmp_mask) &
                ~(_mm_movemask_ps((__m128)state_mask)));

        __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
        __m128i res_s = _mm_shuffle_epi8(state_all, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
        _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            bases_c[size_c] = bases_a[i];
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    return size_c;
}

int bp_intersect_simd4x_count(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b)
{
    int i = 0, j = 0, res = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);
    uint64_t bits[2] __attribute__((aligned(16)));
    
    while (i < qs_a && j < qs_b) {
        __m128i base_a = _mm_lddqu_si128((__m128i*)(bases_a + i));
        __m128i base_b = _mm_lddqu_si128((__m128i*)(bases_b + j));
        __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
        __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));

        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        }

        __m128i and_all = all_zero_si128;

        // shift0:
        __m128i cmp_mask0 = _mm_cmpeq_epi32(base_a, base_b);
        __m128i and_0 = _mm_and_si128(state_a, state_b);
        and_all = _mm_or_si128(and_all, _mm_and_si128(and_0, cmp_mask0));

        // shift1:
        __m128i base_sf1 = _mm_shuffle_epi32(base_b, cyclic_shift1);
        __m128i state_sf1 = _mm_shuffle_epi32(state_b, cyclic_shift1);
        __m128i cmp_mask1 = _mm_cmpeq_epi32(base_a, base_sf1);
        __m128i and_1 = _mm_and_si128(state_a, state_sf1);
        and_all = _mm_or_si128(and_all, _mm_and_si128(and_1, cmp_mask1));

        // shift2:
        __m128i base_sf2 = _mm_shuffle_epi32(base_b, cyclic_shift2);
        __m128i state_sf2 = _mm_shuffle_epi32(state_b, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(base_a, base_sf2);
        __m128i and_2 = _mm_and_si128(state_a, state_sf2);
        and_all = _mm_or_si128(and_all, _mm_and_si128(and_2, cmp_mask2));

        // shift3:
        __m128i base_sf3 = _mm_shuffle_epi32(base_b, cyclic_shift3);
        __m128i state_sf3 = _mm_shuffle_epi32(state_b, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(base_a, base_sf3);
        __m128i and_3 = _mm_and_si128(state_a, state_sf3);
        and_all = _mm_or_si128(and_all, _mm_and_si128(and_3, cmp_mask3));

        // popcnt:
        _mm_store_si128((__m128i*)bits, and_all);
        res += _mm_popcnt_u64(bits[0]);
        res += _mm_popcnt_u64(bits[1]);
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            res += _mm_popcnt_u32(states_a[i] & states_b[j]);
            i++; j++;
        } else if (bases_a[i] < bases_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return res;
}

unsigned long long inter_cnt = 0, no_match_cnt = 0, cmp_cnt = 0;
unsigned long long multimatch_cnt = 0, skew_cnt = 0, low_select_cnt = 0;
unsigned long long byte_check_cnt[4] = {0, 0, 0, 0};

int bp_intersect_filter_simd4x(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b,
            int *bases_c, PackState* states_c)
{
    inter_cnt++;
    int len_a = std::min(size_a, size_b), len_b = std::max(size_a, size_b);
    if (len_a * 32 < len_b) skew_cnt++;
    
    int i = 0, j = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);
    
    while (i < qs_a && j < qs_b) {
        cmp_cnt++;

        __m128i base_a = _mm_lddqu_si128((__m128i*)(bases_a + i));
        __m128i base_b = _mm_lddqu_si128((__m128i*)(bases_b + j));
        __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
        __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));

        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        }
       
        int bn = 0;
        __m128i byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
        __m128i byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];
        if (__builtin_expect(ms_order == -1, 0)) {
            multimatch_cnt++;
            byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[1]);
            byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[1]);
            byte_check_mask = _mm_and_si128(byte_check_mask,
                    _mm_cmpeq_epi8(byte_group_a, byte_group_b));
            bc_mask = _mm_movemask_epi8(byte_check_mask);
            ms_order = byte_check_mask_dict[bc_mask];
            bn++;
            if (__builtin_expect(ms_order == -1, 0)) {
                byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[2]);
                byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[2]);
                byte_check_mask = _mm_and_si128(byte_check_mask,
                        _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                bc_mask = _mm_movemask_epi8(byte_check_mask);
                ms_order = byte_check_mask_dict[bc_mask];
                bn++;
                if (__builtin_expect(ms_order == -1, 0)) {
                    byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[3]);
                    byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[3]);
                    byte_check_mask = _mm_and_si128(byte_check_mask,
                            _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                    bc_mask = _mm_movemask_epi8(byte_check_mask);
                    ms_order = byte_check_mask_dict[bc_mask];
                    bn++;
                }
            }
        }

        // __m128i byte_group_a0 = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
        // __m128i byte_group_b0 = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
        // __m128i byte_group_a1 = _mm_shuffle_epi8(base_a, byte_check_group_a_order[1]);
        // __m128i byte_group_b1 = _mm_shuffle_epi8(base_b, byte_check_group_b_order[1]);
        // __m128i byte_check_mask = _mm_and_si128(
        //         _mm_cmpeq_epi8(byte_group_a0, byte_group_b0),
        //         _mm_cmpeq_epi8(byte_group_a1, byte_group_b1));
        // int bc_mask = _mm_movemask_epi8(byte_check_mask);
        // int ms_order = byte_check_mask_dict[bc_mask];
        // bn += 1;
        // if (__builtin_expect(ms_order == -1, 0)) {
        //     byte_group_a0 = _mm_shuffle_epi8(base_a, byte_check_group_a_order[2]);
        //     byte_group_b0 = _mm_shuffle_epi8(base_b, byte_check_group_b_order[2]);
        //     byte_group_a1 = _mm_shuffle_epi8(base_a, byte_check_group_a_order[3]);
        //     byte_group_b1 = _mm_shuffle_epi8(base_b, byte_check_group_b_order[3]);
        //     __m128i byte_check_mask_p = _mm_and_si128(
        //             _mm_cmpeq_epi8(byte_group_a0, byte_group_b0),
        //             _mm_cmpeq_epi8(byte_group_a1, byte_group_b1));
        //     byte_check_mask = _mm_and_si128(byte_check_mask, byte_check_mask_p);
        //     bc_mask = _mm_movemask_epi8(byte_check_mask);
        //     ms_order = byte_check_mask_dict[bc_mask];
        //     bn += 2;          
        // }

        // if (__builtin_expect(ms_order == -2, 1)) continue; // "no match" in this two block.
        if (ms_order == -2) {no_match_cnt++; continue;}  // "no match" in this two block.
        byte_check_cnt[bn]++;
        // inter_cnt++;        
        // if (ms_order == -2) {no_match_cnt++; continue;} // "no match" in this two block.
        // byte_check_cnt[bn]++;

        __m128i sf_base_b = _mm_shuffle_epi8(base_b, match_shuffle_dict[ms_order]);
        __m128i sf_state_b = _mm_shuffle_epi8(state_b, match_shuffle_dict[ms_order]);
        __m128i cmp_mask = _mm_cmpeq_epi32(base_a, sf_base_b);
        // __m128i and_state = _mm_and_si128(cmp_mask, _mm_and_si128(state_a, sf_state_b));
        __m128i and_state = _mm_and_si128(state_a, sf_state_b);
        __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
        cmp_mask = _mm_andnot_si128(state_mask, cmp_mask);
        int mask = _mm_movemask_ps((__m128)cmp_mask);

        __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
        __m128i res_s = _mm_shuffle_epi8(and_state, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
        _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            bases_c[size_c] = bases_a[i];
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

    double selectivity = (double) size_c / len_a;
    if ((selectivity < 0.3 || len_a < 8) && len_a * 32 >= len_b) low_select_cnt++;

    return size_c;   
}
int bp_intersect_filter_simd4x_count(int* bases_a, PackState* states_a, int size_a,
            int* bases_b, PackState* states_b, int size_b)
{
    inter_cnt++;
    int len_a = std::min(size_a, size_b), len_b = std::max(size_a, size_b);
    if (len_a * 32 < len_b) skew_cnt++;
    int size_c = 0;

    int i = 0, j = 0, res = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);
    uint64_t bits[2] __attribute__((aligned(16)));
    
    while (i < qs_a && j < qs_b) {

        cmp_cnt++;
        __m128i base_a = _mm_lddqu_si128((__m128i*)(bases_a + i));
        __m128i base_b = _mm_lddqu_si128((__m128i*)(bases_b + j));
        __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
        __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));

        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (bases_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (bases_b + j), _MM_HINT_NTA);
            _mm_prefetch((char*) (states_b + j), _MM_HINT_NTA);
        }
       
        int bn = 0;
        __m128i byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
        __m128i byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];
        if (__builtin_expect(ms_order == -1, 0)) {
            multimatch_cnt++;
            byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[1]);
            byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[1]);
            byte_check_mask = _mm_and_si128(byte_check_mask,
                    _mm_cmpeq_epi8(byte_group_a, byte_group_b));
            bc_mask = _mm_movemask_epi8(byte_check_mask);
            ms_order = byte_check_mask_dict[bc_mask];
            bn++;
            if (__builtin_expect(ms_order == -1, 0)) {
                byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[2]);
                byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[2]);
                byte_check_mask = _mm_and_si128(byte_check_mask,
                        _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                bc_mask = _mm_movemask_epi8(byte_check_mask);
                ms_order = byte_check_mask_dict[bc_mask];
                bn++;
                if (__builtin_expect(ms_order == -1, 0)) {
                    byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[3]);
                    byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[3]);
                    byte_check_mask = _mm_and_si128(byte_check_mask,
                            _mm_cmpeq_epi8(byte_group_a, byte_group_b));
                    bc_mask = _mm_movemask_epi8(byte_check_mask);
                    ms_order = byte_check_mask_dict[bc_mask];
                    bn++;
                }
            }
        }

        // __m128i byte_group_a0 = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
        // __m128i byte_group_b0 = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
        // __m128i byte_group_a1 = _mm_shuffle_epi8(base_a, byte_check_group_a_order[1]);
        // __m128i byte_group_b1 = _mm_shuffle_epi8(base_b, byte_check_group_b_order[1]);
        // __m128i byte_check_mask = _mm_and_si128(
        //         _mm_cmpeq_epi8(byte_group_a0, byte_group_b0),
        //         _mm_cmpeq_epi8(byte_group_a1, byte_group_b1));
        // int bc_mask = _mm_movemask_epi8(byte_check_mask);
        // int ms_order = byte_check_mask_dict[bc_mask];
        // bn += 1;
        // if (__builtin_expect(ms_order == -1, 0)) {
        //     byte_group_a0 = _mm_shuffle_epi8(base_a, byte_check_group_a_order[2]);
        //     byte_group_b0 = _mm_shuffle_epi8(base_b, byte_check_group_b_order[2]);
        //     byte_group_a1 = _mm_shuffle_epi8(base_a, byte_check_group_a_order[3]);
        //     byte_group_b1 = _mm_shuffle_epi8(base_b, byte_check_group_b_order[3]);
        //     __m128i byte_check_mask_p = _mm_and_si128(
        //             _mm_cmpeq_epi8(byte_group_a0, byte_group_b0),
        //             _mm_cmpeq_epi8(byte_group_a1, byte_group_b1));
        //     byte_check_mask = _mm_and_si128(byte_check_mask, byte_check_mask_p);
        //     bc_mask = _mm_movemask_epi8(byte_check_mask);
        //     ms_order = byte_check_mask_dict[bc_mask];
        //     bn += 2;          
        // }

        // if (__builtin_expect(ms_order == -2, 1)) continue; // "no match" in this two block.

                
        if (ms_order == -2) {no_match_cnt++; continue;} // "no match" in this two block.
        byte_check_cnt[bn]++;

        __m128i sf_base_b = _mm_shuffle_epi8(base_b, match_shuffle_dict[ms_order]);
        __m128i sf_state_b = _mm_shuffle_epi8(state_b, match_shuffle_dict[ms_order]);
        __m128i cmp_mask = _mm_cmpeq_epi32(base_a, sf_base_b);
        __m128i and_state = _mm_and_si128(cmp_mask, _mm_and_si128(state_a, sf_state_b));

        __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
        cmp_mask = _mm_andnot_si128(state_mask, cmp_mask);
        int mask = _mm_movemask_ps((__m128)cmp_mask);
        size_c += _mm_popcnt_u32(mask);
        
        // popcnt:
        _mm_store_si128((__m128i*)bits, and_state);
        res += _mm_popcnt_u64(bits[0]);
        res += _mm_popcnt_u64(bits[1]);
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            res += _mm_popcnt_u32(states_a[i] & states_b[j]);
            i++; j++;
        } else if (bases_a[i] < bases_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    double selectivity = (double) size_c / len_a;
    if ((selectivity < 0.3 || len_a < 8) && len_a * 32 >= len_b) low_select_cnt++;

    return res;
}

int merge(int *set_a, int size_a, int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    while (i < size_a && j < size_b) {
        if (set_a[i] < set_b[j]) {
            set_c[size_c++] = set_a[i++];
        } else {
            set_c[size_c++] = set_b[j++];
        }
    }
    memcpy(set_c + size_c, set_a + i, (size_a - i) * sizeof(int));
    size_c += (size_a - i);
    memcpy(set_c + size_c, set_b + j, (size_b - j) * sizeof(int));
    size_c += (size_b - j); 

    return size_c;
}

int bp_merge_one(int* bases_a, PackState* states_a, int size_a,
            int v_base, PackState v_bit)
{
    int new_size_a = size_a;
    int i = 0;
    while (i < size_a && bases_a[i] < v_base) ++i;
    if (i == size_a) {
        bases_a[i] = v_base;
        states_a[i] = v_bit;
        new_size_a++;
    } else if (bases_a[i] == v_base) {
        states_a[i] |= v_bit;
    } else {
        // for (int j = size_a; j > i; --j) {
        //     bases_a[j] = bases_a[j - 1];
        //     states_a[j] = states_a[j - 1];
        // }
        memmove(bases_a + i + 1, bases_a + i, (size_a - i) * sizeof(int));
        memmove(states_a + i + 1, states_a + i, (size_a - i) * sizeof(PackState));
        bases_a[i] = v_base;
        states_a[i] = v_bit;
        new_size_a++;
    }

    return new_size_a;   
}

int bp_subtract_visited_simd4x(int* bases_a, PackState* states_a, int size_a,
            PackState* visited, int* bases_c, PackState* states_c)
{
    int i = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    while (i < qs_a) {
        __m128i base_a = _mm_lddqu_si128((__m128i*)(bases_a + i));
        __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
        __m128i state_b = _mm_set_epi32(
            visited[bases_a[i + 3]], visited[bases_a[i + 2]],
            visited[bases_a[i + 1]], visited[bases_a[i]]
            );

        i += 4;
        
        __m128i state_c = _mm_andnot_si128(state_b, state_a);
        __m128i state_mask = _mm_cmpeq_epi32(state_c, all_zero_si128);
        int mask = (15 & ~(_mm_movemask_ps((__m128)state_mask)));
        if (mask != 0) {
            __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
            __m128i res_s = _mm_shuffle_epi8(state_c, shuffle_mask[mask]);
            _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
            _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);
            size_c += _mm_popcnt_u32(mask);            
        }
    }

    while (i < size_a) {
        states_c[size_c] = (states_a[i] & ~(visited[bases_a[i]]));
        if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
        i++;
    }

    return size_c;    
}

int bp_subtract_unvisited_simd4x(int* bases_a, PackState* states_a, int size_a,
            PackState* visited, int* bases_c, PackState* states_c)
{
    int i = 0, size_c = 0;
    int qs_a = size_a - (size_a & 3);
    while (i < qs_a) {
        __m128i base_a = _mm_lddqu_si128((__m128i*)(bases_a + i));
        __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
        __m128i state_b = _mm_set_epi32(
            visited[bases_a[i + 3]], visited[bases_a[i + 2]],
            visited[bases_a[i + 1]], visited[bases_a[i]]
            );

        i += 4;
                
        __m128i state_c = _mm_and_si128(state_b, state_a);
        __m128i state_mask = _mm_cmpeq_epi32(state_c, all_zero_si128);
        int mask = (15 & ~(_mm_movemask_ps((__m128)state_mask)));
        if (mask != 0) {
            __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
            __m128i res_s = _mm_shuffle_epi8(state_c, shuffle_mask[mask]);
            _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
            _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);
            size_c += _mm_popcnt_u32(mask);            
        }
    }

    while (i < size_a) {
        states_c[size_c] = (states_a[i] & visited[bases_a[i]]);
        if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
        i++;
    }

    return size_c;    
}

int bp_subtract_visited(int* bases_a, PackState* states_a, int size_a,
            PackState* visited, int* bases_c, PackState* states_c)
{
    int size_c = 0;
    for (int i = 0; i < size_a; ++i) {
        int u_base = bases_a[i];
        PackState u_state = states_a[i];
        u_state &= (~visited[u_base]);
        if (u_state != 0) {
            bases_c[size_c] = u_base;
            states_c[size_c] = u_state;
            size_c++;
        }
    }

    return size_c;
}

int bp_subtract_unvisited(int* bases_a, PackState* states_a, int size_a,
            PackState* visited, int* bases_c, PackState* states_c)
{
    int size_c = 0;
    for (int i = 0; i < size_a; ++i) {
        int u_base = bases_a[i];
        PackState u_state = states_a[i];
        u_state &= visited[u_base];
        if (u_state != 0) {
            bases_c[size_c] = u_base;
            states_c[size_c] = u_state;
            size_c++;
        }
    }

    return size_c;
}

int subtract(int *set_a, int size_a, int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    while (i < size_a) {
        if (j == size_b) {
            memmove(set_c + size_c, set_a + i, (size_a - i) * sizeof(int));
            return size_c + (size_a - i);
        }
        if (set_a[i] < set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++;
        } else if (set_a[i] > set_b[j]) {
            j++;
        } else {
            i++; j++;
        }
    }

    return size_c;
}

int subtract(int* base_a, PackState* state_a, int size_a,
            int* base_b, PackState* state_b, int size_b,
            int *base_c, PackState* state_c)
{
    int i = 0, j = 0, size_c = 0;
    while (i < size_a) {
        if (j == size_b) {
            memmove(base_c + size_c, base_a + i, (size_a - i) * sizeof(int));
            memmove(state_c + size_c, state_a + i, (size_a - i) * sizeof(PackState));
            return size_c + (size_a - i);
        }
        if (base_a[i] < base_b[j]) {
            base_c[size_c] = base_a[i];
            state_c[size_c] = state_a[i];
            size_c++; i++;
        } else if (base_a[i] > base_b[j]) {
            j++;
        } else {
            state_c[size_c] = (state_a[i] & (~state_b[j]));
            if (state_c[size_c] != 0) base_c[size_c++] = base_a[i];
            i++; j++;
        }
    }

    return size_c;
}