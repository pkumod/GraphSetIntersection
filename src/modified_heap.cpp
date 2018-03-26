#include "modified_heap.hpp"

ModifiedHeap::ModifiedHeap(int max_size)
{
    align_malloc((void**)&heap, 32, sizeof(MHNode) * (max_size + 1));
    align_malloc((void**)&up_info, 32, sizeof(UpdateInfo) * (max_size + 1));
    align_malloc((void**)&up_list, 32, sizeof(int) * (max_size + 1));
    size = max_size;
    reset_label = 0; up_list_cnt = 0;
    for (int i = 0; i < size; ++i) {
        heap[i + 1] = MHNode(i, 0, 0.0);
        up_info[i] = UpdateInfo(i + 1, 0, 0.0);
    }
}

ModifiedHeap::~ModifiedHeap()
{
    free(heap);
    free(up_info);
    free(up_list);
}

const int RESET_LABEL_MASK = 0x7fffffff;
const int UPDATE_LABEL_MASK = 0x80000000;

void ModifiedHeap::inc(int key, double delta)
{
    UpdateInfo& info = up_info[key];
    // assert(info.pos > 0);
    if ((info.state & RESET_LABEL_MASK) != reset_label) {
        info.state = reset_label;
        info.utd_val = delta;
    } else
        info.utd_val += delta;

    if ((info.state & UPDATE_LABEL_MASK) == 0) {
        up_list[up_list_cnt++] = key;
        info.state |= UPDATE_LABEL_MASK;
    }
}

inline bool ModifiedHeap::is_zero(int pos)
{
    return heap[pos].label != reset_label;
}

int ModifiedHeap::get_size()
{
    return size;
}

bool ModifiedHeap::in_heap(int key)
{
    return up_info[key].pos > 0;
}

MHNode ModifiedHeap::top()
{
    adjust();
    return heap[1];
}

int ModifiedHeap::pop()
{
    adjust();
    int key = heap[1].key;    
    swap_mhnode(1, size);
    size--;
    up_info[key].pos = 0;
    down(1);
    
    return key;
}

void ModifiedHeap::del(int key)
{
    adjust();
    int pos = up_info[key].pos;
    if (pos == -1) return;
    swap_mhnode(pos, size);
    size--;
    up_info[key].pos = 0;
    if (is_zero(pos)) down(pos);
    else up(pos);
    
}

void ModifiedHeap::reset()
{
    reset_label++;
    up_list_cnt = 0;
}

void ModifiedHeap::adjust()
{
    for (int i = 0; i < up_list_cnt; ++i) {
        int key = up_list[i];
        UpdateInfo& info = up_info[key];
        heap[info.pos].val = info.utd_val;
        heap[info.pos].label = reset_label;
        up(info.pos);
        info.state -= UPDATE_LABEL_MASK;        
    }
    up_list_cnt = 0;
}

void ModifiedHeap::up(int pos)
{
    // printf("up: %d %d\n", heap[pos].key, pos);
    // if (is_zero(pos)) return;
    while (pos >= 2 &&
        (is_zero(pos >> 1) || heap[pos].val > heap[pos >> 1].val)) {
        // printf("swap: %d %d <-> %d %d\n", heap[pos].key, pos, heap[pos >> 1].key, pos >> 1);
        swap_mhnode(pos, pos >> 1);
        pos >>= 1;
    }
}

void ModifiedHeap::down(int pos)
{
    // printf("down: %d %d\n", heap[pos].key, pos);
    // for (int i = 1; i < 10; ++i) printf("%d %d %d %.2f\n", i, heap[i].key, heap[i].label, heap[i].val);
    if (is_zero(pos)) heap[pos].val = 0.0;
    int j = (pos << 1);
    while (j <= size) {
        if (j < size) {
            char zero_state = is_zero(j) | (is_zero(j + 1) << 1);
            if (zero_state == 0) j += (heap[j + 1].val > heap[j].val);
            else if (zero_state == 1) j++;
        }

        if (is_zero(j) || heap[pos].val >= heap[j].val) break;
        // printf("swap: %d %d <-> %d %d\n", heap[pos].key, pos, heap[j].key, j);
        swap_mhnode(pos, j);
        pos = j; j = pos << 1;
    }
}

inline void ModifiedHeap::swap_mhnode(int pos_x, int pos_y)
{
    std::swap(heap[pos_x], heap[pos_y]);
    up_info[heap[pos_x].key].pos = pos_x;
    up_info[heap[pos_y].key].pos = pos_y;
}
