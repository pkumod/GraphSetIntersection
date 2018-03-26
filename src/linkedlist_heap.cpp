#include "linkedlist_heap.hpp"

LinkedListHeap::LinkedListHeap(int _capacity)
{
    capacity = _capacity;
    align_malloc((void**)&linkedlist, 32, sizeof(LHNode) * capacity);
    align_malloc((void**)&headers, 32, sizeof(LHHeader) * capacity * 2); // directed graph.
    align_malloc((void**)&up_list, 32, sizeof(LHUpdateInfo) * capacity);

    size = capacity;
    for (int i = 0; i < size; ++i) {
        linkedlist[i] = LHNode(0, 0, i - 1, i + 1);
        headers[i] = LHHeader(-1, -1, 0);
    }
    linkedlist[size - 1].next = -1;
    headers[0].first = 0; headers[0].last = size - 1;
    head = 0; tail = size - 1;
    reset_label = 0;
    up_list_cnt = 0;
}

LinkedListHeap::~LinkedListHeap()
{
    free(linkedlist);    
    free(headers);
    free(up_list);
}
void LinkedListHeap::inc(int key)
{
    LHNode& cur_node = linkedlist[key];
    
    // assert(cur_node.val != -1);

    if ((cur_node.state & RESET_LABEL_MASK) != reset_label) {
        cur_node.val = 0;
        cur_node.state = reset_label;
    }
    cur_node.val++;

    if ((cur_node.state & UPDATE_LABEL_MASK) == 0) {
        up_list[up_list_cnt++] = LHUpdateInfo(key, cur_node.val - 1);        
        cur_node.state |= UPDATE_LABEL_MASK;
    }
}

int LinkedListHeap::top()
{
    adjust();
    return head;
}

int LinkedListHeap::pop()
{
    adjust();
    
    LHNode& cur_node = linkedlist[head];    
    int val = cur_node.val;
    if (cur_node.state != reset_label) val = 0;
    int key = head;

    cur_node.val = -1;
    head = cur_node.next;
    if (cur_node.next == -1) tail = -1;
    else linkedlist[cur_node.next].prev = -1;

    LHHeader& top_header = headers[val];
    if (top_header.first == top_header.last) top_header.first = top_header.last = -1;
    else top_header.first = head;

    size--;

    return key;
}

void LinkedListHeap::del(int key)
{
    adjust();

    LHNode& cur_node = linkedlist[key];

    if (cur_node.val == -1) return;

    int val = cur_node.val;
    if (cur_node.state != reset_label) val = 0;

    cur_node.val = -1;
    if (cur_node.prev == -1) head = cur_node.next;
    else linkedlist[cur_node.prev].next = cur_node.next;
    if (cur_node.next == -1) tail = cur_node.prev;
    else linkedlist[cur_node.next].prev = cur_node.prev;

    LHHeader& cur_header = headers[val];
    if (cur_header.first == cur_header.last) cur_header.first = cur_header.last = -1;
    else if (cur_header.first == key) cur_header.first = cur_node.next;
    else if (cur_header.last == key) cur_header.last = cur_node.prev;

    size--;
}

void LinkedListHeap::reset()
{
    reset_label++;
    up_list_cnt = 0;
    headers[0] = LHHeader(head, tail, reset_label);
}

void LinkedListHeap::adjust()
{
    for (int i = 0; i < up_list_cnt; ++i) {        
        up(up_list[i].key, up_list[i].old_val);
    }
    up_list_cnt = 0;
}

bool LinkedListHeap::in_heap(int key)
{
    return linkedlist[key].val != -1;
}

bool LinkedListHeap::is_top_zero()
{
    adjust();
    if (linkedlist[head].state != reset_label || linkedlist[head].val == 0)
        return true;
    return false;
}

int LinkedListHeap::get_top_val()
{
    adjust();
    if (linkedlist[head].state != reset_label) return 0;
    else return linkedlist[head].val;
}

int LinkedListHeap::get_size()
{
    return size;
}

void LinkedListHeap::up(int key, int old_val)
{
    LHNode& cur_node = linkedlist[key];

    // update old header.
    LHHeader& cur_header = headers[old_val];
    if (cur_header.first == cur_header.last) cur_header.first = cur_header.last = -1;
    else if (cur_header.first == key) cur_header.first = cur_node.next;
    else if (cur_header.last == key) cur_header.last = cur_node.prev;

    // split cur_node from linkedlist.
    if (cur_node.prev == -1) head = cur_node.next;
    else linkedlist[cur_node.prev].next = cur_node.next;
    if (cur_node.next == -1) tail = cur_node.prev;
    else linkedlist[cur_node.next].prev = cur_node.prev;

    // find the new header.
    LHHeader& new_header = headers[cur_node.val];
    if (new_header.label != reset_label || new_header.first == -1) { // empty header.
        new_header.first = new_header.last = key;
        new_header.label = reset_label;

        int right_val = cur_node.val - 1;
        while (right_val > old_val && (headers[right_val].label != reset_label ||
                headers[right_val].first == -1)) right_val--;

        const LHHeader& right_header = headers[right_val];
        if (right_header.first == -1) {
            if (cur_node.prev == -1) head = key;
            else linkedlist[cur_node.prev].next = key;
            if (cur_node.next == -1) tail = key;
            else linkedlist[cur_node.next].prev = key;
        } else {
            int right_node_pos = right_header.first;
            int left_node_pos = linkedlist[right_node_pos].prev;
            cur_node.prev = left_node_pos;
            cur_node.next = right_node_pos;
            if (left_node_pos == -1) head = key;
            else linkedlist[left_node_pos].next = key;
            linkedlist[right_node_pos].prev = key;
        }
    } else {
        int left_node_pos = new_header.last;
        int right_node_pos = linkedlist[left_node_pos].next;
        cur_node.prev = left_node_pos;
        cur_node.next = right_node_pos;
        linkedlist[left_node_pos].next = key;
        if (right_node_pos == -1) tail = key;
        else linkedlist[right_node_pos].prev = key;
        new_header.last = key;
    }

    cur_node.state -= UPDATE_LABEL_MASK;
}

void LinkedListHeap::print()
{
    adjust();

    if (head == -1) {
        printf("empty list!\n");
        return;
    }

    printf("list:");
    int cur_key = head;
    while (cur_key != -1) {
        int val = linkedlist[cur_key].val;
        if ((linkedlist[cur_key].state & RESET_LABEL_MASK) != reset_label) val = 0;
        printf("(%d,%d) ", cur_key, val);
        cur_key = linkedlist[cur_key].next;
    }
    
    printf("\n");
}

bool LinkedListHeap::check()
{
    bool flag = true;
    int prev_val = INT32_MAX;
    int cnt = 0;
    int cur_key = head;
    while (cur_key != -1 && flag) {
        int val = linkedlist[cur_key].val;
        if ((linkedlist[cur_key].state & RESET_LABEL_MASK) != reset_label) val = 0;
        
        if (++cnt > size) flag = false;
        if (prev_val < val) {flag = false; printf("val=%d prev_val=%d key=%d\n", val, prev_val, cur_key);}
        cur_key = linkedlist[cur_key].next;
        prev_val = val;
    }

    flag = flag && (cnt == size);
    if (cnt != size) printf("cnt=%d size=%d\n", cnt, size);

    if (!flag) printf("check failed!\n");

    return flag;
}

bool LinkedListHeap::check_size()
{
    bool flag = true;
    int cnt = 0;
    int cur_key = head;
    while (cur_key != -1 && flag) {     
        if (++cnt > size) flag = false;        
        cur_key = linkedlist[cur_key].next;
    }

    flag = flag && (cnt == size);
    if (cnt != size) printf("cnt=%d size=%d\n", cnt, size);
    if (!flag) printf("check_size failed!\n");

    return flag;
}