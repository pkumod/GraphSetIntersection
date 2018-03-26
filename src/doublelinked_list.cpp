#include "doublelinked_list.hpp"

DoubleLinkedList::DoubleLinkedList(int _capacity, int _range): 
        capacity(_capacity), range(_range)
{
    align_malloc((void**)&nodes, 32, sizeof(DoubleLinkedNode) * capacity);
    align_malloc((void**)&key2pos, 32, sizeof(int) * capacity);
    memset(key2pos, -1, sizeof(int) * range);

    head = -1; tail = -1; node_idx = 0;
}

DoubleLinkedList::~DoubleLinkedList()
{
    free(nodes);
    free(key2pos);
}

void DoubleLinkedList::add(int key) // add to list's tail.
{
    if (key2pos[key] != -1) return;
    key2pos[key] = node_idx;

    DoubleLinkedNode& cur_node = nodes[node_idx];
    cur_node.key = key;
    cur_node.next = -1;
    cur_node.prev = tail;
    if (tail != -1) nodes[tail].next = node_idx;
    tail = node_idx;
    if (head == -1) head = node_idx;
    
    node_idx++;
}

void DoubleLinkedList::del(int key)
{
    if (key2pos[key] == -1) return;
    DoubleLinkedNode& cur_node = nodes[key2pos[key]];
    if (cur_node.prev == -1) head = cur_node.next;
    else nodes[cur_node.prev].next = cur_node.next;
    if (cur_node.next == -1) tail = cur_node.prev;
    else nodes[cur_node.next].prev = cur_node.prev;
    key2pos[key] = -1;
}

int DoubleLinkedList::get_head()
{
    if (head != -1) return nodes[head].key;
    else return -1;
}

int DoubleLinkedList::get_tail()
{
    if (tail != -1) return nodes[tail].key;
    else return -1;
}

int DoubleLinkedList::pop_head()
{
    if (head == -1) return -1;
    DoubleLinkedNode& cur_node = nodes[head];
    head = cur_node.next;
    if (cur_node.next == -1) tail = -1;
    else nodes[cur_node.next].prev = -1;   
    key2pos[cur_node.key] = -1;

    return cur_node.key;
}

int DoubleLinkedList::pop_tail()
{
    if (tail == -1) return -1;
    DoubleLinkedNode& cur_node = nodes[tail];
    tail = cur_node.prev;
    if (cur_node.prev == -1) head = -1;
    else nodes[cur_node.prev].next = -1;
    key2pos[cur_node.key] = -1;

    return cur_node.key;
}

void DoubleLinkedList::print()
{
    if (head == -1) {
        printf("empty list!\n");
        return;
    }

    int cur_idx = head;
    printf("list:");
    while (cur_idx != -1) {
        printf(" %d", nodes[cur_idx].key);
        cur_idx = nodes[cur_idx].next;
    }
    printf("\n");
}

void DoubleLinkedList::print_inverse()
{
    if (tail == -1) {
        printf("empty list!\n");
        return;        
    }

    int cur_idx = tail;
    printf("inversed_list:");
    while (cur_idx != -1) {
        printf(" %d", nodes[cur_idx].key);
        cur_idx = nodes[cur_idx].prev;
    }
    printf("\n");
}