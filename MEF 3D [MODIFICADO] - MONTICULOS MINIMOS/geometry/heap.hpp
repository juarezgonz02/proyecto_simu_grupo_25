#include <iostream>
#include <cmath>
#include "node.hpp"
using namespace std;

class MinHeap{
    Node **heap_array; // pointer to array of elements in heap
    int capacity; // maximum possible size of min heap
    int heap_size; // Current number of elements in min heap

public:
    // Constructor: Initialise a capacity and heap_array;
    MinHeap(int capacity){
        this->heap_size = 0;
        this->capacity = capacity;
        this->heap_array = (Node **)malloc(sizeof(Node *) * capacity);
    }

    ~MinHeap(){
        free(heap_array);
    }

    // method to heapify a subtree with the root at given index i
    void MinHeapify(int i){
        /* A recursive method to heapify 'heap_array' */
        int l = left(i);
        int r = right(i);

        int smallest = i;
        if (l < heap_size && heap_array[l]->get_ID() < heap_array[i]->get_ID())
            smallest = l;
        if (r < heap_size && heap_array[r]->get_ID() < heap_array[smallest]->get_ID())
            smallest = r;

        if (smallest != i){
            swap(heap_array[i], heap_array[smallest]);
            MinHeapify(smallest);
        }
    }

    // method to get index of parent of node at index i
    int parent(int i){ return floor((i-1)/2); }

    // method to get index of left child of node at index i
    int left(int i){ return (2*i + 1); }

    // method to get index of right child of node at index i
    int right(int i){ return (2*i + 2); }

      // method to inserts a new key 'k'
    void insert(Node* k){

        if (heap_size == capacity){
            cout << "\nOverflow: Could not insertKey\n";
            return;
        }

        // Inserting the new key at the end
        int i = heap_size;
        heap_array[heap_size++] = k;

        while (i != 0 && heap_array[parent(i)]->get_ID() > heap_array[i]->get_ID()){
            swap(heap_array[i], heap_array[parent(i)]);
            i = parent(i);
        }

        cout << heap_array[heap_size-1]->get_ID() << " " << heap_array[heap_size-1]->get_x_coordinate() << "\n";
    }


    Node* getNodeById(int id){

        int i = 0;
        int found = findFrom(i, &id);

        return heap_array[found];

    }

private: 
    int findFrom(int i, int* find){
        
        int found = -1; 

        if(i > heap_size-1){
            return -1;
        }

        if(heap_array[i]->get_ID()-1 == *find){
            found =  i;
        }
        else{

            int L = findFrom(left(i), find);
            
            int R = findFrom(right(i), find);

            return max(L, R);
        }

        return found; 
    }
}

;