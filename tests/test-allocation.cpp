#include <cassert>
#include <random>

#include "Array.h"

struct Item {
    unsigned id;
};

unsigned N = 5;
unsigned M = 10;
unsigned allocator_called = 0;
unsigned deallocator_called = 0;

Item* myAllocator(size_t n)
{
    ++allocator_called;
    static unsigned next_id = 0;
    Item* result = new Item[n];
    for (unsigned i = 0; i<n; ++i) {
        result[i].id = next_id++;
    }
    return result;
}

void myDeallocator(Item* array)
{
    ++deallocator_called;
    delete[] array;
}

void test_alloc()
{
    {
        auto data = Array2D<Item>(N, M, myAllocator, myDeallocator);
        for (unsigned n = 0; n<N; ++n) {
            for (unsigned m = 0; m<M; ++m) {
                assert(data[n][m].id==n*M+m);
            }
        }
    }
    assert(allocator_called==N);
    assert(deallocator_called==N);
}

int main(void)
{
    test_alloc();
    puts("OK");
}