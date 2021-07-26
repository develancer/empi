#include <cassert>
#include <cstdio>
#include <utility>

class Item
{
public:
    static int copy_called;
    static int copy_ctor_called;
    static int ctor_called;
    static int move_ctor_called;

    static void reset()
    {
        copy_called = copy_ctor_called = ctor_called = move_ctor_called = 0;
    }

    Item() {
        ++ctor_called;
    }

    Item(const Item&) {
        ++copy_ctor_called;
    }

    Item(Item&&) {
        ++move_ctor_called;
    }

    void operator=(const Item&) {
        ++copy_called;
    }
};

int Item::copy_called = 0;
int Item::copy_ctor_called = 0;
int Item::ctor_called = 0;
int Item::move_ctor_called = 0;

class Internal
{
    Item item;

public:
    Internal(Item item) : item(std::move(item))
    { }
};

class External
{
    Internal internal;

public:
    External(Item item) : internal(std::move(item))
    { }
};

int main(void)
{
    Item item;
    Item::reset();

    External external(std::move(item));
    assert(Item::copy_called == 0);
    assert(Item::copy_ctor_called == 0);
    assert(Item::ctor_called == 0);
    assert(Item::move_ctor_called == 3);

    puts("OK");
}
