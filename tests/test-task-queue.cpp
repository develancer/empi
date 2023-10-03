/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstring>
#include <thread>
#include <vector>
#include "TaskQueue.h"
#include "Testing.h"

const int THREAD_COUNT = 4;
const int ITEM_COUNT = 8;

int processed[ITEM_COUNT];

void consumer_main(TaskQueue<int>* task_queue) {
    int number;
    while (task_queue->get(number)) {
        processed[number] = true;
        task_queue->notify();
    }
}

void test_put_twice() {
    TaskQueue<int> task_queue;
    memset(processed, 0, sizeof processed);

    std::vector<std::thread> threads;
    for (int i=0; i<THREAD_COUNT; ++i) {
        threads.emplace_back(&consumer_main, &task_queue);
    }

    task_queue.put(0);
    task_queue.wait_for_tasks();
    for (int i=0; i<ITEM_COUNT; ++i) {
        ASSERT_EQUALS(i == 0, processed[i]);
    }
    task_queue.put(1);
    task_queue.wait_for_tasks();
    for (int i=0; i<ITEM_COUNT; ++i) {
        ASSERT_EQUALS(i <= 1, processed[i]);
    }

    task_queue.terminate();

    for (auto &thread : threads) {
        thread.join();
    }
}

void test_put_list() {
    TaskQueue<int> task_queue;
    memset(processed, 0, sizeof processed);

    std::vector<std::thread> threads;
    for (int i=0; i<THREAD_COUNT; ++i) {
        threads.emplace_back(&consumer_main, &task_queue);
    }

    std::list<int> numbers;
    for (int i=0; i<ITEM_COUNT/2; ++i) {
        numbers.push_back(i);
    }
    task_queue.put(std::move(numbers));
    task_queue.wait_for_tasks();
    for (int i=0; i<ITEM_COUNT; ++i) {
        ASSERT_EQUALS(i < ITEM_COUNT/2, processed[i]);
    }

    numbers.clear();
    for (int i=ITEM_COUNT/2; i<ITEM_COUNT; ++i) {
        numbers.push_back(i);
    }
    task_queue.put(std::move(numbers));
    task_queue.wait_for_tasks();
    for (int i=0; i<ITEM_COUNT; ++i) {
        ASSERT_EQUALS(1, processed[i]);
    }

    task_queue.terminate();

    for (auto &thread : threads) {
        thread.join();
    }
}

int main() {
    test_put_twice();
    test_put_list();
    puts("OK");
}
