/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <cstdio>
#include <vector>
#include "heap.hpp"

typedef HeapItem<double> Item;

int main(void) {
	const int REPEATS = 42;
	const int LENGTH = 1000;

	std::vector<Item> data(LENGTH);
	std::vector<double> values(LENGTH);
	for (int i=0; i<LENGTH; ++i) {
		values[i] = rand() / (double) RAND_MAX;
		data[i].key = i;
		data[i].value = values[i];
	}
	Heap<double> heap(values);

	for (int repeat=0; repeat<REPEATS; ++repeat) {
		std::sort(data.begin(), data.end(), [](const Item& a, const Item& b) { return a.value > b.value; });

		for (int i=0; i<LENGTH; ++i) {
			Item item = data[i], peeked = heap.peek();
			if (item.key != peeked.key || item.value != peeked.value) {
				printf("ERROR at repeat %d item %d: [%zu]%lf instead of [%zu]%lf!\n", repeat, i, peeked.key, peeked.value, item.key, item.value);
				return 1;
			}
			heap.update(item.key, -1);
		}
		if (heap.peek().value != -1) {
			printf("ERROR after repeat %d: %lf!\n", repeat, heap.peek().value);
			return 1;
		}

		std::random_shuffle(data.begin(), data.end());
		for (int i=0; i<LENGTH; ++i) {
			data[i].value = rand() / (double) RAND_MAX;
			heap.update(data[i].key, data[i].value);
		}
	}
}
