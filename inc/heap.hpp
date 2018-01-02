/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_HEAP_HPP
#define	EMPI_HEAP_HPP

#include <algorithm>
#include <map>
#include <stdexcept>

template<typename V>
struct HeapItem {
	size_t key;
	V value;
};

template<typename V>
class Heap {
public:
	typedef V value_t;
	typedef HeapItem<V> item_t;

private:
	const size_t size;
	std::vector<value_t> heap;
	std::vector<size_t> heap_from_key;
	std::vector<size_t> key_from_heap;

	void init(void) {
		for (size_t i=0; i<size; ++i) {
			heap_from_key[i] = key_from_heap[i] = i;
		}
		size_t i = size / 2 - 1;
		do {
			_update(i, heap[i], false);
		} while (i-- > 0);
	}

	void _update(size_t i, value_t value, bool allow_swap_with_parent) {
		size_t a, b, p;

		#define SHOULD_SWAP_WITH_PARENT (i > 0 && value > heap[p=(i-1)/2])
		#define SHOULD_SWAP_WITH_CHILD_A (a < size && value < heap[a])
		#define SHOULD_SWAP_WITH_CHILD_B (b < size && value < heap[b] && heap[b] > heap[a])
		#define SHOULD_SWAP_WITH_CHILD (a=2*i+1, b=a+1, SHOULD_SWAP_WITH_CHILD_A || SHOULD_SWAP_WITH_CHILD_B)
		#define SWAP_WITH(n) std::swap(heap[i], heap[n]), std::swap(key_from_heap[i], key_from_heap[n]), std::swap(heap_from_key[key_from_heap[i]], heap_from_key[key_from_heap[n]]), i = n

		if (allow_swap_with_parent && SHOULD_SWAP_WITH_PARENT) {
			do {
				SWAP_WITH(p);
			} while (SHOULD_SWAP_WITH_PARENT);

		} else if (SHOULD_SWAP_WITH_CHILD) {
			do {
				if (SHOULD_SWAP_WITH_CHILD_B) {
					SWAP_WITH(b);
				} else if (SHOULD_SWAP_WITH_CHILD_A) {
					SWAP_WITH(a);
				}
			} while (SHOULD_SWAP_WITH_CHILD);

		}
		#undef SHOULD_SWAP_WITH_PARENT
		#undef SHOULD_SWAP_WITH_CHILD_A
		#undef SHOULD_SWAP_WITH_CHILD_B
		#undef SHOULD_SWAP_WITH_CHILD
		#undef SWAP_WITH
	}

public:
	Heap(const std::vector<value_t>& values) : size(values.size()), heap(values), heap_from_key(size), key_from_heap(size) {
		init();
	}

	Heap(std::vector<value_t>&& values) : size(values.size()), heap(std::move(values)), heap_from_key(size), key_from_heap(size) {
		init();
	}

	item_t peek(void) const {
		return item_t{key_from_heap.front(), heap.front()};
	}

	void update(size_t key, value_t value) {
		size_t i = heap_from_key[key];
		heap[i] = value;
		_update(i, value, true);
	}

	Heap(const Heap&) =delete;
	Heap(Heap&&) =default;
	void operator=(const Heap&) =delete;
};

#endif	/* EMPI_HEAP_HPP */
