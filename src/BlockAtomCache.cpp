/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <mutex>
#include "BlockAtomCache.h"

//////////////////////////////////////////////////////////////////////////////

std::shared_ptr<BlockExtendedAtom> BlockAtomCache::get(size_t key) {
    std::shared_lock items_lock_guard(items_mutex);
    auto iterator = items.find(key);
    if (iterator == items.end()) {
        return nullptr;
    }
    return iterator->second.atom;
}

void BlockAtomCache::set(size_t key, IndexRange index_range, std::shared_ptr<BlockExtendedAtom> atom) {
    std::unique_lock items_lock_guard(items_mutex);
    items.insert_or_assign(key, BlockAtomCacheItem{index_range, std::move(atom)});
}

void BlockAtomCache::remove_overlapping(IndexRange range_to_overlap) {
    std::unique_lock items_lock_guard(items_mutex);
    for (auto iterator = items.begin(); iterator != items.end();) {
        if (!iterator->second.index_range.overlap(range_to_overlap)) {
            ++iterator;
        } else {
            iterator = items.erase(iterator);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

BlockAtomCacheSlot::BlockAtomCacheSlot(const std::shared_ptr<BlockAtomCache>& cache, size_t key) : cache(cache), key(key) {}

std::shared_ptr<BlockExtendedAtom> BlockAtomCacheSlot::get() const {
    return cache.lock()->get(key);
}

void BlockAtomCacheSlot::set(IndexRange index_range, std::shared_ptr<BlockExtendedAtom> atom) {
    cache.lock()->set(key, index_range, std::move(atom));
}
