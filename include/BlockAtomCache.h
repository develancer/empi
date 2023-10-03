/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_ATOM_CACHE_H
#define EMPI_BLOCK_ATOM_CACHE_H

#include <map>
#include <memory>
#include <shared_mutex>
#include <optional>
#include "IndexRange.h"

// forward declaration to avoid cyclic dependency between BlockAtomCache.h and BlockAtom.h
class BlockExtendedAtom;

/**
 * Internal data object used in BlockAtomCache.
 */
struct BlockAtomCacheItem {
    IndexRange index_range;
    std::shared_ptr<BlockExtendedAtom> atom;
};

/**
 * Cache for BlockExtendedAtom instances used internally in Block class.
 */
//using BlockAtomCache = std::map<size_t, BlockAtomCacheItem>;
class BlockAtomCache
{
    std::map<size_t, BlockAtomCacheItem> items;
    std::shared_mutex items_mutex;

public:
    [[nodiscard]] std::shared_ptr<BlockExtendedAtom> get(size_t key);

    void set(size_t key, IndexRange index_range, std::shared_ptr<BlockExtendedAtom> atom);

    void remove_overlapping(IndexRange range_to_overlap);
};

/**
 * Reference to one of the items in the cache, identified by a unique key.
 */
class BlockAtomCacheSlot {
    std::weak_ptr<BlockAtomCache> cache;
    size_t key;

public:
    BlockAtomCacheSlot(const std::shared_ptr<BlockAtomCache>& cache, size_t key);

    [[nodiscard]] std::shared_ptr<BlockExtendedAtom> get() const;

    void set(IndexRange index_range, std::shared_ptr<BlockExtendedAtom> atom);
};

#endif //EMPI_BLOCK_ATOM_CACHE_H
