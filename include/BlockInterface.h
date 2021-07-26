/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_INTERFACE_H
#define EMPI_BLOCK_INTERFACE_H

/**
 * Simple interface so that Computer instances can notify Block instances about completed computation requests.
 */
class BlockInterface {
public:
    /**
     * This method is called internally by Computer instances
     * to notify that the last computation request has been completed.
     * It is then the block's responsibility to update its internal cache.
     */
    virtual void notify() = 0;
};

#endif //EMPI_BLOCK_INTERFACE_H
