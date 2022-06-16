/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_THREAD_H
#define EMPI_THREAD_H

#include <thread>
#include <utility>

class Thread : public std::thread
{
public:
    template<class Function>
    explicit Thread(Function function) : std::thread(std::move(function))
    {
        affix_to_cpu();
    }

private:
    void affix_to_cpu();
};

#endif //EMPI_THREAD_H
