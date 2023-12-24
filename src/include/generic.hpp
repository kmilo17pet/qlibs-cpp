#ifndef QLIBS_GENERIC
#define QLIBS_GENERIC

#include "include/types.hpp"

namespace qlibs {
    namespace generic {
        using compareFcn_t = int (*)( const void *, const void *, void * );
        using forEachFcn_t = int (* const)( int, void *, void * );

        void swap( void * const x, void * const y, size_t n );
        void sort( void * const pbase, size_t n, size_t size, compareFcn_t cmp, void *arg = nullptr );
        void reverse( void * const pbase, const size_t size, const size_t init, const size_t end );
        void rotate( void * const pbase, const size_t size, const size_t n, const int k );
        void* set( void * const pbase, const size_t size, const size_t n,  const void * const ref );
        void* lSearch( const void *key, const void *pbase, const size_t n, const size_t size, compareFcn_t compar, void *arg = nullptr );
        void* bSearch( const void *key, const void *pbase, const size_t n, const size_t size, compareFcn_t compar, void *arg = nullptr );
        int forEach( void *pbase, const size_t size, const size_t n, forEachFcn_t f, const bool dir = false, void *arg = nullptr );
    }
}

#endif /*QLIBS_GENERIC*/