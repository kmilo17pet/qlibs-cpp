#ifndef QLIBS_GENERIC
#define QLIBS_GENERIC

#include "include/types.hpp"

namespace qlibs {
    namespace generic {
        using compareFcn_t = int (*)( const void *, const void *, void * );
        using forEachFcn_t = int (* const)( int, void *, void * );

        void swap( void * const x,
                   void * const y,
                   size_t n ) noexcept;
        void sort( void * const pbase,
                   size_t n,
                   size_t size,
                   compareFcn_t cmp,
                   void *arg = nullptr ) noexcept;
        void reverse( void * const pbase,
                      const size_t size,
                      const size_t init,
                      const size_t end ) noexcept;
        void rotate( void * const pbase,
                     const size_t size,
                     const size_t n,
                     const int k ) noexcept;
        void* set( void * const pbase,
                   const size_t size,
                   const size_t n,
                   const void * const ref ) noexcept;
        void* lSearch( const void *key,
                       const void *pbase,
                       const size_t n,
                       const size_t size,
                       compareFcn_t compar,
                       void *arg = nullptr ) noexcept;
        void* bSearch( const void *key,
                       const void *pbase,
                       const size_t n,
                       const size_t size,
                       compareFcn_t compar,
                       void *arg = nullptr ) noexcept;
        int forEach( void *pbase,
                     const size_t size,
                     const size_t n,
                     forEachFcn_t f,
                     const bool dir = false,
                     void *arg = nullptr ) noexcept;
    }
}

#endif /*QLIBS_GENERIC*/