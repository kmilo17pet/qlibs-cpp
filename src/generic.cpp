#include <include/generic.hpp>

using namespace qlibs;

struct sortStackNode {
    uint8_t *lo;
    uint8_t *hi;
};

static void sortStackPush( sortStackNode **top,
                           uint8_t *low,
                           uint8_t *high ) noexcept;

/*cstat -CERT-EXP36-C_b*/

/*============================================================================*/
void generic::swap( void * const x,
                    void * const y,
                    size_t n ) noexcept
{
    uint8_t * const a = static_cast<uint8_t *>( x );
    uint8_t * const b = static_cast<uint8_t *>( y );
    size_t i = 0U, j = 0U;
    do {
        const uint8_t tmp = a[ i ];
        /*cstat -CERT-INT30-C_a*/
        a[ i++ ] = b[ j ];
        b[ j++ ] = tmp;
        /*cstat +CERT-INT30-C_a*/
    } while( --n > 0U );
}
/*============================================================================*/
static void sortStackPush( sortStackNode **top,
                           uint8_t *low,
                           uint8_t *high ) noexcept
{
    top[ 0 ]->lo = low;
    top[ 0 ]->hi = high;
    ++top[ 0 ];
}
/*============================================================================*/
void generic::sort( void * const pbase,
                    size_t n,
                    size_t size,
                    generic::compareFcn_t cmp,
                    void *arg ) noexcept
{
    if ( ( nullptr != pbase ) && ( size > 0U ) && ( n > 0U ) && ( nullptr != cmp ) ) {
        const size_t max_thresh = 4U*size;
        uint8_t * const base_ptr = static_cast<uint8_t *>( pbase );
        uint8_t * const end_ptr = &base_ptr[ size*( n - 1U ) ];
        uint8_t * tmp_ptr = base_ptr, *run_ptr;
        uint8_t * const thresh = ( end_ptr < ( base_ptr + max_thresh) ) ? end_ptr : ( base_ptr + max_thresh ) ;

        if ( n > 4U ) {
            uint8_t *lo = base_ptr, *hi = &lo[ size*( n - 1U ) ];
            sortStackNode stack[ 8U*sizeof(size_t) ], *top = stack;

            sortStackPush( &top, nullptr, nullptr );
            while ( stack < top ) {
                uint8_t *left_ptr, *right_ptr;
                uint8_t *mid = &lo[ size*( ( static_cast<size_t>( hi - lo )/size ) >> 1U ) ];
                
                if ( cmp( mid, lo, arg ) < 0 ) {
                    generic::swap( mid, lo, size );
                }
                if ( cmp( hi, mid, arg ) < 0 ) {
                    generic::swap( mid, hi, size );
                }
                else {
                    goto jump_over; /*MISRAC deviation allowed*/
                }
                if ( cmp( mid, lo, arg ) < 0 ) {
                    generic::swap( mid, lo, size );
                }

                jump_over:
                left_ptr  = lo + size;
                right_ptr = hi - size;
                do {
                    while ( cmp( left_ptr, mid, arg ) < 0 ) {
                        left_ptr += size;
                    }
                    while ( cmp( mid, right_ptr, arg ) < 0 ) {
                        right_ptr -= size;
                    }

                    if ( left_ptr < right_ptr ) {
                        generic::swap( left_ptr, right_ptr, size );
                        if ( mid == left_ptr ) {
                            mid = right_ptr;
                        }
                        else if ( mid == right_ptr ) {
                            mid = left_ptr;
                        }
                        else {
                            /*nothing to do here*/
                        }
                        left_ptr += size;
                        right_ptr -= size;
                    }
                    else if ( left_ptr == right_ptr ) {
                        left_ptr += size;
                        right_ptr -= size;
                        break;
                    }
                    else {
                        /*nothing to do here*/
                    }
                } while (left_ptr <= right_ptr);

                if ( static_cast<size_t>( right_ptr - lo ) <= max_thresh ) {
                    if ( static_cast<size_t>( hi - left_ptr ) <= max_thresh ) {
                        --top; /*POP form the stack*/
                        lo = top->lo;
                        hi = top->hi;
                    }
                    else {
                        lo = left_ptr;
                    }
                }
                else if ( static_cast<size_t>( hi - left_ptr ) <= max_thresh ) {
                    hi = right_ptr;
                }
                else if ( ( right_ptr - lo ) > ( hi - left_ptr ) ) {
                    sortStackPush( &top, lo, right_ptr );
                    lo = left_ptr;
                }
                else {
                    sortStackPush( &top, left_ptr, hi );
                    hi = right_ptr;
                }
            }
        }
        for ( run_ptr = tmp_ptr + size ; run_ptr <= thresh ; run_ptr += size ) {
            if ( cmp( run_ptr, tmp_ptr, arg ) < 0 ) {
                tmp_ptr = run_ptr;
            }
        }

        if ( tmp_ptr != base_ptr ) {
            generic::swap( tmp_ptr, base_ptr, size );
        }
        run_ptr = base_ptr + size;
        /*cstat -MISRAC++2008-6-2-1*/
        while ( ( run_ptr += size ) <= end_ptr ) {
        /*cstat +MISRAC++2008-6-2-1*/
            tmp_ptr = run_ptr - size;
            while ( cmp( run_ptr, tmp_ptr, arg ) < 0 ) {
                tmp_ptr -= size;
            }
            tmp_ptr += size;
            if ( tmp_ptr != run_ptr ) {
                uint8_t *tra = run_ptr + size;
                while ( --tra >= run_ptr ) {
                    const uint8_t c = *tra;
                    uint8_t *hi = tra;
                    uint8_t *lo = tra;
                    /*cstat -MISRAC++2008-6-2-1*/
                    while ( (lo -= size) >= tmp_ptr ) {
                    /*cstat +MISRAC++2008-6-2-1*/
                        *hi = *lo;
                        hi = lo;
                    }
                    *hi = c;
                }
            }
        }
    }
}
/*============================================================================*/
void generic::reverse( void * const pbase,
                       const size_t size,
                       const size_t init,
                       const size_t end ) noexcept
{
    if ( ( nullptr != pbase ) && ( size > 0U ) && ( end > init ) ) {
        size_t s = size*init, e = size*end;
        uint8_t * const v = static_cast<uint8_t*>( pbase );
        
        while( s < e ) {
            generic::swap( &v[ s ], &v[ e ], size );
            s += size;
            e -= size;
        }
    }
}
/*============================================================================*/
void generic::rotate( void * const pbase,
                      const size_t size,
                      const size_t n,
                      const int k ) noexcept
{
    if ( ( nullptr != pbase ) && ( 0 != k ) && ( n > 0U ) ) {
        size_t r;

        if ( k > 0 ) {
            r = static_cast<size_t>( k );
            r %= n;
            generic::reverse( pbase, size, n - r, n - 1U );
            generic::reverse( pbase, size, 0U, n - r - 1U );
            generic::reverse( pbase, size, 0U, n - 1U );
        }
        else {
            /*cstat -MISRAC++2008-5-0-9*/
            r = static_cast<size_t>( -k );
            /*cstat +MISRAC++2008-5-0-9*/
            r %= n;
            generic::reverse( pbase, size, 0U, r - 1U );
            generic::reverse( pbase, size, r, n - 1U );
            generic::reverse( pbase, size, 0U, n - 1U );
        }
    }
}
/*============================================================================*/
void* generic::set( void * const pbase,
                    const size_t size,
                    const size_t n,
                    const void * const ref ) noexcept
{
    uint8_t * const p = static_cast<uint8_t*>( pbase );
    void *retVal = nullptr;

    if ( ( nullptr != pbase ) && ( size > 0U ) && ( n > 0U ) && ( nullptr != ref ) ) {
        for ( size_t i = 0U ; i < n ; i++ ) {
            retVal = memcpy( &p[ size*i ], ref, size );
        }
    }

    return retVal;
}
/*============================================================================*/
void* generic::lSearch( const void *key,
                        const void *pbase,
                        const size_t n,
                        const size_t size,
                        generic::compareFcn_t compar,
                        void *arg ) noexcept
{
    const uint8_t * const pb = static_cast<const uint8_t *>( pbase );
    void* retVal = nullptr;

    for ( size_t i = 0U ; i < n ; ++i ) {
        const uint8_t * const element = &pb[ i*size ];

        if ( 0 == compar( key, element, arg ) ) {
            retVal = reinterpret_cast<void*>( const_cast<uint8_t*>( element ) );
            break;
        }
    }
    return retVal;
}
/*============================================================================*/
void* generic::bSearch( const void *key,
                        const void *pbase,
                        const size_t n,
                        const size_t size,
                        generic::compareFcn_t compar,
                        void *arg ) noexcept
{
    const uint8_t *base = static_cast<const uint8_t *>( pbase );
    size_t lim  = n;
    const uint8_t *p;
    void *retVal = nullptr;

    while ( 0U != lim ) {
        int cmp;

        p = &base[ ( lim >> 1U )*size ];
        cmp = compar( key, p, arg );
        if ( 0 == cmp ) {
            retVal = reinterpret_cast<void*>( const_cast<uint8_t*>( p ) );
            break;
        }
        else if ( cmp > 0 ) {
            base = &p[ size ];
            lim--;
        }
        else {
            /*nothing to do here*/
        }
        lim >>= 1U;
    }
    return retVal;
}
/*============================================================================*/
int generic::forEach( void *pbase,
                      const size_t size,
                      const size_t n,
                      generic::forEachFcn_t f,
                      const bool dir,
                      void *arg ) noexcept
{
    int retVal = 0;

    if ( ( nullptr != pbase ) && ( nullptr != f ) && ( n > 0U ) ) {
        uint8_t * const pb = static_cast<uint8_t *>( pbase );
        
        if ( 1 != f( -1, nullptr, arg ) ) {
            size_t i;
            uint8_t *element;

            if ( !dir ) {
                for ( i = 0U ; i < n ; ++i ) {
                    element = &pb[ i*size ];
                    /*cstat -MISRAC++2008-5-0-9*/
                    retVal = f( static_cast<int>( i ), element, arg );
                    /*cstat +MISRAC++2008-5-0-9*/
                    if ( 1 == retVal ) {
                        break;
                    }
                }
            }
            else {
                i = n;
                while ( i-- > 0U ) {
                    element = &pb[ i*size ];
                    /*cstat -MISRAC++2008-5-0-9*/
                    retVal = f( static_cast<int>( i ), element, arg );
                    /*cstat +MISRAC++2008-5-0-9*/
                    if ( 1 == retVal ) {
                        break;
                    }
                }
            }
        }
        /*cstat -MISRAC++2008-5-0-9*/
        (void)f( static_cast<int>( n ), nullptr, arg );
        /*cstat +MISRAC++2008-5-0-9*/
    }

    return retVal;
}
/*============================================================================*/

/*cstat +CERT-EXP36-C_b*/