#include "include/mathex.hpp"
#include <cmath>
#include <algorithm>

namespace qlibs {

/*============================================================================*/
real_t normalize( const real_t x,
                  const real_t xMin,
                  const real_t xMax )
{
    return ( x - xMin )/( xMax - xMin );
}
/*============================================================================*/
real_t mapMinMax( const real_t x,
                  const real_t xMin,
                  const real_t xMax,
                  const real_t yMin,
                  const real_t yMax )
{
    return ( ( yMax - yMin )*normalize( x, xMin, xMax ) ) + yMin;
}
/*============================================================================*/
bool inRangeCoerce( real_t &x,
                    const real_t lowerL,
                    const real_t upperL )
{
    bool retVal = false;

    if ( isnan( x ) ) {
        x = lowerL;
    }
    else {
        if ( x < lowerL ) {
            x = lowerL;
        }
        else if ( x > upperL ) {
            x = upperL;
        }
        else {
            retVal = true;
        }
    }

    return retVal;
}
/*============================================================================*/
bool isEqual( const real_t a, const real_t b, const real_t tol )
{
    return ( fabs( a - b ) <= fabs( tol ) );
}
/*============================================================================*/
bool inPolygon( const real_t x,
                const real_t y,
                const real_t * const px,
                const real_t * const py,
                const size_t p )
{
    size_t i;
    bool retVal = false;
    real_t max_y = py[ 0 ], max_x = px[ 0 ], min_y = py[ 0 ], min_x = px[ 0 ];

    for ( i = 0u ; i < p ; ++i ) {
        max_y = std::max( py[ i ], max_y );
        max_x = std::max( px[ i ], max_x );
        min_y = std::min( py[ i ], min_y );
        min_x = std::min( px[ i ], min_x );
    }

    if ( ( y >= min_y ) && ( y <= max_y ) && ( x >= min_x ) && ( x <= max_x ) ) {
        size_t j = p - 1u;

        for ( i = 0u ; i < p ; ++i ) {
            if ( ( px[ i ] > x ) != ( px[ j ] > x ) ) {
                const real_t dx = px[ j ] - px[ i ];
                const real_t dy = py[ j ] - py[ i ];
                if ( y < ( ( dy*( x - px[ i ] ) )/( dx + py[ i ] ) ) ) {
                    retVal = !retVal;
                }
            }
            j = i;
        }
    }

    return retVal;
}
/*============================================================================*/
bool inCircle( const real_t x,
               const real_t y,
               const real_t cx,
               const real_t cy,
               const real_t r )
{
    const real_t d = ( ( x - cx )*( x - cx ) ) + ( ( y - cy )*( y - cy ) );
    return ( d <= (r*r) );
}
/*============================================================================*/

} /*namespace qlibs*/