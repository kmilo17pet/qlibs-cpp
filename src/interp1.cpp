#include <include/interp1.hpp>
#include <include/ffmath.hpp>
#include <include/mathex.hpp>

using namespace qlibs;

/*cstat -CERT-INT30-C_a*/
/*============================================================================*/
bool interp1::setMethod( const interp1Method m )
{
    bool retValue = false;
    static const interp1Fcn_t im[ ] = {
        &interp1::next,
        &interp1::previous,
        &interp1::nearest,
        &interp1::linear,
        &interp1::sine,
        &interp1::cubic,
        &interp1::hermite,
        &interp1::spline,
        &interp1::cSpline,
    };

    if ( ( m > 0 ) && ( m < INTERP1_MAX ) ) {
        method = im[ m ];
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t interp1::next( const real_t x,
                      const real_t * const tx,
                      const real_t * const ty,
                      const size_t tableSize )
{
    real_t y = ffmath::getNan();

    if ( ( tableSize >= 2U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        size_t nearestIndex = tableSize - 1U;

        for ( size_t i = 0U ; i < tableSize - 1U ; ++i ) {
            if ( x < tx[ i + 1U ] ) {
                nearestIndex = i;
                break;
            }
        }
        
        if ( x >= tx[ tableSize - 1U ] ) {
            y = ty[ tableSize - 1U ];
        }
        else {
            y = ty[ nearestIndex + 1U ];
        }
    }

    return y;
}
/*============================================================================*/
real_t interp1::previous( const real_t x,
                          const real_t * const tx,
                          const real_t * const ty,
                          const size_t tableSize )
{
    real_t y = ffmath::getNan();

    if ( ( tableSize >= 2U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        size_t nearestIndex = 0U;

        if ( x <= tx[ 0 ] ) {
            y = ty[ 0 ];
        }
        else {
            for ( size_t i = 1U ; i < tableSize; ++i ) {
                if ( x < tx [ i ] ) {
                    break;
                }
                nearestIndex = i;
            }
            y = ty[ nearestIndex ];
        }
    }
    return y;
}
/*============================================================================*/
real_t interp1::nearest( const real_t x,
                         const real_t * const tx,
                         const real_t * const ty,
                         const size_t tableSize )
{
    real_t y = ffmath::getNan();

    if ( ( tableSize >= 2U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        size_t nearestIndex = 0U;
        real_t minDistance = ffmath::absf( x - tx[ 0 ] );

        for ( size_t i = 1U ; i < tableSize; ++i ) {
            const real_t distance = ffmath::absf( x - tx[ i ] );

            if ( distance <= minDistance) {
                minDistance = distance;
                nearestIndex = i;
            }
        }
        y = ty[ nearestIndex ];
    }
    return y;
}
/*============================================================================*/
real_t interp1::linear( const real_t x,
                        const real_t * const tx,
                        const real_t * const ty,
                        const size_t tableSize )
{
    real_t y = ffmath::getNan();

    if ( ( tableSize >= 2U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        if ( x < tx[ 0 ] ) {
            const real_t x0 = tx[ 0 ];
            const real_t x1 = tx[ 1 ];
            const real_t y0 = ty[ 0 ];
            const real_t y1 = ty[ 1 ];
            y = y0 + ( ( y1 - y0 )/( x1 - x0 ) )*( x - x0 );
        }
        else if ( x > tx[ tableSize - 1U ] ) {
            const real_t x0 = tx[ tableSize - 2U ];
            const real_t x1 = tx[ tableSize - 1U ];
            const real_t y0 = ty[ tableSize - 2U ];
            const real_t y1 = ty[ tableSize - 1U ];
            y = y1 + ( ( y0 - y1 )/( x0 - x1 ) )*( x - x1 );
        }
        else {
            const int maxIndex = static_cast<int>( tableSize ) - 1;
            for ( int i = 0; i < maxIndex; ++i ) {
                if ( ( x >= tx[ i ] ) && ( x <= tx[ i + 1 ] ) ) {
                    const real_t x0 = tx[ i ];
                    const real_t x1 = tx[ i + 1 ];
                    const real_t y0 = ty[ i ];
                    const real_t y1 = ty[ i + 1 ];
                    y = y0 + ( ( y1 - y0 )/( x1 - x0 ) )*( x - x0 );
                    break;
                }
            }
        }
    }
    return y;
}
/*============================================================================*/
real_t interp1::sine( const real_t x,
                      const real_t * const tx,
                      const real_t * const ty,
                      const size_t tableSize )
{
    real_t y = ffmath::getNan();
    if ( ( tableSize >= 2U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        if ( x < tx[ 0 ] ) {
            const real_t x0 = tx[ 0 ];
            const real_t x1 = tx[ 1 ];
            const real_t y0 = ty[ 0 ];
            const real_t y1 = ty[ 1 ];
            const real_t w = 0.5_re - 0.5_re*ffmath::cos( ffmath::FFP_PI*( x - x0 )/( x1 - x0 ) );
            y = y0 + w*( y1 - y0 );
        }
        else if ( x > tx[ tableSize - 1U ] ) {
            const real_t x0 = tx[ tableSize - 2U ];
            const real_t x1 = tx[ tableSize - 1U ];
            const real_t y0 = ty[ tableSize - 2U ];
            const real_t y1 = ty[ tableSize - 1U ];
            const real_t w = 0.5_re - 0.5_re*ffmath::cos( ffmath::FFP_PI*( x - x1 )/( x0 - x1 ) );
            y = y1 + w*( y0 - y1 );
        }
        else {
            for ( size_t i = 1; i < tableSize; ++i ) {
                if ( x <= tx[ i ] ) {
                    const real_t x0 = tx[ i - 1 ];
                    const real_t x1 = tx[ i ];
                    const real_t y0 = ty[ i - 1 ];
                    const real_t y1 = ty[ i];
                    const real_t w = 0.5_re - 0.5_re*ffmath::cos( ffmath::FFP_PI*( x - x0 )/( x1 - x0 ) );
                    y = y0 + w*( y1 - y0 );
                }
            }
        }
    }
    return y;
}
/*============================================================================*/
real_t interp1::cubic( const real_t x,
                       const real_t * const tx,
                       const real_t * const ty,
                       const size_t tableSize )
{
    real_t y = ffmath::getNan();
    if ( ( tableSize >= 4U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        if ( x < tx[ 0 ] ) {
            const real_t x0 = tx[ 0 ];
            const real_t x1 = tx[ 1 ];
            const real_t y0 = ty[ 0 ];
            const real_t y1 = ty[ 1 ];
            const real_t h = x1 - x0;
            const real_t t = ( x - x0 )/h;
            const real_t t2 = t*t;
            const real_t t3 = t2*t;

            y = ( 2.0_re*t3 - 3.0_re*t2 + 1.0_re )*y0 +
                ( t3 - 2.0_re*t2 + t )*h*( y0 - y1 ) +
                ( -2.0_re*t3 + 3.0_re*t2 )*y1 +
                ( t3 - t2 )*h*( y1 - y0 );
        }
        else if ( x > tx[ tableSize - 1U ] ) {
            const real_t x0 = tx[ tableSize - 2U ];
            const real_t x1 = tx[ tableSize - 1U ];
            const real_t y0 = ty[ tableSize - 2U ];
            const real_t y1 = ty[ tableSize - 1U ];
            const real_t h = x1 - x0;
            const real_t t = ( x - x1 )/h;
            const real_t t2 = t*t;
            const real_t t3 = t2*t;

            y = ( 2.0_re*t3 - 3.0_re*t2 + 1.0_re )*y1 +
                ( t3 - 2.0_re*t2 + t )*h*( y0 - ty[ tableSize - 3U ] ) +
                ( -2.0_re*t3 + 3.0_re*t2 )*y0 +
                ( t3 - t2 )*h*( y1 - y0 );
        }
        else {
            for ( size_t i = 1U; i < tableSize; ++i ) {
                if ( x <= tx[ i ] ) {
                    const real_t x0 = tx[ i - 1U ];
                    const real_t x1 = tx[ i ];
                    const real_t y0 = ty[ i - 1U ];
                    const real_t y1 = ty[ i ];
                    const real_t h = x1 - x0;
                    const real_t t = ( x - x0 )/h;
                    const real_t t2 = t*t;
                    const real_t t3 = t2*t;
                    y = ( 2.0_re*t3 - 3.0_re*t2 + 1.0_re )*y0 +
                        ( t3 - 2.0_re*t2 + t )*h*( y0 - ty[ i - 2U ] ) +
                        ( -2.0_re*t3 + 3.0_re*t2 )*y1 +
                        ( t3 - t2 )*h*( y1 - y0 );
                    break;
                }
            }
        }
    }
    return y;
}
/*============================================================================*/
real_t interp1::hermite( const real_t x,
                         const real_t * const tx,
                         const real_t * const ty,
                         const size_t tableSize )
{
    real_t y = ffmath::getNan();

    if ( ( tableSize >= 2U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        if ( x < tx[ 0 ] ) {
            const real_t x0 = tx[ 0 ];
            const real_t x1 = tx[ 1 ];
            const real_t y0 = ty[ 0 ];
            const real_t y1 = ty[ 1 ];
            y = y0 + ( ( y1 - y0 )/( x1 - x0 ) )*( x - x0 );
        }
        else if ( x > tx[ tableSize -1U ] ) {
            const real_t x0 = tx[ tableSize - 2U ];
            const real_t x1 = tx[ tableSize - 1U ];
            const real_t y0 = ty[ tableSize - 2U ];
            const real_t y1 = ty[ tableSize - 1U ];
            y = y1 + ( ( y0 - y1 )/( x0 - x1 ) )*( x - x1 );
        }
        else {
            y = 0.0_re;
            for ( size_t i = 0U ; i < tableSize; ++i ) {
                real_t term = ty[ i ];

                for ( size_t j = 0U ; j < tableSize; ++j ) {
                    if ( i != j ) {
                        term *= ( x - tx[ j ] )/( tx[ i ] - tx[ j ] );
                    }
                }

                y += term;
            }
        }
    }
    return y;
}
/*============================================================================*/
real_t interp1::spline( const real_t x,
                        const real_t * const tx,
                        const real_t * const ty,
                        const size_t tableSize )
{
    real_t y = ffmath::getNan();

    if ( ( tableSize >= 4U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        size_t i = 0U;
        /* Extrapolation for x beyond the range*/
        if ( x <= tx[ 0 ] ) {
            i = 0U;
        }
        else if ( x >= tx[ tableSize - 1U ] ) {
            i = tableSize - 2U; /* Use the last interval for extrapolation*/
        }
        else {
          while ( x >= tx[ i + 1U ] ) { 
            i++;
          }
        }

        if ( isEqual( x , tx[ i + 1U ] ) ) { 
            y = ty[ i + 1U ];
        }
        else {
            const real_t t = ( x - tx[ i ] )/( tx[ i + 1U ] - tx[ i ] );
            const real_t t_2 = t*t;
            const real_t t_3 = t_2*t;
            const real_t h00 = 2.0_re*t_3 - 3.0_re*t_2 + 1.0_re;
            const real_t h10 = t_3 - 2.0_re*t_2 + t;
            const real_t h01 = 3.0_re * t_2 - 2.0_re*t_3;
            const real_t h11 = t_3 - t_2;
            const real_t x1_x0 = tx[ i + 1U ] - tx[ i ];
            const real_t y0 = ty[ i ];
            const real_t y1 = ty[ i + 1 ];
            real_t m0, m1;

            if ( 0U == i ) {
                m0 = ( ty[ 1 ] - ty[ 0 ] )/( tx[ 1 ] - tx[ 0 ] );
                m1 = ( ty[ 2 ] - ty[ 0 ] )/( tx[ 2 ] - tx[ 0 ] );
            }
            else if ( ( tableSize - 2U ) == i ) {
                m0 = ( ty[ tableSize - 1U ] - ty[ tableSize - 3U ] )/( tx[ tableSize - 1U ] - tx[ tableSize - 3U ] );
                m1 = ( ty[ tableSize - 1U ] - ty[ tableSize - 2U ] )/( tx[ tableSize - 1U ] - tx[ tableSize - 2U ] );
            }
            else {
                m0 = slope( tx, ty, i );
                m1 = slope( tx, ty, i + 1U );
            }
            y = ( h00*y0 ) + ( h01*y1 ) + ( h10*x1_x0*m0 ) + ( h11*x1_x0*m1 );
        }
    }

    return y;
}
/*============================================================================*/
real_t interp1::cSpline( const real_t x,
                         const real_t * const tx,
                         const real_t * const ty,
                         const size_t tableSize )
{
    real_t y = ffmath::getNan();

    if ( ( tableSize >= 4U ) && ( nullptr != tx ) && ( nullptr != ty ) ) {
        size_t i = 0U;
        /* Extrapolation for x beyond the range*/
        if ( x <= tx[ 0 ] ) {
            i = 0U;
        }
        else if ( x >= tx[ tableSize - 1U ] ) {
            i = tableSize - 2U; /* Use the last interval for extrapolation*/
        }
        else {
          while ( x >= tx[ i + 1U ] ) {
            i++;
          }
        }
        
        if ( isEqual( x , tx[ i + 1U ] ) ) {
            y = ty[ i + 1U ];
        }
        else {
            const real_t x0 = tx[ i ];
            const real_t x1 = tx[ i + 1U ];
            const real_t y0 = ty[ i ];
            const real_t y1 = ty[ i + 1U ];
            const real_t fd2i_xl1 = leftSecondDerivate( tx, ty, tableSize - 1, i + 1 );
            const real_t fd2i_x = rightSecondDerivate( tx, ty, tableSize - 1, i + 1 );
            const real_t x0_x1 = x0 - x1;
            const real_t x1_2 = x1*x1;
            const real_t x1_3 = x1_2*x1;
            const real_t x0_2 = x0*x0;
            const real_t x0_3 = x0_2*x0;
            const real_t d = ( fd2i_x - fd2i_xl1 )/( 6.0_re*x0_x1 );
            const real_t c = ( ( x0*fd2i_xl1 ) - ( x1*fd2i_x ) )/( 2.0_re*x0_x1 );
            const real_t b = ( y0 - y1 - c*( x0_2 - x1_2 ) - d*( x0_3 - x1_3 ) )/x0_x1;
            const real_t a = y1 - ( b*x1 ) - ( c*x1_2 ) - ( d*x1_3 );
            y = a + x*( b + x*( c + ( x*d ) ) );
        }
    }

    return y;
}
/*============================================================================*/
real_t interp1::slope( const real_t * const tx,
                       const real_t * const ty,
                       const size_t i )
{
    real_t m;

    if ( isEqual( tx[ i + 1U ], tx[ i - 1U ] ) ) {
        m = 0.0_re;
    }
    else {
        m = ( ty[ i + 1U ] - ty[ i - 1U ] )/( tx[ i + 1U ] - tx[ i - 1U ] );
    }

    return m;
}
/*============================================================================*/
real_t interp1::firstDerivate( const real_t * const tx,
                               const real_t * const ty,
                               const size_t n,
                               const size_t i )
{
    real_t fd1_x;

    if ( 0 == i ) {
        const real_t dx = tx[ 1 ] - tx[ 0 ];
        const real_t dy = ty[ 1 ] - ty[ 0 ];
        fd1_x = 1.5_re*( dy/dx );
        fd1_x -= 1.0_re/( ( tx[ 2 ] - tx[ 0 ] )/( ty[ 2 ] - ty[ 0 ] ) + dx/dy );
    }
    else if ( n == i ) {
        const real_t dx = tx[ n ] - tx[ n - 1 ];
        const real_t dy = ty[ n ] - ty[ n - 1 ];
        fd1_x = 1.5_re*( dy/dx );
        fd1_x -= 1.0_re/( ( tx[ n ] - tx[ n - 2 ] )/( ty[ n ] - ty[ n - 2 ] ) + ( dx/dy ) );
    }
    else {
        const real_t tmp1 = ( tx[ i + 1 ] - tx[ i ] )/( ty[ i + 1 ] - ty[ i ] );
        const real_t tmp2 = ( tx[ i ] - tx[ i - 1 ] )/( ty[ i ] - ty[ i - 1 ] );

        if ( ( tmp1*tmp2 ) < 0.0_re ) {
            fd1_x = 0.0_re;
        }
        else{
            fd1_x = 2.0_re/( tmp1 + tmp2 );
        }
    }
    return fd1_x;
}
/*============================================================================*/
real_t interp1::leftSecondDerivate( const real_t * const tx,
                                    const real_t * const ty,
                                    const size_t n,
                                    const size_t i )
{
    const real_t fdi_x = firstDerivate( tx, ty, n, i );
    const real_t fdi_xl1 = firstDerivate( tx, ty, n, i - 1U );
    const real_t xi_delta = tx[ i ] - tx[ i - 1U ];
    real_t fd2l_x = -2.0_re*( fdi_x + ( 2.0_re*fdi_xl1 ) )/xi_delta;

    fd2l_x += 6.0_re*( ty[ i ] - ty[ i - 1U ] )/( xi_delta*xi_delta );
    return fd2l_x;
}
/*============================================================================*/
real_t interp1::rightSecondDerivate( const real_t * const tx,
                                     const real_t * const ty,
                                     const size_t n,
                                     const size_t i )
{
    const real_t fdi_x = firstDerivate( tx, ty, n, i );
    const real_t fdi_xl1 = firstDerivate( tx, ty, n, i - 1U );
    const real_t xi_delta = tx[ i ] - tx[ i - 1U ];
    real_t fd2r_x = 2.0_re*( ( 2.0_re*fdi_x ) + fdi_xl1 )/xi_delta;

    fd2r_x -= 6.0_re*( ty[ i ] - ty[ i - 1U ] )/( xi_delta*xi_delta );
    return fd2r_x;
}
/*============================================================================*/
/*cstat +CERT-INT30-C_a*/
