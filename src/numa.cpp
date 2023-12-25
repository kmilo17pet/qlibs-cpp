#include "include/numa.hpp"

using namespace qlibs;

/*===========================================================================*/
void state::init( const real_t x0,
                  const real_t sn_1,
                  const real_t sn_2 ) noexcept
{
    x[ 0 ] = x0;
    x[ 1 ] = sn_1;
    x[ 2 ] = sn_2;
}
/*===========================================================================*/
real_t state::integrate( const real_t s,
                         const real_t dt ) noexcept
{
    switch( intMethod ) {
        case INTEGRATION_RECTANGULAR:
            x[ 0 ] += s*dt;
            break;
        case INTEGRATION_TRAPEZOIDAL:
            x[ 0 ] += 0.5*( s + x[ 1 ] )*dt;
            break;
        case INTEGRATION_SIMPSON:
            x[ 0 ] += ( 1.0/6.0 )*( s + ( 4.0*x[ 1 ] ) + x[ 2 ] )*dt;
            break;
        default:
            break;
    }
    update( s );

    return x[ 0 ];
}
/*===========================================================================*/
real_t state::derivative( const real_t s,
                          const real_t dt ) noexcept
{
    x[ 0 ] = ( s - x[ 1 ] )/dt;
    update( s );

    return x[ 0 ];
}
/*===========================================================================*/