#include <include/numa.hpp>

using namespace qlibs;

/*===========================================================================*/
void nState::init( const real_t x0,
                   const real_t sn_1,
                   const real_t sn_2 ) noexcept
{
    x[ 0 ] = x0;
    x[ 1 ] = sn_1;
    x[ 2 ] = sn_2;
}
/*===========================================================================*/
real_t nState::integrate( const real_t s,
                          const real_t dt,
                          const bool bUpdate ) noexcept
{
    switch( iMethod ) {
        case INTEGRATION_RECTANGULAR:
            x[ 0 ] += s*dt;
            break;
        case INTEGRATION_TRAPEZOIDAL:
            x[ 0 ] += 0.5_re*( s + x[ 1 ] )*dt;
            break;
        case INTEGRATION_SIMPSON:
            x[ 0 ] += ( 1.0_re/6.0_re )*( s + ( 4.0_re*x[ 1 ] ) + x[ 2 ] )*dt;
            break;
        case INTEGRATION_QUADRATIC:
            x[ 0 ] += ( 1.0_re/12.0_re )*( ( 5.0_re*s ) + ( 8.0_re*x[ 1 ] ) - x[ 2 ] )*dt;
            break;
        default:
            break;
    }
    if ( bUpdate ) {
        update( s );
    }

    return x[ 0 ];
}
/*===========================================================================*/
real_t nState::derive( const real_t s,
                       const real_t dt,
                       const bool bUpdate ) noexcept
{
    float ds = 0.0F;

    switch( dMethod ) {
        case DERIVATION_2POINTS:
            ds = ( s - x[ 1 ] )/dt;
            break;
        case DERIVATION_BACKWARD:
            ds = ( ( 3.0_re*s ) - ( 4.0_re*x[ 1 ] ) + x[ 2 ] )/( 2.0_re*dt );
            break;
        case DERIVATION_FORWARD:
            ds = ( ( 4.0_re*x[ 1 ] ) - ( 3.0_re*x[ 2 ] ) - s )/( 2.0_re*dt );
            break;
        default:
            break;
    }
    if ( bUpdate ) {
        update( s );
    }

    return ds;
}
/*===========================================================================*/