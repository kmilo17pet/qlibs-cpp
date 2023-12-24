#include <include/rms.hpp>

using namespace qlibs;

/*============================================================================*/
bool rms::setup( real_t * const window, const size_t wsize )
{
    bool retValue = false;

    if ( ( nullptr != window ) && ( wsize > 0U) ) {
        (void)smootherEXPW::setup( 0.99 );
        (void)smootherMWM2::setup( window, wsize );
        (void)smootherLPF1::setup( 0.75 );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t rms::update( const real_t x )
{
    real_t y;

    y = sqrt( smootherEXPW::smooth( x*x ) );
    y = smootherMWM2::smooth( y );
    y = smootherLPF1::smooth( y );

    return y;
}
/*============================================================================*/
bool rms::setParams( const real_t l, const real_t a )
{
    bool retValue = false;

    if ( ( l > 0.0 ) && ( l <= 1.0 ) && ( a > 0.0 ) && ( a <= 1.0 ) ) {
        lambda = l;
        alpha = a;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/