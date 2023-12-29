#include <include/rms.hpp>
#include <include/ffmath.hpp>

using namespace qlibs;

/*============================================================================*/
bool rms::setup( real_t * const window,
                 const size_t wsize ) noexcept
{
    bool retValue = false;

    if ( ( nullptr != window ) && ( wsize > 0U ) ) {
        (void)smootherEXPW::setup( 0.99_re );
        (void)smootherMWM2::setup( window, wsize );
        (void)smootherLPF1::setup( 0.75_re );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t rms::update( const real_t x ) noexcept
{
    real_t y;

    y = ffmath::sqrt( smootherEXPW::smooth( x*x ) );
    y = smootherMWM2::smooth( y );
    y = smootherLPF1::smooth( y );

    return y;
}
/*============================================================================*/
bool rms::setParams( const real_t l,
                     const real_t a ) noexcept
{
    bool retValue = false;

    if ( ( l > 0.0_re ) && ( l <= 1.0_re ) && ( a > 0.0_re ) && ( a <= 1.0_re ) ) {
        lambda = l;
        alpha = a;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/