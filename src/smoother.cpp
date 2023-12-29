#include <include/smoother.hpp>
#include <include/ltisys.hpp>
#include <include/ffmath.hpp>

using namespace qlibs;

/*============================================================================*/
void smoother::windowSet( real_t *w,
                          const size_t wsize,
                          const real_t x )
{
    for ( size_t i = 0U ; i < wsize ; ++i ) {
        w[ i ] = x;
    }
}
/*============================================================================*/
bool smootherLPF1::setup( const real_t a )
{
    bool retValue = false;

    if ( ( a > 0.0_re ) && ( a < 1.0_re ) ) {
        alpha = a;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherLPF1::smooth( const real_t x )
{
    real_t y;

    if ( init ) {
        y1 = x;
        init = false;
    }
    y = x + ( alpha*( y1 - x ) );
    y1 = y;

    return y;
}
/*============================================================================*/
bool smootherLPF2::setup( const real_t a )
{
    bool retValue = false;

    if ( ( a > 0.0_re ) && ( a < 1.0_re ) ) {
        real_t aa, p1, r;
        aa = a*a;
        /*cstat -MISRAC2012-Dir-4.11_b*/
        p1 = ffmath::sqrt( 2.0_re*a ); /*arg always positive*/
        /*cstat +MISRAC2012-Dir-4.11_b*/
        r = 1.0_re + p1 + aa;
        k = aa/r;
        a1 = 2.0_re*( aa - 1.0_re )/r;
        a2 = ( 1.0_re - p1 + aa )/r;
        b1 = 2.0_re*k;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherLPF2::smooth( const real_t x )
{
    real_t y;

    if ( init ) {
        y1 = x;
        y2 = x;
        x1 = x;
        x2 = x;
        init = false;
    }
    y = ( k*x ) + ( b1*x1 ) + ( k*x2 ) - ( a1*y1 ) - ( a2*y2 );
    x2 = x1;
    x1 = x;
    y2 = y1;
    y1 = y;

    return y;
}
/*============================================================================*/
bool smootherMWM1::setup( real_t *window,
                          const size_t w_size )
{
    bool retValue = false;

    if ( ( nullptr != window ) && ( wsize > 0U ) ) {
        w = window;
        wsize = w_size;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherMWM1::smooth( const real_t x )
{
    if ( init ) {
        windowSet( w, wsize, x );
        init = false;
    }

    return discreteSystem::updateFIR( w, wsize, x )/static_cast<real_t>( wsize );
}
/*============================================================================*/
bool smootherMWM2::setup( real_t *window,
                          const size_t w_size )
{
    bool retValue = false;

    if ( ( nullptr != window ) && ( w_size > 0U ) ) {
        tdl::setup( window, w_size );
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherMWM2::smooth( const real_t x )
{
    const real_t wsize = static_cast<real_t>( itemCount );

    if ( init ) {
        flush( x );
        sum = x*wsize;
        init = false;
    }
    sum += x - getOldest();
    insertSample( x );

    return sum/wsize;
}
/*============================================================================*/
bool smootherMOR1::setup( real_t *window,
                          const size_t w_size,
                          const real_t a )
{
    bool retValue = false;

    if ( ( nullptr != window ) && ( w_size > 0U ) && ( a > 0.0_re ) && ( a < 1.0_re ) ) {
        w = window;
        wsize = w_size;
        alpha = a;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherMOR1::smooth( const real_t x )
{
    real_t mc;

    if ( init ) {
        windowSet( w, wsize, x );
        m = x;
        init = false;
    }
    /*shift, sum and compensate*/
    mc = discreteSystem::updateFIR( w, wsize, x ) - x;
    if ( ffmath::absf( m - x ) > ( alpha*ffmath::absf( m ) ) ) {
        w[ 0 ] = m; /*replace the outlier with the dynamic median*/
    }
    /*compute new mean for next iteration*/
    m = ( mc + w[ 0 ] ) / static_cast<real_t>( wsize );
    return w[ 0 ];
}
/*============================================================================*/
bool smootherMOR2::setup( real_t *window,
                          const size_t w_size,
                          const real_t a )
{
    bool retValue = false;

    if ( ( nullptr != window ) && ( w_size > 0U ) && ( a > 0.0_re ) && ( a < 1.0_re ) ) {
        alpha = a;
        tdl::setup( window, w_size );
        sum = 0.0_re;
        m = 0.0_re;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherMOR2::smooth( const real_t x )
{
    real_t xx = x;
    const real_t wsize = static_cast<real_t>( itemCount );

    if ( init ) {
        flush( x );
        sum = wsize*x;
        m = x;
        init = false;
    }
    /*is it an outlier?*/
    if ( ffmath::absf( m - x ) > ( alpha*ffmath::absf( m ) ) ) {
        xx = m; /*replace the outlier with the dynamic median*/
    }
    sum += xx - getOldest();
    m = sum/wsize;
    insertSample( xx );

    return x;
}
/*============================================================================*/
bool smootherGMWF::setup( const real_t sg,
                          const real_t c,
                          real_t *window,
                          real_t *kernel,
                          const size_t wk_size )
{
    bool retValue = false;
    const size_t ws = wk_size;

    if ( ( nullptr != window ) && ( wk_size > 0U ) && ( c < static_cast<real_t>( ws ) ) && ( sg > 0.0_re ) ) {
        real_t r, sum = 0.0_re;
        size_t i;
        real_t l, center;
        /*cstat -MISRAC++2008-5-0-7*/
        l = static_cast<real_t>( wk_size - 1U );
        /*cstat +MISRAC++2008-5-0-7*/
        center = c - l;
        r = 2.0_re*sg*sg;
        for ( i = 0U ; i < ws ; ++i ) {
            /*cstat -MISRAC++2008-5-0-7*/
            real_t d = static_cast<real_t>( i ) - l; /*symmetry*/
            /*cstat +MISRAC++2008-5-0-7*/
            d -= center;
            /*cstat -CERT-FLP32-C_b*/
            kernel[ i ] = ffmath::exp( -( d*d )/r );
            /*cstat +CERT-FLP32-C_b*/
            sum += kernel[ i ];
        }
        for ( i = 0U ; i < ws ; ++i ) {
            kernel[ i ] /= sum;
        }
        w = window;
        k = kernel;
        wsize = ws;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherGMWF::smooth( const real_t x )
{
    if ( init ) {
        windowSet( w, wsize, x );
        init = false;
    }
    return discreteSystem::updateFIR( w, wsize, x, k );
}
/*============================================================================*/
bool smootherEXPW::setup( const real_t lam )
{
    bool retValue = false;
    
    if ( ( lam > 0.0_re ) && ( lam < 1.0_re ) ) {
        lambda = lam;
        m = 0.0_re;
        w = 1.0_re;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherEXPW::smooth( const real_t x )
{
    real_t iw;

    if ( init ) {
        m = x;
        w = 1.0_re;
        init = false;
    }
    w = ( lambda*w ) + 1.0_re;
    iw = 1.0_re/w;
    m = ( m*( 1.0_re - iw ) ) + ( iw*x );

    return m;
}
/*============================================================================*/
bool smootherKLMN::setup( const real_t processNoiseCov,
                          const real_t measureNoiseCov,
                          const real_t estErrorCov )
{
    bool retValue = false;

    if ( ( processNoiseCov > 0.0_re ) && ( measureNoiseCov > 0.0_re ) && ( estErrorCov > 0.0_re ) ) {
        p = estErrorCov;
        q = processNoiseCov;
        r = measureNoiseCov;
        A = 1.0_re;
        H = 1.0_re;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherKLMN::smooth( const real_t x )
{
    real_t pH;

    if ( init ) {
        xS = x;
        init = false;
    }
    /* Predict */
    xS = A*xS;
    p = ( A*A*p ) + q; /* p(n|n-1)=A^2*p(n-1|n-1)+q */
    /* Measurement */
    pH = p*H;
    gain =  pH/( r + ( H*pH ) );
    xS += gain*( x - ( H*xS ) );
    p = ( 1.0_re - ( gain*H ) )*p; /*covariance update*/

    return xS;
}
/*============================================================================*/
bool smootherDESF::setup( const real_t a,
                          const real_t b,
                          const size_t nS )
{
    bool retValue = false;

    if ( ( a > 0.0_re ) && ( a < 1.0_re ) && ( b > 0.0_re ) && ( b < 1.0_re ) ) {
        alpha = a;
        beta = b;
        n = static_cast<real_t>( nS );
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherDESF::smooth( const real_t x )
{
    real_t lt_1;

    if ( init ) {
        lt = x;
        bt = x;
        init = false;
    }
    lt_1 = lt;
    lt = ( ( 1.0_re - alpha)*lt_1 ) + ( alpha*x ); /*level*/
    bt = ( ( 1.0_re - beta )*bt ) + ( beta*( lt - lt_1 ) ); /*trend*/

    return lt + ( n*bt ); /*model/forecast*/
}
/*============================================================================*/
bool smootherALNF::setup( const real_t a,
                          const real_t m,
                          real_t *window,
                          const size_t wsize )
{
    bool retValue = false;

    if ( ( nullptr != window ) && ( wsize > 0U ) && ( a > 0.0_re ) && ( a < 1.0_re ) && ( m > 0.0_re ) && ( m < 1.0_re ) ) {
        alpha = a;
        mu = m;
        xx = window;
        w = &window[ wsize ];
        w_1 = ( mu > 0.0_re ) ? &window[ 2U*wsize ] : nullptr;
        n = wsize;
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
real_t smootherALNF::smooth( const real_t x )
{
    real_t xe;

    if ( init ) {
        const real_t np = 1.0_re/static_cast<real_t>( n );

        windowSet( xx, n, x );
        windowSet( w, n, np );
        if ( nullptr != w_1 ) {
            windowSet( w_1, n, np );
        }
        init = false;
    }
    xe = discreteSystem::updateFIR( xx, n, x, w );
    if ( nullptr != w_1 ) {
        real_t * const ww_1 = &w[ n ];

        for ( size_t i = 0U ; i < n ; ++i ) {
            const real_t w0 = w[ i ];
            const real_t w1 = ww_1[ i ];
 
            w[ i ] += ( alpha*( x - xe )*xx[ i ] ) + ( mu*( w0 - w1 ) );
            w_1[ i ] = w0;
        }
    }
    else {
        for ( size_t i = 0U ; i < n ; ++i ) {
            w[ i ] += alpha*( x - xe )*xx[ i ];
        }
    }

    return xe;
}
/*============================================================================*/
