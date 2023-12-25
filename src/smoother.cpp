#include <include/smoother.hpp>
#include <include/ltisys.hpp>

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

    if ( ( a > 0.0 ) && ( a < 1.0 ) ) {
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

    if ( ( a > 0.0 ) && ( a < 1.0 ) ) {
        real_t aa, p1, r;
        aa = a*a;
        /*cstat -MISRAC2012-Dir-4.11_b*/
        p1 = sqrt( 2.0*a ); /*arg always positive*/
        /*cstat +MISRAC2012-Dir-4.11_b*/
        r = 1.0 + p1 + aa;
        k = aa/r;
        a1 = 2.0*( aa - 1.0 )/r;
        a2 = ( 1.0 - p1 + aa )/r;
        b1 = 2.0*k;
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

    if ( ( nullptr != window ) && ( w_size > 0U ) && ( a > 0.0 ) && ( a < 1.0 ) ) {
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
    if ( fabs( m - x ) > ( alpha*fabs( m ) ) ) {
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

    if ( ( nullptr != window ) && ( w_size > 0U ) && ( a > 0.0 ) && ( a < 1.0 ) ) {
        alpha = a;
        tdl::setup( window, w_size );
        sum = 0.0;
        m = 0.0;
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
    if ( fabs( m - x ) > ( alpha*fabs( m ) ) ) {
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
                          const size_t w_size )
{
    bool retValue = false;
    const size_t ws = wsize/2U;

    if ( ( nullptr != window ) && ( w_size > 0U ) && ( c < static_cast<real_t>( ws ) ) && ( sg > 0.0 ) ) {
        real_t * const kernel = &window[ ws ];
        real_t r, sum = 0.0;
        size_t i;
        real_t l, center;
        /*cstat -MISRAC++2008-5-0-7*/
        l = static_cast<real_t>( wsize - 1U )/2.0;
        /*cstat +MISRAC++2008-5-0-7*/
        center = c - l;
        r = 2.0*sg*sg;
        for ( i = 0U ; i < ws ; ++i ) {
            /*cstat -MISRAC++2008-5-0-7*/
            real_t d = static_cast<real_t>( i ) - l; /*symmetry*/
            /*cstat +MISRAC++2008-5-0-7*/
            d -= center;
            /*cstat -CERT-FLP32-C_b*/
            kernel[ i ] =  exp( -( d*d )/r );
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
    
    if ( ( lam > 0.0 ) && ( lam < 1.0 ) ) {
        lambda = lam;
        m = 0.0;
        w = 1.0;
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
        w = 1.0;
        init = false;
    }
    w = ( lambda*w ) + 1.0;
    iw = 1.0/w;
    m = ( m*( 1.0 - iw ) ) + ( iw*x );

    return m;
}
/*============================================================================*/
bool smootherKLMN::setup( const real_t processNoiseCov,
                          const real_t measureNoiseCov,
                          const real_t estErrorCov )
{
    bool retValue = false;

    if ( ( processNoiseCov > 0.0 ) && ( measureNoiseCov > 0.0 ) && ( estErrorCov > 0.0 ) ) {
        p = estErrorCov;
        q = processNoiseCov;
        r = measureNoiseCov;
        A = 1.0;
        H = 1.0;
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
    p = ( 1.0 - ( gain*H ) )*p; /*covariance update*/

    return xS;
}
/*============================================================================*/
bool smootherDESF::setup( const real_t a,
                          const real_t b,
                          const real_t N )
{
    bool retValue = false;

    if ( ( n >= 0.0 ) && ( a > 0.0 ) && ( a < 1.0 ) && ( b > 0.0 ) && ( b < 1.0 ) ) {
        alpha = a;
        beta = b;
        n = round( N );
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
    lt = ( ( 1.0 - alpha)*lt_1 ) + ( alpha*x ); /*level*/
    bt = ( ( 1.0 - beta )*bt ) + ( beta*( lt - lt_1 ) ); /*trend*/

    return lt + ( n*bt ); /*model/forecast*/
}
/*============================================================================*/
bool smootherALNF::setup( const real_t a,
                          const real_t m,
                          real_t *window,
                          const size_t wsize )
{
    bool retValue = false;

    if ( ( nullptr != window ) && ( wsize > 0U ) && ( a > 0.0 ) && ( a < 1.0 ) && ( m > 0.0 ) && ( m < 1.0 ) ) {
        alpha = a;
        mu = m;
        xx = window;
        w = &window[ wsize ];
        w_1 = ( mu > 0.0 ) ? &window[ 2U*wsize ] : nullptr;
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
        const real_t np = 1.0/static_cast<real_t>( n );

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
