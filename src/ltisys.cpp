#include <include/ltisys.hpp>

using namespace qlibs;

/*============================================================================*/
void ltisys::normalizeTransferFunction( real_t *num,
                                        real_t *den,
                                        size_t n_num,
                                        size_t n_den )
{
    const real_t a0 = den[ 0 ];
    size_t i;

    for ( i = 0U ; i < n_num ; ++i ) {
        num[ i ] /= a0;
    }
    for ( i = 0U ; i < n_den ; ++i ) {
        den[ i ] /= a0;
    }
    b0 = num[ 0 ];
}
/*============================================================================*/
real_t ltisys::saturate( real_t y )
{
    if ( y < min ) {
        y = min;
    }
    else if ( y > max ) {
        y = max;
    }
    else {
        /*do nothing*/
    }
    return y;
}
/*============================================================================*/
bool ltisys::setDelay( real_t * const w,
                       const size_t nD,
                       const real_t initVal ) noexcept
{
    bool retValue = false;

    if ( isInitialized() ) {
        tdl::setup( w, nD, initVal );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool ltisys::setSaturation( const real_t minV,
                            const real_t maxV ) noexcept
{
    bool retValue = false;

    if ( isInitialized() && ( maxV > minV ) ) {
        min = minV;
        max = maxV;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t ltisys::excite( real_t u )
{
    real_t y = 0.0_re;

    if ( isInitialized() ) {
        if ( tdl::isInitialized() ) {
            insertSample( u );
            u = getOldest();
        }
        y = saturate( update( u ) );
    }

    return y;
}
/*============================================================================*/
real_t discreteSystem::updateFIR( real_t *w,
                                  const size_t wsize,
                                  const real_t x,
                                  const real_t * const c )
{
    size_t i;
    real_t y = 0.0_re;

    if ( nullptr != c ) {
        for ( i = ( wsize - 1U ) ; i >= 1U ; --i ) {
            w[ i ] = w[ i - 1U ];
            y += w[ i ]*c[ i ];
        }
        y += c[ 0 ]*x;
    }
    else {
        for ( i = ( wsize - 1U ) ; i >= 1U ; --i ) {
            w[ i ] = w[ i - 1U ];
            y += w[ i ];
        }
        y += x;
    }
    w[ 0 ] = x;

    return y;
}
/*============================================================================*/
bool discreteSystem::setInitStates( const real_t *xi )
{
    bool retValue = false;

    if ( isInitialized() ) {
        if ( nullptr != xi ) {
            for ( size_t i = 0U; i < n ; ++i ) {
                xd[ i ] = xi[ i ];
            }
        }
        else {
            for ( size_t i = 0U; i < n ; ++i ) {
                xd[ i ] = 0.0_re;
            }
        }
       retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool discreteSystem::setup( real_t *num,
                            real_t *den,
                            real_t *x,
                            const size_t n_b,
                            const size_t n_a ) noexcept
{
    bool retValue = false;

    if ( ( nullptr != num ) && ( nullptr != den ) && ( nullptr != x ) && ( n_b > 0U ) ) {
        b = num;
        na = n_a - 1U;
        nb = n_b;
        n = ( na > nb ) ? na : nb;
        xd = x;
        a = &den[ 1 ];
        type = LTISYS_TYPE_DISCRETE;
        normalizeTransferFunction( num, den, n_a, n_b );
        (void)setInitStates();
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t discreteSystem::update( const real_t u )
{
    real_t v = u;

    /*using direct-form 2*/
    for ( size_t i = 0 ; i < na ; ++i ) {
        v -= a[ i ]*xd[ i ];
    }
    return updateFIR( xd, n, v, b );
}
/*============================================================================*/
bool continuousSystem::setup( real_t *num,
                              real_t *den,
                              nState *x,
                              const size_t nD,
                              const real_t dT ) noexcept
{
    bool retValue = false;

    if ( ( nullptr != num ) && ( nullptr != den ) && ( nullptr != x ) && ( nD >= 1U ) ) {
        b = &num[ 1 ];
        n = nD;
        nb = n;
        na = nD + 1U;
        xc = x;
        dt = dT;
        a = &den[ 1 ];
        type = LTISYS_TYPE_CONTINUOUS;
        normalizeTransferFunction( num, den, na, na );
        (void)setInitStates();
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool continuousSystem::setInitStates( const real_t *xi )
{
    bool retValue = false;

    if ( isInitialized() ) {
        if ( nullptr != xi ) {
            for ( size_t i = 0U; i < n ; ++i ) {
                xc[ i ].init( xi[ 0 ], xi[ 0 ], xi[ 0 ] );
            }
        }
        else {
            for ( size_t i = 0U; i < n ; ++i ) {
                xc[ i ].init();
            }
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
#if defined( LTISYS_EVAL_MODEL_CONTROLLABLE )
real_t continuousSystem::update( const real_t u )
{
    real_t y = 0.0_re;
    real_t dx0 = 0.0_re;

    if ( 1U == n ) {
        dx0 = ( u - ( xc[ 0 ]*a[ 0 ] ) );
        (void)xc[ 0 ].integrate( dx0, dt );
    }
    else {
        /*compute states of the system by using the controllable canonical form*/
        for ( size_t i = ( n - 1U ) ; i >= 1U ; --i ) {
            dx0 += a[ i ]*xc[ i ](); /*compute the first derivative*/
            /*integrate to obtain the remaining states*/
            (void)xc[ i ].integrate( xc[ i - 1U ](), dt );
            /*compute the first part of the output*/
            y += ( b[ i ] - ( a[ i ]*b0 ) )*xc[ i ];
        }
        /*compute remaining part of the output that depends of the first state*/
        dx0 = u - ( dx0 + ( a[ 0 ]*xc[ 0 ] ) );
        (void)xc[ 0 ].integrate( dx0, dt ); /*integrate to get the first state*/
        /*compute the remaining part of the output*/
    }
    /*get the output of the system*/
    y += ( b[ 0 ] - ( a[ 0 ]*b0 ) )*xc[ 0 ];

    return y;
}
#elif defined( LTISYS_EVAL_MODEL_OBSERVABLE )
/*============================================================================*/
real_t continuousSystem::update( const real_t u )
{
    real_t y = 0.0_re;
    const real_t x0 = xc[ 0 ](); /*save first state for computation*/
    const size_t N = n - 1U;
    /*compute the last derivative dx_{n}*/
    const real_t dxn = ( -a[ N ]*x0 ) + ( ( b[ N ] - a[ N ]*b0 )*u );

    if ( 1U == n ) {
        xc[ 0 ].integrate( dxn, dt );
    }
    else {
        /*compute states of the system by using the observable canonical form*/
        for ( size_t i = 0U; i < N ; ++i ) {
            const real_t dxi = ( -a[ i ]*x0 ) + xc[ i + 1U ] + ( ( b[ i ] - a[ i ]*b0 )*u );
            (void)xc[ i ].integrate( dxi , dt );
        }
        /*integrate the remaining last state*/
        (void)xc[ N ].integrate( dxn, dt );
    }
    /*get the output of the system*/
    y = xc[ 0 ] + ( b0*u );

    return y;
}
#else
    #error "LTISYS evaluation mode not defined"
#endif
/*============================================================================*/
bool continuousSystem::setIntegrationMethod( integrationMethod m )
{
    bool retValue = false;

    if ( isInitialized() ) {
        /*cstat -MISRAC++2008-0-1-2_a*/
        if ( ( INTEGRATION_RECTANGULAR == m ) || ( INTEGRATION_TRAPEZOIDAL == m ) || ( INTEGRATION_SIMPSON == m ) ) {
        /*cstat +MISRAC++2008-0-1-2_a*/
            for ( size_t i = 0U; i < n ; ++i ) {
                xc[ i ].setIntegrationMethod( m );
            }
            retValue = true;
        }
    }

    return retValue;
}
/*============================================================================*/