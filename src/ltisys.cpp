#include <include/ltisys.hpp>

using namespace qlibs;

/*============================================================================*/
void ltisys::normalizeTransferFunction( real_t *num, real_t *den, size_t n_num, size_t n_den )
{
    const real_t a0 = den[ 0 ];
    size_t i;

    for ( i = 0 ; i < n_num ; ++i ) {
        num[ i ] /= a0;
    }
    for ( i = 0 ; i < n_den ; ++i ) {
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
bool ltisys::setDelay( real_t * const w, const size_t nD, const real_t initVal )
{
    bool retValue = false;

    if ( isInitialized() ) {
        tdl::setup( w, nD, initVal );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool ltisys::setSaturation( const real_t minV, const real_t maxV ) 
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
real_t discreteSystem::updateFIR( real_t *w, const size_t wsize, const real_t x, const real_t * const c )
{
    size_t i;
    real_t y = 0.0f;

    if ( nullptr != c ) {
        for ( i = ( wsize - 1u ) ; i >= 1u ; --i ) {
            w[ i ] = w[ i - 1u ];
            y += w[ i ]*c[ i ];
        }
        y += c[ 0 ]*x;
    }
    else {
        for ( i = ( wsize - 1u ) ; i >= 1u ; --i ) {
            w[ i ] = w[ i - 1u ];
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
        for ( size_t i = 0u; i < n ; ++i ) {
            const real_t zero = 0.0f;
            const real_t * const iv = ( nullptr != xi ) ? &xi[ i ] : &zero;
            xd[ i ] = iv[ i ];
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool discreteSystem::setup( real_t *num, real_t *den, real_t *x, const size_t n_b, const size_t n_a )
{
    bool retValue = false;

    if ( ( nullptr != num ) && ( nullptr != den ) && ( nullptr != x ) && ( n_a > 0U ) ) {
        b = num;
        na = n_a - 1;
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
real_t discreteSystem::excite( real_t u )
{
    real_t y = 0.0;

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
bool continuousSystem::setup( real_t *num, real_t *den, state *x,const size_t n_a, const real_t dT )
{
    bool retValue = false;

    if ( ( nullptr != num ) && ( nullptr != den ) && ( nullptr != x ) && ( n_a > 0U ) ) {
        b = &num[ 1 ];
        n = n_a - 1u;
        nb = n;
        na = n_a;
        xc = x;
        dt = dT;
        a = &den[ 1 ];
        type = LTISYS_TYPE_CONTINUOUS;
        normalizeTransferFunction( num, den, n_a, n_a );
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
        for ( size_t i = 0u; i < n ; ++i ) {
            const real_t zero = 0.0f;
            const real_t * const iv = ( nullptr != xi ) ? &xi[ i ] : &zero;
            xc[ i ].init( iv[ 0 ], iv[ 0 ], iv[ 0 ] );
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t continuousSystem::update( const real_t u )
{
    real_t y = 0.0f;
    real_t dx0 = 0.0f;

    if ( 1u == n ) {
        dx0 = ( u - ( xc[ 0 ]*a[ 0 ] ) );
        (void)xc[ 0 ].integrate( dx0, dt );
        y = ( b[ 0 ] - ( a[ 0 ]*b0 ) )*xc[ 0 ];
    }
    else {
        /*compute states of the system by using the controllable canonical form*/
        for ( size_t i = ( n - 1u ) ; i >= 1u ; --i ) {
            dx0 += a[ i ]*xc[ i ](); /*compute the first derivative*/
            /*integrate to obtain the remaining states*/
            (void)xc[ i ].integrate( xc[ i - 1u ](), dt );
            /*compute the first part of the output*/
            y += ( b[ i ] - ( a[ i ]*b0 ) )*xc[ i ];
        }
        /*compute remaining part of the output that depends of the first state*/
        dx0 = u - ( dx0 + ( a[ 0 ]*xc[ 0 ] ) );
        (void)xc[ 0 ].integrate( dx0, dt ); /*integrate to get the first state*/
        /*compute the remaining part of the output*/
        y += ( b[ 0 ] - ( a[ 0 ]*b0 ) )*xc[ 0 ];
    }

    return y;
}
/*============================================================================*/
real_t continuousSystem::excite( real_t u )
{
    real_t y = 0.0;

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
bool continuousSystem::setIntegrationMethod( integrationMethod m )
{
    bool retValue = false;

    if ( isInitialized() ) {
        /*cstat -MISRAC++2008-0-1-2_a*/
        if ( ( INTEGRATION_RECTANGULAR == m ) || ( INTEGRATION_TRAPEZOIDAL == m ) || ( INTEGRATION_SIMPSON == m ) ) {
        /*cstat +MISRAC++2008-0-1-2_a*/
            for ( size_t i = 0; i < n ; ++i ) {
                xc[ i ].setIntegrationMethod( m );
            }
            retValue = true;
        }
    }

    return retValue;
}
/*============================================================================*/