#include <include/fis.hpp>
#include <include/mathex.hpp>
#include <include/ffmath.hpp>

using namespace qlibs;

/*cstat -CERT-FLP32-C_b*/

/*! @cond  */
/*============================================================================*/
real_t fis::core::TriMF( const fis::ioBase * const in,
                         const real_t *p,
                         const size_t n )
{
    real_t a, b, c, tmp;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    c = p[ 2 ];
    tmp = Min( ( x - a )/( b - a ) , ( c - x )/( c - b ) );

    return Max( tmp , 0.0_re );
}
/*============================================================================*/
real_t fis::core::TrapMF( const fis::ioBase * const in,
                          const real_t *p,
                          const size_t n )
{
    real_t a, b, c, d, tmp;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    c = p[ 2 ];
    d = p[ 3 ];
    tmp = Min( ( x - a )/( b - a ) , 1.0_re );
    tmp = Min( tmp, ( d - x )/( d - c ) ) ;

    return Max( tmp , 0.0_re );
}
/*============================================================================*/
real_t fis::core::GBellMF( const fis::ioBase * const in,
                           const real_t *p,
                           const size_t n )
{
    real_t a, b, c;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    c = p[ 2 ];

    return ( 1.0_re/( 1.0_re + ffmath::pow( ffmath::absf( ( x - c )/a ) , 2.0_re*b ) ) );
}
/*============================================================================*/
real_t fis::core::GaussMF( const fis::ioBase * const in,
                           const real_t *p,
                           const size_t n )
{
    real_t a, c, tmp;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    c = p[ 1 ];
    tmp = ( x - c )/a;

    return ffmath::exp( -0.5_re*tmp*tmp );
}
/*============================================================================*/
real_t fis::core::Gauss2MF( const fis::ioBase * const in,
                            const real_t *p,
                            const size_t n )
{
    real_t c1, c2, f1, f2;
    const real_t x = in[ 0 ].value;

    c1 = p[ 1 ];
    c2 = p[ 3 ];
    f1 = ( x <= c1 ) ? GaussMF( in , p, n ) : 1.0_re;
    f2 = ( x <= c2 ) ? GaussMF( in , &p[ 2 ], n ) : 1.0_re;

    return f1*f2;
}
/*============================================================================*/
real_t fis::core::SigMF( const fis::ioBase * const in,
                         const real_t *p,
                         const size_t n )
{
    real_t a, b;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];

    return 1.0_re/( 1.0_re + ffmath::exp( -a*( x - b ) ) );
}
/*============================================================================*/
real_t fis::core::TSigMF( const fis::ioBase * const in,
                          const real_t *p,
                          const size_t n )
{
    real_t a, b, y;
    const real_t x = in[ 0 ].value;
    const real_t min = in[ 0 ].min;
    const real_t max = in[ 0 ].max;
    (void)n;

    a = p[ 0 ]; /*slope*/
    b = p[ 1 ]; /*inflection*/
    if ( ffmath::isEqual( x, 1.0_re ) ) {
        if ( a >= 0.0_re ) {
            y = max;
        }
        else {
            y = min;
        }
    }
    else if ( ffmath::isEqual( x, 0.0_re ) ) {
        if ( a >= 0.0_re ) {
            y = min;
        }
        else {
            y = max;
        }
    }
    else {
        y = b - ( ffmath::log( ( 1.0_re/x ) - 1.0_re )/a );
    }

    return y;
}
/*============================================================================*/
real_t fis::core::DSigMF( const fis::ioBase * const in,
                          const real_t *p,
                          const size_t n )
{
    return ffmath::absf( SigMF( in , p, n ) - SigMF( in , &p[ 2 ], n ) );
}
/*============================================================================*/
real_t fis::core::PSigMF( const fis::ioBase * const in,
                          const real_t *p,
                          const size_t n )
{
    return ffmath::absf( SigMF( in , p, n )*SigMF( in , &p[ 2 ], n ) );
}
/*============================================================================*/
real_t fis::core::SMF( const fis::ioBase * const in,
                       const real_t *p,
                       const size_t n )
{
    real_t a, b, tmp, y;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    if ( x <= a ) { // skipcq: CXX-W2041
        y =  0.0_re;
    }
    else if ( x >= b ) {
        y = 1.0_re;
    }
    /*cppcheck-suppress knownConditionTrueFalse */
    else if ( ( x >= a ) && ( x <= ( ( a + b )*0.5_re ) ) ) {
        tmp = ( x - a )/( b - a );
        y = 2.0_re*tmp*tmp;
    }
    /*cppcheck-suppress knownConditionTrueFalse */
    else if ( ( x <= b ) && ( x >= ( ( a + b )*0.5_re ) ) ) {
        tmp = ( x - b )/( b - a );
        y = ( 1.0_re - ( 2.0_re*tmp*tmp ) );
    }
    else {
        y = 0.0_re;
    }

    return y;
}
/*============================================================================*/
real_t fis::core::TSMF( const fis::ioBase * const in,
                        const real_t *p,
                        const size_t n )
{
    real_t diff, a, b, ta, tb, ma, mb;
    const real_t x = in[ 0 ].value;
    (void)n;
    fis::ioBase tmp;

    a = p[ 0 ]; /*start*/
    b = p[ 1 ]; /*end*/
    diff = b - a;
    diff = 0.5_re*diff*diff;
    ta = a + ffmath::sqrt( x*diff );
    tmp.value = ta;
    ma = SMF( &tmp, p, n );
    tb = b + ffmath::sqrt( -( x - 1.0_re )*diff );
    tmp.value = tb;
    mb = SMF( &tmp, p, n );
    return  ( ffmath::absf( x - ma ) < ffmath::absf( x - mb ) ) ? ta : tb;
}
/*============================================================================*/
real_t fis::core::ZMF( const fis::ioBase * const in,
                       const real_t *p,
                       const size_t n )
{
    real_t a, b, tmp, y;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    if ( x <= a ) {
        y = 1.0_re;
    }
    else if ( x >= b ) { // skipcq: CXX-W2041
        y = 0.0_re;
    }
    /*cppcheck-suppress knownConditionTrueFalse */
    else if ( ( x >= a ) && ( x <= ( ( a + b )*0.5_re ) ) ) {
        tmp = ( x - a )/( b - a );
        y = 1.0_re - ( 2.0_re*tmp*tmp );
    }
    /*cppcheck-suppress knownConditionTrueFalse */
    else if ( ( x <= b ) && ( x >= ( ( a + b )*0.5_re ) ) ) {
        tmp = ( x - b )/( b - a );
        y = 2.0_re*tmp*tmp;
    }
    else {
        y = 0.0_re;
    }

    return y;
}
/*============================================================================*/
real_t fis::core::LinSMF( const fis::ioBase * const in,
                          const real_t *p,
                          const size_t n )
{
    real_t a, b, y;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    if ( a < b ) {
        if ( x < a ) {
            y = 0.0_re;
        }
        else if ( x > b ) {
            y = 1.0_re;
        }
        else {
            y = ( x - a )/( b - a );
        }
    }
    else if ( ffmath::isEqual( a, b ) ) {
        y = ( x < a ) ? 0.0_re : 1.0_re;
    }
    else {
        y = 0.0_re;
    }

    return y;
}
/*============================================================================*/
real_t fis::core::LinZMF( const fis::ioBase * const in,
                          const real_t *p,
                          const size_t n )
{
    real_t a, b, y;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    if ( a < b ) {
        if ( x < a ) {
            y = 1.0_re;
        }
        else if ( x > b ) {
            y = 0.0_re;
        }
        else {
            y = ( a - x )/( a - b );
        }
    }
    else if ( ffmath::isEqual( a, b ) ) {
        y = ( x < a ) ? 1.0_re : 0.0_re;
    }
    else {
        y = 0.0_re;
    }

    return y;
}
/*============================================================================*/
real_t fis::core::TZMF( const fis::ioBase * const in,
                        const real_t *p,
                        const size_t n )
{
    real_t diff, a, b, ta, tb, ma, mb;
    const real_t x = in[ 0 ].value;
    (void)n;
    fis::ioBase tmp;

    a = p[ 0 ]; /*start*/
    b = p[ 1 ]; /*end*/
    diff = b - a;
    diff = 0.5_re*diff*diff;
    ta = a + ffmath::sqrt( -( x - 1.0_re )*diff );
    tmp.value = ta;
    ma = SMF( &tmp, p, n );
    tb = b + ffmath::sqrt( x*diff );
    tmp.value = tb;
    mb = SMF( &tmp, p, n );
    return  ( ffmath::absf( x - ma ) < ffmath::absf( x - mb ) ) ? ta : tb;
}
/*============================================================================*/
real_t fis::core::PiMF( const fis::ioBase * const in,
                        const real_t *p,
                        const size_t n )
{
    return ffmath::absf( SMF( in , p, n )*ZMF( in , &p[ 2 ], n ) );
}
/*============================================================================*/
real_t fis::core::SingletonMF( const fis::ioBase * const in,
                               const real_t *p,
                               const size_t n )
{
    const real_t x = in[ 0 ].value;
    (void)n;

    return ( ffmath::isEqual( x, p[ 0 ] ) ) ? 1.0_re : 0.0_re;
}
/*============================================================================*/
real_t fis::core::ConcaveMF( const fis::ioBase * const in,
                             const real_t *p,
                             const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t i, e, y;
    (void)n;

    i = p[ 0 ];
    e = p[ 1 ];
    if ( ( i <= e ) && ( x < e ) ) {
        y = ( e - i )/( ( 2.0_re*e ) - i -x );
    }
    else if ( ( i > e ) && ( x > e ) ) {
        y = ( i - e )/( -( 2.0_re*e ) + i +x );
    }
    else {
        y = 1.0_re;
    }

    return y;
}
/*============================================================================*/
real_t fis::core::TConcaveMF( const fis::ioBase * const in,
                              const real_t *p,
                              const size_t n )
{
    real_t i, e;

    i = p[ 0 ];
    e = p[ 1 ];

    return ( ( i - e )/ConcaveMF( in, p, n ) ) + ( 2.0_re*e ) - i;
}
/*============================================================================*/
real_t fis::core::SpikeMF( const fis::ioBase * const in,
                           const real_t *p,
                           const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t w, c;
    (void)n;

    w = p[ 0 ];
    c = p[ 1 ];

    return ffmath::exp( -ffmath::absf( 10.0_re*( x - c )/w ) );
}
/*============================================================================*/
real_t fis::core::TLinSMF( const fis::ioBase * const in,
                           const real_t *p,
                           const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t a, b;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];

    return ( ( b - a )*x ) + a;
}
/*============================================================================*/
real_t fis::core::TLinZMF( const fis::ioBase * const in,
                           const real_t *p,
                           const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t a, b;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];

    return a - ( ( a - b )*x );
}
/*============================================================================*/
real_t fis::core::TRampMF( const ioBase * const in,
                           const real_t *p,
                           const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t a, b;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];

    return ( b - a )*x + a;
}
/*============================================================================*/
real_t fis::core::RectangleMF( const fis::ioBase * const in,
                               const real_t *p,
                               const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t s, e;
    (void)n;

    s = p[ 0 ];
    e = p[ 1 ];

    return ( ( x >= s ) && ( x <= e ) ) ? 1.0_re : 0.0_re;
}
/*============================================================================*/
real_t fis::core::CosineMF( const fis::ioBase * const in,
                            const real_t *p,
                            const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t c, w, y;
    constexpr real_t pi = 3.14159265358979323846_re;
    (void)n;

    c = p[ 0 ];
    w = p[ 1 ];
    if ( ( x < ( c - ( 0.5_re*w ) ) ) || ( x > ( c + ( 0.5_re*w ) ) ) ) {
        y = 0.0_re;
    }
    else {
        y = 0.5_re*( 1.0_re + ffmath::cos( 2.0_re/w*pi*( x - c) ) );
    }

    return y;
}
/*============================================================================*/
real_t fis::core::ConstantMF( const fis::ioBase * const in,
                              const real_t *p,
                              const size_t n )
{
    (void)in;
    (void)n;
    return p[ 0 ];
}
/*============================================================================*/
real_t fis::core::LinearMF( const fis::ioBase * const in,
                            const real_t *p,
                            const size_t n )
{
    real_t px = 0.0_re;
    size_t i;

    for ( i = 0U ; i < n ; ++i ) {
        px += in[ i ].value*p[ i ];
    }
    px += p[ i ];

    return px;
}
/*============================================================================*/
real_t fis::core::bound( real_t y,
                         const real_t minValue,
                         const real_t maxValue )
{
    (void)inRangeCoerce( y, minValue, maxValue );

    return y;
}
/*============================================================================*/
real_t fis::core::Min( const real_t a,
                       const real_t b )
{
    return bound( ( a < b ) ? a : b );
}
/*============================================================================*/
real_t fis::core::Prod( const real_t a,
                        const real_t b )
{
    return bound( a*b );
}
/*============================================================================*/
real_t fis::core::Max( const real_t a,
                       const real_t b )
{
    return bound( ( a > b ) ? a : b );
}
/*============================================================================*/
real_t fis::core::ProbOr( const real_t a,
                          const real_t b )
{
    return bound( a + b - ( a*b ) );
}
/*============================================================================*/
real_t fis::core::Sum( const real_t a,
                       const real_t b )
{
    return bound( a + b );
}
/*============================================================================*/
real_t fis::core::deFuzzCentroid( fis::output * const o,
                                  const fis::deFuzzState stage )
{
    real_t d = 0.0_re;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            o->v[ sum_xy ] += o->x*o->y;
            o->v[ sum_y ] += o->y;
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ sum_xy ] = 0.0_re; /*store sum(x*y)*/
            o->v[ sum_y ]= 0.0_re; /*store sum(y)*/
            break;
        case FIS_DEFUZZ_END:
            d = o->v[ sum_xy ]/o->v[ sum_y ]; /*sum(x*y)/sum(y)*/
            break;
        default:
            break;
    }

    return d;
}
/*============================================================================*/
real_t fis::core::deFuzzBisector( fis::output * const o,
                                  const fis::deFuzzState stage )
{
    size_t k;
    real_t d = 0.0_re;
    fis::instance *f;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            o->v[ sum_y ] += o->y;
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ sum_y ] = 0.0_re; /*store sum(y)*/
            break;
        case FIS_DEFUZZ_END:
            o->v[ currentArea ] = 0.0_re;
            o->v[ halfArea ] = 0.5_re*o->v[ sum_y ];
            f = o->owner;
            for ( k = 0U ; k < f->nPoints ; ++k ) {
                o->y = 0.0_re;
                o->getNextX( k );
                o->value = o->x;
                f->fuzzyAggregate();
                o->v[ currentArea ] += o->y;
                if ( o->v[ currentArea ] >= o->v[ halfArea ] ) {
                    break;
                }
            }
            d = o->x;
            break;
        default:
            break;
    }

    return d;
}
/*============================================================================*/
real_t fis::core::deFuzzLOM( fis::output * const o,
                             const fis::deFuzzState stage )
{
    real_t d = 0.0_re;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            if ( o->y >= o->v[ yMax ] ) {
                o->v[ yMax ] = o->y;
                o->v[ xLargest ] = o->x;
            }
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ yMax ] = -1.0_re;
            o->v[ xLargest ] = o->max;
            break;
        case FIS_DEFUZZ_END:
            d = o->v[ xLargest ];
            break;
        default:
            break;
    }

    return d;
}
/*============================================================================*/
real_t fis::core::deFuzzSOM( fis::output * const o,
                             const fis::deFuzzState stage )
{
    real_t d = 0.0_re;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            if ( o->y >= o->v[ yMax ] ) {
                o->v[ yMax ] = o->y;
                o->v[ xLargest ] = o->x;
            }
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ yMax ] = -1.0_re;
            o->v[ xLargest ] = o->min;
            break;
        case FIS_DEFUZZ_END:
            d = o->v[ xLargest ];
            break;
        default:
            break;
    }

    return d;
}
/*============================================================================*/
real_t fis::core::deFuzzMOM( fis::output * const o,
                             const fis::deFuzzState stage )
{
    real_t d = 0.0_re;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            if ( o->y > o->v[ yMax ] ) {
                o->v[ yMax ] = o->y;
                o->v[ xSmallest ] = o->x;
                o->v[ xLargest ] = o->x;
                o->v[ sp ] = 1.0_re;
            }
            else if ( ffmath::isEqual( o->y , o->v[ yMax ] ) && ( o->v[ sp ]  > 0.0_re ) ) {
                o->v[ xLargest ] = o->x;
            }
            else if ( o->y < o->v[ yMax ] ) {
                o->v[ sp ] = -1.0_re;
            }
            else {
                /*nothing to do*/
            }
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ yMax ] = -1.0_re;
            o->v[ xSmallest ] = o->min;
            o->v[ xLargest ]= o->max;
            o->v[ sp ] = -1.0_re;
            break;
        case FIS_DEFUZZ_END:
            d = 0.5_re*( o->v[ xSmallest ] + o->v[ xLargest ] );
            break;
        default:
            break;
    }

    return d;
}
/*============================================================================*/
real_t fis::core::deFuzzWtAverage( fis::output * const o,
                                   const fis::deFuzzState stage )
{
    real_t d = 0.0_re;

    if ( FIS_DEFUZZ_END == stage ) {
        d = o->v[ sum_wz ]/o->v[ sum_w ];
    }

    return d;
}
/*============================================================================*/
real_t fis::core::deFuzzWtSum( fis::output * const o,
                               const fis::deFuzzState stage )
{
    real_t d = 0.0_re;

    if ( FIS_DEFUZZ_END == stage ) {
        d = o->v[ sum_wz ];
    }

    return d;
}
/*============================================================================*/
/*! @endcond  */

/*cstat +CERT-FLP32-C_b*/