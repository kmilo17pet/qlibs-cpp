#include "include/fis.hpp"
#include "include/mathex.hpp"
#include <cmath>

using namespace qlibs;

/*cstat -CERT-FLP32-C_b*/

/*============================================================================*/
real_t fisCore::TriMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b, c, tmp;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    c = p[ 2 ];
    tmp = Min( ( x - a )/( b - a ) , ( c - x )/( c - b ) );

    return Max( tmp , 0.0 );
}
/*============================================================================*/
real_t fisCore::TrapMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b, c, d, tmp;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    c = p[ 2 ];
    d = p[ 3 ];
    tmp = Min( ( x - a )/( b - a ) , 1.0 );
    tmp = Min( tmp, ( d - x )/( d - c ) ) ;

    return Max( tmp , 0.0 );
}
/*============================================================================*/
real_t fisCore::GBellMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b, c;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    c = p[ 2 ];
    
    return ( 1.0/( 1.0 + pow( fabs( ( x - c )/a ) , 2.0*b ) ) );
}
/*============================================================================*/
real_t fisCore::GaussMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, c, tmp;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    c = p[ 1 ];
    tmp = ( x - c )/a;

    return exp( -0.5*tmp*tmp );
}
/*============================================================================*/
real_t fisCore::Gauss2MF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t c1, c2, f1, f2;
    const real_t x = in[ 0 ].value;

    c1 = p[ 1 ];
    c2 = p[ 3 ];
    f1 = ( x <= c1 ) ? GaussMF( in , p, n ) : 1.0;
    f2 = ( x <= c2 ) ? GaussMF( in , &p[ 2 ], n ) : 1.0;

    return f1*f2;
}
/*============================================================================*/
real_t fisCore::SigMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];

    return 1.0/( 1.0 + exp( -a*( x - b ) ) );
}
/*============================================================================*/
real_t fisCore::TSigMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b, y;
    const real_t x = in[ 0 ].value;
    const real_t min = in[ 0 ].min;
    const real_t max = in[ 0 ].max;
    (void)n;

    a = p[ 0 ]; /*slope*/
    b = p[ 1 ]; /*inflection*/
    if ( isEqual( x, 1.0 ) ) {
        if ( a >= 0.0 ) {
            y = max;
        }
        else {
            y = min;
        }
    }
    else if ( isEqual( x, 0.0 ) ) {
        if ( a >= 0.0 ) {
            y = min;
        }
        else {
            y = max;
        }
    }
    else {
        y = b - ( log( ( 1.0/x ) - 1.0 )/a );
    }

    return y;
}
/*============================================================================*/
real_t fisCore::DSigMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    return fabs( SigMF( in , p, n ) - SigMF( in , &p[ 2 ], n ) );
}
/*============================================================================*/
real_t fisCore::PSigMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    return fabs( SigMF( in , p, n )*SigMF( in , &p[ 2 ], n ) );
}
/*============================================================================*/
real_t fisCore::SMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b, tmp, y;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    if ( x <= a ) { // skipcq: CXX-W2041
        y =  0.0;
    }
    else if ( x >= b ) {
        y = 1.0;
    }
    else if ( ( x >= a ) && ( x <= ( ( a + b )*0.5 ) ) ) {
        tmp = ( x - a )/( b - a );
        y = 2.0*tmp*tmp;
    }
    else if ( ( x <= b ) && ( x >= ( ( a + b )*0.5 ) ) ) {
        tmp = ( x - b )/( b - a );
        y = ( 1.0 - ( 2.0*tmp*tmp ) );
    }
    else {
        y = 0.0;
    }

    return y;
}
/*============================================================================*/
real_t fisCore::TSMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t diff, a, b, ta, tb, ma, mb;
    const real_t x = in[ 0 ].value;
    (void)n;
    fisIOBase tmp;

    a = p[ 0 ]; /*start*/
    b = p[ 1 ]; /*end*/
    diff = b - a;
    diff = 0.5*diff*diff;
    ta = a + sqrt( x*diff );
    tmp.value = ta;
    ma = SMF( &tmp, p, n );
    tb = b + sqrt( -( x - 1.0 )*diff );
    tmp.value = tb;
    mb = SMF( &tmp, p, n );
    return  ( fabs( x - ma ) < fabs( x - mb ) ) ? ta : tb;
}
/*============================================================================*/
real_t fisCore::ZMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b, tmp, y;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    if ( x <= a ) {
        y = 1.0;
    }
    else if ( x >= b ) { // skipcq: CXX-W2041
        y = 0.0;
    }
    else if ( ( x >= a ) && ( x <= ( ( a + b )*0.5 ) ) ) {
        tmp = ( x - a )/( b - a );
        y = 1.0 - ( 2.0*tmp*tmp );
    }
    else if ( ( x <= b ) && ( x >= ( ( a + b )*0.5 ) ) ) {
        tmp = ( x - b )/( b - a );
        y = 2.0*tmp*tmp;
    }
    else {
        y = 0.0;
    }

    return y;
}
/*============================================================================*/
real_t fisCore::LinSMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b, y;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    if ( a < b ) {
        if ( x < a ) {
            y = 0.0;
        }
        else if ( x > b ) {
            y = 1.0;
        }
        else {
            y = ( x - a )/( b - a );
        }
    }
    else if ( isEqual( a, b ) ) {
        y = ( x < a ) ? 0.0 : 1.0;
    }
    else {
        y = 0.0;
    }

    return y;
}
/*============================================================================*/
real_t fisCore::LinZMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t a, b, y;
    const real_t x = in[ 0 ].value;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];
    if ( a < b ) {
        if ( x < a ) {
            y = 1.0;
        }
        else if ( x > b ) {
            y = 0.0;
        }
        else {
            y = ( a - x )/( a - b );
        }
    }
    else if ( isEqual( a, b ) ) {
        y = ( x < a ) ? 1.0 : 0.0;
    }
    else {
        y = 0.0;
    }

    return y;
}
/*============================================================================*/
real_t fisCore::TZMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t diff, a, b, ta, tb, ma, mb;
    const real_t x = in[ 0 ].value;
    (void)n;
    fisIOBase tmp;

    a = p[ 0 ]; /*start*/
    b = p[ 1 ]; /*end*/
    diff = b - a;
    diff = 0.5*diff*diff;
    ta = a + sqrt( -( x - 1.0 )*diff );
    tmp.value = ta;
    ma = SMF( &tmp, p, n );
    tb = b + sqrt( x*diff );
    tmp.value = tb;
    mb = SMF( &tmp, p, n );
    return  ( fabs( x - ma ) < fabs( x - mb ) ) ? ta : tb;
}
/*============================================================================*/
real_t fisCore::PiMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    return fabs( SMF( in , p, n )*ZMF( in , &p[ 2 ], n ) );
}
/*============================================================================*/
real_t fisCore::SingletonMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    const real_t x = in[ 0 ].value;
    (void)n;

    return ( isEqual( x, p[ 0 ] ) ) ? 1.0 : 0.0;
}
/*============================================================================*/
real_t fisCore::ConcaveMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t i, e, y;
    (void)n;

    i = p[ 0 ];
    e = p[ 1 ];
    if ( ( i <= e ) && ( x < e ) ) {
        y = ( e - i )/( ( 2.0*e ) - i -x );
    }
    else if ( ( i > e ) && ( x > e ) ) {
        y = ( i - e )/( -( 2.0*e ) + i +x );
    }
    else {
        y = 1.0;
    }

    return y;
}
/*============================================================================*/
real_t fisCore::TConcaveMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t i, e;

    i = p[ 0 ];
    e = p[ 1 ];

    return ( ( i - e )/ConcaveMF( in, p, n ) ) + ( 2.0*e ) - i;
}
/*============================================================================*/
real_t fisCore::SpikeMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t w, c;
    (void)n;

    w = p[ 0 ];
    c = p[ 1 ];

    return exp( -fabs( 10.0*( x - c )/w ) );
}
/*============================================================================*/
real_t fisCore::TLinSMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t a, b;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];

    return ( ( b - a )*x ) + a;
}
/*============================================================================*/
real_t fisCore::TLinZMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t a, b;
    (void)n;

    a = p[ 0 ];
    b = p[ 1 ];

    return a - ( ( a - b )*x );
}
/*============================================================================*/
real_t fisCore::RectangleMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t s, e;
    (void)n;

    s = p[ 0 ];
    e = p[ 1 ];

    return ( ( x >= s ) && ( x <= e ) ) ? 1.0 : 0.0;
}
/*============================================================================*/
real_t fisCore::CosineMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    const real_t x = in[ 0 ].value;
    real_t c, w, y;
    constexpr real_t pi = 3.14159265358979323846;
    (void)n;

    c = p[ 0 ];
    w = p[ 1 ];
    if ( ( x < ( c - ( 0.5*w ) ) ) || ( x > ( c + ( 0.5*w ) ) ) ) {
        y = 0.0;
    }
    else {
        y = 0.5*( 1.0 + cos( 2.0/w*pi*( x - c) ) );
    }

    return y;
}
/*============================================================================*/
real_t fisCore::ConstantMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    (void)in;
    (void)n;
    return p[ 0 ];
}
/*============================================================================*/
real_t fisCore::LinearMF( const fisIOBase * const in, const real_t *p, const size_t n )
{
    real_t px = 0.0;
    size_t i;

    for ( i = 0U ; i < n ; ++i ) {
        px += in[ i ].value*p[ i ];
    }
    px += p[ i ];

    return px;
}
/*============================================================================*/
real_t fisCore::bound( real_t y, const real_t minValue, const real_t maxValue )
{
    (void)inRangeCoerce( y, minValue, maxValue );

    return y;
}
/*============================================================================*/
real_t fisCore::Min( const real_t a, const real_t b )
{
    return bound( ( a < b ) ? a : b, 0.0, 1.0 );
}
/*============================================================================*/
real_t fisCore::Prod( const real_t a, const real_t b )
{
    return bound( a*b, 0.0, 1.0 );
}
/*============================================================================*/
real_t fisCore::Max( const real_t a, const real_t b )
{
    return bound( ( a > b ) ? a : b, 0.0, 1.0 );
}
/*============================================================================*/
real_t fisCore::ProbOr( const real_t a, const real_t b )
{
    return bound( a + b - ( a*b ), 0.0, 1.0 );
}
/*============================================================================*/
real_t fisCore::Sum( const real_t a, const real_t b )
{
    return bound( a + b, 0.0, 1.0 );
}
/*============================================================================*/
real_t fisCore::deFuzzCentroid( fisOutput * const o, const fisDeFuzzState stage )
{
    real_t d = 0.0;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            o->v[ sum_xy ] += o->x*o->y;
            o->v[ sum_y ] += o->y;
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ sum_xy ] = 0.0; /*store sum(x*y)*/
            o->v[ sum_y ]= 0.0; /*store sum(y)*/
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
real_t fisCore::deFuzzBisector( fisOutput * const o, const fisDeFuzzState stage )
{
    size_t k;
    real_t d = 0.0;
    fis *f;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            o->v[ sum_y ] += o->y;
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ sum_y ] = 0.0; /*store sum(y)*/
            break;
        case FIS_DEFUZZ_END:
            o->v[ currentArea ] = 0.0;
            o->v[ halfArea ] = 0.5*o->v[ sum_y ];
            f = o->owner;
            for ( k = 0U ; k < f->nPoints ; ++k ) {
                o->y = 0.0;
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
real_t fisCore::deFuzzLOM( fisOutput * const o, const fisDeFuzzState stage )
{
    real_t d = 0.0;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            if ( o->y >= o->v[ yMax ] ) {
                o->v[ yMax ] = o->y;
                o->v[ xLargest ] = o->x;
            }
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ yMax ] = -1.0;
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
real_t fisCore::deFuzzSOM( fisOutput * const o, const fisDeFuzzState stage )
{
    real_t d = 0.0;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            if ( o->y >= o->v[ yMax ] ) {
                o->v[ yMax ] = o->y;
                o->v[ xLargest ] = o->x;
            }
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ yMax ] = -1.0;
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
real_t fisCore::deFuzzMOM( fisOutput * const o, const fisDeFuzzState stage )
{
    real_t d = 0.0;

    switch ( stage ) {
        case FIS_DEFUZZ_COMPUTE:
            if ( o->y > o->v[ yMax ] ) {
                o->v[ yMax ] = o->y;
                o->v[ xSmallest ] = o->x;
                o->v[ xLargest ] = o->x;
                o->v[ sp ] = 1.0;
            }
            else if ( isEqual( o->y , o->v[ yMax ] ) && ( o->v[ sp ]  > 0.0 ) ) {
                o->v[ xLargest ] = o->x;
            }
            else if ( o->y < o->v[ yMax ] ) {
                o->v[ sp ] = -1.0;
            }
            else {
                /*nothing to do*/
            }
            break;
        case FIS_DEFUZZ_INIT:
            o->v[ yMax ] = -1.0;
            o->v[ xSmallest ] = o->min;
            o->v[ xLargest ]= o->max;
            o->v[ sp ] = -1.0;
            break;
        case FIS_DEFUZZ_END:
            d = 0.5*( o->v[ xSmallest ] + o->v[ xLargest ] );
            break;
        default:
            break;
    }

    return d;
}
/*============================================================================*/
real_t fisCore::deFuzzWtAverage( fisOutput * const o, const fisDeFuzzState stage )
{
    real_t d = 0.0;

    if ( FIS_DEFUZZ_END == stage ) {
        d = o->v[ sum_wz ]/o->v[ sum_w ];
    }

    return d;
}
/*============================================================================*/
real_t fisCore::deFuzzWtSum( fisOutput * const o, const fisDeFuzzState stage )
{
    real_t d = 0.0;

    if ( FIS_DEFUZZ_END == stage ) {
        d = o->v[ sum_wz ];
    }

    return d;
}
/*============================================================================*/

/*cstat +CERT-FLP32-C_b*/