#include <include/ffmath.hpp>
#include <iostream>
#include <cmath>

using namespace qlibs;

template<typename T1, typename T2>
static inline void cast_reinterpret( T1 &f, const T2 &u )
{
    static_assert( sizeof(T1) == sizeof(T2), "Types must match sizes" );
    (void)memcpy( &f, &u, sizeof(T2) );
}

static float getAbnormal( const int i );
static float compute_cbrt( float x ,
                           bool r );
static float lgamma_positive( float x );
static float poly_laguerre_recursion( unsigned int n,
                                      float alpha,
                                      float x );
static float poly_laguerre_large_n( unsigned int n,
                                    float alpha,
                                    float x );
static float poly_laguerre_hyperg( unsigned int n,
                                   float alpha,
                                   float x );
static float poly_legendre_p( unsigned int l,
                              float x );
static float ellint_rf( float x,
                        float y,
                        float z );
static float ellint_rd( float x,
                        float y,
                        float z );
static float ellint_rc( float x,
                        float y );
static float ellint_rj( float x,
                        float y,
                        float z,
                        float p );

/*============================================================================*/
static inline float absolute( const float x )
{
    return ( x < 0.0F ) ? -x : x;
}
/*============================================================================*/
static float getAbnormal( const int i )
{
    static const uint32_t u_ab[ 2 ] = { 0x7F800000U, 0x7FBFFFFFU };
    static float f_ab[ 2 ] = { 0.0F, 0.0F };
    static bool init = true;

    if ( init ) {
        cast_reinterpret( f_ab, u_ab );
        init = false;
    }

    return f_ab[ i ];
}
/*============================================================================*/
bool ffmath::isEqual( const float a,
                      const float b,
                      const float tol ) noexcept
{
    return ( absolute( a - b ) <= absolute( tol ) );
}
/*============================================================================*/
float ffmath::getInf( void )
{
    return getAbnormal( 0 );
}
/*============================================================================*/
float ffmath::getNan( void )
{
    return getAbnormal( 1 );
}
/*============================================================================*/
ffmath::classification ffmath::classify( const float f )
{
    uint32_t u = 0U;
    ffmath::classification retVal;

    cast_reinterpret( u, f );
    u &= 0x7FFFFFFFU;

    if ( 0U == u ) {
        retVal = ffmath::classification::FFP_ZERO;
    }
    else if ( u < 0x00800000U ) {
        retVal = ffmath::classification::FFP_SUBNORMAL;
    }
    else if ( u < 0x7F800000U ) {
        retVal = ffmath::classification::FFP_NORMAL;
    }
    else if ( 0x7F800000U == u ) {
        retVal = ffmath::classification::FFP_INFINITE;
    }
    else {
        retVal = ffmath::classification::FFP_NAN;
    }

    return retVal;
}
/*============================================================================*/
float ffmath::sign( float x )
{
    float s;

    if ( x > 0.0F ) {
        s = 1.0F;
    }
    else if ( x < 0.0F ) {
        s = -1.0F;
    }
    else if ( ffmath::classification::FFP_ZERO == ffmath::classify( x ) ) {
        s = 0.0F;
    }
    else {
        s = ffmath::getNan();
    }

    return s;
}
/*============================================================================*/
float ffmath::absf( float x )
{
    return absolute( x );
}
/*============================================================================*/
float ffmath::recip( float x )
{
    uint32_t y = 0U;
    float z = 0.0F;

    cast_reinterpret( y, x );
    y = 0x7EF311C7U - y;
    cast_reinterpret( z, y );

    return z*( 2.0F - ( x*z ) );
}
/*============================================================================*/
float ffmath::sqrt( float x )
{
    float retVal;

    if ( x < 0.0F ) {
        retVal = ffmath::getNan();
    }
    else if ( ffmath::classification::FFP_ZERO == ffmath::classify( x ) ) {
        retVal = 0.0F;
    }
    else {
        uint32_t y = 0U;
        float z = 0.0F;

        cast_reinterpret( y, x );
        y = ( ( y - 0x00800000U ) >> 1U ) + 0x20000000U;
        cast_reinterpret( z, y );
        z = 0.5F*( ( x/z ) + z );
        retVal = 0.5F*( ( x/z ) + z );
    }

    return retVal;
}
/*============================================================================*/
float ffmath::rSqrt( float x )
{
    float retVal;

    if ( x < 0.0F ) {
        retVal = ffmath::getNan();
    }
    else if ( ffmath::classification::FFP_ZERO  == ffmath::classify( x ) ) {
        retVal = ffmath::getInf();
    }
    else {
        uint32_t y = 0U;
        const float z = 0.5F*x;

        cast_reinterpret( y, x );
        y = 0x5F375A86U - ( y >> 1U );
        cast_reinterpret( x, y );
        retVal = x*( 1.5F - ( z*x*x ) );
    }

    return retVal;
}
/*============================================================================*/
static float compute_cbrt( float x , bool r )
{
    float retVal, y = 0.0F, c, d;
    const float k[ 3 ] = { 1.752319676F, 1.2509524245F, 0.5093818292F };
    uint32_t i = 0U;
    bool neg = false;

    if ( x < 0.0F ) {
        x = -x;
        neg = true;
    }
    cast_reinterpret( i, x );
    i = 0x548C2B4BU - ( i/3U );
    cast_reinterpret( y, i );
    c = x*y*y*y;
    y = y*( k[ 0 ] - ( c*( k[ 1 ] - ( k[ 2 ]*c ) ) ) );

    d = x*y*y;
    c = 1.0F - ( d*y );
    retVal = 1.0F + ( 0.333333333333F*c );
    retVal *= ( r ) ? y : d;
    return ( neg )? -retVal : retVal;
}
/*============================================================================*/
float ffmath::cbrt( float x )
{
    return compute_cbrt( x, false );
}
/*============================================================================*/
float ffmath::rCbrt( float x )
{
    return ( ffmath::classification::FFP_ZERO  == ffmath::classify( x ) ) ?
             ffmath::getInf() : compute_cbrt( x, true );
}
/*============================================================================*/
float ffmath::rounding( float x )
{
    x += 12582912.0F;
    x -= 12582912.0F;
    return x;
}
/*============================================================================*/
float ffmath::floor( float x )
{
    return ffmath::rounding( x - 0.5F );
}
/*============================================================================*/
float ffmath::ceil( float x )
{
    return ffmath::rounding( x + 0.5F );
}
/*============================================================================*/
float ffmath::trunc( float x )
{
    /*cstat -CERT-FLP34-C -CERT-FLP36-C*/
    return static_cast<float>( static_cast<int32_t>( x ) );
    /*cstat +CERT-FLP34-C +CERT-FLP36-C*/
}
/*============================================================================*/
float ffmath::frac( float x )
{
    return x - ffmath::trunc( x );
}
/*============================================================================*/
float ffmath::rem( float x,
                   float y )
{
    return ( classification::FFP_ZERO  == classify( x ) ) ? ffmath::getNan()
                                                          : ( x - ( y*ffmath::trunc( x/y ) ) );
}
/*============================================================================*/
float ffmath::mod( float x,
                   float y )
{
    float m;

    if ( classification::FFP_ZERO  == classify( y ) ) {
        m = x;
    }
    else {
        m = x - ( y*ffmath::floor( x/y ) );
        if ( y > 0.0F ) {
            if ( m >= y ) {
                m = 0.0F;
            }
            else if ( m < 0.0F ) {
                const float tmp = y + m;
                m = ( absolute( tmp - y ) <= 1.175494351e-38F ) ? 0.0F : tmp;
            }
            else {
                /*nothing to do here*/
            }
        }
        else {
            if ( m <= y ) {
                m = 0.0F;
            }
            else if ( m > 0.0F ) {
                const float tmp = y + m;
                m = ( absolute( tmp - y ) <= 1.175494351e-38F ) ? 0.0F : tmp;
            }
            else {
                /*nothing to do here*/
            }
        }
    }

    return m;
}
/*============================================================================*/
float ffmath::sin( float x )
{
    float y;
    if ( absolute( x ) <= 0.0066F ) {
        y = x;
    }
    else {
        x *= -ffmath::FFP_1_PI;
        y = x + 25165824.0F;
        x -= y - 25165824.0F;
        x *= absolute( x ) - 1.0F;
        y = x*( ( 3.5841304553896F*absolute( x ) ) + 3.1039673861526F );
    }

    return y;
}
/*============================================================================*/
float ffmath::cos( float x )
{
    float y;
    const float abs_x = absolute( x );

    if ( ffmath::isEqual( abs_x, ffmath::FFP_PI_2 ) ) {
        y = 1.0e-12F;
    }
    else {
        y = ffmath::sin( x + ffmath::FFP_PI_2 );
    }

    return y;
}
/*============================================================================*/
float ffmath::tan( float x )
{
    return ffmath::sin( x )/ffmath::cos( x );
}
/*============================================================================*/
float ffmath::asin( float x )
{
    x = ffmath::sqrt( 1.0F + x ) - ffmath::sqrt( 1.0F - x );
    return x*( ( 0.131754508171F*absolute( x ) ) + 0.924391722181F );
}
/*============================================================================*/
float ffmath::acos( float x )
{
    return ffmath::FFP_PI_2 - ffmath::asin( x );
}
/*============================================================================*/
float ffmath::atan( float x )
{
    float abs_x;

    x /= absolute( x ) + 1.0F;
    abs_x = absolute( x );

    return x*( ( abs_x*( ( -1.45667498914F*abs_x ) + 2.18501248371F ) ) + 0.842458832225F );
}
/*============================================================================*/
float ffmath::atan2( float y,
                     float x )
{
    float t, f;

    t = ffmath::FFP_PI - ( ( y < 0.0F ) ? 6.28318530717958647692F : 0.0F );
    f = ( absolute( x ) <= 1.175494351e-38F ) ? 1.0F : 0.0F;
    y = ffmath::atan( y/( x + f ) ) + ( ( x < 0.0F ) ? t : 0.0F );

    return y + ( f*( ( 0.5F*t ) - y ) );
}
/*============================================================================*/
float ffmath::exp2( float x )
{
    float retVal;

    if ( x <= -126.0F ) {
        retVal = 0.0F;
    }
    else if ( x > 128.0F ) {
        retVal = ffmath::getInf();
    }
    else {
        float y = 0.0F;
        uint32_t exponent;
        /*cstat -MISRAC++2008-5-0-7 -CERT-FLP34-C*/
        exponent = static_cast<uint32_t>( x + 127.0F );
        /*cstat +MISRAC++2008-5-0-7 +CERT-FLP34-C*/
        /*cstat -CERT-FLP36-C*/
        x += 127.0F - static_cast<float>( exponent );
        /*cstat +CERT-FLP36-C*/
        exponent <<= 23U;
        cast_reinterpret( y, exponent );
        x *= ( x*0.339766027F ) + 0.660233972F;
        retVal = y*( x + 1.0F );
    }

    return retVal;
}
/*============================================================================*/
float ffmath::log2( float x )
{
    float retVal;

    if ( x < 0.0F ) {
        retVal = ffmath::getNan();
    }
    else if ( ffmath::classification::FFP_ZERO == ffmath::classify( x ) ) {
        retVal = -ffmath::getInf();
    }
    else {
        uint32_t y = 0U, y2;

        cast_reinterpret( y, x );
        y2 = y;
        y >>= 23U;
        /*cstat -CERT-FLP36-C*/
        retVal = static_cast<float>( y );
        /*cstat +CERT-FLP36-C*/
        y = ( y2 & 0x007FFFFFU ) | 0x3F800000U;
        cast_reinterpret( x, y );
        retVal += -128.0F + ( x*( ( -0.333333333F*x ) + 2.0F ) ) - 0.666666666F;
    }

    return retVal;
}
/*============================================================================*/
float ffmath::exp( float x )
{
    return ffmath::exp2( ffmath::FFP_LOG2E*x );
}
/*============================================================================*/
float ffmath::exp10( float x )
{
    return ffmath::exp2( 3.32192809F*x );
}
/*============================================================================*/
float ffmath::log( float x )
{
    return ffmath::FFP_LN2*ffmath::log2(x);
}
/*============================================================================*/
float ffmath::log10( float x )
{
    return 0.301029996F*ffmath::log2(x);
}
/*============================================================================*/
float ffmath::pow( float b,
                   float e )
{
    return ffmath::exp2( e*ffmath::log2( b ) );
}
/*============================================================================*/
float ffmath::sinh( float x )
{
    x = ffmath::exp( x );
    return ( ( x - 1.0F )/x )*0.5F;
}
/*============================================================================*/
float ffmath::cosh( float x )
{
    x = ffmath::exp( x );
    return ( ( x + 1.0F )/x )*0.5F;
}
/*============================================================================*/
float ffmath::tanh( float x )
{
    x = ffmath::exp( -2.0F*x );
    return ( 1.0F - x )/( 1.0F + x );
}
/*============================================================================*/
float ffmath::asinh( float x )
{
    return ffmath::log( x + ffmath::sqrt( ( x*x ) + 1.0F ) );
}
/*============================================================================*/
float ffmath::acosh( float x )
{
    return ( x < 1.0F ) ? ffmath::getNan()
                        : ffmath::log( x + ffmath::sqrt( ( x*x ) - 1.0F ) );
}
/*============================================================================*/
float ffmath::atanh( float x )
{
    return ffmath::log( ( 1.0F + x )/( 1.0F - x ) )*0.5F;
}
/*============================================================================*/
float ffmath::wrapToPi( float x )
{
    return ffmath::mod( x + ffmath::FFP_PI,  6.28318530717958647692F ) - ffmath::FFP_PI;
}
/*============================================================================*/
float ffmath::wrapTo2Pi( float x )
{
    return ffmath::mod( x, 6.28318530717958647692F );
}
/*============================================================================*/
float ffmath::wrapTo180( float x )
{
    return ffmath::mod( x + 180.0F, 360.0F ) - 180.0F;
}
/*============================================================================*/
float ffmath::wrapTo360( float x )
{
    return ffmath::mod( x, 360.0F );
}
/*============================================================================*/
float ffmath::midpoint( float a,
                        float b )
{
    float y;
    constexpr float lo = 2.0F*1.175494351e-38F;
    constexpr float hi = 0.5F*3.402823466e+38F;
    const float abs_a = absolute( a );
    const float abs_b = absolute( b );

    if ( ( abs_a <= hi ) && ( abs_b <= hi ) ) {
        y = 0.5F*( a + b );
    }
    else if ( abs_a < lo ) {
        y = a + ( 0.5F*b );
    }
    else if ( abs_b < lo) {
        y = ( 0.5F*a ) + b;
    }
    else {
        y = ( 0.5F*a ) + ( 0.5F*b );
    }

    return y;
}
/*============================================================================*/
float ffmath::lerp( float a,
                    float b,
                    float t )
{
    float y;

    if ( ( ( a <= 0.0F ) && ( b >= 0.0F ) ) || ( ( a >= 0.0F ) && ( b <= 0.0F ) ) ) {
        y = ( t*b ) + ( a*( 1.0F - t ) );
    }
    else if ( ffmath::isEqual( t, 1.0F ) ) {
        y = b;
    }
    else {
        const float x = a + t*( b - a );
        y = ( ( t > 1.0F ) == ( b > a ) ) ? ( ( b < x ) ? x : b )
                                          : ( ( b > x ) ? x : b );
    }

    return y;
}
/*============================================================================*/
float ffmath::normalize( const float x,
                         const float xMin,
                         const float xMax ) noexcept
{
    return ( x - xMin )/( xMax - xMin );
}
/*============================================================================*/
float ffmath::map(  const float x,
                    const float xMin,
                    const float xMax,
                    const float yMin,
                    const float yMax ) noexcept
{
    return ( ( yMax - yMin )*ffmath::normalize( x, xMin, xMax ) ) + yMin;
}
/*============================================================================*/
bool ffmath::inRangeCoerce( float &x,
                            const float lowerL,
                            const float upperL ) noexcept
{
    bool retVal = false;

    if ( ffmath::isNan( x ) ) {
        x = lowerL;
    }
    else {
        if ( x < lowerL ) {
            x = lowerL;
        }
        else if ( x > upperL ) {
            x = upperL;
        }
        else {
            retVal = true;
        }
    }

    return retVal;
}
/*============================================================================*/
bool ffmath::inPolygon( const float x,
                        const float y,
                        const float * const px,
                        const float * const py,
                        const size_t p ) noexcept
{
    size_t i;
    bool retVal = false;
    float max_y = py[ 0 ], max_x = px[ 0 ], min_y = py[ 0 ], min_x = px[ 0 ];

    for ( i = 0U ; i < p ; ++i ) {
        max_y = ( py[ i ] > max_y ) ? py[ i ] : max_y;
        max_x = ( px[ i ] > max_x ) ? px[ i ] : max_x;
        min_y = ( py[ i ] < min_y ) ? py[ i ] : min_y;
        min_x = ( px[ i ] < min_x ) ? px[ i ] : min_x;
    }

    if ( ( y >= min_y ) && ( y <= max_y ) && ( x >= min_x ) && ( x <= max_x ) ) {
        size_t j = p - 1U;

        for ( i = 0U ; i < p ; ++i ) {
            if ( ( px[ i ] > x ) != ( px[ j ] > x ) ) {
                const float dx = px[ j ] - px[ i ];
                const float dy = py[ j ] - py[ i ];
                if ( y < ( ( dy*( x - px[ i ] ) )/( dx + py[ i ] ) ) ) {
                    retVal = !retVal;
                }
            }
            j = i;
        }
    }

    return retVal;
}
/*============================================================================*/
bool ffmath::inCircle( const float x,
                       const float y,
                       const float cx,
                       const float cy,
                       const float r ) noexcept
{
    const float d = ( ( x - cx )*( x - cx ) ) + ( ( y - cy )*( y - cy ) );
    return ( d <= ( r*r ) );
}
/*============================================================================*/
float ffmath::erf( float x )
{
    float retVal;

    if ( x >= 6.912F ) {
        retVal = 1.0F;
    }
    else {
        x = ffmath::exp( 3.472034176F*x );
        retVal = ( x/( ( absolute( x ) + 1.0F )*2.0F ) ) - 1.0F;
    }

    return retVal;
}
/*============================================================================*/
float ffmath::erfc( float x )
{
    return 1.0F - ffmath::erf( x );
}
/*============================================================================*/
float ffmath::rexp( float x,
                    int32_t *pw2 )
{
    uint32_t lu = 0U, iu;
    int32_t i = 0;

    cast_reinterpret( lu, x );
    iu  = ( lu >> 23U ) & 0x000000FFU;  /* Find the exponent (power of 2) */
    cast_reinterpret( i, iu );
    i -= 0x7E;
    pw2[ 0 ] = static_cast<int>( i );
    lu &= 0x807FFFFFU; /* strip all exponent bits */
    lu |= 0x3F000000U; /* mantissa between 0.5 and 1 */
    cast_reinterpret( x, lu );

    return x;
}
/*============================================================================*/
float ffmath::ldexp( float x,
                     int32_t pw2 )
{
    uint32_t lu = 0U, eu;
    int32_t e = 0;

    cast_reinterpret( lu, x );
    eu = ( lu >> 23U ) & 0x000000FFU;
    cast_reinterpret( e, eu );
    e += pw2;
    cast_reinterpret( eu, e );
    lu = ( ( eu & 0xFFU ) << 23U ) | ( lu & 0x807FFFFFU );
    cast_reinterpret( x, lu );

    return x;
}
/*============================================================================*/
float ffmath::hypot( float x,
                     float y )
{
    float retVal;
    const auto xClass = ffmath::classify( x );
    const auto yClass = ffmath::classify( y );

    if ( ( xClass < classification::FFP_INFINITE ) && ( yClass < classification::FFP_INFINITE ) ) {
        float a, b, an, bn;
        int32_t e = 0;

        if ( x >= y ) {
            a = x;
            b = y;
        }
        else {
            a = y;
            b = x;
        }
        /* Write a = an * 2^e, b = bn * 2^e with 0 <= bn <= an < 1.*/
        an = ffmath::rexp( a, &e );
        bn = ffmath::ldexp( b, -e );
        retVal = ffmath::sqrt( ( an*an ) + ( bn*bn ) );
        retVal = ffmath::ldexp( retVal, e );
    }
    else {
        retVal = ( ( classification::FFP_INFINITE == xClass ) ||
                   ( classification::FFP_INFINITE == yClass )  ) ? ffmath::getInf()
                                                                 : ffmath::getNan();
    }

    return retVal;
}
/*============================================================================*/
float ffmath::nextAfter( float x,
                         float y )
{
    float retVal = 0.0F;
    uint32_t ax, ay, uxi = 0U, uyi = 0U;

    cast_reinterpret( uxi, x );
    cast_reinterpret( uyi, y );
    if ( ffmath::isNan( x ) || ffmath::isNan( y ) ) {
        retVal = ffmath::getNan();
    }
    else if ( uxi == uyi ) {
        retVal = y;
    }
    else {
        ax = uxi & 0x7FFFFFFFU;
        ay = uyi & 0x7FFFFFFFU;
        if ( 0U == ax ) {
            uxi = ( 0U == ay ) ? uyi : ( ( uyi & 0x80000000U ) | 1U );
        }
        else if ( ( ax > ay ) || ( 0U != ( ( uxi^uyi ) & 0x80000000U ) ) ) {
            --uxi;
        }
        else {
            ++uxi;
        }
        cast_reinterpret( retVal, uxi );
    }

    return retVal;
}
/*============================================================================*/
float ffmath::tgamma( float x )
{
    float result;

    const auto fClass = ffmath::classify( x );
    if ( classification::FFP_NAN == fClass ) {
        result = ffmath::getNan();
    }
    else if ( classification::FFP_ZERO == fClass ) {
        result = ffmath::getInf(); /* a huge value */
    }
    else if ( classification::FFP_INFINITE == fClass ) {
        result = ( x > 0.0F ) ? ffmath::getInf() : ffmath::getNan();
    }
    else {
        bool parity = false;
        float fact = 1.0F;
        float y = x;
        float y1;

        if ( y <= 0.0F ) {
            float isItAnInt;

            y = -x;
            y1 = ffmath::trunc( y );
            isItAnInt = y - y1;
            if ( !ffmath::isEqual( 0.0F, isItAnInt ) ) {
                const float tmp = 2.0F*ffmath::trunc( y1*0.5F );

                if ( !ffmath::isEqual( y1, tmp ) ) {
                    parity = true;
                }
                fact = -FFP_PI/ffmath::sin( FFP_PI*isItAnInt );
                y += 1.0F;
            }
            else {
                result = ffmath::getNan();
                goto EXIT_TGAMMA;
            }
        }
        if ( y < 1.19209290E-07F ) { /* y < eps */
            if ( y >= 1.175494351e-38F ) { /* y >= MinimumX */
                result = 1.0F/y;
            }
            else {
                result = ffmath::getInf();
            }
        }
        else if ( y < 12.0F ) {
            float num = 0.0F, den = 1.0F, z;
            int n = 0;

            y1 = y;
            if ( y < 1.0F ) {
                z = y;
                y += 1.0F;
            }
            else {
                n = static_cast<int>( y ) - 1;
                /*cstat -CERT-FLP36-C */
                y -= static_cast<float>( n );
                /*cstat +CERT-FLP36-C */
                z = y - 1.0F;
            }

            num = z*( num + -1.71618513886549492533811e+0F );
            den = ( den*z ) -3.08402300119738975254353e+1F;
            num = z*( num + 2.47656508055759199108314e+1F );
            den = ( den*z ) + 3.15350626979604161529144e+2F;
            num = z*( num - 3.79804256470945635097577e+2F );
            den = ( den*z ) - 1.01515636749021914166146e+3F;
            num = z*( num + 6.29331155312818442661052e+2F );
            den = ( den*z ) - 3.10777167157231109440444e+3F;
            num = z*( num + 8.66966202790413211295064e+2F );
            den = ( den*z ) + 2.25381184209801510330112e+4F;
            num = z*( num - 3.14512729688483675254357e+4F );
            den = ( den*z ) + 4.75584627752788110767815e+3F;
            num = z*( num - 3.61444134186911729807069e+4F );
            den = ( den*z ) - 1.34659959864969306392456e+5F;
            num = z*( num + 6.64561438202405440627855e+4F );
            den = ( den*z ) - 1.15132259675553483497211e+5F;

            result = ( num/den ) + 1.0F;
            if ( y1 < y ) {
                  result /= y1;
            }
            else if ( y1 > y ) {
                for ( int i = 0; i < n ; ++i ) {
                    result *= y;
                    y += 1.0F;
                }
            }
        }
        else {
            if ( x <= 171.624F ) { /* x <= xBig */
                const float yy  = y*y;
                float sum = 5.7083835261e-03F;

                sum = ( sum/yy ) - 1.910444077728e-03F;
                sum = ( sum/yy ) + 8.4171387781295e-04F;
                sum = ( sum/yy ) - 5.952379913043012e-04F;
                sum = ( sum/yy ) + 7.93650793500350248e-04F;
                sum = ( sum/yy ) - 2.777777777777681622553e-03F;
                sum = ( sum/yy ) + 8.333333333333333331554247e-02F;

                sum = ( sum /y ) - y + FFP_LN_SQRT_2PI;
                sum += ( y - 0.5F )*ffmath::log( y );
                result = ffmath::exp( sum );
            }
            else {
                result = ffmath::getInf();
            }
        }
        if ( parity ) {
            result = -result;
        }
        if ( !ffmath::isEqual( fact, 1.0F ) ) {
            result = fact/result;
        }
    }

    EXIT_TGAMMA:
    return result;
}
/*============================================================================*/
static float lgamma_positive( float x )
{
    constexpr float d1 = -5.772156649015328605195174e-1F;
    constexpr float d2 = 4.227843350984671393993777e-1F;
    constexpr float d4 = 1.791759469228055000094023e+0F;
    constexpr float p1[ 8 ] = { 4.945235359296727046734888e+0F,
                                2.018112620856775083915565e+2F,
                                2.290838373831346393026739e+3F,
                                1.131967205903380828685045e+4F,
                                2.855724635671635335736389e+4F,
                                3.848496228443793359990269e+4F,
                                2.637748787624195437963534e+4F,
                                7.225813979700288197698961e+3F };
    constexpr float q1[ 8 ] = { 6.748212550303777196073036e+1F,
                                1.113332393857199323513008e+3F,
                                7.738757056935398733233834e+3F,
                                2.763987074403340708898585e+4F,
                                5.499310206226157329794414e+4F,
                                6.161122180066002127833352e+4F,
                                3.635127591501940507276287e+4F,
                                8.785536302431013170870835e+3F };
    constexpr float p2[ 8 ] = { 4.974607845568932035012064e+0F,
                                5.424138599891070494101986e+2F,
                                1.550693864978364947665077e+4F,
                                1.847932904445632425417223e+5F,
                                1.088204769468828767498470e+6F,
                                3.338152967987029735917223e+6F,
                                5.106661678927352456275255e+6F,
                                3.074109054850539556250927e+6F };
    constexpr float q2[ 8 ] = { 1.830328399370592604055942e+2F,
                                7.765049321445005871323047e+3F,
                                1.331903827966074194402448e+5F,
                                1.136705821321969608938755e+6F,
                                5.267964117437946917577538e+6F,
                                1.346701454311101692290052e+7F,
                                1.782736530353274213975932e+7F,
                                9.533095591844353613395747e+6F };
    constexpr float p4[ 8 ] = { 1.474502166059939948905062e+04F,
                                2.426813369486704502836312e+06F,
                                1.214755574045093227939592e+08F,
                                2.663432449630976949898078e+09F,
                                2.940378956634553899906876e+10F,
                                1.702665737765398868392998e+11F,
                                4.926125793377430887588120e+11F,
                                5.606251856223951465078242e+11F };
    constexpr float q4[ 8 ] = { 2.690530175870899333379843e+03F,
                                6.393885654300092398984238e+05F,
                                4.135599930241388052042842e+07F,
                                1.120872109616147941376570e+09F,
                                1.488613728678813811542398e+10F,
                                1.016803586272438228077304e+11F,
                                3.417476345507377132798597e+11F,
                                4.463158187419713286462081e+11F };
    constexpr float pnt68 = 0.6796875F;
    float result;

    if ( x > 171.624F ) {
        result = ffmath::getInf();
    }
    else {
        float y, corrector, num, den;

        y = x;
        if ( y <= 1.19209290E-07F ) { /* y < eps */
            result = -ffmath::log( y );
        }
        else if ( y <= 1.5F ) {
            float xMinus;

            if ( y < pnt68 ) {
                corrector = -ffmath::log( y );
                xMinus = y;
            }
            else {
                corrector = 0.0F;
                xMinus = ( y - 0.5F ) - 0.5F;
            }
            if ( ( y <= 0.5F ) || ( y >= pnt68 ) ) {
                den = 1.0F;
                num = 0.0F;
                for ( size_t i = 0U ; i < 8U ; ++i ) {
                    num = ( num*xMinus ) + p1[ i ];
                    den = ( den*xMinus ) + q1[ i ];
                }
                result = corrector + ( xMinus*( d1 + xMinus*( num/den ) ) );
            }
            else {
                xMinus = ( y - 0.5F ) - 0.5F;
                den = 1.0F;
                num = 0.0F;
                for ( size_t i = 0U; i < 8U ; ++i ) {
                    num = num*xMinus + p2[ i ];
                    den = den*xMinus + q2[ i ];
                }
                result = corrector + ( xMinus*( d2 + xMinus*( num/den ) ) );
            }
        }
        else if ( y <= 4.0F ) {
            const float xMinus = y - 2.0F;
            den = 1.0F;
            num = 0.0F;
            for ( size_t i = 0U; i < 8U ; ++i ) {
                num = num*xMinus + p2[ i ];
                den = den*xMinus + q2[ i ];
            }
            result = xMinus*( d2 + xMinus*( num/den ) );
        }
        else if ( y <= 12.0F ) {
            const float xMinus = y - 4.0F;
            den = -1.0F;
            num = 0.0F;
            for ( size_t i = 0U; i < 8U ; ++i ) {
                num = num*xMinus + p4[ i ];
                den = den*xMinus + q4[ i ];
            }
            result = d4 + xMinus*( num/den );
        }
        else {
            result = 0.0F;
            if ( y <= 4294967296.87842273712158203125F ) { /* y < xBig^(1/4)*/
                const float yy = y*y;
                result = 5.7083835261e-03F;
                result = ( result/yy ) - 1.910444077728e-03F;
                result = ( result/yy ) + 8.4171387781295e-04F;
                result = ( result/yy ) - 5.952379913043012e-04F;
                result = ( result/yy ) + 7.93650793500350248e-04F;
                result = ( result/yy ) - 2.777777777777681622553e-03F;
                result = ( result/yy ) + 8.333333333333333331554247e-02F;
            }
            result /= y;
            corrector = ffmath::log( y );
            result += ffmath::FFP_LN_SQRT_2PI - ( 0.5F*corrector );
            result += y*( corrector - 1.0F );
        }
    }

    return result;
}
/*============================================================================*/
float ffmath::lgamma( float x )
{
    float result;

    const auto fClass = ffmath::classify( x );
    if ( classification::FFP_NAN == fClass ) {
        result = ffmath::getNan();
    }
    else if ( ( classification::FFP_ZERO == fClass ) || ( classification::FFP_INFINITE == fClass ) ) {
        result = ffmath::getInf();
    }
    else {
        if ( x < 0.0F ) {
            if ( x <= -4503599627370496.0F ) { /* x < 2^52 */
                result = ffmath::getInf();
            }
            else {
                float y, y1, isItAnInt;

                y = -x;
                y1 = ffmath::trunc( y );
                isItAnInt = y - y1;
                if ( ffmath::isEqual( 0.0F, isItAnInt ) ) {
                    result = ffmath::getInf();
                }
                else {
                    float a;

                    a = ffmath::sin( FFP_PI*isItAnInt );
                    result = ffmath::log( FFP_PI/absolute( a*x ) ) - lgamma_positive( -x );
                }
            }
        }
        else {
            result = lgamma_positive( x );
        }
    }

    return result;
}
/*============================================================================*/
float ffmath::factorial( float x )
{
    constexpr float ft[ 35 ] = { 1.0F, 1.0F, 2.0F, 6.0F, 24.0F, 120.0F, 720.0F,
                                 5040.0F, 40320.0F, 362880.0F, 3628800.0F,
                                 39916800.0F, 479001600.0F, 6227020800.0F,
                                 87178291200.0F, 1307674368000.0F,
                                 20922789888000.0F, 355687428096000.0F,
                                 6402373705728001.0F, 121645100408832000.0F,
                                 2432902008176640000.0F, 51090942171709440000.0F,
                                 1124000727777607680000.0F,
                                 25852016738884978212864.0F,
                                 620448401733239544217600.0F,
                                 15511210043330988202786816.0F,
                                 403291461126605719042260992.0F,
                                 10888869450418351940239884288.0F,
                                 304888344611713836734530715648.0F,
                                 8841761993739700772720181510144.0F,
                                 265252859812191104246398737973248.0F,
                                 8222838654177921277277005322125312.0F,
                                 263130836933693591553328612565319680.0F,
                                 8683317618811885938715673895318323200.0F,
                                 295232799039604119555149671006000381952.0F,
                                 };
    float y;

    if ( x > 34.0F ) {
        y = ffmath::getInf();
    }
    else if ( x >= 0.0F ) {
        const auto i = static_cast<size_t>( x );
        y = ft[ i ];
    }
    else {
        y = ffmath::getNan();
    }

    return y;
}
/*============================================================================*/
static float poly_laguerre_recursion( unsigned int n,
                                      float alpha,
                                      float x )
{
    const float l0 = 1.0F;
    float y;

    if ( 0U == n ) {
        y = l0;
    }
    else {
        const float l1 = -x + 1.0F + alpha;
        if ( 1U == n ) {
            y = l1;
        }
        else {
            float ln2 = l0;
            float ln1 = l1;
            float ln = 0.0F;
            for ( size_t i = 2U ; i <= n ; ++i ) {
                /*cstat -CERT-FLP36-C*/
                const auto nn = static_cast<float>( i );
                /*cstat +CERT-FLP36-C*/
                ln = ( ( ( 2.0F*nn ) - 1.0F ) + alpha - x )*( ln1/nn )
                     - ( ( nn - 1.0F ) + alpha )*( ln2/nn );
                ln2 = ln1;
                ln1 = ln;
            }
            y = ln;
        }
    }

    return y;
}
/*============================================================================*/
static float poly_laguerre_large_n( unsigned int n,
                                    float alpha,
                                    float x )
{
    constexpr float PI_2_SQ = 2.467401100272339498076235031476244330406188964F;
    /*cstat -CERT-FLP36-C*/
    const float m = static_cast<float>( n );
    /*cstat +CERT-FLP36-C*/
    const float a = -m;
    const float b = alpha + 1.0F;
    const float eta = ( 2.0F*b ) - ( 4.0F*a );
    const float cos2th = x/eta;
    const float sin2th = 1.0F - cos2th;
    const float th = ffmath::acos( ffmath::sqrt( cos2th ) );
    const float pre_h = PI_2_SQ*eta*eta*cos2th*sin2th;
    const float lg_b = ffmath::lgamma( b + m );
    const float ln_fact = ffmath::lgamma( m + 1.0F );

    const float preTerm1 = 0.5F*( 1.0F - b )*ffmath::log( 0.25F*x*eta );
    const float preTerm2 = 0.25F*ffmath::log( pre_h );
    const float lnPre = lg_b - ln_fact + ( 0.5F*x ) + preTerm1 - preTerm2;
    const float serTerm1 = ffmath::sin( ffmath::FFP_PI*a );
    const float serTerm2 = ffmath::sin( ( 0.25F*eta*( ( 2.0F*th ) - ffmath::sin( 2.0F*th ) ) + ffmath::FFP_PI_4 ) );

    return ffmath::exp( lnPre )*( serTerm1 + serTerm2 );
}
/*============================================================================*/
static float poly_laguerre_hyperg( unsigned int n,
                                   float alpha,
                                   float x )
{
    const float b = alpha + 1.0F;
    const float mx = -x;
    const float tc_sgn = ( x < 0.0F ) ? 1.0F : ( ( 1 == ( n % 2 ) ) ? -1.0F : 1.0F );
    const float ax = absolute( x );
    float tc = 1.0F;

    for ( size_t i = 1U ; i <= n ; ++i ) {
        /*cstat -CERT-FLP36-C*/
        const float k = static_cast<float>( i );
        /*cstat +CERT-FLP36-C*/
        tc *= ax/k;
    }

    float term = tc*tc_sgn;
    float sum = term;
    const int N = static_cast<int>( n );
    for ( int i = ( N - 1 ) ; i >= 0; --i ) {
        /*cstat -CERT-FLP36-C -MISRAC++2008-5-0-7*/
        const float k = static_cast<float>( i );
        term *= ( b + k )/static_cast<float>( N - i )*( k + 1.0F )/mx;
        /*cstat +CERT-FLP36-C +MISRAC++2008-5-0-7*/
        sum += term;
    }

    return sum;
}
/*============================================================================*/
float ffmath::assoc_laguerre( unsigned int n,
                              unsigned int m,
                              float x )
{
    // include/tr1/poly_laguerre.tcc
    float y;
    /*cstat -CERT-FLP36-C*/
    const float alpha = static_cast<float>( m );
    const float N = static_cast<float>( n );
    /*cstat +CERT-FLP36-C*/
    if ( ( x < 0.0f ) || ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else if ( 0U == n ) {
        y = 1.0F;
    }
    else if ( 1U == n ) {
        y = 1.0F + alpha - x;
    }
    else if ( ffmath::isEqual( 0.0F, x ) ) {
        float prod = alpha + 1.0F;
        for ( size_t i = 2U ; i < n ; ++i ) {
            /*cstat -CERT-FLP36-C*/
            const float k = static_cast<float>( i );
            /*cstat +CERT-FLP36-C*/
            prod *= ( alpha + k )/k;
        }
        y = prod;
    }
    else if ( ( n > 10000000 ) && ( alpha > -1.0F ) && ( x < ( ( 2.0F*( alpha + 1.0F ) ) + ( 4.0F*N ) ) ) ) {
        y = poly_laguerre_large_n( n, alpha, x );
    }
    else if ( ( alpha >= 0.0F ) || ( ( x > 0.0F) && ( alpha < -( N + 1.0F ) ) ) ) {
        y = poly_laguerre_recursion( n, alpha, x );
    }
    else {
        y = poly_laguerre_hyperg( n, alpha, x );
    }

    return y;
}
/*============================================================================*/
static float poly_legendre_p( unsigned int l,
                              float x )
{
    float y;

    if ( ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else if ( ffmath::isEqual( -1.0F, x ) ) {
        y = ( 1 == ( l % 2 ) ) ? -1.0F : 1.0F;
    }
    else {
        float p_lm2 = 1.0F;

        if ( 0U == l ) {
            y = p_lm2;
        }
        else {
            float p_lm1 = x;

            if ( 1U == l ) {
                y = p_lm1;
            }
            else {
                float p_l = 0.0F;

                for ( size_t i = 2U ; i <= l ; ++i ) {
                    /*cstat -CERT-FLP36-C*/
                    const float ll = static_cast<float>( i );
                    /*cstat +CERT-FLP36-C*/
                    p_l = ( 2.0F*x* p_lm1 ) - p_lm2 - ( ( x*p_lm1 ) - p_lm2)/ll;
                    p_lm2 = p_lm1;
                    p_lm1 = p_l;
                }
                y = p_l;
            }
        }
    }

    return y;
}
/*============================================================================*/
float ffmath::assoc_legendre( unsigned int n,
                              unsigned int m,
                              float x )
{
    float y;
    const float phase = 1.0F;

    if ( m > n ) {
        y = 0.0F;
    }
    else if ( ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else if ( 0U == m ) {
        y = poly_legendre_p( n, x );
    }
    else {
        float p_mm = 1.0F;
        const float root = ffmath::sqrt( 1.0F - x )*ffmath::sqrt( 1.0F + x );
        float fact = 1.0F;

        for ( size_t i = 1U ; i <= m ; ++i ) {
            p_mm *= phase*fact*root;
            fact += 2.0F;
        }

        if ( n == m ) {
            y = p_mm;
        }
        else {
            /*cstat -CERT-FLP36-C*/
            const float p_mp1m = ( 2.0F*static_cast<float>( m ) + 1.0F )*x*p_mm;
            /*cstat +CERT-FLP36-C*/
            if ( n == ( m + 1U ) ) {
                y = p_mp1m;
            }
            else {
                float p_lm2m = p_mm;
                float p_lm1m = p_mp1m;
                float p_lm = 0.0F;
                /*cstat -CERT-FLP36-C*/
                const float M = static_cast<float>( m );
                for ( size_t i = ( m + 2U ) ; i <= n ; ++i ) {
                    const float j = static_cast<float>( i );
                    /*cstat +CERT-FLP36-C*/
                    p_lm = ( ( ( 2.0F*j ) - 1.0F )*x*p_lm1m ) - ( ( j + M - 1.0F )*p_lm2m/( j - M ) );
                    p_lm2m = p_lm1m;
                    p_lm1m = p_lm;
                }
                y = p_lm;
            }
        }
    }

    return y;
}
/*============================================================================*/
float ffmath::beta( float x,
                    float y )
{
    float result;

    if ( ffmath::isNan( x ) || ffmath::isNan( y ) ) {
        result = ffmath::getNan();
    }
    else {
        const float bet = ffmath::lgamma( x ) + ffmath::lgamma( y ) - ffmath::lgamma( x + y );
        result = ffmath::exp( bet );
    }

    return result;
}
/*============================================================================*/
static float ellint_rf( float x,
                        float y,
                        float z )
{
    constexpr float min = 1.175494351e-38F;
    const float loLim = 5.0F*min;
    float result;


    if ( ( x < 0.0F ) || ( y < 0.0F ) || ( z < 0.0F ) ) {
        result = ffmath::getNan();
    }
    else if ( ( ( x + y ) < loLim ) || ( ( x + z ) < loLim ) || ( ( y + z) < loLim ) ) {
        result = ffmath::getNan();
    }
    else {
        constexpr float c0 = 0.25F;
        constexpr float c1 = 0.04166666666666666666666666666667F;
        constexpr float c2 = 0.1F;
        constexpr float c3 = 0.06818181818181818181818181818182F;
        constexpr float c4 = 0.07142857142857142857142857142857F;
        constexpr float errTol = 0.0024607833005759250678823324F;

        float xn = x;
        float yn = y;
        float zn = z;
        float mu, xnDev, ynDev, znDev;
        const size_t maxIter = 100U;
        float epsilon;
        float e2, e3, s;

        for ( size_t iter = 0U ; iter < maxIter ; ++iter ) {
            float abs_xnDev, abs_ynDev, abs_znDev, lambda;
            float xRoot, yRoot, zRoot;

            mu = ( xn + yn + zn )*0.33333333333333333333333333333333F;
            xnDev = 2.0F - ( mu + xn )/mu;
            ynDev = 2.0F - ( mu + yn )/mu;
            znDev = 2.0F - ( mu + zn )/mu;
            abs_xnDev = absolute( xnDev );
            abs_ynDev = absolute( ynDev );
            abs_znDev = absolute( znDev );
            epsilon = ( abs_xnDev > abs_ynDev ) ? abs_xnDev : abs_ynDev;
            epsilon = ( abs_znDev > epsilon ) ? abs_znDev : epsilon;
            if ( epsilon < errTol ) {
                break;
            }
            xRoot = ffmath::sqrt( xn );
            yRoot = ffmath::sqrt( yn );
            zRoot = ffmath::sqrt( zn );
            lambda = xRoot*( yRoot + zRoot ) + ( yRoot*zRoot );
            xn = c0*( xn + lambda );
            yn = c0*( yn + lambda );
            zn = c0*( zn + lambda );
        }
        e2 = xnDev*ynDev;
        e3 = ( e2*znDev );
        e2 = e2 - ( znDev*znDev );
        s  = 1.0F + ( ( c1*e2 ) - c2 - ( c3*e3 ) )*e2  + ( c4*e3 );
        result = s/ffmath::sqrt( mu );
    }

    return result;
}
/*============================================================================*/
static float ellint_rd( float x,
                        float y,
                        float z )
{
    constexpr float errTol = 0.0017400365588678507952624663346341549F;
    constexpr float loLim = 4.103335708781587555782386855921935e-26F;
    float result;

    if ( ( x < 0.0F ) || ( y < 0.0F ) ) {
        result = ffmath::getNan();
    }
    else if ( ( ( x + y ) < loLim ) || ( z < loLim ) ) {
        result = ffmath::getNan();
    }
    else {
        constexpr float c0 = 0.25F;
        constexpr float c1 = 0.214285714285714273819039021873322781175F;
        constexpr float c2 = 0.16666666666666665741480812812369549646F;
        constexpr float c3 = 0.409090909090909116141432377844466827809F;
        constexpr float c4 = 0.1153846153846153910205174497605185024875F;
        float xn = x;
        float yn = y;
        float zn = z;
        float sigma = 0.0F;
        float power4 = 1.0F;
        float epsilon;
        float mu, xnDev, ynDev, znDev;
        float ea, eb, ec, ed, ef, s1, s2;
        const size_t maxIter = 100U;

        for ( size_t iter = 0U ; iter < maxIter ; ++iter ) {
            float abs_xnDev, abs_ynDev, abs_znDev, lambda;
            float xRoot, yRoot, zRoot;

            mu = ( xn + yn + ( 3.0F*zn ) )*0.2F;
            xnDev = ( mu - xn )/mu;
            ynDev = ( mu - yn )/mu;
            znDev = ( mu - zn )/mu;
            abs_xnDev = absolute( xnDev );
            abs_ynDev = absolute( ynDev );
            abs_znDev = absolute( znDev );
            epsilon = ( abs_xnDev > abs_ynDev ) ? abs_xnDev : abs_ynDev;
            epsilon = ( abs_znDev > epsilon ) ? abs_znDev : epsilon;
            if ( epsilon < errTol ) {
                break;
            }
            xRoot = ffmath::sqrt( xn );
            yRoot = ffmath::sqrt( yn );
            zRoot = ffmath::sqrt( zn );
            lambda = xRoot*( yRoot + zRoot ) + ( yRoot*zRoot );
            sigma += power4/( zRoot*( zn + lambda ) );
            power4 *= c0;
            xn = c0*( xn + lambda );
            yn = c0*( yn + lambda );
            zn = c0*( zn + lambda );
        }
        ea = xnDev*ynDev;
        eb = znDev*znDev;
        ec = ea - eb;
        ed = ea - ( 6.0F*eb );
        ef = ed + ec + ec;
        s1 = ed*( -c1 + ( c3*ed/3.0F ) - ( 1.5F*c4*znDev*ef ) );
        s2 = znDev*( ( c2*ef ) + znDev*( -( c3*ec ) - ( znDev*c4 ) - ea ) );
        result = ( 3.0F*sigma ) + power4*ffmath::rSqrt( mu )*( 1.0F + s1 + s2)/mu;
    }

    return result;
}
/*============================================================================*/
static float ellint_rc( float x,
                        float y )
{
    float result;
    constexpr float loLim = 5.8774717550000002558112628881984982848919e-38F;
    constexpr float errTol = 0.049606282877419791144113503378321F;

    if ( ( x < 0.0F ) || ( y < 0.0F ) || ( y < loLim ) ) {
        result = ffmath::getNan();
    }
    else {
        constexpr float c0 = 0.25F;
        constexpr float c1 = 0.142857142857142849212692681248881854116F;
        constexpr float c2 = 0.409090909090909116141432377844466827809F;
        constexpr float c3 = 0.3F;
        constexpr float c4 = 0.375F;
        float xn = x;
        float yn = y;
        const size_t maxIter = 100;
        float mu, s, sn;

        for ( size_t iter = 0U ; iter < maxIter ; ++iter ) {
            float lambda;

            mu = ( xn + 2.0F*yn )*0.333333333333333333333333333333333F;
            sn = ( yn + mu )/mu - 2.0F;
            if ( absolute( sn) < errTol ){
                break;
            }
            lambda = ( 2.0F*ffmath::sqrt( xn )*ffmath::sqrt( yn ) ) + yn;
            xn = c0*( xn + lambda );
            yn = c0*( yn + lambda );
        }
        s = sn*sn*( c3 + sn*( c1 + sn*( c4 + ( sn*c2 ) ) ) );
        result = ( 1.0F + s )*ffmath::rSqrt( mu );
    }

    return result;
}
/*============================================================================*/
static float ellint_rj( float x,
                        float y,
                        float z,
                        float p )
{
    float result;
    constexpr float loLim = 4.103335708781587555782386855921935e-26F;
    constexpr float errTol = 0.049606282877419791144113503378321F;

    if ( ( x < 0.0F ) || ( y < 0.0F ) || ( z < 0.0F ) ) {
        result = ffmath::getNan();
    }
    else if ( ( ( x + y ) < loLim ) || ( ( x + z ) < loLim ) || ( ( y + z) < loLim )|| ( p < loLim ) ) {
        result = ffmath::getNan();
    }
    else {
        constexpr float c0 = 0.25F;
        constexpr float c1 = 0.214285714285714273819039021873322781175F;
        constexpr float c2 = 0.333333333333333333333333333333333333333F;
        constexpr float c3 = 0.136363636363636353543427048862213268876F;
        constexpr float c4 = 0.1153846153846153910205174497605185024875F;
        float xn = x;
        float yn = y;
        float zn = z;
        float pn = p;
        float sigma = 0.0F;
        float power4 = 1.0F;
        float epsilon;
        float mu, xnDev, ynDev, znDev, pnDev;
        float ea, eb, ec, e2, e3, s1, s2, s3;
        const size_t maxIter = 100;


        for ( size_t iter = 0U ; iter < maxIter ; ++iter ) {
            float abs_xnDev, abs_ynDev, abs_znDev, abs_pnDev, lambda, alpha1, alpha2, beta;
            float xRoot, yRoot, zRoot;
            mu = 0.2F*( xn + yn + zn + ( 2.0F*pn ) );
            xnDev = ( mu - xn )/mu;
            ynDev = ( mu - yn )/mu;
            znDev = ( mu - zn )/mu;
            pnDev = ( mu - pn )/mu;
            abs_xnDev = absolute( xnDev );
            abs_ynDev = absolute( ynDev );
            abs_znDev = absolute( znDev );
            abs_pnDev = absolute( pnDev );
            epsilon = ( abs_xnDev > abs_ynDev ) ? abs_xnDev : abs_ynDev;
            epsilon = ( abs_znDev > epsilon ) ? abs_znDev : epsilon;
            epsilon = ( abs_pnDev > epsilon ) ? abs_pnDev : epsilon;
            if ( epsilon < errTol ) {
                break;
            }
            xRoot = ffmath::sqrt( xn );
            yRoot = ffmath::sqrt( yn );
            zRoot = ffmath::sqrt( zn );
            lambda = xRoot*( yRoot + zRoot ) + ( yRoot*zRoot );
            alpha1 = pn*( xRoot + yRoot + zRoot ) + xRoot*yRoot*zRoot;
            alpha2 = alpha1*alpha1;
            beta = pn*( pn + lambda )*( pn + lambda );
            sigma += power4*ellint_rc( alpha2, beta );
            power4 *= c0;
            xn = c0*( xn + lambda );
            yn = c0*( yn + lambda );
            zn = c0*( zn + lambda );
            pn = c0*( pn + lambda );
        }
        ea = xnDev*( ynDev + znDev ) + ynDev*znDev;
        eb = xnDev*ynDev*znDev;
        ec = pnDev*pnDev;
        e2 = ea - 3.0F*ec;
        e3 = eb + 2.0F*pnDev*( ea - ec );
        s1 = 1.0F + e2*( -c1 + 0.75F*c3*e2 - 1.5F*c4*e3 );
        s2 = eb*( 0.5F*c2 + pnDev*( -c3 - c3 + pnDev*c4 ) );
        s3 = pnDev*ea*( c2 - ( pnDev*c3 ) ) - ( c2*pnDev*ec );
        result = 3.0F*sigma + power4*( s1 + s2 + s3)/( mu * ffmath::sqrt( mu ) );
    }

    return result;
}
/*============================================================================*/
float ffmath::comp_ellint_1( float k )
{
    float y;

    if ( ffmath::isNan( k ) || ( absolute( k ) >= 1.0F ) ) {
        y = ffmath::getNan();
    }
    else {
        y = ellint_rf( 0.0F, 1.0F - ( k*k ), 1.0F );
    }

    return y;
}
/*============================================================================*/
float ffmath::comp_ellint_2( float k )
{
    float y;
    const float abs_k = absolute( k );

    if ( ffmath::isNan( k ) || ( abs_k > 1.0F ) ) {
        y = ffmath::getNan();
    }
    else if ( ffmath::isEqual( 1.0F, abs_k ) ) {
        y = 1.0F;
    }
    else {
        const float kk = k*k;
        const float one_m_kk = 1.0F - kk;

        y = ellint_rf( 0.0F, one_m_kk, 1.0F ) - 0.333333333333333F*kk*ellint_rd( 0.0F, one_m_kk, 1.0F );
    }

    return y;
}
/*============================================================================*/
float ffmath::comp_ellint_3( float k,
                             float nu )
{
    float y;
    const float abs_k = absolute( k );

    if ( ffmath::isNan( k ) || ffmath::isNan( nu ) ) {
        y = ffmath::getNan();
    }
    else if ( ffmath::isEqual( 1.0F, nu ) ) {
        y = ffmath::getInf();
    }
    else if ( abs_k > 1.0F ) {
        y = ffmath::getNan();
    }
    else {
        const float kk = k*k;
        const float one_m_kk = 1.0F - kk;

        y = ellint_rf( 0.0F, one_m_kk, 1.0F ) +
            0.333333333333333F*nu*ellint_rj( 0.0F, one_m_kk, 1.0F, 1.0F - nu );
    }

    return y;
}
/*============================================================================*/
float ffmath::ellint_1( float k,
                        float phi )
{
    float y;

    if ( ffmath::isNan( k ) || ffmath::isNan( phi ) || ( absolute(k) > 1.0F ) ) {
        y = ffmath::getNan();
    }
    else {
        const float n = ffmath::floor( phi/ffmath::FFP_PI + 0.5F );
        const float phi_red = phi - n*ffmath::FFP_PI;
        const float s = ffmath::sin( phi_red );
        const float c = ffmath::cos( phi_red );
        const float f = s*ellint_rf( c*c, 1.0F - k*k*s*s, 1.0F );
        if ( ffmath::classification::FFP_ZERO == ffmath::classify( n ) ) {
            y = f;
        }
        else {
            y = f + 2.0F*n*ffmath::comp_ellint_1( k );
        }
    }

    return y;
}
/*============================================================================*/
float ffmath::ellint_2( float k,
                        float phi )
{
    float y;

    if ( ffmath::isNan( k ) || ffmath::isNan( phi ) || ( absolute( k ) > 1.0F ) ) {
        y = ffmath::getNan();
    }
    else {
        const float n = ffmath::floor( phi/ffmath::FFP_PI + 0.5F );
        const float phi_red = phi - n*ffmath::FFP_PI;
        const float kk = k*k;
        const float s = ffmath::sin( phi_red );
        const float ss = s*s;
        const float sss = ss*s;
        const float c = ffmath::cos( phi_red );
        const float cc = c*c;
        const float tmp = 1.0F - kk*ss;
        const float e = s*ellint_rf( cc, tmp, 1.0F )
                        - 0.333333333333F*kk*sss*ellint_rd( cc, tmp, 1.0F );

        if ( ffmath::classification::FFP_ZERO == ffmath::classify( n ) ) {
            y = e;
        }
        else {
            y = e + 2.0F*n*ffmath::comp_ellint_2( k );
        }
    }

    return y;
}
/*============================================================================*/
float ffmath::ellint_3( float k,
                        float nu,
                        float phi )
{
    float y;

    if ( ffmath::isNan( k ) || ffmath::isNan( nu ) || ffmath::isNan( phi ) || ( absolute( k ) > 1.0F ) ) {
        y = ffmath::getNan();
    }
    else {
        const float n = ffmath::floor( phi/ffmath::FFP_PI + 0.5F );
        const float phi_red = phi - n*ffmath::FFP_PI;
        const float kk = k*k;
        const float s = ffmath::sin( phi_red );
        const float ss = s*s;
        const float sss = ss*s;
        const float c = ffmath::cos( phi_red );
        const float cc = c*c;
        const float tmp = 1.0F - kk*ss;
        const float pi = s*ellint_rf( cc, tmp, 1.0F )
                        + 0.333333333333F*nu*sss*ellint_rj( cc, tmp, 1.0F, 1.0F - nu*ss );

        if ( ffmath::classification::FFP_ZERO == ffmath::classify( n ) ) {
            y = pi;
        }
        else {
            y = pi + 2.0F*n*ffmath::comp_ellint_3( k, nu );
        }
    }

    return y;
}
/*============================================================================*/