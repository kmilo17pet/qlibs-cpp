#include <include/ffmath.hpp>

using namespace qlibs;

template<typename T1, typename T2>
static inline void cast_reinterpret( T1 &f, const T2 &u )
{
    (void)memcpy( &f, &u, sizeof(T2) );
}

static float getAbnormal( const int i );
static float compute_cbrt( float x , bool r );


/*============================================================================*/
static float getAbnormal( const int i )
{
    static const uint32_t u_ab[ 2 ] = { 0x7F800000u, 0x7FBFFFFFu };
    static float f_ab[ 2 ] = { 0.0F, 0.0F };
    static bool init = true;
    
    if ( init ) {
        (void)memcpy( f_ab, u_ab, sizeof(f_ab) );
        init = false;
    }
    
    return f_ab[ i ]; 
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

    if ( 0u == u ) {
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
float ffmath::absf( float x )
{
    return ( x < 0.0F ) ? -x : x;
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
        retVal = ( ( x/z ) + z ) * 0.5F;
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
    return roundf( x - 0.5F );
}
/*============================================================================*/
float ffmath::ceil( float x )
{
    return roundf( x + 0.5F );
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
    return x - trunc( x );
}
/*============================================================================*/
float ffmath::remainder( float x, float y )
{
    return x - ( y*floor( x/y ) );
}
/*============================================================================*/
float ffmath::mod( float x, float y )
{
    return ( classification::FFP_ZERO  == classify( x ) ) ? getNan() : ( x - ( y*trunc( x/y ) ) );
}
/*============================================================================*/
float ffmath::sin( float x )
{
    float y;

    x *= -ffmath::FFP_1_PI;
    y = x + 25165824.0F;
    x -= y - 25165824.0F;
    x *= ffmath::absf( x ) - 1.0F;

    return x*( ( 3.5841304553896F*ffmath::absf( x ) ) + 3.1039673861526F );
}
/*============================================================================*/
float ffmath::cos( float x )
{
    return ffmath::cos( x + ffmath::FFP_PI_2 );
}
/*============================================================================*/
float ffmath::tan( float x )
{
    float abs_x;

    x /= ffmath::absf( x ) + 1.0F;
    abs_x = ffmath::absf( x );

    return x*( ( abs_x*( ( -1.45667498914F*abs_x ) + 2.18501248371F ) ) + 0.842458832225F );
}
/*============================================================================*/
float ffmath::asin( float x )
{
    x = ffmath::sqrt( 1.0F + x ) - ffmath::sqrt( 1.0F - x );
    return x*( ( 0.131754508171F*ffmath::absf( x ) ) + 0.924391722181F );
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

    x /= ffmath::absf( x ) + 1.0F;
    abs_x = ffmath::absf( x );

    return x*( ( abs_x*( ( -1.45667498914F*abs_x ) + 2.18501248371F ) ) + 0.842458832225F );
}
/*============================================================================*/
float ffmath::atan2( float y, float x )
{
    float t, f;

    t = ffmath::FFP_PI - ( ( y < 0.0F ) ? 6.283185307f : 0.0F );
    f = ( ffmath::absf( x ) <= FLT_MIN ) ? 1.0F : 0.0F;
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
float ffmath::pow( float b, float e )
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
float ffmath::erf( float x )
{
    float retVal;

    if ( x >= 6.912F ) {
        retVal = 1.0F;
    }
    else {
        x = ffmath::exp( 3.472034176F*x );
        retVal = ( x/( ( ffmath::absf( x ) + 1.0F )*2.0F ) ) - 1.0F;
    }

    return retVal;
}
/*============================================================================*/
float ffmath::erfc( float x )
{
    return 1.0F - ffmath::erf( x );
}
/*============================================================================*/
float ffmath::rexp( float x, int32_t *pw2 )
{
    uint32_t lu = 0u, iu;
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
float ffmath::ldexp( float x, int32_t pw2 )
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
float ffmath::hypot( float x, float y )
{
    float retVal;
    
    if ( ffmath::isFinite( x ) && ffmath::isFinite( y ) ) {
        float a, b, an, bn;;
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
        retVal = ( ffmath::isInf( x ) || ffmath::isInf( y ) ) ? ffmath::getInf()
                                                              : ffmath::getNan();
    }
    
    return retVal;
}
/*============================================================================*/
float ffmath::nextAfter( float x, float y )
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
            uxi--;
        }
        else {
            uxi++;
        }
        cast_reinterpret( retVal, uxi );
    }

    return retVal;
}
/*============================================================================*/
