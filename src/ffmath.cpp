#include <include/ffmath.hpp>
/*
references :
https://github.com/akohlmey/fastermath
*/

struct limits {
    static constexpr float Min() { return FLT_MIN; }
    static constexpr float Max() { return FLT_MAX; }
    static constexpr float epsilon() { return FLT_EPSILON; }
};

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
static float poly_laguerre_recursion( size_t n,
                                      float alpha,
                                      float x );
static float poly_laguerre_large_n( size_t n,
                                    float alpha,
                                    float x );
static float poly_laguerre_hyperg( size_t n,
                                   float alpha,
                                   float x );
static float poly_legendre_p( size_t l,
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
static float expint_E1_series( float x );
static float expint_Ei( float x );
static float expint_E1_asymp( float x );
static float expint_Ei_asymp( float x );
static float expint_E1( float x );
static float expint_En_cont_frac( size_t n,
                                  float x );
static float expint_Ei_series( float x );
static float riemann_zeta_glob( float s );
static float riemann_zeta_product( float s );
static void gamma_temme( float mu,
                         float& gam1,
                         float& gam2,
                         float& gam_pl,
                         float& gam_mi);
static void bessel_jn( float nu,
                       float x,
                       float& j_nu,
                       float& n_nu,
                       float& j_pnu,
                       float& n_pnu );
static void sph_bessel_jn( size_t n,
                           float x,
                           float& j_n,
                           float& n_n,
                           float& jp_n,
                           float& np_n );
static void bessel_ik( float nu,
                       float x,
                       float& i_nu,
                       float& k_nu,
                       float& i_pnu,
                       float& k_pnu );
static float cyl_bessel_ij_series( float nu,
                                   float x,
                                   float sgn,
                                   size_t max_iter );
static void cyl_bessel_jn_asymp( float nu,
                                 float x,
                                 float& Jnu,
                                 float& Nnu );

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
float ffmath::copysign( float mag,
                        float sgn )
{
    uint32_t u_mag = 0U, u_sgn = 0U;

    cast_reinterpret( u_mag, mag );
    cast_reinterpret( u_sgn, sgn );
    u_mag &= 0x7FFFFFFFU;
    u_mag |= u_sgn & 0x80000000U;
    cast_reinterpret( mag, u_mag );

    return mag;
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
    constexpr float c13 = 1.0F/3.0F;
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
    retVal = 1.0F + ( c13*c );
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
                m = ( absolute( tmp - y ) <= limits::Min() ) ? 0.0F : tmp;
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
                m = ( absolute( tmp - y ) <= limits::Min() ) ? 0.0F : tmp;
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

    t = ffmath::FFP_PI - ( ( y < 0.0F ) ? ffmath::FFP_2PI : 0.0F );
    f = ( absolute( x ) <= limits::Min() ) ? 1.0F : 0.0F;
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
        float ip, fp;
        float ep_f = 0.0F;
        int32_t ep_i;

        ip = ffmath::floor( x + 0.5F );
        fp = x - ip;
        /*cstat -CERT-FLP34-C -MISRAC++2008-5-0-21*/
        ep_i = ( static_cast<int32_t>( ip ) + 127 ) << 23;
        /*cstat +CERT-FLP34-C +MISRAC++2008-5-0-21*/
        x = 1.535336188319500e-4F;
        x = ( x*fp ) + 1.339887440266574e-3F;
        x = ( x*fp ) + 9.618437357674640e-3F;
        x = ( x*fp ) + 5.550332471162809e-2F;
        x = ( x*fp ) + 2.402264791363012e-1F;
        x = ( x*fp ) + 6.931472028550421e-1F;
        x = ( x*fp ) + 1.0F;
        cast_reinterpret( ep_f, ep_i );
        retVal = ep_f*x;
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
        float z, px;
        int32_t ip, fp;
        int32_t val_i = 0;

        cast_reinterpret( val_i, x );
        /*cstat -MISRAC++2008-5-0-21*/
        fp = val_i & 8388607;
        ip = val_i & 2139095040;
        fp |= 1065353216;
        cast_reinterpret( x, fp );
        ip >>= 23;
        ip -= 127;
        /*cstat +MISRAC++2008-5-0-21*/
        if ( x > ffmath::FFP_SQRT2 ) {
            x *= 0.5F;
            ++ip;
        }
        x -= 1.0F;
        px = 7.0376836292e-2F;
        px = ( px*x ) - 1.1514610310e-1F;
        px = ( px*x ) + 1.1676998740e-1F;
        px = ( px*x ) - 1.2420140846e-1F;
        px = ( px*x ) + 1.4249322787e-1F;
        px = ( px*x ) - 1.6668057665e-1F;
        px = ( px*x ) + 2.0000714765e-1F;
        px = ( px*x ) - 2.4999993993e-1F;
        px = ( px*x ) + 3.3333331174e-1F;
        z = x*x;
        z = ( x*z*px ) - ( 0.5F*z ) + x;
        z *= ffmath::FFP_LOG2E;
        /*cstat -CERT-FLP36-C*/
        retVal = static_cast<float>( ip ) + z;
        /*cstat +CERT-FLP36-C*/
    }

    return retVal;
}
/*============================================================================*/
float ffmath::exp( float x )
{
    return ffmath::exp2( ffmath::FFP_LOG2E*x );
}
/*============================================================================*/
float ffmath::expm1( float x )
{
    return ffmath::exp2( ffmath::FFP_LOG2E*x ) - 1.0F;
}
/*============================================================================*/
float ffmath::exp10( float x )
{
    return ffmath::exp2( 3.32192809F*x );
}
/*============================================================================*/
float ffmath::log( float x )
{
    return ffmath::FFP_LN2*ffmath::log2( x );
}
/*============================================================================*/
float ffmath::log1p( float x )
{
    return ffmath::FFP_LN2*ffmath::log2( 1.0F + x );
}
/*============================================================================*/
float ffmath::log10( float x )
{
    return ffmath::FFP_LOG10_2*ffmath::log2( x );
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
    const float epx = ffmath::exp( x );
    const float enx = 1.0F/epx;

    return 0.5F*( epx - enx );
}
/*============================================================================*/
float ffmath::cosh( float x )
{
    const float epx = ffmath::exp( x );
    const float enx = 1.0F/epx;

    return 0.5F*( epx + enx );
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
    return ffmath::mod( x + ffmath::FFP_PI, ffmath::FFP_2PI ) - ffmath::FFP_PI;
}
/*============================================================================*/
float ffmath::wrapTo2Pi( float x )
{
    return ffmath::mod( x, ffmath::FFP_2PI );
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
    constexpr float lo = 2.0F*limits::Min();
    constexpr float hi = 0.5F*limits::Max();
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
        if ( y < limits::epsilon() ) {
            if ( y >= limits::Min() ) {
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
    constexpr float pnt68 = 0.6796875F;
    float result;

    if ( x > 171.624F ) {
        result = ffmath::getInf();
    }
    else {
        float y, corrector, num, den;

        y = x;
        if ( y <= limits::epsilon() ) {
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
                num = ( num*xMinus ) + 4.945235359296727046734888e+0F;
                num = ( num*xMinus ) + 2.018112620856775083915565e+2F;
                num = ( num*xMinus ) + 2.290838373831346393026739e+3F;
                num = ( num*xMinus ) + 1.131967205903380828685045e+4F;
                num = ( num*xMinus ) + 2.855724635671635335736389e+4F;
                num = ( num*xMinus ) + 3.848496228443793359990269e+4F;
                num = ( num*xMinus ) + 2.637748787624195437963534e+4F;
                num = ( num*xMinus ) + 7.225813979700288197698961e+3F;
                den = ( den*xMinus ) + 6.748212550303777196073036e+1F;
                den = ( den*xMinus ) + 1.113332393857199323513008e+3F;
                den = ( den*xMinus ) + 7.738757056935398733233834e+3F;
                den = ( den*xMinus ) + 2.763987074403340708898585e+4F;
                den = ( den*xMinus ) + 5.499310206226157329794414e+4F;
                den = ( den*xMinus ) + 6.161122180066002127833352e+4F;
                den = ( den*xMinus ) + 3.635127591501940507276287e+4F;
                den = ( den*xMinus ) + 8.785536302431013170870835e+3F;
                result = corrector + ( xMinus*( d1 + xMinus*( num/den ) ) );
            }
            else {
                xMinus = ( y - 0.5F ) - 0.5F;
                den = 1.0F;
                num = 0.0F;
                num = ( num*xMinus ) + 4.974607845568932035012064e+0F;
                num = ( num*xMinus ) + 5.424138599891070494101986e+2F;
                num = ( num*xMinus ) + 1.550693864978364947665077e+4F;
                num = ( num*xMinus ) + 1.847932904445632425417223e+5F;
                num = ( num*xMinus ) + 1.088204769468828767498470e+6F;
                num = ( num*xMinus ) + 3.338152967987029735917223e+6F;
                num = ( num*xMinus ) + 5.106661678927352456275255e+6F;
                num = ( num*xMinus ) + 3.074109054850539556250927e+6F;
                den = ( den*xMinus ) + 1.830328399370592604055942e+2F;
                den = ( den*xMinus ) + 7.765049321445005871323047e+3F;
                den = ( den*xMinus ) + 1.331903827966074194402448e+5F;
                den = ( den*xMinus ) + 1.136705821321969608938755e+6F;
                den = ( den*xMinus ) + 5.267964117437946917577538e+6F;
                den = ( den*xMinus ) + 1.346701454311101692290052e+7F;
                den = ( den*xMinus ) + 1.782736530353274213975932e+7F;
                den = ( den*xMinus ) + 9.533095591844353613395747e+6F;
                result = corrector + ( xMinus*( d2 + xMinus*( num/den ) ) );
            }
        }
        else if ( y <= 4.0F ) {
            const float xMinus = y - 2.0F;
            den = 1.0F;
            num = 0.0F;
            num = ( num*xMinus ) + 4.974607845568932035012064e+0F;
            num = ( num*xMinus ) + 5.424138599891070494101986e+2F;
            num = ( num*xMinus ) + 1.550693864978364947665077e+4F;
            num = ( num*xMinus ) + 1.847932904445632425417223e+5F;
            num = ( num*xMinus ) + 1.088204769468828767498470e+6F;
            num = ( num*xMinus ) + 3.338152967987029735917223e+6F;
            num = ( num*xMinus ) + 5.106661678927352456275255e+6F;
            num = ( num*xMinus ) + 3.074109054850539556250927e+6F;
            den = ( den*xMinus ) + 1.830328399370592604055942e+2F;
            den = ( den*xMinus ) + 7.765049321445005871323047e+3F;
            den = ( den*xMinus ) + 1.331903827966074194402448e+5F;
            den = ( den*xMinus ) + 1.136705821321969608938755e+6F;
            den = ( den*xMinus ) + 5.267964117437946917577538e+6F;
            den = ( den*xMinus ) + 1.346701454311101692290052e+7F;
            den = ( den*xMinus ) + 1.782736530353274213975932e+7F;
            den = ( den*xMinus ) + 9.533095591844353613395747e+6F;
            result = xMinus*( d2 + xMinus*( num/den ) );
        }
        else if ( y <= 12.0F ) {
            const float xMinus = y - 4.0F;
            den = -1.0F;
            num = 0.0F;
            num = ( num*xMinus ) + 1.474502166059939948905062e+04F;
            num = ( num*xMinus ) + 2.426813369486704502836312e+06F;
            num = ( num*xMinus ) + 1.214755574045093227939592e+08F;
            num = ( num*xMinus ) + 2.663432449630976949898078e+09F;
            num = ( num*xMinus ) + 2.940378956634553899906876e+10F;
            num = ( num*xMinus ) + 1.702665737765398868392998e+11F;
            num = ( num*xMinus ) + 4.926125793377430887588120e+11F;
            num = ( num*xMinus ) + 5.606251856223951465078242e+11F;
            den = ( den*xMinus ) + 2.690530175870899333379843e+03F;
            den = ( den*xMinus ) + 6.393885654300092398984238e+05F;
            den = ( den*xMinus ) + 4.135599930241388052042842e+07F;
            den = ( den*xMinus ) + 1.120872109616147941376570e+09F;
            den = ( den*xMinus ) + 1.488613728678813811542398e+10F;
            den = ( den*xMinus ) + 1.016803586272438228077304e+11F;
            den = ( den*xMinus ) + 3.417476345507377132798597e+11F;
            den = ( den*xMinus ) + 4.463158187419713286462081e+11F;
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
    static const float ft[ 35 ] = { 1.0F, 1.0F, 2.0F, 6.0F, 24.0F, 120.0F, 720.0F,
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
        y = ft[ static_cast<size_t>( x ) ];
    }
    else {
        y = ffmath::getNan();
    }

    return y;
}
/*============================================================================*/
static float poly_laguerre_recursion( size_t n,
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
static float poly_laguerre_large_n( size_t n,
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
    const float th2 = 2.0F*th;
    const float serTerm2 = ffmath::sin( ( 0.25F*eta*( ( th2 ) - ffmath::sin( th2 ) ) + ffmath::FFP_PI_4 ) );

    return ffmath::exp( lnPre )*( serTerm1 + serTerm2 );
}
/*============================================================================*/
static float poly_laguerre_hyperg( size_t n,
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
float ffmath::assoc_laguerre( size_t n,
                              size_t m,
                              float x )
{
    // include/tr1/poly_laguerre.tcc
    float y;
    /*cstat -CERT-FLP36-C*/
    const float alpha = static_cast<float>( m );
    const float N = static_cast<float>( n );
    /*cstat +CERT-FLP36-C*/
    if ( ( x < 0.0F ) || ffmath::isNan( x ) ) {
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
    else if ( ( alpha > -1.0F ) && ( x < ( ( 2.0F*( alpha + 1.0F ) ) + ( 4.0F*N ) ) ) ) {
        y = poly_laguerre_large_n( n, alpha, x );
    }
    else if ( ( ( x > 0.0F ) && ( alpha < -( N + 1.0F ) ) ) ) {
        y = poly_laguerre_recursion( n, alpha, x );
    }
    else {
        y = poly_laguerre_hyperg( n, alpha, x );
    }

    return y;
}
/*============================================================================*/
static float poly_legendre_p( size_t l,
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
float ffmath::assoc_legendre( size_t n,
                              size_t m,
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
    constexpr float loLim = 5.0F*limits::Min();
    float result;


    if ( ( x < 0.0F ) || ( y < 0.0F ) || ( z < 0.0F ) ||
        ( ( x + y ) < loLim ) || ( ( x + z ) < loLim ) || ( ( y + z) < loLim )
       ) {
        result = ffmath::getNan();
    }
    else {
        constexpr float c0 = 1.0F/4.0F;
        constexpr float c1 = 1.0F/24.0F;
        constexpr float c2 = 1.0F/10.0F;
        constexpr float c3 = 3.0F/44.0F;
        constexpr float c4 = 1.0F/14.0F;
        constexpr float errTol = 0.0024607833005759250678823324F;
        constexpr float c13 = 1.0F/3.0F;

        float xn = x;
        float yn = y;
        float zn = z;
        float mu, xnDev, ynDev, znDev;
        const size_t maxIter = 100U;
        float e2, e3, s;

        for ( size_t iter = 0U ; iter < maxIter ; ++iter ) {
            float abs_xnDev, abs_ynDev, abs_znDev, lambda, epsilon;
            float xRoot, yRoot, zRoot;

            mu = ( xn + yn + zn )*c13;
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

    if ( ( x < 0.0F ) || ( y < 0.0F ) || ( ( x + y ) < loLim ) || ( z < loLim ) ) {
        result = ffmath::getNan();
    }
    else {
        constexpr float c0 = 1.0F/4.0F;
        constexpr float c1 = 3.0F/14.0F;
        constexpr float c2 = 1.0F/6.0F;
        constexpr float c3 = 9.0F/22.0F;
        constexpr float c4 = 3.0F/26.0F;
        float xn = x;
        float yn = y;
        float zn = z;
        float sigma = 0.0F;
        float power4 = 1.0F;
        float mu, xnDev, ynDev, znDev;
        float ea, eb, ec, ed, ef, s1, s2;
        const size_t maxIter = 100U;

        for ( size_t iter = 0U ; iter < maxIter ; ++iter ) {
            float abs_xnDev, abs_ynDev, abs_znDev, lambda, epsilon;
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
        constexpr float c0 = 1.0F/4.0F;
        constexpr float c1 = 1.0F/7.0F;
        constexpr float c2 = 9.0F/22.0F;
        constexpr float c3 = 3.0F/10.0F;
        constexpr float c4 = 3.0F/8.0F;
        constexpr float c13 = 1.0F/3.0F;
        float xn = x;
        float yn = y;
        const size_t maxIter = 100;
        float mu, s, sn;

        for ( size_t iter = 0U ; iter < maxIter ; ++iter ) {
            float lambda;


            mu = ( xn + 2.0F*yn )*c13;
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

    if ( ( x < 0.0F ) || ( y < 0.0F ) || ( z < 0.0F ) || ( ( x + y ) < loLim ) ||
         ( ( x + z ) < loLim ) || ( ( y + z) < loLim )|| ( p < loLim ) ) {
        result = ffmath::getNan();
    }
    else {
        constexpr float c0 = 1.0F/4.0F;
        constexpr float c1 = 3.0F/14.0F;
        constexpr float c2 = 1.0F/3.0F;
        constexpr float c3 = 3.0F/22.0F;
        constexpr float c4 = 3.0F/26.0F;
        float xn = x;
        float yn = y;
        float zn = z;
        float pn = p;
        float sigma = 0.0F;
        float power4 = 1.0F;
        float mu, xnDev, ynDev, znDev, pnDev;
        float ea, eb, ec, e2, e3, s1, s2, s3;
        const size_t maxIter = 100;

        for ( size_t iter = 0U ; iter < maxIter ; ++iter ) {
            float abs_xnDev, abs_ynDev, abs_znDev, abs_pnDev, lambda, alpha1;
            float alpha2, beta, epsilon;
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
        constexpr float c13 = 1.0F/3.0F;

        y = ellint_rf( 0.0F, one_m_kk, 1.0F ) - c13*kk*ellint_rd( 0.0F, one_m_kk, 1.0F );
    }

    return y;
}
/*============================================================================*/
float ffmath::comp_ellint_3( float k,
                             float nu )
{
    float y;
    const float abs_k = absolute( k );

    if ( ffmath::isNan( k ) || ffmath::isNan( nu ) || ( abs_k > 1.0F ) ) {
        y = ffmath::getNan();
    }
    else if ( ffmath::isEqual( 1.0F, nu ) ) {
        y = ffmath::getInf();
    }
    else {
        const float kk = k*k;
        const float one_m_kk = 1.0F - kk;
        constexpr float c13 = 1.0F/3.0F;

        y = ellint_rf( 0.0F, one_m_kk, 1.0F ) + c13*nu*ellint_rj( 0.0F, one_m_kk, 1.0F, 1.0F - nu );
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
            y = f + ( 2.0F*n*ffmath::comp_ellint_1( k ) );
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
        constexpr float c13 = 1.0F/3.0F;
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
                        - c13*kk*sss*ellint_rd( cc, tmp, 1.0F );

        if ( ffmath::classification::FFP_ZERO == ffmath::classify( n ) ) {
            y = e;
        }
        else {
            y = e + ( 2.0F*n*ffmath::comp_ellint_2( k ) );
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
        constexpr float c13 = 1.0F/3.0F;
        const float pi = s*ellint_rf( cc, tmp, 1.0F ) + c13*nu*sss*ellint_rj( cc, tmp, 1.0F, 1.0F - nu*ss );

        if ( ffmath::classification::FFP_ZERO == ffmath::classify( n ) ) {
            y = pi;
        }
        else {
            y = pi + ( 2.0F*n*ffmath::comp_ellint_3( k, nu ) );
        }
    }

    return y;
}
/*============================================================================*/
static float expint_E1_series( float x )
{
    float term = 1.0F;
    float eSum = 0.0F;
    float oSum = 0.0F;
    const size_t maxIter = 1000U;

    for ( size_t i = 1U; i < maxIter; ++i ) {
        /*cstat -CERT-FLP36-C*/
        const auto j = static_cast<float>( i );
        /*cstat +CERT-FLP36-C*/
        term *= -x/j;
        if ( absolute( term ) < limits::epsilon() ) {
            break;
        }
        if ( term >= 0.0F ) {
            eSum += term/j;
        }
        else {
            oSum += term/j;
        }
    }
    return - eSum - oSum - ffmath::FFP_GAMMA_E - ffmath::log( x );
}
/*============================================================================*/
static float expint_E1_asymp( float x )
{
    float term = 1.0F;
    float eSum = 1.0F;
    float oSum = 0.0F;
    const size_t maxIter = 1000U;

    for ( size_t i = 1U; i < maxIter; ++i ) {
        const float prev = term;
        /*cstat -CERT-FLP36-C*/
        term *= -static_cast<float>( i )/x;
        /*cstat +CERT-FLP36-C*/
        if ( absolute( term ) > absolute( prev ) ) {
            break;
        }
        if ( term >= 0.0F ) {
            eSum += term;
        }
        else {
            oSum += term;
        }
    }
    return ffmath::exp( -x )*( eSum + oSum )/x;
}
/*============================================================================*/
static float expint_Ei_asymp( float x )
{
    float term = 1.0F;
    float sum = 1.0F;
    const size_t maxIter = 1000U;

    for ( size_t i = 1U; i < maxIter; ++i ) {
        const float prev = term;
        /*cstat -CERT-FLP36-C*/
        term *= -static_cast<float>( i )/x;
        /*cstat +CERT-FLP36-C*/
        if ( ( term < limits::epsilon() ) || ( term >= prev ) ) {
            break;
        }
        sum += term;
    }

    return ffmath::exp( x )*sum/x;
}
/*============================================================================*/
static float expint_E1( float x )
{
    float y;

    if ( x < 0.0F ) {
        y = -expint_Ei( -x );
    }
    else if ( x < 1.0F ) {
        y = expint_E1_series( x );
    }
    else if ( x < 100.F ) {
        y = expint_En_cont_frac( 1, x );
    }
    else {
        y = expint_E1_asymp( x );
    }

    return y;
}
/*============================================================================*/
static float expint_En_cont_frac( size_t n,
                                  float x )
{
    float y = ffmath::getNan();
    const int maxIter = 1000;
    const int nm1 = static_cast<int>( n ) - 1;
    /*cstat -CERT-FLP36-C*/
    float b = x + static_cast<float>( n );
    /*cstat +CERT-FLP36-C*/
    float c = 1.0F/limits::Min() ;
    float d = 1.0F/b;
    float h = d;

    for ( int i = 1; i <= maxIter; ++i ) {
        /*cstat -MISRAC++2008-5-0-7 -CERT-FLP36-C*/
        const float a = -static_cast<float>( i*( nm1 + i ) );
        /*cstat MISRAC++2008-5-0-7 +CERT-FLP36-C*/
        b += 2.0F;
        d = 1.0F/( ( a*d ) + b );
        c = b + ( a/c );
        const float del = c*d;
        h *= del;
        if ( absolute( del - 1.0F ) < limits::epsilon() ) {
            y = h*ffmath::exp( -x );
            break;
        }
    }

    return y;
}
/*============================================================================*/
static float expint_Ei_series( float x )
{
    float term = 1.0F;
    float sum = 0.0F;
    const size_t maxIter = 1000U;

    for ( size_t i = 1U; i < maxIter; ++i ) {
        /*cstat -CERT-FLP36-C*/
        const auto j = static_cast<float>( i );
        /*cstat +CERT-FLP36-C*/
        term *= x/j;
        sum += term/j;
        if ( term < ( limits::epsilon()*sum ) ) {
            break;
        }
      }
    return ffmath::FFP_GAMMA_E + sum + ffmath::log( x );
}
/*============================================================================*/
static float expint_Ei( float x )
{
    constexpr float logEps = 36.044F;
    float y;

    if ( x < 0.0F ) {
        y = -expint_E1( -x );
    }
    else if ( x < logEps ) {
        y = expint_Ei_series( x );
    }
    else {
        y = expint_Ei_asymp( x );
    }

    return y;
}
/*============================================================================*/
float ffmath::expint( float num )
{
    return ( ffmath::isNan( num ) ) ? ffmath::getNan() : expint_Ei( num );
}
/*============================================================================*/
float ffmath::hermite( size_t n,
                       float x )
{
    float y = 0.0F;

    if ( ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else {
        const float H_0 = 1.0F;

        if ( 0U == n ) {
            y= H_0;
        }
        else {
            const float H_1 = 2.0F*x;
            if ( 1U == n ) {
                y = H_1;
            }
            else {
                float H_nm1, H_nm2;

                H_nm2 = H_0;
                H_nm1 = H_1;
                for ( size_t i = 2U; i <= n; ++i ) {
                    /*cstat -CERT-FLP36-C*/
                    const auto j = static_cast<float>( i - 1U );
                    /*cstat +CERT-FLP36-C*/
                    y = 2.0F*( ( x*H_nm1 ) - ( j*H_nm2 ) );
                    H_nm2 = H_nm1;
                    H_nm1 = y;
                }
            }
        }
    }

    return y;
}
/*============================================================================*/
float ffmath::laguerre( size_t n,
                        float x )
{
    return assoc_laguerre( n, 0, x );
}
/*============================================================================*/
float ffmath::legendre( size_t n,
                        float x )
{
    float y = 0.0F;

    if ( ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else if ( ffmath::isEqual( 1.0F, x ) ) {
        y = 1.0F;
    }
    else if ( ffmath::isEqual( -1.0F, x ) ) {
        y =  ( ( n % 2 ) == 1 ) ? -1.0F : 1.0F;
    }
    else {
        float p_lm2 = 1.0F;

        if ( 0 == n ) {
            y = p_lm2;
        }
        else {
            float p_lm1 = x;

            if ( 1 == n ) {
                y = p_lm1;
            }
            else {
                for ( size_t ll = 2U ; ll <= n; ++ll ) {
                    /*cstat -CERT-FLP36-C*/
                    const auto ll_f = static_cast<float>( ll );
                    /*cstat +CERT-FLP36-C*/
                    y = ( 2.0F*x*p_lm1 ) - p_lm2 - ( x*p_lm1 - p_lm2)/ll_f;
                    p_lm2 = p_lm1;
                    p_lm1 = y;
                }
            }
        }
    }

    return y;
}
/*============================================================================*/
static float riemann_zeta_glob( float s )
{
    constexpr float maxBinCoeff = 86.4982335337F;
    float zeta = 0.0F;
    const float ss = s;
    bool neg = false;

    if ( s < 0.0F ) {
        s = 1.0F - s;
        neg = true;
    }
    float num = 0.5F;
    const size_t maxIt = 10000U;
    /*cstat -MISRAC++2008-6-6-4*/
    for ( size_t i = 0U ; i < maxIt; ++i ) {
        bool punt = false;
        float sgn = 1.0F;
        float term = 0.0F;
        for ( size_t j = 0U ; j <= i ; ++j ) {
            /*cstat -CERT-FLP36-C*/
            const float ii = static_cast<float>( i );
            const float jj = static_cast<float>( j );
            /*cstat +CERT-FLP36-C*/
            float bin_coeff = ffmath::lgamma( 1.0F + ii ) -
                              ffmath::lgamma( 1.0F + jj ) -
                              ffmath::lgamma( 1.0F + ii - jj );

            if ( bin_coeff > maxBinCoeff ) {
                punt = true;
                break;
            }
            bin_coeff = ffmath::exp( bin_coeff );
            term += sgn*bin_coeff*ffmath::pow( 1.0F + jj, -s );
            sgn *= -1.0F;
        }
        if ( punt ) {
            break;
        }
        term *= num;
        zeta += term;
        if ( absolute( term/zeta ) < limits::epsilon() ) {
            break;
        }
        num *= 0.5F;
    }
    /*cstat +MISRAC++2008-6-6-4*/
    zeta /= 1.0F - ffmath::pow( 2.0F, 1.0F - s );

    if ( neg ) {
        zeta *= ffmath::pow( 2.0F*ffmath::FFP_PI, ss )*
                ffmath::sin( ffmath::FFP_PI_2*ss )*
                ffmath::exp( ffmath::lgamma( s ) )/ffmath::FFP_PI;
    }
    return zeta;
}
/*============================================================================*/
static float riemann_zeta_product( float s )
{
    static const uint8_t prime[ 29 ] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
        71, 73, 79, 83, 89, 97, 101, 103, 107, 109
    };
    static constexpr size_t num_primes = sizeof(prime)/sizeof(float);
    float zeta = 1.0F;

    for ( size_t i = 0U ; i < num_primes; ++i ) {
        const float f = 1.0F - ffmath::pow( static_cast<float>( prime[ i ] ), -s );

        zeta *= f;
        if ( ( 1.0F - f ) < limits::epsilon() ) {
            break;
        }
    }
    zeta = 1.0F/zeta;

    return zeta;
}
/*============================================================================*/
float ffmath::riemann_zeta( float s )
{
    float z;

    if ( ffmath::isNan( s ) ) {
        z = ffmath::getNan();
    }
    else if ( ffmath::isEqual( 1.0F, s ) ) {
        z = ffmath::getInf();
    }
    else if ( s < -19.0F ) {
        z = riemann_zeta_product( 1.0F - s );
        z *= ffmath::pow( 2.0F*ffmath::FFP_PI, s )*
             ffmath::sin( ffmath::FFP_PI_2*s )*
             ffmath::exp( ffmath::lgamma( 1.0F - s ) )/
             ffmath::FFP_PI;
    }
    else if ( s < 20.0F ) {
        z = riemann_zeta_glob( s );
    }
    else {
        z = riemann_zeta_product( s );
    }

    return z;
}
/*============================================================================*/
static void gamma_temme( float mu,
                         float& gam1,
                         float& gam2,
                         float& gam_pl,
                         float& gam_mi)
{

    gam_pl = 1.0F/ffmath::tgamma( 1.0F + mu );
    gam_mi = 1.0F/ffmath::tgamma( 1.0F - mu );

    if ( absolute( mu ) < limits::epsilon() ) {
        gam1 = -ffmath::FFP_GAMMA_E;
    }
    else {
        gam1 = ( gam_mi - gam_pl )/( 2.0F*mu );
    }
    gam2 = 0.5F*( gam_mi + gam_pl );
}
/*============================================================================*/
static void bessel_jn( float nu,
                       float x,
                       float& j_nu,
                       float& n_nu,
                       float& j_pnu,
                       float& n_pnu )
{
    if ( ffmath::isEqual( 0.0F, x ) ) {
        if ( ffmath::isEqual( 0.0F, nu ) ) {
            j_nu = 1.0F;
            j_pnu = 0.0F;
        }
        else if ( ffmath::isEqual( 1.0F, nu )  ) {
            j_nu = 0.0F;
            j_pnu = 0.5F;
        }
        else {
            j_nu = 0.0F;
            j_pnu = 0.0F;
        }
        n_nu = -ffmath::getInf();
        n_pnu = ffmath::getInf();
    }
    else {
        const float eps = limits::epsilon();
        const float fp_min = 1.08420217256745973463717809E-19F;
        const int max_iter = 15000;
        const float x_min = 2.0F;
        /*cstat -CERT-FLP34-C*/
        const int nl = ( x < x_min ) ? static_cast<int>( nu + 0.5F )
                                     : ffmath::Max( 0, static_cast<int>( nu - x + 1.5F ) ) ;
        /*cstat +CERT-FLP34-C -CERT-FLP36-C*/
        const float mu = nu - static_cast<float>( nl );
        /*cstat +CERT-FLP36-C*/
        const float mu2 = mu*mu;
        const float xi = 1.0F/x;
        const float xi2 = 2.0F*xi;
        const float w = xi2/ffmath::FFP_PI;
        float iSign = 1.0F;
        float h = nu*xi;
        if ( h < fp_min ) {
            h = fp_min;
        }
        float b = xi2*nu;
        float d = 0.0F;
        float c = h;
        for ( int i = 1; i <= max_iter; ++i ) {
            b += xi2;
            d = b - d;
            if ( absolute( d ) < fp_min ) {
                d = fp_min;
            }
            c = b - ( 1.0F/c );
            if ( absolute( c ) < fp_min ) {
                c = fp_min;
            }
            d = 1.0F / d;

            const float del = c * d;
            h *= del;
            if ( d < 0.0F ) {
                iSign = -iSign;
            }
            if ( absolute( del - 1.0F ) < eps ) {
                break;
            }
        }
        float j_nul = iSign*fp_min;
        float j_pnu_l = h*j_nul;
        const float j_nul1 = j_nul;
        const float j_pnu1 = j_pnu_l;
        float fact = nu*xi;

        for ( int l = nl; l >= 1; --l ) {
            const float j_nu_temp = ( fact*j_nul ) + j_pnu_l;

            fact -= xi;
            j_pnu_l = ( fact*j_nu_temp ) - j_nul;
            j_nul = j_nu_temp;
        }
        if ( ffmath::isEqual( 0.0F, j_nul ) ) {
            j_nul = eps;
        }
        const float f = j_pnu_l/j_nul;
        float n_mu, n_nu1, n_pmu, Jmu;
        if ( x < x_min ) {
            const float x2 = 0.5F*x;
            const float pi_mu = ffmath::FFP_PI*mu;
            const float fact_l = ( absolute( pi_mu ) < eps ) ? 1.0F : pi_mu/ffmath::sin( pi_mu );
            d = -ffmath::log( x2 );
            float e = mu*d;
            const float fact2 = ( absolute( e ) < eps ) ? 1.0F : ffmath::sinh(e)/e;
            float gam1, gam2, gam_pl, gam_mi;

            gamma_temme( mu, gam1, gam2, gam_pl, gam_mi );
            float ff = ( 2.0F/ffmath::FFP_PI )*fact_l*( gam1*ffmath::cosh( e ) + ( gam2*fact2*d ) );
            e = ffmath::exp( e );
            float p = e/( ffmath::FFP_PI*gam_pl );
            float q = 1.0F/( e*ffmath::FFP_PI*gam_mi );
            const float pi_mu2 = pi_mu / 2.0F;
            const float fact3 = ( absolute( pi_mu2 ) < eps ) ? 1.0F : ffmath::sin( pi_mu2 )/pi_mu2;
            const float r = ffmath::FFP_PI*pi_mu2*fact3*fact3;
            float sum = ff + ( r*q );
            float sum1 = p;
            c = 1.0F;
            d = -x2*x2;
            for ( int i = 1; i <= max_iter ; ++i ) {
                /*cstat -CERT-FLP36-C*/
                const float j = static_cast<float>( i );
                /*cstat +CERT-FLP36-C*/
                ff = ( j*ff + p + q )/( j*j - mu2 );
                c *= d/j;
                p /= j - mu;
                q /= j + mu;
                const float del = c*( ff + ( r*q ) );
                sum += del;
                const float del1 = ( c*p ) - ( j*del );
                sum1 += del1;
                if ( absolute( del ) < ( eps*( 1.0F + absolute( sum ) ) ) ) {
                    break;
                }
            }

            n_mu = -sum;
            n_nu1 = -sum1*xi2;
            n_pmu = ( mu*xi*n_mu ) - n_nu1;
            Jmu = w/( n_pmu - ( f*n_mu ) );
        }
        else {
            float a = 0.25F - mu2;
            float q = 1.0F;
            float p = -xi*0.5F;
            const float br = 2.0F*x;
            float bi = 2.0F;
            float fact_g = a*xi/( ( p*p ) + ( q*q ) );
            float cr = br + ( q*fact_g );
            float ci = bi + ( p*fact_g );
            float den = ( br*br ) + ( bi*bi );
            float dr = br/den;
            float di = -bi/den;
            float dlr = ( cr*dr ) - ( ci*di );
            float dli = ( cr*di ) + ( ci*dr );
            float temp = ( p*dlr ) - ( q*dli );
            q = ( p*dli ) + ( q * dlr );
            p = temp;

            for ( int i = 2; i <= max_iter; ++i ) {
                a += static_cast<float>( 2*( i - 1 ) );
                bi += 2.0F;
                dr = ( a*dr ) + br;
                di = ( a*di ) + bi;
                if ( ( absolute( dr ) + absolute( di ) ) < fp_min ) {
                    dr = fp_min;
                }
                fact_g = a/( ( cr*cr ) + ( ci*ci ) );
                cr = br + ( cr*fact_g );
                ci = bi - ( ci*fact_g );
                if ( ( absolute( cr ) + absolute( ci ) ) < fp_min ) {
                    cr = fp_min;
                }
                den = ( dr*dr ) + ( di*di );
                dr /= den;
                di /= -den;
                dlr = ( cr*dr ) - ( ci*di );
                dli = ( cr*di ) + ( ci*dr );
                temp = ( p*dlr ) - ( q*dli );
                q = ( p*dli ) + ( q*dlr );
                p = temp;
                if ( absolute( dlr - 1.0F ) + absolute( dli ) < eps ) {
                    break;
                }
            }
            const float gam = ( p - f )/q;

            Jmu = ffmath::sqrt( w/( ( p - f )*gam + q ) );
            if ( ( Jmu*j_nul ) < 0.0F) {
                Jmu = -Jmu;
            }
            n_mu = gam*Jmu;
            n_pmu = ( p + ( q/gam ) )*n_mu;
            n_nu1 = ( mu*xi*n_mu ) - n_pmu;
        }
        fact = Jmu/j_nul;
        j_nu = fact*j_nul1;
        j_pnu = fact*j_pnu1;
        for ( int i = 1; i <= nl; ++i ) {
            /*cstat -CERT-FLP36-C*/
            const float n_nu_temp = ( mu + static_cast<float>( i ) )*xi2*n_nu1 - n_mu;
            /*cstat +CERT-FLP36-C*/
            n_mu = n_nu1;
            n_nu1 = n_nu_temp;
        }
        n_nu = n_mu;
        n_pnu = ( nu*xi*n_mu ) - n_nu1;
    }
}
/*============================================================================*/
static void sph_bessel_jn( size_t n,
                           float x,
                           float& j_n,
                           float& n_n,
                           float& jp_n,
                           float& np_n )
{
    /*cstat -CERT-FLP36-C*/
    const float nu = static_cast<float>( n ) + 0.5F;
    /*cstat +CERT-FLP36-C*/
    constexpr float sqrtpi2 = 1.25331413731550012080617761967005208134651184F;
    float j_nu, n_nu, jp_nu, np_nu;
    const float factor = sqrtpi2*ffmath::rSqrt( x );
    const float inv_2x = 1.0F/( 2.0F*x );

    bessel_jn( nu, x, j_nu, n_nu, jp_nu, np_nu );
    j_n = factor*j_nu;
    n_n = factor*n_nu;
    jp_n = ( factor*jp_nu ) - ( j_n*inv_2x );
    np_n = ( factor*np_nu ) - ( n_n*inv_2x );
}
/*============================================================================*/
float ffmath::sph_bessel( size_t n,
                          float x )
{
    float y;

    if ( ( x < 0.0F ) || ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else if ( ffmath::isEqual( 0.0F, x ) ) {
        y = ( 0U == n ) ? 1.0F : 0.0F;
    }
    else {
        float j_n, n_n, jp_n, np_n;

        sph_bessel_jn( n, x, j_n, n_n, jp_n, np_n );
        y = j_n;
    }

    return y;
}
/*============================================================================*/
float ffmath::sph_neumann( size_t n,
                           float x )
{
    float y;

    if ( ( x < 0.0F ) || ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else if ( ffmath::isEqual( 0.0F, x ) ) {
        y = -ffmath::getInf();
    }
    else {
        float j_n, n_n, jp_n, np_n;

        sph_bessel_jn( n, x, j_n, n_n, jp_n, np_n );
        y = n_n;
    }
    return y;
}
/*============================================================================*/
static void bessel_ik( float nu,
                       float x,
                       float& i_nu,
                       float& k_nu,
                       float& i_pnu,
                       float& k_pnu )
{
    if ( ffmath::isEqual( 0.0F, x ) ) {
        if ( ffmath::isEqual( 0.0F, nu ) ) {
            i_nu = 1.0F;
            i_pnu = 0.0F;
        }
        else if ( ffmath::isEqual( 1.0F, x ) ) {
            i_nu = 0.0F;
            i_pnu = 0.5F;
        }
        else {
            i_nu = 0.0F;
            i_pnu = 0.0F;
        }
        k_nu = ffmath::getInf();
        k_pnu = -ffmath::getInf();
    }
    else {
        const float eps = limits::epsilon();
        const float fp_min = 10.0F*limits::epsilon();
        const int max_iter = 15000;
        const float x_min = 2.0F;
        /*cstat -CERT-FLP34-C -CERT-FLP36-C*/
        const int nl = static_cast<int>( nu + 0.5F );
        const float mu = nu - static_cast<float>( nl );
        /*cstat +CERT-FLP34-C +CERT-FLP36-C*/
        const float mu2 = mu*mu;
        const float xi = 1.0F/x;
        const float xi2 = 2.0F*xi;
        float h = nu*xi;

        if ( h < fp_min ) {
            h = fp_min;
        }
        float b = xi2*nu;
        float d = 0.0F;
        float c = h;

        for ( int i = 1; i <= max_iter; ++i ) {
            b += xi2;
            d = 1.0F/( b + d );
            c = b + ( 1.0F/c );
            const float del = c*d;
            h *= del;
            if ( absolute( del - 1.0F ) < eps ) {
                break;
            }
        }

        float i_nul = fp_min;
        float i_pnu_l = h*i_nul;
        const float i_nul1 = i_nul;
        const float i_pnu_1 = i_pnu_l;
        float fact_m = nu*xi;

        for ( int l = nl ; l >= 1 ; --l ) {
            const float i_nu_temp = ( fact_m*i_nul )  + i_pnu_l;

            fact_m -= xi;
            i_pnu_l = ( fact_m*i_nu_temp ) + i_nul;
            i_nul = i_nu_temp;
        }
        const float f = i_pnu_l/i_nul;
        float Kmu, k_nu1;

        if ( x < x_min ) {
            const float x2 = 0.5F*x;
            const float pi_mu = ffmath::FFP_PI*mu;
            const float fact = ( absolute( pi_mu ) < eps ) ? 1.0F : pi_mu/ffmath::sin( pi_mu );
            d = -ffmath::log( x2 );
            float e = mu*d;
            const float fact2 = ( absolute( e ) < eps ) ? 1.0F : ffmath::sinh( e )/e ;
            float gam1, gam2, gam_pl, gam_mi;
            gamma_temme(mu, gam1, gam2, gam_pl, gam_mi);
            float ff = fact*( ( gam1*ffmath::cosh( e ) ) + ( gam2*fact2*d ) );
            float sum = ff;
            e = ffmath::exp( e );
            float p = e/( 2.0F*gam_pl );
            float q = 1.0F/( 2.0F*e*gam_mi );
            float sum1 = p;
            c = 1.0F;
            d = x2*x2;

            for ( int i = 1; i <= max_iter; ++i ) {
                const float j = static_cast<float>( i );
                ff = ( j*ff + p + q )/( j*j - mu2 );
                c *= d/j;
                p /= j - mu;
                q /= j + mu;
                const float del = c*ff;
                sum += del;
                sum1 += c*( p - ( j*ff ) );
                if ( absolute( del ) < ( eps*absolute( sum ) ) ) {
                    break;
                }
            }

            Kmu = sum;
            k_nu1 = sum1*xi2;
        }
        else {
            float del_h = d;
            float q1 = 0.0F;
            float q2 = 1.0F;
            const float a1 = 0.25F - mu2;
            float q = a1;
            float a = -a1;
            float s = 1.0F + ( q*del_h );

            b = 2.0F*( 1.0F + x );
            d = 1.0F/b;
            h = d;
            c = a1;
            for ( int i = 2 ; i <= max_iter; ++i) {
                a -= static_cast<float>( 2*( i - 1 ) );
                c = -a*c/static_cast<float>( i );
                const float q_new = ( q1 - ( b*q2 ) )/a;
                q1 = q2;
                q2 = q_new;
                q += c*q_new;
                b += 2.0F;
                d = 1.0F/( b + ( a*d ) );
                del_h = ( ( b*d ) - 1.0F )*del_h;
                h += del_h;
                const float del_s = q*del_h;
                s += del_s;
                if ( absolute( del_s/s ) < eps ) {
                    break;
                }
            }
            h = a1*h;
            Kmu = ffmath::sqrt( ffmath::FFP_PI/( 2.0F*x ) )*ffmath::exp(-x)/s;
            k_nu1 = Kmu*( mu + x + 0.5F - h )*xi;
        }
        const float k_pmu = ( mu*xi*Kmu ) - k_nu1;
        const float i_num_u = xi/( ( f*Kmu ) - k_pmu );

        i_nu = i_num_u*i_nul1/i_nul;
        i_pnu = i_num_u*i_pnu_1/i_nul;
        for ( int i = 1 ; i <= nl ; ++i ) {
            /*cstat -CERT-FLP36-C*/
            const float k_nu_temp = ( ( mu + static_cast<float>( i ) )*xi2*k_nu1 ) + Kmu;
            /*cstat +CERT-FLP36-C*/
            Kmu = k_nu1;
            k_nu1 = k_nu_temp;
        }
        k_nu = Kmu;
        k_pnu = ( nu*xi*Kmu ) - k_nu1;
    }
}
/*============================================================================*/
static float cyl_bessel_ij_series( float nu,
                                   float x,
                                   float sgn,
                                   size_t max_iter )
{
    float y;

    if ( ffmath::isEqual( 0.0F, x ) ) {
        y = ( ffmath::isEqual( 0.0F, nu ) ) ? 1.0F : 0.0F;
    }
    else {
        const float x2 = 0.5F*x;
        float fact = nu*ffmath::log( x2 );
        float Jn = 1.0F;
        float term = 1.0F;
        const float xx4 = sgn*x2*x2;

        fact -= ffmath::lgamma( nu + 1.0F );
        fact = ffmath::exp( fact );
        for ( size_t i = 1U ; i < max_iter; ++i ) {
            /*cstat -CERT-FLP36-C*/
            const float j = static_cast<float>( i );
            /*cstat +CERT-FLP36-C*/
            term *= xx4/( j*( nu + j ) );
            Jn += term;
            if ( absolute( term/Jn ) < limits::epsilon() ) {
                break;
            }
        }
        y = fact*Jn;
    }

    return y;
}
/*============================================================================*/
float ffmath::cyl_bessel_i( float nu,
                            float x )
{
    float y;

    if ( ( nu < 0.0F ) || ( x < 0.0F ) || ffmath::isNan( nu ) || ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else if ( ( x*x ) < ( 10.0F*( nu + 1.0F ) ) ) {
        y = cyl_bessel_ij_series( nu, x, 1.0F, 200U );
    }
    else {
        float I_nu, K_nu, Ip_nu, Kp_nu;

        bessel_ik( nu, x, I_nu, K_nu, Ip_nu, Kp_nu );
        y = I_nu;
    }
    return y;
}
/*============================================================================*/
static void cyl_bessel_jn_asymp( float nu,
                                 float x,
                                 float & Jnu,
                                 float & Nnu)
{
    const float mu = 4.0F*nu*nu;
    const float x8 = 8.0F*x;
    float P = 0.0F;
    float Q = 0.0F;
    float term = 1.0F;
    const float eps = limits::epsilon();
    size_t i = 0U;

    do {
        bool epsP, epsQ;
        float k2_1;
        /*cstat -MISRAC++2008-0-1-2_b*/
        float k = static_cast<float>( i );
        /*cstat +MISRAC++2008-0-1-2_b*/
        k2_1 = 2.0F*k - 1.0F;
        term *= ( i == 0U ) ? 1.0F : -( mu - ( k2_1*k2_1 ) )/( k*x8 );
        epsP = absolute( term ) < ( eps*absolute( P ) );
        P += term;
        ++i;
        k = static_cast<float>( i );
        k2_1 = 2.0F*k - 1.0F;
        term *= ( mu - ( k2_1*k2_1 ) )/( k*x8 );
        epsQ = absolute( term ) < ( eps*absolute( Q ) );
        Q += term;
        if ( epsP && epsQ && ( k > ( 0.5F*nu ) ) ) {
            break;
        }
        ++i;
    } while ( i < 1000U );
    const float chi = x - ( nu + 0.5F )*ffmath::FFP_PI_2;
    const float c = ffmath::cos( chi );
    const float s = ffmath::sin( chi );
    const float coeff = ffmath::sqrt( 2.0F/( ffmath::FFP_PI*x ) );
    Jnu = coeff*( ( c*P ) - ( s*Q ) );
    Nnu = coeff*( ( s*P ) + ( c*Q ) );;
}
/*============================================================================*/
float ffmath::cyl_bessel_j( float nu,
                            float x )
{
    float y;

    if ( ( nu < 0.0F ) || ( x < 0.0F ) || ffmath::isNan( nu ) || ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else if ( ( x*x ) < ( 10.0F*( nu + 1.0F ) ) ) {
        y = cyl_bessel_ij_series( nu, x, -1.0F, 200U );
    }
    else if ( x > 1000.0F ) {
        float j_nu, n_nu;

        cyl_bessel_jn_asymp( nu, x, j_nu, n_nu );
        y = j_nu;
    }
    else {
        float J_nu, N_nu, Jp_nu, Np_nu;

        bessel_jn( nu, x, J_nu, N_nu, Jp_nu, Np_nu );
        y = J_nu;
    }
    return y;
}
/*============================================================================*/
float ffmath::cyl_bessel_k( float nu,
                            float x )
{
    float y;

    if ( ( nu < 0.0F ) || ( x < 0.0F ) || ffmath::isNan( nu ) || ffmath::isNan( x ) ) {
        y = ffmath::getNan();
    }
    else {
        float I_nu, K_nu, Ip_nu, Kp_nu;

        bessel_ik( nu, x, I_nu, K_nu, Ip_nu, Kp_nu );
        y = K_nu;
    }
    return y;
}
/*============================================================================*/
float ffmath::sph_legendre( size_t l,
                            size_t m,
                            float theta )
{
    float y;

    if ( ffmath::isNan( theta ) ) {
        y = ffmath::getNan();
    }
    else {
        const float x = ffmath::cos( theta );
        constexpr float pi4 = 4.0F*ffmath::FFP_PI;

        if ( m > l ) {
            y = 0.0F;
        }
        else if ( 0U == m ) {
            float P = ffmath::legendre( l, x );
            /*cstat -CERT-FLP36-C*/
            const float fact = ffmath::sqrt( static_cast<float>( 2U*l + 1U )/pi4 );
            /*cstat +CERT-FLP36-C*/
            P *= fact;
            y = P;
        }
        else if ( ffmath::isEqual( 1.0F, x ) || ffmath::isEqual( -1.0F, x ) ) {
            y = 0.0F;
        }
        else {
            /*cstat -CERT-FLP36-C*/
            const float mf = static_cast<float>( m );
            const float y_mp1m_factor = x*ffmath::sqrt( static_cast<float>( 2U*m + 3U ) );
            /*cstat +CERT-FLP36-C*/
            const float sgn = ( 1U == ( m % 2U ) ) ? -1.0F : 1.0F;
            const float ln_circ = ffmath::log( 1.0F - ( x*x ) );
            const float ln_poc_h = ffmath::lgamma( mf + 0.5F ) - ffmath::lgamma( mf );
            const float ln_pre_val = ( -0.25F*ffmath::FFP_LN_PI ) + 0.5F*( ln_poc_h + ( mf*ln_circ ) );
            const float sr = ffmath::sqrt( ( 2.0F + ( 1.0F/mf ) )/pi4);
            float y_mm = sgn*sr*ffmath::exp( ln_pre_val );
            float y_mp1m = y_mp1m_factor*y_mm;

            if ( l == m ) {
                y = y_mm;
            }
            else if ( l == ( m + 1U ) ) {
                y = y_mp1m;
            }
            else {
                float y_lm = 0.0F;

                for ( size_t ll = ( m + 2U ) ; ll <= l; ++ll ) {
                    /*cstat -CERT-FLP36-C*/
                    const float ll_m_m = static_cast<float>( ll - m );
                    const float ll_p_m = static_cast<float>( ll + m );
                    const float ll2_p_1 = static_cast<float>( ( 2U*ll ) + 1U );
                    const float ll2_m_1 = static_cast<float>( ( 2U*ll ) - 1U );
                    const float ll_pm_m1 = static_cast<float>( ll + m - 1U );
                    const float ll_mm_m1 = static_cast<float>( ll - m - 1U );
                    /*cstat +CERT-FLP36-C*/
                    const float rat1 = ll_m_m/ll_p_m;
                    const float fact1 = ffmath::sqrt( rat1*ll2_p_1*ll2_m_1 );
                    /*cstat -CERT-FLP36-C*/
                    const float fact2 = ffmath::sqrt( rat1*( ll_mm_m1/ll_pm_m1 )*ll2_p_1/static_cast<float>( 2U*ll - 3U ) );
                    /*cstat -CERT-FLP36-C*/
                    y_lm = ( x*y_mp1m*fact1 - ll_pm_m1*y_mm*fact2 )/ll_m_m;
                    y_mm = y_mp1m;
                    y_mp1m = y_lm;
                }
                y = y_lm;
            }
        }
    }
    return y;
}
/*============================================================================*/