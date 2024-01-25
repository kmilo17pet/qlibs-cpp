#include <include/ffmath.hpp>

using namespace qlibs;

template<typename T1, typename T2>
static inline void cast_reinterpret( T1 &f, const T2 &u )
{
    static_assert( sizeof(T1)==sizeof(T2), "Types must match sizes" );
    (void)memcpy( &f, &u, sizeof(T2) );
}

static float getAbnormal( const int i );
static float compute_cbrt( float x , bool r );
static float lgamma_positive( float x );

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
    return ( ffmath::absf( a - b ) <= ffmath::absf( tol ) );
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
    return x - trunc( x );
}
/*============================================================================*/
float ffmath::rem( float x, float y )
{
    return ( classification::FFP_ZERO  == classify( x ) ) ? getNan() : ( x - ( y*trunc( x/y ) ) );
}
/*============================================================================*/
float ffmath::mod( float x, float y )
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
    return ffmath::sin( x + ffmath::FFP_PI_2 );
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
float ffmath::atan2( float y, float x )
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
float ffmath::rexp( float x, int32_t *pw2 )
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
float ffmath::tgamma( float x )
{
    float result;

    const auto fClass = ffmath::classify( x );
    if ( classification::FFP_NAN == fClass ) {
        result = getNan();
    }
    else if ( classification::FFP_ZERO == fClass ) {
        result = getInf(); /* a huge value */
    }
    else if ( classification::FFP_INFINITE == fClass ) {
        if ( x > 0.0F ) {
            result = getInf(); /* a huge value */
        }
        else {
            result = getNan();
        }
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
                result = getInf();
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
        result = ffmath::getInf(); /* a huge value */
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
        result = getNan();
    }
    else if ( classification::FFP_ZERO == fClass ) {
        result = getInf();
    }
    else if ( classification::FFP_INFINITE == fClass ) {
        result = getInf();
    }
    else {
        if ( x < 0.0F ) {
            if ( x <= -4503599627370496.0F ) { /* x < 2^52 */
                result = getInf();
            }
            else {
                float y, y1, isItAnInt;

                y = -x;
                y1 = ffmath::trunc( y );
                isItAnInt = y - y1;
                if ( ffmath::isEqual( 0.0F, isItAnInt ) ) {
                    result = getInf();
                }
                else {
                    float a;

                    a = sin( FFP_PI*isItAnInt );
                    result = ffmath::log( FFP_PI/ffmath::absf( a*x ) ) - lgamma_positive( -x );
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