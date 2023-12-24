#include <include/fp16.hpp>


using namespace qlibs;

fp16_t fp16::Min = -2147483647;
fp16_t fp16::Max = 2147483647;
bool fp16::rounding = true;
bool fp16::saturation = false;

const fp16_t fp16::exp_max = 681391;
const fp16_t fp16::f2 = 131072;
const fp16_t fp16::f3 = 196608;
const fp16_t fp16::f16 = 1048576;
const fp16_t fp16::f100 = 6553600;
const fp16_t fp16::f6_5 = 425984;
const float fp16::one_fp16_f = 0.0000152587890625f;
const double fp16::one_fp16_d = 0.0000152587890625;
const uint32_t fp16::overflow_mask = 0x80000000U;
const uint32_t fp16::fraction_mask = 0x0000FFFFU;
const uint32_t fp16::integer_mask = 0xFFFF0000U;

const fp16_t fp16::f_e = 178145;
const fp16_t fp16::f_log2e = 94548;
const fp16_t fp16::f_log10e = 28462;
const fp16_t fp16::f_ln2 = 45426;
const fp16_t fp16::f_ln10 = 150902;
const fp16_t fp16::f_pi = 205887;
const fp16_t fp16::f_pi_2 = 102944;
const fp16_t fp16::f_2pi = 411775;
const fp16_t fp16::f_pi_4 = 51471;
const fp16_t fp16::f_1_pi = 20861;
const fp16_t fp16::f_2_pi = 41722;
const fp16_t fp16::f_2_sqrtpi = 73949;
const fp16_t fp16::f_sqrt2 =92682 ;
const fp16_t fp16::f_sqrt1_2 = 46341;
const fp16_t fp16::epsilon = 1;
const fp16_t fp16::MaxValue = 2147483647;
const fp16_t fp16::overflow = -2147483647 - 1;
const fp16_t fp16::one = 65536;
const fp16_t fp16::one_half = 32768;
const fp16_t fp16::f_180_pi = 3754936;
const fp16_t fp16::f_pi_180 = 1144;
const fp16_t fp16::f_180 = 11796480;
const fp16_t fp16::f_360 = 23592960;


fp16 fp16Constant::operator+(const fp16Constant& other) const {
    fp16 ret;
    ret.value = fp16::add( value, other.value );
    return ret;
}

fp16 fp16Constant::operator-(const fp16Constant& other) const {
    fp16 ret;
    ret.value = fp16::sub( value, other.value );
    return ret;
}

fp16 fp16Constant::operator*(const fp16Constant& other) const {
    fp16 ret;
    ret.value = fp16::mul( value, other.value );
    return ret;
}

fp16 fp16Constant::operator/(const fp16Constant& other) const {
    fp16 ret;
    ret.value = fp16::div( value, other.value );
    return ret;
}

/*cstat -MISRAC++2008-5-0-21 -MISRAC++2008-5-0-9 -ATH-shift-neg -CERT-INT34-C_c -MISRAC++2008-5-0-10*/

/*============================================================================*/
int fp16::toInt( const fp16 &x )
{
    int retValue;

    if ( rounding ) {
        if ( x.value >= 0 ) {
            retValue = ( x.value + ( one >> 1 ) ) / one;
        }
        else {
            retValue = ( x.value - ( one >> 1 ) ) / one;
        }
    }
    else {
        retValue = static_cast<int>( static_cast<uint32_t>( x.value ) >> 16 );
    }

    return retValue;
}
/*============================================================================*/
float fp16::toFloat( const fp16 &x )
{
    /*cstat -CERT-FLP36-C*/
    return static_cast<float>( x.value )*one_fp16_f;
    /*cstat +CERT-FLP36-C*/
}
/*============================================================================*/
double fp16::toDouble( const fp16 &x )
{
    return static_cast<double>( x.value )*one_fp16_d;
}
/*============================================================================*/
fp16_t fp16::fromInt( const int x )
{
    return static_cast<fp16_t>( static_cast<uint32_t>( x ) << 16U );
}
/*============================================================================*/
fp16_t fp16::fromFloat( const float x )
{
    float f_value;
    /*cstat -CERT-FLP36-C -CERT-FLP34-C*/
    f_value = x * static_cast<float>( one );
    if ( rounding ) {
        f_value += ( f_value >= 0.0f ) ? 0.5f : -0.5f;
    }
    return static_cast<fp16_t>( f_value );
    /*cstat +CERT-FLP36-C +CERT-FLP34-C*/
}
/*============================================================================*/
fp16_t fp16::fromDouble( const double x )
{
    double d_value;
    /*cstat -CERT-FLP36-C -CERT-FLP34-C*/
    d_value = x * static_cast<double>( one );
    if ( rounding ) {
        d_value += ( d_value >= 0.0 ) ? 0.5 : -0.5;
    }
    return static_cast<fp16_t>( d_value );
    /*cstat +CERT-FLP36-C +CERT-FLP34-C*/
}
/*============================================================================*/
fp16_t fp16::saturate( const fp16_t nsInput, const fp16_t x, const fp16_t y )
{
    fp16_t retValue = nsInput;

    if ( saturation ) {
        if ( overflow == nsInput ) {
            retValue = ( ( x >= 0 ) == ( y >= 0 ) ) ? Max : Min;
        }
    }

    return retValue;
}
/*============================================================================*/
fp16_t fp16::add( const fp16_t X, const fp16_t Y )
{
    const uint32_t x = static_cast<uint32_t>( X );
    const uint32_t y = static_cast<uint32_t>( Y );
    uint32_t retValue;

    retValue =  x + y;
    if ( ( 0U == ( ( x ^ y ) & overflow_mask ) ) && ( 0U != ( ( x ^ retValue ) & overflow_mask ) ) ) {
        retValue = static_cast<uint32_t>( overflow );
    }

    return saturate( static_cast<fp16_t>( retValue ), X, X );
}
/*============================================================================*/
fp16_t fp16::sub( const fp16_t X, const fp16_t Y )
{
    const uint32_t x = static_cast<uint32_t>( X );
    const uint32_t y = static_cast<uint32_t>( Y );
    uint32_t retValue;

    retValue =  x - y;
    if ( ( 0U != ( ( x ^ y ) & overflow_mask ) ) && ( 0U != ( ( x ^ retValue ) & overflow_mask ) ) ) {
        retValue = static_cast<uint32_t>( overflow );
    }

    return saturate( static_cast<fp16_t>( retValue ), X, X );
}
/*============================================================================*/
fp16_t fp16::mul( const fp16_t x, const fp16_t y )
{
    fp16_t retValue = overflow;
    fp16_t a, c, ac, ad_cb, mulH;
    uint32_t b, d, bd, tmp, mulL;

    a = ( x >> 16 );
    c = ( y >> 16 );
    b = static_cast<uint32_t>( x & 0xFFFF );
    d = static_cast<uint32_t>( y & 0xFFFF );
    ac = a*c;
    ad_cb = static_cast<fp16_t>( ( static_cast<uint32_t>( a )*d ) + ( static_cast<uint32_t>( c )*b ) );
    bd = b*d;
    mulH = ac + ( ad_cb >> 16 );
    tmp = static_cast<uint32_t>( ad_cb ) << 16;
    mulL = bd + tmp;
    if ( mulL < bd ) {
        ++mulH;
    }
    /*cstat -MISRAC++2008-5-0-3*/
    a = ( mulH < 0 ) ? -1 : 0;
    /*cstat +MISRAC++2008-5-0-3*/
    if ( a == ( mulH >> 15 ) ) {
        if ( rounding ) {
            uint32_t tmp2;

            tmp2 = mulL;
            mulL -= static_cast<uint32_t>( one_half );
            mulL -= static_cast<uint32_t>( mulH ) >> 31;
            if ( mulL > tmp2 ) {
                --mulH;
            }
            retValue = static_cast<fp16_t>( mulH << 16 ) | static_cast<fp16_t>( mulL >> 16 );
            retValue += 1;
        }
        else {
            retValue = static_cast<fp16_t>( mulH << 16 ) | static_cast<fp16_t>( mulL >> 16 );
        }
    }

    return saturate( retValue, x, y );
}
/*============================================================================*/
fp16_t fp16::div( const fp16_t x, const fp16_t y )
{
    fp16_t retValue = Min;

    if ( 0 != y ) {
        uint32_t xRem, xDiv, bit = 0x10000U;

        xRem = static_cast<uint32_t>( ( x >= 0 ) ? x : -x );
        xDiv = static_cast<uint32_t>( ( y >= 0 ) ? y : -y );

        while ( xDiv < xRem ) {
            xDiv <<= 1;
            bit <<= 1;
        }
        retValue = overflow;
        /*cstat -MISRAC++2008-0-1-2_a*/
        if ( 0U != bit ) {
        /*cstat +MISRAC++2008-0-1-2_a*/
            uint32_t quotient = 0U;

            if ( 0U != ( xDiv & 0x80000000U ) ) {
                if ( xRem >= xDiv ) {
                    quotient |= bit;
                    xRem -= xDiv;
                }
                xDiv >>= 1;
                bit >>= 1;
            }

            while ( ( 0U != bit ) && ( 0U != xRem ) ) {
                if ( xRem >= xDiv ) {
                    quotient |= bit;
                    xRem -= xDiv;
                }
                xRem <<= 1;
                bit >>= 1;
            }
            if ( rounding ) {
                if ( xRem >= xDiv ) {
                    ++quotient;
                }
            }

            retValue = static_cast<fp16_t>( quotient );
            if ( 0U != ( static_cast<uint32_t>( x ^ y ) & overflow_mask ) ) {
                if ( quotient == static_cast<uint32_t>( Min ) ) {
                    retValue = overflow;
                }
                else {
                    retValue = -retValue;
                }
            }
        }
    }

    return saturate( retValue, x, y );
}
/*============================================================================*/
fp16_t fp16::abs( fp16_t x )
{
    fp16_t retValue;

    if ( x == Min ) {
        retValue = overflow;
    }
    else {
        retValue = ( x >= 0 ) ? x : -x;
    }

    return retValue;
}
/*============================================================================*/
fp16 fp16::abs( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::abs( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::abs( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::abs( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::sqrt( fp16_t x )
{
    fp16_t retValue = overflow;

    if ( x > 0 ) {
        uint32_t bit;
        uint8_t n;

        retValue = 0;
        /*cstat -MISRAC++2008-5-0-3*/
        bit = ( 0 != ( x & static_cast<fp16_t>( 4293918720 ) ) ) ? ( 1U << 30U ) : ( 1U << 18U );
        /*cstat +MISRAC++2008-5-0-3*/
        while ( bit > static_cast<uint32_t>( x ) ) {
            bit >>= 2U;
        }

        for ( n = 0U ; n < 2U ; ++n ) {
            while ( 0U != bit ) {
                if ( x >= static_cast<fp16_t>( static_cast<uint32_t>( retValue ) + bit ) ) {
                    x -= static_cast<fp16_t>( static_cast<uint32_t>( retValue ) + bit );
                    retValue = static_cast<fp16_t>( ( static_cast<uint32_t>( retValue ) >> 1U ) + bit );
                }
                else {
                    retValue = ( retValue >> 1 );
                }
                bit >>= 2U;
            }

            if ( 0U == n ) {
                if ( x > 65535 ) {
                    x -= retValue;
                    x = ( x << 16 ) - one_half;
                    retValue = ( retValue << 16 ) + one_half;
                }
                else {
                    x <<= 16;
                    retValue <<= 16;
                }
                bit = 1U << 14U;
            }
        }
    }
    if ( ( 1U == rounding ) && ( x > retValue ) ) {
        ++retValue;
    }

    return static_cast<fp16_t>( retValue );
}
/*============================================================================*/
fp16 fp16::sqrt( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::sqrt( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::sqrt( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::sqrt( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::exp( fp16_t x )
{
    fp16_t retValue, term;
    bool isNegative;
    int i;

    if ( 0 == x ) {
        retValue = one;
    }
    else if ( x == one ) {
        retValue = f_e;
    }
    else if ( x >= exp_max ) {
        retValue = Max;
    }
    else if ( x <= -exp_max ) {
        retValue = 0;
    }
    else {
        isNegative = ( x < 0 );
        if ( isNegative ) {
            x = -x;
        }

        retValue = x + one;
        term = x;

        for ( i = 2 ; i < 30 ; ++i ) {
            term = mul( term, div( x, fromInt( i ) ) );
            retValue += term;

            if ( ( term < 500 ) && ( ( i > 15 ) || ( term < 20 ) ) ) {
                break;
            }
        }

        if ( isNegative ) {
            retValue = div( one, retValue );
        }
    }

    return retValue;
}
/*============================================================================*/
fp16 fp16::exp( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::exp( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::exp( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::exp( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::log( fp16_t x )
{
    fp16_t retValue = overflow;
    static const fp16_t e4 = 3578144; /*e^4*/

    if ( x > 0 ) {
        fp16_t guess = f2, delta;
        int scaling = 0, count = 0;

        while ( x > f100 ) {
            x = div( x, e4 );
            scaling += 4;
        }

        while ( x < one ) {
            x = mul( x, e4 );
            scaling -= 4;
        }

        do {
            const fp16_t e = exp( guess );

            delta = div( x - e , e );

            if ( delta > f3 ) {
                delta = f3;
            }
            guess += delta;
        } while ( ( count++ < 10 ) && ( ( delta > 1 ) || ( delta < -1 ) ) );

        retValue = guess + fromInt( scaling );
    }

    return retValue;
}
/*============================================================================*/
fp16 fp16::log( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::log( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::log( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::log( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::rs( fp16_t x )
{
    fp16_t retValue;

    if ( rounding ) {
        retValue = ( x >> 1U ) + ( x & 1 );
    }
    else {
        retValue = x >> 1;
    }

    return retValue;
}
/*============================================================================*/
fp16_t fp16::log2i( fp16_t x )
{
    fp16_t retValue = 0;

    while ( x >= f2 ) {
        ++retValue;
        x = rs( x );
    }

    if ( 0 == x ) {
        retValue = retValue << 16;
    }
    else {
        int i;
        for ( i = 16 ; i > 0 ; --i ) {
            x = mul( x, x );
            retValue <<= 1;
            if ( x >= f2 ) {
                retValue |= 1;
                x = rs( x );
            }
        }
        if ( rounding ) {
            x = mul( x, x );
            if ( x >= f2 ) {
                ++retValue;
            }
        }
    }

    return retValue;
}
/*============================================================================*/
fp16_t fp16::log2( fp16_t x )
{
    fp16_t retValue = overflow;

    if ( x > 0 ) {
        if ( x < one ) {
            if ( 1 == x ) {
                retValue = -f16;
            }
            else {
                fp16_t inv;
                inv = div( one, x );
                retValue = -log2i( inv );
            }
        }
        else {
            retValue = log2i( x );
        }
    }
    if ( saturation ) {
        if ( overflow == retValue ) {
            retValue = Min;
        }
    }

    return retValue;
}
/*============================================================================*/
fp16 fp16::log2( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::log2( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::log2( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::log2( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::wrapToPi( fp16_t x )
{
    if ( ( x < -f_pi ) || ( x > f_pi ) ) {
        while ( x > f_pi ) {
            x -= f_2pi;
        }
        while ( x <= -f_pi ) {
            x += f_2pi;
        }
    }

    return x;
}
/*============================================================================*/
fp16 fp16::wrapToPi( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::wrapToPi( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::wrapToPi( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::wrapToPi( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::wrapTo180( fp16_t x )
{
    if ( ( x < -f_180 ) || ( x > f_180 ) ) {
        while ( x > f_180 ) {
            x -= f_360;
        }
        while ( x <= -f_pi ) {
            x += f_360;
        }
    }

    return x;
}
/*============================================================================*/
fp16 fp16::wrapTo180( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::wrapTo180( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::wrapTo180( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::wrapTo180( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::sin( fp16_t x )
{
    fp16_t retValue, x2;

    x = wrapToPi( x );
    x2 = mul( x ,x );
    retValue = x;
    x = mul( x, x2 );
    retValue -= ( x / 6 ); /*x^3/3!*/
    x = mul( x, x2 );
    retValue += ( x / 120 ); /*x^5/5!*/
    x = mul( x, x2 );
    retValue -= ( x / 5040 ); /*x^7/7!*/
    x = mul( x, x2 );
    retValue += ( x / 362880 ); /*x^9/9!*/
    x = mul( x, x2);
    retValue -= ( x / 39916800 ); /*x^11/11!*/

    return retValue;
}
/*============================================================================*/
fp16 fp16::sin( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::sin( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::sin( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::sin( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::cos( fp16_t x )
{
    return sin( x + f_pi_2 );
}
/*============================================================================*/
fp16 fp16::cos( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::cos( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::cos( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::cos( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::tan( fp16_t x )
{
    fp16_t a ,b;

    a = sin( x );
    b = cos( x );

    return div( a, b );
}
/*============================================================================*/
fp16 fp16::tan( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::tan( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::tan( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::tan( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::atan2( fp16_t y, fp16_t x )
{
    fp16_t absY, mask, angle, r, r_3;
    const fp16_t QFP16_0_981689 = 0x0000FB50;
    const fp16_t QFP16_0_196289 = 0x00003240;
    static const fp16_t f_3pi_div_4 = 154415; /*3*pi/4*/

    mask = ( y >> ( sizeof(fp16_t)*7U ) );
    absY = ( y + mask ) ^ mask;
    if ( x >= 0 ) {
        r = div( ( x - absY ), ( x + absY ) );
        angle = f_pi_4;
    }
    else {
        r = div( ( x + absY ), ( absY - x ) );
        angle = f_3pi_div_4;
    }
    r_3 = mul( mul( r, r ), r );
    angle += mul( QFP16_0_196289, r_3 ) - mul( QFP16_0_981689, r );
    /*cstat -ATH-neg-check-nonneg*/
    if ( y < 0 ) {
        angle = -angle;
    }
    /*cstat +ATH-neg-check-nonneg*/
    return angle;
}
/*============================================================================*/
/*cstat -CERT-EXP30-C_b*/
fp16 fp16::atan2( const fp16 &y, const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::atan2( y.raw(), x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::atan2( const fp16 &y, const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::atan2( y.raw(), x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::atan2( const fp16Constant &y, const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::atan2( y.raw(), x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::atan2( const fp16Constant &y, const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::atan2( y.raw(), x.raw() );
    return ret;
}
/*cstat +CERT-EXP30-C_b*/
/*============================================================================*/
fp16_t fp16::atan( fp16_t x )
{
    return atan2( x, one );
}
/*============================================================================*/
fp16 fp16::atan( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::atan( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::atan( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::atan( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::asin( fp16_t x )
{
    fp16_t retValue = 0;

    if ( ( x <= one ) && ( x >= -one ) ) {
        retValue = one - mul( x, x );
        retValue = div( x, sqrt( retValue ) );
        retValue = atan( retValue );
    }

    return retValue;
}
/*============================================================================*/
fp16 fp16::asin( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::asin( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::asin( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::asin( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::acos( fp16_t x )
{
    return ( f_pi_2 - asin( x ) );
}
/*============================================================================*/
fp16 fp16::acos( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::acos( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::acos( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::acos( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::cosh( fp16_t x )
{
    fp16_t retValue = overflow;
    fp16_t epx, enx;

    if ( 0 == x ) {
        retValue = one;
    }
    else if ( ( x >= exp_max ) || ( x <= -exp_max ) ) {
        retValue = Max;
    }
    else {
        epx = exp( x );
        enx = exp( -x );
        if ( ( overflow != epx ) && ( overflow != enx ) ) {
            retValue = epx + enx;
            retValue = ( retValue >> 1 );
        }
    }

    return retValue;
}
/*============================================================================*/
fp16 fp16::cosh( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::cosh( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::cosh( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::cosh( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::sinh( fp16_t x )
{
    fp16_t retValue = overflow;
    fp16_t epx, enx;

    if ( 0 == x ) {
        retValue = one;
    }
    else if ( x >= exp_max ) {
        retValue = Max;
    }
    else if ( x <= -exp_max ) {
        retValue = -Max;
    }
    else {
        epx = exp( x );
        enx = exp( -x );
        if ( ( overflow != epx ) && ( overflow != enx ) ) {
            retValue = epx - enx;
            retValue = ( retValue >> 1 );
        }
    }

    return retValue;
}
/*============================================================================*/
fp16 fp16::sinh( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::sinh( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::sinh( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::sinh( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::tanh( fp16_t x )
{
    fp16_t retValue, epx, enx;

    if ( 0 == x ) {
        retValue = 0;
    }
    else if ( x >  f6_5 ) { /* tanh for any x>6.5 ~= 1*/
        retValue = one;
    }
    else if ( x < -f6_5 ) { /* tanh for any x<6.5 ~= -1*/
        retValue = -one;
    }
    else {
        retValue = abs( x );
        epx = exp( retValue );
        enx = exp( -retValue );
        retValue = div( epx - enx, epx + enx );
        retValue = ( x > 0 ) ? retValue : -retValue;
    }

    return retValue;
}
/*============================================================================*/
fp16 fp16::tanh( const fp16 &x )
{
    fp16 ret;

    ret.value = fp16::tanh( x.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::tanh( const fp16Constant &x )
{
    fp16 ret;

    ret.value = fp16::tanh( x.raw() );
    return ret;
}
/*============================================================================*/
fp16_t fp16::powi( fp16_t x, fp16_t y )
{
    fp16_t retValue;
    fp16_t n;
    int32_t i;

    retValue = one;
    n = y >> 16;
    if ( 0 == n ) {
        retValue = one;
    }
    else if ( one == n ) {
        retValue = x;
    }
    else {
        for ( i = 0 ; i < n ; ++i ) {
            retValue = mul( x, retValue );
            if ( overflow == retValue ) {
                break;
            }
        }
    }

    return retValue;
}
/*============================================================================*/
fp16_t fp16::pow( fp16_t x, fp16_t y )
{
    fp16_t retValue = overflow;

    if ( ( 0U == ( static_cast<uint32_t>( y ) & fraction_mask ) ) && ( y > 0 ) ) {
        retValue = powi( x, y );
    }
    else {
        fp16_t tmp;
        tmp = mul( y, log( abs( x ) ) );
        if ( overflow != tmp ) {
            retValue = exp( tmp );
            if ( x < 0 ) {
                retValue = -retValue;
            }
        }
    }

    return retValue;
}
/*============================================================================*/
/*cstat -CERT-EXP30-C_b*/
fp16 fp16::pow( const fp16 &x, const fp16 &y )
{
    fp16 ret;

    ret.value = fp16::pow( x.raw(), y.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::pow( const fp16 &x, const fp16Constant &y )
{
    fp16 ret;

    ret.value = fp16::pow( x.raw(), y.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::pow( const fp16Constant &x, const fp16 &y )
{
    fp16 ret;

    ret.value = fp16::pow( x.raw(), y.raw() );
    return ret;
}
/*============================================================================*/
fp16 fp16::pow( const fp16Constant &x, const fp16Constant &y )
{
    fp16 ret;

    ret.value = fp16::pow( x.raw(), y.raw() );
    return ret;
}
/*cstat +CERT-EXP30-C_b*/
/*============================================================================*/
char* fp16::itoa( char *buf, uint32_t scale, uint32_t val, uint8_t skip )
{
    while ( 0U != scale ) {
        const uint32_t digit = ( val / scale );
        if ( ( 0U == skip ) || ( 0U != digit ) || ( 1U == scale ) ) {
            skip = 0U;
            /*cstat -MISRAC++2008-5-0-3*/
            *buf++ = static_cast<char>( '0' ) + static_cast<char>( digit );
            /*cstat +MISRAC++2008-5-0-3*/
            val %= scale;
        }
        scale /= 10U;
    }

    return buf;
}
/*============================================================================*/
char* fp16::toASCII( const fp16_t num, char *str, int decimals )
{
    char * const retValue = str;

    if ( overflow == num ) {
        str[ 0 ] = 'o';
        str[ 1 ] = 'v';
        str[ 2 ] = 'e';
        str[ 3 ] = 'r';
        str[ 4 ] = 'f';
        str[ 5 ] = 'l';
        str[ 6 ] = 'o';
        str[ 7 ] = 'w';
        str[ 8 ] = '\0';
    }
    else {
        const uint32_t iScales[ 6 ] = { 1U, 10U, 100U, 1000U, 10000U, 100000U };
        uint32_t uValue, fPart, scale;
        fp16_t iPart;

        uValue = static_cast<uint32_t>( ( num >= 0 ) ? num : -num );
        if ( num < 0 ) {
            *str++ = '-';
        }

        iPart = static_cast<fp16_t>( uValue >> 16 );
        fPart = uValue & fraction_mask;
        if ( decimals > 5 ) {
            decimals = 5;
        }
        if ( decimals < 0 ) {
            decimals = 0;
        }
        scale = iScales[ decimals ];
        fPart = static_cast<uint32_t>( mul( static_cast<fp16_t>( fPart ), static_cast<fp16_t>( scale ) ) );

        if ( fPart >= scale ) {
            iPart++;
            fPart -= scale;
        }
        str = itoa( str, 10000U, static_cast<uint32_t>( iPart ), 1U );

        if ( 1U != scale ) {
            *str++ = '.';
            str = itoa( str, scale/10U, fPart, 0U );
        }
        *str = '\0';
    }

    return retValue;
}
/*============================================================================*/
char* fp16::toASCII( const fp16 &x, char *str, int decimals )
{
    return toASCII( x.raw(), str, decimals );
}
/*============================================================================*/
char* fp16::toASCII( const fp16Constant &x, char *str, int decimals )
{
    return toASCII( x.raw(), str, decimals );
}
/*============================================================================*/

/*cstat +MISRAC++2008-5-0-21 +MISRAC++2008-5-0-9 +ATH-shift-neg +CERT-INT34-C_c +MISRAC++2008-5-0-10*/