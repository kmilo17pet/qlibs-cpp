#ifndef QLIBS_FP16
#define QLIBS_FP16

#include "include/types.hpp" 
#include <iostream>

namespace qlibs {

    using fp16Raw_t = int32_t;

    class fp16 {
        struct fp16Hidden {
            fp16Raw_t x;
        };
        private:
            fp16Raw_t value{ overflow };
            static fp16Raw_t Min; // skipcq: CXX-W2009
            static fp16Raw_t Max; // skipcq: CXX-W2009
            static bool rounding; // skipcq: CXX-W2009
            static bool saturation; // skipcq: CXX-W2009
            
            static const fp16Raw_t exp_max;
            static const fp16Raw_t f2;
            static const fp16Raw_t f3;
            static const fp16Raw_t f16;
            static const fp16Raw_t f100;
            static const fp16Raw_t f6_5;
            static const float one_fp16_f;
            static const double one_fp16_d;
            static const double one_dbl;
            static const uint32_t overflow_mask;
            static const uint32_t fraction_mask;
            static const uint32_t integer_mask;

            static uint32_t overflowCheck( uint32_t res, const uint32_t x, const uint32_t y );
            static fp16Raw_t saturate( const fp16Raw_t nsInput, const fp16Raw_t x, const fp16Raw_t y );
            static fp16Raw_t fromInt( const int x );
            static fp16Raw_t fromFloat( const float x );
            static fp16Raw_t fromDouble( const double x );

            static fp16Raw_t add( const fp16Raw_t X, const fp16Raw_t Y );
            static fp16Raw_t sub( const fp16Raw_t X, const fp16Raw_t Y );
            static fp16Raw_t mul( const fp16Raw_t x, const fp16Raw_t y );
            static fp16Raw_t div( const fp16Raw_t x, const fp16Raw_t y );

            static fp16Raw_t abs( fp16Raw_t x );
            static fp16Raw_t sqrt( fp16Raw_t x );
            static fp16Raw_t exp( fp16Raw_t x );
            static fp16Raw_t log( fp16Raw_t x );
            static fp16Raw_t log2( fp16Raw_t x );
            static fp16Raw_t wrapToPi( fp16Raw_t x );
            static fp16Raw_t wrapTo180( fp16Raw_t x );
            static fp16Raw_t sin( fp16Raw_t x );
            static fp16Raw_t cos( fp16Raw_t x );
            static fp16Raw_t tan( fp16Raw_t x );
            static fp16Raw_t atan2( fp16Raw_t y, fp16Raw_t x );
            static fp16Raw_t atan( fp16Raw_t x );
            static fp16Raw_t asin( fp16Raw_t x );
            static fp16Raw_t acos( fp16Raw_t x );
            static fp16Raw_t cosh( fp16Raw_t x );
            static fp16Raw_t sinh( fp16Raw_t x );
            static fp16Raw_t tanh( fp16Raw_t x );
            static fp16Raw_t powi( fp16Raw_t x, fp16Raw_t y );
            static fp16Raw_t pow( fp16Raw_t x, fp16Raw_t y );

            static char* itoa( char *buf, uint32_t scale, uint32_t val, uint8_t skip );
            static char* toASCII( const fp16Raw_t num, char *str, int decimals );

            static fp16Raw_t rs( fp16Raw_t x );
            static fp16Raw_t log2i( fp16Raw_t x );

            constexpr fp16( fp16Hidden val ) : value( val.x ) {}
            friend constexpr fp16 operator"" _fp( long double val );
            friend constexpr fp16 operator"" _fp( unsigned long long val );

        public:
            constexpr fp16() : value( 0 ) {}
            fp16( const fp16& other) : value( other.value ) {}

            static const fp16Raw_t f_e;
            static const fp16Raw_t f_log2e;
            static const fp16Raw_t f_log10e;
            static const fp16Raw_t f_ln2;
            static const fp16Raw_t f_ln10;
            static const fp16Raw_t f_pi;
            static const fp16Raw_t f_pi_2;
            static const fp16Raw_t f_2pi;
            static const fp16Raw_t f_pi_4;
            static const fp16Raw_t f_1_pi;
            static const fp16Raw_t f_2_pi;
            static const fp16Raw_t f_2_sqrtpi;
            static const fp16Raw_t f_sqrt2;
            static const fp16Raw_t f_sqrt1_2;
            static const fp16Raw_t epsilon;
            static const fp16Raw_t MaxValue;
            static const fp16Raw_t overflow;
            static const fp16Raw_t one;
            static const fp16Raw_t one_half;
            static const fp16Raw_t f_180_pi;
            static const fp16Raw_t f_pi_180;
            static const fp16Raw_t f_180;
            static const fp16Raw_t f_360;

            inline fp16Raw_t raw( void ) const
            {
                return value;
            }
            inline fp16 operator+( const fp16 &other )
            {
                return fp16( { add( value, other.value ) } );
            }
            inline fp16& operator+=( const fp16 &other ) {
                value = add( value, other.value );
                return *this;
            }
            inline fp16 operator-() const {
                return fp16( { -value } );
            }
            inline fp16& operator-=( const fp16 &other ) {
                value = sub( value, other.value );
                return *this;
            }
            inline fp16 operator-( const fp16 &other )
            {
                return fp16( { sub( value, other.value ) } );
            }
            inline fp16 operator*( const fp16 &other )
            {
                return fp16( { mul( value, other.value ) } );
            }
            inline fp16& operator*=( const fp16  &other ) {
                value = mul( value, other.value );
                return *this;
            }
            inline fp16 operator/( const fp16 &other )
            {
                return fp16( { div( value, other.value ) } );
            }
            inline fp16& operator/=( const fp16  &other ) {
                value = div( value, other.value );
                return *this;
            }
            inline fp16& operator++()
            {
                value = add( value, one );
                return *this;
            }
            inline fp16 operator++(int)
            {
                fp16 temp = *this;
                value = add( value, one );;
                return temp;
            }
            inline fp16& operator--()
            {
                value = sub( value, one );
                return *this;
            }
            inline fp16 operator--(int)
            {
                fp16 temp = *this;
                value = sub( value, one );;
                return temp;
            }
            inline bool operator>(const fp16 &other) const
            {
                return value > other.value;
            }
            inline bool operator>=(const fp16 &other) const
            {
                return value >= other.value;
            }
            inline bool operator<(const fp16 &other) const
            {
                return value < other.value;
            }
            inline bool operator<=(const fp16 &other) const
            {
                return value <= other.value;
            }
            inline bool operator==(const fp16 &other) const
            {
                return value == other.value;
            }
            inline bool operator!=(const fp16 &other) const
            {
                return value != other.value;
            }
            inline fp16& operator=( int x )
            {
                value = fromInt( x );
                return *this;
            }
            inline fp16& operator=( float x )
            {
                value = fromFloat( x );
                return *this;
            }
            inline fp16& operator=( double x )
            {
                value = fromDouble( x );
                return *this;
            }
            inline fp16& operator=( const fp16 &other )
            {
                value = other.value;
                return *this;
            }
            static int toInt( const fp16 &x );
            static float toFloat( const fp16 &x );
            static double toDouble( const fp16 &x );

            static inline fp16 from( const int x )
            {
                return fp16( { fromInt( x ) } );
            }
            static inline fp16 from( const float x )
            {
                return fp16( { fromFloat( x ) } );
            }
            static inline fp16 from( const double x )
            {
                return fp16( { fromDouble( x ) } );
            }
            static inline fp16 abs( const fp16 &x )
            {
                return fp16( { abs( x.raw() ) } );
            }
            static inline fp16 sqrt( const fp16 &x )
            {
                return fp16( { sqrt( x.raw() ) } );
            }
            static inline fp16 exp( const fp16 &x )
            {
                return fp16( { exp( x.raw() ) } );
            }
            static inline fp16 log( const fp16 &x )
            {
                return fp16( { log( x.raw() ) } );
            }
            static inline fp16 log2( const fp16 &x )
            {
                return fp16( { log2( x.raw() ) } );
            }
            static inline fp16 wrapToPi( const fp16 &x )
            {
                return fp16( { wrapToPi( x.raw() ) } );
            }
            static inline fp16 wrapTo180( const fp16 &x )
            {
                return fp16( { wrapTo180( x.raw() ) } );
            }
            static inline fp16 sin( const fp16 &x )
            {
                return fp16( { sin( x.raw() ) } );
            }
            static inline fp16 cos( const fp16 &x )
            {
                return fp16( { cos( x.raw() ) } );
            }
            static inline fp16 tan( const fp16 &x )
            {
                return fp16( { tan( x.raw() ) } );
            }
            static inline fp16 atan2( const fp16 &y, const fp16 &x )
            {
                return fp16( { atan2( y.raw(), x.raw() ) } );
            }
            static inline fp16 atan( const fp16 &x )
            {
                return fp16( { atan( x.raw() ) } );
            }
            static inline fp16 asin( const fp16 &x )
            {
                return fp16( { asin( x.raw() ) } );
            }
            static inline fp16 acos( const fp16 &x )
            {
                return fp16( { acos( x.raw() ) } );
            }
            static inline fp16 cosh( const fp16 &x )
            {
                return fp16( { cosh( x.raw() ) } );
            }
            static inline fp16 sinh( const fp16 &x )
            {
                return fp16( { sinh( x.raw() ) } );
            }
            static inline fp16 tanh( const fp16 &x )
            {
                return fp16( { tanh( x.raw() ) } );
            }
            static inline fp16 pow( const fp16 &x, const fp16 &y )
            {
                return fp16( { pow( x.raw(), y.raw() ) } );
            }
            static inline char* toASCII( const fp16 &x, char *str, int decimals )
            {
                return toASCII( x.raw(), str, decimals );
            }

    };

    inline std::ostream& operator<<( std::ostream& os, const fp16& obj )
    {
        char buff[ 64 ] = { 0 };
        os << fp16::toASCII( obj, buff, static_cast<int>( os.precision() ) );
        return os;
    }
    /*cstat -MISRAC++2008-5-0-7 -CERT-FLP34-C -MISRAC++2008-5-0-9*/
    constexpr fp16 operator"" _fp( long double val )
    {
        return { { static_cast<fp16Raw_t>( ( ( static_cast<double>( val )*static_cast<double>( fp16::one ) ) >= 0.0 ) ? ( static_cast<double>( val )*static_cast<double>( fp16::one ) ) + 0.5 :  ( static_cast<double>( val )*static_cast<double>( fp16::one ) ) - 0.5 ) } };
    }
    constexpr fp16 operator"" _fp(unsigned long long val)
    {
        return { { static_cast<fp16Raw_t>( static_cast<uint32_t>( val ) << 16 ) } };
    }
    /*cstat +MISRAC++2008-5-0-7 -CERT-FLP34-C +MISRAC++2008-5-0-9*/


}

#endif /*QLIBS_FP16*/