#ifndef QLIBS_FP16
#define QLIBS_FP16

#include "include/types.hpp" 
#include <iostream>

namespace qlibs {

    using fp16_t = int32_t;

    class fp16;
    class fp16Constant {
        private:
            fp16_t value{ 0 };
        public:
            constexpr fp16Constant() : value(0) {}
            constexpr fp16Constant( fp16_t x ) : value( x ) {}
            fp16_t raw() const {
                return value;
            }
            fp16Constant operator-() const {
                return fp16Constant( -value );
            }

            fp16 operator+( const fp16Constant& other ) const;
            fp16 operator-( const fp16Constant& other ) const;
            fp16 operator*( const fp16Constant& other ) const;
            fp16 operator/( const fp16Constant& other ) const;
    };
    class fp16 {
        friend class fp16Constant;
        private:
            fp16_t value{ overflow };
            static fp16_t Min; // skipcq: CXX-W2009
            static fp16_t Max; // skipcq: CXX-W2009
            static bool rounding; // skipcq: CXX-W2009
            static bool saturation; // skipcq: CXX-W2009
            
            static const fp16_t exp_max;
            static const fp16_t f2;
            static const fp16_t f3;
            static const fp16_t f16;
            static const fp16_t f100;
            static const fp16_t f6_5;
            static const float one_fp16_f;
            static const double one_fp16_d;
            static const uint32_t overflow_mask;
            static const uint32_t fraction_mask;
            static const uint32_t integer_mask;

            static uint32_t overflowCheck( uint32_t res, const uint32_t x, const uint32_t y );
            static fp16_t saturate( const fp16_t nsInput, const fp16_t x, const fp16_t y );
            static fp16_t fromInt( const int x );
            static fp16_t fromFloat( const float x );
            static fp16_t fromDouble( const double x );

            static fp16_t add( const fp16_t X, const fp16_t Y );
            static fp16_t sub( const fp16_t X, const fp16_t Y );
            static fp16_t mul( const fp16_t x, const fp16_t y );
            static fp16_t div( const fp16_t x, const fp16_t y );

            static fp16_t abs( fp16_t x );
            static fp16_t sqrt( fp16_t x );
            static fp16_t exp( fp16_t x );
            static fp16_t log( fp16_t x );
            static fp16_t log2( fp16_t x );
            static fp16_t wrapToPi( fp16_t x );
            static fp16_t wrapTo180( fp16_t x );
            static fp16_t sin( fp16_t x );
            static fp16_t cos( fp16_t x );
            static fp16_t tan( fp16_t x );
            static fp16_t atan2( fp16_t y, fp16_t x );
            static fp16_t atan( fp16_t x );
            static fp16_t asin( fp16_t x );
            static fp16_t acos( fp16_t x );
            static fp16_t cosh( fp16_t x );
            static fp16_t sinh( fp16_t x );
            static fp16_t tanh( fp16_t x );
            static fp16_t powi( fp16_t x, fp16_t y );
            static fp16_t pow( fp16_t x, fp16_t y );

            static char* itoa( char *buf, uint32_t scale, uint32_t val, uint8_t skip );
            static char* toASCII( const fp16_t num, char *str, int decimals );

            static fp16_t rs( fp16_t x );
            static fp16_t log2i( fp16_t x );

        public:
            virtual ~fp16() {}
            fp16() = default;
            fp16( const fp16Constant &x )
            {
                value = x.raw();
            }
            explicit fp16( const int x )
            {
                value = fromInt( x );
            }
            explicit fp16( const float x )
            {
                value = fromFloat( x );
            }
            explicit fp16( const double x )
            {
                value = fromDouble( x );
            }
            static const fp16_t f_e;
            static const fp16_t f_log2e;
            static const fp16_t f_log10e;
            static const fp16_t f_ln2;
            static const fp16_t f_ln10;
            static const fp16_t f_pi;
            static const fp16_t f_pi_2;
            static const fp16_t f_2pi;
            static const fp16_t f_pi_4;
            static const fp16_t f_1_pi;
            static const fp16_t f_2_pi;
            static const fp16_t f_2_sqrtpi;
            static const fp16_t f_sqrt2;
            static const fp16_t f_sqrt1_2;
            static const fp16_t epsilon;
            static const fp16_t MaxValue;
            static const fp16_t overflow;
            static const fp16_t one;
            static const fp16_t one_half;
            static const fp16_t f_180_pi;
            static const fp16_t f_pi_180;
            static const fp16_t f_180;
            static const fp16_t f_360;

            fp16_t raw( void ) const
            {
                return value;
            }

            fp16 operator+( const fp16 &other )
            {
                fp16 result;
                result.value = add( value, other.value );
                return result;
            }
            fp16& operator+=( const fp16 &other ) {
                value = add( value, other.value );
                return *this;
            }
            fp16& operator+=( const fp16Constant &other ) {
                value = add( value, other.raw() );
                return *this;
            }
            fp16 operator-() const {
                fp16 result;
                result.value = -value;
                return result;
            }
            fp16& operator-=( const fp16 &other ) {
                value = sub( value, other.value );
                return *this;
            }
            fp16& operator-=( const fp16Constant &other ) {
                value = sub( value, other.raw() );
                return *this;
            }
            fp16 operator-( const fp16 &other )
            {
                fp16 result;
                result.value = sub( value, other.value );
                return result;
            }
            fp16 operator*( const fp16 &other )
            {
                fp16 result;
                result.value = mul( value, other.value );
                return result;
            }
            fp16& operator*=( const fp16  &other ) {
                value = mul( value, other.value );
                return *this;
            }
            fp16& operator*=( const fp16Constant &other ) {
                value = mul( value, other.raw() );
                return *this;
            }
            fp16 operator/( const fp16 &other )
            {
                fp16 result;
                result.value = div( value, other.value );
                return result;
            }
            fp16 operator/( const fp16Constant &other )
            {
                fp16 result;
                result.value = div( value, other.raw() );
                return result;
            }
            fp16& operator/=( const fp16  &other ) {
                value = div( value, other.value );
                return *this;
            }
            fp16& operator/=( const fp16Constant &other ) {
                value = div( value, other.raw() );
                return *this;
            }
            fp16& operator++()
            {
                value = add( value, one );
                return *this;
            }
            fp16 operator++(int)
            {
                fp16 temp = *this;
                value = add( value, one );;
                return temp;
            }
            fp16& operator--()
            {
                value = sub( value, one );
                return *this;
            }
            fp16 operator--(int)
            {
                fp16 temp = *this;
                value = sub( value, one );;
                return temp;
            }
            bool operator>(const fp16 &other) const
            {
                return value > other.value;
            }
            bool operator>=(const fp16 &other) const
            {
                return value >= other.value;
            }
            bool operator<(const fp16 &other) const
            {
                return value < other.value;
            }
            bool operator<=(const fp16 &other) const
            {
                return value <= other.value;
            }
            bool operator==(const fp16 &other) const
            {
                return value == other.value;
            }
            bool operator!=(const fp16 &other) const
            {
                return value != other.value;
            }

            fp16& operator=( int x )
            {
                value = fromInt( x );
                return *this;
            }
            fp16& operator=( float x )
            {
                value = fromFloat( x );
                return *this;
            }
            fp16& operator=( double x )
            {
                value = fromDouble( x );
                return *this;
            }
            fp16& operator=( const fp16Constant &x )
            {
                value = x.raw();
                return *this;
            }

            fp16( const fp16& other) : value( other.value ) {}
            fp16& operator=( const fp16 &other )
            {
                value = other.value;
                return *this;
            }
            static int toInt( const fp16 &x );
            static float toFloat( const fp16 &x );
            static double toDouble( const fp16 &x );

            inline static fp16 from( const int x )
            {
                return fp16( x );
            }
            inline static fp16 from( const float x )
            {
                return fp16( x );
            }
            inline static fp16 from( const double x )
            {
                return fp16( x );
            }
            static fp16 abs( const fp16 &x );
            static fp16 abs( const fp16Constant &x );
            static fp16 sqrt( const fp16 &x );
            static fp16 sqrt( const fp16Constant &x );
            static fp16 exp( const fp16 &x );
            static fp16 exp( const fp16Constant &x );
            static fp16 log( const fp16 &x );
            static fp16 log( const fp16Constant &x );
            static fp16 log2( const fp16 &x );
            static fp16 log2( const fp16Constant &x );
            static fp16 wrapToPi( const fp16 &x );
            static fp16 wrapToPi( const fp16Constant &x );
            static fp16 wrapTo180( const fp16 &x );
            static fp16 wrapTo180( const fp16Constant &x );
            static fp16 sin( const fp16 &x );
            static fp16 sin( const fp16Constant &x );
            static fp16 cos( const fp16 &x );
            static fp16 cos( const fp16Constant &x );
            static fp16 tan( const fp16 &x );
            static fp16 tan( const fp16Constant &x );
            static fp16 atan2( const fp16 &y, const fp16 &x );
            static fp16 atan2( const fp16 &y, const fp16Constant &x );
            static fp16 atan2( const fp16Constant &y, const fp16Constant &x );
            static fp16 atan2( const fp16Constant &y, const fp16 &x );
            static fp16 atan( const fp16 &x );
            static fp16 atan( const fp16Constant &x );
            static fp16 asin( const fp16 &x );
            static fp16 asin( const fp16Constant &x );
            static fp16 acos( const fp16 &x );
            static fp16 acos( const fp16Constant &x );
            static fp16 cosh( const fp16 &x );
            static fp16 cosh( const fp16Constant &x );
            static fp16 sinh( const fp16 &x );
            static fp16 sinh( const fp16Constant &x );
            static fp16 tanh( const fp16 &x );
            static fp16 tanh( const fp16Constant &x );
            static fp16 pow( const fp16 &x, const fp16 &y );
            static fp16 pow( const fp16 &x, const fp16Constant &y );
            static fp16 pow( const fp16Constant &x, const fp16Constant &y );
            static fp16 pow( const fp16Constant &x, const fp16 &y );

            static char* toASCII( const fp16 &x, char *str, int decimals );
            static char* toASCII( const fp16Constant &x, char *str, int decimals );

            friend fp16 operator+( const fp16 &obj, const fp16Constant &val );
            friend fp16 operator-( const fp16 &obj, const fp16Constant &val );
            friend fp16 operator-( const fp16Constant &val, const fp16 &obj );
            friend fp16 operator*( const fp16 &obj, const fp16Constant &val );
            friend fp16 operator*( const fp16Constant &val, const fp16 &obj );
            friend fp16 operator/( const fp16 &obj, const fp16Constant &val );
            friend fp16 operator/( const fp16Constant &val, const fp16 &obj );
    };

    /*cstat -CERT-EXP30-C_b*/
    inline fp16 operator+( const fp16 &obj, const fp16Constant &val )
    {
        fp16 result;
        result.value = fp16::add( obj.raw(), val.raw() );
        return result;
    }
    inline fp16 operator-( const fp16 &obj, const fp16Constant &val )
    {
        fp16 result;
        result.value = fp16::sub( obj.raw(), val.raw() );
        return result;
    }
    inline fp16 operator-( const fp16Constant &val, const fp16 &obj )
    {
        fp16 result;
        result.value = fp16::sub( val.raw(), obj.raw() );
        return result;
    }
    inline fp16 operator*( const fp16 &obj, const fp16Constant &val )
    {
        fp16 result;
        result.value = fp16::mul( obj.raw(), val.raw() );
        return result;
    }
    inline fp16 operator*( const fp16Constant &val, const fp16 &obj )
    {
        fp16 result;
        result.value = fp16::mul( obj.raw(), val.raw() );
        return result;
    }
    inline fp16 operator/( const fp16 &obj, const fp16Constant &val )
    {
        fp16 result;
        result.value = fp16::div( obj.raw(), val.raw() );
        return result;
    }
    inline fp16 operator/( const fp16Constant &val, const fp16 &obj )
    {
        fp16 result;
        result.value = fp16::div( val.raw(), obj.raw() );
        return result;
    }
    /*cstat +CERT-EXP30-C_b*/
    inline bool operator>( const fp16& x, fp16Constant val )
    {
        return x.raw() > val.raw();
    }
    inline bool operator>=( const fp16& x, fp16Constant val )
    {
        return x.raw() >= val.raw();
    }
    inline bool operator<( const fp16& x, fp16Constant val )
    {
        return x.raw() < val.raw();
    }
    inline bool operator<=( const fp16& x, fp16Constant val )
    {
        return x.raw() <= val.raw();
    }
    inline bool operator==( const fp16& x, fp16Constant val )
    {
        return x.raw() == val.raw();
    }
    inline bool operator!=( const fp16& x, fp16Constant val )
    {
        return x.raw() != val.raw();
    }

    inline std::ostream& operator<<( std::ostream& os, const fp16& obj )
    {
        char buff[ 64 ] = { 0 };
        os << fp16::toASCII( obj, buff, static_cast<int>( os.precision() ) );
        return os;
    }
    inline std::ostream& operator<<( std::ostream& os, const fp16Constant& obj )
    {
        char buff[ 64 ] = { 0 };
        os << fp16::toASCII( obj, buff, static_cast<int>( os.precision() ) );
        return os;
    }

    constexpr fp16Constant operator"" _fp( long double val )
    {
        /*cstat -CERT-EXP30-C_b -CERT-FLP34-C -MISRAC++2008-5-0-7*/
        return fp16Constant( static_cast<fp16_t>( ( ( static_cast<double>( val )*static_cast<double>( fp16::one ) ) >= 0.0 ) ? ( static_cast<double>( val )*static_cast<double>( fp16::one ) ) + 0.5 :  ( static_cast<double>( val )*static_cast<double>( fp16::one ) ) - 0.5 ) );
        /*cstat +CERT-EXP30-C_b +CERT-FLP34-C +MISRAC++2008-5-0-7*/
    }
    constexpr fp16Constant operator"" _fp(unsigned long long val)
    {
        /*cstat -MISRAC++2008-5-0-9*/
        return fp16Constant( static_cast<fp16_t>( static_cast<uint32_t>( val ) << 16 ) );
        /*cstat +MISRAC++2008-5-0-9*/
    }
}

#endif /*QLIBS_FP16*/