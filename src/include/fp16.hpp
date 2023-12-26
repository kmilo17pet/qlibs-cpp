/*!
 * @file fp16.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Fixed-Point math Q16.16 with rounding and saturated arithmetic.
 **/

#ifndef QLIBS_FP16
#define QLIBS_FP16

#include "include/types.hpp"
#include <iostream>

namespace qlibs {

    /** @addtogroup qfp16 Fixed-Point Q16.16 math
    * @brief Fixed-point Q16.16 math library
    *  @{
    */

    /*! @cond  */
    using fp16Raw_t = int32_t;
    /*! @endcond  */

    /**  @brief Fixed-point Q16.16 type with width of exactly 32 bits.*/
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
            static const uint32_t overflow_mask;
            static const uint32_t fraction_mask;
            static const uint32_t integer_mask;
            static const fp16Raw_t f_pi_2;
            static const fp16Raw_t overflow;
            static const fp16Raw_t one;
            static const fp16Raw_t one_half;

            static uint32_t overflowCheck( uint32_t res,
                                           const uint32_t x,
                                           const uint32_t y ) noexcept;
            static fp16Raw_t saturate( const fp16Raw_t nsInput,
                                       const fp16Raw_t x,
                                       const fp16Raw_t y ) noexcept;
            static fp16Raw_t fromInt( const int x ) noexcept;
            static fp16Raw_t fromFloat( const float x ) noexcept;
            static fp16Raw_t fromDouble( const double x ) noexcept;
            static fp16Raw_t add( const fp16Raw_t X,
                                  const fp16Raw_t Y ) noexcept;
            static fp16Raw_t sub( const fp16Raw_t X,
                                  const fp16Raw_t Y ) noexcept;
            static fp16Raw_t mul( const fp16Raw_t x,
                                  const fp16Raw_t y ) noexcept;
            static fp16Raw_t div( const fp16Raw_t x,
                                  const fp16Raw_t y ) noexcept;
            static fp16Raw_t abs( fp16Raw_t x ) noexcept;
            static fp16Raw_t ceil( fp16Raw_t x ) noexcept;
            static fp16Raw_t sqrt( fp16Raw_t x ) noexcept;
            static fp16Raw_t exp( fp16Raw_t x ) noexcept;
            static fp16Raw_t log( fp16Raw_t x ) noexcept;
            static fp16Raw_t log2( fp16Raw_t x ) noexcept;
            static fp16Raw_t radToDeg( const fp16Raw_t x );
            static fp16Raw_t degToRad( const fp16Raw_t x );
            static fp16Raw_t wrapToPi( fp16Raw_t x ) noexcept;
            static fp16Raw_t wrapTo180( fp16Raw_t x ) noexcept;
            static fp16Raw_t sin( fp16Raw_t x ) noexcept;
            static fp16Raw_t cos( fp16Raw_t x ) noexcept;
            static fp16Raw_t tan( fp16Raw_t x ) noexcept;
            static fp16Raw_t atan2( fp16Raw_t y,
                                    fp16Raw_t x ) noexcept;
            static fp16Raw_t atan( fp16Raw_t x ) noexcept;
            static fp16Raw_t asin( fp16Raw_t x ) noexcept;
            static fp16Raw_t acos( fp16Raw_t x ) noexcept;
            static fp16Raw_t cosh( fp16Raw_t x ) noexcept;
            static fp16Raw_t sinh( fp16Raw_t x ) noexcept;
            static fp16Raw_t tanh( fp16Raw_t x ) noexcept;
            static fp16Raw_t powi( fp16Raw_t x,
                                   fp16Raw_t y ) noexcept;
            static fp16Raw_t pow( fp16Raw_t x,
                                  fp16Raw_t y ) noexcept;
            static char* itoa( char *buf, 
                               uint32_t scale,
                               uint32_t val,
                               uint8_t skip ) noexcept;
            static char* toASCII( const fp16Raw_t num,
                                  char *str,
                                  int decimals ) noexcept;
            static fp16Raw_t rs( fp16Raw_t x ) noexcept;
            static fp16Raw_t log2i( fp16Raw_t x ) noexcept;

            constexpr fp16( fp16Hidden val ) : value( val.x ) {}
            friend constexpr fp16 operator"" _fp( long double val );
            friend constexpr fp16 operator"" _fp( unsigned long long val );
        public:
            constexpr fp16() : value( 0 ) {}
            fp16( const fp16& other) : value( other.value ) {}

            /**
            * @brief Check for (q16.16) fixed-point overflow
            * @return @c true if the fixed-point has overflowed otherwise @c false.
            */
            inline bool isOverflow( void ) const noexcept
            {
                return overflow == value;
            }

            /**
            * @brief Check for (q16.16) fixed-point fp16::exp() operation reaches
            * the @c EXP_MAX value 
            * @return @c true if the fixed-point has reached @c EXP_MAX,
            * otherwise @c false.
            */
            inline bool isExpMax( void ) const noexcept
            {
                return ( exp_max == value ) || ( -exp_max == value );
            }

            /**
            * @brief Get the (q16.16) raw integer value from x
            * @return The raw integer value that represents the fixed-point.
            */
            inline fp16Raw_t raw( void ) const noexcept
            {
                return value;
            }

            inline fp16 operator+( const fp16 &other ) noexcept
            {
                return fp16( { add( value, other.value ) } );
            }
            inline fp16& operator+=( const fp16 &other ) noexcept
            {
                value = add( value, other.value );
                return *this;
            }
            inline fp16 operator-() const noexcept
            {
                return fp16( { -value } );
            }
            inline fp16& operator-=( const fp16 &other ) noexcept
            {
                value = sub( value, other.value );
                return *this;
            }
            inline fp16 operator-( const fp16 &other ) noexcept
            {
                return fp16( { sub( value, other.value ) } );
            }
            inline fp16 operator*( const fp16 &other ) noexcept
            {
                return fp16( { mul( value, other.value ) } );
            }
            inline fp16& operator*=( const fp16 &other ) noexcept
            {
                value = mul( value, other.value );
                return *this;
            }
            inline fp16 operator/( const fp16 &other ) noexcept
            {
                return fp16( { div( value, other.value ) } );
            }
            inline fp16& operator/=( const fp16  &other ) noexcept
            {
                value = div( value, other.value );
                return *this;
            }
            inline fp16& operator++() noexcept
            {
                value = add( value, one );
                return *this;
            }
            inline fp16 operator++(int) noexcept
            {
                fp16 temp = *this;
                value = add( value, one );;
                return temp;
            }
            inline fp16& operator--() noexcept
            {
                value = sub( value, one );
                return *this;
            }
            inline fp16 operator--(int) noexcept
            {
                fp16 temp = *this;
                value = sub( value, one );;
                return temp;
            }
            inline bool operator>(const fp16 &other) const noexcept
            {
                return value > other.value;
            }
            inline bool operator>=(const fp16 &other) const noexcept
            {
                return value >= other.value;
            }
            inline bool operator<(const fp16 &other) const noexcept
            {
                return value < other.value;
            }
            inline bool operator<=(const fp16 &other) const noexcept
            {
                return value <= other.value;
            }
            inline bool operator==(const fp16 &other) const noexcept
            {
                return value == other.value;
            }
            inline bool operator!=(const fp16 &other) const noexcept
            {
                return value != other.value;
            }
            inline fp16& operator=( int x ) noexcept
            {
                value = fromInt( x );
                return *this;
            }
            inline fp16& operator=( float x ) noexcept
            {
                value = fromFloat( x );
                return *this;
            }
            inline fp16& operator=( double x ) noexcept
            {
                value = fromDouble( x );
                return *this;
            }
            inline fp16& operator=( const fp16 &other ) noexcept
            {
                value = other.value;
                return *this;
            }

            /**
            * @brief Returns the fixed-point value @a x converted to int.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns @a x converted to int.
            */
            static int toInt( const fp16 &x ) noexcept;

             /**
            * @brief Returns the fixed-point value @a x converted to float.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns @a x converted to float.
            */
            static float toFloat( const fp16 &x ) noexcept;

            /**
            * @brief Returns the fixed-point value @a x converted to double.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns @a x converted to double.
            */
            static double toDouble( const fp16 &x ) noexcept;

            /**
            * @brief Returns the int value @a x converted to fixed-point q16.16.
            * @param[in] x The integer value.
            * @return This function returns @a x converted to fixed-point q16.16.
            */
            static inline fp16 from( const int x ) noexcept
            {
                return fp16( { fromInt( x ) } );
            }

            /**
            * @brief Returns the float value @a x converted to fixed-point q16.16.
            * @param[in] x The floating-point value.
            * @return This function returns @a x converted to fixed-point q16.16.
            */
            static inline fp16 from( const float x ) noexcept
            {
                return fp16( { fromFloat( x ) } );
            }

            /**
            * @brief Returns the double value @a x converted to fixed-point q16.16.
            * @param[in] x The double precision floating-point value.
            * @return This function returns @a x converted to fixed-point q16.16.
            */
            static inline fp16 from( const double x ) noexcept
            {
                return fp16( { fromDouble( x ) } );
            }

            /**
            * @brief Returns the largest integer value less than or equal to @a x.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns the largest integral value not greater
            * than @a x.
            */
            static inline fp16 floor( const fp16 &x ) noexcept
            {
                /*cstat -MISRAC++2008-5-0-9*/
                return fp16( { static_cast<fp16Raw_t>( static_cast<uint32_t>( x.raw() ) & integer_mask ) } );
                /*cstat +MISRAC++2008-5-0-9*/
            }

            static inline fp16 ceil( const fp16 &x ) noexcept
            {
                return fp16( { ceil( x.raw() ) } );
            }

            /**
            * @brief Returns the nearest integer value of the fixed-point argument @a x
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns the nearest integral value of @a x.
            */
            static inline fp16 round( const fp16 &x ) noexcept
            {
                return fp16( { x.raw() + one } );
            }

            /**
            * @brief Returns the absolute value of @a x.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns the absolute value of x.
            */
            static inline fp16 abs( const fp16 &x ) noexcept
            {
                return fp16( { abs( x.raw() ) } );
            }

            /**
            * @brief Returns the fixed-point square root of @a x.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns the square root of @a x. For negative
            * numbers, returns @c overflow.
            */
            static inline fp16 sqrt( const fp16 &x ) noexcept
            {
                return fp16( { sqrt( x.raw() ) } );
            }

            /**
            * @brief Returns the fixed-point value of e raised to the xth power.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns the exponential value of x. @c overflow
            * when an operation overflow is detected.
            */
            static inline fp16 exp( const fp16 &x ) noexcept
            {
                return fp16( { exp( x.raw() ) } );
            }

            /**
            * @brief Returns the fixed-point natural logarithm (base-e logarithm) of @a x.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns natural logarithm of @a x. For negative
            * values returns @c overflow
            */
            static inline fp16 log( const fp16 &x ) noexcept
            {
                return fp16( { log( x.raw() ) } );
            }

            /**
            * @brief Returns the fixed-point log base 2 of @a x.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns log base 2 of @a x. For negative values
            * returns @c overflow
            */
            static inline fp16 log2( const fp16 &x ) noexcept
            {
                return fp16( { log2( x.raw() ) } );
            }

            /**
            * @brief Converts angle units from radians to degrees.
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in radians.
            * @return This function returns the angle converted in degrees.
            */
            static inline fp16 radToDeg( const fp16 &x ) noexcept
            {
                return fp16( { radToDeg( x.raw() ) } );
            }

            /**
            * @brief Converts angle units from degrees to radians.
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in degrees.
            * @return This function returns the angle converted in radians.
            */
            static inline fp16 degToRad( const fp16 &x ) noexcept
            {
                return fp16( { degToRad( x.raw() ) } );
            }

            /**
            * @brief Wrap the fixed-point angle in radians to [−pi pi]
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in radians.
            * @return This function returns the wrapped angle in the range [−pi, pi]
            * of @a x.
            */
            static inline fp16 wrapToPi( const fp16 &x ) noexcept
            {
                return fp16( { wrapToPi( x.raw() ) } );
            }

            /**
            * @brief Wrap the fixed-point angle in degrees  to [−180 180]
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in degrees.
            * @return This function returns the wrapped angle in the range [−180, 180]
            * of @a x.
            */
            static inline fp16 wrapTo180( const fp16 &x ) noexcept
            {
                return fp16( { wrapTo180( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point sine of the radian angle @a x.
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in radians.
            * @return This function returns sine of @a x.
            */
            static inline fp16 sin( const fp16 &x ) noexcept
            {
                return fp16( { sin( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point cosine of the radian angle @a x.
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in radians.
            * @return This function returns cosine of @a x.
            */
            static inline fp16 cos( const fp16 &x ) noexcept
            {
                return fp16( { cos( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point tangent  of the radian angle @a x.
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in radians.
            * @return This function returns tangent of @a x.
            */
            static inline fp16 tan( const fp16 &x ) noexcept
            {
                return fp16( { tan( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point arc tangent in radians of @a y / @a x
            * based on the signs of both values to determine the correct quadrant.
            * @param[in] y The fixed-point(q16.16) value representing an x-coordinate.
            * @param[in] x The fixed-point(q16.16) value representing an y-coordinate.
            * @return This function returns the principal arc tangent of y/x, in the
            * interval [-pi,+pi] radians.
            */
            static inline fp16 atan2( const fp16 &y,
                                      const fp16 &x ) noexcept
            {
                /*cstat -CERT-EXP30-C_b*/
                return fp16( { atan2( y.raw(), x.raw() ) } );
                /*cstat +CERT-EXP30-C_b*/
            }

            /**
            * @brief Computes the fixed-point arc tangent of @a x in radians.
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in radians.
            * @return This function returns arc tangent of @a x.
            */
            static inline fp16 atan( const fp16 &x ) noexcept
            {
                return fp16( { atan( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point arc sine of @a x in radians.
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in radians.
            * @return This function returns arc sine of @a x.
            */
            static inline fp16 asin( const fp16 &x ) noexcept
            {
                return fp16( { asin( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point arc cosine of @a x in radians.
            * @param[in] x The fixed-point(q16.16) value representing an angle expressed
            * in radians.
            * @return This function returns arc cosine of @a x.
            */
            static inline fp16 acos( const fp16 &x ) noexcept
            {
                return fp16( { acos( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point hyperbolic cosine of @a x.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns hyperbolic cosine of @a x. If overflow
            * detected returns @c overflow. If the function saturates, returns
            * @c EXP_MAX or @c EXP_MIN.
            */
            static inline fp16 cosh( const fp16 &x ) noexcept
            {
                return fp16( { cosh( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point hyperbolic sine of @a x.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns hyperbolic sine of @a x. If overflow
            * detected returns @c overflow. If the function saturates, returns
            * @c EXP_MAX or @c EXP_MIN.
            */
            static inline fp16 sinh( const fp16 &x ) noexcept
            {
                return fp16( { sinh( x.raw() ) } );
            }

            /**
            * @brief Computes the fixed-point hyperbolic tangent  of @a x.
            * @param[in] x The fixed-point(q16.16) value.
            * @return This function returns hyperbolic tangent of @a x. If overflow
            * detected returns @c overflow. If the function saturates, returns
            * @c EXP_MAX or @c EXP_MIN.
            */
            static inline fp16 tanh( const fp16 &x ) noexcept
            {
                return fp16( { tanh( x.raw() ) } );
            }

            /**
            * @brief Returns @a x raised to the power of @a y. (x^y)
            * @param[in] x The fixed-point(q16.16) base value.
            * @param[in] y The fixed-point(q16.16) power value.
            * @return This function returns the result of raising @a x to the power @a y.
            * @c overflow when an operation overflow is detected.
            */
            static inline fp16 pow( const fp16 &x,
                                    const fp16 &y ) noexcept
            {
                /*cstat -CERT-EXP30-C_b*/
                return fp16( { pow( x.raw(), y.raw() ) } );
                /*cstat +CERT-EXP30-C_b*/
            }

            /**
            * @brief Converts the fixed-point value to a formatted string.
            * @param[in] num The fixed-point(q16.16) value to be converted to string.
            * @param[in] str Array in memory where to store the resulting null-terminated
            * string.
            * @param[in] decimals Number of decimals to show in the string representation.
            * @note: Max decimal allowed = 5
            * @return A pointer to the resulting null-terminated string, same as
            * parameter @a str
            */
            static inline char* toASCII( const fp16 &x,
                                         char *str,
                                         int decimals ) noexcept
            {
                return toASCII( x.raw(), str, decimals );
            }

    };

    /*! @cond  */
    inline std::ostream& operator<<( std::ostream& os,
                                     const fp16& obj )
    {
        char buff[ 64 ] = { 0 };
        os << fp16::toASCII( obj, buff, static_cast<int>( os.precision() ) );
        return os;
    }
    /*cstat -MISRAC++2008-5-0-7 -CERT-FLP34-C -MISRAC++2008-5-0-9*/
    constexpr fp16 operator"" _fp( long double val )
    {
        return { { static_cast<fp16Raw_t>( ( ( static_cast<double>( val )*65536.0 ) >= 0.0 ) ? ( static_cast<double>( val )*65536.0 ) + 0.5 :  ( static_cast<double>( val )*65536.0 ) - 0.5 ) } };
    }
    constexpr fp16 operator"" _fp(unsigned long long val)
    {
        return { { static_cast<fp16Raw_t>( static_cast<uint32_t>( val ) << 16 ) } };
    }
    /*cstat +MISRAC++2008-5-0-7 -CERT-FLP34-C +MISRAC++2008-5-0-9*/
    /*! @endcond  */

    /**  @brief @c e The base of natural logarithms as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_E        = 2.718281828459045235360_fp;

    /**  @brief @c log2(e) The base @c 2 logarithm of @c e as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_LOG2E    = 1.442695040888963407360_fp;

    /**  @brief @c log10(e) The base @c 10 logarithm of @c e as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_LOG10E   = 0.434294481903251827651_fp;

    /**  @brief @c ln(2) The natural logarithm of @c 2 as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_LN2      = 0.693147180559945309417_fp;

    /**  @brief @c ln(10) The natural logarithm of @c 10 as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_LN10     = 2.302585092994045684020_fp;

    /**  @brief @c pi The circumference of a circle with diameter @c 1 as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_PI       = 3.141592653589793238460_fp;

    /**  @brief @c pi/2 Half of @c pi as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_PI_2     = 1.570796326794896619230_fp;

    /**  @brief @c pi/4 A quarter of @c pi as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_PI_4     = 0.785398163397448309616_fp;

    /**  @brief @c 1/pi The inverse of @c pi as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_1_PI     = 0.318309886183790671538_fp;

    /**  @brief @c 2/pi Twice the inverse of @c pi as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_2_PI     = 0.636619772367581343076_fp;

    /**  @brief @c 2/sqrt(pi) The inverse of the square root of @c pi as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_2_SQRTPI = 1.128379167095512573900_fp;

    /**  @brief @c sqrt(2) The square root of @c 2 as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_SQRT2    = 1.414213562373095048800_fp;

    /**  @brief @c 1/sqrt(2) The inverse of the square root of @c 2 as a Fixed-point Q16.16 value.*/
    constexpr fp16 FP_SQRT1_2  = 0.707106781186547524401_fp;

     /** @}*/

}

#endif /*QLIBS_FP16*/