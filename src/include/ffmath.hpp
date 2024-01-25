/*!
 * @file ffmath.hpp
 * @author J. Camilo Gomez C.
 * @version 1.06
 * @note This file is part of the qLibs++ distribution.
 * @brief Fast floating-point math library for applications where speed is more
 * important than accuracy
 **/

#ifndef QLIBS_FFMATH
#define QLIBS_FFMATH

#include <include/qlibs_types.hpp>

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {

    /** @brief Namespace for the Fast floating-point math library*/
    namespace ffmath {


        /** @addtogroup  qffmath Float Fast-Math
        * @brief Fast floating-point math library for applications where speed is
        * more important than accuracy
        *  @{
        */

       /**
        * @brief  Enum with the possible categorizations of a 32-bit floating-point number
        */
        enum class classification {
            FFP_ZERO = 0,   /*!< Indicates that the value is positive or negative zero */
            FFP_SUBNORMAL,  /*!< Indicates that the value is subnormal  */
            FFP_NORMAL,     /*!< Indicates that the value is normal, i.e. not an infinity, subnormal, not-a-number or zero*/
            FFP_INFINITE,   /*!< Indicates that the value is not representable by the underlying type, positive or negative infinity*/
            FFP_NAN,        /*!< Indicates that the value is not-a-number @c nan (NaN)*/
        };


        /**
        * @brief Determines if the parameters given as floating-point values are
        * approximately equal.
        * @param[in] a Input to be compared.
        * @param[in] b Input to be compared.
        * @param[in] tol Tolerance
        * @return @c true when both values are approximately equal.
        */
        bool isEqual( const float a,
                      const float b,
                      const float tol = 1.175494351e-38F ) noexcept;

        /**
        * @brief Returns positive infinity @c inf as a 32-bit floating point number
        * @return The @c +inf value
        */
        float getInf( void );

        /**
        * @brief Returns Not a Number (NaN) @c nan as a 32-bit floating point number
        * @return The @c nan value
        */
        float getNan( void );

        /**
        * @brief Categorizes the floating-point number @a x. This function
        * determines whether its argument is a normal floating-point number, or one
        * of several special categories of values, including NaN (not a number),
        * infinity, subnormal floating-point values or zero. To determine what
        * category the argument belongs to, compare the return value with the
        * any of the following number classification macros:
        * - classification::FFP_ZERO
        * - classification::FFP_SUBNORMAL
        * - classification::FFP_NORMAL
        * - classification::FFP_INFINITE
        * - classification::FFP_NAN
        * @param[in] f The number you want to test.
        * @return One of the items in the ffmath::classification enumeration,
        * specifying the category of @a x.
        */
        classification classify( const float f );

        /**
        * @brief Determine if @a x is Not-A-Number (NaN)
        * @param[in] x The number you want to test.
        * @return @c true if the value of @a x is (NaN), otherwise
        * returns @c false.
        */
        inline bool isNan( const float x )
        {
            return ( classification::FFP_NAN == classify( x ) );
        }

        /**
        * @brief Determine if @a x is Infinity.
        * @param[in] x The number you want to test.
        * @return @c true if the value of @a x is ±Infinity, otherwise returns @c false.
        */
        inline bool isInf( const float x )
        {
            return ( classification::FFP_INFINITE == classify( x ) );
        }

        /**
        * @brief Determines if the given floating point number @a x has finite
        * value i.e. it is normal, subnormal or zero, but not infinite or NaN.
        * @param[in] x The number you want to test.
        * @return @c true if @a x has a finite value, @c false otherwise
        */
        inline bool isFinite( const float x )
        {
            return ( classify( x ) < classification::FFP_INFINITE );
        }

        /**
        * @brief Determines if the given floating point number @a x is normal, i.e.
        * is neither zero, subnormal, infinite, or NaN.
        * @param[in] x The number you want to test.
        * @return @c true if @a x has a normal value, @c false otherwise
        */
        inline bool isNormal( const float x )
        {
            return ( classification::FFP_NORMAL == classify( x ) );
        }

        /**
        * @brief Computes the sign function ( signum function).
        * @param[in] x The floating point value
        * @return The sign function of @a x
        */
        float sign( float x );

        /**
        * @brief Computes the absolute value of a floating point value @a x.
        * @param[in] x The floating point value
        * @return The absolute value of @a x
        */
        float absf( float x );

        /**
        * @brief Computes the multiplicative inverse or reciprocal for the value
        * @a x, denoted by @c 1/x or <tt>x^(−1)</tt>
        * @param[in] x The floating point value
        * @return The reciprocal value of @a x
        */
        float recip( float x );

        /**
        * @brief Computes the square-root of @a x
        * @param[in] x The floating point value
        * @return If all validations are ok, square root of @a x, is returned.
        * If the domain validation fails, @c nan is returned
        */
        float sqrt( float x );

        /**
        * @brief Computes the reciprocal square-root of @a x denoted as
        * <tt>1/sqrt(x)</tt>
        * @param[in] x The floating point value
        * @return If all validations are ok, the reciprocal square root of @a x, is
        * returned. If the domain validation fails @c nan is returned
        */
        float rSqrt( float x );

        /**
        * @brief Computes the cubic-root of @a x
        * @param[in] x The floating point value
        * @return If all validations are ok, cubic root of @a x, is returned. If
        * the domain validation fails, @c nan is returned
        */
        float cbrt( float x );

        /**
        * @brief Computes the reciprocal cubic-root of @a x denoted as
        * <tt>1/cbrt(x)</tt>
        * @param[in] x The floating point value
        * @return If all validations are ok, the reciprocal cubic root of @a x, is
        * returned. If the domain validation fails @c nan is returned
        */
        float rCbrt( float x );

        /**
        * @brief Computes the nearest integer value to @a x (in floating-point
        * format), rounding halfway cases away from zero.
        * @param[in] x The floating point value
        * @return The nearest integer value to @a x, rounding halfway cases away
        * from zero
        */
        float rounding( float x );

        /**
        * @brief Computes the largest integer value not greater than @a x.
        * @param[in] x The floating point value
        * @return The largest integer value not greater than @a x
        */
        float floor( float x );

        /**
        * @brief Computes the smallest integer value not less than @a x.
        * @param[in] x The floating point value
        * @return The smallest integer value not less than @a x
        */
        float ceil( float x );

        /**
        * @brief  Computes the nearest integer not greater in magnitude than @a x.
        * @param[in] x The floating point value
        * @return The nearest integer value not greater in magnitude than @a x
        * (in other words, @a x rounded towards zero)
        */
        float trunc( float x );

        /**
        * @brief Obtain the fractional part of @a x.
        * @param[in] x The floating point value
        * @return The fractional part of @a x
        */
        float frac( float x );

        /**
        * @brief Computes the floating point remainder after division of @a x
        * by @a y, where  @a x is the dividend and @a y is the divisor. This
        * function is often called the remainder operation, which can be
        * expressed as @c r=a-(b*trunc(a/b)) . This function follows the
        * convention that @c rem(x,0) is @c nan.
        * @note The concept of remainder after division is not uniquely defined,
        * and the two functions ffmath::mod() and ffmath::rem() each compute a
        * different variation.
        * The ffmath::mod() function produces a result that is either zero or
        * has the same sign as the divisor. The ffmath::rem() function produces
        * a result that is either zero or has the same sign as the dividend.
        * Another difference is the convention when the divisor is zero. The
        * ffmath::mod() function follows the convention that @c mod(x,0)
        * returns @c x, whereas the rem function follows the convention that
        * @c rem(x,0) returns @c nan.
        * @param[in] x The floating point value
        * @param[in] y The floating point value
        * @return If successful, returns the IEEE floating-point remainder of the
        * division @c x/y. If the domain validation fails, a @c nan value is returned.
        */
        float rem( float x, float y );

        /**
        * @brief Computes the floating point remainder after division of @a x
        * by @a y, where  @a x is the dividend and @a y is the divisor. This
        * function is often called the modulo operation, which can be expressed
        * as @c b=a-m*floor(a/m) . This function follows the convention that
        * @c mod(x,0) returns x.
        * @note The concept of remainder after division is not uniquely defined,
        * and the two functions ffmath::mod() and ffmath::rem() each compute a
        * different variation.
        * The ffmath::mod() function produces a result that is either zero or
        * has the same sign as the divisor. The ffmath::rem() function produces
        * a result that is either zero or has the same sign as the dividend.
        * Another difference is the convention when the divisor is zero. The
        * ffmath::mod() function follows the convention that @c mod(x,0)
        * returns @c x, whereas the rem function follows the convention that
        * @c rem(x,0) returns @c nan.
        * @param[in] x The floating point value
        * @param[in] y The floating point value
        * @return If successful, returns the IEEE floating-point remainder of the
        * division @c x/y. If the domain validation fails, a @c nan value is
        * returned.
        */
        float mod( float x, float y );

        /**
        * @brief Computes the sine of @a x (measured in radians).
        * @param[in] x The floating point value
        * @return If all validations are ok, the sine of @a x @c sin(x) in the
        * range  [-1 ; +1], is returned. If the domain validation fails, a @c nan
        * value is returned.
        */
        float sin( float x );

        /**
        * @brief Computes the cosine of @a x (measured in radians).
        * @param[in] x The floating point value
        * @return If all validations are ok, the cosine of @a x @c cos(x) in the range
        * [-1 ; +1], is returned. If the domain validation fails, a @c nan value is
        * returned.
        */
        float cos( float x );

        /**
        * @brief Computes the tangent of @a x (measured in radians).
        * @param[in] x The floating point value
        * @return If all validations are ok, the tangent of @a x @c tan(x) is
        * returned. If the domain validation fails, a @c nan value is returned.
        */
        float tan( float x );

        /**
        * @brief Computes the principal value of the arc sine of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the arc sine of @a x @c arcsin(x) in the range
        * [-pi/2 ; pi/2]. is returned.If the domain validation fails, a @c nan value is
        * returned.
        */
        float asin( float x );

        /**
        * @brief Computes the principal value of the arc cosine of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the arc cosine of @a x @c arccos(x) in the
        * range [0 ; pi]. is returned.If the domain validation fails, a @c nan value is
        * returned.
        */
        float acos( float x );

        /**
        * @brief Computes the principal value of the arc tangent of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the arc sine of @a x @c arctan(x) in the range
        * [-pi/2 ; pi/2]. is returned.If the domain validation fails, a @c nan value is
        * returned.
        */
        float atan( float x );

        /**
        * @brief Computes the arc tangent of y/x using the signs of arguments to
        * determine the correct quadrant.
        * @param[in] y The floating point value
        * @param[in] x The floating point value
        * @return If all validations are ok, the arc tangent of @c y/x <tt>arctan(y/x)</tt>
        * in the range [-pi ; +pi] radians, is returned. If the domain validation fails,
        * a @c nan value is returned.
        */
        float atan2( float y, float x );

        /**
        * @brief Computes 2 raised to the given power @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the base-2 exponential of @a x <tt>2^x</tt> is
        * returned. If the range validation fails due to overflow, @c inf is
        * returned.
        */
        float exp2( float x );

        /**
        * @brief Computes the base 2 logarithm of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the base-2 logarithm of @a x @c log_2(x) is
        * returned. If the domain validation fails, a @c nan value is returned.
        * If the pole validation fails, @c -inf is returned.
        */
        float log2( float x );

        /**
        * @brief Computes the e (Euler's number, 2.7182818) raised to the given
        * power @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the base-e exponential of @a x <tt>e^(x)</tt>
        * is returned. If the range validation fails due to overflow, @c +inf is
        * returned.
        */
        float exp( float x );

        /**
        * @brief Computes the value of 10 raised to the power of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the base-e exponential of @a x <tt>10^(x)</tt>
        * is returned. If the range validation fails due to overflow, @c ±inf or
        * @c nan is returned.
        */
        float exp10( float x );

        /**
        * @brief Computes the natural (base e) logarithm of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the natural (base-e) logarithm of @a x
        * @c ln(x) is returned. If the domain validation fails, a @c nan value is
        * returned. If the pole validation fails, @c -inf is returned.
        */
        float log( float x );

        /**
        * @brief Computes the common (base-10) logarithm of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the common (base-10) logarithm of @a x
        * @c log_10(x) is returned. If the domain validation fails, a @c nan value is
        * returned. If the pole validation fails, @c -inf is returned.
        */
        float log10( float x );

        /**
        * @brief Computes the value of @a b raised to the power @a e.
        * @param[in] b Base as floating point value
        * @param[in] e Exponent as floating point value
        * @return If all validations are ok, @a b raised to the power of @a e @c b^e is
        * returned. If the domain validation fails, a @c nan value is returned. If
        * the pole or range validation fails due to overflow, @c ±inf is
        * returned.
        */
        float pow( float b, float e );

        /**
        * @brief Computes hyperbolic sine of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the hyperbolic sine of @a x @c sinh(x) is
        * returned. If the range validation fails, a @c ±inf is value is
        * returned.
        */
        float sinh( float x );

        /**
        * @brief Computes hyperbolic cosine of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the hyperbolic cosine of @a x @c cosh(x) is
        * returned. If the range validation fails, a @c inf value is returned.
        */
        float cosh( float x );

        /**
        * @brief Computes hyperbolic tangent of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the hyperbolic tangent of @a x
        * @c tanh(x) is returned.
        */
        float tanh( float x );

        /**
        * @brief Computes the inverse hyperbolic sine of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the inverse hyperbolic sine of @a x
        * <tt>sinh^-1(x)</tt> is returned.
        */
        float asinh( float x );

        /**
        * @brief Computes the inverse hyperbolic cosine of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the inverse hyperbolic cosine of @a x
        * <tt>cosh^-1(x)</tt> is returned.
        */
        float acosh( float x );

        /**
        * @brief Computes the inverse hyperbolic tangent of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the inverse hyperbolic tangent of @a x
        * <tt>tanh^-1(x)</tt> is returned. If the domain validation fails, a @c nan
        * value is returned. If the pole validation fails, @c ±inf is returned.
        */
        float atanh( float x );

        /**
        * @brief Wraps angle @a x, in radians, to the interval [−pi, pi] such that
        * pi maps to pi and −pi maps to −pi. In general, odd, positive multiples
        * of pi map to pi and odd, negative multiples of pi map to −pi.
        * @param x The angle in radians
        * @return The angle @a x wrapped to the [-pi, pi] interval
        */
        float wrapToPi( float x );

        /**
        * @brief Wraps angle @a x, in radians, to the interval [0, 2*pi] such that
        * 0 maps to 0 and 2*pi and 2*pi maps to 2*pi. In general, positive multiples
        * of 2*pi map to 2*pi and negative multiples of 2*pi map to 0.
        * @param x The angle in radians
        * @return The angle @a x wrapped to the [0, 2*pi] interval
        */
        float wrapTo2Pi( float x );

        /**
        * @brief Wraps angle @a x, in degrees, to the interval [–180, 180] such
        * that 180 maps to 180 and –180 maps to –180. In general, odd, positive
        * multiples of 180 map to 180 and odd, negative multiples of 180 map to –180.
        * @param x The angle in degrees
        * @return The angle @a x wrapped to the [-pi, pi] interval
        */
        float wrapTo180( float x );

        /**
        * @brief Wraps angle @a x, in degrees, to the interval [0, 360] such
        * that 0 maps to 0 and 360 maps to 360. In general, positive multiples
        * of 360 map to 360 and negative multiples of 360 map to zero.
        * @param x The angle in degrees
        * @return The angle @a x wrapped to the [0, 360] interval
        */
        float wrapTo360( float x );

        /**
        * @brief Computes the error function of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the value of the error function is
        * returned.
        */
        float erf( float x );

        /**
        * @brief Computes the complementary error function of @a x.
        * @param[in] x The floating point value
        * @return If all validations are ok, the value of the complementary error
        * function is returned.
        */
        float erfc( float x );

        /**
        * @brief Decomposes given floating point value @a x into a normalized
        * fraction and an integral power of two.
        * @param[in] x The floating point value
        * @param[in] pw2 Pointer to integer value to store the exponent to
        * @return If @a x is zero, returns zero and stores zero in @a pw2. Otherwise
        * (if @a x is not zero), if all validations are ok, returns the value @a y in the
        * range (-1;-0.5], [0.5; 1) and stores an integer value in @a pw2 such that
        * <tt>y×2^(pw2) = x</tt> . If the value to be stored in @a pw2 is outside
        * the range of an @c int, the behavior is unspecified. If @a x is not a
        * floating-point number, the behavior is unspecified.
        */
        float rexp( float x, int32_t *pw2 );

        /**
        * @brief Multiplies a floating point value @a x by the number 2 raised to
        * the @a pw2 power.
        * @param[in] x The floating point value
        * @param[in] pw2 Integer value
        * @return @a x multiplied by 2 to the power of @a pw2
        * <tt>x×2^pwd</tt> is returned. If the range validation fails due to
        * overflow, @c ±inf is returned. If the range validation fails
        * due to underflow, the correct result (after rounding) is returned.
        */
        float ldexp( float x, int32_t pw2 );

        /**
        * @brief Computes the square root of the sum of the squares of @a x and @a y,
        * without undue overflow or underflow at intermediate stages of the
        * computation.
        * @param[in] x The floating point value
        * @param[in] y The floating point value
        * @return The hypotenuse of a right-angled triangle,
        * <tt>sqrt(x^2 + y^2)</tt>, is returned. If the range validation fails
        * due to overflow, @c +inf is returned. If the range validation
        * fails due to underflow, the correct result
        * (after rounding) is returned.
        */
        float hypot( float x, float y );

        /**
        * @brief Returns the next representable value of @a x in the direction of
        * @a y. If @a x equals to @a y, @a y is returned.
        * @param[in] x The floating point value
        * @param[in] y The floating point value
        * @return The next representable value of @a x in the
        * direction of @a y is returned. If @a x equals @a y, then @a y is returned.
        * If the range validation fails due to overflow, @c ±inf is returned
        * (with the same sign as @a x). If the range validation fails due to underflow,
        * the correct result is returned.
        */
        float nextAfter( float x, float y );

        /**
        * @brief Computes the gamma function of @a x
        * @param[in] x The floating point value
        * @return Upon successful completion, this function shall return Gamma(x).
        * If @a x is a negative integer, a Inf() value shall be returned. If the
        * correct value would cause overflow, tgamma() shall return ±Inf,
        * with the same sign as the correct value of the function.
        * If x is NaN, a NaN shall be returned.
        * If x is +Inf, x shall be returned.
        * If x is ±0, tgamma() shall return ±Inf.
        * If x is -Inf, a NaN  value shall be returned.
        * For IEEE Std 754-1985 float, overflow happens when 0 < x < 1/FLT_MAX,
        * and 171.7 < x.
        */
        float tgamma( float x );

        /**
        * @brief Computes the natural logarithm of the absolute value of the
        * gamma function of @a x
        * @note The argument @a x need not be a non-positive integer ( is defined
        * over the reals, except the non-positive integers).
        * @param[in] x The floating point value
        * @return Upon successful completion, this function shall return the
        * logarithmic gamma of x. If x is a non-positive integer, lgamma() shall
        * return +Inf. If the correct value would cause overflow, lgamma() shall
        * return ±Inf having the same sign as the correct value.
        * If x is NaN, a NaN shall be returned.
        * If x is 1 or 2, +0 shall be returned.
        * If x is ±Inf, +Inf shall be returned.
        */
        float lgamma( float x );

        /**
        * @brief Return the factorial of the integer part of @a x.
        * @note The argument @a x needs to be positive
        * @note This function overflows when @a x > 34
        * @param[in] x The floating point value
        * @return Upon successful completion, this function shall return the
        * factorial of the integer part of x. If x is non-positive, factorial() shall
        * return NaN. If the correct value would cause overflow, lgamma() shall
        * return +Inf.
        */
        float factorial( float x );

        /*cstat -MISRAC++2008-0-1-4_b*/

        /** @brief The base of natural logarithms ( e ) given as a single-precision floating-point number*/
        constexpr float FFP_E               = ( 2.71828182845904523540F );
        /** @brief The base 2 logarithm of e ( log_2 e ) given as a single-precision floating-point number */
        constexpr float FFP_LOG2E           = ( 1.44269504088896340740F );
        /** @brief The base 10 logarithm of e ( log_10 e ) given as a single-precision floating-point number */
        constexpr float FFP_LOG10E          = ( 0.43429448190325182765F );
        /** @brief The natural logarithm of 2 ( ln 2 ) given as a single-precision floating-point number */
        constexpr float FFP_LN2             = ( 0.69314718055994530942F );
        /** @brief The natural logarithm of 10 ( ln 10 ) given as a single-precision floating-point number */
        constexpr float FFP_LN10            = ( 2.30258509299404568402F );
        /** @brief The circumference of a circle with diameter 1, ( π ) given as a single-precision floating-point number */
        constexpr float FFP_PI              = ( 3.14159265358979323846F );
        /** @brief Half of π ( π/2 ) given as a single-precision floating-point number */
        constexpr float FFP_PI_2            = ( 1.57079632679489661923F );
        /** @brief A quarter of π ( π/4 ) given as a single-precision floating-point number */
        constexpr float FFP_PI_4            = ( 0.78539816339744830962F );
        /** @brief The inverse of π  ( 1/π ) given as a single-precision floating-point number */
        constexpr float FFP_1_PI            = ( 0.31830988618379067154F );
        /** @brief Twice the inverse of π  (  2/π ) given as a single-precision floating-point number */
        constexpr float FFP_2_PI            = ( 0.63661977236758134308F );
        /** @brief The inverse of the square root of π ( 2/√π ) given as a single-precision floating-point number */
        constexpr float FFP_2_SQRTPI        = ( 1.12837916709551257390F );
        /** @brief The square root of 2 ( √2 ) given as a single-precision floating-point number */
        constexpr float FFP_SQRT2           = ( 1.41421356237309504880F );
        /** @brief The inverse of square root of 2 ( 1/√2 ) given as a single-precision floating-point number */
        constexpr float FFP_SQRT1_2         = ( 0.70710678118654752440F );
        /** @brief The natural logarithm of the square root of 2π given as a single-precision floating-point number */
        constexpr float FFP_LN_SQRT_2PI     = ( 0.9189385332046727417803297F );



        /*cstat +MISRAC++2008-0-1-4_b*/

        /** @}*/
    }
}


#endif /*QLIBS_FFMATH*/