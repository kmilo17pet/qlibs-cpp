/*!
 * @file ffmath.hpp
 * @author J. Camilo Gomez C.
 * @version 1.07
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


        /** @addtogroup qffmath Float Fast-Math
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
        * @brief Returns positive @a infinity @c inf as a 32-bit floating point number
        * @return The @c +inf value
        */
        float getInf( void );

        /**
        * @brief Returns Not a Number ( @a NaN ) @c nan as a 32-bit floating point number
        * @return The @c nan value
        */
        float getNan( void );

        /**
        * @brief Categorizes the floating-point number @a x. This function
        * determines whether its argument is a normal floating-point number, or one
        * of several special categories of values, including @a NaN (not a number),
        * @a infinity, subnormal floating-point values or zero. To determine what
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
        * @brief Returns the greater of the given values.
        * @param[in] x Value to compare.
        * @param[in] y Value to compare.
        * @return The greater of @a x and @a y. If they are equivalent,
        * returns @a y
        */
        template<typename T>
        inline T Max( const T x, const T y )
        {
            return ( x > y ) ? x : y;
        }

        /**
        * @brief Returns the smaller  of the given values.
        * @param[in] x Value to compare.
        * @param[in] y Value to compare.
        * @return The smaller of @a x and @a y. If they are equivalent,
        * returns @a y
        */
        template<typename T>
        inline T Min( const T x, const T y )
        {
            return ( x < y ) ? x : y;
        }

        /**
        * @brief Determine if @a x is Not-A-Number @a NaN
        * @param[in] x The number you want to test.
        * @return @c true if the value of @a x is @a NaN, otherwise
        * returns @c false.
        */
        inline bool isNan( const float x )
        {
            return ( classification::FFP_NAN == classify( x ) );
        }

        /**
        * @brief Determine if @a x is @a Infinity.
        * @param[in] x The number you want to test.
        * @return @c true if the value of @a x is ±Infinity, otherwise returns @c false.
        */
        inline bool isInf( const float x )
        {
            return ( classification::FFP_INFINITE == classify( x ) );
        }

        /**
        * @brief Determines if the given floating point number @a x has finite
        * value i.e. it is normal, subnormal or zero, but not @a infinite or @a NaN.
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
        * @brief  Composes a floating point value with the magnitude of @a mag
        * and the sign of @a sgn
        * @param[in] mag floating-point value.
        * @param[in] sgn floating-point value
        * @return The floating point value with the magnitude of @a mag and the
        * sign of @a sgn is returned. If @a mag is @c nan, then @c nan with the
        * sign of @a sgn is returned. if @c sgn is @c -0, the result is only
        * negative if the implementation supports the signed zero consistently
        * in arithmetic operations.
        */
        float copysign( float mag,
                        float sgn );

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
        * @return Upon successful completion, square root of @a x, is returned.
        * If the domain validation fails, @c nan is returned
        */
        float sqrt( float x );

        /**
        * @brief Computes the reciprocal square-root of @a x denoted as
        * <tt>1/sqrt(x)</tt>
        * @param[in] x The floating point value
        * @return Upon successful completion, the reciprocal square root of @a x, is
        * returned. If the domain validation fails @c nan is returned
        */
        float rSqrt( float x );

        /**
        * @brief Computes the cubic-root of @a x
        * @param[in] x The floating point value
        * @return Upon successful completion, cubic root of @a x, is returned. If
        * the domain validation fails, @c nan is returned
        */
        float cbrt( float x );

        /**
        * @brief Computes the reciprocal cubic-root of @a x denoted as
        * <tt>1/cbrt(x)</tt>
        * @param[in] x The floating point value
        * @return Upon successful completion, the reciprocal cubic root of @a x, is
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
        * expressed as <tt>r=a-(b*trunc(a/b))</tt> . This function follows the
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
        float rem( float x,
                   float y );

        /**
        * @brief Computes the floating point remainder after division of @a x
        * by @a y, where  @a x is the dividend and @a y is the divisor. This
        * function is often called the modulo operation, which can be expressed
        * as <tt>b=a-m*floor(a/m)</tt> . This function follows the convention that
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
        float mod( float x,
                   float y );

        /**
        * @brief Computes the sine of @a x (measured in radians).
        * @param[in] x The floating point value
        * @return Upon successful completion, the sine of @a x @c sin(x) in the
        * range <tt>[-1 ; +1]</tt>, is returned. If the domain validation fails, a @c nan
        * value is returned.
        */
        float sin( float x );

        /**
        * @brief Computes the cosine of @a x (measured in radians).
        * @param[in] x The floating point value
        * @return Upon successful completion, the cosine of @a x @c cos(x) in the range
        * <tt>[-1 ; +1]</tt>, is returned. If the domain validation fails, a @c nan value is
        * returned.
        */
        float cos( float x );

        /**
        * @brief Computes the tangent of @a x (measured in radians).
        * @param[in] x The floating point value
        * @return Upon successful completion, the tangent of @a x @c tan(x) is
        * returned. If the domain validation fails, a @c nan value is returned.
        */
        float tan( float x );

        /**
        * @brief Computes the principal value of the arc sine of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the arc sine of @a x @c arcsin(x) in the range
        * <tt>[-pi/2 ; pi/2]</tt>. is returned.If the domain validation fails, a @c nan value is
        * returned.
        */
        float asin( float x );

        /**
        * @brief Computes the principal value of the arc cosine of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the arc cosine of @a x @c arccos(x) in the
        * range <tt>[0 ; pi]</tt>. is returned.If the domain validation fails, a @c nan value is
        * returned.
        */
        float acos( float x );

        /**
        * @brief Computes the principal value of the arc tangent of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the arc sine of @a x @c arctan(x) in the range
        * <tt>[-pi/2 ; pi/2]</tt>. is returned.If the domain validation fails, a @c nan value is
        * returned.
        */
        float atan( float x );

        /**
        * @brief Computes the arc tangent of @c y/x using the signs of arguments to
        * determine the correct quadrant.
        * @param[in] y The floating point value
        * @param[in] x The floating point value
        * @return Upon successful completion, the arc tangent of @c y/x <tt>arctan(y/x)</tt>
        * in the range <tt>[-pi ; +pi]</tt> radians, is returned. If the domain validation fails,
        * a @c nan value is returned.
        */
        float atan2( float y,
                     float x );

        /**
        * @brief Computes 2 raised to the given power @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the base-2 exponential of @a x <tt>2^x</tt> is
        * returned. If the range validation fails due to overflow, @c inf is
        * returned.
        */
        float exp2( float x );

        /**
        * @brief Computes the base 2 logarithm of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the base-2 logarithm of @a x @c log_2(x) is
        * returned. If the domain validation fails, a @c nan value is returned.
        * If the pole validation fails, @c -inf is returned.
        */
        float log2( float x );

        /**
        * @brief Computes the @a e ( Euler's number, @c 2.7182818 ) raised to the given
        * power @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the base-e exponential of @a x <tt>e^(x)</tt>
        * is returned. If the range validation fails due to overflow, @c +inf is
        * returned.
        */
        float exp( float x );

        /**
        * @brief Returns @a e raised to the given power minus one <tt>e^x-1</tt>
        * power @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the base-e exponential of @a x minus
        * one <tt>e^(x)-1</tt> is returned. If the range validation fails due
        * to overflow, @c +inf is
        * returned.
        */
        float expm1( float x );

        /**
        * @brief Computes the value of 10 raised to the power of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the base-e exponential of @a x <tt>10^(x)</tt>
        * is returned. If the range validation fails due to overflow, @c ±inf or
        * @c nan is returned.
        */
        float exp10( float x );

        /**
        * @brief Computes the natural (base e) logarithm of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the natural (base-e) logarithm of @a x
        * @c ln(x) is returned. If the domain validation fails, a @c nan value is
        * returned. If the pole validation fails, @c -inf is returned.
        */
        float log( float x );

        /**
        * @brief Computes the natural (base e) logarithm of 1 plus the given
        * number @a x @c ln(1+x) .
        * @param[in] x The floating point value
        * @return Upon successful completion, the natural (base-e) logarithm
        * of 1 plus the given number @a x @c ln(1+x) is returned. If the domain
        * validation fails, a @c nan value is returned. If the pole validation
        * fails, @c -inf is returned.
        */
        float log1p( float x );

        /**
        * @brief Computes the common (base-10) logarithm of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the common (base-10) logarithm of @a x
        * @c log_10(x) is returned. If the domain validation fails, a @c nan value is
        * returned. If the pole validation fails, @c -inf is returned.
        */
        float log10( float x );

        /**
        * @brief Computes the value of @a b raised to the power @a e.
        * @param[in] b Base as floating point value
        * @param[in] e Exponent as floating point value
        * @return Upon successful completion, @a b raised to the power of @a e @c b^e is
        * returned. If the domain validation fails, a @c nan value is returned. If
        * the pole or range validation fails due to overflow, @c ±inf is
        * returned.
        */
        float pow( float b,
                   float e );

        /**
        * @brief Computes hyperbolic sine of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the hyperbolic sine of @a x @c sinh(x) is
        * returned. If the range validation fails, a @c ±inf is value is
        * returned.
        */
        float sinh( float x );

        /**
        * @brief Computes hyperbolic cosine of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the hyperbolic cosine of @a x @c cosh(x) is
        * returned. If the range validation fails, a @c inf value is returned.
        */
        float cosh( float x );

        /**
        * @brief Computes hyperbolic tangent of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the hyperbolic tangent of @a x
        * @c tanh(x) is returned.
        */
        float tanh( float x );

        /**
        * @brief Computes the inverse hyperbolic sine of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the inverse hyperbolic sine of @a x
        * <tt>sinh^-1(x)</tt> is returned.
        */
        float asinh( float x );

        /**
        * @brief Computes the inverse hyperbolic cosine of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the inverse hyperbolic cosine of @a x
        * <tt>cosh^-1(x)</tt> is returned.
        */
        float acosh( float x );

        /**
        * @brief Computes the inverse hyperbolic tangent of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the inverse hyperbolic tangent of @a x
        * <tt>tanh^-1(x)</tt> is returned. If the domain validation fails, a @c nan
        * value is returned. If the pole validation fails, @c ±inf is returned.
        */
        float atanh( float x );

        /**
        * @brief Wraps angle @a x, in radians, to the interval <tt>[−pi, pi]</tt> such that
        * @c pi maps to @c pi and @c −pi maps to @c −pi. In general, odd, positive multiples
        * of @c pi map to @c pi and odd, negative multiples of @c pi map to @c −pi.
        * @param x The angle in radians
        * @return The angle @a x wrapped to the <tt>[-pi, pi]</tt> interval
        */
        float wrapToPi( float x );

        /**
        * @brief Wraps angle @a x, in radians, to the interval <tt>[0, 2*pi]</tt> such that
        * @c 0 maps to @c 0 and @c 2*pi and @c 2*pi maps to @c 2*pi. In general, positive multiples
        * of @c 2*pi map to @c 2*pi and negative multiples of @c 2*pi map to @c 0.
        * @param x The angle in radians
        * @return The angle @a x wrapped to the <tt>[0, 2*pi]</tt> interval
        */
        float wrapTo2Pi( float x );

        /**
        * @brief Wraps angle @a x, in degrees, to the interval <tt>[–180, 180]</tt> such
        * that @c 180 maps to @c 180 and @c –180 maps to @c –180. In general, odd, positive
        * multiples of @c 180 map to @c 180 and odd, negative multiples of @c 180 map to @c –180.
        * @param x The angle in degrees
        * @return The angle @a x wrapped to the <tt>[-180, 180]</tt> interval
        */
        float wrapTo180( float x );

        /**
        * @brief Wraps angle @a x, in degrees, to the interval <tt>[0, 360]</tt> such
        * that @c 0 maps to @c 0 and @c 360 maps to @c 360. In general, positive multiples
        * of @c 360 map to @c 360 and negative multiples of @c 360 map to zero.
        * @param x The angle in degrees
        * @return The angle @a x wrapped to the <tt>[0, 360]</tt> interval
        */
        float wrapTo360( float x );

        /**
        * @brief Computes the midpoint of the floating-points @a  a and @a b.
        * @param a A floating point value
        * @param b A floating point value
        * @return Half the sum of @a a and @a b. No overflow occurs. A at most
        * one inexact operation occurs.
        */
        float midpoint( float a,
                        float b );

        /**
        * @brief  Computes the linear interpolation between @a a and @a b, if
        * the parameter t is inside [0, 1] (the linear extrapolation otherwise),
        * i.e. the result of <tt> a+t*(b-a) </tt> with accounting for floating-point
        * calculation imprecision.
        * @param a A floating point value
        * @param b A floating point value
        * @param t A floating point value
        * @return Return <tt> a+t*(b-a) </tt>. When both @a a and @a b are finite, the
        * following properties are guaranteed:
        * If @c t==0, the result is equal to @a a.
        * If @c t==1, the result is equal to @a b.
        * If @c t>=0 and @c t<=1, the result is finite.
        * If @a t is finite and @c a==b, the result is equal to @a a.
        * If @a t is finite or <tt> (a-b)!=0</tt> with @a t infinity, the result
        * is not @c nan.
        * Let @c CMP(x,y) be @c 1 if @c x>y, @c -1 if @c x<y, and @c 0 otherwise.
        * For any @c t1 and @c t2, the product of: <tt> CMP(lerp(a,b,t2)</tt>,
        * <tt>lerp(a,b,t1)</tt>, <tt> CMP(t2,t1)</tt> , and @c CMP(b,a) is non-negative.
        * (That is, lerp() is monotonic.).
        */
        float lerp( float a,
                    float b,
                    float t );

        /**
        * @brief Scales the given input @a x in value range given by  @a xMin and
        * @a xMax to value range specified by the @a yMin and @a yMax.
        * @param[in] x Input
        * @param[in] xMin Input minimum value for range
        * @param[in] xMax Input maximum  value for range
        * @param[in] yMin Output minimum value for range
        * @param[in] yMax Output maximum value for range
        * @return The scaled value in range @a yMin and @a yMax.
        */
        float map( const float x,
                   const float xMin,
                   const float xMax,
                   const float yMin,
                   const float yMax ) noexcept;
        /**
        * @brief Normalize the given input @a x in value range given by @a xMin and
        * @a xMax to value range between 0 and 1.
        * @param[in] x Input
        * @param[in] xMin Input minimum value for range
        * @param[in] xMax Input maximum  value for range
        * @return The scaled value in range [0 - 1].
        */
        real_t normalize( const float x,
                          const float xMin,
                          const float xMax ) noexcept;
        /**
        * @brief Determines if the value pointed by @a x falls within a range
        * specified by the upper limit and lower limit inputs and coerces the value
        * to fall within the range
        * @param[in,out] x Input
        * @param[in] lowerL Lower limit.
        * @param[in] upperL Upper limit.
        * @return @c true when the value falls within the specified range, otherwise
        * false
        */
        bool inRangeCoerce( float &x,
                            const float lowerL,
                            const float upperL ) noexcept;

        /**
        * @brief Determines if the point at ( @a x, @a y ) is inside the polygon
        * given by the set of points on  @a px and @a py.
        * @param[in] x Point x-coordinate
        * @param[in] y Point y-coordinate
        * @param[in] px x-coordinate points of the polygon
        * @param[in] py y-coordinate points of the polygon
        * @param[in] p Number of points that represent the polygon
        * @return @c true when the given point is inside the polygon
        */
        bool inPolygon( const float x,
                        const float y,
                        const float * const px,
                        const float * const py,
                        const size_t p ) noexcept;

        /**
        * @brief Determines if the point at ( @a x, @a y) is inside the circle
        * with radius @a r located at @a cx and @a cy.
        * @param[in] x Point x-coordinate
        * @param[in] y Point y-coordinate
        * @param[in] cx X-location of the circle
        * @param[in] cy Y-location of the circle
        * @param[in] r Radio of the circle
        * @return @c true when the given point is inside the circle
        */
        bool inCircle( const float x,
                       const float y,
                       const float cx,
                       const float cy,
                       const float r ) noexcept;

        /**
        * @brief Computes the error function of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the value of the error function is
        * returned.
        */
        float erf( float x );

        /**
        * @brief Computes the complementary error function of @a x.
        * @param[in] x The floating point value
        * @return Upon successful completion, the value of the complementary error
        * function is returned.
        */
        float erfc( float x );

        /**
        * @brief Decomposes given floating point value @a x into a normalized
        * fraction and an integral power of two.
        * @param[in] x The floating point value
        * @param[in] pw2 Pointer to integer value to store the exponent to
        * @return If @a x is zero, returns zero and stores zero in @a pw2. Otherwise
        * (if @a x is not zero), returns the value @a y in the
        * range <tt>(-1;-0.5]</tt>, <tt>[0.5; 1)</tt> and stores an integer value in @a pw2 such that
        * <tt>y×2^(pw2) = x</tt> . If the value to be stored in @a pw2 is outside
        * the range of an @c int, the behavior is unspecified. If @a x is not a
        * floating-point number, the behavior is unspecified.
        */
        float rexp( float x,
                    int32_t *pw2 );

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
        float ldexp( float x,
                     int32_t pw2 );

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
        float hypot( float x,
                     float y );

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
        float nextAfter( float x,
                         float y );

        /**
        * @brief Computes the gamma function of @a x
        * @param[in] x The floating point value
        * @return Upon successful completion, this function shall return @c Gamma(x).
        * If @a x is a negative integer, a @c inf value shall be returned. If the
        * correct value would cause overflow, tgamma() shall return @c ±Inf,
        * with the same sign as the correct value of the function.
        * If @a x is @c nan, a @c nan shall be returned.
        * If @a x is @c +inf, @a x shall be returned.
        * If @a x is @c ±0, tgamma() shall return @c ±Inf.
        * If @a x is @c -inf, a @c nan  value shall be returned.
        * For IEEE Std 754-1985 float, overflow happens when <tt> 0 < x < 1/FLT_MAX</tt>,
        * and <tt>171.7 < x</tt>.
        */
        float tgamma( float x );

        /**
        * @brief Computes the natural logarithm of the absolute value of the
        * gamma function of @a x
        * @note The argument @a x need not be a non-positive integer ( is defined
        * over the reals, except the non-positive integers).
        * @param[in] x The floating point value
        * @return Upon successful completion, this function shall return the
        * logarithmic gamma of @a x. If @a x is a non-positive integer, lgamma() shall
        * return @c +inf. If the correct value would cause overflow, lgamma() shall
        * return @c ±inf having the same sign as the correct value.
        * If @a x is @c nan, a @c nan shall be returned.
        * If @a x is @c 1 or @c 2, @c +0 shall be returned.
        * If @a x is @c ±inf, @c +inf shall be returned.
        */
        float lgamma( float x );

        /**
        * @brief Return the factorial of the integer part of @a x.
        * @note The argument @a x needs to be positive
        * @warning For @a x values greater than @c 14, result is imprecise because of
        * the limited precision of the 32-bit floating point data-type.
        * With @a x values greater than @c 35, this function overflows.
        * @param[in] x The floating point value
        * @return Upon successful completion, this function shall return the
        * factorial of the integer part of @a x. If @a x is non-positive, factorial() shall
        * return @c nan. If the correct value would cause overflow, factorial() shall
        * return @c +inf.
        */
        float factorial( float x );

        /**
        * @brief Computes the associated Laguerre polynomials of the degree @a n,
        * order @a m, and argument @a x.
        * @param[in] n The degree of the polynomial, an unsigned integer value
        * @param[in] m The order of the polynomial, an unsigned integer value
        * @param[in] x The argument, a floating-point or integer value
        * @return Upon successful completion, the value of the associated Laguerre
        * polynomial of @a x shall be returned. If the argument is @c nan, a @c nan is
        * returned. If @a x is negative, @c nan is returned. If @a n or @a m is greater or
        * equal to @c 128, the behavior is implementation-defined.
        */
        float assoc_laguerre( size_t n,
                              size_t m,
                              float x );

        /**
        * @brief Computes the associated Legendre polynomials of the degree @a n,
        * order @a m, and argument @a x.
        * @param[in] n The degree of the polynomial, an unsigned integer value
        * @param[in] m The order of the polynomial, an unsigned integer value
        * @param[in] x The argument, a floating-point or integer value
        * @return Upon successful completion, the value of the associated Legendre
        * polynomial of x shall be returned. If the argument is @c nan, a @c nan is
        * returned. If <tt>|x| > 1</tt>, @c nan is returned due the domain error.
        * If @a n is greater or equal to @c 128, the behavior is implementation-defined.
        */
        float assoc_legendre( size_t n,
                              size_t m,
                              float x );

        /**
        * @brief Computes the Beta function of @a x and @a y.
        * @param[in] x The argument, a floating-point or integer value
        * @param[in] y The argument, a floating-point or integer value
        * @return Upon successful completion, the value of the Beta function of
        * @a x and @a y. If the argument is @c nan, @c nan is returned. The function
        * is only required to be defined where both @a x and @a y are greater
        * than zero, and is allowed to return @c nan otherwise.
        */
        float beta( float x,
                    float y );

        /**
        * @brief Computes the complete elliptic integral of the first kind of @a k
        * @param[in] k Elliptic modulus or eccentricity as a floating-point value
        * @return Upon successful completion, the value of the complete elliptic
        * integral of the first kind of @a k. If the argument is @c nan, @c nan is
        * returned. If <tt>|k| > 1</tt>, NaN is returned due the domain error
        */
        float comp_ellint_1( float k );

        /**
        * @brief Computes the complete elliptic integral of the second kind of @a k
        * @param[in] k Elliptic modulus or eccentricity as a floating-point value
        * @return Upon successful completion, the value of the complete elliptic
        * integral of the second kind of @a k. If the argument is @c nan, @c nan is
        * returned. If <tt>|k| > 1</tt>, @c nan is returned due the domain error
        */
        float comp_ellint_2( float k );

        /**
        * @brief Computes the complete elliptic integral of the third kind of
        * @a k and @a nu.
        * @param[in] k Elliptic modulus or eccentricity as a floating-point value
        * @param[in] nu Elliptic characteristic as a floating-point value
        * @return Upon successful completion, the value of the complete elliptic
        * integral of the third kind of @a k and @a nu. If the argument is @c nan,
        * @c nan is returned. If <tt>|k| > 1</tt>, @c nan is returned due the domain error
        */
        float comp_ellint_3( float k,
                             float nu );

        /**
        * @brief Computes the incomplete elliptic integral of the first kind of
        * @a k and @a phi
        * @param[in] k Elliptic modulus or eccentricity as a floating-point value
        * @param[in] phi Jacobi amplitude as a floating-point value given in radians
        * @return Upon successful completion, the value of the incomplete elliptic
        * integral of the first kind of @a k and @a phi. If the argument is @c nan,
        * @c nan is returned. If <tt>|k| > 1</tt>, @c nan is returned due the domain error
        */
        float ellint_1( float k,
                        float phi );

        /**
        * @brief Computes the incomplete elliptic integral of the second kind of
        * @a k and @a phi
        * @param[in] k Elliptic modulus or eccentricity as a floating-point value
        * @param[in] phi Jacobi amplitude as a floating-point value given in radians
        * @return Upon successful completion, the value of the incomplete elliptic
        * integral of the second kind of @a k and @a phi. If the argument is @c nan,
        * @c nan is returned. If <tt>|k| > 1</tt>, @c nan is returned due the domain error
        */
        float ellint_2( float k,
                        float phi );

        /**
        * @brief Computes the incomplete elliptic integral of the third kind of
        * @a k, @a nu and @a phi
        * @param[in] k Elliptic modulus or eccentricity as a floating-point value
        * @param[in] nu Elliptic characteristic as a floating-point value
        * @param[in] phi Jacobi amplitude as a floating-point value given in radians
        * @return Upon successful completion, the value of the incomplete elliptic
        * integral of the third kind of @a k, @a nu and @a phi. If the argument is @c nan,
        * @c nan is returned. If <tt>|k| > 1</tt>, @c nan is returned due the domain error
        */
        float ellint_3( float k,
                        float nu,
                        float phi );

        /**
        * @brief Computes the Exponential integral of @a num
        * @param[in] num A floating-point value
        * @return Upon successful completion, the value of the exponential integral
        * of @a num. If the argument is @c nan, @c nan is returned. If the argument
        * is @c ±0, @c -inf is returned.
        */
        float expint( float num );

        /**
        * @brief Computes the (physicist's) Hermite polynomials of the degree
        * @a n and argument @a x
        * @param[in] n The degree of the polynomial
        * @param[in] x The argument, a floating-point value
        * @return Upon successful completion, the value of order-n Hermite polynomial
        * of @a x. If the argument is @c nan, @c nan is returned. If <tt>n>=128</tt>,
        * the behavior is implementation-defined.
        */
        float hermite( size_t n,
                       float x );

        /**
        * @brief Computes the non-associated Laguerre polynomials of the degree @a n,
        * and argument @a x.
        * @param[in] n The degree of the polynomial, an unsigned integer value
        * @param[in] x The argument, a floating-point or integer value
        * @return Upon successful completion, the value of the non-associated Laguerre
        * polynomial of @a x shall be returned. If the argument is @c nan, a @c nan is
        * returned. If @a x is negative, @c nan is returned. If @a n is greater or
        * equal than @c 128, the behavior is implementation-defined.
        */
        float laguerre( size_t n,
                        float x );

        /**
        * @brief Computes the unassociated Legendre polynomials of the degree @a n,
        * and argument @a x.
        * @param[in] n The degree of the polynomial, an unsigned integer value
        * @param[in] x The argument, a floating-point or integer value
        * @return Upon successful completion, the value of the unassociated Legendre
        * polynomial of @a x shall be returned. If the argument is @c nan, a @c nan is
        * returned. The function is not required to be defined for <tt>|x| > 1 </tt>.
        * If @a n is greater or equal than @c 128, the behavior is implementation-defined.
        */
        float legendre( size_t n,
                        float x );

        /**
        * @brief Computes the Riemann zeta function of @a s
        * @param[in] s A floating-point value
        * @return Upon successful completion, the value of the Riemann zeta function
        * of @a s. If the argument is @c nan, @c nan is returned.
        */
        float riemann_zeta( float s );

        /**
        * @brief Computes the spherical Bessel function of the first kind @a n,
        * and @a x.
        * @param[in] n The order of the function
        * @param[in] x The argument to the function, a floating-point or integer value
        * @return Upon successful completion, the value of the spherical Bessel
        * function of the first kind of @a n and @a x. If the argument is @c nan, a @c nan is
        * returned. If @a n is greater or equal than @c 128, the behavior is
        * implementation-defined.
        */
        float sph_bessel( size_t n,
                          float x );

        /**
        * @brief Computes the spherical Bessel function of the second kind also
        * known as the spherical Neumann function of @a n and @a x.
        * @param[in] n The order of the function
        * @param[in] x The argument to the function, a floating-point or integer value
        * @return Upon successful completion, the value of the spherical Bessel
        * function of the second kind (spherical Neumann function) of  @a n and
        * @a x. If the argument is @c nan, a @c nan is
        * returned. If @a n is greater or equal than @c 128, the behavior is
        * implementation-defined.
        */
        float sph_neumann( size_t n,
                           float x );

        /**
        * @brief Computes the regular modified cylindrical Bessel function of
        * @a nu and @a x
        * @param[in] nu The order of the function
        * @param[in] x The argument to the function, a floating-point or integer value
        * @return Upon successful completion, the value of the regular modified
        * cylindrical Bessel function of @a nu and @a x. If the argument is
        * @c nan, a @c nan is returned. If @a nu is greater or equal than @c 128,
        * the behavior is implementation-defined.
        */
        float cyl_bessel_i( float nu,
                            float x );

        /**
        * @brief Computes the cylindrical Bessel function of the first kind of @a nu
        * and @a x
        * @param[in] nu The order of the function
        * @param[in] x The argument to the function, a floating-point or integer value
        * @return Upon successful completion, the value of the irregular modified cylindrical Bessel function
        * (also known as modified Bessel function of the second kind) of @a nu
        * and @a x. If the argument is @c nan, a @c nan is returned. If @a nu is
        * greater or equal than @c 128, the behavior is implementation-defined.
        */
        float cyl_bessel_j( float nu,
                            float x );

        /**
        * @brief Computes the irregular modified cylindrical Bessel function
        * (also known as modified Bessel function of the second kind) of @a nu
        * and @a x
        * @param[in] nu The order of the function
        * @param[in] x The argument to the function, a floating-point or integer value
        * @return Upon successful completion, the value of the irregular modified cylindrical Bessel function
        * (also known as modified Bessel function of the second kind) of @a nu
        * and @a x. If the argument is @c nan, a @c nan is returned. If @a nu is
        * greater or equal than @c 128, the behavior is implementation-defined.
        */
        float cyl_bessel_k( float nu,
                            float x );

        /**
        * @brief Computes the cylindrical Neumann function ( also known as Bessel
        * function of the second kind or Weber function) of @a nu and @a x.
        * @param[in] nu The order of the function
        * @param[in] x The argument to the function, a floating-point or integer value
        * @return Upon successful completion, the value of the cylindrical Neumann
        * function ( Bessel function of the second kind) of @a nu
        * and @a x. If the argument is @c nan, a @c nan is returned. If @a nu is
        * greater or equal than @c 128, the behavior is implementation-defined.
        */
        float cyl_neumann( float nu,
                           float x );

        /**
        * @brief 1-3) Computes the spherical associated Legendre function of
        * degree @a l, order @a m, and polar angle @a theta
        * @param[in] l The degree
        * @param[in] m The order
        * @param[in] theta Polar angle, measured in radians
        * @return Upon successful completion, the value of the spherical associated
        * Legendre function (that is, spherical harmonic with sigma = 0) of @a l,
        * @a m, and @a theta. If the argument @a theta is @c nan, a @c nan is returned.
        * If @a l is greater or equal than @c 128, the behavior is implementation-defined.
        */
        float sph_legendre( size_t l,
                            size_t m,
                            float theta );

        /*cstat -MISRAC++2008-0-1-4_b*/

        /** @brief The base of natural logarithms ( e ) given as a single-precision floating-point number*/
        constexpr float FFP_E               = ( 2.718281828459045235360287471352662498F );
        /** @brief The base 2 logarithm of e ( log_2 e ) given as a single-precision floating-point number */
        constexpr float FFP_LOG2E           = ( 1.442695040888963407359924681001892137F );
        /** @brief The base 10 logarithm of e ( log_10 e ) given as a single-precision floating-point number */
        constexpr float FFP_LOG10E          = ( 0.434294481903251827651128918916605082F );
        /** @brief The natural logarithm of 2 ( ln 2 ) given as a single-precision floating-point number */
        constexpr float FFP_LN2             = ( 0.693147180559945309417232121458176568F );
        /** @brief The natural logarithm of 10 ( ln 10 ) given as a single-precision floating-point number */
        constexpr float FFP_LN10            = ( 2.302585092994045684017991454684364208F );
        /** @brief The circumference of a circle with diameter 1, ( π ) given as a single-precision floating-point number */
        constexpr float FFP_PI              = ( 3.141592653589793238462643383279502884F );
        /** @brief Twice circumference of a circle with diameter 1, ( 2π ) given as a single-precision floating-point number */
        constexpr float FFP_2PI             = ( 6.283185307179586231995926937088370323F );
        /** @brief Half of π ( π/2 ) given as a single-precision floating-point number */
        constexpr float FFP_PI_2            = ( 1.570796326794896557998981734272092580F );
        /** @brief A quarter of π ( π/4 ) given as a single-precision floating-point number */
        constexpr float FFP_PI_4            = ( 0.785398163397448278999490867136046290F );
        /** @brief The inverse of π  ( 1/π ) given as a single-precision floating-point number */
        constexpr float FFP_1_PI            = ( 0.318309886183790671537767526745028724F );
        /** @brief Twice the inverse of π  ( 2/π ) given as a single-precision floating-point number */
        constexpr float FFP_2_PI            = ( 0.636619772367581382432888403855031356F );
        /** @brief The inverse of the square root of π ( 1/√π ) given as a single-precision floating-point number */
        constexpr float FFP_1_SQRTPI        = ( 0.564189583547756286948079451560772586F );
        /** @brief Twice the inverse of the square root of π ( 1/√π ) given as a single-precision floating-point number */
        constexpr float FFP_2_SQRTPI        = ( 1.128379167095512558560699289955664426F );
        /** @brief The square root of 2 ( √2 ) given as a single-precision floating-point number */
        constexpr float FFP_SQRT2           = ( 1.414213562373095048801688724209698079F );
        /** @brief The square root of 3 ( √3 ) given as a single-precision floating-point number */
        constexpr float FFP_SQRT3           = ( 1.732050807568877293527446341505872367F );
        /** @brief The inverse of square root of 2 ( 1/√2 ) given as a single-precision floating-point number */
        constexpr float FFP_SQRT1_2         = ( 0.707106781186547461715008466853760182F );
        /** @brief The inverse of square root of 3 ( 1/√3 ) given as a single-precision floating-point number */
        constexpr float FFP_SQRT1_3         = ( 0.577350269189625764509148780501957456F );
        /** @brief The natural logarithm of the square root of 2π given as a single-precision floating-point number */
        constexpr float FFP_LN_SQRT_2PI     = ( 0.918938533204672669540968854562379419F );
        /** @brief Constant Euler-Mascheroni given as a single-precision floating-point number*/
        constexpr float FFP_GAMMA_E         = ( 0.577215664901532860606512090082402431F );
        /** @brief The golden ratio, (1+√5)/2 given as a single-precision floating-point number*/
        constexpr float FFP_PHI             = ( 1.618033988749894848204586834365638118F );
        /** @brief Radian, value of ( 180/π ) given as a single-precision floating-point number*/
        constexpr float FFP_RADIAN          = ( 57.29577951308232286464772187173366546F );
        /** @brief The base 10 logarithm of 2 ( log_10 2 ) given as a single-precision floating-point number */
        constexpr float FFP_LOG10_2         = ( 0.301029995663981198017467022509663366F );
        /** @brief The natural logarithm of π ( ln(π) ) given as a single-precision floating-point number */
        constexpr float FFP_LN_PI           = ( 1.144729885849400163877476188645232468F );
        /*cstat +MISRAC++2008-0-1-4_b*/

        /**
        * @brief Provides several mathematical constants as single-precision floating-point numbers
        */
        struct numbers {
            /** @brief The base of natural logarithms ( e ) given as a single-precision floating-point number*/
            static constexpr float e()              { return FFP_E; }
            /** @brief The base 2 logarithm of e ( log_2 e ) given as a single-precision floating-point number */
            static constexpr float log2e()          { return FFP_LOG2E; }
            /** @brief The base 10 logarithm of e ( log_10 e ) given as a single-precision floating-point number */
            static constexpr float log10e()         { return FFP_LOG10E; }
            /** @brief The natural logarithm of 2 ( ln 2 ) given as a single-precision floating-point number */
            static constexpr float ln2()            { return FFP_LN2; }
            /** @brief The natural logarithm of 10 ( ln 10 ) given as a single-precision floating-point number */
            static constexpr float ln10()           { return FFP_LN10; }
            /** @brief The circumference of a circle with diameter 1, ( π ) given as a single-precision floating-point number */
            static constexpr float pi()             { return FFP_PI; }
            /** @brief The inverse of π  ( 1/π ) given as a single-precision floating-point number */
            static constexpr float inv_pi()         { return FFP_1_PI; }
            /** @brief The inverse of the square root of π ( 1/√π ) given as a single-precision floating-point number */
            static constexpr float inv_sqrtpi()     { return FFP_1_SQRTPI; }
            /** @brief The square root of 2 ( √2 ) given as a single-precision floating-point number */
            static constexpr float sqrt2()          { return FFP_SQRT2; }
            /** @brief The square root of 4 ( √4 ) given as a single-precision floating-point number */
            static constexpr float sqrt3()          { return FFP_SQRT3; }
            /** @brief The inverse of square root of 2 ( 1/√2 ) given as a single-precision floating-point number */
            static constexpr float inv_sqrt2()      { return FFP_SQRT1_2; }
            /** @brief The inverse of square root of 2 ( 1/√2 ) given as a single-precision floating-point number */
            static constexpr float inv_sqrt3()      { return FFP_SQRT1_3; }
            /** @brief Constant Euler-Mascheroni given as a single-precision floating-point number*/
            static constexpr float egamma()         { return FFP_GAMMA_E; }
            /** @brief The golden ratio, (1+√5)/2 given as a single-precision floating-point number*/
            static constexpr float phi()            { return FFP_PHI; }
            /** @brief Twice circumference of a circle with diameter 1, ( 2π ) given as a single-precision floating-point number */
            static constexpr float twice_pi()       { return FFP_2PI; }
            /** @brief Half of π ( π/2 ) given as a single-precision floating-point number */
            static constexpr float half_pi()        { return FFP_PI_2; }
            /** @brief A quarter of π ( π/4 ) given as a single-precision floating-point number */
            static constexpr float quarter_pi()     { return FFP_PI_4; }
            /** @brief Radian, value of ( 180/π ) given as a single-precision floating-point number*/
            static constexpr float radian()         { return FFP_RADIAN; }
            /** @brief The natural logarithm of π ( ln(π) ) given as a single-precision floating-point number */
            static constexpr float lnpi()           { return FFP_LN_PI; }
            /** @brief The natural logarithm of the square root of 2π ( ln( √2π ) ) given as a single-precision floating-point number */
            static constexpr float ln_sqrt_2pi()    { return FFP_LN_SQRT_2PI; }
        };

        /** @}*/
    }
}


#endif /*QLIBS_FFMATH*/