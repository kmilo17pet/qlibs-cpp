/*!
 * @file mathex.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Extra math and analysis functions for double precision floating point numbers.
 **/

#ifndef QLIBS_MATHEX
#define QLIBS_MATHEX

#include <include/qlibs_types.hpp>


/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /** @addtogroup  qfmathex Math extra
    * @brief Extra floating-point math and analysis functions
    * @{
    */


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
    real_t mapMinMax( const real_t x,
                      const real_t xMin,
                      const real_t xMax,
                      const real_t yMin,
                      const real_t yMax ) noexcept;
     /**
    * @brief Normalize the given input @a x in value range given by @a xMin and
    * @a xMax to value range between 0 and 1.
    * @param[in] x Input
    * @param[in] xMin Input minimum value for range
    * @param[in] xMax Input maximum  value for range
    * @return The scaled value in range [0 - 1].
    */
    real_t normalize( const real_t x,
                      const real_t xMin,
                      const real_t xMax ) noexcept;
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
    bool inRangeCoerce( real_t &x,
                        const real_t lowerL,
                        const real_t upperL ) noexcept;

    /**
    * @brief Determines if the parameters given as floating-point values are
    * approximately equal.
    * @param[in] a Input to be compared.
    * @param[in] b Input to be compared.
    * @param[in] tol Tolerance
    * @return @c true when both values are approximately equal.
    */
    bool isEqual( const real_t a,
                  const real_t b,
                  const real_t tol = 1.175494351e-38F ) noexcept;

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
    bool inPolygon( const real_t x,
                    const real_t y,
                    const real_t * const px,
                    const real_t * const py,
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
    bool inCircle( const real_t x,
                   const real_t y,
                   const real_t cx,
                   const real_t cy,
                   const real_t r ) noexcept;

    /** @}*/
}


#endif /*QLIBS_MATHEX*/