/*!
 * @file interp1.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Helper class that implements a Tapped Delay Line (TDL). A TDL is a
 * delay line that provides access to its contents at arbitrary intermediate
 * delay length values.
 * This class runs in constant time O(1), so it becomes useful when you need to
 * work with long delayed lines.
 **/

#ifndef QLIBS_INTERP1
#define QLIBS_INTERP1

#include <include/qlibs_types.hpp>

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /** @addtogroup  qinterp1 1D Interpolation
    * @brief One-dimensional interpolation class.
    *  @{
    */

        /**
        * @brief An enum with all the available interpolation methods.
        */
    enum interp1Method {
        INTERP1_NEXT = 0,               /*!< Return the next neighbor.*/
        INTERP1_PREVIOUS,               /*!< Return the previous neighbor.*/
        INTERP1_NEAREST,                /*!< Return the nearest neighbor.*/
        INTERP1_LINEAR,                 /*!< Linear interpolation from nearest neighbors.*/
        INTERP1_SINE,                   /*!< Sine interpolation.*/
        INTERP1_CUBIC,                  /*!< Cubic interpolation.*/
        INTERP1_HERMITE,                /*!< Piecewise cubic Hermite interpolation.*/
        INTERP1_SPLINE,                 /*!< Catmull spline interpolation.*/
        INTERP1_CONSTRAINED_SPLINE,     /*!< A special kind of spline that doesn't overshoot.*/
        INTERP1_MAX,
    };

    /**
    * @brief A 1D interpolation object.
    */
    class interp1 {
        private:
            using interp1Fcn_t = real_t (*)( const real_t x,
                                             const real_t * const tx,
                                             const real_t * const ty,
                                             const size_t tableSize );
            const real_t *xData{ nullptr };
            const real_t *yData{ nullptr };
            const size_t dataSize{ 0U };
            interp1Fcn_t method{ &linear };
            static real_t slope( const real_t * const tx,
                                 const real_t * const ty,
                                 const size_t i );
            static real_t firstDerivate( const real_t * const tx,
                                         const real_t * const ty,
                                         const size_t n,
                                         const size_t i );
            static real_t leftSecondDerivate( const real_t * const tx,
                                              const real_t * const ty,
                                              const size_t n,
                                              const size_t i );
            static real_t rightSecondDerivate( const real_t * const tx,
                                               const real_t * const ty,
                                               const size_t n,
                                               const size_t i );
            static real_t next( const real_t x,
                                const real_t * const tx,
                                const real_t * const ty,
                                const size_t tableSize );
            static real_t previous( const real_t x,
                                    const real_t * const tx,
                                    const real_t * const ty,
                                    const size_t tableSize );
            static real_t nearest( const real_t x,
                                   const real_t * const tx,
                                   const real_t * const ty,
                                   const size_t tableSize );
            static real_t linear( const real_t x,
                                  const real_t * const tx,
                                  const real_t * const ty,
                                  const size_t tableSize );
            static real_t sine( const real_t x,
                                const real_t * const tx,
                                const real_t * const ty,
                                const size_t tableSize );
            static real_t cubic( const real_t x,
                                 const real_t * const tx,
                                 const real_t * const ty,
                                 const size_t tableSize );
            static real_t hermite( const real_t x,
                                   const real_t * const tx,
                                   const real_t * const ty,
                                   const size_t tableSize );
            static real_t spline( const real_t x,
                                  const real_t * const tx,
                                  const real_t * const ty,
                                  const size_t tableSize );
            static real_t cSpline( const real_t x,
                                   const real_t * const tx,
                                   const real_t * const ty,
                                   const size_t tableSize );
        public:
            virtual ~interp1() {}
            /**
            * @brief Constructor for the 1D interpolation instance.
            * @param[in] xTable An array of size @a sizeTable with the x points sorted in ascending order.
            * @param[in] yTabLe An array of size @a sizeTable with the y points.
            * @param[in] sizeTable The number of points in @a xTable @a yTable
            */
            interp1( const real_t * const xTable,
                     const real_t * const yTabLe,
                     const size_t sizeTable ) : xData( xTable ), yData( yTabLe ), dataSize( sizeTable ) {}
            
            /**
            * @brief Constructor for the 1D interpolation instance.
            * @param[in] xTable An array of size @a sizeTable with the x points sorted in ascending order.
            * @param[in] yTabLe An array of size @a sizeTable with the y points.
            */
            template <size_t sizeTable>
            interp1( real_t (&xTable)[ sizeTable ],
                     real_t (&yTable)[ sizeTable ] ) : interp1( xTable, yTable, sizeTable ) {}

            /**
            * @brief Specify the interpolation method to use.
            * @param[in] m The interpolation method.
            * @return @c true on success otherwise @c false.
            */
            bool setMethod( const interp1Method m );
            template <typename T>

            /**
            * @brief Interpolate input point @a x to determine the value of @a y
            * at the points @a xi using the current method. If value if beyond
            * the endpoints, extrapolation is performed using the current method.
            * @return @c The interpolated-extrapolated @a y value.
            */
            inline real_t get( const T x )
            {
                return method( static_cast<real_t>( x ), xData, yData, dataSize );
            }
    };

    /** @}*/
}


#endif /*QLIBS_INTERP1*/