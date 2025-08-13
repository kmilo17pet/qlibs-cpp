/*!
 * @file rms.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Computes the RMS (Root Mean Square) of a signal using a 2-step
 * recursive average specially designed for micro-controllers with FPU.
 **/

#ifndef QLIBS_RMS
#define QLIBS_RMS

#include <include/qlibs_types.hpp>
#include <include/smoother.hpp>


/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {

    /** @addtogroup qrms RMS estimator
    * @brief Class for Recursive Root Mean Square RMS estimation
    *  @{
    */

    /**
    * @brief RMS calculator object
    */
    class rms : public smootherEXPW, public smootherMWM2, public smootherLPF1 {
        public:
            virtual ~rms() {}
            rms() = default;

            /**
            * @brief Computes the moving root mean square (RMS) of the input signal.
            * The object uses both the exponential weighting method and the sliding
            * window method to compute the moving RMS.
            * @note For a 50-60Hz sine waveform leave default parameters, use a window
            * of 8 samples and make sure you call the this method at a rate of 1mS.
            * In this configuration you will obtain an optimal estimate that converges
            * in a cycle of the wave
            * @param[in] x The raw signal.
            * @return A recursive estimation of the RMS for the incoming signal @a x.
            */
            real_t update( const real_t x ) noexcept;

            /**
            * @brief Initialize the RMS instance by setting the default optimal
            * parameters.
            * @param[in] window A pointer to the window storage, an array of
            * @a wsize elements.
            * @param[in] wsize The window size.
            * @return @c true on success, otherwise returns @c false.
            */
            bool setup( real_t * const window,
                        const size_t wsize ) noexcept;

            /**
            * @brief Initialize the RMS instance by setting the default optimal
            * parameters.
            * @param[in] win The array for the window storage
            * @return @c true on success, otherwise returns @c false.
            */
            template <size_t windowSize>
            bool setup( real_t (&win)[ windowSize ] ) noexcept
            {
                return setup( win, windowSize );
            }

            /**
            * @brief Check if the RMS instance has been initialized using setup().
            * @return @c true if the RMS instance has been initialized, otherwise
            * return @c false.
            */
            explicit operator bool() const noexcept {
                return smootherLPF1::isInitialized(); // controls truthiness
            }

            /**
            * @brief Change the recursive parameters for the moving RMS estimator.
            * @param[in] l Exponential weighting factor, specified as a positive
            * real scalar in the range [0,1]. A forgetting factor of 0.9 gives more
            * weight to the older data than does a forgetting factor of 0.1. A forgetting
            * factor of 1.0 indicates infinite memory. All the past samples are given an
            * equal weight.
            * @param[in] a A parameter to tune the 2nd stage filter. Should be
            * a value between  [ 0 < a < 1 ]. A higher value will result in a
            * smoother output but also increasing the convergence time.
            * @return @c true on success, otherwise returns @c false.
            */
            bool setParams( const real_t l,
                            const real_t a ) noexcept;
    };

    /** @}*/
}


#endif /*QLIBS_RMS*/