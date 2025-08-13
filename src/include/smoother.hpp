/*!
 * @file smoother.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief API to smooth noisy signals.
 **/

#ifndef QLIBS_SMOOTHER
#define QLIBS_SMOOTHER

#include <include/qlibs_types.hpp>
#include <include/tdl.hpp>

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /** @addtogroup qssmoother Smoothing filters
    * @brief Filters to smooth noisy signals
    *  @{
    */

    /**
    * @brief The smoother base abstract class
    */
    class smoother {
        protected:
            /*! @cond  */
            bool init{ true };
            bool isSetup{ false };
            static void windowSet( real_t *w,
                                   const size_t wsize,
                                   const real_t x );
            /*! @endcond  */
        public:
            virtual ~smoother() {}

            /**
            * @brief Check if the smoother filter has been initialized using setup().
            * @return @c true if the smoother has been initialized, otherwise
            * return @c false.
            */
            bool isInitialized( void ) const
            {
                return isSetup;
            }

            /**
            * @brief Check if the smoother filter has been initialized using setup().
            * @return @c true if the smoother has been initialized, otherwise
            * return @c false.
            */
            explicit operator bool() const noexcept
            {
                return isInitialized();
            }

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            virtual real_t smooth( const real_t x ) = 0;

            /**
            * @brief Reset the the smoother filter.
            * @return @c true on success, otherwise return false.
            */
            inline bool reset( void )
            {
                init = true;
                return true;
            }
    };

    /**
    * @brief A 1st order Low-Pass Filter
    */
    class smootherLPF1 : public smoother {
        protected:
            /*! @cond  */
            real_t alpha{ 0.9_re };
            real_t y1{ 0.0_re };
            /*! @endcond  */
        public:
            virtual ~smootherLPF1() {}

            /**
            * @brief Setup an initialize the 1st order Low-Pass Filter.
            * @param[in] a The filter adjustment parameter.
            * A value between  [ 0 < @a alpha < 1 ]
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const real_t a = 0.9_re );

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief A 2nd order Low-Pass Filter
    */
    class smootherLPF2 : public smoother {
        protected:
            /*! @cond  */
            real_t y1{ 0.0_re };
            real_t y2{ 0.0_re };
            real_t x1{ 0.0_re };
            real_t x2{ 0.0_re };
            real_t k{ 0.0_re };
            real_t a1{ 0.0_re };
            real_t a2{ 0.0_re };
            real_t b1{ 0.0_re };
            /*! @endcond  */
        public:
            virtual ~smootherLPF2() {}

            /**
            * @brief Setup an initialize the 2nd order Low-Pass Filter.
            * @param[in] a The filter adjustment parameter.
            * A value between  [ 0 < @a a < 1 ]
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const real_t a = 0.9_re );

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief A Moving Window Median filter
    * @note Time complexity is O(n)
    */
    class smootherMWM1 : public smoother, private nonCopyable {
        protected:
            /*! @cond  */
            real_t *w{ nullptr };
            size_t wsize{ 0U };
            /*! @endcond  */
        public:
            virtual ~smootherMWM1() {}

            /**
            * @brief Setup an initialize the Moving Window Median filter
            * @param[in] window An array to hold the samples of the window
            * @param[in] w_size The number of elements in @a window
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( real_t *window,
                        const size_t w_size );

            /**
            * @brief Setup an initialize the Moving Window Median filter
            * @param[in] window An array to hold the samples of the window
            * @return @c true on success, otherwise return @c false.
            */
            template <size_t windowSize>
            bool setup( real_t (&window)[ windowSize ] )
            {
                return setup( window, windowSize );
            }

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief A Moving Window Median filter
    * @note Time complexity is O(1)
    */
    class smootherMWM2 : public smoother, public tdl {
        protected:
            /*! @cond  */
            real_t sum{ 0.0_re };
            /*! @endcond  */
        public:
            virtual ~smootherMWM2() {}

            /**
            * @brief Setup an initialize the Moving Window Median filter
            * @param[in] window An array to hold the samples of the window
            * @param[in] w_size The number of elements in @a window
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( real_t *window,
                        const size_t w_size );

            /**
            * @brief Setup an initialize the Moving Window Median filter
            * @param[in] window An array to hold the samples of the window
            * @return @c true on success, otherwise return @c false.
            */
            template <size_t windowSize>
            bool setup( real_t (&window)[ windowSize ] )
            {
                return setup( window, windowSize );
            }

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief A Moving Outlier Removal filter
    * @note Time complexity is O(n)
    */
    class smootherMOR1 : public smoother, private nonCopyable {
        protected:
            /*! @cond  */
            real_t *w{ nullptr };
            real_t m{ 0.0_re };
            real_t alpha{ 0.8_re };
            size_t wsize{ 0U };
            /*! @endcond  */
        public:
            virtual ~smootherMOR1() {}

            /**
            * @brief Setup an initialize the Moving Outlier Removal filter
            * @param[in] window An array to hold the samples of the window
            * @param[in] w_size The number of elements in @a window
            * @param[in] a A value to adjust the filter behavior
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( real_t *window,
                        const size_t w_size,
                        const real_t a = 0.9_re );

            /**
            * @brief Setup an initialize the Moving Outlier Removal filter
            * @param[in] window An array to hold the samples of the window
            * @param[in] a A value to adjust the filter behavior
            * @return @c true on success, otherwise return @c false.
            */
            template <size_t windowSize>
            bool setup( real_t (&window)[ windowSize ],
                        const real_t a = 0.9_re )
            {
                return setup( window, windowSize, a );
            }

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief A Moving Outlier Removal filter
    * @note Time complexity is O(1)
    */
    class smootherMOR2 : public smoother, public tdl {
        protected:
            /*! @cond  */
            real_t sum{ 0.0_re };
            real_t m{ 0.0_re };
            real_t alpha{ 0.8_re };
            /*! @endcond  */
        public:
            virtual ~smootherMOR2() {}

            /**
            * @brief Setup an initialize the Moving Outlier Removal filter
            * @param[in] window An array to hold the samples of the window
            * @param[in] w_size The number of elements in @a window
            * @param[in] a A value to adjust the filter behavior
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( real_t *window,
                        const size_t w_size,
                        const real_t a = 0.9_re );

            /**
            * @brief Setup an initialize the Moving Outlier Removal filter
            * @param[in] window An array to hold the samples of the window
            * @param[in] a A value to adjust the filter behavior
            * @return @c true on success, otherwise return @c false.
            */
            template <size_t windowSize>
            bool setup( real_t (&window)[ windowSize ],
                        const real_t a = 0.9_re )
            {
                return setup( window, windowSize, a );
            }

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief A Gaussian filter
    */
    class smootherGMWF : public smoother, private nonCopyable {
        protected:
            /*! @cond  */
            real_t *w{ nullptr };
            real_t *k{ nullptr };
            size_t wsize;
            /*! @endcond  */
        public:
            virtual ~smootherGMWF() {}

            /**
            * @brief Setup an initialize the Gaussian filter
            * @param[in] sg Standard deviation [ @a sigma > 0 ]
            * @param[in] c Offset of the gaussian center [ 0 < @a offset < ( @a wsize - 1 ) ]
            * @param[in] window An array to hold the samples of the window
            * @param[in] kernel An array to hold the gaussian kernel coefficients.
            * @param[in] wk_size The number of elements in @a window or @a kernel.
            * @note Size of both @a window and @a kernel should be equal
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const real_t sg,
                        const real_t c,
                        real_t *window,
                        real_t *kernel,
                        const size_t wk_size );

            /**
            * @brief Setup an initialize the Gaussian filter
            * @param[in] sg Standard deviation [ @a sigma > 0 ]
            * @param[in] c Offset of the gaussian center [ 0 < @a offset < ( @a wsize - 1 ) ]
            * @param[in] window An array to hold the samples of the window
            * @param[in] kernel An array to hold the gaussian kernel coefficients.
            * @note Size of both @a window and @a kernel should be equal
            * @return @c true on success, otherwise return @c false.
            */
            template <size_t windowAndKernelSize>
            bool setup( const real_t sg,
                        const real_t c,
                        real_t (&window)[ windowAndKernelSize ],
                        real_t (&kernel)[ windowAndKernelSize ] )
            {
                return setup( sg, c, window, kernel, windowAndKernelSize );
            }

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief An Exponential weighting filter
    */
    class smootherEXPW : public smoother {
        protected:
            /*! @cond  */
            real_t lambda{ 0.8_re };
            real_t m{ 0_re };
            real_t w{ 1.0_re };
            /*! @endcond  */
        public:
            virtual ~smootherEXPW() {}

            /**
            * @brief Setup an initialize the Exponential weighting filter
            * @param[in] lam Forgetting factor, a value between [ 0 < @a lam < 1 ]
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const real_t lam = 0.8_re );

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief A scalar Kalman filter
    */
    class smootherKLMN : public smoother {
        protected:
            /*! @cond  */
            real_t xS{ 0.0_re };  /* state */
            real_t A{ 1.0_re };  /* x(n)=A*x(n-1)+u(n),u(n)~N(0,q) */
            real_t H{ 1.0_re };  /* z(n)=H*x(n)+w(n),w(n)~N(0,r) */
            real_t q{ 100.0_re };  /* process(predict) noise covariance */
            real_t r{ 0.9_re };  /* measure noise covariance */
            real_t p{ 100.0_re };  /* estimated error covariance */
            real_t gain;
            /*! @endcond  */
        public:
            virtual ~smootherKLMN() {}

            /**
            * @brief Setup an initialize the scalar Kalman filter
            * @param[in] processNoiseCov Process(predict) noise covariance
            * @param[in] measureNoiseCov Measure noise covariance
            * @param[in] estErrorCov Initial estimated error covariance
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const real_t processNoiseCov,
                        const real_t measureNoiseCov,
                        const real_t estErrorCov );

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief Double exponential smoothing (Holt’s Method)
    */
    class smootherDESF : public smoother {
        protected:
            /*! @cond  */
            real_t alpha;
            real_t beta;
            real_t lt;
            real_t bt;
            real_t n;
            /*! @endcond  */
        public:
            virtual ~smootherDESF() {}

            /**
            * @brief Setup an initialize the Double exponential smoothing instance
            * @param[in] a Weight for the level [ 0 < @a a < 1 ]
            * @param[in] b Weight for the trend [ 0 < @a b < 1 ]
            * @param[in] nS Number of steps for the forecast
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const real_t a,
                        const real_t b,
                        const size_t nS );

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /**
    * @brief Adaptive Filter LMS
    */
    class smootherALNF : public smoother, private nonCopyable {
        protected:
            /*! @cond  */
            real_t alpha;
            real_t mu;
            real_t *w{ nullptr };
            real_t *w_1{ nullptr };
            real_t *xx{ nullptr };
            size_t n{ 0U };
            /*! @endcond  */
        public:
            virtual ~smootherALNF() {}

            /**
            * @brief Setup an initialize the Adaptive Filter LMS
            * @param[in] a Learning rate [ 0 < @a a < 1 ]
            * @param[in] m Momentum [ 0 < @a m < 1 ]
            * @param[in] wsize The size of the @a window array
            * @param[in] window An array to store the window samples.
            * @param[in] weights An array to store the filter weights
            * @param[in] w1 An array to keep previous estimated wights. To ignore
            * pass @c nullptr as argument.
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const real_t a,
                        const real_t m,
                        const size_t wsize,
                        real_t *window,
                        real_t *weights,
                        real_t *w1 = nullptr );

            /**
            * @brief Setup an initialize the Adaptive Filter LMS
            * @param[in] a Learning rate [ 0 < @a a < 1 ]
            * @param[in] m Momentum [ 0 < @a m < 1 ]
            * @param[in] window An array to store the window samples.
            * @param[in] weights An array to store the filter weights
            * @param[in] w1 An array to keep previous estimated wights.
            * @return @c true on success, otherwise return @c false.
            */
            template <size_t windowSize>
            bool setup( const real_t a,
                        const real_t m,
                        real_t (&window)[ windowSize ],
                        real_t (&weights)[ windowSize ],
                        real_t (&w1)[ windowSize ] )
            {
                return setup( a, m, windowSize, window, weights, w1 );
            }

            /**
            * @brief Setup an initialize the Adaptive Filter LMS
            * @param[in] a Learning rate [ 0 < @a a < 1 ]
            * @param[in] window An array to store the window samples.
            * @param[in] weights An array to store the filter weights.
            * @return @c true on success, otherwise return @c false.
            */
            template <size_t windowSize>
            bool setup( const real_t a,
                        real_t (&window)[ windowSize ],
                        real_t (&weights)[ windowSize ] )
            {
                return setup( a, 0.0F, windowSize, window, weights, nullptr);
            }

            /**
            * @brief Perform the smooth operation recursively for the input signal @a x.
            * @pre Instance must be previously initialized
            * @param[in] x A sample of the input signal.
            * @return The smoothed output.
            */
            real_t smooth( const real_t x ) override;
    };

    /** @}*/
}


#endif /*QLIBS_SMOOTHER*/