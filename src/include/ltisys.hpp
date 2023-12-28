/*!
 * @file ltisys.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief API to evaluation continuous and discrete LTI systems.
 **/

#ifndef QLIBS_LTISYS
#define QLIBS_LTISYS

#include "include/types.hpp"
#include "include/tdl.hpp"
#include "include/numa.hpp"


/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /** @addtogroup  qltisys LTI Systems
    * @brief Classes for recursive evaluation of LTI systems
    * defined by transfer functions
    * @{
    */

    enum ltisysType {
        LTISYS_TYPE_UNKNOWN = 0,
        LTISYS_TYPE_CONTINUOUS,
        LTISYS_TYPE_DISCRETE,
    };

    /**
    * @brief Type to specify continuous states
    */
    template<size_t order>
    using continuousStates = state[ order ];

    /**
    * @brief Type to specify discrete states
    */
    template<size_t order>
    using discreteStates = real_t[ order ];

    /**
    * @brief Continuous transfer function for easy LTI system definition
    * @note Initial conditions are zero by default.
    */
    template<size_t order>
    struct continuousTF
    {
        real_t num[ order+1 ];
        real_t den[ order+1 ];
        continuousStates<order> states = {};
        continuousTF( const real_t ( &numerator )[ order + 1 ], const real_t ( &denominator )[ order + 1 ] ) {
            for ( size_t i = 0; i <= order; ++i ) {
                num[ i ] = numerator[ i ];
                den[ i ] = denominator[ i ];
            }
        }
    };

    /**
    * @brief Discrete transfer function for easy LTI system definition
    * @note Initial conditions are zero by default.
    */
    template<size_t NB, size_t NA>
    struct discreteTF
    {
        real_t num[ NB ];
        real_t den[ NA ];
        discreteStates<(NA>NB)? NA:NB> states = {};
        discreteTF( const real_t ( &numerator )[ NB ], const real_t ( &denominator)[ NA ] ) {
            for ( size_t i = 0; i <= NB; ++i ) {
                num[ i ] = numerator[ i ];
            }
            for ( size_t i = 0; i <= NA; ++i) {
                den[ i ] = denominator[ i ];
            }
        }
    };


    /**
    * @brief A LTI system base class
    */
    class ltisys : public tdl {
        protected:
            real_t *a{ nullptr };
            real_t *b{ nullptr };
            size_t n{ 0U };
            size_t na{ 0U };
            size_t nb{ 0U };
            real_t b0{ 0.0 };
            real_t min{ DBL_MIN };
            real_t max{ DBL_MAX };
            void normalizeTransferFunction( real_t *num,
                                            real_t *den,
                                            size_t n_num,
                                            size_t n_den );
            real_t saturate( real_t y );
            ltisysType type{ LTISYS_TYPE_UNKNOWN };
        public:
            virtual ~ltisys() {}
            ltisys() = default;

            /**
            * @brief Drives the LTI system recursively using the input signal provided
            * @pre Instance must be previously initialized.
            * @param[in] u A sample of the input signal that excites the system
            * @return The system response.
            */
            virtual real_t excite( real_t u ) = 0;

            /**
            * @brief Check if the LTI system is initialized.
            * @return @c true if the system has been initialized, otherwise 
            * return @c false.
            */
            virtual bool isInitialized( void ) const = 0;

            /**
            * @brief Set the initial states for the given system
            * @pre System should be previously initialized by using ::setup() method
            * @param[in] xi An array of n-elements with the initial state values. User
            * can pass @c nullptr as argument to set initial conditions equal to zero.
            * @return @c true on success, otherwise return @c false.
            */
            virtual bool setInitStates( const real_t *xi = nullptr ) = 0;

            /**
            * @brief Get the LTI system type.
            * @return @c true if the system has been initialized, otherwise 
            * return @c false.
            */
            ltisysType getType( void ) const
            {
                return type;
            }

            /**
            * @brief Set the input delay for LTI system.
            * @param[in] w A array of n elements with the delay window for the
            *  input channel.
            * @param[in] nD The number of elements of @a w.
            * @param[in] initVal The initial value of the input channel.
            * @return ::LTISYS_TYPE_CONTINUOUS or LTISYS_TYPE_DISCRETE. If the
            * system is not initialized, ::LTISYS_TYPE_UNKNOWN 
            */
            bool setDelay( real_t * const w,
                           const size_t nD,
                           const real_t initVal = 0.0 ) noexcept;

            /**
            * @brief Setup the output saturation for the LTI system.
            * @param[in] minV The minimal value allowed for the output.
            * @param[in] maxV The maximal value allowed for the output.
            * @return @c true on success, otherwise return @c false.
            */
            bool setSaturation( const real_t minV,
                                const real_t maxV ) noexcept;
    }; 

    /**
    * @brief A LTI discrete system object
    * @details The instance should be initialized using the
    * discreteSystem::setup() method.
    */
    class discreteSystem : public ltisys {
        private:
            real_t *xd{ nullptr };
            real_t update( const real_t u );
        public:
            virtual ~discreteSystem() {}

            /**
            * @brief Constructor for a the discrete LTI system.
            * @param[in,out] num : An array of @c nb elements with the numerator
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] den An array of or @c na elements with the 
            * denominator coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] x Initial conditions of the system. Should be an,
            * an array of type discreteStates with max(na,nb) elements
            * The supplied array will be updated on every invocation of
            * discreteSystem::excite().
            * @param[in] n_b Number of elements in @a num
            * 
            * example: \f$ b_{0}+b_{1}z^{-1}+b_{2}z^{-2}+b_{3}z^{-3}, nb = 4 \f$
            * @param[in] n_a Number of elements in @a den.
            * 
            * example 1: \f$ a_{0}+a_{1}z^{-1}+a_{2}z^{-2}+a_{3}z^{-3}, na = 4 \f$
            * @return @c true on success, otherwise return @c false.
            * @note By default, initial conditions are set to zero. To change the initial
            * conditions to the desired values, use the 
            * discreteSystem::setInitStates() method.
            */
            discreteSystem( real_t *num,
                            real_t *den,
                            real_t *x,
                            const size_t n_b,
                            const size_t n_a ) noexcept
            {
                (void)setup( num, den, x, n_b, n_a );
            }

            /**
            * @brief Constructor for a the discrete LTI system.
            * @param[in,out] num : An array of @c nb elements with the numerator
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] den An array of or @c na elements with the 
            * denominator coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] x Initial conditions of the system. Should be an,
            * an array of type discreteStates with max(na,nb) elements
            * The supplied array will be updated on every invocation of
            * discreteSystem::excite().
            * @return @c true on success, otherwise return @c false.
            * @note By default, initial conditions are set to zero. To change the initial
            * conditions to the desired values, use the 
            * discreteSystem::setInitStates() method.
            */
            template<size_t NB, size_t NA>
            discreteSystem( real_t (&num)[ NB ],
                            real_t (&den)[ NA ],
                            real_t *x ) noexcept
            {
                (void)setup( num, den, x, NB, NA );
            }

            /**
            * @brief Constructor for a the discrete LTI system from a transfer
            * function definition
            * @param[in,out] dtf : The transfer function definition
            */
            template<size_t NB, size_t NA>
            discreteSystem( discreteTF<NB, NA>& dtf ) noexcept
            {
                (void)setup( dtf.num, dtf.den, dtf.states, NB, NA );
            }

            /**
            * @brief Setup and initialize an instance of the discrete LTI system.
            * @param[in,out] num : An array of @c nb elements with the numerator
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] den An array of or @c na elements with the 
            * denominator coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] x Initial conditions of the system. Should be an,
            * an array of type discreteStates with max(na,nb) elements
            * The supplied array will be updated on every invocation of
            * discreteSystem::excite().
            * @param[in] n_b Number of elements in @a num
            * 
            * example: \f$ b_{0}+b_{1}z^{-1}+b_{2}z^{-2}+b_{3}z^{-3}, nb = 4 \f$
            * @param[in] n_a Number of elements in @a den.
            * 
            * example 1: \f$ a_{0}+a_{1}z^{-1}+a_{2}z^{-2}+a_{3}z^{-3}, na = 4 \f$
            * @return @c true on success, otherwise return @c false.
            * @note By default, initial conditions are set to zero. To change the initial
            * conditions to the desired values, use the 
            * discreteSystem::setInitStates() method.
            */
            bool setup( real_t *num,
                        real_t *den,
                        real_t *x,
                        const size_t n_b,
                        const size_t n_a ) noexcept;

            /**
            * @brief Setup and initialize an instance of the discrete LTI system.
            * @param[in,out] num : An array of @c nb elements with the numerator
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] den An array of or @c na elements with the 
            * denominator coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] x Initial conditions of the system. Should be an,
            * an array of type discreteStates with max(na,nb) elements
            * The supplied array will be updated on every invocation of
            * discreteSystem::excite().
            * @return @c true on success, otherwise return @c false.
            * @note By default, initial conditions are set to zero. To change the initial
            * conditions to the desired values, use the 
            * discreteSystem::setInitStates() method.
            */
            template<size_t NB, size_t NA>
            bool setup( real_t (&num)[ NB ],
                        real_t (&den)[ NA ],
                        real_t *x )
            {
                return setup( num, den, x, NB, NA );
            }

            /**
            * @brief Setup and initialize an instance of the discrete LTI system
            * from a transfer function definition
            * @param[in,out] dtf : The transfer function definition
            * @return @c true on success, otherwise return @c false.
            */
            template<size_t NB, size_t NA>
            bool setup( discreteTF<NB, NA>& dtf ) noexcept
            {
                return setup( dtf.num, dtf.den, dtf.states, NB, NA );
            }

            /**
            * @brief Check if the LTI discrete system is initialized.
            * @return @c true if the system has been initialized, otherwise 
            * return @c false.
            */
            bool isInitialized( void ) const override
            {
                return ( nullptr != xd ) && ( LTISYS_TYPE_DISCRETE == type );
            }

            /**
            * @brief Set the initial states for the discrete system
            * @pre System should be previously initialized by using the
            * discreteSystem::setup() method
            * @param[in] xi An array of n-elements with the initial state values. User
            * can pass @c nullptr as argument to set initial conditions equal to zero.
            * @return @c true on success, otherwise return @c false.
            */
            bool setInitStates( const real_t *xi = nullptr ) override;

            /**
            * @brief Evaluate the discrete FIR filter by updating the delay lines of
            * @a x inside the window @a w of size @a wsize with the coefficients given
            * in @a c. If @a c it's not supplied, this function just perform the window
            * update.
            * @param[in,out] w An array of @a wsize elements that holds the window with
            * the delay lines of @a x.
            * @param[in] wsize The number of elements of @a w.
            * @param[in] x A sample of the input signal.
            * @param[in] c An array of @a wsize elements with the coefficients of the
            * FIR filter. Coefficients should be given in descending powers of the
            * nth-degree polynomial. To ignore pass @c nullptr.
            * @return If @a c is provided, returns the evaluation of the FIR filter.
            * otherwise return the sum of the updated window @a w.
            */
            static real_t updateFIR( real_t *w,
                                     const size_t wsize,
                                     const real_t x,
                                     const real_t * const c = nullptr );

            /**
            * @brief Drives the discrete LTI system recursively using the input 
            * signal provided
            * @pre Instance must be previously initialized using the 
            * discreteSystem::setup() method
            * @note The user must ensure that this function is executed in the 
            * required sample time @a T either by using a HW or SW timer, a 
            * real time task or a timing service.
            * @param[in] u A sample of the input signal that excites the system
            * @return The system response.
            */
            real_t excite( real_t u ) override;

    };

    /**
    * @brief A LTI continuous system object
    * @details The instance should be initialized using the 
    * continuousSystem::setup() method.
    */
    class continuousSystem : public ltisys {
        private:
            real_t dt{ 1.0 };
            state *xc{ nullptr };
            real_t update( const real_t u );
        public:
            virtual ~continuousSystem() {}

            /**
            * @brief Constructor for an instance of a LTI continuous system.
            * @param[in,out] num : An array of n+1 elements with the numerator 
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] den An array of n+1 elements with the denominator 
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] x Initial conditions of the system. For a continuos system,
            * an array of type continuousStates with n elements.
            * The supplied array will be updated on every invocation of
            * continuousSystem::excite().
            * @param[in] N The system order ( n ).
            * 
            * example 2: \f$ \frac{ b_{0}s^{2}+b_{1}s+b_{2} }{ a_{0}s^{2} + a_{1}s + a_{2} }, na = 3 \f$
            * @note Size of @a num and @a den should be equal.
            * @param[in] dT The time-step of the continuos system.
            * @note By default, initial conditions are set to zero. To change the initial
            * conditions to the desired values, use the 
            * continuousSystem::setInitStates() method.
            */
            continuousSystem( real_t *num,
                              real_t *den,
                              state *x,
                              const size_t N,
                              const real_t dT ) noexcept
            {
                (void)setup( num, den, x, N, dT );
            }

            /**
            * @brief Constructor for an instance of a LTI continuous system
            * from a transfer function definition
            * @param[in] ctf : The continuous transfer-function
            * @param[in] dT The time-step of the continuos system.
            */
            template<size_t order>
            continuousSystem( continuousTF<order>& ctf, const real_t dT )
            {
                (void)setup( ctf.num, ctf.den, ctf.states, order, dT );
            }

            /**
            * @brief Constructor for an instance of a LTI continuous system.
            * @param[in,out] num : An array of n+1 elements with the numerator 
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] den An array of n+1 elements with the denominator 
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] x Initial conditions of the system. For a continuos system,
            * an array of type continuousStates with n elements.
            * The supplied array will be updated on every invocation of
            * continuousSystem::excite().
            * @note Size of @a num and @a den should be equal.
            * @param[in] dT The time-step of the continuos system.
            * @note By default, initial conditions are set to zero. To change the initial
            * conditions to the desired values, use the 
            * continuousSystem::setInitStates() method.
            */
            template<size_t systemOrder>
            continuousSystem( real_t (&num)[ systemOrder+1 ],
                              real_t (&den)[ systemOrder+1 ],
                              state (&x)[ systemOrder ],
                              const real_t dT ) noexcept
            {
                (void)setup( num, den, x, systemOrder, dT );
            }

            /**
            * @brief Setup and initialize an instance of a LTI continuous system.
            * @param[in,out] num : An array of n+1 elements with the numerator 
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] den An array of n+1 elements with the denominator 
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] x Initial conditions of the system. For a continuos system,
            * an array of type continuousStates with n elements.
            * The supplied array will be updated on every invocation of
            * continuousSystem::excite().
            * @param[in] N The system order ( n ).
            * 
            * example 2: \f$ \frac{ b_{0}s^{2}+b_{1}s+b_{2} }{ a_{0}s^{2} + a_{1}s + a_{2} }, na = 3 \f$
            * @note Size of @a num and @a den should be equal.
            * @param[in] dT The time-step of the continuos system.
            * @return @c true on success, otherwise return @c false.
            * @note By default, initial conditions are set to zero. To change the initial
            * conditions to the desired values, use the 
            * continuousSystem::setInitStates() method.
            */
            bool setup( real_t *num,
                        real_t *den,
                        state *x,
                        const size_t N,
                        const real_t dT ) noexcept;

            /**
            * @brief Setup and initialize an instance of a LTI continuous system.
            * @param[in,out] num : An array of n+1 elements with the numerator 
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] den An array of n+1 elements with the denominator 
            * coefficients of the transfer function.
            * Coefficients should be given in descending powers of the n or nb-degree
            * polynomial. Coefficients will be normalized internally.
            * @param[in,out] x Initial conditions of the system. For a continuos system,
            * an array of type continuousStates with n elements.
            * The supplied array will be updated on every invocation of
            * continuousSystem::excite().
            * @note Size of @a num and @a den should be equal.
            * @param[in] dT The time-step of the continuos system.
            * @return @c true on success, otherwise return @c false.
            * @note By default, initial conditions are set to zero. To change the initial
            * conditions to the desired values, use the 
            * continuousSystem::setInitStates() method.
            */
            template<size_t systemOrder>
            bool setup( real_t (&num)[ systemOrder+1 ],
                        real_t (&den)[ systemOrder+1 ],
                        state (&x)[ systemOrder ],
                        const real_t dT ) noexcept
            {
                return setup( num, den, x, systemOrder, dT );
            }

            /**
            * @brief Setup and initialize an instance of a LTI continuous system
            * from a transfer function definition
            * @param[in] ctf : The continuous transfer-function
            * @param[in] dT The time-step of the continuos system.
            * @return @c true on success, otherwise return @c false.
            */
            template<size_t order>
            bool setup( continuousTF<order>& ctf, const real_t dT )
            {
                return setup( ctf.num, ctf.den, ctf.states, order, dT );
            }

            /**
            * @brief Check if the LTI continuous system is initialized.
            * @return @c true if the system has been initialized, otherwise 
            * return @c false.
            */
            bool isInitialized( void ) const override
            {
                return ( nullptr != xc ) && ( LTISYS_TYPE_CONTINUOUS == type );
            }

            /**
            * @brief Set the initial states for the continuous system
            * @pre System should be previously initialized by using the
            * discreteSystem::setup() method
            * @param[in] xi An array of n-elements with the initial state values. User
            * can pass @c nullptr as argument to set initial conditions equal to zero.
            * @return @c true on success, otherwise return @c false.
            */
            bool setInitStates( const real_t *xi = nullptr ) override;

            /**
            * @brief Set integration method of the continuous system.
            * @param[in] m The desired integration method. Use one of the following:
            *
            * @c INTEGRATION_RECTANGULAR : Integrate using the Rectangular rule.
            *
            * @c INTEGRATION_TRAPEZOIDAL : (default) Integrate using the Trapezoidal rule.
            *
            * @c INTEGRATION_SIMPSON : Integrate using the Simpson's 1/3 rule.
            *
            * @return @c true on success, otherwise return @c false.
            */
            bool setIntegrationMethod( integrationMethod m );

            /**
            * @brief Drives the continuous LTI system recursively using the input 
            * signal provided
            * @pre Instance must be previously initialized using the 
            * continuousSystem::setup() method
            * @note The user must ensure that this function is executed in the time
            * specified in @a dt either by using a HW or SW timer, a real time task,
            * or a timing service.
            * @param[in] u A sample of the input signal that excites the system
            * @return The system response.
            */
            real_t excite( real_t u ) override;
    };

    /** @}*/
}

#endif /*QLIBS_LTISYS*/