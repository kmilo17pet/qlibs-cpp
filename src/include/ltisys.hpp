/*!
 * @file ltisys.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief API to evaluation continuous and discrete LTI systems.
 **/

#ifndef QLIBS_LTISYS
#define QLIBS_LTISYS

#include <include/qlibs_types.hpp>
#include <include/tdl.hpp>
#include <include/numa.hpp>


/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /** @addtogroup  qltisys LTI Systems
    * @brief Classes for recursive evaluation of LTI systems
    * defined by transfer functions
    * @{
    */

    /** @cond **/
    /*only for continuousSystem*/
    #if !defined( LTISYS_EVAL_MODEL_CONTROLLABLE ) || !defined( LTISYS_EVAL_MODEL_OBSERVABLE )
        #define LTISYS_EVAL_MODEL_CONTROLLABLE
    #endif
    /** @endcond **/

    /**
    * @brief All the possible natures of a LTI system.
    */
    enum ltisysType {
        LTISYS_TYPE_UNKNOWN = 0,    /*!< Unknown type of system */
        LTISYS_TYPE_CONTINUOUS,     /*!< Continuous-time system*/
        LTISYS_TYPE_DISCRETE,       /*!< Discrete-time controller */
    };

    /**
    * @brief Type to specify continuous states
    */
    template<size_t order>
    using continuousStates = nState[ order ];

    /**
    * @brief Type to specify discrete states
    */
    template<size_t order>
    using discreteStates = real_t[ order ];

    /**
    * @brief Continuous transfer function for easy LTI system definition
    * @note Initial conditions are zero by default.
    * @tparam order The continuous system order
    */
    template<size_t order>
    struct continuousTF {
        /** @cond **/
        real_t num[ order+1 ];
        real_t den[ order+1 ];
        continuousStates<order> states = {};
        /** @endcond **/

        /**
        * @brief Constructor for the continuousTF class
        * @param[in] numerator An array of size @c order+1 with the coefficients
        * of the numerator in descending powers of s.
        * @param[in] denominator An array of size @c order+1 with the coefficients
        * of the denominator in descending powers of s.
        */
        continuousTF( const real_t ( &numerator )[ order + 1 ], const real_t ( &denominator )[ order + 1 ] )
        {
            static_assert( order >= 1 , "Order should be greater than 0" );
            for ( size_t i = 0; i <= order; ++i ) {
                num[ i ] = numerator[ i ];
                den[ i ] = denominator[ i ];
            }
        }
    };


    /**
     * @brief Represents a time delay value for use in transportDelay constructor.
     *
     * This utility provides a convenient way to convert a continuous-time delay
     * (in seconds) into a discrete-time delay step count based on a given time step `dt`.
     *
     * Example usage:
     * @code
     * transportDelay< 2.5_td(dt) > processDelay;
     * @endcode
     */
    struct timeDelay {
        /// Delay duration in seconds.
        real_t value;

        /**
        * @brief Construct a new timeDelay object.
        * @param v The delay value in seconds.
        */
        constexpr explicit timeDelay(real_t v) : value(v) {}

        /**
        * @brief Computes the number of discrete steps equivalent to the delay.
        * @param dt The time step used in the simulation.
        * @return The delay expressed in number of steps, rounded to the nearest integer.
        */
        constexpr size_t operator()(const real_t dt) const {
            return static_cast<size_t>( ( value/dt )  + 0.5_re);
        }

        /**
        * @brief Alternate syntax to compute delay in steps using indexing operator.
        * @param dt The time step used in the simulation.
        * @return The delay expressed in number of steps, rounded to the nearest integer.
        */
        constexpr size_t operator[](const real_t dt) const {
            return static_cast<size_t>( ( value/dt ) + 0.5_re);
        }
    };

    /**
    * @brief Literal for creating a timeDelay from a floating-point value.
    *
    * Example:
    * @code
    * auto d = 0.2_td; // same as timeDelay(0.2_re)
    * @endcode
    *
    * @param v The delay value in seconds.
    * @return A timeDelay instance.
    */
    constexpr timeDelay operator"" _td(long double v) {
        return timeDelay(static_cast<real_t>(v));
    }

    /**
    * @brief Computes the delay in discrete steps using the comma operator.
    *
    * This allows concise syntax like:
    * @code
    * constexpr real_t dt = 0.01_re;
    * size_t steps = 0.2_td, dt; // same as timeDelay(0.2_re)(dt)
    * @endcode
    *
    * @param td A timeDelay object.
    * @param dt The time step.
    * @return The delay in steps.
    */
    constexpr size_t operator,(const timeDelay td, const real_t dt) {
        return static_cast<size_t>( ( td.value/dt ) + 0.5_re );
    }

    /**
    * @brief Computes the number of discrete delays required for a specified
    * amount of time using a defined time-step.
    * @see transportDelay
    * @param[in] Time The amount of time to delay
    * @param[in] dt The time step
    * @return The number of discrete delays required to delay @a Time seconds
    * using the time step @a dt
    */
    constexpr size_t delayFromTime( const real_t Time, const real_t dt )
    {
        return static_cast<size_t>( ( Time/dt ) + 0.5_re );
    }

    /**
    * @brief Computes the number of discrete delays required for a specified
    * amount of time using a defined time-step.
    * @see transportDelay
    * @param[in] Time The amount of time to delay
    * @param[in] dt The time step
    * @return The number of discrete delays required to delay @a Time seconds
    * using the time step @a dt
    */
    constexpr size_t delayFromTime( const timeDelay Time, const real_t dt )
    {
        return static_cast<size_t>( ( Time.value/dt ) + 0.5_re );
    }

    /** @cond **/
    class ITransportDelay {
    public:
        virtual ~ITransportDelay() = default;
        virtual real_t delay(const real_t xInput) noexcept = 0;
        virtual real_t operator()(const real_t xInput) noexcept = 0;
        virtual size_t getNumberOfDelays() const noexcept = 0;
    };
    /** @endcond **/

    /**
    * @brief Delays the input by a specified amount of time. You can use this
    * class to simulate a time delay.
    * @see delayFromTime
    * @tparam numberOfDelay The number of discrete delays to be used. Use the
    * timeDelay facility or  delayFromTime() function to determine the number
    * of discrete delays required for a specified amount of time.
    * Example :
    * @code{.cpp}
    * constexpr real_t dt = 0.1_re;
    * transportDelay< 2.5_td(dt) > myDelay1;
    * transportDelay< delayFromTime(5.2, dt) > myDelay2;
    * transportDelay< 4.3_td[dt] > myDelay2;
    * @endcode
    */
    template<size_t numberOfDelays>
    class transportDelay : public ITransportDelay {
        private:
            real_t buf[ numberOfDelays + 1 ];
            tdl dl;
        public:
            virtual ~transportDelay() {}
            /**
            * @brief Constructor for the transportDelay class
            * @param[in] initValue The output generated by the block between the
            * start of the simulation and the Time delay.
            */
            transportDelay( const real_t initValue = 0.0_re )
            {
                static_assert( numberOfDelays >= 1 , "Delay taps should be greater than 0" );
                dl.setup( buf, initValue);
            }

            /**
            * @brief Delays the input by a specified amount of time.
            * @param[in] xInput The signal to be delayed.
            * @return The delayed input signal
            */
            real_t delay( const real_t xInput ) noexcept override
            {
                dl.insertSample( xInput );
                return dl.getOldest();
            }

            /**
            * @brief Delays the input by a specified amount of time.
            * @param[in] xInput The signal to be delayed.
            * @return The delayed input signal
            */
            real_t operator()( const real_t xInput ) noexcept override
            {
                return delay( xInput );
            }

            /**
            * @brief Returns the number of delay steps configured for this instance.
            *
            * @return The number of delays (equal to the template parameter).
            */
            size_t getNumberOfDelays() const noexcept {
                return numberOfDelays;
            }
    };

    /**
    * @brief Delays the input by a specified amount of samples. You can use this
    * class to simulate a discrete time delay.
    * @tparam Delay The number of samples to delay the signal.
    */
    template<size_t delay>
    using discreteDelay = transportDelay<delay>;


    /**
    * @brief Discrete transfer function for easy LTI system definition
    * @note Initial conditions are zero by default.
    * @tparam NB Order of the numerator of the transfer function
    * @tparam NA Order of the denominator of the transfer function
    */
    template<size_t NB, size_t NA>
    struct discreteTF {
        /** @cond **/
        real_t num[ NB ];
        real_t den[ NA ];
        discreteStates<(NA>NB)? NA:NB> states = {};
        /** @endcond **/

        /**
        * @brief Constructor for the discreteTF class
        * @param[in] numerator An array of size @c NB with the coefficients
        * of the numerator in descending powers of z.
        * @param[in] denominator An array of size @c NA with the coefficients
        * of the denominator in descending powers of z.
        */
        discreteTF( const real_t ( &numerator )[ NB ], const real_t ( &denominator)[ NA ] ) {
            for ( size_t i = 0; i < NB; ++i ) {
                num[ i ] = numerator[ i ];
            }
            for ( size_t i = 0; i < NA; ++i ) {
                den[ i ] = denominator[ i ];
            }
        }
    };


    /**
    * @brief A LTI system base class
    */
    class ltisys : public tdl {
        protected:
            /** @cond **/
            real_t *a{ nullptr };
            real_t *b{ nullptr };
            size_t n{ 0U };
            size_t na{ 0U };
            size_t nb{ 0U };
            real_t b0{ 0.0_re };
            real_t min{ -REAL_MAX };
            real_t max{ +REAL_MAX };
            void normalizeTransferFunction( real_t *num,
                                            real_t *den,
                                            size_t n_num,
                                            size_t n_den );
            real_t saturate( real_t y );
            ltisysType type{ LTISYS_TYPE_UNKNOWN };
            virtual real_t update( const real_t u ) = 0;
            /** @endcond **/
        public:
            virtual ~ltisys() {}
            ltisys() = default;

            /**
            * @brief Drives the LTI system recursively using the provided input
            * sample.
            * @details This function evaluates the system response based on the
            * given input signal and the internal state of the system. It must
            * be called periodically with a fixed time step to maintain correct
            * system behavior.
            *
            * @pre The instance must be properly initialized before calling this method.
            *
            * @note The user is responsible for ensuring that this function is invoked at consistent
            * intervals equal to the system's time step @a dt. This timing should
            * be enforced using a hardware timer, software timer, real-time task,
            * or another reliable timing mechanism.
            *
            * @param[in] u A new input sample that excites (drives) the system.
            * @return The system's output (response) at the current time step.
            */
            real_t excite( real_t u ) noexcept;

            /**
            * @brief Drives the LTI system recursively using the provided input
            * sample.
            * @details This function evaluates the system response based on the
            * given input signal and the internal state of the system. It must
            * be called periodically with a fixed time step to maintain correct
            * system behavior.
            *
            * @pre The instance must be properly initialized before calling this method.
            *
            * @note The user is responsible for ensuring that this function is invoked at consistent
            * intervals equal to the system's time step @a dt. This timing should
            * be enforced using a hardware timer, software timer, real-time task,
            * or another reliable timing mechanism.
            *
            * @param[in] u A new input sample that excites (drives) the system.
            * @return The system's output (response) at the current time step.
            */
            real_t operator()( const real_t u ) {
                return excite( u );
            }

            /**
            * @brief Check if the LTI system is initialized.
            * @return @c true if the system has been initialized, otherwise
            * return @c false.
            */
            virtual bool isInitialized( void ) const noexcept = 0;


            /**
            * @brief Check if the LTI system is initialized.
            * @return @c true if the system has been initialized, otherwise
            * return @c false.
            */
            explicit operator bool() const noexcept {
                return isInitialized();
            }

            /**
            * @brief Set the initial states for the given system
            * @pre System should be previously initialized by using the @c setup() method
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
                           const real_t initVal = 0.0_re ) noexcept;

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
            real_t update( const real_t u ) override;
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
            bool isInitialized( void ) const noexcept override
            {
                return ( nullptr != xd ) && ( LTISYS_TYPE_DISCRETE == type );
            }

            /**
            * @brief Check if the LTI discrete system is initialized.
            * @return @c true if the system has been initialized, otherwise
            * return @c false.
            */
            explicit operator bool() const noexcept {
                return isInitialized();
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
    };

    /**
    * @brief A LTI continuous system object
    * @details The instance should be initialized using the
    * continuousSystem::setup() method.
    */
    class continuousSystem : public ltisys {
        private:
            real_t dt{ 1.0_re };
            nState *xc{ nullptr };
            real_t update( const real_t u ) override;
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
            * @param[in] nD The system order ( n ).
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
                              nState *x,
                              const size_t nD,
                              const real_t dT ) noexcept
            {
                (void)setup( num, den, x, nD, dT );
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
                              nState (&x)[ systemOrder ],
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
            * @param[in] nD The system order ( n ).
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
                        nState *x,
                        const size_t nD,
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
                        nState (&x)[ systemOrder ],
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
            bool isInitialized( void ) const noexcept override
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
    };


    using customProcessModel = real_t(*)(real_t, void*);


    /**
     * @brief A Smith Predictor implementation for compensating time delays in
     * control systems.
     * @details The Smith Predictor is a model-based feedforward control strategy
     * designed to improve performance in systems with significant dead time.
     * It estimates the delay-free output using an internal model of the process,
     * a delay block, and optionally an output filter model.
     */
    class smithPredictor {
    private:
        ltisys *model{nullptr};
        customProcessModel modelAlternate{ nullptr };
        ITransportDelay *modelDelay;
        ltisys *filter{ nullptr };
        real_t yp_hat;
        void *alternateData{ nullptr };
    public:
        virtual ~smithPredictor() {}
        /**
         * @brief Constructs a Smith Predictor with a plant model, delay model,
         * and optional initial output estimate.
         * @param[in] modelTf Reference to the LTI system model representing the
         * delay-free plant.
         * @param[in] mDelay Reference to the transport delay block modeling the
         * plant’s dead time.
         * @param[in] initialCondition Initial value for the internal output
         * prediction (@c yp_hat). Default is 0.0.
         *
         * @note Both the model and delay block are passed by reference and stored internally as pointers.
         */
        smithPredictor( ltisys& modelTf,
                        ITransportDelay& mDelay,
                        const real_t initialCondition = 0.0_re )
                        : model(&modelTf), modelDelay(&mDelay), yp_hat(initialCondition) {}

        /**
         * @brief Constructs a Smith Predictor with a plant model, delay model,
         * and optional initial output estimate.
         * @param[in] modelCustom Reference to the system model representing the
         * delay-free plant. User should define a custom function that recreates
         * the system dynamics.
         * @param[in] mDelay Reference to the transport delay block modeling the
         * plant’s dead time.
         * @param[in] initialCondition Initial value for the internal output
         * prediction (@c yp_hat). Default is 0.0.
         */
        smithPredictor( customProcessModel modelCustom,
                        ITransportDelay& mDelay,
                        const real_t initialCondition = 0.0_re )
                        : modelAlternate(modelCustom), modelDelay(&mDelay), yp_hat(initialCondition) {}

        /**
         * @brief Constructs a Smith Predictor with a plant model, delay model,
         * and optional initial output estimate.
         * @param[in] modelTf Reference to the LTI system model representing the
         * delay-free plant.
         * @param[in] mDelay Reference to the transport delay block modeling the
         * plant’s dead time.
         * @param[in] filterTf Reference to the LTI system used as the robustness
         * filter.
         * @param[in] initialCondition Initial value for the internal output
         * prediction (@c yp_hat). Default is 0.0.
         *
         * @note Both the model and delay block are passed by reference and stored internally as pointers.
         */
        smithPredictor( ltisys& modelTf,
                        ITransportDelay& mDelay,
                        ltisys& filterTf,
                        const real_t initialCondition = 0.0_re )
                        : model(&modelTf), modelDelay(&mDelay), filter(&filterTf), yp_hat(initialCondition) {}

        /**
         * @brief Constructs a Smith Predictor with a plant model, delay model,
         * and optional initial output estimate.
         * @param[in] modelCustom Reference to the system model representing the
         * delay-free plant. User should define a custom function that recreates
         * the system dynamics.
         * @param[in] mDelay Reference to the transport delay block modeling the
         * plant’s dead time.
         * @param[in] filterTf Reference to the LTI system used as the robustness
         * filter.
         * @param[in] initialCondition Initial value for the internal output
         * prediction (@c yp_hat). Default is 0.0.
         *
         * @note Both the model and delay block are passed by reference and stored internally as pointers.
         */
        smithPredictor( customProcessModel modelCustom,
                        ITransportDelay& mDelay,
                        ltisys& filterTf,
                        const real_t initialCondition = 0.0_re )
                        : modelAlternate(modelCustom), modelDelay(&mDelay), filter(&filterTf), yp_hat(initialCondition) {}

        /**
         * @brief Updates the internal prediction based on control input and
         * measured plant output.
         * @details This method should be called at each control step. It updates
         * the internal predicted output by propagating the input through the
         * internal delay-free model and adjusting the prediction using the
         * actual plant output.
         *
         * @param[in] ut The control input applied to the real plant.
         * @param[in] yt The actual measured output from the delayed plant.
         * @return @c true if the prediction was successfully updated, @c false otherwise.
         *
         * @note The time-step used to call this function must match the time
         * base used by both the plant model and the delay block.
         */
        bool updatePrediction( const real_t ut,
                               const real_t yt ) noexcept;

        /**
         * @brief Retrieves the current delay-free predicted output of the system.
         * @return The internally computed predicted output.
         */
        real_t getPrediction() const noexcept {
            return yp_hat;
        }

        /**
         * @brief Sets an optional filter for the internal Smith Predictor model.
         * @details The filter is used to improve the system's robustness, attenuate
         * measurement noise, and ensure internal stability, particularly in the case
         * of unstable plant dynamics.
         *
         * @param[in] filterTf Reference to the LTI system used as the robustness
         * filter.
         * @return @c true if the filter was set successfully.
         *
         * @note The filter is applied internally to the model-based prediction
         * and does not affect the plant or delay block directly.
         */
        bool setFilter( ltisys& filterTf ) noexcept {
            bool retValue = filter != &filterTf;
            if ( retValue ) {
                filter = &filterTf;
            }
            return retValue;
        }

        /**
         * @brief Sets the model data when alternate when a custom delay-free
         * plant model is used
         *
         * @param[in] data The data passed to the user-defined model at the moment
         * of evaluation
         * @return @c true on success. False otherwise
         */
        bool setModelData( void* data ) noexcept {
            bool retValue = data != alternateData;
            if ( retValue ) {
                alternateData = data;
            }
            return retValue;
        }
    };

    /** @}*/
}

#endif /*QLIBS_LTISYS*/