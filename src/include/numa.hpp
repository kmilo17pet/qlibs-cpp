/*!
 * @file numa.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Numerical approximations for real-time signals.
 **/

#ifndef QLIBS_NUMA
#define QLIBS_NUMA

#include <include/qlibs_types.hpp>

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {

    /** @addtogroup qnuma Numerical state class
    * @brief Class to compute in real-time the numerical approximations of
    * integral and derivative operations for data values sampled periodically.
    *  @{
    */

    /**
    * @brief An enum with all the available integration methods.
    */
    enum integrationMethod {
        INTEGRATION_RECTANGULAR,    /*!< Numerical integration using the rectangular rule*/
        INTEGRATION_TRAPEZOIDAL,    /*!< Numerical integration using the trapezoidal rectangular rule*/
        INTEGRATION_SIMPSON,        /*!< Numerical integration using the Simpson's 1/3 rule.*/
        INTEGRATION_QUADRATIC,      /*!< Numerical integration using a parabola fit to three points.*/
    };

    /**
    * @brief An enum with all the available derivation methods.
    */
    enum derivationMethod {
        DERIVATION_2POINTS,         /*!< Numerical derivation using two points*/
        DERIVATION_BACKWARD,        /*!< Numerical derivation using the three-point backward difference.*/
        DERIVATION_FORWARD,         /*!< Numerical derivation using the three-point forward difference.*/
    };

    /**
    * @brief A numerical state object.
    * @details A numerical state object can be used to compute in real-time the
    * numerical approximations of integral and derivative operations for data
    * values sampled periodically.
    */
    class nState {
        private:
            real_t x[ 3 ] = { 0.0_re, 0.0_re, 0.0_re };
            void inline update( const real_t &s )
            {
                x[ 2 ] = x[ 1 ];
                x[ 1 ] = s;
            }
            integrationMethod iMethod{ INTEGRATION_TRAPEZOIDAL };
            derivationMethod dMethod{ DERIVATION_2POINTS };
        public:
            virtual ~nState() {}

            /**
            * @brief Constructor for the state object
            * @param[in] x0 initial condition at time t(0)
            * @param[in] sn_1 initial condition at time (t-1)
            * @param[in] sn_2 initial condition at time (t-2)
            */
            nState( const real_t x0 = 0.0_re,
                    const real_t sn_1 = 0.0_re,
                    const real_t sn_2 = 0.0_re ) noexcept
            {
                init( x0, sn_1, sn_2 );
            }

            /**
            * @brief Initialize the state object
            * @param[in] x0 initial condition at time t(0)
            * @param[in] sn_1 initial condition at time (t-1)
            * @param[in] sn_2 initial condition at time (t-2)
            */
            void init( const real_t x0 = 0.0_re,
                       const real_t sn_1 = 0.0_re,
                       const real_t sn_2 = 0.0_re ) noexcept;

            /**
            * @brief Perform a numerical integration step by using the specified
            * integration method.
            * @param[in] s The input signal
            * @param[in] dt The time-step given in seconds.
            * @param[in] bUpdate Flag to update the states ( @c true by default).
            * @return The current value of the integration step.
            */
            real_t integrate( const real_t s,
                              const real_t dt,
                              const bool bUpdate = true ) noexcept;

             /**
            * @brief Perform a numerical derivation step by using the specified
            * derivation method.
            * @param[in] s The input signal
            * @param[in] dt The time-step given in seconds.
            * @param[in] bUpdate Flag to update the states ( @c true by default).
            * @return The current value of the derivation step.
            */
            real_t derive( const real_t s,
                           const real_t dt,
                           const bool bUpdate = true ) noexcept;

            /**
            * @brief Sets the numerical integration method.
            * @details Configures the method used to approximate the integral of
            * input values. This allows the user to select from several common
            * numerical integration techniques.
            * @param[in] m The desired integration method. Supported options:
            *
            * @c ::INTEGRATION_RECTANGULAR : Integrate using the Rectangular rule.
            *
            * @c ::INTEGRATION_TRAPEZOIDAL : (default) Integrate using the Trapezoidal rule.
            *
            * @c ::INTEGRATION_SIMPSON : Integrate using the Simpson's 1/3 rule.
            *
            * @c ::INTEGRATION_QUADRATIC : Integrate using a parabola fit to three points.
            *
            * @note The effectiveness and accuracy of each method depend on the
            * signal characteristics and time step.
            */
            inline void setIntegrationMethod( integrationMethod m ) noexcept
            {
                iMethod = m;
            }

            /**
            * @brief Sets the numerical derivation method.
            * @details Configures the method used to compute the numerical
            * derivative of input values. Different methods offer varying accuracy
            * and responsiveness depending on signal characteristics.
            *
            * @param[in] m The desired derivation method. Supported options:
            *
            * @c ::DERIVATION_2POINTS : (default) Uses a simple two-point
            * (first-order backward) difference.
            *
            * @c ::DERIVATION_BACKWARD : Derivative using the three-point backward-difference.
            *
            * @c ::DERIVATION_FORWARD : Derivative using the three-point forward-difference.
            *
            * @note Choose the method based on the required balance between
            * accuracy and latency.
            */
            inline void setDerivationMethod( derivationMethod m ) noexcept
            {
                dMethod = m;
            }

             /**
            * @brief Get the value of the state.
            * @return The current value of the state.
            */
            real_t operator()( void ) const noexcept
            {
                return x[ 0 ];
            }

            /*! @cond  */
            real_t operator*( real_t rValue) const noexcept
            {
                return x[ 0 ]*rValue;
            }
            real_t operator/( real_t rValue) const noexcept
            {
                return x[ 0 ]/rValue;
            }
            real_t operator+( real_t rValue) const noexcept
            {
                return x[ 0 ] + rValue;
            }
            real_t operator-( real_t rValue) const noexcept
            {
                return x[ 0 ] - rValue;
            }
            friend real_t operator*( real_t rValue,
                                     const nState& s ) noexcept;
            friend real_t operator/( real_t rValue,
                                     const nState& s ) noexcept;
            friend real_t operator+( real_t rValue,
                                     const nState& s ) noexcept;
            friend real_t operator-( real_t rValue,
                                     const nState& s ) noexcept;
            /*! @endcond  */
    };

    /*! @cond  */
    inline real_t operator*( real_t rValue,
                             const nState& s ) noexcept
    {
        return rValue*s.x[ 0 ];
    }
    inline real_t operator/( real_t rValue,
                             const nState& s ) noexcept
    {
        return rValue/s.x[ 0 ];
    }
    inline real_t operator+( real_t rValue,
                             const nState& s ) noexcept
    {
        return rValue + s.x[ 0 ];
    }
    inline real_t operator-( real_t rValue,
                             const nState& s ) noexcept
    {
        return rValue - s.x[ 0 ];
    }
    /*! @endcond  */

    /**
    * @brief A numerical integration class.
    * @details A numerical integration class that can be used to compute
    * in real-time the numerical approximation of an integral for data values
    * sampled periodically. It supports optional output saturation limits.
    */
    class integrator : public nState, private nonCopyable {
        private:
            real_t dt;
            real_t min{ -REAL_MAX };
            real_t max{ +REAL_MAX };
        public:
            virtual ~integrator() {}

            /**
            * @brief Constructs an integrator block with a given @a timeStep time
            * and optional initial condition.
            * @param[in] timeStep The fixed time step (dt) used to compute the
            * integration.
            * @param[in] initialCondition The initial output value of the
            * integrator. Default is 0.0.
            * @note It is assumed that input samples will be provided at regular
            * intervals of timeStep.
            */
            integrator( const real_t timeStep,
                        const real_t initialCondition = 0.0_re );

            /**
            * @brief Sets the saturation limits for the integrator output.
            * @param[in] minV The minimum value the output can reach.
            * @param[in] maxV The maximum value the output can reach.
            * @return @c true if the limits are valid and applied; @c false
            * otherwise (e.g., minV > maxV).
            * @note If not set, the output is unbounded.
            */
            bool setSaturation( const real_t minV, const real_t maxV ) noexcept;

            /**
            * @brief Performs one step of numerical integration.
            * @param[in] xDot The input value to be integrated
            * @return The integrated value (i.e., the output of the integrator)
            * after applying saturation.
            * @note This should be called at intervals equal to the time step provided in the constructor.
            */
            real_t operator()( const real_t xDot );
    };

    /**
    * @brief A numerical derivative class.
    * @details A numerical derivative class that can be used to compute
    * in real-time the numerical approximations of derivatives for data values
    * sampled periodically.
    */
    class derivative : public nState, private nonCopyable {
        private:
            real_t dt;
        public:
            virtual ~derivative() {}

            /**
            * @brief Constructs a derivative block with a given @a timeStep and
            * optional initial condition.
            * @param[in] timeStep The fixed time step (dt) used to compute the
            * derivative.
            * @param[in] initialCondition The initial input value. Default is @c 0.0
            * @note It is assumed that input samples will be provided at regular
            * intervals of timeStep.
            */
            derivative( const real_t timeStep,
                        const real_t initialCondition = 0.0_re );

            /**
            * @brief Computes the numerical derivative based on the current
            * input value.
            * @param[in] xt The current input value.
            * @return The estimated derivative of the input signal.
            * @note This should be called at intervals equal to the time step
            * provided in the constructor.
            */
            real_t operator()( const real_t xt );
    };

    /** @}*/
}

#endif /*QLIBS_NUMA*/