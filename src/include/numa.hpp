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
            * @return none
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
            * @return none
            */
            void init( const real_t x0 = 0.0_re,
                       const real_t sn_1 = 0.0_re,
                       const real_t sn_2 = 0.0_re ) noexcept;

            /**
            * @brief Perform a numerical integration step.
            * @param[in] s The input signal
            * @param[in] dt The time-step given in seconds.
            * @param[in] bUpdate Flag to update the states ( @c true by default).
            * @return The current value of the integration step.
            */
            real_t integrate( const real_t s,
                              const real_t dt,
                              const bool bUpdate = true ) noexcept;

             /**
            * @brief Perform a numerical derivation step by using the delta rule.
            * @param[in] s The input signal
            * @param[in] dt The time-step given in seconds.
            * @param[in] bUpdate Flag to update the states ( @c true by default).
            * @return The current value of the derivation step.
            */
            real_t derive( const real_t s,
                           const real_t dt,
                           const bool bUpdate = true ) noexcept;

            /**
            * @brief Set integration method .
            * @param[in] m The desired integration method. Use one of the following:
            *
            * @c ::INTEGRATION_RECTANGULAR : Integrate using the Rectangular rule.
            *
            * @c ::INTEGRATION_TRAPEZOIDAL : (default) Integrate using the Trapezoidal rule.
            *
            * @c ::INTEGRATION_SIMPSON : Integrate using the Simpson's 1/3 rule.
            *
            * @c ::INTEGRATION_QUADRATIC : Integrate using a parabola fit to three points.
            *
            * @return @c true on success, otherwise return @c false.
            */
            inline void setIntegrationMethod( integrationMethod m ) noexcept
            {
                iMethod = m;
            }

            /**
            * @brief Set derivation method.
            * @param[in] m The desired derivation method. Use one of the following:
            *
            * @c ::DERIVATION_2POINTS : (default) Derivative using two points.
            *
            * @c ::DERIVATION_BACKWARD : Derivative using the three-point backward-difference.
            *
            * @c ::DERIVATION_FORWARD : Derivative using the three-point forward-difference.
            *
            * @return @c true on success, otherwise return @c false.
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

    /** @}*/
}

#endif /*QLIBS_NUMA*/