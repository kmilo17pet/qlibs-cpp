/*!
 * @file numa.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Numerical approximations for real-time signals.
 **/

#ifndef QLIBS_NUMA
#define QLIBS_NUMA

#include "include/types.hpp" 

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {

    enum integrationMethod {
        INTEGRATION_RECTANGULAR,    /*!< Numerical integration using the rectangular rule*/
        INTEGRATION_TRAPEZOIDAL,    /*!< Numerical integration using the trapezoidal rectangular rule*/
        INTEGRATION_SIMPSON,        /*!< Numerical integration using the Simpson's 1/3 rule.*/
    };

    class state {
        private:
            real_t x[ 3 ] = { 0.0, 0.0, 0.0 };
            void inline update( const real_t &s )
            {
                x[ 2 ] = x[ 1 ];
                x[ 1 ] = s;
            }
            integrationMethod intMethod{ INTEGRATION_TRAPEZOIDAL };
        public:
            virtual ~state() {}

            /**
            * @brief Constructor for the state object
            * @param[in] x0 initial condition at time t(0)
            * @param[in] sn_1 initial condition at time (t-1)
            * @param[in] sn_2 initial condition at time (t-2)
            * @return none
            */
            state( const real_t x0 = 0.0,
                   const real_t sn_1 = 0.0,
                   const real_t sn_2 = 0.0 ) noexcept
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
            void init( const real_t x0 = 0.0,
                       const real_t sn_1 = 0.0,
                       const real_t sn_2 = 0.0 ) noexcept;

            /**
            * @brief Perform a numerical integration step.
            * @param[in] s The input signal
            * @param[in] dt The time-step given in seconds.
            * @return The current value of the integration step.
            */
            real_t integrate( const real_t s,
                              const real_t dt ) noexcept;

             /**
            * @brief Perform a numerical derivation step by using the delta rule.
            * @param[in] s The input signal
            * @param[in] dt The time-step given in seconds.
            * @return The current value of the derivation step.
            */
            real_t derivative( const real_t s,
                               const real_t dt ) noexcept;

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
            * @return @c true on success, otherwise return @c false.
            */
            inline void setIntegrationMethod( integrationMethod m ) noexcept
            {
                intMethod = m;
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
            real_t operator+( real_t rValue) const noexcept
            {
                return x[ 0 ]+rValue;
            }
            friend real_t operator*( real_t rValue,
                                     const state& s ) noexcept;
            friend real_t operator+( real_t rValue,
                                     const state& s ) noexcept;
            /*! @endcond  */
    };

    /*! @cond  */
    inline real_t operator*( real_t rValue,
                             const state& s ) noexcept
    {
        return rValue*s.x[ 0 ];
    }
    inline real_t operator+( real_t rValue,
                             const state& s ) noexcept
    {
        return rValue+s.x[ 0 ];
    }
    /*! @endcond  */

}

#endif /*QLIBS_NUMA*/