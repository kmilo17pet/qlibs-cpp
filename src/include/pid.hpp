/*!
 * @file pid.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief API to control systems using the PID algorithm. This controller
 * features anti-windup, tracking mode, and derivative filter.
 **/

#ifndef QLIBS_PID
#define QLIBS_PID

#include <include/qlibs_types.hpp>
#include <include/numa.hpp>

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /** @addtogroup qpid PID Controller
    * @brief Closed Loop PID Controller class
    * @{
    */

    /**
    * @brief Enumeration class with the operational modes for the PID controller
    */
    enum class pidMode {
        PID_AUTOMATIC,  /*!< Fully operational closed-loop PID controller */
        PID_MANUAL,     /*!< Open-loop with manual input*/
    };

    /**
    * @brief Enumeration class with the operational modes for the PID controller
    */
    enum class pidType {
        PID_TYPE_P,     /*!< Proportional controller */
        PID_TYPE_PI,    /*!< Proportional-Integral controller */
        PID_TYPE_PD,    /*!< Proportional-Derivative controller */
        PID_TYPE_PID,   /*!< Proportional-Integral-Derivative controller */
    };

    /**
    * @brief Direction modes of the PID controller
    */
    enum class pidDirection {
        PID_FORWARD,    /*!< Forward control action */
        PID_BACKWARD,   /*!< Reverse the control action*/
    };

    /*! @cond  */
    class pidController;
    /*! @endcond  */

    /**
    * @brief PID Gains structure
    */
    struct pidGains {
        real_t Kc{ 0.0_re };  /*!< Proportional gain */
        real_t Ki{ 0.0_re };  /*!< Integral gain */
        real_t Kd{ 0.0_re };  /*!< Derivative gain */

        /*! @cond  */
        constexpr pidGains() = default;
        constexpr pidGains(real_t kc, real_t ki, real_t kd)
            : Kc(kc), Ki(ki), Kd(kd) {}

        constexpr pidGains operator+(const pidGains& other) const {
            return pidGains{
                this->Kc + other.Kc,
                this->Ki + other.Ki,
                this->Kd + other.Kd
            };
        }
        /*! @endcond  */
    };

    /*! @cond  */
    constexpr pidGains operator"" _kc(long double v) {
        return {static_cast<real_t>(v), 0.0_re, 0.0_re};
    }
    constexpr pidGains operator"" _ki(long double v) {
        return {0.0_re, static_cast<real_t>(v), 0.0_re};
    }
    constexpr pidGains operator"" _kd(long double v) {
        return {0.0_re, 0.0_re, static_cast<real_t>(v)};
    }
    /*! @endcond  */

    /**
    * @brief A PID Auto-tuning object
    * @details The instance should be bound to a configured PID controller by
    * using the pidController::bindAutoTuning() method
    */
    class pidAutoTuning {
        friend class pidController;
        public:
            /**
            * @brief Constant to keep the auto-tuner enabled indefinitely
            */
            static const uint32_t UNDEFINED;
        protected:
            /*! @cond  */
            real_t p00{ 1.0_re };       /*covariance value*/
            real_t p01{ 0.0_re };       /*covariance value*/
            real_t p10{ 0.0_re };       /*covariance value*/
            real_t p11{ 1.0_re };       /*covariance value*/
            real_t b1{ 0.1_re };        /*estimation value*/
            real_t a1{ 0.9_re };        /*estimation value*/
            real_t uk{ 0.0_re };        /*process input*/
            real_t yk{ 0.0_re };        /*process output*/
            real_t l{ 0.9898_re };      /*memory factor [ 0.9 < l < 1 ]*/
            real_t k{ 1.0_re };         /*process static gain*/
            real_t tao{ 1.0_re };       /*process time constant*/
            real_t mu{ 0.95_re };       /*variation attenuation*/
            real_t speed{ 0.25_re };    /*final controller speed*/
            uint32_t it{ UNDEFINED };   /*enable time*/
            static bool isValidValue( const real_t x ) noexcept;
            pidType type{ pidType::PID_TYPE_PI };
            void initialize( const pidGains current,
                             const real_t dt ) noexcept;
            inline void enable( const uint32_t tEnable ) noexcept
            {
                it = ( 0UL == tEnable ) ? pidAutoTuning::UNDEFINED : tEnable;
            }
            inline bool isComplete( void ) const noexcept
            {
                return ( ( 0UL == it ) && ( it != pidAutoTuning::UNDEFINED ) );
            }
            inline void setMemoryFactor( const real_t lambda ) noexcept
            {
                l = lambda;
            }
            inline void setMomentum( const real_t Mu ) noexcept
            {
                mu = Mu;
            }
            inline void setEstimatedControllerSpeed( const real_t alpha ) noexcept
            {
                speed = alpha;
            }
            inline void setEstimatedControllerType( const pidType t ) noexcept
            {
                type = t;
            }
            static inline bool isValidParam( const real_t p ) noexcept
            {
                return ( p > 0.0_re ) && ( p <= 1.0_re );
            }
            bool step( const real_t u,
                       const real_t y,
                       const real_t dt ) noexcept;
            /*! @endcond  */
        public:
            pidAutoTuning() = default;
            pidGains getEstimates( void ) const noexcept;
    };

    /**
    * @brief A PID controller object
    * @details The instance should be initialized using the pid::setup() method.
    */
    class pidController : public pidGains, public nState, private nonCopyable {
        private:
            real_t b, c, sat_Min, sat_Max, epsilon, kw, kt, D, u1, beta, uSat, gainBlend;
            real_t dt{ 1.0_re };
            pidGains nextGains;
            real_t m, mInput;
            const real_t *yr{ nullptr };
            real_t alpha, gamma; /*MRAC additive controller parameters*/
            nState m_state; /*MRAC additive controller state*/
            nState b_state; /*Bumpless-transfer state*/
            pidAutoTuning *adapt{ nullptr };
            pidMode mode{ pidMode::PID_AUTOMATIC };
            pidDirection dir{ pidDirection::PID_FORWARD };
            bool isInitialized{ false };
            real_t error( real_t w, real_t y, real_t k = 1.0_re ) noexcept;
            static real_t saturate( real_t x,
                                    const real_t vMin,
                                    const real_t vMax ) noexcept;
            void adaptGains( const real_t u,
                             const real_t y ) noexcept;
        public:
            virtual ~pidController() {}
            pidController() = default;

            /**
            * @brief Setup and initialize the PID controller instance.
            * @param[in] kc Proportional Gain.
            * @param[in] ki Integral Gain.
            * @param[in] kd Derivative Gain.
            * @param[in] dT Time step in seconds.
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const real_t kc,
                        const real_t ki,
                        const real_t kd,
                        const real_t dT ) noexcept;

            /**
            * @brief Setup and initialize the PID controller instance.
            * @param[in] g The structure with the controller gains
            * @param[in] dT Time step in seconds.
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const pidGains &g,
                        const real_t dT ) noexcept
            {
                return setup( g.Kc, g.Ki, g.Kd, dT );
            }

            /**
            * @brief Set the PID control action direction.
            * @param[in] d Desired output direction.
            * @return @c true on success, otherwise return @c false.
            */
            bool setDirection( const pidDirection d ) noexcept;

            /**
            * @brief Set/Change the PID controller gains by using the [Kc, Ti Td ]
            * triplet.
            * @param[in] kc Proportional Gain.
            * @param[in] ti Integral time.
            * @param[in] td Derivative time.
            * @return @c true on success, otherwise return @c false.
            */
            bool setParams( const real_t kc,
                            const real_t ti,
                            const real_t td ) noexcept;

            /**
            * @brief Set/Change the PID controller gains.
            * @param[in] kc Proportional Gain.
            * @param[in] ki Integral Gain.
            * @param[in] kd Derivative Gain.
            * @return @c true on success, otherwise return @c false.
            */
            bool setGains( const real_t kc,
                           const real_t ki,
                           const real_t kd ) noexcept;

            /**
            * @brief Set/Change the PID controller gains.
            * @param[in] g The structure with the controller gains
            * @return @c true on success, otherwise return @c false.
            */
            bool setGains( const pidGains &g ) noexcept;

            /**
            * @brief Set/Change extra PID controller gains.
            * @param[in] Kw Saturation feedback gain. Used for antiWindup and bumpless
            * transfer. A zero value disables these features.
            * @param[in] Kt Manual input gain.
            * @return @c true on success, otherwise return @c false.
            */
            bool setExtraGains( const real_t Kw,
                                const real_t Kt ) noexcept;

            /**
            * @brief Setup the output saturation for the PID controller.
            * @param[in] Min The minimal value allowed for the output.
            * @param[in] Max The maximal value allowed for the output.
            * @return @c true on success, otherwise return @c false.
            */
            bool setSaturation( const real_t Min,
                                const real_t Max ) noexcept;

            /**
            * @brief Convert the controller gains to conform the series or
            * interacting form.
            * @return @c true on success, otherwise return @c false.
            */
            bool setSeries( void ) noexcept;

            /**
            * @brief Set the minimum value considered as error.
            * @param[in] eps The minimal error value ( eps > 0 ).
            * @return @c true on success, otherwise return @c false.
            */
            bool setEpsilon( const real_t eps ) noexcept;

            /**
            * @brief Set the tuning parameter for the derivative filter.
            * @param[in] Beta The tuning parameter. [ 0 <= Beta < 1 ]
            * @return @c true on success, otherwise return @c false.
            */
            bool setDerivativeFilter( const real_t Beta ) noexcept;

            /**
            * @brief Set the time constant for the derivative filter.
            * @note
            * @param[in] Tf Derivative filter time constant [ Tf >= 0 ]
            * @return @c true on success, otherwise return @c false.
            */
            bool setDerivativeFilterTimeConstant( const real_t Tf ) noexcept;

            /**
            * @brief Change the controller operational mode.
            * In pidMode::PID_AUTOMATIC, the computed output of the PID controller
            * will be used as the control signal for the process. In
            * pidMode::PID_MANUAL mode, the manual input will be used as the control
            * signal for the process and the PID controller loop will continue
            * operating to guarantee the bumpless-transfer when a switch to
            * the pidMode::PID_AUTOMATIC is performed;
            * @param[in] Mode The desired operational mode.
            * @return @c true on success, otherwise return @c false.
            */
            bool setMode( const pidMode Mode ) noexcept;

            /**
            * @brief Set the PID Reference(Set-Point) Weighting. This value is used in
            * order to avoid the increase of the rise time due to the smoothing of the
            * reference signal applied to the closed-loop system
            * @note A value close to zero en @a gc can be used to eliminate the
            * reduce the effect of the phenomena called "derivative-kick"
            * @param[in] gb The reference weight value for the proportional element.
            * @param[in] gc The reference weight value for the derivative element.
            * @return @c true on success, otherwise return @c false.
            */
            bool setReferenceWeighting( const real_t gb,
                                        const real_t gc ) noexcept;

            /**
            * @brief Set the PID manual input mode. This value will be used
            * as the manual input when the controller it set into the pidMode::PID_MANUAL
            * mode. Bumpless-transfer is guaranteed.
            * @param[in] manualInput The value of the manual input.
            * @return @c true on success, otherwise return @c false.
            */
            bool setManualInput( const real_t manualInput ) noexcept;

            /**
            * @brief Enable the additive MRAC(Model Reference Adaptive Control) feature.
            * @param[in] modelRef A pointer to the output of the model reference.
            * @param[in] Gamma Adjustable parameter to indicate the adaptation speed.
            * @param[in] Alpha Adjustable parameter for the adaptation momentum.
            * @return @c true on success, otherwise return @c false.
            */
            bool setModelReferenceControl( const real_t &modelRef,
                                           const real_t Gamma = 0.5_re,
                                           const real_t Alpha = 0.01_re ) noexcept;

            /**
            * @brief Removes the additive MRAC(Model Reference Adaptive Control)
            * feature from the control loop.
            * MRAC(Model Reference Adaptive Control) feature from the control loop.
            * @return @c true on success, otherwise return @c false.
            */
            bool removeModelReferenceControl( void ) noexcept;

            /**
            * @brief Computes the control action for given PID controller instance.
            * @pre Instance must be previously initialized by pidController::setup()
            * @note The user must ensure that this function is executed in the time
            * specified in @a dt either by using a HW or SW timer, a real time task,
            * or a timing service.
            * @param[in] w The reference value aka SetPoint.
            * @param[in] y The controlled variable aka Process-variable.
            * @return The control action.
            */
            real_t control( const real_t w,
                            const real_t y ) noexcept;


            /**
            * @brief Computes the control action for given PID controller instance.
            * @pre Instance must be previously initialized by pidController::setup()
            * @note The user must ensure that this function is executed in the time
            * specified in @a dt either by using a HW or SW timer, a real time task,
            * or a timing service.
            * @param[in] w The reference value aka SetPoint.
            * @param[in] y The controlled variable aka Process-variable.
            * @return The control action.
            */
            real_t operator()( const real_t w,
                               const real_t y )
            {
                return control( w, y );
            }

            /**
            * @brief Binds the specified instance to enable the PID controller auto
            * tuning algorithm.
            * @param[in] at The PID auto-tuning instance.
            * @return @c true on success, otherwise return @c false.
            */
            bool bindAutoTuning( pidAutoTuning &at ) noexcept;

            /**
            * @brief Unbind and disable the PID controller auto-tuning algorithm.
            * @return @c true on success, otherwise return @c false.
            */
            inline bool unbindAutoTuning( void ) noexcept
            {
                const bool retValue = nullptr != adapt;
                adapt = nullptr;
                return retValue;
            }

            /**
            * @brief Set the number of time steps where the auto tuner algorithm will
            * modify the controller gains.
            * @pre Controller must have a pidAutoTuning object already binded with
            * pidController::bindAutoTuning()
            * @note To disable auto-tuning pass a 0uL value to the @a tEnable argument.
            * @param[in] tEnable The number of time steps. To keep the auto tuner
            * enabled indefinitely pass pidAutoTuning::UNDEFINED as argument.
            * @return @c true on success, otherwise return @c false.
            */
            bool enableAutoTuning( const uint32_t tEnable ) noexcept;

            /**
            * @brief Verifies that the auto tuning process has finished with new
            * gains set on the controller
            * @return @c true if auto-tuning its complete, otherwise return @c false.
            */
            bool isAutoTuningComplete( void ) const noexcept;

            /**
            * @brief Change parameters of the auto-tuning algorithm.
            * @param[in] Mu Algorithm momentum. [ 0 <= Mu <= 1 ].
            * @param[in] Alpha Final controller speed adjustment. [ 0 < Alpha <= 1 ].
            * @param[in] lambda Algorithm forgetting factor [ 0.8 <= lambda <= 1 ].
            * @param[in] Tb Gains blend-time[  Tb > 0  ].
            * @return @c true on success, @c false on failure.
            */
            bool setAutoTuningParameters( const real_t Mu,
                                          const real_t Alpha,
                                          const real_t lambda,
                                          const real_t Tb = 3.0_re ) noexcept;

            /**
            * @brief Select the PID type to tune.
            * @param[in] t The type of controller to tune.
            * @return @c true on success, @c false on failure.
            */
            bool setAutoTuningControllerType( const pidType t ) noexcept;

            /**
            * @brief Reset the internal PID controller calculations.
            * @return @c true on success, otherwise return @c false.
            */
            bool reset( void ) noexcept;

            /**
            * @brief Retrieve the current PID gains.
            * @return A struct with the pid gains @a Kc @a Ki and @a Kd
            */
            inline pidGains getGains( void ) const noexcept
            {
                return { Kc, Ki, Kd };
            }

            /**
            * @brief Check if the PID instance has been initialized using setup().
            * @return @c true if the PID instance has been initialized, otherwise
            * return @c false.
            */
            explicit operator bool() const noexcept {
                return isInitialized;
            }
    };

    /** @}*/
}

#endif /*QLIBS_PID*/