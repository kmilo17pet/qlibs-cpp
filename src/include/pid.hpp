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

#include "include/types.hpp"
#include "include/numa.hpp"
#include "include/ltisys.hpp"

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
        real_t Kc;  /*!< Proportional gain */
        real_t Ki;  /*!< Integral gain */
        real_t Kd;  /*!< Derivative gain */
    };

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
            real_t p00{ 1.0 };
            real_t p01{ 0.0 };
            real_t p10{ 0.0 };
            real_t p11{ 1.0 };   /*covariance values*/
            real_t b1{ 0.1 };
            real_t a1{ 0.9 };   /*estimation  values*/
            real_t uk{ 0.0 };
            real_t yk{ 0.0 };      /*process I/O measurements*/
            real_t l{ 0.9898 }; /*memory factor [ 0.9 < l < 1 ]*/
            real_t il{ 1.0 };
            real_t k{ 1.0 };
            real_t tao{ 1.0 };      /*process metrics*/
            real_t mu{ 0.95};
            real_t speed{ 0.25 };   /*fine adjustments  [ 0 < mu < speed ] [ 0 < speed < 1 ]*/
            uint32_t it{ UNDEFINED };/*enable time*/
            static bool isValidValue( const real_t x ) noexcept;
        public:
            pidAutoTuning() = default;
            bool step( const real_t u,
                       const real_t y,
                       const real_t dt ) noexcept;
            pidGains getEstimates( const real_t dt ) const noexcept;
    };

    /**
    * @brief A PID controller object
    * @details The instance should be initialized using the pid::setup() method.
    */
    class pidController : public pidGains, private nonCopyable {
        private:
            //real_t Kc, Ki, Kd;
            real_t b, c, min, max, epsilon, kw, kt, D, u1, beta;
            real_t dt{ 1.0 };
            real_t m, mInput;
            const real_t *yr{ nullptr };
            real_t alpha, gamma; /*MRAC additive controller parameters*/
            state c_state; /*controller integral & derivative state*/
            state m_state; /*MRAC additive controller state*/
            state b_state; /*Bumpless-transfer state*/
            pidAutoTuning *adapt{ nullptr };
            pidMode mode{ pidMode::PID_AUTOMATIC };
            pidDirection dir{ pidDirection::PID_FORWARD };
            bool init{ false };
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
            * @param[in] eps The minimal error value.
            * @return @c true on success, otherwise return @c false.
            */
            bool setEpsilon( const real_t eps ) noexcept;

            /**
            * @brief Set the tuning parameter for the derivative filter.
            * @param[in] Beta The tuning parameter. [ 0 < Beta < 1 ]
            * @return @c true on success, otherwise return @c false.
            */
            bool setDerivativeFilter( const real_t Beta ) noexcept;

            /**
            * @brief Change the controller operational mode.
            * In pidMode::PID_AUTOMATIC, the computed output of the PID controller
            * will be used as the control signal for the process. In
            * pidMode::PID_MANUAL mode, the manual input will be used as the control
            * signal for the process and the PID controller loop will continue
            * operating to guarantee the bumpless-transfer when a switch to 
            * the pidMode::PID_AUTOMATIC its performed;
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
            * @param[in] Alpha Adjustable parameter the adaptation momentum.
            * @return @c true on success, otherwise return @c false.
            */
            bool setModelReferenceControl( const real_t &modelRef,
                                           const real_t Gamma = 0.5,
                                           const real_t Alpha = 0.01 ) noexcept;

            /**
            * @brief Removes the Enable the additive 
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
                bool retValue = nullptr != adapt;
                adapt = nullptr;
                return retValue;
            }

            /**
            * @brief Set the number of time steps where the auto tuner algorithm will
            * modify the controller gains.
            * @pre Controller must have an pidAutoTuning object already binded wih
            * pidController::bindAutoTuning()
            * @note To disable auto-tuning pass a 0uL value to the @a tEnable argument.
            * @param[in] tEnable The number of time steps. To keep the auto tuner
            * enabled indefinitely pass pidAutoTuning::UNDEFINED as argument.
            * @return @c true on success, otherwise return @c false.
            */
            bool enableAutoTuning( const uint32_t tEnable ) noexcept;

            /**
            * @brief Verifies that the auto tuning process has finished with new
            * gains  set on the controller
            * @return @c true if auto-tuning its complete, otherwise return @c false.
            */
            bool isAutoTuningComplete( void ) const noexcept;

            /**
            * @brief Change parameters of the auto-tuning algorithm.
            * @param[in] Mu Algorithm momentum. [ 0 <= Mu <= 1 ].
            * @param[in] Alpha Final controller speed adjustment. [ 0 < Alpha <= 1 ].
            * @param[in] lambda Algorithm forgetting factor [ 0.8 <= lambda <= 1 ].
            * @return @c true on success, @c false on failure.
            */
            bool setAutoTuningParameters( const real_t Mu,
                                          const real_t Alpha,
                                          const real_t lambda ) noexcept;

            /**
            * @brief Reset the internal PID controller calculations.
            * @return @c true on success, otherwise return @c false.
            */
            bool reset( void ) noexcept;
    };

    /** @}*/
}

#endif /*QLIBS_PID*/