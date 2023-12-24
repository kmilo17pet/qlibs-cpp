#ifndef QLIBS_PID
#define QLIBS_PID

#include "include/types.hpp"
#include "include/numa.hpp"
#include "include/ltisys.hpp"
#include <math.h>

namespace qlibs {

    enum class pidMode {
        PID_AUTOMATIC,
        PID_MANUAL,
    };

    enum class pidDirection {
        PID_FORWARD,
        PID_BACKWARD,
    };

    class pidController;

    struct pidGains {
        real_t Kc, Ki, Kd;
    };

    class pidAutoTuning {
        friend class pidController;
        public:
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
            static bool isValidValue( const real_t x );
        public:
            
            pidAutoTuning() = default;
            bool step( const real_t u, const real_t y, const real_t dt );
            pidGains getEstimates( const real_t dt );
    };

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
            static real_t saturate( real_t x, const real_t vMin, const real_t vMax );
            void adaptGains( const real_t u, const real_t y );
        public:
            virtual ~pidController() {}
            pidController() = default;
            bool setup( const real_t kc, const real_t ki, const real_t kd, const real_t dT );
            bool setDirection( const pidDirection d );
            bool setParams( const real_t kc, const real_t ti, const real_t td );
            bool setGains( const real_t kc, const real_t ki, const real_t kd );
            bool setExtraGains( const real_t Kw, const real_t Kt );
            bool setSaturation( const real_t Min, const real_t Max );
            bool setSeries( void );
            bool setEpsilon( const real_t eps );
            bool setDerivativeFilter( const real_t Beta );
            bool setMode( const pidMode Mode );
            bool setReferenceWeighting( const real_t gb, const real_t gc );
            bool setManualInput( const real_t manualInput );
            bool setModelReferenceControl( const real_t &modelRef, const real_t Gamma = 0.5, const real_t Alpha = 0.01 );
            bool removeModelReferenceControl( void );
            real_t control( const real_t w, const real_t y );
            bool bindAutoTuning( pidAutoTuning &at );
            bool enableAutoTuning( const uint32_t tEnable );
            bool isAutoTuningComplete( void ) const;
            bool setAutoTuningParameters( const real_t Mu, const real_t Alpha, const real_t lambda );
            bool reset( void );
    };

}

#endif /*QLIBS_PID*/