#ifndef QLIBS_LTISYS
#define QLIBS_LTISYS

#include "include/types.hpp"
#include "include/tdl.hpp"
#include "include/numa.hpp"
#include <cmath>
#include <cfloat>

namespace qlibs {

    enum ltisysType {
        LTISYS_TYPE_UNKNOWN = 0,
        LTISYS_TYPE_CONTINUOUS,
        LTISYS_TYPE_DISCRETE,
    };

    using continuousStates = state[];
    using discreteStates = real_t[];

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
            void normalizeTransferFunction( real_t *num, real_t *den, size_t n_num, size_t n_den );
            real_t saturate( real_t y );
            ltisysType type{ LTISYS_TYPE_UNKNOWN };
        public:
            virtual ~ltisys() {}
            ltisys() = default;

            virtual real_t excite( real_t u ) = 0;
            virtual bool isInitialized( void ) const = 0;
            virtual bool setInitStates( const real_t *xi = nullptr ) = 0;
            ltisysType getType( void ) const
            {
                return type;
            }
            bool setDelay( real_t * const w, const size_t nD, const real_t initVal = 0.0 );
            bool setSaturation( const real_t minV, const real_t maxV );
    }; 

    class discreteSystem : public ltisys {
        private:
            real_t *xd;
            real_t update( const real_t u );
        public:
            virtual ~discreteSystem() {}
            discreteSystem( real_t *num, real_t *den, real_t *x, const size_t n_b, const size_t n_a )
            {
                (void)setup( num, den, x, n_b, n_a );
            }
            bool setup( real_t *num, real_t *den, real_t *x, const size_t n_b, const size_t n_a );
            bool isInitialized( void ) const override
            {
                return ( nullptr != xd );
            }
            bool setInitStates( const real_t *xi = nullptr ) override;
            static real_t updateFIR( real_t *w, const size_t wsize, const real_t x, const real_t * const c = nullptr );
            real_t excite( real_t u ) override;

    };

    class continuousSystem : public ltisys {
        private:
            real_t dt;
            state *xc;
            real_t update( const real_t u );
        public:
            virtual ~continuousSystem() {}
            continuousSystem( real_t *num, real_t *den, state *x, const size_t n_a, const real_t dT )
            {
                (void)setup( num, den, x, n_a, dT );
            }
            bool setup( real_t *num, real_t *den, state *x, const size_t n_a, const real_t dT );
            bool isInitialized( void ) const override
            {
                return ( nullptr != xc );
            }
            bool setInitStates( const real_t *xi = nullptr ) override;
            
            bool setIntegrationMethod( integrationMethod m );
            real_t excite( real_t u ) override;
    };
}

#endif /*QLIBS_LTISYS*/