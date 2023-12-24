#ifndef QLIBS_RMS
#define QLIBS_RMS

#include "include/types.hpp" 
#include "include/smoother.hpp"

namespace qlibs {

    class rms : public smootherEXPW, public smootherMWM2, public smootherLPF1 {
        public:
            virtual ~rms() {}
            rms() = default;
            real_t update( const real_t x );
            bool setup( real_t * const window, const size_t wsize );
            template <size_t windowSize>
            bool setup( real_t (&win)[ windowSize ] )
            {
                return setup( win, windowSize );
            }
            bool setParams( const real_t lambda, const real_t alpha );
    };

}

#endif /*QLIBS_RMS*/