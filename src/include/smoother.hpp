#ifndef QLIBS_SMOOTHER
#define QLIBS_SMOOTHER

#include "include/types.hpp" 
#include "include/tdl.hpp" 
#include <cmath>

namespace qlibs {

    class smoother {
        protected:
            bool init{ true };
            static void windowSet( real_t *w, const size_t wsize, const real_t x );
        public:
            virtual ~smoother() {}
            bool isInitialized( void ) const;
            virtual real_t smooth( const real_t x ) = 0;
            bool reset( void )
            {
                init = true;
                return true;
            }
    };

    class smootherLPF1 : public smoother {
        protected:
            real_t alpha{ 0.9 };
            real_t y1{ 0.0 };
        public:
            virtual ~smootherLPF1() {}
            bool setup( const real_t a = 0.9 );
            real_t smooth( const real_t x ) override;
    };

    class smootherLPF2 : public smoother {
        protected:
            real_t y1{ 0.0 };
            real_t y2{ 0.0 };
            real_t x1{ 0.0 };
            real_t x2{ 0.0 };
            real_t k{ 0.0 };
            real_t a1{ 0.0 };
            real_t a2{ 0.0 };
            real_t b1{ 0.0 };
        public:
            virtual ~smootherLPF2() {}
            bool setup( const real_t a = 0.9 );
            real_t smooth( const real_t x ) override;
    };

    class smootherMWM1 : public smoother, private nonCopyable {
        protected:
            real_t *w{ nullptr };
            size_t wsize;
        public:
            virtual ~smootherMWM1() {}
            bool setup( real_t *window, const size_t w_size );
            real_t smooth( const real_t x ) override;
    };

    class smootherMWM2 : public smoother, public tdl {
        protected:
            real_t sum{ 0.0 };
        public:
            virtual ~smootherMWM2() {}
            bool setup( real_t *window, const size_t w_size );
            real_t smooth( const real_t x ) override;
    };

    class smootherMOR1 : public smoother, private nonCopyable {
        protected:
            real_t *w{ nullptr };
            real_t m{ 0.0 };
            real_t alpha{ 0.0 };
            size_t wsize;
        public:
            virtual ~smootherMOR1() {}
            bool setup( real_t *window, const size_t w_size, const real_t a = 0.9 );
            real_t smooth( const real_t x ) override;
    };

    class smootherMOR2 : public smoother, public tdl {
        protected:
            real_t sum{ 0.0 };
            real_t m{ 0.0 };
            real_t alpha;
        public:
            virtual ~smootherMOR2() {}
            bool setup( real_t *window, const size_t w_size, const real_t a = 0.9 );
            real_t smooth( const real_t x ) override;
    };

    class smootherGMWF : public smoother, private nonCopyable {
        protected:
            real_t *w{ nullptr };
            real_t *k{ nullptr };
            size_t wsize;
        public:
            virtual ~smootherGMWF() {}
            bool setup( const real_t sg, const size_t c, real_t *window, const size_t w_size );
            real_t smooth( const real_t x ) override;
    };

    class smootherEXPW : public smoother {
        protected:
            real_t lambda{ 0U };
            real_t m{ 0U };
            real_t w{ 0U };
        public:
            virtual ~smootherEXPW() {}
            bool setup( const real_t lam );
            real_t smooth( const real_t x ) override;
    };

    class smootherKLMN : public smoother {
        protected:
            real_t xS;  /* state */
            real_t A{ 1.0 };  /* x(n)=A*x(n-1)+u(n),u(n)~N(0,q) */
            real_t H{ 1.0 };  /* z(n)=H*x(n)+w(n),w(n)~N(0,r) */
            real_t q;  /* process(predict) noise covariance */
            real_t r;  /* measure noise covariance */
            real_t p;  /* estimated error covariance */
            real_t gain;
        public:
            virtual ~smootherKLMN() {}
            bool setup( const real_t processNoiseCov, const real_t measureNoiseCov, const real_t estErrorCov );
            real_t smooth( const real_t x ) override;
    };

    class smootherDESF : public smoother {
        protected:
            real_t alpha;
            real_t beta;
            real_t lt;
            real_t bt;
            real_t n;
        public:
            virtual ~smootherDESF() {}
            bool setup( const real_t a, const real_t b, const real_t N );
            real_t smooth( const real_t x ) override;
    };

    class smootherALNF : public smoother, private nonCopyable {
        protected:
            real_t alpha;
            real_t mu;
            real_t *w{ nullptr };
            real_t *w_1{ nullptr };
            real_t *xx{ nullptr };
            size_t n{ 0U };
        public:
            virtual ~smootherALNF() {}
            bool setup( const real_t a, const real_t m, real_t *window, const size_t wsize );
            real_t smooth( const real_t x ) override;
    };


}

#endif /*QLIBS_SMOOTHER*/