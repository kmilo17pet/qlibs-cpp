#ifndef QLIBS_TDL
#define QLIBS_TDL

#include "include/types.hpp" 
#include <cmath>

namespace qlibs {

    class tdl : private nonCopyable {
        protected:
            real_t *head{ nullptr };
            real_t *tail{ nullptr };
            real_t *rd{ nullptr };
            real_t *wr{ nullptr };
            size_t itemCount{ 0U };
            const real_t undefined{ nan("" ) }; // skipcq: CXX-W2010

            void insertNewest( const real_t sample );
            void removeOldest( void );
        public:
            virtual ~tdl() {}
            tdl() = default;
            tdl( real_t * const area, const size_t n, const real_t initVal = 0.0 )
            {
                setup( area, n, initVal );
            }
            template <size_t numberOfDelays>
            tdl( real_t (&area)[ numberOfDelays ], const real_t initVal = 0.0 )
            {
                setup( area, numberOfDelays, initVal );
            }
            void setup( real_t * const area, const size_t n, const real_t initVal = 0.0 ) noexcept;
            template <size_t numberOfDelays>
            void setup( real_t (&area)[ numberOfDelays ], const real_t initVal = 0.0 ) noexcept
            {
                setup( area, numberOfDelays, initVal );
            }
            void flush( const real_t initVal = 0.0 );
            real_t getOldest( void ) const;
            real_t getRecent( void ) const;
            real_t getAtIndex( const size_t i ) const;
            void insertSample( const real_t sample );
            const real_t& operator[]( int index );
            void operator()( const real_t sample )
            {
                insertSample( sample );
            }
            bool isInitialized( void ) const {
                return ( nullptr != head );
            }
    };
}

#endif /*QLIBS_TDL*/