#ifndef QLIBS_MATHEX
#define QLIBS_MATHEX

#include "include/types.hpp"
#include <cfloat>

namespace qlibs {

    real_t mapMinMax( const real_t x,
                      const real_t xMin,
                      const real_t xMax,
                      const real_t yMin,
                      const real_t yMax );

    real_t normalize( const real_t x,
                      const real_t xMin,
                      const real_t xMax );

    bool inRangeCoerce( real_t &x,
                        const real_t lowerL,
                        const real_t upperL );

    bool isEqual( const real_t a, const real_t b, const real_t tol = DBL_MIN );

    bool inPolygon( const real_t x,
                    const real_t y,
                    const real_t * const px,
                    const real_t * const py,
                    const size_t p );

    bool inCircle( const real_t x,
                   const real_t y,
                   const real_t cx,
                   const real_t cy,
                   const real_t r );

}

#endif /*QLIBS_MATHEX*/