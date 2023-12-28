/*!
 * @file qlibs_types.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Common types for qLibs-cpp.
 **/

#ifndef QLIBS_TYPES
#define QLIBS_TYPES

#if defined( ARDUINO ) || defined( ARDUINO_ARCH_AVR) || defined( ARDUINO_ARCH_SAMD ) || defined( ENERGIA_ARCH_MSP430ELF )
    #define ARDUINO_PLATFORM
#endif


#if ( __cplusplus < 201103L || defined ( __AVR_ARCH__ ) || defined( ARDUINO_PLATFORM ) || defined ( STM8Sxx ) )
    #include <stddef.h>
    #include <stdint.h>
    #include <stdlib.h>
    #include <string.h>
    #include <ctype.h>
    #include <math.h>
    #include <limits.h>
    #include <float.h>
#else 
    #include <cstddef>
    #include <cstdint>
    #include <cstdlib>
    #include <cstring>
    #include <cctype>
    #include <cmath>
    #include <climits>
    #include <cfloat>
    using namespace std;
#endif

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /**
    * @brief A type to instantiate a real variable double-precision of 64-bits IEEE 754.
    */
    using real_t = double;

    /*! @cond */
    class nonCopyable {
        protected:
            nonCopyable() {}
            ~nonCopyable() {}
        private:
            nonCopyable( const nonCopyable & );
            nonCopyable& operator=( const nonCopyable & );
    };
    /*! @endcond */
}

/*! @cond */
#define     Q_UNUSED(arg)     (void)(arg)
#define     Q_NONE            /*EMPTY MACRO*/
/*! @endcond */



#endif /*QOS_CPP_TYPES*/
