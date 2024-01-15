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
#else
    #include <cstddef>
    #include <cstdint>
    #include <cstdlib>
    #include <cstring>
    #include <cctype>
    using namespace std;
#endif

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /**
    * @brief A type to instantiate a real variable double-precision of 64-bits IEEE 754.
    */
    using real_t = float;

    /*! @cond */
    /*cstat -CERT-FLP36-C -CERT-FLP34-C*/
    constexpr real_t operator "" _re( unsigned long long int x )
    {
        return static_cast<real_t>( x );
    }
    constexpr real_t operator "" _re( long double x )
    {
        return static_cast<real_t>( x );
    }
    /*cstat +CERT-FLP36-C +CERT-FLP34-C*/

    /*cstat -MISRAC++2008-0-1-4_b*/
    constexpr real_t REAL_MAX = 3.402823466e+38F;        // max value
    constexpr real_t REAL_MIN = 1.175494351e-38F;        // min normalized positive value
    /*cstat +MISRAC++2008-0-1-4_b*/

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
