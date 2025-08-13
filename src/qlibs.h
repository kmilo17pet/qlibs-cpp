/*!
 * @file qlibs.h
 * @author J. Camilo Gomez C.
 * @version 1.4.0
 * @note This file is part of the qlibs++ distribution.
 * @brief Global inclusion header
 **/


/*
qlibs++ - A collection of useful C++ libraries for embedded systems.
MIT License
C++11 and MISRA C++ 2008 / CERT Compliant

Copyright (C) 2012 Eng. Juan Camilo Gomez Cadavid MSc. All Rights Reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

VISIT https://github.com/kmilo17pet/QuarkTS-cpp TO ENSURE YOU ARE USING THE LATEST
VERSION.


This file is part of the QuarkTS++ OS distribution.
*/

#ifndef QLIBS_CPP_H
#define QLIBS_CPP_H

    #define QLIBS_CPP_VERSION         "1.4.0"
    #define QLIBS_CPP_VERNUM          ( 140 )
    #define QLIBS_CPP_CAPTION         "qLibs++ " QLIBS_CPP_VERSION

    #include <include/qlibs_types.hpp>
    #include <include/tdl.hpp>
    #include <include/numa.hpp>
    #include <include/ltisys.hpp>
    #include <include/smoother.hpp>
    #include <include/fis.hpp>
    #include <include/pid.hpp>
    #include <include/crc.hpp>
    #include <include/rms.hpp>
    #include <include/fp16.hpp>
    #include <include/ffmath.hpp>
    #include <include/bitfield.hpp>
    #include <include/interp1.hpp>
    #include <include/algorithm.hpp>

    namespace qlibs {
        namespace build {
            constexpr const uint32_t number = 2458;
            constexpr const char* date = __DATE__;
            constexpr const char* time = __TIME__;
            constexpr const char* std = "c++11";
        }
        namespace version {
            constexpr const char* str = QLIBS_CPP_VERSION;
            constexpr const uint8_t number = QLIBS_CPP_VERNUM;
            constexpr const uint8_t mayor = 1U;
            constexpr const uint8_t minor = 4U;
            constexpr const uint8_t rev = 0U;
        }
        namespace product {
            constexpr const char* author = "J. Camilo Gomez C.";
            constexpr const char* copyright = "Copyright (C) 2012 J. Camilo Gomez C. All Rights Reserved.";
            constexpr const char* name = "qLibs++";
            constexpr const char* category = "Library";
            constexpr const char* caption = QLIBS_CPP_CAPTION;
            constexpr const char* compliance = "MISRAC++2008,SEI-CERT";
            constexpr const char* license = "MIT";
            constexpr const char* source_model = "Open Source";
        }
    }

#endif /*QOS_CPP_H*/
