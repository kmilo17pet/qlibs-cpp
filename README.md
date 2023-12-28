[![Built for](https://img.shields.io/badge/built%20for-microcontrollers-lightgrey?logo=WhiteSource)](https://github.com/kmilo17pet/qlibs-cpp)
[![DeepSource](https://app.deepsource.com/gh/kmilo17pet/qlibs-cpp.svg/?label=active+issues&show_trend=true&token=v2gmYWIv1qjuUS9q01v_ncon)](https://app.deepsource.com/gh/kmilo17pet/qlibs-cpp/)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2773da9d554f4f8abc31e76b7bcd33c3)](https://app.codacy.com/gh/kmilo17pet/qlibs-cpp/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![CodeFactor](https://www.codefactor.io/repository/github/kmilo17pet/qlibs-cpp/badge)](https://www.codefactor.io/repository/github/kmilo17pet/qlibs-cpp)
[![MISRAC++2008](https://img.shields.io/badge/MISRAC++2008-Compliant-blue.svg?logo=verizon)](https://www.misra.org.uk/)
[![CERT](https://img.shields.io/badge/CERT-Compliant-blue.svg?logo=cplusplus)](https://wiki.sei.cmu.edu/confluence/display/seccode/SEI+CERT+Coding+Standards)
[![C++ Standard](https://img.shields.io/badge/STD-C++11-green.svg?logo=cplusplus)](https://en.cppreference.com/w/cpp/11)
[![arduino-library-badge](https://www.ardu-badge.com/badge/qlibs.svg?)](https://www.ardu-badge.com/qlibs)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg?logo=textpattern)](https://github.com/kmilo17pet/QuarkTS-cpp/graphs/commit-activity)
[![License](https://img.shields.io/github/license/kmilo17pet/QuarkTS-cpp?logo=livejournal)](https://github.com/kmilo17pet/QuarkTS-cpp/blob/master/LICENSE)

# qlibs++
qlibs++ is the [qlibs](https://github.com/kmilo17pet/qlibs)  port for C++.

![qlibs++logo](https://github.com/kmilo17pet/qlibs-cpp/assets/11412210/76b07eb7-b522-449e-be1b-c4349bb02dd6)
# qlibs++ : A collection of useful libraries for embedded systems


* Download the latest release [here](https://github.com/kmilo17pet/qlibs-cpp/releases)
* Documentation and API Reference [here](https://kmilo17pet.github.io/qlibs-cpp/)


Below is the list of the classes and modules provided and their features:

- smoother : Filters to smooth noisy signals
   - `LPF1`: _Low Pass Filter Order 1_
   - `LPF2`: _Low Pass Filter Order 2_
   - `MWM1`: _Moving Window Median O(n)_
   - `MWM2`: _Moving Window Median O(1): With TDL(works efficient for large windows)_
   - `MOR1`: _Moving Outlier Removal O(n)_
   - `MOR2`: _Moving Outlier Removal O(1): With TDL(works efficient for large windows)_
   - `GMWF`: _Gaussian filter_
   - `KLMN`: _Scalar Kalman filter_
   - `EXPW`: _Exponential weighting filter_
- pidController : Closed Loop PID Controller
  - Derivative filter
  - Anti-Windup
  - Tracking Mode
  - Auto-tunning 
  - Additive MRAC
- ltisys : Recursive LTI systems evaluation by transfer functions
  - Continuous
  - Discrete
- fis : Fuzzy Inference System Engine
  - Mamdani
  - Sugeno
  - Tsukamoto
- fp16 : Q16.16 Fixed-point math
  - Basic operations
  - Trigonometric functions
  - Exponential functions
- crc : Generic Cyclic Redundancy Check (CRC) calculator
  - CRC8
  - CRC16
  - CRC32
- bitfield: A bit-field manipulation library
- tdl : Tapped Delay Line in O(1). 
- rms : Recursive Root Mean Square(RMS) calculation of a signal.
- Type-generic array operations
- Single precision floating-point vector(1D-Array) operations
- Fast single precision floating-point math
