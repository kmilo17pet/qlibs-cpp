cmake_minimum_required( VERSION 3.2 )
project( qlibs-cpp
         VERSION 1.0.0
         DESCRIPTION "A collection of useful libraries for embedded systems"
         LANGUAGES CXX )

add_library( ${PROJECT_NAME} INTERFACE )

target_sources( ${PROJECT_NAME} INTERFACE
             "${CMAKE_CURRENT_LIST_DIR}/numa.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/tdl.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/ltisys.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/smoother.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/pid.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/fisCore.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/fis.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/rms.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/crc.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/fp16.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/bitfield.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/generic.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/ffmath.cpp"
             "${CMAKE_CURRENT_LIST_DIR}/interp1.cpp"
           )

target_compile_options( ${PROJECT_NAME} INTERFACE
           $<$<COMPILE_LANGUAGE:CXX>:
               -std=c++11 -Wall -Wextra -Wconversion -pedantic
               -flto -Wimplicit-fallthrough=0 -Wformat-security
               -Wduplicated-cond -Wfloat-equal -Wshadow -Wsign-conversion
               -Wlogical-not-parentheses -Wnull-dereference  -Wstringop-overflow
               -fdump-rtl-expand -fstrict-aliasing  -Wstrict-aliasing
               -Wparentheses -Wsequence-point -Wreturn-type -Wswitch -Wtrigraphs
               -Wunused -Wuninitialized -Wunknown-pragmas -Wfloat-equal -Wshadow
               -Wwrite-strings -Wsign-compare  -Wmissing-declarations -Wformat
               -Wmissing-format-attribute -Wno-deprecated-declarations -Wpacked
               -Wredundant-decls -Wunreachable-code
            >
)
set_target_properties( ${PROJECT_NAME}  PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories( ${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR} )
