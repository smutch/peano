cmake_minimum_required(VERSION 3.8)
project(peano)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES
    src/peano.c
    src/peano.h)

add_library(peano SHARED ${SOURCE_FILES})

if(!DISABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
