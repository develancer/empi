cmake_minimum_required(VERSION 3.15)
project(empi CXX C)
set(CMAKE_CXX_STANDARD 17)

enable_testing()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, defaults to Release." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

if(APPLE)
    execute_process(
        COMMAND brew --prefix fftw
        OUTPUT_VARIABLE FFTW_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    include_directories("${FFTW_PREFIX}/include")
    link_directories("${FFTW_PREFIX}/lib")
endif()

if(NOT UNIX)
    find_package(FFTW3 CONFIG REQUIRED)
    find_path(FFTW_INCLUDE_DIR NAMES fftw3.h)
    add_library(fftw3 ALIAS FFTW3::fftw3)
    include_directories(${FFTW_INCLUDE_DIR})
endif()

find_package(Git)
set(VERSION_FILE "${CMAKE_SOURCE_DIR}/version")

if(EXISTS ${VERSION_FILE})
    file(READ ${VERSION_FILE} FILE_VERSION)
    string(STRIP ${FILE_VERSION} APP_VERSION)
endif()
if(NOT DEFINED APP_VERSION AND GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_TAG
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE GIT_TAG_RESULT
    )
    if(GIT_TAG_RESULT EQUAL "0")
        set(APP_VERSION ${GIT_TAG})
    endif()
endif()
if(NOT DEFINED APP_VERSION)
    set(APP_VERSION "1.?.?")
endif()
add_definitions("-DAPP_VERSION=\"${APP_VERSION}\"")
add_definitions("-D_USE_MATH_DEFINES")

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    set(WITH_CUDA TRUE CACHE BOOL "Compile with CUDA support?")
elseif(WITH_CUDA)
    message(SEND_ERROR "CUDA is not available on this system")
endif()
if(WITH_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 11)
    add_compile_definitions(HAVE_CUDA)
    # link all executables with nvcc
    string(REPLACE <CMAKE_CXX_COMPILER> nvcc CMAKE_CXX_LINK_EXECUTABLE ${CMAKE_CXX_LINK_EXECUTABLE})
endif()

file(GLOB EMPI_CPU_SOURCES src/*.cpp)
file(GLOB EMPI_GPU_SOURCES src/*.cu)
file(GLOB EMPI_TEST_SOURCES tests/*.cpp)

add_library(empi-cpu STATIC ${EMPI_CPU_SOURCES} vendor/sqlite3.c)
target_include_directories(empi-cpu PUBLIC include vendor/include)
set_property(SOURCE vendor/sqlite3.c PROPERTY COMPILE_FLAGS -DSQLITE_OMIT_LOAD_EXTENSION)

add_library(empi-fake-gpu STATIC src/special/alloc.cpp)
target_include_directories(empi-fake-gpu PUBLIC include)

if(UNIX)
    # on UNIX, you can compile statically with cmake -DSTANDALONE=1
    set(STANDALONE FALSE CACHE BOOL "Build a standalone binary?")
elseif(STANDALONE)
    message(SEND_ERROR "Standalone compilation is only available on UNIX")
endif()

if(STANDALONE AND WITH_CUDA)
    message(SEND_ERROR "Standalone compilation cannot use CUDA")
endif()

if(WITH_CUDA)
    add_library(empi-gpu STATIC ${EMPI_GPU_SOURCES} src/special/alloc.cu)
    target_include_directories(empi-gpu PUBLIC include)
    set_property(SOURCE src/SpectrogramCalculatorCUDACallback.cu PROPERTY COMPILE_FLAGS --relocatable-device-code=true)
elseif(STANDALONE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--whole-archive -lrt -lpthread -Wl,--no-whole-archive")
else()
    set(CMAKE_EXE_LINKER_FLAGS "-pthread")
endif()

add_executable(empi src/special/empi.cpp)
install(TARGETS empi)
target_link_libraries(empi empi-cpu fftw3)
if(WITH_CUDA)
    target_link_libraries(empi empi-gpu cufft_static culibos)
else()
    target_link_libraries(empi empi-fake-gpu)
endif()

foreach(TEST_SOURCE ${EMPI_TEST_SOURCES})
    get_filename_component(TEST ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST} tests/${TEST}.cpp)
    target_link_libraries(${TEST} empi-cpu)
    add_test(${TEST} ${TEST})
endforeach()

target_link_libraries(test-best-match empi-fake-gpu fftw3)

target_link_libraries(test-block empi-fake-gpu)

target_link_libraries(test-block-atom empi-fake-gpu)

if(WITH_CUDA)
    add_executable(test-cuda tests/special/test-cuda.cpp)
    target_link_libraries(test-cuda empi-cpu empi-gpu cufft_static culibos)
    add_test(test-cuda test-cuda)
endif()

target_link_libraries(test-dictionary empi-fake-gpu)

target_link_libraries(test-fftw fftw3)

target_link_libraries(test-measure empi-fake-gpu fftw3)

target_link_libraries(test-subsample empi-fake-gpu fftw3)
