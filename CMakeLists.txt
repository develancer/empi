# University of Warsaw, Department of Biomedical Physics
#   Enhanced Matching Pursuit Implementation (empi)
# See README.md and LICENCE for details.

cmake_minimum_required(VERSION 2.6)
project(empi CXX C)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# Basic configuration
file(GLOB EMPI_SOURCES src/*.cpp)
add_executable(empi ${EMPI_SOURCES})
set_property(TARGET empi PROPERTY CXX_STANDARD 11)

# SQLite source code is included
add_library(sqlite src/sqlite3.c)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSQLITE_OMIT_LOAD_EXTENSION")

# You can compile statically with cmake -DSTANDALONE=1
set(STANDALONE FALSE CACHE BOOL "Build a standalone binary?")
if(STANDALONE AND UNIX)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
	set(CMAKE_EXE_LINKER_FLAGS "-static")
	set_property(TARGET empi PROPERTY LINK_SEARCH_END_STATIC 1)
endif()

# FFTW is required
find_path(FFTW_INCLUDES fftw3.h)
find_library(FFTW_LIBRARIES NAMES fftw3 libfftw3-3 fftw)
if(NOT FFTW_LIBRARIES)
	message(FATAL_ERROR "Required FFTW library has not been found!")
endif()
include_directories(inc ${FFTW_INCLUDES})

# OpenMP is optional
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(MSVC)
	# Force to always compile with W4
	if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
		string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
	endif()

	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
	if(NOT STANDALONE)
		# Enable AVX2 optimizations
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /arch:AVX2")
	endif()

elseif(CMAKE_COMPILER_IS_GNUCXX)
	# Add GCC warning level
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
endif()

# linking it all
target_link_libraries(empi ${FFTW_LIBRARIES})
target_link_libraries(empi sqlite)
install(TARGETS empi RUNTIME DESTINATION bin)

if(UNIX)
	add_custom_target(distclean)
	add_custom_command(
		DEPENDS clean
		COMMENT "Cleaning for source distribution"
		COMMAND rm
		ARGS    -rf CMakeCache.txt CMakeFiles CMakeScripts Makefile *.cmake
		TARGET  distclean
	)

	# Cross-compiling
	set(CROSS_CXX_FLAGS --std=c++11 -Wall -Wextra -Wpedantic -DNDEBUG -O2 -I${CMAKE_CURRENT_SOURCE_DIR}/inc)
	set(WIN_LD_FLAGS -static -Wl,-Bdynamic -lfftw3 -Wl,-Bstatic)
	set(OSX_LD_FLAGS -lfftw3)

	add_custom_target(empi-osx64)
	add_custom_command(
		COMMENT "Cross-compiling version for 64-bit OS X"
		COMMAND x86_64-apple-darwin15-clang -c -O3 -o sqlite3.o src/sqlite3.c
		COMMAND x86_64-apple-darwin15-clang++-libc++ ${CROSS_CXX_FLAGS} -o empi-osx64 ${EMPI_SOURCES} sqlite3.o ${OSX_LD_FLAGS}
		TARGET  empi-osx64
	)
	add_custom_target(empi-win32.exe)
	add_custom_command(
		COMMENT "Cross-compiling static version for 32-bit Windows"
		COMMAND i686-w64-mingw32-gcc -c -fopenmp -O3 -o sqlite3.o src/sqlite3.c
		COMMAND i686-w64-mingw32-g++ -fopenmp ${CROSS_CXX_FLAGS} -o empi-win32.exe ${EMPI_SOURCES} sqlite3.o ${WIN_LD_FLAGS}
		TARGET  empi-win32.exe
	)
	add_custom_target(empi-win64.exe)
	add_custom_command(
		COMMENT "Cross-compiling static version for 64-bit Windows"
		COMMAND x86_64-w64-mingw32-gcc -c -fopenmp -O3 -o sqlite3.o src/sqlite3.c
		COMMAND x86_64-w64-mingw32-g++ -fopenmp ${CROSS_CXX_FLAGS} -o empi-win64.exe ${EMPI_SOURCES} sqlite3.o ${WIN_LD_FLAGS}
		TARGET  empi-win64.exe
	)

	add_custom_target(cross)
	add_dependencies(cross empi-osx64 empi-win32.exe empi-win64.exe)
	set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "empi-osx64;empi-win32.exe;empi-win64.exe;sqlite3.o")
endif()
