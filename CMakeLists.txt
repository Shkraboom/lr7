cmake_minimum_required(VERSION 3.16)
project(lr7_profiling LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# FFTW3
set(FFTW3_USE_OK OFF)
find_library(FFTW3_LIB NAMES fftw3 libfftw3 fftw3-3
  PATHS /usr/local/lib /opt/homebrew/lib
  HINTS $ENV{HOME}/brew/lib $ENV{FFTW3_ROOT} "$ENV{FFTW3_ROOT}/lib")
find_path(FFTW3_INC NAMES fftw3.h
  PATHS /usr/local/include /opt/homebrew/include
  HINTS $ENV{HOME}/brew/include $ENV{FFTW3_ROOT} "$ENV{FFTW3_ROOT}/include")
if(FFTW3_LIB AND FFTW3_INC)
  set(FFTW3_USE_OK ON)
  set(FFTW3_LIBRARIES ${FFTW3_LIB})
  set(FFTW3_INCLUDE_DIRS ${FFTW3_INC})
endif()
if(NOT FFTW3_USE_OK)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(FFTW3 QUIET fftw3)
    if(FFTW3_FOUND)
      set(FFTW3_USE_OK ON)
    endif()
  endif()
endif()

if(NOT FFTW3_USE_OK)
  message(FATAL_ERROR "FFTW3 not found. Install fftw (e.g. brew install fftw) or set FFTW3_ROOT.")
endif()

add_executable(lr7_benchmark
    main.cpp
    param.cpp
    param.h
    calcfft_original.cpp
    calcfft_original.h
    calcfft_optimized.cpp
    calcfft_optimized.h
)

target_include_directories(lr7_benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${FFTW3_INCLUDE_DIRS}
)
if(FFTW3_LIBRARY_DIRS)
  target_link_directories(lr7_benchmark PRIVATE ${FFTW3_LIBRARY_DIRS})
endif()

target_link_libraries(lr7_benchmark PRIVATE ${FFTW3_LIBRARIES})

target_compile_options(lr7_benchmark PRIVATE -O2)
