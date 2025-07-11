cmake_minimum_required(VERSION 3.10)

project(MyProject C)
add_compile_options(-Wall -Werror -Wextra -g -O0)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(LIBM m)
set(EXE_NAME song_recognition)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PC_FFTW3 REQUIRED fftw3)

file(GLOB SRCS PUBLIC src/*.c)
add_executable(${EXE_NAME} main.c ${SRCS})

# project headers
target_include_directories(${EXE_NAME} PUBLIC include)

# fftw3 headers
target_include_directories(${EXE_NAME} PUBLIC ${PC_FFTW3_INCLUDE_DIRS})

# add fftw3 lib folder
target_link_directories(${EXE_NAME} PUBLIC ${PC_FFTW3_LIBRARY_DIRS})

# link external libs
target_link_libraries(${EXE_NAME} PUBLIC ${PC_FFTW3_LIBRARIES} ${LIBM})
