cmake_minimum_required(VERSION 3.10)
project(MyProject C)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PC_FFTW3 REQUIRED fftw3)

add_compile_options(-Wall -Wextra -g)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# 1. Update the source paths to include the "src/" prefix
add_executable(main 
    src/main.c 
    src/func.c
    src/database.c
)

# 2. Add your local "include" folder to the search path
# This allows you to use #include "func.h" without relative paths
target_include_directories(main PUBLIC 
    include
    ${PC_FFTW3_INCLUDE_DIRS}
)

target_link_directories(main PUBLIC ${PC_FFTW3_LIBRARY_DIRS})
target_link_libraries(main PUBLIC ${PC_FFTW3_LIBRARIES} m)