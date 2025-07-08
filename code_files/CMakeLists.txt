# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Define the project name and specify the language (C)
project(MyProject C)

# Find the pkg-config utility, which helps locate libraries
find_package(PkgConfig REQUIRED)

# Use pkg-config to find the FFTW3 library. This populates the PC_FFTW3_* variables.
pkg_check_modules(PC_FFTW3 REQUIRED fftw3)

# Add compiler flags that will be used for all targets
add_compile_options(-Wall -Wextra -g)

# Set the C standard for the project
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Define the executable target and list its source files
add_executable(main main.c func.c)

# Add the FFTW include directories to our target.
# This is how the compiler finds "fftw3.h".
target_include_directories(main PUBLIC ${PC_FFTW3_INCLUDE_DIRS})

# --- ADDED THIS LINE ---
# Add the FFTW library search path to our target.
# This tells the linker where to find the actual library file (e.g., libfftw3.a).
target_link_directories(main PUBLIC ${PC_FFTW3_LIBRARY_DIRS})

# Link the required libraries to our target by name.
# This links the FFTW library (-lfftw3) and the standard math library (m).
target_link_libraries(main PUBLIC ${PC_FFTW3_LIBRARIES} m)
