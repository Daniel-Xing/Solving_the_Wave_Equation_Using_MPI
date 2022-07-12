cmake_minimum_required(VERSION 3.0.0)
project(MPI-CourseWork VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
    set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")
    set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++17")
elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1)
    set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++1z")
    set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++1z")
endif()

# add_subdirectory(Tests)
#add_compile_options(-D_Float128=__float128)
# target_compile_features(${TARGET_NAME} PRIVATE cxx_std_17)

find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(MPI-CourseWork main.cpp)

include_directories("src")

target_sources(MPI-CourseWork
PRIVATE
    "src/Mpi_CourseWork.h"
    "src/Mpi_CourseWork.cpp"
    "src/Serial_Wave_Equation.h"
    "src/Serial_Wave_Equation.cpp"
)

find_package(MPI REQUIRED)
target_link_libraries(MPI-CourseWork PUBLIC MPI::MPI_CXX)

#target_link_libraries(MPI-CourseWork ${MPI_INCLUDE_PATH})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)