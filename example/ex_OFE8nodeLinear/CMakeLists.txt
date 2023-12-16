cmake_minimum_required(VERSION 3.10)
project(ex_OFE8nodeLinear_GEVP VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the default build type to Release if it is not already set
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find the Threads package
find_package(Threads REQUIRED)

# Source files
file(GLOB_RECURSE SRC_FILES ../../src/OFE8nodeLinear.cpp ../../src/ConstitutiveLaw.cpp)
file(GLOB_RECURSE OTHER_FILES ../../src/SkyLine.cpp)

# Add executable
add_executable(main main.cpp ${SRC_FILES} ${OTHER_FILES})

# Include directories
target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_link_libraries(main Threads::Threads)

# Print the build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")