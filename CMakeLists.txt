cmake_minimum_required(VERSION 3.16)
project(rovaca VERSION 1.0.0) # Change version number here when releasing
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.27")
    cmake_policy(SET CMP0003 NEW)
    cmake_policy(SET CMP0144 NEW)
endif()

include(GNUInstallDirs)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(DEFINED LIBDEFLATE_PATH)
    message(STATUS "Using custom libdeflate path")
    link_directories(${CMAKE_SOURCE_DIR})
    find_library(LIBDEFLATE libdeflate.so.0
                PATHS ${CMAKE_SOURCE_DIR}
                NO_DEFAULT_PATH)
else()
    message(STATUS "Using system libdeflate path")
    find_library(LIBDEFLATE libdeflate.so.0 
                PATHS /usr/local/lib)
endif()

if(NOT LIBDEFLATE)
    message(FATAL_ERROR "libdeflate not found")
endif()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN")
# SET(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)


set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

if(CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    set(CMAKE_C_FLAGS "-Wall -Wextra -O0 -g -fPIC")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -O0 -g -fPIC")
else()
    set(CMAKE_C_FLAGS "-Wall -Wextra -O3  -g -fPIC")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3  -g -fPIC")
endif()

# set(MIN_BOOST_VERSION 1.69)
# find_package(Boost ${MIN_BOOST_VERSION} REQUIRED)
# if(Boost_FOUND)
#     message(STATUS "* * * * * * * * * * * boost version: ${Boost_LIB_VERSION} * * * * * * * * * * *")
# endif()

option(ENABLE_GCOV "Enable Gcov" OFF)
option(ENABLE_ASAN "Enable Address Sanitizer" OFF)
option(ENABLE_GOOGLE_TEST "Enable Google Test" OFF)

if(ENABLE_GCOV)
    include(gcov)
    message(STATUS "* * * * * * * * * * * Gcov Enabled * * * * * * * * * * *")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
endif()

if(ENABLE_ASAN)
    include(asan)
    message(STATUS "* * * * * * * * * * * Address Sanitizer Enabled * * * * * * * * * * *")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
endif()

if(ENABLE_GOOGLE_TEST)
    message(STATUS "* * * * * * * * * * * Google Test Enabled * * * * * * * * * * *")
    set(GTest_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_lib/googletest-1.14.0/lib64/cmake/GTest)
    find_package(GTest REQUIRED COMPONENTS)
    if(GTest_FOUND)
        message(STATUS "* * * * * * * * * * * google test version: ${GTest_VERSION} * * * * * * * * * * *")
    endif()
    enable_testing()
endif()

set(HTSLIB_DIR ${CMAKE_SOURCE_DIR}/third_lib/htslib)
set(HTSLIB_INCLUDE_DIR ${HTSLIB_DIR}/include)
set(HTSLIB_LIBRARY_DIR ${HTSLIB_DIR}/lib)

set(SPDLOG_DIR ${CMAKE_SOURCE_DIR}/third_lib/spdlog)
set(SPDLOG_INCLUDE_DIR ${SPDLOG_DIR}/include)
set(SPDLOG_LIBRARY_DIR ${SPDLOG_DIR}/lib64)

if(NOT EXISTS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhts.so.3)
    file(COPY ${HTSLIB_LIBRARY_DIR}/libhts.so.1.18-1-g6dbcddf DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
         FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
    file(RENAME ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhts.so.1.18-1-g6dbcddf ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhts.so.3)
endif()

add_library(hts SHARED IMPORTED)
set_target_properties(hts PROPERTIES IMPORTED_LOCATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhts.so.3)

add_subdirectory(src)

install(FILES "${LIBDEFLATE}" DESTINATION lib)
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhts.so.3 DESTINATION lib)
