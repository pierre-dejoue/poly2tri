cmake_minimum_required(VERSION 3.12)

project(poly2tri LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

option(P2T_BUILD_TESTS "Build tests" OFF)
option(P2T_BUILD_TESTBED "Build the testbed application" OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/third_party/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/config/")

add_subdirectory(poly2tri)

if(P2T_BUILD_TESTS)
    enable_testing()
    add_subdirectory(unittest)
endif()

if(P2T_BUILD_TESTBED)
    add_subdirectory(testbed)
endif()

install(EXPORT poly2tri_project DESTINATION lib/poly2tri)
install(FILES ${PROJECT_SOURCE_DIR}/AUTHORS DESTINATION .)
install(FILES ${PROJECT_SOURCE_DIR}/LICENSE DESTINATION .)
