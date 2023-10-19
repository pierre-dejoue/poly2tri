# Dependencies
if (WIN32)
    set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost 1.69 REQUIRED COMPONENTS
    filesystem
    unit_test_framework
)

#
# Build unit tests of the library
#
add_executable(test_poly2tri
    test_cdt.cpp
    test_triangle.cpp
    utility.cpp
)

target_include_directories(test_poly2tri
    PRIVATE
    ${Boost_INCLUDE_DIRS}
)

target_compile_definitions(test_poly2tri
    PRIVATE
    P2T_BASE_DIR="${PROJECT_SOURCE_DIR}"
)

target_link_libraries(test_poly2tri
    PRIVATE
    poly2tri
    ${Boost_LIBRARIES}
)

add_test(NAME poly2tri COMMAND test_poly2tri)

#
# Build unit tests of the geometric predicates
#
add_executable(test_predicates
    test_predicates.cpp
)

target_include_directories(test_predicates
    PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${poly2tri_SOURCE_DIR}/poly2tri/include
    ${poly2tri_SOURCE_DIR}/poly2tri/src
)

# The geometric predicates are header-only, hence the test suite does not need to be linked with the poly2tri library
target_link_libraries(test_predicates
    PRIVATE
    ${Boost_LIBRARIES}
)

add_test(NAME predicates COMMAND test_predicates)
