enable_testing()

# Dependencies
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.69.0 REQUIRED COMPONENTS
    filesystem
    unit_test_framework
)

if (WIN32)
  # disable autolinking in boost
  add_definitions(-DBOOST_ALL_NO_LIB)
endif()

# Build Unit Tests
add_executable(test_poly2tri
    main.cpp
    TriangleTest.cpp
)

target_include_directories(test_poly2tri
    PRIVATE
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(test_poly2tri
    PRIVATE
    poly2tri
    ${Boost_LIBRARIES}
)

add_test(NAME poly2tri COMMAND test_poly2tri)
