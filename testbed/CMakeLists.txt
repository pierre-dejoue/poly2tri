find_package(OpenGL REQUIRED)

add_executable(testbed
    main.cc
)

target_include_directories(testbed
    PRIVATE
    ../third_party/glfw/include
)

target_link_libraries(testbed
    PRIVATE
    glfw
    OpenGL::GL
    poly2tri
)
