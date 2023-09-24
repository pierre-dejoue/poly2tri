# Dependencies
find_package(OpenGL REQUIRED)
include(glfw)

# Build testbed
add_executable(p2t
    main.cpp
)

target_link_libraries(p2t
    PRIVATE
    glfw
    OpenGL::GL
    poly2tri
)
