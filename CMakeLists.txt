cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)

# OpenGL Dependencies
find_package(OpenGL REQUIRED)
# GLFW Dependencies
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_SOURCE_DIR}/glfw)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable( 
    Main
        src/Main.cpp
        src/PointSet.cpp
        src/glad.c
        src/imgui/imgui.cpp
        src/imgui/imgui_draw.cpp
        src/imgui/imgui_impl_glfw_gl3.cpp
)
target_link_libraries(
    Main
        glfw 
        ${OPENGL_gl_LIBRARY}
 )