cmake_minimum_required(VERSION 3.1)

project(ParallelMulticlassBlueNoise)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")

# OpenGL Dependencies
find_package(OpenGL REQUIRED)

# GLFW Dependencies
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_SOURCE_DIR}/glfw)

# Rainbow Noise Library(From Li-Yi Wei's SIGGRAPH code)
add_subdirectory(${CMAKE_SOURCE_DIR}/RainbowNoiseLibrary)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/RainbowNoiseLibrary)

# Organize source file
set(IMGUI_SRC
        src/imgui/imgui.cpp
        src/imgui/imgui_draw.cpp
        src/imgui/imgui_impl_glfw_gl3.cpp
        )

set(LIBS ${LIBS} glfw ${OPENGL_gl_LIBRARY} RainbowNoiseLibrary)

add_executable(
        Main src/Main.cpp src/PointSet.cpp src/glad.c ${IMGUI_SRC}
)
target_link_libraries(Main ${LIBS})

add_executable(DartThrowing src/DartThrowing.cpp)
target_link_libraries(DartThrowing ${LIBS})