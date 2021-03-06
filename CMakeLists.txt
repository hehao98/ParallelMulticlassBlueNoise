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
add_subdirectory(${CMAKE_SOURCE_DIR}/RainbowNoiseSpectrum)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/RainbowNoiseLibrary)
include_directories(${CMAKE_SOURCE_DIR}/RainbowNoiseSpectrum)

# Organize source file
set(IMGUI_SRC
        src/imgui/imgui.cpp
        src/imgui/imgui_draw.cpp
        src/imgui/imgui_impl_glfw_gl3.cpp
        )

set(LIBS ${LIBS} glfw ${OPENGL_gl_LIBRARY} RainbowNoiseLibrary RainbowNoiseSpectrum)

add_executable(
        Main src/Main.cpp src/PointSet.cpp src/SpectrumImage.cpp src/glad.c ${IMGUI_SRC}
)
target_link_libraries(Main ${LIBS})

add_executable(DartThrowing src/DartThrowing.cpp)
target_link_libraries(DartThrowing ${LIBS})

add_executable(ParallelDartThrowing src/ParallelDartThrowing.cpp)
target_link_libraries(ParallelDartThrowing ${LIBS})

add_executable(ExDartThrowing src/ExperimentalDartThrowing.cpp)
target_link_libraries(ExDartThrowing ${LIBS})

add_executable(SFT src/SFT.cpp)
target_link_libraries(SFT ${LIBS})

add_executable(PFM2PPM src/PFM2PPM.cpp)
target_link_libraries(PFM2PPM ${LIBS})