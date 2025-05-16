cmake_minimum_required(VERSION 3.20)
project(FluidSim VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(glad REQUIRED)

# Find GLEW
find_package(GLEW REQUIRED)

# Find GLFW
find_package(glfw3 REQUIRED)

# Find GLM
find_package(glm REQUIRED)

add_executable(FluidSim src/main.cpp)

target_link_directories(FluidSim PRIVATE
    "C:/vcpkg/installed/x64-windows/debug/lib" # Explicitly add the debug library path
)

target_link_libraries(FluidSim PRIVATE
    GLEW::GLEW
    glfw3dll # Use the namespaced target for linking
    glm::glm
    glad
    opengl32 # Add explicit OpenGL linking
    glu32    # Add explicit GLU linking
)