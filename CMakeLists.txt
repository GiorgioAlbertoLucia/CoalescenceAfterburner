cmake_minimum_required(VERSION 3.14)
project(li4_coalescence VERSION 1.0 LANGUAGES CXX)

# ── C++ standard ─────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Build type default ────────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ── Find ROOT via root-config ─────────────────────────────────────────────────
# This approach works regardless of how ROOT was installed (CVMFS, conda,
# manual build) as long as root-config is in PATH or ROOTSYS is set.

find_program(ROOT_CONFIG root-config
    HINTS $ENV{ROOTSYS}/bin
    REQUIRED
)

# Retrieve flags from root-config
execute_process(COMMAND ${ROOT_CONFIG} --version
    OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${ROOT_CONFIG} --incdir
    OUTPUT_VARIABLE ROOT_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${ROOT_CONFIG} --libs
    OUTPUT_VARIABLE ROOT_LIBS_RAW OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${ROOT_CONFIG} --cflags
    OUTPUT_VARIABLE ROOT_CXX_FLAGS_RAW OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${ROOT_CONFIG} --libdir
    OUTPUT_VARIABLE ROOT_LIBRARY_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

# Convert space-separated flag strings into CMake lists
separate_arguments(ROOT_LIBS     UNIX_COMMAND ${ROOT_LIBS_RAW})
separate_arguments(ROOT_CXX_FLAGS UNIX_COMMAND ${ROOT_CXX_FLAGS_RAW})

message(STATUS "ROOT version : ${ROOT_VERSION}")
message(STATUS "ROOT incdir  : ${ROOT_INCLUDE_DIRS}")
message(STATUS "ROOT libdir  : ${ROOT_LIBRARY_DIR}")
message(STATUS "ROOT libs    : ${ROOT_LIBS}")

# ── Source files ──────────────────────────────────────────────────────────────
set(SOURCES
    src/EventReader.cxx
    src/CoalescenceEngine.cxx
    src/ToyMCEngine.cxx
    src/CoalescenceEngineLi4.cxx
    src/CoalescenceEngineHe4.cxx
    src/CoalescenceEngineHe3.cxx
    src/IntegratorA4.cxx
    src/IntegratorA3.cxx
    src/ToyMCEngineA4.cxx
    src/ToyMCEngineA3.cxx
)

# ── Include directories ───────────────────────────────────────────────────────
set(HEADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ── Library (optional: makes unit testing easier later) ───────────────────────
add_library(li4_core STATIC ${SOURCES})

target_include_directories(li4_core
    PUBLIC
        ${HEADERS_DIR}
        ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(li4_core
    PUBLIC
        ${ROOT_LIBS}
)

target_compile_options(li4_core
    PRIVATE
        ${ROOT_CXX_FLAGS}
        -Wall
        -Wextra
        -Wpedantic
        $<$<CONFIG:Release>:-O2>
        $<$<CONFIG:Debug>:-g -O0>
)

# ── Executable ────────────────────────────────────────────────────────────────
add_executable(li4_coalescence main.cxx)

target_include_directories(li4_coalescence
    PRIVATE
        ${HEADERS_DIR}
        ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(li4_coalescence
    PRIVATE
        li4_core
        ${ROOT_LIBS}
)

# ── ToyMC executable ─────────────────────────────────────────────────────────
add_executable(li4_toymc main_toymc.cxx)

target_include_directories(li4_toymc
    PRIVATE
        ${HEADERS_DIR}
        ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(li4_toymc
    PRIVATE
        li4_core
)

# ── A=4 integrator executable ─────────────────────────────────────────────────
add_executable(li4_a4 main_a4.cxx)


target_include_directories(li4_a4 PRIVATE ${HEADERS_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(li4_a4 PRIVATE li4_core)

# ── A=3 integrator executable ─────────────────────────────────────────────────
add_executable(he3_a3 main_a3.cxx)

target_include_directories(he3_a3 PRIVATE ${HEADERS_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(he3_a3 PRIVATE li4_core)

# ── Install ───────────────────────────────────────────────────────────────────
install(TARGETS li4_coalescence DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# ── Summary ───────────────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
message(STATUS " li4_coalescence configuration")
message(STATUS "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
message(STATUS " Install prefix : ${CMAKE_INSTALL_PREFIX}")
message(STATUS " Build type     : ${CMAKE_BUILD_TYPE}")
message(STATUS " C++ standard   : ${CMAKE_CXX_STANDARD}")
message(STATUS " ROOT version   : ${ROOT_VERSION}")
message(STATUS "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
message(STATUS "")
