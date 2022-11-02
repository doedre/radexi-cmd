cmake_minimum_required(VERSION 3.16)

project(radexi LANGUAGES C)

option(RADEXI_BUILD_TESTS "Build and execute tests" OFF)

# Load dependencies
include(Dependencies.cmake)

find_package(yaml REQUIRED)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(RADEXI_SOURCES
  ${CMAKE_SOURCE_DIR}/src/input/range.c
  ${CMAKE_SOURCE_DIR}/src/lamda/doc.c
  ${CMAKE_SOURCE_DIR}/src/lamda/node.c
  ${CMAKE_SOURCE_DIR}/src/lamda/parse.c
)

add_library(radexi SHARED ${RADEXI_SOURCES})

target_include_directories(radexi
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_compile_options(radexi
  PRIVATE
    -Wall
    -Wextra
    -pedantic
    $<$<CONFIG:DEBUG>:-g3>
    $<$<CONFIG:DEBUG>:-Og>
    $<$<CONFIG:RELEASE>:-O3>
)

target_compile_definitions(radexi
  PRIVATE
    #  $<$<CONFIG:DEBUG>:USE_SMTH>
)

target_compile_features(radexi
  PRIVATE
    c_std_99
    c_static_assert
    c_variadic_macros
)

if (RADEXI_BUILD_TESTS)
  find_package(cmocka REQUIRED)
  add_subdirectory(test)
endif()

