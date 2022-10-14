cmake_minimum_required(VERSION 3.16)

project(radexi LANGUAGES C)

option(RADEXI_BUILD_TESTS "Build and execute tests" OFF)
option(RADEXI_BUILD_3RDPARTY "Build 3rdparty utilities" ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_compile_options(
  -Wall
  -Wextra
  -pedantic
  -Wno-format-truncation
  $<$<CONFIG:DEBUG>:-g3>
  $<$<CONFIG:DEBUG>:-Og>
  $<$<CONFIG:RELEASE>:-O3>
)

add_compile_definitions(
  #  $<$<CONFIG:DEBUG>:USE_SMTH>
)

configure_file(CMakeLists.txt.in 3rdparty/libyaml/CMakeLists.txt)
execute_process(
  COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/3rdparty/libyaml
)
execute_process(
  COMMAND "${CMAKE_COMMAND}" --build .
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/3rdparty/libyaml
)

add_subdirectory(3rdparty/libyaml)

set(SOURCES
  ${CMAKE_SOURCE_DIR}/src/input/range.c
  ${CMAKE_SOURCE_DIR}/src/lamda/doc.c
  ${CMAKE_SOURCE_DIR}/src/lamda/node.c
  ${CMAKE_SOURCE_DIR}/src/lamda/parse.c
)

set(INCLUDES
  ${CMAKE_SOURCE_DIR}/include
)

set(TESTS
  ${CMAKE_SOURCE_DIR}/test/range_test.c
  ${CMAKE_SOURCE_DIR}/test/lamda_test.c
)

include_directories(${INCLUDES})

add_library(radexi SHARED ${SOURCES})

add_executable(radexi-cmd
  ${CMAKE_SOURCE_DIR}/src/main.c
  ${SOURCES}
)

if(RADEXI_BUILD_TESTS)
  configure_file(CMakeLists.txt.in 3rdparty/cmocka/CMakeLists.txt)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/3rdparty/cmocka
  )
  execute_process(
    COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/3rdparty/cmocka
  )

  add_subdirectory(3rdparty/cmocka)

  add_subdirectory(test)
endif()
