add_executable(radexi_tests
  main.c
)

target_include_directories(radexi_tests
  PRIVATE
  ${INCLUDES}
  ${CMAKE_SOURCE_DIR}/3rdparty/cmocka/include
)

target_link_libraries(radexi_tests
  radexi
  cmocka
)

add_custom_command(
  TARGET radexi_tests POST_BUILD
  COMMAND radexi_tests
  VERBATIM
)
