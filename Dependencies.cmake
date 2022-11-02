cmake_minimum_required(VERSION 3.16)

include(cmake/3rdpartyInstall.cmake)

build_thirdparty_project(
  NAME            yaml
  PATH            3rdparty/libyaml
  GIT_REPOSITORY  https://github.com/yaml/libyaml.git
  GIT_TAG         master
)

if (RADEXI_BUILD_TESTS)
  build_thirdparty_project(
    NAME            cmocka
    PATH            3rdparty/cmocka
    GIT_REPOSITORY  https://gitlab.com/cmocka/cmocka.git
    GIT_TAG         master
  )
endif()


