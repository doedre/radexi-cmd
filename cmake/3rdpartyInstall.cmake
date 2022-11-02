cmake_minimum_required(VERSION 3.16)

function(build_thirdparty_project)
  set(ONE_VALUE_KEYWORDS NAME PATH GIT_REPOSITORY GIT_TAG)
  set(MULTI_VALUE_KEYWORDS FLAGS)
  cmake_parse_arguments(_BUILD "" "${ONE_VALUE_KEYWORDS}" "${MULTI_VALUE_KEYWORDS}" ${ARGN})
  
  if (NOT DEFINED _BUILD_PATH)
      message(FATAL_ERROR "PATH is a required keyword to `build_thirdparty_project`.")
  endif()

  if (NOT DEFINED _BUILD_NAME)
      message(FATAL_ERROR "NAME is a required keyword to 'build_thirdparty_project`.")
  endif()

  include(FetchContent)
  FetchContent_Declare(${_BUILD_NAME}
    GIT_REPOSITORY "${_BUILD_GIT_REPOSITORY}"
    GIT_TAG "${_BUILD_GIT_TAG}"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/${_BUILD_PATH}"
    OVERRIDE_FIND_PACKAGE
  )
endfunction()
