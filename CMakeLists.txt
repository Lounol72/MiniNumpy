cmake_minimum_required(VERSION 4.0)
project(MiniNumpy LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Type de build : Debug ou Release" FORCE)
endif()

enable_testing()
# OPTIONS de compilations
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(
  -Wall -Wextra -Werror -Wpedantic
  "$<$<CONFIG:Debug>:-fsanitize=address,undefined>"
)
add_link_options(
  "$<$<CONFIG:Debug>:-fsanitize=address,undefined>"
)

#FICHIERS
file(GLOB SRCS CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB EXECS CONFIGURE_DEPENDS "apps/*.cpp")
file(GLOB TESTS CONFIGURE_DEPENDS "tests/*.cpp")

add_library(mnp STATIC ${SRCS})
target_include_directories(mnp PUBLIC include)

add_library(doctest INTERFACE)
target_include_directories(doctest INTERFACE third_party)

foreach(FICHIER IN LISTS EXECS)
  get_filename_component(EXEC ${FICHIER} NAME_WE)
  add_executable(${EXEC} ${FICHIER})
  message(STATUS "fichier=${FICHIER}  exec=${EXEC}")
  target_link_libraries(${EXEC} PRIVATE mnp)
endforeach()

foreach(FICHIER IN LISTS TESTS)
  get_filename_component(EXEC ${FICHIER} NAME_WE)
  add_executable(${EXEC} ${FICHIER})
  message(STATUS "fichier=${FICHIER}  exec=${EXEC}")
  target_link_libraries(${EXEC} PRIVATE mnp)
  target_link_libraries(${EXEC} PRIVATE doctest)
  add_test(NAME ${EXEC} COMMAND ${EXEC})
endforeach()

#==========================================================
# Macros
#==========================================================
macro(create_shell_task TARGET_NAME)
  add_custom_target(${TARGET_NAME}
    COMMAND ${ARGN}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
COMMENT "Exécution de la tâche : ${TARGET_NAME}"
VERBATIM
  )
endmacro()
