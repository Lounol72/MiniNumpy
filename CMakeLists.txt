cmake_minimum_required(VERSION 4.0)
project(MiniNumpy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB SRCS CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB EXECS CONFIGURE_DEPENDS "apps/*.cpp")

add_library(mnp STATIC ${SRCS})
target_include_directories(mnp PUBLIC include)

foreach(FICHIER IN LISTS EXECS)
  get_filename_component(EXEC ${FICHIER} NAME_WE)
  add_executable(${EXEC} ${FICHIER})
  message(STATUS "fichier=${FICHIER}  exec=${EXEC}")
  target_link_libraries(${EXEC} PRIVATE mnp)
endforeach()
