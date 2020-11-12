cmake_minimum_required(VERSION 3.4)

add_library(iMath Vector.h Matrix.h Matrix.cpp)
add_executable(main main.cpp Material.h Material.cpp )
target_link_libraries(main iMath  )
