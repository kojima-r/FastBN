cmake_minimum_required(VERSION 3.10)

project(FastBN_GPU)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(fast_bn
  fast_bn_main.cpp
  fast_bn_dag.cpp
  fast_bn_score.cpp
)
target_compile_options(fast_bn PRIVATE 
  -O3
)


target_compile_options(fast_bn PRIVATE -pg)
target_link_options(fast_bn PRIVATE -pg)