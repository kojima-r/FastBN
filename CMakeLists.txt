#
# GPU 版の切り替えは
# -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++
# で行う
#
cmake_minimum_required(VERSION 3.10)

project(FastBN_GPU)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS ${CMAKE_CXX_COMPILER_ID})

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_executable(fast_bn
        fast_bn_main.cpp
        fast_bn_score.cpp
        fast_bn_dag.cpp
    )
    target_compile_options(fast_bn PRIVATE 
        -O3
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    find_package(OpenACC REQUIRED)
    add_executable(fast_bn
        fast_bn_main.cpp
        fast_bn_score.cpp
        fast_bn_dag.cpp
    )
    if(OpenACC_CXX_FOUND)
        target_link_libraries(fast_bn PRIVATE OpenACC::OpenACC_CXX)
    endif()
    target_compile_options(fast_bn PRIVATE 
        -acc
        -target=gpu
        -gpu=ccnative,cuda$ENV{CUDA_VERSION}
        -Minfo=accel
        -fast
    )
    target_link_options(fast_bn PRIVATE 
        -acc
        -target=gpu
        -gpu=ccnative,cuda$ENV{CUDA_VERSION}
        -Minfo=accel
        -fast
    )
endif()

#target_compile_options(fast_bn PRIVATE -pg)
#target_link_options(fast_bn PRIVATE -pg)