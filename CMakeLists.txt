#
# GPU 版の切り替えは
# -DCMAKE_C_COMPILER=nvc -DCMAKE_CXX_COMPILER=nvc++
# で行う
#
# CUDA のバージョンを環境変数 CUDA_VERSION に記述する（例：12.6, 12.9, 13.1）
#
cmake_minimum_required(VERSION 3.10)

project(FastBN_GPU LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS ${CMAKE_CXX_COMPILER_ID})

execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_executable(fast_bn
        fast_bn_main.cpp
        fast_bn_score.cpp
        fast_bn_counts.cpp
    )
    target_compile_definitions(fast_bn PRIVATE GIT_HASH="${GIT_HASH}")
    target_compile_options(fast_bn PRIVATE 
        -O3
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    find_package(OpenACC REQUIRED)
    add_executable(fast_bn
        fast_bn_main.cpp
        fast_bn_score.cpp
        fast_bn_counts.cpp
        fast_bn_utils.cpp
    )
    target_compile_definitions(fast_bn PRIVATE GIT_HASH="${GIT_HASH}")
    if(OpenACC_CXX_FOUND)
        target_link_libraries(fast_bn PRIVATE OpenACC::OpenACC_CXX)
    endif()
    target_compile_options(fast_bn PRIVATE
        -acc
        -target=gpu
        -gpu=ccnative,cuda$ENV{CUDA_VERSION}
        -Minfo=accel
        -fast
        -cuda
    )
    target_link_options(fast_bn PRIVATE 
        -acc
        -cuda
    )
endif()

#target_compile_options(fast_bn PRIVATE -pg)
#target_link_options(fast_bn PRIVATE -pg)