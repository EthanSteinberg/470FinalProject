cmake_minimum_required(VERSION 3.0)
project(networkprob)

add_executable(main src/main)
add_executable(tests src/test)
add_library(networkprob SHARED src/matlabffi.cpp)

target_compile_options(main PUBLIC -std=c++14 -Wall -Wextra -O0 -g -march=native)
target_compile_options(tests PUBLIC -std=c++14 -Wall -Wextra -O0 -g -march=native)
target_compile_options(networkprob PUBLIC -std=c++14 -Wall -Wextra -O0 -g -march=native)

target_include_directories(main PUBLIC src)
target_include_directories(tests PUBLIC src)
target_include_directories(networkprob PUBLIC src)