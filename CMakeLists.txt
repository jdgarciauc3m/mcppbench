cmake_minimum_required(VERSION 3.16)
project(mcppbench)

set(CMAKE_CXX_STANDARD 17)

find_package(fmt)

enable_testing()

add_subdirectory(blackscholes)