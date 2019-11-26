cmake_minimum_required(VERSION 3.15)
project(DDSketch)

set(CMAKE_CXX_STANDARD 14)

add_executable(DDSKetch main.cpp ddsketch.cc ddsketch.h error.cpp error.h)