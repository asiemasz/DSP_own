cmake_minimum_required(VERSION 3.16)

file(GLOB C_SOURCES
    "src/*.c"
)

add_library(DSP_own STATIC ${C_SOURCES})
target_include_directories(DSP_own PUBLIC src)