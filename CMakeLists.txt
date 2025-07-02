cmake_minimum_required(VERSION 3.10)
project(wc_iff_loader CXX)

add_library(wc_iff
    src/IFFReader.cpp
    src/WC3Decoder.cpp
)

target_include_directories(wc_iff PUBLIC include)
set_target_properties(wc_iff PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

add_executable(test_wc_iff test/test_wc_iff.cpp)
target_link_libraries(test_wc_iff wc_iff)

enable_testing()
add_test(NAME test_wc_iff COMMAND test_wc_iff)
