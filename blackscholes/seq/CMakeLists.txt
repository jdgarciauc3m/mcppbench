add_compile_options(-Wall -Wextra -Wpedantic -Werror)
add_compile_options(
    -Wcast-qual
    -Wconversion-null
    -Wmissing-declarations
    -Woverlength-strings
    -Wpointer-arith
    -Wunused-local-typedefs
    -Wunused-result
    -Wvarargs
    -Wvla
    -Wwrite-strings)


add_executable(blackscholes-seq main.cpp)

add_test(
    NAME blackscholes-seq-small
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/blackscholes-seq 1
    ${CMAKE_CURRENT_SOURCE_DIR}/../input/in_4.txt
    ${CMAKE_CURRENT_BINARY_DIR}/out_4.txt
)

add_test(
    NAME blackscholes-seq-regresion
    COMMAND ${CMAKE_COMMAND} -E compare_files
    ${CMAKE_CURRENT_BINARY_DIR}/out_4.txt
    ${CMAKE_CURRENT_BINARY_DIR}/../seq-orig/out_4.txt
)

add_test(
    NAME blackscholes-seq-large
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/blackscholes-seq 1
    ${CMAKE_CURRENT_SOURCE_DIR}/../input/in_10M.txt
    ${CMAKE_CURRENT_BINARY_DIR}/out_10M.txt
)