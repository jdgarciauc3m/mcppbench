add_executable(blackscholes-seq-orig main.cpp)

add_test(
    NAME blackscholes-seq-orig-small
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/blackscholes-seq-orig 1
      ${CMAKE_CURRENT_SOURCE_DIR}/../input/in_4.txt
      ${CMAKE_CURRENT_BINARY_DIR}/out_4.txt
)

add_test(
    NAME blackscholes-seq-orig-large
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/blackscholes-seq-orig 1
    ${CMAKE_CURRENT_SOURCE_DIR}/../input/in_10M.txt
    ${CMAKE_CURRENT_BINARY_DIR}/out_10M.txt
)
