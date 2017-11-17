add_executable(test_peano test_peano.c)
target_link_libraries(test_peano peano criterion)
add_test(NAME test COMMAND test_peano)
