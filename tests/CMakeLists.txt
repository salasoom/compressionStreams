add_executable(testBasics testBasics.cpp)
target_link_libraries(testBasics compressedStream)
add_test(NAME testBasics COMMAND testBasics)
