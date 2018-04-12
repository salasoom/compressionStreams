add_library(compressedStream compressedStream.cpp)

target_include_directories(compressedStream
  PUBLIC ${CMAKE_SOURCE_DIR}/include
)
