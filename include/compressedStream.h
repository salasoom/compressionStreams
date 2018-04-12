#ifndef COMPRESSED_STREAM_H
#define COMPRESSED_STREAM_H

#include <sstream>
using namespace std;

class CompressedStream : stringstream {
public:
  CompressedStream() {};
};

#endif
