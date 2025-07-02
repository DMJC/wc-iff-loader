# wc-iff-loader

A small C++ library for decoding EA/Origin IFF files as used by Wing Commander 3.
The library provides a generic IFF parser and a simple decoder that extracts
3D model data (vertices and faces) from WC3 model files.

## Building

This project uses CMake:

```bash
mkdir build && cd build
cmake ..
make
ctest
```

The example test generates a minimal WC3 style IFF file and verifies that the
library decodes it correctly.
