#ifndef WC_IFF_LOADER_WC3DECODER_H
#define WC_IFF_LOADER_WC3DECODER_H

#include "wc_iff_loader/IFFReader.h"
#include <vector>

namespace wc_iff_loader {

struct Vertex { float x, y, z; };
struct Face { uint16_t v1, v2, v3; };

struct Model {
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
};

class WC3Decoder {
public:
    // Decode WC3 model from root FORM chunk
    static Model decode(const IffChunk& root);
};

} // namespace wc_iff_loader

#endif // WC_IFF_LOADER_WC3DECODER_H
