#ifndef WC_IFF_LOADER_IFFREADER_H
#define WC_IFF_LOADER_IFFREADER_H

#include <cstdint>
#include <string>
#include <vector>
#include <istream>

namespace wc_iff_loader {

struct IffChunk {
    std::string id;                     // Chunk ID (e.g. "FORM", "VERT")
    uint32_t size = 0;                  // Size of chunk data
    std::string type;                   // FORM/LIST type
    std::vector<uint8_t> data;          // Raw data for leaf chunks
    std::vector<IffChunk> children;     // Child chunks for FORM/LIST/CAT
};

class IFFReader {
public:
    // Load an IFF file from disk. Returns true on success.
    bool load(const std::string& path);

    // Access parsed root chunk
    const IffChunk& root() const { return rootChunk; }

private:
    IffChunk rootChunk;
    bool parseChunk(std::istream& is, IffChunk& parent, uint32_t size);
    static uint32_t readU32BE(std::istream& is);
    static float readF32BE(std::istream& is);
};

} // namespace wc_iff_loader

#endif // WC_IFF_LOADER_IFFREADER_H
