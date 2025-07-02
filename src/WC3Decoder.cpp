#include "wc_iff_loader/WC3Decoder.h"
#include <cstring>

namespace wc_iff_loader {

static const IffChunk* findChunk(const IffChunk& parent, const std::string& id) {
    for (const auto& c : parent.children) {
        if (c.id == id)
            return &c;
    }
    return nullptr;
}

Model WC3Decoder::decode(const IffChunk& root) {
    Model model;
    const IffChunk* vert = findChunk(root, "VERT");
    const IffChunk* face = findChunk(root, "FACE");
    if (!vert || !face)
        return model;

    std::istream vstream(nullptr);
    // parse vertices
    {
        const std::vector<uint8_t>& data = vert->data;
        if (data.size() < 4)
            return model;
        size_t offset = 0;
        uint32_t count = (data[offset] << 24) | (data[offset+1] << 16) | (data[offset+2] << 8) | data[offset+3];
        offset += 4;
        model.vertices.resize(count);
        for (uint32_t i=0; i<count; ++i) {
            uint32_t bx = (data[offset] << 24) | (data[offset+1] << 16) | (data[offset+2] << 8) | data[offset+3]; offset += 4;
            uint32_t by = (data[offset] << 24) | (data[offset+1] << 16) | (data[offset+2] << 8) | data[offset+3]; offset += 4;
            uint32_t bz = (data[offset] << 24) | (data[offset+1] << 16) | (data[offset+2] << 8) | data[offset+3]; offset += 4;
            float fx, fy, fz;
            std::memcpy(&fx, &bx, sizeof(float));
            std::memcpy(&fy, &by, sizeof(float));
            std::memcpy(&fz, &bz, sizeof(float));
            model.vertices[i] = {fx, fy, fz};
        }
    }

    // parse faces
    {
        const std::vector<uint8_t>& data = face->data;
        if (data.size() < 4)
            return model;
        size_t offset = 0;
        uint32_t count = (data[offset] << 24) | (data[offset+1] << 16) | (data[offset+2] << 8) | data[offset+3];
        offset += 4;
        model.faces.resize(count);
        for (uint32_t i=0; i<count; ++i) {
            uint16_t i1 = (data[offset] << 8) | data[offset+1]; offset += 2;
            uint16_t i2 = (data[offset] << 8) | data[offset+1]; offset += 2;
            uint16_t i3 = (data[offset] << 8) | data[offset+1]; offset += 2;
            model.faces[i] = {i1, i2, i3};
        }
    }

    return model;
}

} // namespace wc_iff_loader
