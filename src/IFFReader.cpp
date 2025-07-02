#include "wc_iff_loader/IFFReader.h"
#include <fstream>
#include <cstring>
#include <stdexcept>

namespace wc_iff_loader {

static bool isContainer(const std::string& id) {
    return id == "FORM" || id == "LIST" || id == "CAT ";
}

uint32_t IFFReader::readU32BE(std::istream& is) {
    uint8_t b[4];
    is.read(reinterpret_cast<char*>(b), 4);
    return (static_cast<uint32_t>(b[0]) << 24) |
           (static_cast<uint32_t>(b[1]) << 16) |
           (static_cast<uint32_t>(b[2]) << 8)  |
           static_cast<uint32_t>(b[3]);
}

float IFFReader::readF32BE(std::istream& is) {
    uint32_t u = readU32BE(is);
    float f;
    std::memcpy(&f, &u, sizeof(f));
    return f;
}

bool IFFReader::load(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file)
        return false;

    char id[4];
    file.read(id, 4);
    if (file.gcount() != 4 || std::string(id,4) != "FORM")
        return false;

    rootChunk.id = "FORM";
    rootChunk.size = readU32BE(file);
    char type[4];
    file.read(type, 4);
    rootChunk.type.assign(type, 4);

    if (!parseChunk(file, rootChunk, rootChunk.size - 4))
        return false;

    return true;
}

bool IFFReader::parseChunk(std::istream& is, IffChunk& parent, uint32_t size) {
    uint32_t bytesRead = 0;
    while (bytesRead < size) {
        char id[4];
        is.read(id, 4);
        if (is.gcount() != 4)
            return false;
        uint32_t chunkSize = readU32BE(is);
        IffChunk chunk;
        chunk.id.assign(id, 4);
        chunk.size = chunkSize;
        bytesRead += 8 + chunkSize + (chunkSize & 1);

        if (isContainer(chunk.id)) {
            char type[4];
            is.read(type, 4);
            chunk.type.assign(type, 4);
            if (!parseChunk(is, chunk, chunkSize - 4))
                return false;
        } else {
            chunk.data.resize(chunkSize);
            is.read(reinterpret_cast<char*>(chunk.data.data()), chunkSize);
        }

        if (chunkSize & 1)
            is.get(); // skip pad

        parent.children.push_back(std::move(chunk));
    }
    return true;
}

} // namespace wc_iff_loader

