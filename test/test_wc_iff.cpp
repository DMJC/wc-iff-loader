#include "wc_iff_loader/IFFReader.h"
#include "wc_iff_loader/WC3Decoder.h"
#include <cassert>
#include <fstream>
#include <cstring>
#include <vector>

using namespace wc_iff_loader;

static void writeU32BE(std::vector<uint8_t>& buf, uint32_t v) {
    buf.push_back((v >> 24) & 0xFF);
    buf.push_back((v >> 16) & 0xFF);
    buf.push_back((v >> 8) & 0xFF);
    buf.push_back(v & 0xFF);
}

static void writeF32BE(std::vector<uint8_t>& buf, float f) {
    uint32_t u;
    std::memcpy(&u, &f, sizeof(u));
    writeU32BE(buf, u);
}

static void writeU16BE(std::vector<uint8_t>& buf, uint16_t v) {
    buf.push_back((v >> 8) & 0xFF);
    buf.push_back(v & 0xFF);
}

int main() {
    std::vector<uint8_t> file;
    // Root FORM header
    file.insert(file.end(), {'F','O','R','M'});
    // Placeholder size
    writeU32BE(file, 0);
    // Type
    file.insert(file.end(), {'W','C','3','M'});

    std::vector<uint8_t> vert;
    writeU32BE(vert, 3); // num vertices
    writeF32BE(vert, 1.0f); writeF32BE(vert, 0.0f); writeF32BE(vert, 0.0f);
    writeF32BE(vert, 0.0f); writeF32BE(vert, 1.0f); writeF32BE(vert, 0.0f);
    writeF32BE(vert, 0.0f); writeF32BE(vert, 0.0f); writeF32BE(vert, 1.0f);
    std::vector<uint8_t> face;
    writeU32BE(face, 1); // num faces
    writeU16BE(face, 0); writeU16BE(face, 1); writeU16BE(face, 2);

    auto startSize = file.size();
    file.insert(file.end(), {'V','E','R','T'});
    writeU32BE(file, vert.size());
    file.insert(file.end(), vert.begin(), vert.end());
    if (vert.size() & 1) file.push_back(0);

    file.insert(file.end(), {'F','A','C','E'});
    writeU32BE(file, face.size());
    file.insert(file.end(), face.begin(), face.end());
    if (face.size() & 1) file.push_back(0);

    uint32_t formSize = file.size() - 8; // exclude 'FORM' and size field
    file[4] = (formSize >> 24) & 0xFF;
    file[5] = (formSize >> 16) & 0xFF;
    file[6] = (formSize >> 8) & 0xFF;
    file[7] = formSize & 0xFF;

    std::ofstream ofs("test.iff", std::ios::binary);
    ofs.write(reinterpret_cast<const char*>(file.data()), file.size());

    IFFReader reader;
    ofs.close();
    assert(reader.load("test.iff"));

    Model model = WC3Decoder::decode(reader.root());
    assert(model.vertices.size() == 3);
    assert(model.faces.size() == 1);
    return 0;
}
