// wc3_viewer.cpp — WC3 REAL IFF viewer with textures + OBJ/MTL export + auto-fit draw
// Build:  g++ wc3_viewer.cpp -std=c++17 -O2 -lSDL2 -lGLEW -lGL -o wc3_viewer

#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <optional>
#include <cctype>

#if defined(__has_include)
#  if __has_include(<openvr.h>)
#    define WC_HAVE_OPENVR 1
#    include <openvr.h>
#  else
#    define WC_HAVE_OPENVR 0
#  endif
#else
#  define WC_HAVE_OPENVR 0
#endif

using std::string;
using std::vector;

// ----------------------------- small utils -----------------------------
static uint32_t be32(const uint8_t* p){ return (uint32_t(p[0])<<24)|(uint32_t(p[1])<<16)|(uint32_t(p[2])<<8)|uint32_t(p[3]); }
static uint16_t le16u(const uint8_t* p){ return uint16_t(p[0]) | (uint16_t(p[1])<<8); }
static int32_t  le32s(const uint8_t* p){ return int32_t(uint32_t(p[0]) | (uint32_t(p[1])<<8) | (uint32_t(p[2])<<16) | (uint32_t(p[3])<<24)); }
static float    le32f(const uint8_t* p){ float f; std::memcpy(&f,p,4); return f; }

static std::string stemFromPath(const std::string& p){
    size_t s = p.find_last_of("/\\");
    std::string fn = (s==std::string::npos) ? p : p.substr(s+1);
    size_t d = fn.find_last_of('.');
    return (d==std::string::npos) ? fn : fn.substr(0,d);
}

struct Vec3{ float x=0,y=0,z=0; };
static Vec3 vecAdd(const Vec3& a, const Vec3& b) { return Vec3{ a.x + b.x, a.y + b.y, a.z + b.z }; }
static Vec3 vecScale(const Vec3& v, float s) { return Vec3{ v.x * s, v.y * s, v.z * s }; }
static float vecLength(const Vec3& v) { return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }
static Vec3 vecNormalize(const Vec3& v) { float l = vecLength(v); if (l <= 1e-6f) return Vec3{ 0,0,0 }; return Vec3{ v.x / l, v.y / l, v.z / l }; }
static Vec3 vecCross(const Vec3& a, const Vec3& b) {
    return Vec3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}
static bool isBackOrFrontLabel(const std::string& name) {
    if (name.empty()) return false;
    std::string upper = name;
    std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) { return (char)std::toupper(c); });
    return upper == "BACK" || upper == "FRONT" || upper == "BMD9" || upper == "BSD1" || upper == "BMDX5";
}
struct Tri{
    uint32_t v[3]{};
    uint16_t tex=0;       // texture index
    float uv[6]{0,0,0,0,0,0}; // pixel UVs (u0,v0,u1,v1,u2,v2)
    bool hasTex=false;
};
struct Texture{
    int w=0,h=0;
    vector<uint8_t> rgba; // RGBA8
    vector<uint8_t> idx;  // original palette indices
    GLuint gl=0;
    std::string name;     // TXMP name (up to 8 chars)
    bool skipRender = false; // true when texture should not be rendered
    bool valid() const { return w>0 && h>0 && rgba.size()==(size_t)w*h*4; }
};

struct ObjGroup {
    std::string name;
    size_t firstTri = 0;  // index into M.tris where this group starts
    size_t triCount = 0;  // how many triangles belong to this group
};

struct SubModel;
struct Model{
    vector<Vec3> verts;
    vector<Tri>  tris;
    vector<Texture> textures;
    string name;
    vector<SubModel> submodels;
    vector<ObjGroup> groups;
    std::string afterburnerModelName;
    std::vector<glm::vec3> afterburnerOffsets;
    std::vector<Model> afterburnerFrames;
};
struct SubModel {
    std::string name;
    std::string weapon;
    glm::vec3   translate;   // local placement
    Model       mesh;
    int         parentIndex; // -1 = parent is the ship
};
struct TurretDef {
    std::string baseModel;   // "<prefix>_TRT" -> file "<prefix>_TRT.IFF"
    std::string gunModel;    // "<prefix>_GUN" -> file "<prefix>_GUN.IFF"
    std::string weaponName;  // "TURLASER" / "ANTIGUN" etc
    // Limits (degrees)
    int yawMin   = -180;
    int yawMax   =  180;
    int pitchMin =  -90;
    int pitchMax =   90;
    int elevCap  =   90;     // Often shows as 0x5A (90)
    int guns = 1;
    // Mount position in ship-local units (heuristic; see notes)
    Vec3 pos{0,0,0};
};

// ----------------------------- global palette -----------------------------
static std::array<uint8_t,3> gPalette[256];
static int gPalBankCount = 1;

static void makeDefaultPalette(){
    for(int i=0;i<256;i++){
        gPalette[i][0]=gPalette[i][1]=gPalette[i][2]=(uint8_t)i;
    }
}

static bool loadPaletteJSON(const std::string& path, int palOffset){
    std::ifstream f(path);
    if(!f){ std::cerr<<"[pal] can't open "<<path<<"\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(f)), {});
    std::vector<int> nums; nums.reserve(256*3);
    int v=0, sign=1; bool in=false;
    auto flush=[&]{ if(in){ nums.push_back(sign*v); v=0; sign=1; in=false; } };
    for(char c: s){
        if(c=='-'){ sign=-1; in=true; }
        else if(c>='0' && c<='9'){ v=v*10+(c-'0'); in=true; }
        else flush();
    }
    flush();
    if(nums.size() < 256*3){
        std::cerr<<"[pal] not enough values ("<<nums.size()<<")\n"; return false;
    }
    gPalBankCount = (int)(nums.size() / (256*3));
    if(gPalBankCount < 1) gPalBankCount = 1;
    int sel = std::clamp(palOffset, 0, gPalBankCount-1);
    size_t base = (size_t)sel * 256u * 3u;
    for(int i=0;i<256;i++){
        gPalette[i][0] = (uint8_t)std::clamp(nums[base+i*3+0],0,255);
        gPalette[i][1] = (uint8_t)std::clamp(nums[base+i*3+1],0,255);
        gPalette[i][2] = (uint8_t)std::clamp(nums[base+i*3+2],0,255);
    }
    std::cerr<<"[pal] loaded "<<path<<" banks="<<gPalBankCount<<" using bank "<<sel<<"\n";
    return true;
}

// ----------------------------- IFF parser -----------------------------
struct Chunk{
    std::string id;       // "FORM" or actual id
    std::string formType; // when id=="FORM"
    uint32_t size=0;      // payload size (EA IFF big-endian)
    size_t   start=0;     // file offset at chunk header
    size_t   payload=0;   // offset to payload (start+8 or +12 for FORM)
    std::vector<Chunk> children; // if FORM
};
struct IFF{
    std::vector<uint8_t> buf;
    Chunk root;

    Chunk parse(size_t& off, size_t end){
        Chunk c;
        if(off+8 > end){ c.id=""; return c; }
        c.start = off;
        c.id.assign((const char*)&buf[off], 4);
        c.size = be32(&buf[off+4]);
        c.payload = off+8;
        off += 8;
        if(c.id=="FORM"){
            if(off+4 > end){ c.id=""; return c; }
            c.formType.assign((const char*)&buf[off],4);
            off += 4;
            size_t formEnd = c.payload + c.size;
            while(off+8 <= formEnd){
                c.children.push_back(parse(off, formEnd));
                if(off & 1) off++;
            }
            off = formEnd;
        }else{
            off += c.size;
        }
        if(off & 1) off++;
        return c;
    }
    bool load(const std::string& path){
        std::ifstream f(path, std::ios::binary);
        if(!f) return false;
        buf.assign(std::istreambuf_iterator<char>(f), {});
        if(buf.size()<12) return false;
        size_t off=0;
        root = parse(off, buf.size());
        return true;
    }
};

static const Chunk* findFirst(const Chunk* c, const char* id){
    if(!c) return nullptr;
    if(c->id == id) return c;
    for(const auto& ch : c->children){
        if(const Chunk* r = findFirst(&ch, id)) return r;
    }
    return nullptr;
}
static const Chunk* findFirst(const Chunk* c, const char* form, const char* formType){
    if(!c) return nullptr;
    if(c->id==form && c->formType==formType) return c;
    for(const auto& ch : c->children){
        if(const Chunk* r = findFirst(&ch, form, formType)) return r;
    }
    return nullptr;
}
static const Chunk* childOf(const Chunk* c, const char* id){
    if(!c) return nullptr;
    for(const auto& ch : c->children) if(ch.id==id) return &ch;
    return nullptr;
}

// ----------------------------- TXMP decode (WC3) -----------------------------
// TXMP payload layout (WC3):
// +0..7  : name (ASCII, 8 bytes)
// +8..9  : width  (u16 LE)
// +10..11: height (u16 LE)
// +12..  : image data; if len == w*h -> RAW pal8, else -> one whole-image RLE stream.
static bool rle_expand(const uint8_t* src, size_t n, size_t expected, vector<uint8_t>& out){
    out.clear(); out.reserve(expected);
    size_t i = 0;
    while (i < n && out.size() < expected){
        uint8_t c = src[i++];
        if (c & 0x80){
            size_t cnt = (c & 0x7F) + 1;
            if (i >= n) return false;
            uint8_t val = src[i++];
            size_t need = expected - out.size();
            if (cnt > need) cnt = need;
            out.insert(out.end(), cnt, val);
        } else {
            size_t cnt = (size_t)c + 1;
            size_t need = expected - out.size();
            size_t take = std::min(cnt, need);
            if (i + take > n) return false;
            out.insert(out.end(), src + i, src + i + take);
            i += cnt;
        }
    }
    return out.size() == expected;
}

static void palToRGBA(const vector<uint8_t>& idx, int w, int h, vector<uint8_t>& rgba){
    rgba.resize((size_t)w*h*4);
    auto clamp255=[](int x){ return (uint8_t)(x<0?0:(x>255?255:x)); };
    for(size_t i=0;i<idx.size();++i){
        uint8_t pi = idx[i];
        const auto& c = gPalette[pi];
        rgba[i*4+0] = clamp255((int)c[0] + 30);
        rgba[i*4+1] = clamp255((int)c[1] + 30);
        rgba[i*4+2] = clamp255((int)c[2] + 30);
        rgba[i*4+3] = (pi==255) ? 0 : 255;
    }
}

static bool decode_TXMP_WC3(const uint8_t* data, size_t len, int& w, int& h,
                            std::vector<uint8_t>& idx, std::vector<uint8_t>& rgba){
    auto sane_dims = [](int W,int H)->bool{
        return (W>0 && H>0 && W<=4096 && H<=4096 && (int64_t)W*(int64_t)H <= (1<<26));
    };

    struct Variant { int offW, offH, offP; };
    const Variant variants[] = {
        {  8, 10, 12 }, // 12-byte header (name[0..7], w@+8,  h@+10, pixels@+12)
        { 16, 18, 20 }, // 20-byte header (name[8..15], w@+16, h@+18, pixels@+20) - out-min.js style
    };

    for (auto v : variants){
        if ((size_t)v.offH + 2 > len) continue;
        int W = (int)le16u(data + v.offW);
        int H = (int)le16u(data + v.offH);
        if (!sane_dims(W,H)) continue;

        const uint8_t* img = (len > (size_t)v.offP) ? data + v.offP : nullptr;
        size_t ilen = img ? (len - (size_t)v.offP) : 0;
        size_t need = (size_t)W * (size_t)H;

        std::vector<uint8_t> idxLocal;
        bool ok = false;

        // --- RAW (with optional trailing pad) ---
        // Some WC3 TXMPs store RAW pal8 plus 1 padding byte (observed on EXCAL: 01,12,15).
        if (ilen >= need && ilen - need <= 2) {
            idxLocal.assign(img, img + need);   // ignore pad
            ok = true;
        }

        // --- RLE fallback (whole-image RLE stream) ---
        if (!ok && img) {
            if (rle_expand(img, ilen, need, idxLocal)) ok = true;
        }

        if (ok) {
            idx = std::move(idxLocal);
            palToRGBA(idx, W, H, rgba);
            w = W; h = H;
            // Optional debug:
            // if (ilen > need) std::cerr << "[TXMP] RAW+pad (" << (ilen - need) << " byte) " << W << "x" << H << "\n";
            return true;
        }
    }
    return false;
}

// ----------------------------- TGA writer -----------------------------
static bool writeTGA(const std::string& path, int w, int h, const uint8_t* rgba){
    if(w<=0 || h<=0 || !rgba) return false;
    std::ofstream f(path, std::ios::binary);
    if(!f) return false;

    uint8_t hdr[18] = {};
    hdr[2] = 2; // uncompressed truecolor
    hdr[12] = w & 0xFF; hdr[13] = (w>>8)&0xFF;
    hdr[14] = h & 0xFF; hdr[15] = (h>>8)&0xFF;
    hdr[16] = 32;       // BGRA32
    hdr[17] = 0x20;     // top-left origin
    f.write((char*)hdr,18);

    // RGBA -> BGRA per pixel
    vector<uint8_t> line((size_t)w*4);
    for(int y=0;y<h;y++){
        const uint8_t* s = rgba + (size_t)y*w*4;
        uint8_t* d = line.data();
        for(int x=0;x<w;x++){
            d[0]=s[2]; d[1]=s[1]; d[2]=s[0]; d[3]=s[3];
            s+=4; d+=4;
        }
        f.write((char*)line.data(), line.size());
    }
    return true;
}

static void appendModelInto(Model& dst,
                            const Model& src,
                            const glm::vec3& translate,
                            const std::string& groupName)
{
    const size_t vBase = dst.verts.size();
    const size_t tBase = dst.textures.size();
    const size_t triStart = dst.tris.size();

    // translate & copy verts
    dst.verts.reserve(dst.verts.size() + src.verts.size());
    for (const auto& v : src.verts) {
        dst.verts.push_back({ v.x + translate.x,
                              v.y + translate.y,
                              v.z + translate.z });
    }

    // append textures
    dst.textures.reserve(dst.textures.size() + src.textures.size());
    for (const auto& tex : src.textures) dst.textures.push_back(tex);

    // reindex tris and remap texture indices
    dst.tris.reserve(dst.tris.size() + src.tris.size());
    for (const auto& s : src.tris) {
        Tri d = s;
        d.v[0] = s.v[0] + (uint32_t)vBase;
        d.v[1] = s.v[1] + (uint32_t)vBase;
        d.v[2] = s.v[2] + (uint32_t)vBase;
        if (s.hasTex) { d.tex = (uint16_t)(s.tex + tBase); d.hasTex = true; }
        dst.tris.push_back(d);
    }

    // optional group for clean OBJ sub-objects
    if constexpr (true) { // flip to false if you don’t want groups
        if (!std::is_sorted(groupName.begin(), groupName.end())) { /*no-op, just avoid unused warning*/ }
        if constexpr (true) {
            // If Model already has groups vector:
            // dst.groups.push_back({groupName, triStart, dst.tris.size()-triStart});
        }
    }
}

static void flattenSubmodelsInto(Model& dst)
{
    if (dst.submodels.empty()) return;

    // Gather for stable iteration (in case we move-from)
    std::vector<SubModel> subs = dst.submodels;
    dst.submodels.clear();

    for (const auto& sm : subs) {
        const Model& src = sm.mesh;

        const size_t vBase = dst.verts.size();
        const size_t tBase = dst.textures.size();
        const size_t triStart = dst.tris.size();

        // 1) append translated vertices
        dst.verts.reserve(dst.verts.size() + src.verts.size());
        for (const auto& v : src.verts) {
            dst.verts.push_back({ v.x + sm.translate.x,
                                  v.y + sm.translate.y,
                                  v.z + sm.translate.z });
        }

        // 2) append textures (keep indices stable by offsetting)
        dst.textures.reserve(dst.textures.size() + src.textures.size());
        for (const auto& tex : src.textures) dst.textures.push_back(tex);

        // 3) append triangles with reindex + texture remap
        dst.tris.reserve(dst.tris.size() + src.tris.size());
        for (const auto& s : src.tris) {
            Tri d = s;
            d.v[0] = s.v[0] + (uint32_t)vBase;
            d.v[1] = s.v[1] + (uint32_t)vBase;
            d.v[2] = s.v[2] + (uint32_t)vBase;
            if (s.hasTex) {
                d.tex = (uint16_t)(s.tex + tBase);
                d.hasTex = true;
            }
            dst.tris.push_back(d);
        }

        // 4) record group range for OBJ export
        ObjGroup g;
        g.name     = sm.name;           // e.g. "turret_base_0" / "turret_gun_0"
        g.firstTri = triStart;
        g.triCount = dst.tris.size() - triStart;
        dst.groups.push_back(std::move(g));
    }
}

static void applyRotation(Model& m, float yawDeg, float pitchDeg, float rollDeg)
{
    glm::mat4 R(1.0f);
    if (std::fabs(yawDeg)   > 1e-6f) R = glm::rotate(R, glm::radians(yawDeg),   glm::vec3(0,1,0));
    if (std::fabs(pitchDeg) > 1e-6f) R = glm::rotate(R, glm::radians(pitchDeg), glm::vec3(1,0,0));
    if (std::fabs(rollDeg)  > 1e-6f) R = glm::rotate(R, glm::radians(rollDeg),  glm::vec3(0,0,1));
    for (auto& v : m.verts) {
        glm::vec4 p(v.x, v.y, v.z, 1.0f);
        p = R * p;
        v.x = p.x; v.y = p.y; v.z = p.z;
    }
}

// ----------------------------- Loader (HCl geometry + textures) -----------------------------
static bool load_wc3_model_hcl_textured(const string& path, Model& M, std::vector<Model>* lodFrames=nullptr){
    // Minimal IFF loader (recurses already in IFF::load)
    IFF iff; if(!iff.load(path)){ std::cerr<<"Not a REAL/FORM IFF\n"; return false; }
    M.name = stemFromPath(path);

    std::string baseDir = path;
    {
        auto pos = baseDir.find_last_of("/\\");
        baseDir = (pos == std::string::npos) ? std::string() : baseDir.substr(0, pos + 1);
    }

    auto loadSubModel = [&](const std::string& rawName, Model& out) -> bool {
        if (rawName.empty()) return false;
        std::vector<std::string> candidates;
        auto addCandidate = [&](std::string name) {
            if (name.empty()) return;
            if (std::find(candidates.begin(), candidates.end(), name) == candidates.end()) {
                candidates.push_back(std::move(name));
            }
        };
        addCandidate(rawName);
        std::string upper = rawName;
        std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c){ return (char)std::toupper(c); });
        addCandidate(upper);
        std::string lower = rawName;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c){ return (char)std::tolower(c); });
        addCandidate(lower);

        for (const auto& candidate : candidates) {
            std::string file = baseDir + candidate + ".IFF";
            if (file == path) continue; // avoid immediate self-recursion
            std::ifstream f(file, std::ios::binary);
            if (!f) continue;
            f.close();
            if (load_wc3_model_hcl_textured(file, out)) return true;
        }
        return false;
    };

    const Chunk* TURT = findFirst(&iff.root,"TURT");
    const Chunk* VERT = findFirst(&iff.root,"VERT");
    const Chunk* TRIS = findFirst(&iff.root,"FORM","TRIS");
    const Chunk* QUAD = findFirst(&iff.root,"FORM","QUAD");
    std::array<std::vector<uint16_t>, 8> lvlFaces;
    for (int lvl = 0; lvl < 8; ++lvl) {
        char id[5];
        std::snprintf(id, sizeof(id), "LVL%d", lvl);
        if (const Chunk* chunk = findFirst(&iff.root, id)) {
            const uint8_t* p = &iff.buf[chunk->start + 8];
            uint32_t sz = be32(&iff.buf[chunk->start + 4]);
            size_t n = sz / 2;
            auto& dst = lvlFaces[lvl];
            dst.resize(n);
            for (size_t i = 0; i < n; ++i) dst[i] = le16u(p + i * 2);
        }
    }
    const Chunk* TXMS = findFirst(&iff.root,"FORM","TXMS");
    if(!VERT){ std::cerr<<"VERT missing\n"; return false; }

    // VERT: prefer LE int32/256; fallback to float32 if insane
    {
        const uint8_t* p=&iff.buf[VERT->start+8];
        uint32_t sz = be32(&iff.buf[VERT->start+4]);
        size_t n=sz/12; M.verts.resize(n);
        bool insane=false;
        for(size_t i=0;i<n;i++){
            float x = le32s(p+i*12+0)/256.0f;
            float y = le32s(p+i*12+4)/256.0f;
            float z = le32s(p+i*12+8)/256.0f;
            if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z) ||
               std::fabs(x)>1e6f || std::fabs(y)>1e6f || std::fabs(z)>1e6f){
                insane=true; break;
            }
            M.verts[i]={x,y,z};
        }
        if(insane){
            for(size_t i=0;i<n;i++){
                float x = le32f(p+i*12+0);
                float y = le32f(p+i*12+4);
                float z = le32f(p+i*12+8);
                M.verts[i]={x,y,z};
            }
        }
    }

    const Chunk* T_FACE = TRIS ? childOf(TRIS,"FACE") : nullptr;
    const Chunk* T_MAPS = TRIS ? childOf(TRIS,"MAPS") : nullptr;
    const Chunk* Q_FACE = QUAD ? childOf(QUAD,"FACE") : nullptr;
    const Chunk* Q_MAPS = QUAD ? childOf(QUAD,"MAPS") : nullptr;

    struct FaceT { uint16_t flag=0; uint32_t v[3]{}; bool hasTex=false; uint16_t tex=0; uint16_t uvpx[6]{}; };
    struct FaceQ { uint16_t flag=0; uint32_t v[4]{}; bool hasTex=false; uint16_t tex=0; uint16_t uvpx[8]{}; };
    vector<FaceT> tris; vector<FaceQ> quads;

    if(T_FACE){
        const uint8_t* f=&iff.buf[T_FACE->start+8];
        uint32_t sz = be32(&iff.buf[T_FACE->start+4]);
        size_t rows=sz/(4*2);
        tris.resize(rows);
        for(size_t i=0;i<rows;i++){
            const uint8_t* r=f+i*8;
            tris[i].flag = le16u(r+0);
            tris[i].v[0] = le16u(r+2);
            tris[i].v[1] = le16u(r+4);
            tris[i].v[2] = le16u(r+6);
        }
    }
    if(Q_FACE){
        const uint8_t* f=&iff.buf[Q_FACE->start+8];
        uint32_t sz = be32(&iff.buf[Q_FACE->start+4]);
        size_t rows=sz/(5*2);
        quads.resize(rows);
        for(size_t i=0;i<rows;i++){
            const uint8_t* r=f+i*10;
            quads[i].flag = le16u(r+0);
            quads[i].v[0] = le16u(r+2);
            quads[i].v[1] = le16u(r+4);
            quads[i].v[2] = le16u(r+6);
            quads[i].v[3] = le16u(r+8);
        }
    }

    if(T_MAPS){
        const uint8_t* m=&iff.buf[T_MAPS->start+8];
        uint32_t sz = be32(&iff.buf[T_MAPS->start+4]);
        size_t rows=sz/(8*2);
        for(size_t i=0;i<rows;i++){
            const uint8_t* r=m+i*16;
            uint16_t fi  = le16u(r+0);
            uint16_t tex = le16u(r+2);
            if(fi < tris.size()){
                auto& ft = tris[fi]; ft.hasTex=true; ft.tex=tex;
                for(int k=0;k<6;k++) ft.uvpx[k]=le16u(r+4+k*2);
            }
        }
    }
    if(Q_MAPS){
        const uint8_t* m=&iff.buf[Q_MAPS->start+8];
        uint32_t sz = be32(&iff.buf[Q_MAPS->start+4]);
        size_t rows=sz/(10*2);
        for(size_t i=0;i<rows;i++){
            const uint8_t* r=m+i*20;
            uint16_t fi  = le16u(r+0);
            uint16_t tex = le16u(r+2);
            if(fi < quads.size()){
                auto& fq = quads[fi]; fq.hasTex=true; fq.tex=tex;
                for(int k=0;k<8;k++) fq.uvpx[k]=le16u(r+4+k*2);
            }
        }
    }

    // Fallback for TRIS without FACE (INDX/MAPS)
    if(!T_FACE && TRIS){
        const Chunk* INDX=nullptr; const Chunk* MAPS=nullptr;
        for(const auto& ch: TRIS->children){
            if(ch.id=="INDX") INDX=&ch;
            else if(ch.id=="MAPS") MAPS=&ch;
        }
        if(INDX){
            size_t ntri = be32(&iff.buf[INDX->start+4]) / (sizeof(uint16_t)*3);
            const uint8_t* p = &iff.buf[INDX->start+8];
            tris.resize(ntri);
            for(size_t i=0;i<ntri;i++){
                tris[i].v[0] = le16u(p + i*6 + 0);
                tris[i].v[1] = le16u(p + i*6 + 2);
                tris[i].v[2] = le16u(p + i*6 + 4);
            }
        }
        if(MAPS){
            const uint8_t* m = &iff.buf[MAPS->start+8];
            size_t rows = be32(&iff.buf[MAPS->start+4]) / 16;
            for(size_t r=0;r<rows;r++){
                const uint8_t* row = m + r*16;
                uint16_t fi  = le16u(row+0);
                uint16_t tex = le16u(row+2);
                if(fi < tris.size()){
                    auto& ft = tris[fi]; ft.hasTex=true; ft.tex=tex;
                    for(int k=0;k<6;k++) ft.uvpx[k]=le16u(row+4+k*2);
                }
            }
        }
    }

    auto inRange=[&](uint32_t i){ return i<M.verts.size(); };
    auto passFlag=[](uint16_t flag){ return (flag & 0xFF) == 1; };

    auto addTriTo = [&](Model& dest, const FaceT& ft){
        if(!(inRange(ft.v[0])&&inRange(ft.v[1])&&inRange(ft.v[2]))) return;
        Tri t{}; t.tex = ft.tex;
        // HCl viewer draws triangles as (v2,v1,v0)
        t.v[0]=ft.v[2]; t.v[1]=ft.v[1]; t.v[2]=ft.v[0];
        if(ft.hasTex){
            // uv for (v0,v1,v2) -> remap to (v2,v1,v0)
            t.uv[0]=(float)ft.uvpx[4]; t.uv[1]=(float)ft.uvpx[5];
            t.uv[2]=(float)ft.uvpx[2]; t.uv[3]=(float)ft.uvpx[3];
            t.uv[4]=(float)ft.uvpx[0]; t.uv[5]=(float)ft.uvpx[1];
            t.hasTex=true;
        }
        dest.tris.push_back(t);
    };

    auto addQuadAsTwo = [&](Model& dest, const FaceQ& fq){
        if(!(inRange(fq.v[0])&&inRange(fq.v[1])&&inRange(fq.v[2])&&inRange(fq.v[3]))) return;
        // tri A: (v2,v1,v0)
        {
            Tri t{}; t.tex=fq.tex;
            t.v[0]=fq.v[2]; t.v[1]=fq.v[1]; t.v[2]=fq.v[0];
            if(fq.hasTex){
                t.uv[0]=(float)fq.uvpx[4]; t.uv[1]=(float)fq.uvpx[5];
                t.uv[2]=(float)fq.uvpx[2]; t.uv[3]=(float)fq.uvpx[3];
                t.uv[4]=(float)fq.uvpx[0]; t.uv[5]=(float)fq.uvpx[1];
                t.hasTex=true;
            }
            dest.tris.push_back(t);
        }
        // tri B: (v3,v2,v0)
        {
            Tri t{}; t.tex=fq.tex;
            t.v[0]=fq.v[3]; t.v[1]=fq.v[2]; t.v[2]=fq.v[0];
            if(fq.hasTex){
                t.uv[0]=(float)fq.uvpx[6]; t.uv[1]=(float)fq.uvpx[7];
                t.uv[2]=(float)fq.uvpx[4]; t.uv[3]=(float)fq.uvpx[5];
                t.uv[4]=(float)fq.uvpx[0]; t.uv[5]=(float)fq.uvpx[1];
                t.hasTex=true;
            }
            dest.tris.push_back(t);
        }
    };
    auto populateModel = [&](Model& dest, const std::vector<uint16_t>& order, bool fallback){
        dest.tris.clear();
        if(!order.empty()){
            for(uint16_t idx : order){
                if(idx < tris.size()){
                    const auto& ft=tris[idx];
                    if(ft.hasTex || passFlag(ft.flag)) addTriTo(dest, ft);
                }else{
                    size_t qi = idx - (uint16_t)tris.size();
                    if(qi<quads.size()){
                        const auto& fq=quads[qi];
                        if(fq.hasTex || passFlag(fq.flag)) addQuadAsTwo(dest, fq);
                    }
                }
            }
        }else if(fallback){
            for(const auto& ft: tris) addTriTo(dest, ft);
            for(const auto& fq: quads) addQuadAsTwo(dest, fq);
        }
    };

    populateModel(M, lvlFaces[0], true);

    // Textures
    M.textures.clear();
    if(TXMS){
        size_t ok=0, rawCnt=0, rleCnt=0, total=0;
        for (const auto& ch : TXMS->children) {
            if (ch.id != "TXMP") continue;
            ++total;
            const uint8_t* p = &iff.buf[ch.payload];
            size_t len = ch.size;
            int w = 0, h = 0; vector<uint8_t> rgba, idx;
            Texture T;
            if (len >= 8) {
                char nbuf[9];
                std::memcpy(nbuf, p, 8);
                nbuf[8] = 0;
                T.name = nbuf;
                while (!T.name.empty() && (T.name.back() == '\0' || T.name.back() == ' ')) T.name.pop_back();
                if (isBackOrFrontLabel(T.name)) {
                    T.skipRender = true;
                    std::cerr << "[tex] suppress " << T.name << " (BACK/FRONT)\n";
                }
            }
            if(decode_TXMP_WC3(p, len, w, h, idx, rgba)){
                if(12 + (size_t)w*h == len) ++rawCnt; else ++rleCnt;
                T.w=w; T.h=h; T.idx.swap(idx); T.rgba.swap(rgba); ++ok;
            }else{
                // placeholder to keep indices aligned
                T.w=2; T.h=2; T.rgba={255,0,255,255, 0,0,0,0, 0,0,0,0, 255,0,255,255};
            }
            M.textures.push_back(std::move(T));
        }
        std::cerr << "[tex] decoded " << ok << "/" << total
                  << " (RAW="<<rawCnt<<", RLE="<<rleCnt<<")\n";
    }

    if (lodFrames) {
        lodFrames->clear();
        lodFrames->resize(8);
        for (int lvl = 0; lvl < 8; ++lvl) {
            Model frame;
            frame.name = M.name + "_frame" + std::to_string(lvl);
            frame.verts = M.verts;
            frame.textures = M.textures;
            populateModel(frame, lvlFaces[lvl], lvl == 0);
            (*lodFrames)[lvl] = std::move(frame);
        }
    }
    std::cerr << "[geom] verts="<<M.verts.size()<<" tris="<<M.tris.size()<<" tex="<<M.textures.size()<<"\n";

    // Helpers for BE/LE (you already have le32s/le16u etc.)
    auto be32 = [&](const uint8_t* p){ return (uint32_t(p[0])<<24)|(uint32_t(p[1])<<16)|(uint32_t(p[2])<<8)|uint32_t(p[3]); };
    auto le32u = [&](const uint8_t* p){ return (uint32_t(p[3])<<24)|(uint32_t(p[2])<<16)|(uint32_t(p[1])<<8)|uint32_t(p[0]); };
    auto le32s = [&](const uint8_t* p){ return (int32_t)le32u(p); };
    
    // Find TURT bytes
    if (TURT) {
        const uint8_t* turt = &iff.buf[TURT->start+8];
        uint32_t turtSize   = be32(&iff.buf[TURT->start+4]);
        const uint8_t* tbeg = turt;
        const uint8_t* tend = turt + turtSize;
    
        struct TurretRec {
            std::string baseModel;  // "CRUL_TRT"
            std::string gunModel;   // "CRUL_GUN"
            std::string weapon;     // "TURLASER" / "ANTIGUN"
            glm::ivec3  pos{0,0,0}; // mount (int units); viewer verts use /256
            float       pitch = 0;  // degrees
            float       yaw   = 0;  // degrees
            int         elevCap = 90; // often 0x5A
            int         guns    = 1;  // often 1
        };
        std::vector<TurretRec> turrets;
    
        // Scan for "<NAME>_TRT\0<NAME>_GUN\0<WEAPON>\0" blocks terminated by 00000001
        const uint8_t* p = tbeg;
        while (p + 16 < tend) {
            // seek zero-terminated string
            if (!(*p >= 'A' && *p <= 'Z')) { ++p; continue; }
            const uint8_t* s1 = p;
            while (p < tend && *p != 0) ++p;
            if (p >= tend) break;
            const uint8_t* e1 = p;
            auto s1len = size_t(e1 - s1);
            if (s1len < 5 || s1len > 64) { ++p; continue; }
            if (!(s1len >= 4)) { ++p; continue; }
//            if (!(s1len >= 4 && (e1[-4]=='_' || e1[-4]=='T' ) && (e1[-3]=='T' || e1[-3]=='U') && e1[-2]=='R' && e1[-1]=='T')) { ++p; continue; }

            std::string baseModel((const char*)s1, s1len); // NAME_TRT

            // next "<NAME>_GUN\0"
            const uint8_t* s2 = e1 + 1; if (s2 >= tend) break;
            const uint8_t* e2 = s2; while (e2 < tend && *e2 != 0) ++e2;
            if (e2 >= tend) break;
            std::string gunModel((const char*)s2, size_t(e2 - s2));

            // next "<WEAPON>\0"
            const uint8_t* s3 = e2 + 1; if (s3 >= tend) break;
            const uint8_t* e3 = s3; while (e3 < tend && *e3 != 0) ++e3;
            if (e3 >= tend) break;
            std::string weapon((const char*)s3, size_t(e3 - s3));

            TurretRec rec; rec.baseModel = baseModel; rec.gunModel = gunModel; rec.weapon = weapon;

            // Find 00000001 terminator following the mount block
            const uint8_t* term = nullptr;
            for (const uint8_t* q = e3 + 1; q + 4 <= tend; ++q) {
                if (q[0]==0x00 && q[1]==0x00 && q[2]==0x00 && q[3]==0x01) { term = q; break; }
            }

            bool gotPos = false;
            if (term && (term - tbeg) >= 32) {
                const uint8_t* r = term - 1; // expect elevation marker
                    const uint8_t* w = term - 25; // "SER\0" + mount data
                        rec.pos.x = le32s(w + 4);
                        rec.pos.y = le32s(w + 8);
                        rec.pos.z = le32s(w + 12);
                        rec.pitch  = le32s(w + 16) / 256.0f;
                        rec.yaw    = le32s(w + 20) / 256.0f;
                        rec.elevCap = le32s(w + 24);
                        rec.guns    = 1; // not used for placement
                        gotPos = true;
            }

            if (!gotPos) {
                std::cerr << "[TURT] " << rec.baseModel << " invalid mount block; using (0,0,0)\n";
            }
            std::cerr << "[TURT] " << rec.baseModel << " " << rec.gunModel << " " << rec.weapon << " pos=("
                      << rec.pos.x << "," << rec.pos.y << "," << rec.pos.z << ")"
                      << " yaw=" << rec.yaw << " pitch=" << rec.pitch << "\n";
            turrets.push_back(std::move(rec));
            if (term) p = term + 4; else break; // continue after terminator
        }
	std::cerr << turrets.size() << " Turrets Found" << std::endl;
        // Load turret parts and attach as submodels
        if (!turrets.empty()) {
            for (size_t i = 0; i < turrets.size(); ++i) {
                const auto& T = turrets[i];

                // Base
                Model baseM;
                if (!loadSubModel(T.baseModel, baseM)) {
                    std::cerr << "[TURT] Missing " << T.baseModel << ".IFF\n";
                }

                // Gun
                Model gunM;
                if (!loadSubModel(T.gunModel, gunM)) {
                    std::cerr << "[TURT] Missing " << T.gunModel << ".IFF\n";
                }

                // Apply mount rotation
                applyRotation(baseM, T.yaw, T.pitch, 0.0f);
                applyRotation(gunM,  T.yaw, T.pitch, 0.0f);

                // Place base at mount (viewer verts are /256.0f)
                glm::vec3 basePos = glm::vec3(T.pos.x/256.0f, T.pos.y/256.0f, T.pos.z/256.0f);

                // Parent indices: ship = -1, base = idx of pushed base
                int baseIdx = (int)M.submodels.size();
                M.submodels.push_back(SubModel{
                    "turret_base_" + std::to_string(i),
                    T.weapon,
                    basePos,
                    std::move(baseM),
                    -1 // parent is ship
                });

                // Gun at mount, shares base rotation
                if (!gunM.tris.empty() || !gunM.verts.empty()) {
                    M.submodels.push_back(SubModel{
                        "turret_gun_" + std::to_string(i),
                        T.weapon,
                        basePos,
                        std::move(gunM),
                        baseIdx // parent is the base we just pushed
                    });
                }
                std::cerr << "[TURT] " << T.baseModel << " pos=("
                          << T.pos.x << "," << T.pos.y << "," << T.pos.z << ")"
                          << " yaw=" << T.yaw << " pitch=" << T.pitch << "\n";
            }
        }
    }


    if (const Chunk* APPR = findFirst(&iff.root, "FORM", "APPR")) {
        if (const Chunk* POLY = findFirst(APPR, "FORM", "POLY")) {
            if (const Chunk* SUPR = childOf(POLY, "SUPR")) {
                const uint8_t* suprData = &iff.buf[SUPR->start + 8];
                uint32_t suprSize = be32(&iff.buf[SUPR->start + 4]);
                if (suprSize < 3) {
                    std::cerr << "[SUPR] payload too small (" << suprSize << ")\n";
                } else {
                    uint16_t entryIndex = le16u(suprData + 0);
                    size_t nameOffset = 2;
                    size_t nameLen = 0;
                    while (nameOffset + nameLen < suprSize && suprData[nameOffset + nameLen] != 0) {
                        ++nameLen;
                    }
                    std::string hangarName;
                    if (nameLen > 0) {
                        hangarName.assign(reinterpret_cast<const char*>(suprData + nameOffset), nameLen);
                        while (!hangarName.empty() && hangarName.back() == ' ') hangarName.pop_back();
                    }
                    size_t tailStart = nameOffset + nameLen;
                    if (tailStart < suprSize && suprData[tailStart] == 0) ++tailStart;
                    int32_t parentIndex = -1;
                    if (tailStart < suprSize) {
                        size_t remain = suprSize - tailStart;
                        if (remain >= 4) {
                            parentIndex = le32s(suprData + tailStart);
                        } else if (remain == 3) {
                            uint32_t raw = uint32_t(suprData[tailStart + 0]) |
                                           (uint32_t(suprData[tailStart + 1]) << 8) |
                                           (uint32_t(suprData[tailStart + 2]) << 16);
                            if (raw & 0x800000) raw |= 0xFF000000u;
                            parentIndex = (int32_t)raw;
                        } else if (remain == 2) {
                            parentIndex = (int32_t)(int16_t)le16u(suprData + tailStart);
                        } else if (remain == 1) {
                            parentIndex = (int32_t)(int8_t)suprData[tailStart];
                        }
                    }

                    if (!hangarName.empty()) {
                        Model hangarModel;
                        if (!loadSubModel(hangarName, hangarModel)) {
                            std::cerr << "[SUPR] Missing hangar IFF " << hangarName << ".IFF\n";
                        } else {
                            std::string subName = "hangar_" + hangarName;
                            M.submodels.push_back(SubModel{
                                subName,
                                hangarName,
                                glm::vec3(0.0f),
                                std::move(hangarModel),
                                -1
                            });
                            std::cerr << "[SUPR] hangar=" << hangarName
                                      << " entry=" << entryIndex
                                      << " parent=" << parentIndex
                                      << "\n";
                        }
                    } else {
                        std::cerr << "[SUPR] empty hangar name\n";
                    }
                }
            }
        }
    }

    if (const Chunk* CRGO = findFirst(&iff.root, "CRGO")) {
        const uint8_t* cargoData = &iff.buf[CRGO->start + 8];
        uint32_t cargoSize = be32(&iff.buf[CRGO->start + 4]);
        const size_t entrySize = 37;
        size_t entryCount = cargoSize / entrySize;
        if (cargoSize % entrySize != 0) {
            std::cerr << "[CRGO] payload size " << cargoSize << " not aligned to " << entrySize << "\n";
        }
        auto normalizeDegrees = [](float deg) {
            if (!std::isfinite(deg)) return 0.0f;
            deg = std::fmod(deg, 360.0f);
            if (deg > 180.0f) deg -= 360.0f;
            if (deg < -180.0f) deg += 360.0f;
            return deg;
        };
        struct CargoRec {
            std::string name;
            glm::vec3   shipPosition{0.0f}; // decoded ship-space placement (X/Z 16.16, Y 8.8)
            float       yawDeg = 0.0f;
            float       pitchDeg = 0.0f;
            float       rollDeg = 0.0f;
            int16_t     sentinel = 0;
            int32_t     rawYaw = 0;
            int32_t     rawX = 0;
            int32_t     rawY = 0;
            int32_t     rawZ = 0;
        };
        std::vector<CargoRec> cargos;
        cargos.reserve(entryCount);
        for (size_t i = 0; i < entryCount; ++i) {
            const uint8_t* e = cargoData + i * entrySize;
            size_t len = 0;
            while (len < 16 && e[len] != 0) ++len;
            std::string name((const char*)e, len);
            while (!name.empty() && name.back() == ' ') name.pop_back();
            int32_t rawX     = le32s(e + 16);
            int32_t rawY     = le32s(e + 20);
            int32_t rawZ     = le32s(e + 24);
            int32_t rawYaw   = le32s(e + 28);
            int8_t  rawRoll  = (int8_t)e[32];
            int16_t rawPitch = (int16_t)le16u(e + 33);
            int16_t sentinel = (int16_t)le16u(e + 35);
            if (!name.empty()) {
                CargoRec rec;
                rec.name = std::move(name);
                rec.rawX = rawX;
                rec.rawY = rawY;
                rec.rawZ = rawZ;
                rec.rawYaw = rawYaw;
                rec.shipPosition = glm::vec3(
                    rawX / 65536.0f,
                    rawY / 256.0f,
                    rawZ / 65536.0f
                );
                rec.yawDeg   = normalizeDegrees(rawYaw / 65536.0f);
                rec.pitchDeg = rawPitch / 256.0f;
                rec.rollDeg  = static_cast<float>(rawRoll);
                rec.sentinel = sentinel;
                cargos.push_back(rec);
            }
        }
        if (!cargos.empty()) {
            std::cerr << "[CRGO] " << cargos.size() << " cargo mounts\n";
        }
        for (size_t i = 0; i < cargos.size(); ++i) {
            const auto& C = cargos[i];
            Model cargoModel;
            if (!loadSubModel(C.name, cargoModel)) {
                std::cerr << "[CRGO] Missing " << C.name << ".IFF\n";
                continue;
            }
            applyRotation(cargoModel, C.yawDeg, C.pitchDeg, C.rollDeg);
            glm::vec3 pos = C.shipPosition;
            std::string subName = "cargo_" + std::to_string(i) + "_" + C.name;
            M.submodels.push_back(SubModel{
                subName,
                C.name,
                pos,
                std::move(cargoModel),
                -1
            });
            std::cerr << "[CRGO] " << C.name
                      << " posRaw=(" << C.rawX << "," << C.rawY << "," << C.rawZ << ")"
                      << " pos=(" << pos.x << "," << pos.y << "," << pos.z << ")"
                      << " yawRaw=" << C.rawYaw
                      << " yaw=" << C.yawDeg
                      << " pitch=" << C.pitchDeg
                      << " roll=" << C.rollDeg
                      << " flag=" << C.sentinel
                      << "\n";
        }
    }

    if (const Chunk* AFTB = findFirst(&iff.root, "AFTB")) {
        const uint8_t* data = &iff.buf[AFTB->start + 8];
        uint32_t payloadSize = be32(&iff.buf[AFTB->start + 4]);
        if (payloadSize < 16) {
            std::cerr << "[AFTB] payload too small (" << payloadSize << ")\n";
        } else {
            uint32_t engineCount = (uint32_t)le32s(data + 0);
            (void)le32s(data + 4); // reserved/unused field
            size_t nameBytes = std::min<size_t>(8, payloadSize - 8);
            std::string abName(reinterpret_cast<const char*>(data + 8), nameBytes);
            while (!abName.empty() && (abName.back() == '\0' || abName.back() == ' ')) abName.pop_back();
            size_t pos = 16;
            size_t available = (payloadSize > pos) ? (payloadSize - pos) : 0;
            size_t entries = available / 12;
            size_t engines = std::min<size_t>(engineCount, entries);
            M.afterburnerModelName = abName;
            M.afterburnerOffsets.clear();
            for (size_t i = 0; i < engines; ++i) {
                const uint8_t* e = data + pos + i * 12;
                glm::vec3 offset(
                    le32s(e + 0) / 256.0f,
                    le32s(e + 4) / 256.0f,
                    le32s(e + 8) / 256.0f
                );
                M.afterburnerOffsets.push_back(offset);
            }
            if (engines < engineCount) {
                std::cerr << "[AFTB] expected " << engineCount << " engine mounts, got " << engines << "\n";
            }
            if (abName.empty()) {
                std::cerr << "[AFTB] empty model name\n";
            } else {
                std::string abPath = baseDir + abName + ".IFF";
                if (abPath == path) {
                    std::cerr << "[AFTB] refusing to load self-reference " << abPath << "\n";
                } else {
                    Model abModel;
                    std::vector<Model> abFrames;
                    if (!load_wc3_model_hcl_textured(abPath, abModel, &abFrames)) {
                        std::cerr << "[AFTB] failed to load " << abPath << "\n";
                    } else {
                        if (abFrames.empty()) {
                            std::cerr << "[AFTB] no frames in " << abPath << "\n";
                        } else {
                            if (abFrames.size() < 8) {
                                abFrames.resize(8);
                            }
                            M.afterburnerFrames.clear();
                            M.afterburnerFrames.reserve(abFrames.size());
                            for (size_t i = 0; i < abFrames.size(); ++i) {
                                const Model& srcFrame = abFrames[i];
                                Model combined;
                                combined.name = srcFrame.name;
                                combined.textures = srcFrame.textures;
                                if (M.afterburnerOffsets.empty()) {
                                combined.verts = srcFrame.verts;
                                combined.tris = srcFrame.tris;
                            } else {
                                for (const auto& offset : M.afterburnerOffsets) {
                                    size_t vBase = combined.verts.size();
                                    combined.verts.reserve(combined.verts.size() + srcFrame.verts.size());
                                    for (const auto& v : srcFrame.verts) {
                                        combined.verts.push_back({ v.x + offset.x, v.y + offset.y, v.z + offset.z });
                                    }
                                    combined.tris.reserve(combined.tris.size() + srcFrame.tris.size());
                                    for (const auto& tri : srcFrame.tris) {
                                        Tri t = tri;
                                        t.v[0] += (uint32_t)vBase;
                                        t.v[1] += (uint32_t)vBase;
                                        t.v[2] += (uint32_t)vBase;
                                        combined.tris.push_back(t);
                                    }
                                }
                            }
                            combined.afterburnerModelName.clear();
                            combined.afterburnerOffsets.clear();
                            combined.afterburnerFrames.clear();
                            M.afterburnerFrames.push_back(std::move(combined));
                        }
                        std::cerr << "[AFTB] model=" << abName
                                  << " frames=" << M.afterburnerFrames.size()
                                  << " engines=" << M.afterburnerOffsets.size()
                                  << "\n";
                    }
                }
            }
        }
    }

    flattenSubmodelsInto(M);
    return true;
}

// ----------------------------- GL batching -----------------------------
struct Vtx { float x,y,z,u,v; };
struct Batch { GLuint vao=0,vbo=0,ebo=0; GLsizei idxCount=0; uint16_t tex=0; };

static vector<Batch> buildBatches(const Model& M){
    std::unordered_map<uint16_t, vector<Vtx>> vbPerTex;
    std::unordered_map<uint16_t, vector<uint32_t>> ibPerTex;

    for(const auto& t: M.tris){
        uint16_t tx = 65535;
        if (t.hasTex && t.tex < M.textures.size()) {
            const Texture& tex = M.textures[t.tex];
            if (tex.skipRender) {
                continue; // omit BACK/FRONT polygons entirely
            }
            if (tex.valid()) {
                tx = t.tex;
            }
        }
        auto& vb = vbPerTex[tx];
        auto& ib = ibPerTex[tx];

        float invW=1.f, invH=1.f;
        if(tx != 65535 && tx < M.textures.size() && M.textures[tx].valid()){
            invW = 1.f / std::max(1, M.textures[tx].w);
            invH = 1.f / std::max(1, M.textures[tx].h);
        }

        uint32_t base = (uint32_t)vb.size();
        for(int i=0;i<3;i++){
            const Vec3& P = M.verts[t.v[i]];
            float u=0.f, v=0.f;
            if(t.hasTex){
                u =  t.uv[i*2+0] * invW;
		v = t.uv[i*2+1] * invH;
            }
            vb.push_back({P.x,P.y,P.z,u,v});
        }
        ib.push_back(base+0); ib.push_back(base+1); ib.push_back(base+2);
    }

    vector<Batch> out;
    for(auto& kv : vbPerTex){
        uint16_t tx = kv.first;
        auto& vb = kv.second; auto& ib = ibPerTex[tx];
        GLuint vao,vbo,ebo; glGenVertexArrays(1,&vao); glBindVertexArray(vao);
        glGenBuffers(1,&vbo); glBindBuffer(GL_ARRAY_BUFFER,vbo);
        glBufferData(GL_ARRAY_BUFFER, vb.size()*sizeof(Vtx), vb.data(), GL_STATIC_DRAW);
        glGenBuffers(1,&ebo); glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, ib.size()*sizeof(uint32_t), ib.data(), GL_STATIC_DRAW);
        glEnableVertexAttribArray(0); glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(Vtx),(void*)0);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,sizeof(Vtx),(void*)(3*sizeof(float)));
        glBindVertexArray(0);
        out.push_back(Batch{vao,vbo,ebo,(GLsizei)ib.size(), tx});
    }
    return out;
}

// ----------------------------- GL helpers -----------------------------
static GLuint makeShader(GLenum type, const char* src){
    GLuint s = glCreateShader(type);
    glShaderSource(s,1,&src,nullptr);
    glCompileShader(s);
    GLint ok=0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if(!ok){
        GLint len=0; glGetShaderiv(s, GL_INFO_LOG_LENGTH, &len);
        std::string log(len,'\0'); glGetShaderInfoLog(s, len, nullptr, log.data());
        std::cerr << "[glsl] " << (type==GL_VERTEX_SHADER?"VS":"FS") << ":\n" << log << "\n";
    }
    return s;
}
static GLuint makeProgram(const char* vs, const char* fs){
    GLuint v=makeShader(GL_VERTEX_SHADER, vs);
    GLuint f=makeShader(GL_FRAGMENT_SHADER, fs);
    GLuint p=glCreateProgram();
    glAttachShader(p,v); glAttachShader(p,f);
    glBindAttribLocation(p,0,"aPos");
    glBindAttribLocation(p,1,"aUV");
    glLinkProgram(p);
    glDeleteShader(v); glDeleteShader(f);
    return p;
}

static void uploadTexture(Texture& T){
    if(!T.valid()) return;
    glGenTextures(1,&T.gl);
    glBindTexture(GL_TEXTURE_2D, T.gl);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, T.w, T.h, 0, GL_RGBA, GL_UNSIGNED_BYTE, T.rgba.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

// ----------------------------- OBJ/MTL export -----------------------------
static bool exportOBJ(const Model& M, const std::string& outBase){
    if(M.verts.empty() || M.tris.empty()) { std::cerr<<"[export] nothing to export\n"; return false; }

    std::string objPath = outBase + ".obj";
    std::string mtlPath = outBase + ".mtl";

    std::ofstream obj(objPath);
    if(!obj){ std::cerr<<"[export] cannot write "<<objPath<<"\n"; return false; }
    std::ofstream mtl(mtlPath);
    if(!mtl){ std::cerr<<"[export] cannot write "<<mtlPath<<"\n"; return false; }

    obj << "mtllib " << mtlPath.substr(mtlPath.find_last_of("/\\")+1) << "\n";
    for(const auto& v : M.verts) obj << "v " << v.x << " " << v.y << " " << v.z << "\n";

    // write textures + MTL
    size_t tgaCount=0;
    std::string texDir;
    size_t slash = mtlPath.find_last_of("/\\");
    if(slash != std::string::npos) texDir = mtlPath.substr(0, slash+1);
    for(size_t i=0;i<M.textures.size();++i){
        const auto& T = M.textures[i];
        char defName[16]; std::snprintf(defName, sizeof(defName), "tex_%02zu", i);
        std::string texName = T.name.empty() ? std::string(defName) : T.name;
        std::string fileName = texName + ".tga";
        std::string fullPath = texDir + fileName;
        if(T.valid()){ writeTGA(fullPath, T.w, T.h, T.rgba.data()); ++tgaCount; }
        mtl << "newmtl m_tex" << i << "\n";
        mtl << "Kd 1 1 1\n";
        mtl << "map_Kd " << fileName << "\n\n";
    }
    mtl << "newmtl m_tex65535\nKd 1 1 1\n\n";

    // Sequentially write faces preserving group ranges for sub-objects (e.g. turrets)
    std::vector<ObjGroup> groups = M.groups;
    std::sort(groups.begin(), groups.end(),
              [](const ObjGroup& a, const ObjGroup& b){ return a.firstTri < b.firstTri; });

    // Start with main model name (hull) until first group
    std::string hullName = M.name.empty() ? std::string("hull") : M.name;
    obj << "o " << hullName << "\n";

    size_t groupIdx = 0;
    size_t nextGroupStart = groupIdx < groups.size() ? groups[groupIdx].firstTri : size_t(-1);
    size_t currentGroupEnd = nextGroupStart; // hull until first group
    std::string currentObj = hullName;

    uint16_t currentMtl = 65535;
    size_t vtBase = 1;
    for(size_t i=0;i<M.tris.size();++i){
        // Enter next group when reaching its first triangle
        if(i == nextGroupStart){
            currentObj = groups[groupIdx].name;
            obj << "o " << currentObj << "\n";
            currentGroupEnd = groups[groupIdx].firstTri + groups[groupIdx].triCount;
            ++groupIdx;
            nextGroupStart = groupIdx < groups.size() ? groups[groupIdx].firstTri : size_t(-1);
        }else if(i == currentGroupEnd){
            // Exit group back to hull until next group
            currentObj = hullName;
            obj << "o " << currentObj << "\n";
            currentGroupEnd = nextGroupStart;
        }

        uint16_t t = M.tris[i].tex;
        if(!(t < M.textures.size() && M.textures[t].valid())) t = 65535;
        if(t != currentMtl){
            obj << "usemtl m_tex" << t << "\n";
            currentMtl = t;
        }

        float invW=1.f, invH=1.f;
        if(t<M.textures.size() && M.textures[t].valid()){
            invW = 1.f / M.textures[t].w;
            invH = 1.f / M.textures[t].h;
        }

        const Tri& tr = M.tris[i];
        for(int k=0;k<3;k++){
            float u=0.f, v=0.f;
            if(tr.hasTex){
                u = tr.uv[k*2+0] * invW;
                v = 1.0f - tr.uv[k*2+1] * invH; // OBJ v-up
            }
            obj << "vt " << u << " " << v << "\n";
        }

        size_t vt0 = vtBase; vtBase += 3;
        obj << "f "
            << (tr.v[0]+1) << "/" << vt0   << " "
            << (tr.v[1]+1) << "/" << vt0+1 << " "
            << (tr.v[2]+1) << "/" << vt0+2 << "\n";
    }
    std::cerr<<"[export] wrote "<<tgaCount<<" TGAs + OBJ/MTL\n";
    return true;
}

// ----------------------------- draw helpers & matrices -----------------------------
struct Mat4 { float m[16]; static Mat4 identity(){ Mat4 r{}; for(int i=0;i<16;i++) r.m[i]=(i%5==0)?1.f:0.f; return r; } };
static Mat4 perspective(float fovy,float aspect,float zn,float zf){
    float f=1.f/std::tan(fovy*0.5f); Mat4 P{}; for(int i=0;i<16;i++) P.m[i]=0;
    P.m[0]=f/aspect; P.m[5]=f; P.m[10]=(zf+zn)/(zn-zf); P.m[11]=-1.f; P.m[14]=(2*zf*zn)/(zn-zf);
    return P;
}
static Mat4 lookAt(float ex,float ey,float ez, float cx,float cy,float cz, float ux,float uy,float uz){
    float fx=cx-ex, fy=cy-ey, fz=cz-ez; float fl=std::sqrt(fx*fx+fy*fy+fz*fz); fx/=fl; fy/=fl; fz/=fl;
    float sx = fy*uz - fz*uy; float sy = fz*ux - fx*uz; float sz = fx*uy - fy*ux;
    float sl=std::sqrt(sx*sx+sy*sy+sz*sz); sx/=sl; sy/=sl; sz/=sl;
    float ux2 = sy*fz - sz*fy; float uy2 = sz*fx - sx*fz; float uz2 = sx*fy - sy*fx;
    Mat4 M = Mat4::identity();
    M.m[0]=sx; M.m[4]=sy; M.m[8]=sz;
    M.m[1]=ux2; M.m[5]=uy2; M.m[9]=uz2;
    M.m[2]=-fx; M.m[6]=-fy; M.m[10]=-fz;
    M.m[12]=-(sx*ex + sy*ey + sz*ez);
    M.m[13]=-(ux2*ex + uy2*ey + uz2*ez);
    M.m[14]=(fx*ex + fy*ey + fz*ez);
    return M;
}
static Mat4 mul(const Mat4& A,const Mat4& B){
    Mat4 R{}; for(int r=0;r<4;r++)for(int c=0;c<4;c++){ float s=0; for(int k=0;k<4;k++) s+=A.m[k*4+r]*B.m[c*4+k]; R.m[c*4+r]=s; } return R;
}
static Mat4 translate(float x,float y,float z){
    Mat4 M = Mat4::identity();
    M.m[12] = x; M.m[13] = y; M.m[14] = z;
    return M;
}
static Mat4 scale1(float s){
    Mat4 M = Mat4::identity();
    M.m[0]=M.m[5]=M.m[10]=s;
    return M;
}

static Mat4 invertRigid(const Mat4& M){
    Mat4 R = Mat4::identity();
    float r00 = M.m[0],  r01 = M.m[4],  r02 = M.m[8];
    float r10 = M.m[1],  r11 = M.m[5],  r12 = M.m[9];
    float r20 = M.m[2],  r21 = M.m[6],  r22 = M.m[10];
    float tx  = M.m[12], ty  = M.m[13], tz  = M.m[14];

    // Transpose rotation
    R.m[0] = r00; R.m[1] = r01; R.m[2] = r02;
    R.m[4] = r10; R.m[5] = r11; R.m[6] = r12;
    R.m[8] = r20; R.m[9] = r21; R.m[10] = r22;

    // Inverted translation
    R.m[12] = -(r00 * tx + r10 * ty + r20 * tz);
    R.m[13] = -(r01 * tx + r11 * ty + r21 * tz);
    R.m[14] = -(r02 * tx + r12 * ty + r22 * tz);
    return R;
}

#if WC_HAVE_OPENVR
static Mat4 fromHmdMatrix34(const vr::HmdMatrix34_t& mat){
    Mat4 M = Mat4::identity();
    for(int r=0;r<3;r++){
        for(int c=0;c<4;c++){
            M.m[c*4 + r] = mat.m[r][c];
        }
    }
    return M;
}

static Mat4 fromHmdMatrix44(const vr::HmdMatrix44_t& mat){
    Mat4 M = Mat4::identity();
    for(int r=0;r<4;r++){
        for(int c=0;c<4;c++){
            M.m[c*4 + r] = mat.m[r][c];
        }
    }
    return M;
}

struct VRContext {
    vr::IVRSystem* system = nullptr;
    vr::TrackedDevicePose_t poses[vr::k_unMaxTrackedDeviceCount]{};
    Mat4 eyeProj[2];
    Mat4 eyeToHead[2];
    Mat4 hmdPose = Mat4::identity();
    GLuint eyeFBO[2]{};
    GLuint eyeColor[2]{};
    GLuint eyeDepth[2]{};
    uint32_t eyeWidth = 0;
    uint32_t eyeHeight = 0;

    bool init(){
        vr::EVRInitError err = vr::VRInitError_None;
        system = vr::VR_Init(&err, vr::VRApplication_Scene);
        if(err != vr::VRInitError_None){
            std::cerr << "[vr] SteamVR init failed: " << vr::VR_GetVRInitErrorAsEnglishDescription(err) << "\n";
            system = nullptr;
            return false;
        }
        if(!vr::VRCompositor()){
            std::cerr << "[vr] No VR compositor available\n";
            vr::VR_Shutdown();
            system = nullptr;
            return false;
        }
        system->GetRecommendedRenderTargetSize(&eyeWidth, &eyeHeight);
        updateEyeData();
        if(!createTargets()){
            shutdown();
            return false;
        }
        return true;
    }

    void shutdown(){
        destroyTargets();
        if(system){
            vr::VR_Shutdown();
            system = nullptr;
        }
    }

    void updateEyeData(){
        if(!system) return;
        for(int eye=0; eye<2; ++eye){
            eyeProj[eye] = fromHmdMatrix44(system->GetProjectionMatrix(static_cast<vr::Hmd_Eye>(eye), 0.05f, 100.0f));
            eyeToHead[eye] = fromHmdMatrix34(system->GetEyeToHeadTransform(static_cast<vr::Hmd_Eye>(eye)));
        }
    }

    bool beginFrame(){
        if(!system) return false;
        vr::VRCompositor()->WaitGetPoses(poses, vr::k_unMaxTrackedDeviceCount, nullptr, 0);
        if(poses[vr::k_unTrackedDeviceIndex_Hmd].bPoseIsValid){
            hmdPose = fromHmdMatrix34(poses[vr::k_unTrackedDeviceIndex_Hmd].mDeviceToAbsoluteTracking);
        }else{
            hmdPose = Mat4::identity();
        }
        return true;
    }

    Mat4 viewForEye(int eyeIndex, const Mat4& baseView) const{
        Mat4 headInv = invertRigid(hmdPose);
        Mat4 eyeInv = invertRigid(eyeToHead[eyeIndex]);
        return mul(eyeInv, mul(headInv, baseView));
    }

    const Mat4& projectionForEye(int eyeIndex) const{ return eyeProj[eyeIndex]; }

    void submit(){
        for(int eye=0; eye<2; ++eye){
            vr::Texture_t tex{};
            tex.handle = reinterpret_cast<void*>(static_cast<uintptr_t>(eyeColor[eye]));
            tex.eType = vr::TextureType_OpenGL;
            tex.eColorSpace = vr::ColorSpace_Gamma;
            vr::VRCompositor()->Submit(static_cast<vr::Hmd_Eye>(eye), &tex);
        }
        vr::VRCompositor()->PostPresentHandoff();
    }

private:
    bool createTargets(){
        glGenFramebuffers(2, eyeFBO);
        glGenTextures(2, eyeColor);
        glGenRenderbuffers(2, eyeDepth);
        for(int eye=0; eye<2; ++eye){
            glBindTexture(GL_TEXTURE_2D, eyeColor[eye]);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, eyeWidth, eyeHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

            glBindFramebuffer(GL_FRAMEBUFFER, eyeFBO[eye]);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, eyeColor[eye], 0);

            glBindRenderbuffer(GL_RENDERBUFFER, eyeDepth[eye]);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, eyeWidth, eyeHeight);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, eyeDepth[eye]);

            GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
            if(status != GL_FRAMEBUFFER_COMPLETE){
                std::cerr << "[vr] framebuffer incomplete for eye " << eye << " status=" << std::hex << status << std::dec << "\n";
                glBindFramebuffer(GL_FRAMEBUFFER, 0);
                glBindRenderbuffer(GL_RENDERBUFFER, 0);
                glBindTexture(GL_TEXTURE_2D, 0);
                return false;
            }
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        glBindTexture(GL_TEXTURE_2D, 0);
        return true;
    }

    void destroyTargets(){
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        glBindTexture(GL_TEXTURE_2D, 0);
        glDeleteFramebuffers(2, eyeFBO);
        glDeleteTextures(2, eyeColor);
        glDeleteRenderbuffers(2, eyeDepth);
        for(int eye=0; eye<2; ++eye){
            eyeFBO[eye] = eyeColor[eye] = eyeDepth[eye] = 0;
        }
    }
};
#endif

static void drawBBox(const Vec3& bmin, const Vec3& bmax, GLuint prog, GLint uMVP, const Mat4& MVP, const std::array<float,3>& color){
    struct BBoxVertex{ float x,y,z,r,g,b; };
    const unsigned idx[]={0,1,1,2,2,3,3,0, 4,5,5,6,6,7,7,4, 0,4,1,5,2,6,3,7};
    BBoxVertex verts[8];
    const float x0=bmin.x, y0=bmin.y, z0=bmin.z;
    const float x1=bmax.x, y1=bmax.y, z1=bmax.z;
    const float cr=color[0], cg=color[1], cb=color[2];
    verts[0]={x0,y0,z0, cr,cg,cb}; verts[1]={x1,y0,z0, cr,cg,cb};
    verts[2]={x1,y1,z0, cr,cg,cb}; verts[3]={x0,y1,z0, cr,cg,cb};
    verts[4]={x0,y0,z1, cr,cg,cb}; verts[5]={x1,y0,z1, cr,cg,cb};
    verts[6]={x1,y1,z1, cr,cg,cb}; verts[7]={x0,y1,z1, cr,cg,cb};

    GLuint vao=0,vbo=0,ebo=0;
    glGenVertexArrays(1,&vao);
    glGenBuffers(1,&vbo);
    glGenBuffers(1,&ebo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER,vbo);
    glBufferData(GL_ARRAY_BUFFER,sizeof(verts),verts,GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(BBoxVertex),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(BBoxVertex),(void*)(3*sizeof(float)));
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(idx),idx,GL_DYNAMIC_DRAW);
    glUseProgram(prog);
    glUniformMatrix4fv(uMVP,1,GL_FALSE,MVP.m);
    glDrawElements(GL_LINES,(GLsizei)(sizeof(idx)/sizeof(idx[0])),GL_UNSIGNED_INT,0);
    glBindVertexArray(0);
    glDeleteBuffers(1,&ebo);
    glDeleteBuffers(1,&vbo);
    glDeleteVertexArrays(1,&vao);
}

// ----------------------------- main -----------------------------
int main(int argc, char** argv){
    if(argc<2){
        std::cerr<<"Usage: "<<argv[0]<<" FILE.IFF [--out BASE] [--export] [--export-only] [--no-fit] [--vr]\n";
        return 1;
    }
    string path = argv[1];
    string outBase = stemFromPath(path);

    // Parse flags up-front so palette is ready BEFORE decoding textures
    bool exportOnStart=false; bool exportOnly=false; bool noFit=false; bool useVR=false;
    for(int i=2;i<argc;i++){
        string a=argv[i];
        if(a=="--out" && i+1<argc) outBase = argv[++i];
        else if(a=="--export") exportOnStart = true;
        else if(a=="--export-only") exportOnly = true;
        else if(a=="--no-fit") noFit = true;
        else if(a=="--vr") useVR = true;
    }

#if !WC_HAVE_OPENVR
    if(useVR){
        std::cerr << "[vr] SteamVR support not available in this build\n";
        return 4;
    }
#endif

    std::vector<std::string> paletteFiles = {"wc3pal.json","wc4pal.json","armpal.json"};
    std::vector<std::array<std::array<uint8_t,3>,256>> loadedPalettes;
    for(const auto& f : paletteFiles){
        if(loadPaletteJSON(f,0)){
            std::array<std::array<uint8_t,3>,256> pal{};
            for(int i=0;i<256;i++) pal[i]=gPalette[i];
            loadedPalettes.push_back(pal);
        }else{
            std::array<std::array<uint8_t,3>,256> pal{};
            for(int i=0;i<256;i++){ pal[i][0]=pal[i][1]=pal[i][2]=(uint8_t)i; }
            loadedPalettes.push_back(pal);
        }
    }
    int currentPalette = 0;
    if(!loadedPalettes.empty()){
        for(int i=0;i<256;i++) gPalette[i] = loadedPalettes[currentPalette][i];
    }else{
        makeDefaultPalette();
    }

    Model M;
    if(!load_wc3_model_hcl_textured(path, M)){ std::cerr<<"Failed to load model\n"; return 2; }

    if(exportOnStart || exportOnly){
        std::cerr<<"[export] writing "<<outBase<<".obj/.mtl + TGAs\n";
        exportOBJ(M, outBase);
        if(exportOnly) return 0;
    }

    if(SDL_Init(SDL_INIT_VIDEO)!=0){ std::cerr<<"SDL init failed: "<<SDL_GetError()<<"\n"; return 3; }
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION,3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION,3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_Window* win = SDL_CreateWindow(("WC3 Viewer - "+M.name).c_str(),
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 800,
        SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
    SDL_GLContext glc = SDL_GL_CreateContext(win);
    glewExperimental=GL_TRUE; glewInit();
#if WC_HAVE_OPENVR
    VRContext vr;
    if(useVR){
        if(!vr.init()){
            std::cerr << "[vr] initialization failed, continuing without VR\n";
            useVR = false;
        }
    }
#endif
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    auto shutdownAndExit = [&](int code){
#if WC_HAVE_OPENVR
        if(useVR) vr.shutdown();
#endif
        SDL_GL_DeleteContext(glc);
        SDL_DestroyWindow(win);
        SDL_Quit();
        return code;
    };

    // Upload textures (+ fallback white)
    for (auto& T : M.textures) {
        if (T.skipRender) continue;
        uploadTexture(T);
    }
    Texture white; white.w=1; white.h=1; white.rgba={255,255,255,255}; uploadTexture(white);

    // Build batches
    auto batches = buildBatches(M);
    std::vector<std::vector<Batch>> afterburnerBatches;
    if (!M.afterburnerFrames.empty()) {
        afterburnerBatches.resize(M.afterburnerFrames.size());
        for (size_t i = 0; i < M.afterburnerFrames.size(); ++i) {
            Model& frameModel = M.afterburnerFrames[i];
            for (auto& T : frameModel.textures) {
                if (T.skipRender) continue;
                uploadTexture(T);
            }
            afterburnerBatches[i] = buildBatches(frameModel);
        }
    }
    int currentAfterburnerFrame = 0;

    // Shaders
    const char* VS =
        "#version 330 core\n"
        "layout(location=0) in vec3 aPos;\n"
        "layout(location=1) in vec2 aUV;\n"
        "uniform mat4 uMVP;\n"
        "out vec2 vUV;\n"
        "void main(){ vUV=aUV; gl_Position = uMVP * vec4(aPos,1.0); }\n";
    const char* FS =
        "#version 330 core\n"
        "in vec2 vUV; uniform sampler2D uTex; out vec4 o;\n"
        "void main(){ vec4 c = texture(uTex, vUV); if(c.a < 0.5) discard; o = c; }\n";
    GLuint prog = makeProgram(VS,FS);
    GLint uMVP = glGetUniformLocation(prog,"uMVP");
    GLint uTex = glGetUniformLocation(prog,"uTex");

    const char* axisVS =
        "#version 330 core\n"
        "layout(location=0) in vec3 aPos;\n"
        "layout(location=1) in vec3 aColor;\n"
        "uniform mat4 uMVP;\n"
        "out vec3 vColor;\n"
        "void main(){ vColor=aColor; gl_Position = uMVP * vec4(aPos,1.0); }\n";
    const char* axisFS =
        "#version 330 core\n"
        "in vec3 vColor; out vec4 o;\n"
        "void main(){ o = vec4(vColor, 1.0); }\n";
    GLuint axisV = makeShader(GL_VERTEX_SHADER, axisVS);
    GLuint axisF = makeShader(GL_FRAGMENT_SHADER, axisFS);
    GLuint axisProg = glCreateProgram();
    glAttachShader(axisProg, axisV);
    glAttachShader(axisProg, axisF);
    glBindAttribLocation(axisProg, 0, "aPos");
    glBindAttribLocation(axisProg, 1, "aColor");
    glLinkProgram(axisProg);
    glDeleteShader(axisV);
    glDeleteShader(axisF);
    GLint axisUMVP = glGetUniformLocation(axisProg, "uMVP");

    struct AxisVertex { float x,y,z,r,g,b; };
    const float axisTip = 1.0f;
    const float axisBase = 0.82f;
    const float head = 0.12f;
    const AxisVertex axisVerts[] = {
        {0.f,0.f,0.f, 1.f,0.f,0.f}, {axisTip,0.f,0.f, 1.f,0.f,0.f},
        {axisTip,0.f,0.f, 1.f,0.f,0.f}, {axisBase, head, head, 1.f,0.f,0.f},
        {axisTip,0.f,0.f, 1.f,0.f,0.f}, {axisBase,-head, head, 1.f,0.f,0.f},
        {axisTip,0.f,0.f, 1.f,0.f,0.f}, {axisBase, head,-head, 1.f,0.f,0.f},
        {axisTip,0.f,0.f, 1.f,0.f,0.f}, {axisBase,-head,-head, 1.f,0.f,0.f},
        {0.f,0.f,0.f, 0.f,1.f,0.f}, {0.f,axisTip,0.f, 0.f,1.f,0.f},
        {0.f,axisTip,0.f, 0.f,1.f,0.f}, { head,axisBase, head, 0.f,1.f,0.f},
        {0.f,axisTip,0.f, 0.f,1.f,0.f}, {-head,axisBase, head, 0.f,1.f,0.f},
        {0.f,axisTip,0.f, 0.f,1.f,0.f}, { head,axisBase,-head, 0.f,1.f,0.f},
        {0.f,axisTip,0.f, 0.f,1.f,0.f}, {-head,axisBase,-head, 0.f,1.f,0.f},
        {0.f,0.f,0.f, 0.f,0.f,1.f}, {0.f,0.f,axisTip, 0.f,0.f,1.f},
        {0.f,0.f,axisTip, 0.f,0.f,1.f}, { head, head,axisBase, 0.f,0.f,1.f},
        {0.f,0.f,axisTip, 0.f,0.f,1.f}, {-head, head,axisBase, 0.f,0.f,1.f},
        {0.f,0.f,axisTip, 0.f,0.f,1.f}, { head,-head,axisBase, 0.f,0.f,1.f},
        {0.f,0.f,axisTip, 0.f,0.f,1.f}, {-head,-head,axisBase, 0.f,0.f,1.f},
    };
    GLuint axisVAO=0, axisVBO=0;
    glGenVertexArrays(1,&axisVAO);
    glGenBuffers(1,&axisVBO);
    glBindVertexArray(axisVAO);
    glBindBuffer(GL_ARRAY_BUFFER, axisVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(axisVerts), axisVerts, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(AxisVertex),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(AxisVertex),(void*)(3*sizeof(float)));
    glBindVertexArray(0);
    const GLsizei axisVertCount = (GLsizei)(sizeof(axisVerts)/sizeof(axisVerts[0]));

    const char* reticleVS =
        "#version 330 core\n"
        "layout(location=0) in vec2 aPos;\n"
        "void main(){ gl_Position = vec4(aPos, 0.0, 1.0); }\n";
    const char* reticleFS =
        "#version 330 core\n"
        "uniform vec3 uColor; out vec4 o;\n"
        "void main(){ o = vec4(uColor, 1.0); }\n";
    GLuint reticleProg = makeProgram(reticleVS, reticleFS);
    GLint reticleColorLoc = glGetUniformLocation(reticleProg, "uColor");

    std::vector<float> reticleCircle;
    const int circleSegments = 64;
    const float reticleRadius = 0.04f;
    const float pi = 3.14159265358979323846f;
    for (int i = 0;i < circleSegments;i++) {
        float angle = (float)i / (float)circleSegments * 2.0f * pi;
        reticleCircle.push_back(std::cos(angle) * reticleRadius);
        reticleCircle.push_back(std::sin(angle) * reticleRadius);
    }
    const GLsizei reticleCircleCount = (GLsizei)(reticleCircle.size() / 2);
    GLuint reticleCircleVAO = 0, reticleCircleVBO = 0;
    glGenVertexArrays(1, &reticleCircleVAO);
    glGenBuffers(1, &reticleCircleVBO);
    glBindVertexArray(reticleCircleVAO);
    glBindBuffer(GL_ARRAY_BUFFER, reticleCircleVBO);
    glBufferData(GL_ARRAY_BUFFER, reticleCircle.size() * sizeof(float), reticleCircle.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);

    const float crossHalf = reticleRadius * 0.8f;
    const float crossVerts[] = {
        -crossHalf, 0.0f,  crossHalf, 0.0f,
         0.0f, -crossHalf, 0.0f,  crossHalf
    };
    GLuint reticleCrossVAO = 0, reticleCrossVBO = 0;
    glGenVertexArrays(1, &reticleCrossVAO);
    glGenBuffers(1, &reticleCrossVBO);
    glBindVertexArray(reticleCrossVAO);
    glBindBuffer(GL_ARRAY_BUFFER, reticleCrossVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(crossVerts), crossVerts, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glBindVertexArray(0);

    // Compute bbox
    Vec3 bmin{+1e9f,+1e9f,+1e9f}, bmax{-1e9f,-1e9f,-1e9f};
    auto expandBounds = [&](const Vec3& v){
        bmin.x = std::min(bmin.x, v.x);
        bmin.y = std::min(bmin.y, v.y);
        bmin.z = std::min(bmin.z, v.z);
        bmax.x = std::max(bmax.x, v.x);
        bmax.y = std::max(bmax.y, v.y);
        bmax.z = std::max(bmax.z, v.z);
    };
    for(const auto& v: M.verts){ expandBounds(v); }
    for (const auto& frame : M.afterburnerFrames) {
        for (const auto& v : frame.verts) {
            expandBounds(v);
        }
    }

    // Auto-fit parameters (draw only)
    Vec3 center{
        (bmin.x + bmax.x)*0.5f,
        (bmin.y + bmax.y)*0.5f,
        (bmin.z + bmax.z)*0.5f
    };
    float dx = (bmax.x - bmin.x);
    float dy = (bmax.y - bmin.y);
    float dz = (bmax.z - bmin.z);
    float longest = std::max({dx, dy, dz});
    if (longest <= 1e-6f) longest = 1.0f;
    float fitScale = noFit ? 1.0f : (2.0f / longest); // longest side -> ~2 units

    float sceneSpan = std::max({ bmax.x - bmin.x, bmax.y - bmin.y, bmax.z - bmin.z });
    float yaw   = 0.8f;          // radians
    float pitch = 0.35f;         // radians
    bool  rotating = false;

    // Initial distance (fit vs no-fit)
    float camDist = noFit ? (sceneSpan * 2.2f + 1.0f) : 4.0f;
    // Reasonable clamps
    float minDist = noFit ? std::max(0.1f, sceneSpan * 0.05f) : 0.2f;
    float maxDist = noFit ? sceneSpan * 50.0f : 50.0f;

    Vec3 camPos{};
    {
        float cp = std::cos(pitch), sp = std::sin(pitch);
        float cy = std::cos(yaw), sy = std::sin(yaw);
        camPos.x = cp * cy * camDist;
        camPos.y = sp * camDist;
        camPos.z = cp * sy * camDist;
    }
    bool showReticle = false;

    // Simple orbit camera + bbox
    bool drawBox=false;
    bool wireframe=false;
    bool showAxis=false;
    uint32_t last = SDL_GetTicks();
    float t=0.f; int winW=1280, winH=800;

    auto drawModelBatches = [&](const Model& model, const std::vector<Batch>& bs){
        for(const auto& b : bs){
            GLuint tex = (b.tex==65535) ? white.gl :
                         (b.tex<model.textures.size() && model.textures[b.tex].gl ? model.textures[b.tex].gl : white.gl);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, tex);
            glBindVertexArray(b.vao);
            glDrawElements(GL_TRIANGLES, b.idxCount, GL_UNSIGNED_INT, 0);
        }
    };

    auto renderScene = [&](const Mat4& MVP){
        glUseProgram(prog);
        glUniformMatrix4fv(uMVP,1,GL_FALSE,MVP.m);
        glUniform1i(uTex, 0);

        drawModelBatches(M, batches);
        if (!afterburnerBatches.empty()) {
            size_t idx = (size_t)std::clamp(currentAfterburnerFrame, 0, (int)afterburnerBatches.size() - 1);
            if (idx < afterburnerBatches.size()) {
                drawModelBatches(M.afterburnerFrames[idx], afterburnerBatches[idx]);
            }
        }
        glBindVertexArray(0);

        if(drawBox){
            drawBBox(bmin, bmax, axisProg, axisUMVP, MVP, std::array<float,3>{1.f, 1.f, 0.f});
            glUseProgram(prog);
        }
        if (showReticle) {
            glDisable(GL_DEPTH_TEST);
            glUseProgram(reticleProg);
            glUniform3f(reticleColorLoc, 1.0f, 1.0f, 1.0f);
            glBindVertexArray(reticleCircleVAO);
            glLineWidth(1.5f);
            glDrawArrays(GL_LINE_LOOP, 0, reticleCircleCount);
            glBindVertexArray(reticleCrossVAO);
            glDrawArrays(GL_LINES, 0, 4);
            glLineWidth(1.0f);
            glBindVertexArray(0);
            glEnable(GL_DEPTH_TEST);
            glUseProgram(prog);
        }
    };

    auto renderAxisOverlay = [&](const Mat4& view){
        if(!showAxis) return;
        GLint prevViewport[4];
        glGetIntegerv(GL_VIEWPORT, prevViewport);
        int axisBase = std::min(winW, winH);
        int axisSize = std::max(80, axisBase / 5);
        axisSize = std::min(axisSize, axisBase);
        glViewport(0, 0, axisSize, axisSize);
        glDisable(GL_DEPTH_TEST);

        Mat4 axisRot = view;
        axisRot.m[12] = axisRot.m[13] = axisRot.m[14] = 0.0f;
        Mat4 axisProj = perspective(0.7f, 1.0f, 0.01f, 10.0f);
        Mat4 axisMV = mul(translate(0.0f, 0.0f, -1.5f), axisRot);
        Mat4 axisMVP = mul(axisProj, axisMV);

        glUseProgram(axisProg);
        glUniformMatrix4fv(axisUMVP, 1, GL_FALSE, axisMVP.m);
        glBindVertexArray(axisVAO);
        glLineWidth(2.0f);
        glDrawArrays(GL_LINES, 0, axisVertCount);
        glLineWidth(1.0f);
        glBindVertexArray(0);
        glEnable(GL_DEPTH_TEST);
        glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);
        glUseProgram(prog);
    };

    while(true){
        uint32_t now = SDL_GetTicks();
        float dt = (now - last) * 0.001f;
        last = now;
        t += dt;

        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) { return shutdownAndExit(0); }
            if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
                winW = e.window.data1; winH = e.window.data2; glViewport(0, 0, winW, winH);
            }
            if (e.type == SDL_KEYDOWN) {
                SDL_Scancode sc = e.key.keysym.scancode;
                if (!M.afterburnerFrames.empty()) {
                    int newFrame = -1;
                    if (sc >= SDL_SCANCODE_0 && sc <= SDL_SCANCODE_7) {
                        newFrame = sc - SDL_SCANCODE_0;
                    } else if (sc >= SDL_SCANCODE_KP_0 && sc <= SDL_SCANCODE_KP_7) {
                        newFrame = sc - SDL_SCANCODE_KP_0;
                    }
                    if (newFrame >= 0 && newFrame < (int)M.afterburnerFrames.size()) {
                        currentAfterburnerFrame = newFrame;
                        std::cerr << "[AFTB] frame=" << currentAfterburnerFrame << "\n";
                    }
                }
                if (sc == SDL_SCANCODE_ESCAPE) { return shutdownAndExit(0); }
                if (sc == SDL_SCANCODE_O || sc == SDL_SCANCODE_E) { std::cerr << "[export] key -> OBJ/MTL/TGA\n"; exportOBJ(M, outBase); }
                if (sc == SDL_SCANCODE_P) {
                    currentPalette = (currentPalette + 1) % loadedPalettes.size();
                    std::cerr << "[pal] switching to " << paletteFiles[currentPalette] << "\n";
                    for (int i = 0;i < 256;i++) gPalette[i] = loadedPalettes[currentPalette][i];
                    for (auto& T : M.textures) {
                        if (!T.idx.empty()) {
                            palToRGBA(T.idx, T.w, T.h, T.rgba);
                            glBindTexture(GL_TEXTURE_2D, T.gl);
                            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, T.w, T.h, GL_RGBA, GL_UNSIGNED_BYTE, T.rgba.data());
                        }
                    }
                    glBindTexture(GL_TEXTURE_2D, 0);
                }
                if (sc == SDL_SCANCODE_B) { drawBox = !drawBox; }
                if (sc == SDL_SCANCODE_W) {
                    wireframe = !wireframe;
                    glPolygonMode(GL_FRONT_AND_BACK, wireframe ? GL_LINE : GL_FILL);
                }
                if (sc == SDL_SCANCODE_X) {
                    showAxis = !showAxis;
                }
                if (sc == SDL_SCANCODE_I) {
                    showReticle = !showReticle;
                }
            }
            if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                rotating = true;
                SDL_SetRelativeMouseMode(SDL_TRUE); // grab mouse while rotating
            }
            if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
                rotating = false;
                SDL_SetRelativeMouseMode(SDL_FALSE);
            }
            if (e.type == SDL_MOUSEMOTION && rotating) {
                const float sens = 0.005f; // radians per pixel
                yaw += e.motion.xrel * sens;
                pitch += -e.motion.yrel * sens;
                pitch = std::clamp(pitch, -1.5f, 1.5f); // avoid flipping over
            }
            if (e.type == SDL_MOUSEWHEEL) {
                if (e.wheel.y != 0) {
                    // zoom in on positive y, out on negative; exponential feels nicer
                    float oldDist = vecLength(camPos);
                    if (oldDist < 1e-6f) oldDist = 1e-6f;
                    float newDist = std::clamp(oldDist * std::pow(0.9f, (float)e.wheel.y), minDist, maxDist);
                    float scale = newDist / oldDist;
                    camPos = vecScale(camPos, scale);
                    camDist = newDist;
                }
            }
            // optional: R to reset view
            if (e.type == SDL_KEYDOWN && e.key.keysym.scancode == SDL_SCANCODE_R) {
                yaw = 0.8f; pitch = 0.35f;
                camDist = noFit ? (sceneSpan * 2.2f + 1.0f) : 4.0f;
            }
        }
        pitch = std::clamp(pitch, -1.5f, 1.5f);

        const Uint8* keys = SDL_GetKeyboardState(nullptr);
        const float turnSpeed = 1.5f;
        if (keys[SDL_SCANCODE_LEFT]) yaw -= turnSpeed * dt;
        if (keys[SDL_SCANCODE_RIGHT]) yaw += turnSpeed * dt;
        if (keys[SDL_SCANCODE_UP]) pitch -= turnSpeed * dt;
        if (keys[SDL_SCANCODE_DOWN]) pitch += turnSpeed * dt;

        float cp = std::cos(pitch), sp = std::sin(pitch);
        float cy = std::cos(yaw), sy = std::sin(yaw);
        Vec3 forward{ cp * cy, sp, cp * sy };
        Vec3 worldUp{ 0.0f, 1.0f, 0.0f };
        Vec3 right = vecCross(worldUp, forward);
        float rightLen = vecLength(right);
        if (rightLen > 1e-6f) {
            right = vecScale(right, 1.0f / rightLen);
        }

        const float moveSpeed = std::max(sceneSpan * 0.5f, 0.5f);
        float moveInput = 0.0f;
        if (keys[SDL_SCANCODE_A]) moveInput += 0.01f;
        if (keys[SDL_SCANCODE_Z]) moveInput -= 0.01f;
        if (std::abs(moveInput) > 0.0f) {
            camPos = vecAdd(camPos, vecScale(forward, moveInput * moveSpeed * dt));
        }

        float truckInput = 0.0f;
        if (keys[SDL_SCANCODE_KP_4]) truckInput -= 0.01f;
        if (keys[SDL_SCANCODE_KP_6]) truckInput += 0.01f;
        if (std::abs(truckInput) > 0.0f && rightLen > 1e-6f) {
            camPos = vecAdd(camPos, vecScale(right, truckInput * moveSpeed * dt));
        }

        float jibInput = 0.0f;
        if (keys[SDL_SCANCODE_KP_8]) jibInput += 0.01f;
        if (keys[SDL_SCANCODE_KP_2]) jibInput -= 0.01f;
        if (std::abs(jibInput) > 0.0f) {
            camPos = vecAdd(camPos, vecScale(worldUp, jibInput * moveSpeed * dt));
        }

        camDist = vecLength(camPos);
        if (camDist < minDist || camDist > maxDist) {
            Vec3 dir = vecNormalize(camPos);
            if (vecLength(dir) <= 1e-6f) {
                dir = vecNormalize(forward);
            }
            float target = std::clamp(camDist, minDist, maxDist);
            camPos = vecScale(dir, target);
            camDist = target;
        }

        Mat4 V = lookAt(camPos.x, camPos.y, camPos.z,
            camPos.x + forward.x, camPos.y + forward.y, camPos.z + forward.z,
            0, 1, 0);
        // Model matrix: translate to origin then scale to fit
        Mat4 Mdl = noFit ? translate(0,0,0) : mul(scale1(fitScale), translate(-center.x, -center.y, -center.z));

        #if WC_HAVE_OPENVR
        if(useVR){
            vr.beginFrame();
            for(int eye=0; eye<2; ++eye){
                glBindFramebuffer(GL_FRAMEBUFFER, vr.eyeFBO[eye]);
                glViewport(0, 0, vr.eyeWidth, vr.eyeHeight);
                glClearColor(0.08f,0.08f,0.10f,1);
                glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
                Mat4 eyeView = vr.viewForEye(eye, V);
                Mat4 MVP = mul(vr.projectionForEye(eye), mul(eyeView, Mdl));
                renderScene(MVP);
            }
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            glViewport(0,0,winW,winH);
            glBindFramebuffer(GL_READ_FRAMEBUFFER, vr.eyeFBO[0]);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glBlitFramebuffer(0,0,vr.eyeWidth,vr.eyeHeight, 0,0, winW, winH, GL_COLOR_BUFFER_BIT, GL_LINEAR);
            glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glEnable(GL_DEPTH_TEST);
            vr.submit();
            SDL_GL_SwapWindow(win);
            continue;
        }
        #endif

        float aspect = (float)winW/std::max(1,winH);
        float fovy = 1.0f; // ~57 deg
        Mat4 P = perspective(fovy, aspect, 0.05f, 100.0f);

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0,0,winW,winH);
        glClearColor(0.08f,0.08f,0.10f,1);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

        Mat4 MVP = mul(P, mul(V, Mdl));
        renderScene(MVP);
        renderAxisOverlay(V);
        SDL_GL_SwapWindow(win);
    }
    return shutdownAndExit(0);
}
