#include "tinyply.h"

#include <vector>
#include <string>
#include <fstream>
#include <cstring>

struct SimpleMesh{
    /* A mesh where we only care about the vertices and the faces.
     * All the faces are triangles
     * */
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<size_t>> faces;
};

SimpleMesh readMesh(std::string fname){

    std::ifstream ss(fname);

    tinyply::PlyFile file;
    file.parse_header(ss);

    SimpleMesh m;
    m.vertices.clear();
    m.faces.clear();

    std::shared_ptr<tinyply::PlyData> vertices,faces;

    const size_t numVerticesBytes = vertices->buffer.size_bytes();
    struct float3 { float x, y, z; };
    struct double3 { double x, y, z; };

    std::vector<float3> verts_floats(vertices->count);
    std::vector<double3> verts_doubles(vertices->count);

    vertices = file.request_properties_from_element("vertex", { "x", "y", "z" });

    if (vertices->t == tinyply::Type::FLOAT32) { 
        std::memcpy(verts_floats.data(), vertices->buffer.get(), numVerticesBytes);
        for (float3 v : verts_floats)
            m.vertices.push_back({(double)v.x,(double)v.y,(double)v.z});
    }
    if (vertices->t == tinyply::Type::FLOAT64) {
        std::memcpy(verts_doubles.data(), vertices->buffer.get(), numVerticesBytes);
        for (double3 v : verts_doubles)
            m.vertices.push_back({v.x,v.y,v.z});
    }


    faces = file.request_properties_from_element("face", { "vertex_indices" });
    const size_t numFacesBytes = faces-> buffer.size_bytes();
    struct ind3 { int32_t a,b,c; };
    std::vector<ind3> faceVerts(faces->count);
    std::memcpy(faceVerts.data(),faces->buffer.get(),numFacesBytes);

    for (ind3 f : faceVerts){
        m.faces.push_back({(size_t)f.a,(size_t)f.b,(size_t)f.c});
    }

    return m;
};
