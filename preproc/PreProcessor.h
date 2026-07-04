#pragma once
#include <string>
#include "../Mesh.h"    // Contains your SoA Mesh struct

class PreProcessor {
public:
    PreProcessor();
    ~PreProcessor();

    bool CGNSReader(const std::string& filename, Mesh& mesh);
    void GenerateFaceData(Mesh& mesh);
    void ComputeProperties(Mesh& mesh);
    void VerifyMesh(const Mesh& mesh);
};
