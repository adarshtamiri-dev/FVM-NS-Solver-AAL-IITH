
//#include "Mesh.h" // Base declarations and object blueprints
#include <iostream> // Output stream for logging execution milestones

#include "preproc/PreProcessor.h"
#include "Mesh.h"

int main() {
    Mesh mesh;
    PreProcessor test_pre;

    // Instantiating the class executes the entire pipeline automatically
    test_pre.CGNSReader("../imports/file2.cgns", mesh);
    test_pre.GenerateFaceData(mesh);
    test_pre.ComputeProperties(mesh);
    test_pre.VerifyMesh(mesh);

    return 0;
}