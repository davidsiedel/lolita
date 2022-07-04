#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_0.hxx"
#include "lolita/lolita_core_n_2.hxx"

int main(int argc, char** argv) {
    // Domain
    auto constexpr domain = lolita2::Domain(2, lolita2::Domain::Frame::Cartesian);

    // Mesh
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    auto msh2 = lolita2::geometry::MeshParser<lolita2::MeshData::Format::Gmsh, domain>(file_path);

    std::cout << msh2.mesh_data_ << std::endl;
}