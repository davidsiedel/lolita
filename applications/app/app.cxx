#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_2.hxx"

int main(int argc, char** argv) {

    // Domain
    auto constexpr domain = lolita2::Domain(2, lolita2::Domain::Frame::Cartesian);

    // Mesh
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    auto msh = lolita2::geometry::MeshParser<lolita2::MeshData::Format::Gmsh, domain>::makeMesh(file_path);

    for (auto const & element : msh.sets_[* msh.domains_[5]].getElements<domain.dim_, 0>())
    {
        std::cout << element.second->tag_ << std::endl;
    }
    
    std::cout << msh << std::endl;

}