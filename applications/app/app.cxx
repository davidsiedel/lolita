#include "lolita/lolita_core_n_3.hxx"

int main(int argc, char** argv)
{

    // Domain
    auto constexpr domain = lolita2::Domain(2, lolita2::Domain::Cartesian);

    // Mesh
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    // auto msh = lolita2::geometry::MeshParser<lolita2::MeshFileFormat::Gmsh, domain>::makeMesh(file_path);

    // // for (auto const & element : msh.sets_[* msh.domains_[5]].getElements<domain.dim_, 0>())
    // // {
    // //     std::cout << element.second->tag_ << std::endl;
    // // }

    // auto constexpr field = lolita2::Field(2);

    // auto constexpr basis = lolita2::Basis(lolita2::Basis::Monomial, 2);
    
    // std::cout << msh << std::endl;

    // auto faces_unknown = std::make_shared<lolita2::geometry::Unknown>(lolita2::geometry::Unknown());
    
    // auto faces_unknowns = msh.template make<lolita2::geometry::FiniteElementUnknown, field, basis>(faces_unknown);

    // auto value = lolita::integer(5);

    // auto HDGf = msh.template make<1, lolita2::geometry::HybridDiscontinuousGalerkin, field, basis, basis>("SQUARE", faces_unknown);
    // auto HDGc = msh.template make<2, lolita2::geometry::HybridDiscontinuousGalerkin, field, basis, basis>("SQUARE", 1);

    // std::cout << msh.template getNumElements<1, 0>("SQUARE") << std::endl;
    // std::cout << faces_unknown->coefficients_.size() << std::endl;

}