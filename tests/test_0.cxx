#include "gtest/gtest.h"
#include "lolita/lolita_core_n_005.hxx"

TEST(t0, t0)
{
    // constants
    auto constexpr domain = lolita2::Domain::cartesian(2);
    auto constexpr basis = lolita2::Basis::monomial(1);
    auto constexpr quadrature = lolita2::Quadrature(lolita2::Quadrature::Gauss, 2);
    // fields
    auto constexpr displacement_field = lolita2::Field::vector("Displacement");
    auto constexpr damage_field = lolita2::Field::scalar("Damage");
    // generalized strains
    auto constexpr displacement_generalized_strain = lolita2::GeneralizedStrain(displacement_field, lolita2::Mapping::smallStrain());
    auto constexpr damage_generalized_strain = lolita2::GeneralizedStrain(damage_field, lolita2::Mapping::gradient(), lolita2::Mapping::identity());
    // behaviors
    auto constexpr displacement_behavior = lolita2::Behavior(displacement_generalized_strain);
    auto constexpr damage_behavior = lolita2::Behavior(damage_generalized_strain);
    // discretization
    auto constexpr hdg = lolita2::HybridDiscontinuousGalerkin(basis, basis, lolita2::HybridDiscontinuousGalerkin::Stabilization::Hdg);
    // finite elements
    auto constexpr displacement_element =  lolita2::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    auto constexpr damage_element =  lolita2::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);
    // mesh    
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    // dofs
    auto degree_of_freedom = std::make_shared<lolita2::geometry::DegreeOfFreedom>("FaceDisplacement");
    // mesh build
    auto elements = lolita2::geometry::MeshFileParser(file_path).template makeFiniteElementSet<domain, displacement_element, damage_element>();
    // problem build
    elements->activate<domain.dim_, 0, 1>("SQUARE");
    elements->activate<domain.dim_ - 1, 0, 1>("SQUARE");
    elements->setDegreeOfFreedom<displacement_field, basis, 1, 0>("SQUARE", degree_of_freedom);
    // show mesh
    std::cout << * elements << std::endl;
    std::cout << degree_of_freedom->coefficients_.size() << std::endl;
    
}