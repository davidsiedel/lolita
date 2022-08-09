#include "gtest/gtest.h"
#include "lolita/lolita_core_n_005.hxx"

TEST(t0, t0)
{

    auto constexpr jjj = std::tuple<char>{'A'};
    // constants
    auto constexpr domain = lolita2::Domain::cartesian(2);
    auto constexpr cell_basis = lolita2::Basis::monomial(1);
    auto constexpr face_basis = lolita2::Basis::monomial(1);
    auto constexpr quadrature = lolita2::Quadrature::gauss(2);
    auto constexpr cells = lolita2::geometry::ElementType::cells(domain);
    auto constexpr faces = lolita2::geometry::ElementType::faces(domain);
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
    auto constexpr hdg = lolita2::HybridDiscontinuousGalerkin(cell_basis, face_basis, lolita2::HybridDiscontinuousGalerkin::Stabilization::Hdg);
    // finite elements
    auto constexpr displacement_element =  lolita2::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    auto constexpr damage_element =  lolita2::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);
    // mesh    
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    // dofs
    auto degree_of_freedom = std::make_shared<lolita2::geometry::DegreeOfFreedom>("FaceDisplacement");
    auto load_f = std::make_shared<lolita2::geometry::Load>([](lolita2::Point const &p, lolita::real const &t) { return 1.0; }, 0, 0);
    // mesh build
    auto elements = lolita2::geometry::MeshFileParser(file_path).template makeFiniteElementSet<domain, displacement_element, damage_element>();
    // problem build
    elements->make<displacement_element, cells, hdg, quadrature>("SQUARE");
    elements->activate<displacement_element, faces>("SQUARE");
    elements->activate<damage_element, cells>("SQUARE");
    elements->activate<damage_element, faces>("SQUARE");
    elements->setDegreeOfFreedom<displacement_element, faces, displacement_field, face_basis>("SQUARE", degree_of_freedom);
    elements->setLoad<displacement_element, cells>("SQUARE", 0, 0, [](lolita2::Point const &p, lolita::real const &t) { return 1.0; });
    elements->setBehavior<displacement_behavior, cells, quadrature>("SQUARE");
    // show mesh
    // std::cout << * elements << std::endl;
    std::cout << degree_of_freedom->coefficients_.size() << std::endl;
    degree_of_freedom->coefficients_.setZero();
    //
    for (auto const & element : elements->getElements<1, 0>())
    {
        element.second->getFiniteElement<0>()->isActivated();
        std::cout
        <<
        element.second->getFiniteElement<0>()->getDegreeOfFreedom("FaceDisplacement").getCoefficients<displacement_field, face_basis>() << std::endl;
    }

    // auto constexpr hhjk = lolita2::expand<lolita::utility::Holder, displacement_behavior>();
    auto constexpr kkkm = lolita::utility::Aggregate<int, char, double>(1, 'A', 2.0);
    auto constexpr kkkm2 = lolita::utility::aggregate_template_t<lolita::utility::Holder, kkkm>();
    // using RESS = std::tuple<lolita::utility::Holder<1>, lolita::utility::Holder<'A'>, lolita::utility::Holder<1>>;
    lolita::utility::TD<decltype(kkkm2)>();
    using SLICE = lolita::utility::tuple_slice_t<std::tuple<int, char, double, float, bool, long>, 2, 4>;
    lolita::utility::TD<SLICE>();
    using UNIQUE = lolita::utility::tuple_unique_t<std::tuple<int, char, int, float, bool, float>>;
    using CAT = lolita::utility::tuple_cat_t<std::tuple<>, std::tuple<int, float>>;
    lolita::utility::TD<UNIQUE>();
    
    
}