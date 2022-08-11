#include "gtest/gtest.h"
#include "lolita/lolita_core_n_005.hxx"

TEST(t0, t0)
{
    // constants
    auto constexpr domain = lolita2::Domain::cartesian(2);
    auto constexpr cell_basis = lolita2::Basis::monomial(1);
    auto constexpr face_basis = lolita2::Basis::monomial(1);
    auto constexpr quadrature = lolita2::Quadrature::gauss(2);
    auto constexpr cells = lolita2::geometry::ElementType::cells(domain);
    auto constexpr faces = lolita2::geometry::ElementType::faces(domain);
    // fields
    auto constexpr displacement_field = lolita2::Field::vector();
    auto constexpr damage_field = lolita2::Field::scalar();
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
    auto elements = lolita2::geometry::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    elements->addElement<cells>("Displacement", "SQUARE");
    elements->addElement<faces>("Displacement", "SQUARE");
    elements->addElement<cells>("Damage", "SQUARE");
    elements->addElement<faces>("Damage", "SQUARE");
    // dofs
    elements->addDegreeOfFreedom<faces, displacement_field, face_basis>("Displacement", "SQUARE", degree_of_freedom);
    // load
    elements->addLoad<faces>("Displacement", "SQUARE", load_f);
    // bhv
    auto lib_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/behaviour/src/libBehaviour.so";
    auto lib_name = "Voce";
    auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
        mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
        mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    };
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;

    auto bhvv = std::make_shared<mgis::behaviour::Behaviour>(mgis::behaviour::load(opts, lib_path, lib_name, hyp));
    elements->addBehavior<cells, quadrature>("SQUARE", bhvv);
    elements->addBehavior<cells>("Displacement", "SQUARE", lib_name);
    // elements->makeQuadrature<cells, displacement_field, quadrature, hdg>("SQUARE", "Displacement");
    //
    // problem build
    // elements->activate<cells>("Displacement", "SQUARE");
    // elements->activate<faces>("Displacement", "SQUARE");
    // elements->activate<cells>("Damage", "SQUARE");
    // elements->activate<faces>("Damage", "SQUARE");
    // elements->activate<displacement_element, lolita2::geometry::ElementType::faces(domain)>("SQUARE");
    // elements->activate<damage_element, lolita2::geometry::ElementType::cells(domain)>("SQUARE");
    // elements->activate<damage_element, lolita2::geometry::ElementType::faces(domain)>("SQUARE");
    // elements->setLoad<displacement_element, lolita2::geometry::ElementType::cells(domain)>(
    //     "SQUARE",
    //     0, 0, [](lolita2::Point const &p, lolita::real const &t) { return 1.0; }
    // );
    // // show mesh
    // std::cout << degree_of_freedom->coefficients_.size() << std::endl;
    // degree_of_freedom->coefficients_.setZero();
    // //
    // for (auto const & element : elements->getElements<1, 0>())
    // {
    //     element.second->getFiniteElement<0>()->isActivated();
    //     std::cout
    //     <<
    //     element.second->getFiniteElement<0>()->getDegreeOfFreedom("FaceDisplacement")->getCoefficients<displacement_field, face_basis>() << std::endl;
    // }
    //
    std::cout << * elements << std::endl;

    
    
}