#include "gtest/gtest.h"
#include "lolita_lolita/lolita_core/lolita_core_n_5000.hxx"

TEST(t0, t0)
{
    // constants
    auto constexpr domain = lolita::Domain::cartesian(2);
    auto constexpr cell_basis = lolita::Basis::monomial(1);
    auto constexpr face_basis = lolita::Basis::monomial(1);
    auto constexpr quadrature = lolita::Quadrature::gauss(2);
    auto constexpr cells = lolita::ElementType::cells(domain);
    auto constexpr faces = lolita::ElementType::faces(domain);
    // fields
    auto constexpr displacement_field = lolita::Field::vector();
    auto constexpr damage_field = lolita::Field::scalar();
    // generalized strains
    auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(displacement_field, lolita::Mapping::smallStrain());
    auto constexpr damage_generalized_strain = lolita::GeneralizedStrain(damage_field, lolita::Mapping::gradient(), lolita::Mapping::identity());
    // behaviors
    auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    auto constexpr damage_behavior = lolita::Behavior(damage_generalized_strain);
    // discretization
    auto constexpr hdg = lolita::HybridDiscontinuousGalerkin(cell_basis, face_basis, lolita::HybridDiscontinuousGalerkin::Stabilization::Hdg);
    // finite elements
    auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    auto constexpr damage_element =  lolita::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);
    // mesh    
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    // load
    auto load_f = std::make_shared<lolita::Load>([](lolita::Point const &p, lolita::Real const &t) { return 1.0; }, 0, 0);
    // mesh build
    auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    elements->addElement<cells>("Displacement", "SQUARE");
    elements->addElement<faces>("Displacement", "SQUARE");
    elements->addElement<cells>("Damage", "SQUARE");
    elements->addElement<faces>("Damage", "SQUARE");
    // dofs
    auto degree_of_freedom = elements->addDegreeOfFreedom<faces, displacement_field, face_basis>("Displacement", "SQUARE", "Displacement");
    // load
    auto load = elements->addLoad<faces>("Displacement", "SQUARE", [](lolita::Point const &p, lolita::Real const &t) { return 1.0; }, 0, 0);
    // bhv
    auto lib_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/bhv_micromorphic/src/libBehaviour.so";
    auto lib_name = "MicromorphicDamageII";
    auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
        mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
        mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    };
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;

    // auto bhvv = std::make_shared<mgis::behaviour::Behaviour>(mgis::behaviour::load(opts, lib_path, lib_name, hyp));
    auto bhvv = std::make_shared<mgis::behaviour::Behaviour>(mgis::behaviour::load(lib_path, lib_name, hyp));
    //
    auto behavior_data = std::make_shared<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* bhvv));
    mgis::behaviour::setMaterialProperty(behavior_data->s0, "FractureEnergy", 1.0);
    auto behaviour_view = mgis::behaviour::make_view(* behavior_data);
    auto result = mgis::behaviour::integrate(* std::make_unique<mgis::behaviour::BehaviourDataView>(mgis::behaviour::make_view(* behavior_data)), * bhvv);
    std::cout << "result : " << result << std::endl;
    for (auto const val : behavior_data->s1.gradients)
    {
        std::cout << val << std::endl;
    }
    //
    elements->addBehavior<cells, quadrature>("SQUARE", bhvv);
    elements->addBehavior<cells>("Displacement", "SQUARE", lib_name);
    // std::cout << lolita::numerics::sqrt_2 << std::endl;
    // elements->makeQuadrature<cells, displacement_field, quadrature, hdg>("SQUARE", "Displacement");
    //
    // problem build
    // elements->activate<cells>("Displacement", "SQUARE");
    // elements->activate<faces>("Displacement", "SQUARE");
    // elements->activate<cells>("Damage", "SQUARE");
    // elements->activate<faces>("Damage", "SQUARE");
    // elements->activate<displacement_element, lolita::core::ElementType::faces(domain)>("SQUARE");
    // elements->activate<damage_element, lolita::core::ElementType::cells(domain)>("SQUARE");
    // elements->activate<damage_element, lolita::core::ElementType::faces(domain)>("SQUARE");
    // elements->setLoad<displacement_element, lolita::core::ElementType::cells(domain)>(
    //     "SQUARE",
    //     0, 0, [](lolita::Point const &p, Real const &t) { return 1.0; }
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
    // std::cout << * elements << std::endl;

    
    
}