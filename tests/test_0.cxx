#include "gtest/gtest.h"
#include "lolita_lolita/lolita_core/lolita_core_n_5000.hxx"

TEST(t0, t0)
{
    // declaring behavior
    auto lib_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/bhv_micromorphic/src/libBehaviour.so";
    auto lib_name = "MicromorphicDamageII";
    auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
        mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
        mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    };
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    //
    //
    //
    // constants
    auto constexpr domain = lolita::Domain::cartesian(2);
    auto constexpr quadrature = lolita::Quadrature::gauss(2);
    auto constexpr cells = lolita::ElementType::cells(domain);
    auto constexpr faces = lolita::ElementType::faces(domain);
    // fields
    auto constexpr displacement_field = lolita::Field::vector();
    auto constexpr damage_field = lolita::Field::scalar();
    // generalized strains
    auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(lolita::Field::vector(), lolita::Mapping::smallStrain());
    auto constexpr damage_generalized_strain = lolita::GeneralizedStrain(lolita::Field::scalar(), lolita::Mapping::gradient(), lolita::Mapping::identity());
    // behaviors
    auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    auto constexpr damage_behavior = lolita::Behavior(damage_generalized_strain);
    // discretization
    auto constexpr hdg = lolita::HybridDiscontinuousGalerkin::hybridDiscontinuousGalerkin(1, 1);
    // finite elements
    auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    auto constexpr damage_element =  lolita::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);
    // mesh    
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    // mesh build
    auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    // dofs
    auto face_displacement = elements->setDegreeOfFreedom<faces, displacement_field, hdg.getFaceBasis()>("SQUARE", "Displacement");
    auto cell_displacement = elements->setDegreeOfFreedom<cells, displacement_field, hdg.getCellBasis()>("SQUARE", "Displacement");
    auto face_damage = elements->setDegreeOfFreedom<faces, damage_field, hdg.getFaceBasis()>("SQUARE", "Damage");
    auto cell_damage = elements->setDegreeOfFreedom<cells, damage_field, hdg.getCellBasis()>("SQUARE", "Damage");
    face_displacement->coefficients_.setZero();
    cell_displacement->coefficients_.setZero();
    face_damage->coefficients_.setZero();
    cell_damage->coefficients_.setZero();
    // systems
    auto displacement_system = lolita::System(face_displacement->coefficients_.size());
    auto damage_system = lolita::System(face_damage->coefficients_.size());
    // load
    auto load = elements->setConstraint<faces>("TOP", "Pull", [](lolita::Point const & p, lolita::Real const & t) { return t; }, 0, 0);
    // adding behavior
    auto micromorphic_damage = elements->setBehavior<cells, quadrature>("SQUARE", lib_path, lib_name, hyp);
    //making operators
    elements->setStrainOperators<cells, damage_element, hdg>("SQUARE", "MicromorphicDamageII", "Damage");
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamageII", "FractureEnergy", [](lolita::Point const & p) { return 1.0; });
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamageII", "CharacteristicLength", [](lolita::Point const & p) { return 1.0; });
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamageII", "PenalisationFactor", [](lolita::Point const & p) { return 1.0; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDamageII", "Temperature", [](lolita::Point const & p) { return 1.0; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDamageII", "EnergyReleaseRate", [](lolita::Point const & p) { return 1.0; });
    elements->setElementOperators<cells, damage_element, hdg>("SQUARE", "Stabilization");
    elements->setStrainValues<cells, damage_element, hdg>("SQUARE", "MicromorphicDamageII", "Damage");
    elements->integrate<cells>("SQUARE", "MicromorphicDamageII");
    elements->assemble<cells, damage_element, hdg>("SQUARE", "MicromorphicDamageII", "Damage", damage_system);
    damage_system.getCorrection();



    auto my_cells = elements->getSubSet<cells>("SQUARE");
    auto lll = my_cells->template setDegreeOfFreedom<displacement_field, hdg.getFaceBasis()>("Displacement");
    auto llm = my_cells->setLoad("Pull", [](lolita::Point const & p, lolita::Real const & t) { return t; }, 0, 0);

    // lolita::setConstraint

    // auto matt = lolita::RealMatrix<>();
    // auto mat_map = std::map<std::string, lolita::RealMatrix<>>();
    // auto matt_t = lolita::RealMatrix<2, 2>();
    // mat_map["Hey"] = matt_t;





    // //
    // // auto behavior_data = std::make_shared<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* bhvv));
    // auto behavior_data = mgis::behaviour::BehaviourData(* micromorphic_damage);
    // // mgis::behaviour::setMaterialProperty(behavior_data->s0, "FractureEnergy", 1.0);
    // mgis::behaviour::setMaterialProperty(behavior_data.s0, "FractureEnergy", 1.0);
    // // auto behaviour_view = mgis::behaviour::make_view(* behavior_data);
    // auto behaviour_view = mgis::behaviour::make_view(behavior_data);
    // // auto result = mgis::behaviour::integrate(* std::make_unique<mgis::behaviour::BehaviourDataView>(mgis::behaviour::make_view(behavior_data)), * bhvv);
    // auto result = mgis::behaviour::integrate(* std::make_unique<mgis::behaviour::BehaviourDataView>(mgis::behaviour::make_view(behavior_data)), * micromorphic_damage);
    // std::cout << "result : " << result << std::endl;
    // for (auto const val : behavior_data.s1.gradients)
    // {
    //     std::cout << val << std::endl;
    // }

    // auto opss = lolita::utility::Holderr<int>();
    // opss.setItem("1", 1);
    // opss.setItem("2", 2);
    // opss.setItem("3", 3);
    // opss.getValue("3") = 4;
    // opss.setItem("3", 5);
    // std::cout << opss.getValue("3") << std::endl;


    // auto map = std::unordered_map<std::basic_string_view<lolita::Character>, int>();
    // auto map_set_tick = std::chrono::high_resolution_clock::now();
    // map["1"] = 1;
    // // map["2"] = 2;
    // // map["3"] = 3;
    // auto map_set_tock = std::chrono::high_resolution_clock::now();
    // auto map_set_time = std::chrono::duration<double>(map_set_tock - map_set_tick);
    // std::cout << "map_set_time : " << map_set_time.count() << std::endl;

    // auto hld = lolita::utility::Holderr<int>();
    // auto hld_set_tick = std::chrono::high_resolution_clock::now();
    // hld.setItem("1", 1);
    // // hld.setItem("2", 2);
    // // hld.setItem("3", 3);
    // auto hld_set_tock = std::chrono::high_resolution_clock::now();
    // auto hld_set_time = std::chrono::duration<double>(hld_set_tock - hld_set_tick);
    // std::cout << "hld_set_time : " << hld_set_time.count() << std::endl;
    
    // auto map_get_tick = std::chrono::high_resolution_clock::now();
    // map.at("1");
    // auto map_get_tock = std::chrono::high_resolution_clock::now();
    // auto map_get_time = std::chrono::duration<double>(map_get_tock - map_get_tick);
    // std::cout << "map_get_time : " << map_get_time.count() << std::endl;
    
    // auto vec = std::vector<int>();
    // vec.push_back(1);
    // auto hld_get_tick = std::chrono::high_resolution_clock::now();
    // // hld.getValue("1");
    // auto val = vec[0];
    // auto hld_get_tock = std::chrono::high_resolution_clock::now();
    // auto hld_get_time = std::chrono::duration<double>(hld_get_tock - hld_get_tick);
    // std::cout << "hld_get_time : " << hld_get_time.count() << std::endl;

    // std::cout << "map size : " << sizeof(hld) << std::endl;
    // std::cout << "hld size : " << sizeof(map) << std::endl;


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