//
// Created by dsiedel on 28/05/22.
//

#include <span>
#include <ranges>
#include "gtest/gtest.h"

#include "lolita/lolita_user.hxx"
#include "lolita/new_element.hxx"
#include "lolita/new_finite_element.hxx"
#include "lolita/new_mesh.hxx"
//#include "lolita/lolita_core_mesh2.hxx"

TEST(test_lolita_22, test_lolita_22) {

//    std::cout << 1 << std::endl;


//    std::cout << "dynamic matrix size : " << sizeof(lolita::matrix::Vector<lolita::real>) << std::endl;
//    std::cout << "empty matrix size : " << sizeof(lolita::matrix::Vector<lolita::real, 0>) << std::endl;
//    std::cout << "2x2 array of empty matrix size : " << sizeof(std::array<std::array<lolita::matrix::Vector<lolita::real, 0>, 2>, 2>) << std::endl;

    auto constexpr domain = lolita::geometry::Domain("Middle", 2, lolita::geometry::Frame::Cartesian());
    auto constexpr u_unknown = lolita::field::Unknown("Displacement", 1, lolita::field::Mapping::LargeStrain());
    auto constexpr p_unknown = lolita::field::Unknown("Pressure", 1, lolita::field::Mapping::Identity());
    auto constexpr d_unknown = lolita::field::Unknown("Damage", 0, lolita::field::Mapping::Gradient(), lolita::field::Mapping::Identity());
    auto constexpr u_discretization = lolita::finite_element::HybridHighOrder(1, 1);
    auto constexpr d_discretization = lolita::finite_element::HybridHighOrder(1, 1);
    auto constexpr u_voce = lolita::behaviour::MgisBehaviour3<u_unknown, p_unknown>();
    auto constexpr u_elasticity = lolita::behaviour::MgisBehaviour3<u_unknown>();
    auto constexpr d_phase_field = lolita::behaviour::MgisBehaviour3<d_unknown>();
    auto constexpr hho_u0 = lolita::finite_element::FiniteElement<u_unknown, u_voce, u_discretization>(lolita::finite_element::Quadrature::Gauss(), 2);
    auto constexpr hho_p = lolita::finite_element::FiniteElement<p_unknown, u_voce, u_discretization>(lolita::finite_element::Quadrature::Gauss(), 2);
    auto constexpr hho_u1 = lolita::finite_element::FiniteElement<u_unknown, u_elasticity, u_discretization>(lolita::finite_element::Quadrature::Gauss(), 2);
    auto constexpr hho_d = lolita::finite_element::FiniteElement<d_unknown, d_phase_field, d_discretization>(lolita::finite_element::Quadrature::Gauss(), 2);
    auto constexpr mixed_0 = lolita::finite_element::ElementGroup<hho_u0>();
    auto constexpr mixed_1 = lolita::finite_element::ElementGroup<hho_d>();
    auto constexpr mixed_2 = lolita::finite_element::ElementGroup<hho_u0, hho_d>();
    auto constexpr coupled_0 = lolita::finite_element::ElementGroup<lolita::finite_element::ElementGroup<hho_u0>{}, lolita::finite_element::ElementGroup<hho_d>{}>();


    auto constexpr conttest = lolita::finite_element::ElementGroup<
            lolita::finite_element::ElementGroup<hho_u0, hho_u0>{},
            lolita::finite_element::ElementGroup<hho_d>{},
            lolita::finite_element::ElementGroup<lolita::finite_element::ElementGroup<hho_u0, hho_u0>{}>{}
    >();

    std::cout << "count : " << conttest.count() << std::endl;
    std::cout << "index : " << conttest.getFiniteElementIndex<hho_d>() << std::endl;

    Eigen::Matrix<std::shared_ptr<lolita::index>, 2, 2>();

//    std::cout << lolita::utility::EnumA("Hello") << std::endl;
//
//    for (auto const & c : lolita::utility::EnumA("Hello").tag_) {
//        std::cout << c << std::endl;
//    }

    static_assert(lolita::geometry::Domain("Middle", 2, lolita::geometry::Frame::Cartesian()) == domain);

    static_assert(lolita::finite_element::HybridHighOrderFiniteElementConcept<hho_u0>);

    auto path = "";
    path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/behaviour/src/libBehaviour.so";
    path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/bhv_micromorphic/src/libBehaviour.so";
    auto name = "";
    name = "Voce";
    name = "MicromorphicDamageII";
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    auto opt = mgis::behaviour::FiniteStrainBehaviourOptions{
            mgis::behaviour::FiniteStrainBehaviourOptions::StressMeasure::PK1,
            mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator::DPK1_DF
    };

    auto mgipara = lolita::behaviour::MgisParameter("Temperature", [](lolita::geometry::Point const & p, lolita::real const & t) { return 1.0; });
    auto bhv_u = lolita::behaviour::MgisBehaviour("Displacement", "SQUARE", path, name, hyp, {mgipara});
    auto bhv_d = lolita::behaviour::MgisBehaviour("Damage", "SQUARE", path, name, hyp, {mgipara});

    auto load = lolita::finite_element::Load("Displacement", "SQUARE", 2, 0, 0, [](lolita::geometry::Point const &p, lolita::real const &t) { return 1.0; },
                                             lolita::finite_element::Loading::Natural());
    auto load_top_x = lolita::finite_element::Load("Displacement", "TOP", 1, 0, 0, [](lolita::geometry::Point const &p, lolita::real const &t) { return 1.0; },
                                                   lolita::finite_element::Loading::Constraint());
    auto load_top_y = lolita::finite_element::Load("Displacement", "TOP", 1, 1, 0, [](lolita::geometry::Point const &p, lolita::real const &t) { return 1.0; },
                                                   lolita::finite_element::Loading::Constraint());

    auto bhvvv = lolita::behaviour::Behaviour<u_voce>();

//    auto bhvholder = lolita::behaviour::BehaviourHolder(bhvvv);

    auto file_path = "";
    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";
//    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/perforated_strip_huge.msh";

    auto testttt = lolita::core::finite_element::FEObject<lolita::core::Element::LinearTriangle(), domain, mixed_2>();

//    std::cout << "end_ 0 : " << testttt.end_ << std::endl;
//    std::cout << "end_ 1 : " << std::get<0>(testttt.elements_)->end_ << std::endl;
//    testttt.getElement<0>()->

//    auto msh2 = lolita::core::mesh::MeshParser<lolita::mesh::Format::Gmsh, domain, hho_u0, hho_d>(file_path, {load, load_top_x, load_top_y}, {bhv_u, bhv_d});
    auto msh2 = lolita::core::mesh::MeshParser<lolita::mesh::Format::Gmsh, domain, mixed_2>(file_path, {load, load_top_x, load_top_y}, {bhv_u, bhv_d});
    auto msh1 = lolita::core::mesh::MeshParser<lolita::mesh::Format::Gmsh, domain, coupled_0>(file_path, {load, load_top_x, load_top_y}, {bhv_u, bhv_d});

    std::cout << msh2.mesh_data_ << std::endl;

    std::cout << lolita::finite_element::Quadrature::Gauss() << std::endl;


    std::cout << lolita::core::ElementDescription<lolita::core::Element::LinearSegment(), domain>::getCoordinates().dim_ << std::endl;

    std::cout << lolita::core::FiniteElementDescription<lolita::core::Element::LinearSegment(), domain, hho_u0>::getNumUnknowns<lolita::core::finite_element::unknown::UnknownType::Structural()>() << std::endl;

    std::cout << lolita::core::FiniteElementDescription<lolita::core::Element::LinearTriangle(), domain, hho_u0>::getNumUnknowns<lolita::core::finite_element::unknown::UnknownType::Subsidiary()>() << std::endl;

    std::cout << lolita::core::FiniteElementDescription<lolita::core::Element::LinearTriangle(), domain, hho_u0>::getNumUnknowns<lolita::core::finite_element::unknown::UnknownType::Structural()>() << std::endl;

    std::cout << lolita::core::FiniteElementDescription<lolita::core::Element::LinearTriangle(), domain, hho_u0>::getNumUnknowns() << std::endl;

    std::cout << lolita::core::FiniteElementDescription<lolita::core::Element::LinearTriangle(), domain, hho_u0>::getNumOwnedUnknowns() << std::endl;

}
