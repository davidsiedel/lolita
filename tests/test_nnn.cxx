//
// Created by dsiedel on 28/05/22.
//

#include <span>
#include <ranges>
#include "gtest/gtest.h"

//#include "lolita/lolita.hxx"
//#include "lolita/lolita_utility.hxx"
//#include "lolita/lolita_algebra.hxx"
//#include "lolita/lolita_user.hxx"
//#include "lolita/lolita_core_1.hxx"
//#include "lolita/lolita_core_2.hxx"
//#include "lolita/lolita_core_4.hxx"
//#include "lolita/lolita_core_5.hxx"

#include "lolita/lolita_core_5_final.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_0.hxx"
#include "lolita/lolita_core_n_2.hxx"

template<auto arg>
struct Hello
{
    
};

TEST(test_lolita_nnn, test_lolita_nnn) {

    // auto constexpr d_domain = lolita2::Domain(2, lolita2::Domain::Frame::Cartesian);
    // auto constexpr cells_unknowns = lolita2::Unknown("Cell Field");
    // auto constexpr faces_unknowns = lolita2::Unknown("Face Field");
    // auto constexpr cell_basis = lolita2::FieldG(1, 2, lolita2::FieldG::Basis::Monomial, lolita2::FieldG::Frame::Cell);
    // auto constexpr face_basis = lolita2::FieldG(1, 1, lolita2::FieldG::Basis::Monomial, lolita2::FieldG::Frame::Face);
    // auto constexpr grad_basis = lolita2::FieldG(2, 1, lolita2::FieldG::Basis::Monomial, lolita2::FieldG::Frame::Cell);
    // auto constexpr ghost = lolita2::Field<>("Cell Field", cell_basis);
    // auto constexpr cell_field = lolita2::Field<cells_unknowns>("Cell Field", cell_basis);
    // auto constexpr face_field = lolita2::Field<faces_unknowns>("Face Field", face_basis);
    // auto constexpr grad_field = lolita2::Field<cell_field, face_field>("Gradient", grad_basis);
    // auto constexpr stab_field = lolita2::Field<grad_field>("Stabilization", face_basis);

    auto load_new = lolita2::Load([](lolita::domain::Point const &p, lolita::real const &t) { return 1.0; });

    // std::cout << lolita::core::finite_element::basis::FiniteElementBasisTraits<lolita::core::geometry::Element::LinearTriangle(), lolita::core::finite_element::basis::Basis::Monomial(), 4>::dim_ << std::endl;

    // lolita::core::finite_element::unknown::ABCD abcd{lolita::core::finite_element::unknown::ABCD::A, 2};

    auto constexpr domain = lolita::domain::Domain("Middle", 2, lolita::domain::Frame::Cartesian());
    auto constexpr u_unknown = lolita::field::Unknown("Displacement", 1, lolita::field::Mapping::LargeStrain());
    auto constexpr d_unknown = lolita::field::Unknown("Damage", 0, lolita::field::Mapping::Gradient(), lolita::field::Mapping::Identity());
    auto constexpr u_discretization = lolita::finite_element::HybridHighOrder(1, 1);
    auto constexpr d_discretization = lolita::finite_element::HybridHighOrder(1, 1);
    auto constexpr hho_u = lolita::finite_element::FiniteElement<u_unknown, 'a', u_discretization>(lolita::finite_element::Quadrature::Gauss(), 2);
    auto constexpr hho_d = lolita::finite_element::FiniteElement<d_unknown, 'a', d_discretization>(lolita::finite_element::Quadrature::Gauss(), 2);
    auto constexpr mixed_u = lolita::finite_element::ElementGroup<hho_u>();
    auto constexpr mixed_d = lolita::finite_element::ElementGroup<hho_d>();
    auto constexpr mixed_ud = lolita::finite_element::ElementGroup<hho_u, hho_d>();
    auto constexpr coupled_ud = lolita::finite_element::ElementGroup<mixed_u, mixed_d>();

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

    auto mgipara = lolita::behaviour::MgisParameter("Temperature", [](lolita::domain::Point const & p, lolita::real const & t) { return 1.0; });
    auto bhv_u = lolita::behaviour::MgisBehaviour("Displacement", "SQUARE", path, name, hyp, {mgipara});
    auto bhv_d = lolita::behaviour::MgisBehaviour("Damage", "SQUARE", path, name, hyp, {mgipara});

    auto load = lolita::finite_element::LoadEntry("Displacement", "SQUARE", 2, 0, 0, [](lolita::domain::Point const &p, lolita::real const &t) { return 1.0; },
                                                  lolita::finite_element::Load::Natural());
    auto load_top_x = lolita::finite_element::LoadEntry("Displacement", "TOP", 1, 0, 0, [](lolita::domain::Point const &p, lolita::real const &t) { return 1.0; },
                                                        lolita::finite_element::Load::Constraint());
    auto load_top_y = lolita::finite_element::LoadEntry("Displacement", "TOP", 1, 1, 0, [](lolita::domain::Point const &p, lolita::real const &t) { return 1.0; },
                                                        lolita::finite_element::Load::Constraint());
    auto file_path = "";
    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";

    auto msh2 = lolita::core::mesh::MeshParser<lolita::mesh::Format::Gmsh, domain, mixed_ud>(file_path, {load, load_top_x, load_top_y}, {bhv_u, bhv_d});

}
