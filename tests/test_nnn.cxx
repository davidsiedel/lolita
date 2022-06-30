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

TEST(test_lolita_nnn, test_lolita_nnn) {

    // std::cout << lolita::core2::finite_element::basis::FiniteElementBasisTraits<lolita::core2::geometry::Element::LinearTriangle(), lolita::core2::finite_element::basis::Basis::Monomial(), 4>::dim_ << std::endl;

    // lolita::core2::finite_element::unknown::ABCD abcd{lolita::core2::finite_element::unknown::ABCD::A, 2};

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

    //    auto bhvholder = lolita::behaviour::BehaviourHolder(bhvvv);


    /*
     * 
     */
    auto file_path = "";
    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";

    auto constexpr node = lolita::core2::geometry::Element::Node();
    auto constexpr segment = lolita::core2::geometry::Element::LinearSegment();
    auto constexpr triangle = lolita::core2::geometry::Element::LinearQuadrangle();
    //    lolita::core2::geometry::ElementGeometryTraits<node, domain>::getNeighbours<arg0, arg1>();
    std::cout << "node with dim " << node.dim() << std::endl;
    std::cout << lolita::core2::geometry::DomainTraits<domain>::getElementCoordinates<node>() << std::endl;
    std::cout << "segment with dim " << segment.dim() << std::endl;
    std::cout << lolita::core2::geometry::DomainTraits<domain>::getElementCoordinates<segment>() << std::endl;
    std::cout << "triangle with dim " << triangle.dim() << std::endl;
    std::cout << lolita::core2::geometry::DomainTraits<domain>::getElementCoordinates<triangle>() << std::endl;

    std::cout << "elements size : " << std::tuple_size_v<lolita::core2::geometry::Elements<lolita::core2::geometry::detail::ElementView, domain>> << std::endl;

    auto msh2 = lolita::core2::mesh::MeshParser<lolita::mesh::Format::Gmsh, domain, mixed_ud>(file_path, {load, load_top_x, load_top_y}, {bhv_u, bhv_d});
    //    auto msh1 = lolita::core2::mesh::MeshParser<lolita::mesh::Format::Gmsh, domain, coupled_ud>(file_path, {load, load_top_x, load_top_y}, {bhv_u, bhv_d});

//    std::cout << "num unknowns displacement : " << msh2.mesh_data_.systems_[0].num_unknowns_ << std::endl;
//    std::cout << "num bindings displacement : " << msh2.mesh_data_.systems_[0].num_bindings_ << std::endl;
//    std::cout << "num unknowns damage : " << msh2.mesh_data_.systems_[1].num_unknowns_ << std::endl;
//    std::cout << "num bindings damage : " << msh2.mesh_data_.systems_[1].num_bindings_ << std::endl;
//
//    std::cout << "here : " << lolita::core2::finite_element::FiniteElementTraits<lolita::core2::geometry::Element::LinearTriangle(), domain, hho_u>::template getNumUnknowns<lolita::core2::finite_element::unknown::Unknown::Subsidiary()>() << std::endl;
//    std::cout << "here : " << lolita::core2::finite_element::FiniteElementTraits<lolita::core2::geometry::Element::LinearTriangle(), domain, hho_u>::getNumUnknowns() << std::endl;
//    std::cout << "here cell structural : " << lolita::core2::finite_element::FiniteElementTraits<lolita::core2::geometry::Element::LinearTriangle(), domain, hho_u>::getDimUnknowns<lolita::core2::finite_element::unknown::Unknown::Structural()>() << std::endl;
//    std::cout << "here face structural : " << lolita::core2::finite_element::FiniteElementTraits<lolita::core2::geometry::Element::LinearSegment(), domain, hho_u>::getDimUnknowns<lolita::core2::finite_element::unknown::Unknown::Structural()>() << std::endl;
//    std::cout << "here cell subsidiary : " << lolita::core2::finite_element::FiniteElementTraits<lolita::core2::geometry::Element::LinearTriangle(), domain, hho_u>::getDimUnknowns<lolita::core2::finite_element::unknown::Unknown::Subsidiary()>() << std::endl;
//    std::cout << "here face subsidiary : " << lolita::core2::finite_element::FiniteElementTraits<lolita::core2::geometry::Element::LinearSegment(), domain, hho_u>::getDimUnknowns<lolita::core2::finite_element::unknown::Unknown::Subsidiary()>() << std::endl;

//    using HHH = lolita::matrix::Vector<lolita::integer, 0>;
//    auto hhh_vec = HHH::LinSpaced(0, 0 + 0 - 1);
//    std::cout << "hhh_vec : " << std::endl;
//    std::cout << hhh_vec << std::endl;
//    using HHH2 = lolita::matrix::Vector<lolita::integer>;
//    auto hhh_vec2 = HHH2::LinSpaced(2, 0, 0 + 2 - 1);
//    std::cout << "hhh_vec2 : " << std::endl;
//    std::cout << hhh_vec2 << std::endl;
//
//    static_assert(lolita::core2::finite_element::unknown::Unknown("Coucou", false) == lolita::core2::finite_element::unknown::Unknown("Coucou", true));
//
   std::cout << msh2.mesh_data_ << std::endl;

}
