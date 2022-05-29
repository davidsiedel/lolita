//
// Created by dsiedel on 03/05/22.
//

#include <span>
#include <ranges>
#include "gtest/gtest.h"

#include "lolita/lolita_user.hxx"
#include "lolita/lolita_element.hxx"
#include "lolita/lolita_core_finite_element.hxx"
//#include "lolita/lolita_core_mesh2.hxx"

TEST(test_lolita_2, test_lolita_2)
{

//    using namespace lolita;
//    using namespace lolita::core;

    std::cout << "dynamic matrix size : " << sizeof(lolita::matrix::Vector<lolita::real>) << std::endl;
    std::cout << "empty matrix size : " << sizeof(lolita::matrix::Vector<lolita::real, 0>) << std::endl;
    std::cout << "2x2 array of empty matrix size : " << sizeof(std::array<std::array<lolita::matrix::Vector<lolita::real, 0>, 2>, 2>) << std::endl;
//    std::cout << "empty unknowns : " << sizeof(lolita::core::element::Unknowns<0, 0, 0>) << std::endl;

    auto constexpr domain = lolita::geometry::Domain("Middle", 2, lolita::geometry::Frame::Cartesian);
//    auto constexpr displacement = lolita::field::DegreeOfFreedom("Displacement", 4, lolita::field::Mapping::Gradient);
//    auto constexpr displacement_field = lolita::field::Tensor("Displacement", 1);
    auto constexpr u_unknown = lolita::field::Unknown("Displacement", 1, lolita::field::Mapping::LargeStrainPlane);
    auto constexpr p_unknown = lolita::field::Unknown("Pressure", 1, lolita::field::Mapping::Identity);
    auto constexpr d_unknown = lolita::field::Unknown("Damage", 0, lolita::field::Mapping::Gradient, lolita::field::Mapping::Identity);
    auto constexpr u_discretization = lolita::finite_element::HybridHighOrder(1, 1);
    auto constexpr d_discretization = lolita::finite_element::HybridHighOrder(1, 1);
    auto constexpr u_voce = lolita::behaviour::MgisBehaviour3<u_unknown, p_unknown>();
    auto constexpr u_elasticity = lolita::behaviour::MgisBehaviour3<u_unknown>();
    auto constexpr d_phase_field = lolita::behaviour::MgisBehaviour3<d_unknown>();
    auto constexpr hho_u0 = lolita::finite_element::FiniteElement<u_unknown, u_voce, u_discretization>(lolita::finite_element::Quadrature::Gauss, 2);
    auto constexpr hho_p = lolita::finite_element::FiniteElement<p_unknown, u_voce, u_discretization>(lolita::finite_element::Quadrature::Gauss, 2);
    auto constexpr hho_u1 = lolita::finite_element::FiniteElement<u_unknown, u_elasticity, u_discretization>(lolita::finite_element::Quadrature::Gauss, 2);
    auto constexpr hho_d = lolita::finite_element::FiniteElement<d_unknown, d_phase_field, d_discretization>(lolita::finite_element::Quadrature::Gauss, 2);

    static_assert(lolita::finite_element::HybridHighOrderFiniteElementConcept<hho_u0>);

//    std::cout << "num struct ukns : " << lolita::core::element::FEObjectUnknownsPol<lolita::core::geometry::tri_03, domain, hho_u0>::_dim_unknowns<lolita::core::unknown::UnknownType::Structural>() << std::endl;
//    std::cout << "num subsid ukns : " << lolita::core::element::FEObjectUnknownsPol<lolita::core::geometry::tri_03, domain, hho_u0>::_dim_unknowns<lolita::core::unknown::UnknownType::Subsidiary>() << std::endl;
//    std::cout << "num struct ukns : " << lolita::core::element::FEObjectUnknownsPol<lolita::core::geometry::seg_02, domain, hho_u0>::_dim_unknowns<lolita::core::unknown::UnknownType::Structural>() << std::endl;
//    std::cout << "num subsid ukns : " << lolita::core::element::FEObjectUnknownsPol<lolita::core::geometry::seg_02, domain, hho_u0>::_dim_unknowns<lolita::core::unknown::UnknownType::Subsidiary>() << std::endl;
//
//
//    std::cout << "num subsid ukns : " << lolita::core::unknown::numUnknowns<lolita::core::element::FiniteElementTriplet(lolita::core::geometry::tri_03, domain, hho_u0), lolita::core::unknown::UnknownType::Structural>() << std::endl;

    std::tuple_size_v<std::tuple<>>;

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

    lolita::matrix::Vector<lolita::integer, 3> aaa;
    aaa = lolita::matrix::Vector<lolita::integer, 3>::LinSpaced(0, 3);
    std::cout << "lin spaces d :" << aaa << std::endl;
    lolita::matrix::Vector<lolita::integer> bbb;
    bbb = lolita::matrix::Vector<lolita::integer>::LinSpaced(3, 0, 3);
    std::cout << "lin spaces d :" << bbb << std::endl;

//    opt.stress_measure = mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
//    opt.tangent_operator = mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator::DPK1_DF;

//    auto bhv = mgis::behaviour::load(opt, path, name, hyp);
    auto mgipara = lolita::behaviour::MgisParameter("Temperature", [] (lolita::geometry::Point const & p, lolita::real const & t) {return 1.0;});
//    auto bhvs = lolita::behaviour::MgisBehaviours{
//            {"Middle", lolita::behaviour::MgisBehaviour("Displacement", path, name, hyp, {mgipara})}
//    };
    auto bhv_u = lolita::behaviour::MgisBehaviour("Displacement", "SQUARE", path, name, hyp, {mgipara});
    auto bhv_d = lolita::behaviour::MgisBehaviour("Damage", "SQUARE", path, name, hyp, {mgipara});

    auto load = lolita::finite_element::Load("Displacement", "SQUARE", 0, 0, [] (lolita::geometry::Point const & p, lolita::real const & t) {return 1.0;}, lolita::finite_element::Loading::Natural);
    auto load_top_x = lolita::finite_element::Load("Displacement", "TOP", 0, 0, [] (lolita::geometry::Point const & p, lolita::real const & t) {return 1.0;}, lolita::finite_element::Loading::Constraint);
    auto load_top_y = lolita::finite_element::Load("Displacement", "TOP", 1, 0, [] (lolita::geometry::Point const & p, lolita::real const & t) {return 1.0;}, lolita::finite_element::Loading::Constraint);

//    auto loads = lolita::finite_element::Loads{
//            {"Middle", lolita::finite_element::Load2("Displacement", 0, 0, [] (lolita::geometry::Point const & p, lolita::real const & t) {return 1.0;}, lolita::finite_element::Loading::Natural)}
//    };

//    std::map<lolita::mesh::Domain, lolita::real> mymap = {
//            {lolita::mesh::Domain(1, "SQUARE"), 1.0}
//    };

    auto bhvvv = lolita::behaviour::Behaviour<u_voce>();

//    auto bhvholder = lolita::behaviour::BehaviourHolder(bhvvv);

    auto file_path = "";
    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";
//    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/perforated_strip_huge.msh";

    auto msh2 = lolita::core::mesh::MeshParser<lolita::mesh::Format::Gmsh, domain, hho_u0, hho_d>(file_path, {load, load_top_x, load_top_y}, {bhv_u, bhv_d});

//    std::cout << msh2.mesh_data_.elements_.getElements<2, 0>()["345"]->getElement<0>()->operators_[0] << std::endl;
//    std::cout << msh2.mesh_data_.elements_.getElements<2, 0>()["345"]->getComponentIndex<0, 0>(2) << std::endl;

//    std::cout << "dom :res : ->> " << * msh2.mesh_data_.elements_.getElements<2, 0>()["12"]->getElement<0>()->domains_[0] << std::endl;

    if (file_path == "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh") {
        auto & element = msh2.mesh_data_.elements_.getElements<2, 0>()["345"];
        std::cout << element->getElement<0>()->operators_[0] << std::endl;
        element->getElement<0>()->integrateBehaviour();
        auto & seg = msh2.mesh_data_.elements_.getElements<1, 0>()["35"];
        seg->getElement<0>()->getCurrentDiameters();
        seg->getElement<0>()->getNeighbours<0, 0>()[0]->getCurrentDiameters();
        seg->getElement<0>()->getNeighbours<1, 0>()[0]->getCurrentDiameters();
        seg->getElement<0>()->getComponents<0, 0>()[0]->getCurrentDiameters();
        std::cout << "u_i : 0, 0 :" << seg->getElement<0>()->degrees_of_freedom_[0][0].unknown_index_ << std::endl;
        std::cout << "u_i : 0, 1 :" << seg->getElement<0>()->degrees_of_freedom_[1][0].unknown_index_ << std::endl;
        for (auto const & dom : seg->getElement<0>()->domains_) {
            std::cout << "dom : " << * dom << std::endl;
        }
//        for (int i = 0; i < 2; ++i) {
//            seg->getElement<0>()->mesh_domains_.push_back(std::make_shared<std::basic_string<lolita::character>>("Hello"));
////            seg->getElement<0>()->domains_.push_back(std::make_shared<lolita::index>(1));
//            seg->getElement<0>()->neighbours_;
//        }
//        std::cout << "dom : " << seg->getElement<0>()->domains_[0] << std::endl;
    }



    std::cout << msh2 << std::endl;

    std::cout << msh2.mesh_data_.dof_indices_[0].num_unknowns_ << std::endl;
    std::cout << msh2.mesh_data_.dof_indices_[0].num_bindings_ << std::endl;



}

