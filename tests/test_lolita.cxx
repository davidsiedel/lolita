//
// Created by dsiedel on 13/04/2022.
//

#include "gtest/gtest.h"

#include "lolita/lolita_containers.hxx"
#include "lolita/lolita.hxx"
//#include "lolita/lolita_unordered_map.hxx"
#include "lolita/lolita_core_element_geometry.hxx"
#include "lolita/lolita_core.hxx"
#include "lolita/lolita_core_finite_element.hxx"
#include "lolita/lolita_core_mesh.hxx"

TEST(test_lolita, test_lolita_1)
{

    using namespace lolita;
    using namespace lolita::core;

    auto file_path = "";
//    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/perforated_strip_huge.msh";
    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";
    auto msh_file = lolita::file::File(file_path);

    auto const constexpr dom = Domain(2, EuclideanFrame::Cartesian);
    auto const constexpr fet = FiniteElement('D', 1, HHO(1, 1), MappingOperator::Gradient, MappingOperator::Identity);
    auto const constexpr met = MixedElement<fet>(Model::Solid);
    auto const constexpr cet = CoupledElement<met, met>(QuadratureRule::Gauss);

    LoadL<dom>('D', "MIDDLE", 0, 0, LoadComponent<dom>([] (auto const &, auto const &) { return Real(0); }, LoadType::Constraint));

    Behaviour bhv = Behaviour(
            "/home/dsiedel/projetcs/lolita/lolita/tests/data/behaviour/src/libBehaviour.so",
            "Voce",
            HypothesisType::PlaneStrain,
            StrainType::LargeStrain
    );
    mgis::behaviour::MaterialDataManager material_data_manager(bhv.ptr_behaviour.get(), 1);
    matrix::MatMap<Vector<Real, 5>>(material_data_manager.s0.gradients.data()) = Vector<Real, 5>::Zero();
//    material_data_manager.s0.gradients = Vector<Real, 5>::Zero();
//    material_data_manager.s0.set = Vector<Real, 5>::Zero();
//    print(material_data_manager.s0.gradients);
//    mgis::behaviour::integrate()
//    Intg res = mgis::behaviour::integrate(
//            material_data_manager,
//            mgis::behaviour::IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR,
//            0.0,
//            0,
//            1
//    );


    auto msh = mesh::MeshBase<MeshFormatType::Gmsh, dom, cet>(msh_file, {});
    auto & elem = msh.element_collection.get<2>().get<0>().get("345").get();
//    print("tag :", elem.tag);
//    print("unknowns :", elem.get<0>().get<0>().unknowns);
//    print(sizeof(elem.get<0>().get<0>().unknowns));
//    print(sizeof(elem.get<0>().get<0>().operators));
//    print(elem.get<0>().get<0>().operators.get(0).get<0>());
//    print(elem.get<0>().get<0>().module.stabilization);
//    print("");
//    print(elem.get<0>().get<0>().operators.get(0).get<1>());
//    print(sizeof(1.0));
    //
    //
    //
    elem.initialize();
    elem.getGrads();
    //
    //
    //
//    print(elem.get<0>().get<0>().module.stabilization);
//

    auto vec = Vector<Real, -1>::Zero(0);
    print(vec);
    print("vec.size() :", vec.size());
    print(sizeof(elem));
    printMesh(msh);


}

