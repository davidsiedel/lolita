//
// Created by dsiedel on 13/04/2022.
//

#include <span>
#include "gtest/gtest.h"

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

    auto const constexpr dom = Domain(EuclideanFrame::Cartesian, 2);
    auto const constexpr qad = Quadrature(QuadratureRule::Gauss, 2);
    auto const constexpr fet = FiniteElement('D', 1, HHO(1, 1), MappingOperator::Gradient, MappingOperator::Identity);
    auto const constexpr met = MixedElement<fet>(Model::Solid);
    auto const constexpr ets = Elements<fet>();
    auto const constexpr cet = CoupledElement<met, met>(QuadratureRule::Gauss);
    constexpr auto str = "hello world abc programming";
//    auto constexpr ndarray = std::array<std::array<Real, 2>, 2>{{{1, 2}, {3, 4}}};
    auto constexpr ndarray = std::array<std::array<Real, 2>, 2>{1, 2, 3, 4};
    auto spn = std::span(ndarray);
    ndarray[0][0];

    auto set_components = [&] <auto K, auto I, auto J> (auto & set_components_impl)
    mutable
    {
        auto constexpr szeK = 2;
        auto constexpr szeI = 1;
        auto constexpr szeJ = 1;
        print("hello test :", K, I, J);
        if constexpr (J < szeJ) {
            set_components_impl.template operator()<K, I, J + 1>(set_components_impl);
        }
//        if constexpr (I < sze && J == sze) {
        else if constexpr (I < szeI) {
            set_components_impl.template operator()<K, I + 1, 0>(set_components_impl);
        }
//        if constexpr (K < sze && I == sze && J == sze) {
        else if constexpr (K < szeK) {
            set_components_impl.template operator()<K + 1, 0, 0>(set_components_impl);
        }
    };
    set_components.template operator()<0, 0, 0>(set_components);


    LoadL<dom>('D', "MIDDLE", 0, 0, LoadComponent<dom>([] (auto const &, auto const &) { return Real(0); }, LoadType::Constraint));

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

//    auto elem_test = element::FiniteElement<element::pnt_00, dom, met>();
//    elem_test.coordinates;
//    auto elem_test = SharedPointer<element::FiniteElementUnknownMod<element::tri_03, dom, fet>>(element::FiniteElementUnknownMod<element::tri_03, dom, fet>());
    auto elem_test = SharedPointer<element::FiniteElementUnknownMod<element::tri_03, dom, fet>>(element::FiniteElementUnknownMod<element::tri_03, dom, fet>());
    elem_test.get().components.get<0>().get<0>().get(0).ptr.get().unknowns;
    auto elem_cmp = SharedPointer<element::FiniteElementUnknownMod<element::seg_02, dom, fet>>(element::FiniteElementUnknownMod<element::seg_02, dom, fet>());
    elem_test.get().components.get<0>().get<0>().get(0) = element::FiniteElementUnknownMod<element::tri_03, dom, fet>::Component<element::seg_02>(
            elem_cmp, 1
    );

    auto ptr_struct_test = SharedPointer<element::ATest<2>>(
            element::ATest<2>()
    );
//    print("here 1:", ptr_struct_test.get().ptr.exists());
//    print("here 2:", ptr_struct_test.get().ptr.get().ptr.exists());
//    print("here 3 :", ptr_struct_test.get().ptr.get().ptr.get().ptr.exists());

    auto ptr = SharedPointer<Indx>::template make<Indx>();
    auto ptr2 = SharedPointer<Indx>::template make<Indx>(3);
//    auto ptr1 = SharedPointer<Indx>::make();
//    auto ptr12 = SharedPointer<Indx>::make(3);

    print("ptr :", ptr.get());
    print("ptr2 :", ptr2.get());

    auto ptrtest = std::make_shared<Indx>(2);


    auto msh = mesh::MeshBase<MeshFormatType::Gmsh, dom, ets>(msh_file, {});
    auto & elem = msh.elements_.get<2>().get<0>().get("345").get();
    auto & elem2 = msh.elements_.get<2>().get<0>().get("345");
//    elem.initialize();
    auto & aa = elem2.get().components.template get<0>().template get<0>().get(0).ptr.get().get<0>();
//    print("aa read");
    auto & bb = elem2.get().template get<0>().get().components.template get<0>().template get<0>().get(0).ptr;
//    print("bb read");
//    print("bb exists : ", elem.template get<0>().components.template get<0>().template get<0>().get(0).exists());
//    aa.itemp = 2;
//    aa.unknowns.get(0, 0);
//    print("aa itemp read");
//    bb.itemp;
////    bb.unknowns.get(0, 0);
//    print("bb itemp read");
//    bb = aa;

    print(elem2.get().get<0>().get().components.get<0>().get<0>().get(0).ptr.get().getCurrentCoordinates());
    print(elem2.get().components.get<0>().get<0>().get(0).ptr.get().getCurrentCoordinates());
    print(elem2.get().get<0>().get().getCurrentCoordinates());
    print(elem2.get().getCurrentCoordinates());


//    print(elem.getCurrentCoordinates());
//    print(elem.getBasisEvaluation<Basis(BasisName::Monomial, 2)>(Vector<Real, 2>::Zero()));
//    elem.components.template get<0>().template get<0>().get(0).orientation;
//    elem.components.template get<0>().template get<0>().get(0).get();
//    elem.components.template get<0>().template get<0>().get(0).get().template get<0>().itemp;
////    print("here coords :");
////    elem.template get<0>().getBasisEvaluation<Basis(BasisName::Monomial, 2)>(Vector<Real, 2>::Zero());
//    print("fe compo exists :", elem.template get<0>().components.template get<0>().template get<0>().get(0).exists());
//    print("itemp :", elem.template get<0>().components.template get<0>().template get<0>().get(0).get().itemp);
//    elem.template get<0>().components.template get<0>().template get<0>().get(0).get();
//    print("here coords :");
//    elem.template get<0>().components.template get<0>().template get<0>().get(0).get().loads;
//    print("here coords :");
//    print(elem.template get<0>().components.template get<1>().template get<0>().get(0).get().coordinates.get());


//    elem.initialize();
//    print(elem.get<0>().getCurrentCoordinates());
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
////    elem.getGrads();
//    //
//    //
//    //
////    print(elem.get<0>().get<0>().module.stabilization);
////
//
//    auto vec = Vector<Real, -1>::Zero(0);
//    print(vec);
//    print("vec.size() :", vec.size());
    print(sizeof(elem));
    printMesh(msh);


}

