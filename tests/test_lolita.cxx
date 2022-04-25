//
// Created by dsiedel on 13/04/2022.
//

#include "gtest/gtest.h"

#include "lolita/lolita_containers.hxx"
#include "lolita/lolita_lolita.hxx"
//#include "lolita/lolita_unordered_map.hxx"
#include "lolita/lolita_core_element_geometry.hxx"
#include "lolita/lolita_core.hxx"
#include "lolita/lolita_core_finite_element.hxx"
#include "lolita/lolita_core_mesh.hxx"

TEST(test_lolita, test_lolita_1)
{

    using namespace lolita;
    using namespace lolita::core;

    auto file_path = "/home/dsiedel/projetcs/feta_9/feta/tests/data/meshes/perforated_strip_huge.msh";
    file_path = "/home/dsiedel/projetcs/feta_7/feta/tests/data/meshes/unit_square_3_cpp.msh";
    auto msh_file = lolita::file::File(file_path);

    auto const constexpr dom = Domain(2, EuclideanFrame::Cartesian);
    auto const constexpr fet = FiniteElement('D', 1, HHO(1, 1), MappingOperator::Gradient, MappingOperator::Identity);
    auto const constexpr met = MixedElement<fet>(Model::Solid);
    auto const constexpr cet = CoupledElement<met, met>(QuadratureRule::Gauss);

//    print("try it :", core::element::FiniteElementUnknown<element::tri_03, fet, met, cet, dom>::numComponentsUnknowns());
//    print("try it :", core::element::FiniteElementUnknown<element::tri_03, fet, met, cet, dom>::numComponentsUnknowns<0>());
//    print("try it :", core::element::FiniteElementUnknown<element::tri_03, fet, met, cet, dom>::numComponentsUnknowns<1>());
//    print("try it :", core::element::FiniteElementUnknown<element::tri_03, fet, met, cet, dom>::numComponentsUnknowns<0, 0>());

    auto const constexpr a = collection::index<Collection<Indx, Real, Intg>, Intg>();
    print("a :", a);
    auto const constexpr tpl = Collection<Indx, Real, Intg>(1, 2, 3);
    auto in = 0;
    auto set_num_components_t = [& in] <typename T>() constexpr mutable {
        if (in > -1) {
            in += T(3);
        }
    };
    collection::apply<Collection<Indx, Real, Intg>>(set_num_components_t);
    print(in);
    auto set_num_components = [& in] (auto const & x) constexpr mutable {
        if (in > -1) {
            in += x;
        }
    };
    auto set_num_components2 = [& in] (auto & x) constexpr mutable {
        x = 0;
    };
    auto tpl2 = Collection<Indx, Real, Intg>(1, 2, 3);
    collection::apply(tpl2, set_num_components2);
    print(tpl2.get<0>(), tpl2.get<1>(), tpl2.get<2>());
    static_assert(isIn(1, 2, 3, 4, 1, 5));

//    core::element::ElemCom<core::element::tri_03, core::element::ElementGeometry>::template Lay<>;
//    core::element::ElemCom<core::element::tri_03, core::element::ElementGeometry>::template Lay<0>;
//    core::element::ElemCom<core::element::tri_03, core::element::ElementGeometry>::template Lay<0, 0>;



//    std::apply(t_add_lambda, std::make_tuple(1, 2, 3));
//    std::for_each(values.begin(), values.end(), add_lambda);
//
//    print("sumres :", sumres);
//
//    print("ordIntegration : ", cet.ordIntegration(dom));
//
////    auto const constexpr coltest = Collection(3, 'A', 1, 3.0);
//    static auto const constexpr ag = Aggregate<Char, Intg, Real>{'A', 1, 3.0};
//    print(aggregate::index(ag, 'A'));
//
//    auto const static constexpr settt = Set(1, 2, 3, 4);
//
////    element::detail::FiniteElementReal<element::pnt_00, 2, FiniteElementMethod::Lagrange, fet, met, dom> telem;
//    auto const constexpr fetlag = FiniteElement('D', 1, Discretization<FiniteElementMethod::Lagrange>(1), MappingOperator::Gradient, MappingOperator::Identity);
////    template<Element E, auto F, auto D, auto M>
//    auto telem = element::FiniteElementOperators<element::tri_03, fet, met, cet, dom>();
//    telem.heyCell();
//    print(telem.getOperator<MappingOperator::Gradient>(0));
//    auto tefac = element::FiniteElementOperators<element::seg_02, fet, met, cet, dom>();
//    tefac.heyFace();
//    struct FiniteElementReal<pnt_00, d, FiniteElementMethod::Lagrange, F, D, M>
//    auto const constexpr col = ag.toCollection();
//    auto const constexpr idx = collection::index(coltest, 'B');
//    static_assert(collection::has<Intg>(ag));
//    static_assert(collection::has(ag, 'A'));
//    print(collection::has<Intg>(ag));

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
//    elem.initialize();
    elem.initialize2();
    //
    //
    //
    print(elem.get<0>().get<0>().module.stabilization);

    using CC = core::element::CC;
    using AA = core::element::AA;
    using BB = core::element::BB;

    auto cc = CC();

//    auto fun = [&] <auto I> (auto const & x) constexpr mutable {
//        print("x is :", x);
//    };

//    auto arrtest = std::array<Indx, 3>{1, 2, 3};

//    std::apply([&](auto const &... x){(..., fun(x));}, arrtest);


//    core::element::CRTPTESTUsed<2>().initialize();

//    using GeomType = core::element::FiniteElementGeometry<core::element::tri_03, dom, cet>;
//    GeomType & geoptype = reinterpret_cast<GeomType &>(elem);
//    print(std::addressof(geoptype));
//    print(std::addressof(elem));
//    print("getCurrentCoordinates() -->");
//    print(geoptype.getCurrentCentroid());
////    elem.init();
////    elem.sayCurrentCentroid();
//    auto ptr_elem = elem.toGeom();
//    print("getCurrentCoordinates()");
//    print(ptr_elem.get().getCurrentCoordinates());
//    print(& elem);
//    print(& ptr_elem.get());
//    Indx idx_test = 9;
//    Intg itg_test = 9;
//    Indx * ptr0;
//    UniquePointer<Indx> ptr1;
//    ptr1.data = & idx_test;
////    ptr0 = reinterpret_cast<Indx *>(&itg_test);
//    ptr0 = & idx_test;
//    print(std::addressof(idx_test));
//    print(std::addressof(* ptr0));
//    idx_test = 3;
//    print(itg_test);
//    print(* ptr0);

//    elem.sayCurrentCentroid();
//    auto ref = elem.toGeometry();
//    print(ref.getCurrentCentroid());

    printMesh(msh);


}

