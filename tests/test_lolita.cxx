//
// Created by dsiedel on 13/04/2022.
//

#include <span>
#include <ranges>
#include "gtest/gtest.h"

#include "lolita/lolita_core_mesh.hxx"

TEST(test_lolita, test_lolita_1)
{

    using namespace lolita;
    using namespace lolita::core;

    /*
     *
     */
//    static_assert(Machine2<SSS>);
//    static_assert(Machine2<SSS>);
    auto const static constexpr label = Label();
    constexpr std::string_view my_str = "hello, world";
    auto const constexpr lbl2 = Label2("hello, world", 1, 2);
    print(lbl2.sv, lbl2.i, lbl2.j);
    static_assert(lbl2.sv == "hello, world");
//    auto const constexpr strvv = StrView("hello, world");
    StrView const constexpr strvv = "hello, world";
    auto const constexpr strv22 = StrView2Char("hello");
    static_assert(strv22 == StrView2Char("hello"));
    print(strv22);
//    auto constexpr strv2 = StrView2("hello", 2);
//    print(strv2.i, strv2._M_str);
//    static_assert(strv2._M_str == "hello");
    StrTest<strv22> mytest("hello", 1);
    print(mytest.v);
//    static_assert(strvv == "hello, world");
//    print(strvv);
    /*
     *
     */
    auto file_path = "";
//    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/perforated_strip_huge.msh";
    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";
    auto msh_file = lolita::file::File(file_path);
    /*
     *
     */
    auto const constexpr energy_field = AuxiliaryField('E', 1, FLDTyp::Internal);
    auto const constexpr dom = Domain(EuclideanFrame::Cartesian, 2);
    auto const constexpr displacement_field = StructuralField('U', 1, MappingOperator::Gradient);
    auto const constexpr pressure_field = StructuralField('P', 0, MappingOperator::Identity);
    auto const constexpr damage_field = StructuralField('D', 0, MappingOperator::Gradient, MappingOperator::Identity);
    auto const constexpr displacement_bhv = BHV2(energy_field, displacement_field, pressure_field);
    auto const constexpr damage_bhv = BHV2(energy_field, damage_field);
    auto const constexpr fet3 = FiniteElement2(damage_field, HHO(1, 1), QuadratureRule::Gauss);
    auto const constexpr ets = Elements<fet3>();
    static_assert(FieldType<decltype(displacement_field)>);
    static_assert(FieldType<decltype(energy_field)>);
    /*
     *
     */
    auto loads = Array<Load<dom>>{
//        Load<dom>('D', "SQUARE", 0, 0, LoadComponent<dom>([] (auto const &, auto const &) { return Real(0); }, LoadType::Constraint)),
            Load<dom>('D', "TOP", 0, 0, LoadComponent<dom>([] (auto const &, auto const &) { return Real(0); }, LoadType::Constraint)),
            Load<dom>('D', "BOTTOM", 1, 0, LoadComponent<dom>([] (auto const &, auto const &) { return Real(0); }, LoadType::Constraint))
    };
    /*
     *
     */
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
//    auto bhv = Behaviour('D', path, name, HypothesisType::PlaneStrain, StrainType::SmallStrain, "SQUARE");
    auto bhv = SharedPointer<mgis::behaviour::Behaviour>(mgis::behaviour::load(path, name, hyp));
//    mgis::behaviour::MaterialDataManager(bhv.get(), 2);
//    auto gfhjk = UniquePointer<mgis::behaviour::MaterialDataManager>(mgis::behaviour::MaterialDataManager(bhv.get(), 2));
    auto bhvs1 = Behaviour('D', 1, bhv, "SQUARE");
    auto bhvs2 = Behaviour('U', 1, bhv, "SQUARE");

//    auto const constexpr qad = Quadrature(QuadratureRule::Gauss, 2);
////    auto const constexpr fet = FiniteElement('D', 1, HHO(1, 1), qad, MappingOperator::Gradient, MappingOperator::Identity);
////    auto const constexpr fet = FiniteElement('D', 1, HHO(1, 1), QuadratureRule::Gauss, 2, MappingOperator::Gradient, MappingOperator::Identity);
//    auto const constexpr fet = FiniteElement('D', 1, HHO(1, 1), QuadratureRule::Gauss, MappingOperator::Gradient, MappingOperator::Identity);
//    auto const constexpr met = MixedElement<fet>(Model::Solid);
//    auto const constexpr cet = CoupledElement<met, met>(QuadratureRule::Gauss);
//    auto const constexpr fes = Collection(fet);
//    auto const constexpr fes2 = std::tuple(fet);




    auto msh = mesh::MeshBase<MeshFormatType::Gmsh, dom, ets>(msh_file, loads, {bhvs1});
    auto & elem = msh.elements_.get<2>().get<0>().get("345").get();

//    print("load res :", elem.get<0>().get().loads.get(0, 0).get().getImposedValue(Vector<Real, 2>::Ones(), 1.0));

    print("unknown_index : ", msh.unknown_indices.get(0).unknown_index);
    print("binding_index : ", msh.unknown_indices.get(0).binding_index);

    print("ns :");
//    print(elem.get<0>().get().ns.get(0,0));

    elem.get<0>().get().operators.get(0);
    print(elem.get<0>().get().operators.get(0));

    auto myvec = Vector<Real>();

    myvec.resize(4);
    myvec.setZero();
    print(myvec);


     print(msh);

}

