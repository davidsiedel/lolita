//
// Created by dsiedel on 14/03/2022.
//

#include "gtest/gtest.h"

#include "new/feta_core_element_element_connectivity.hxx"
#include "new/feta_core_field_description.hxx"
#include "new/_feta_shared_pointer.hxx"

#include "new/feta_core_mesh_parser_implementation.hxx"
//#include "new/00_finite_element_description_base.hxx"
//#include "new/mixed_element_description.hxx"
#include <chrono>
#include <span>
//#include "new/00_finite_element_data.hxx"
#include "new/_feta_matrix.hxx"
#include "new/feta_core_element_fe_degree_of_freedom.hxx"
//#include "new/_feta_tensor.hxx"
#include "new/feta_core_element_element_reference.hxx"
#include "new/feta_core_fe_description.hxx"
//#include <sys/resource.h>
//#include <sys/time.h>

TEST(test_gmsh_new, test_miscs) {

    using namespace feta::core::element;
    using namespace feta::core;
    using namespace feta;
//
//
////    auto const mappings = Mappings(MappingDescription2(Mapping::Gradient, 1), MappingDescription2(Mapping::Gradient, 1));
////    auto const mappings2 = Mappings();
////    auto const feoptshho = DiscretizationDescription<Discretization::HybridHighOrder>(FieldMapping(Field::Scalar, Mapping::Gradient), 1, 2);
//
////    auto const disc = FieldMapping(Field::Scalar, Mapping::Identity);
////    auto const disc2 = FieldMapping(Field::Scalar, Mapping::Gradient);
////    auto const mappings2 = Mappings(core::MappingDescription(Mapping::Gradient, 1), core::MappingDescription(Mapping::Identity, 1));
//
////    DiscreteFieldDescription();
//
//    auto const constexpr field_desc = FieldDescription::getMappingFieldDescription(FieldDescription(Field::Vector, 2), Mapping::Identity);
//
//    DiscreteFieldDescription(Field::Vector, DiscreteFieldComponent(Body::Cell, 1), DiscreteFieldComponent(Body::Cell, 2));
//
//    auto const hhoset = Setter<FiniteElementType::HybridHighOrder>();
//
////    Matrix<Real, 2, 2> mymatrix = Matrix<Real, 2, 2>::Zero();
//    auto mymatrix = Matrix<Real, 2, 5>();
//    auto bary = geometry::getBarycenter(mymatrix);
//    mymatrix.setZero();
//
//
//
//    Load2<2, field_desc>();
//    Load2<2, field_desc>(LoadComponent<2>(), LoadComponent<2>());
//
//    auto exponents_values = Array<Indx, 2, 2>();
//
////    auto const constexpr gh = FiniteElementDescription2('D', FiniteElementType::HybridHighOrder, DiscreteFieldDescription(Field::Scalar, 2, DiscreteFieldComponent(Structure::Cell, 1), DiscreteFieldComponent(Structure::Cell, 2)), MappingDescription(Mapping::Gradient, 1));
//
//    Indx const constexpr DIM = 2;
//
//    DavidDerived::Type klopm;
//    DavidDerived::Type3 ksjhsh;
//
//    using DamageField = FiniteElementDescription<
//            FiniteElementType::HybridHighOrder,
//            Field::Scalar,
//            Mapping::Gradient,
//            Mapping::Identity
//    >;
//
//    using DisplacementField = FiniteElementDescription<
//            FiniteElementType::HybridHighOrder,
//            Field::Vector,
//            Mapping::SmallStrain
//    >;
//
//    DisplacementField constexpr displacement_field('U', 2, 1, 1, 1);
//
//    DamageField constexpr damage_field('D', 2, 1, 1, 1, 1);
//
//    static_assert(Body::Cell == static_cast<Body>(0));
//
//    using MixedField = MixedElementDescription3<
//            damage_field,
//            displacement_field
//    >;
//
////    ElementReference<tri_3>::node_connectivity.get<0>().get<0>().get(0, 0);
////    ElementReference<qua_4>::node_connectivity.get<1>().get<0>().get(3, 1);
////    ElementReference<seg_2>::node_connectivity.get<1>().get<0>();
////
////    auto test_tuple_mono = Collection(Collection(Collection(Real(0))));
////
////    test_tuple_mono.get<0>().get<0>().get<0>();
//
////    auto const constexpr mixed_new = MixedElementDescriptionNew(2, FrameType::Cartesian, Quadrature::Gauss, damage_field);
//
//    MixedField constexpr mixed_element_description(2, FrameType::Cartesian, Quadrature::Gauss);
//
//    using InteriorDomain = MeshInteriorDomain<mixed_element_description>;
//
//    auto zero_load = [] (Matrix<Real, DIM> const & point_arg, Real const & time_arg) { return 0.0; };
//
//    auto const lll = Load(
//            LoadComponent<DIM>(zero_load, 0, LoadType::Volumetric),
//            LoadComponent<DIM>(zero_load, 0, LoadType::Volumetric),
//            LoadComponent<DIM>(zero_load, 0, LoadType::Volumetric)
//    );
//
//    InteriorDomain mesh_interior_domain_component = InteriorDomain(
//            "MIDDLE",
//            "/home/dsiedel/projetcs/feta_9/feta/tests/data/behaviour/src/libBehaviour.so",
//            "Voce",
//            HypothesisType::PlaneStrain,
//            StrainType::LargeStrain,
//            Load1<damage_field>(
//                    LoadComponent1<damage_field>(zero_load, 0, LoadType::Volumetric)
//            ),
//            Load1<displacement_field>(
//                    LoadComponent1<displacement_field>(zero_load, 0, LoadType::Volumetric),
//                    LoadComponent1<displacement_field>(zero_load, 1, LoadType::Volumetric)
//            )
//    );
//
//    InteriorDomain mesh_boundary_domain = InteriorDomain(
//            "TOP",
//            Body::Face,
//            Load1<damage_field>(
//                    LoadComponent1<damage_field>(zero_load, 0, LoadType::Neumann)
//            ),
//            Load1<displacement_field>(
//                    LoadComponent1<displacement_field>(zero_load, 0, LoadType::Dirichlet),
//                    LoadComponent1<displacement_field>(zero_load, 1, LoadType::Neumann)
//            )
//    );
//
//    InteriorDomain mesh_boundary_domain2 = InteriorDomain(
//            "BOTTOM",
//            Body::Face,
//            Load1<damage_field>(
//                    LoadComponent1<damage_field>(zero_load, 0, LoadType::Neumann)
//            ),
//            Load1<displacement_field>(
//                    LoadComponent1<displacement_field>(zero_load, 0, LoadType::Dirichlet),
//                    LoadComponent1<displacement_field>(zero_load, 1, LoadType::Neumann)
//            )
//    );
//
//    MeshInteriorDomains<mixed_element_description> mesh_interior_domains = {
//            MeshInteriorDomains<mixed_element_description>::Item{
//                    "MIDDLE",
//                    SharedPointer<MeshInteriorDomain<mixed_element_description>>(mesh_interior_domain_component)
//            }
//    };

    Matrix<Real, 3> matrix_segtet;
    print("matrix_segtet.segment<0>(0) :", matrix_segtet.segment<0>(0));

    std::string file_path = "/home/dsiedel/projetcs/feta4/feta/tests/test_mesh/test_gmsh/data/cook_01_quadrangles_structured.msh";
    file_path = "/home/dsiedel/projetcs/feta5/feta/tests/test_mesh/test_gmsh/data/unit_square_2_cpp.msh";
    file_path = "/home/dsiedel/projetcs/feta5/feta/tests/test_mesh/test_gmsh/data/unit_square_cpp.msh";
    file_path = "/home/dsiedel/projetcs/feta_9/feta/tests/data/meshes/triangles_3.msh";
    file_path = "/home/dsiedel/projetcs/feta_9/feta/tests/data/meshes/unit_square_cpp.msh";
//    file_path = "/home/dsiedel/projetcs/feta_9/feta/tests/data/meshes/unit_square_3_cpp.msh";
    file_path = "/home/dsiedel/projetcs/feta_9/feta/tests/data/meshes/perforated_strip_huge.msh";
    file_path = "/home/dsiedel/projetcs/feta_7/feta/tests/data/meshes/unit_square_3_cpp.msh";
//    file_path = "/home/dsiedel/projetcs/feta_9/feta/tests/data/meshes/cook_01_quadrangles_structured.msh";
//    file_path = "/home/dsiedel/projetcs/feta_6/feta/tests/test_mesh/test_parsers/test_gmsh/data/unit_square_cpp_500.msh";
    feta::file::File<file::FileType::Input> mshfile(file_path);


    auto const constexpr fld_d = DiscretizationDescription('D', Field::Scalar, FiniteElementMethod::HybridHighOrder, Mapping::Gradient, Mapping::Identity);
    auto const constexpr fem_d = FiniteElementDescription<fld_d>(elt::Model::Solid, 1, 2);

    auto const constexpr fld_u = DiscretizationDescription('U', Field::Vector, FiniteElementMethod::HybridHighOrder, Mapping::SmallStrain);
    auto const constexpr fem_u = FiniteElementDescription<fld_u>(elt::Model::Solid, 1, 1, 1);

    auto const constexpr me_ud = MixedElementDescription<fem_u, fem_d>(elt::Model::Solid, 2, EuclideanFrame::Cartesian, Quadrature::Gauss);

    Element<tri_3, me_ud> em;
//    em.components.get<0>().get<0>().get(0).get().get<0>().unknowns.get(0, 0);
//    Nesting nesting;
//    nesting.a;

//    auto unknown_index = 0;
//    auto binding_index = 0;
//    UnknownComponent<tri_3, BasisDescription(Basis::Monomial, 1)> uknc(unknown_index);
//    Unknowns<FieldDescription(Field::Scalar, 2), tri_3, BasisDescription(Basis::Monomial, 1)> ukn(unknown_index, binding_index);
////    ukn.get(0, 0).isBound();
//    print("index now is :", unknown_index);


    ParserBase<me_ud> parser_base;

//    parser_base.element_sets.operator()("J").get<0>().get<0>().data.insert({"J"}, {});

    using ParserType = ParserInterface<MeshFormatType::Gmsh, me_ud>;

    auto start = std::chrono::high_resolution_clock::now();
//    GmshParser<2> ggg(mshfile);

    ParserHelper<feta::core::MeshFormatType::Gmsh> parser_helper(mshfile);
    ParserType ggg(mshfile, {});

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = duration_cast<std::chrono::milliseconds>(stop - start);

    struct PrintIt {

        void
        operator () (
                ParserType const & ggg
        )
        {

            auto elsets = ggg.element_sets;
            auto elems = ggg.element_collection;
            print("START ELEMS");
            print("NODES");
            for (auto const & [key, val] : elems.get<0>().get<0>().data) {
                print(key, val.get().tag);
                for (auto const & ghj: val.get().neighbours.get<0>().get<0>().data) {
                    print(" --> ", ghj.get().getHash());
                }
                for (auto const & ghj: val.get().neighbours.get<1>().get<0>().data) {
                    print(" ----> ", ghj.get().getHash());
                }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
            }
            print("NODES DONE");

            print("SEGS");
            for (auto const & [key, val] : elems.get<1>().get<0>().data) {
                print(key, val.get().tag);
                for (int i = 0; i < val.get().components.get<0>().get<0>().getSize(); ++i) {
                    print("--> ", val.get().components.get<0>().get<0>()(i).get().getHash(), val.get().components.get<0>().get<0>()(i).orientation);
                }
                for (auto const & ghj: val.get().neighbours.get<1>().get<0>().data) {
                    print(" ----> ", ghj.get().getHash());
                }
                for (auto const & ghj: val.get().neighbours.get<0>().get<0>().data) {
                    print(" <-> ", ghj.get().getHash());
                }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
            }
            print("SEGS DONE");

            print("QUAS");
            for (auto const & [key, val] : elems.get<2>().get<1>().data) {
                print(key, val.get().tag);
                for (int i = 0; i < val.get().components.get<0>().get<0>().getSize(); ++i) {
                    print("--> ", val.get().components.get<0>().get<0>()(i).get().getHash(), val.get().components.get<0>().get<0>()(i).orientation);
                }
                for (int i = 0; i < val.get().components.get<1>().get<0>().getSize(); ++i) {
                    print("----> ", val.get().components.get<1>().get<0>()(i).get().getHash(), val.get().components.get<1>().get<0>()(i).orientation);
                }
                for (auto const & ghj: val.get().neighbours.get<0>().get<0>().data) {
                    print(" <-> ", ghj.get().getHash());
                }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
            }
            print("QUAS DONE");

            print("TRIS");
            for (auto const & [key, val] : elems.get<2>().get<0>().data) {
                print(key, val.get().tag);
                for (int i = 0; i < val.get().components.get<0>().get<0>().getSize(); ++i) {
                    print("--> ", val.get().components.get<0>().get<0>()(i).get().getHash(), val.get().components.get<0>().get<0>()(i).orientation);
                }
                for (int i = 0; i < val.get().components.get<1>().get<0>().getSize(); ++i) {
                    print("----> ", val.get().components.get<1>().get<0>()(i).get().getHash(), val.get().components.get<1>().get<0>()(i).orientation);
                }
                for (auto const & ghj: val.get().neighbours.get<0>().get<0>().data) {
                    print(" <-> ", ghj.get().getHash());
                }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
            }
            print("TRIS DONE");

            print("END ELEMS");

            print("START ELSETS");
            for (auto const & [key, val] : elsets.data) {
                print(key);
                for (auto const & [elhash, elptr] : val.get<0>().get<0>().data) {
                    print(elhash);
                }
            }
            for (auto const & [key, val] : elsets.data) {
                print(key);
                for (auto const & [elhash, elptr] : val.get<1>().get<0>().data) {
                    print(elhash);
                }
            }
            for (auto const & [key, val] : elsets.data) {
                print(key);
                for (auto const & [elhash, elptr] : val.get<2>().get<0>().data) {
                    print(elhash);
                }
            }
            print("END ELSETS");
        }

//            auto elsets = ggg.element_sets;
//            auto elems = ggg.elements;
//            print("START ELEMS");
//            print("NODES");
//            for (auto const & [key, val] : elems.get<0>().get<0>().data) {
//                print(key, val.get().tag);
//                for (auto const & ghj: val.get().neighbours.get<0>().get<0>().data) {
//                    print(" --> ", ghj.ptr_neighbour().getHash());
//                }
//                for (auto const & ghj: val.get().neighbours.get<1>().get<0>().data) {
//                    print(" ----> ", ghj.ptr_neighbour().getHash());
//                }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
//            }
//            print("NODES DONE");
//
//            print("SEGS");
//            for (auto const & [key, val] : elems.get<1>().get<0>().data) {
//                print(key, val.get().tag);
//                for (int i = 0; i < val.get().components.get<0>().get<0>().getSize(); ++i) {
//                    print("--> ", val.get().components.get<0>().get<0>()(i).get().getHash(), val.get().components.get<0>().get<0>()(i).orientation);
//                }
//                for (auto const & ghj: val.get().neighbours.get<1>().get<0>().data) {
//                    print(" ----> ", ghj.ptr_neighbour().getHash());
//                }
//                for (auto const & ghj: val.get().neighbours.get<0>().get<0>().data) {
//                    print(" <-> ", ghj.ptr_neighbour().getHash());
//                }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
//            }
//            print("SEGS DONE");
//
//            print("QUAS");
//            for (auto const & [key, val] : elems.get<2>().get<1>().data) {
//                print(key, val.get().tag);
//                for (int i = 0; i < val.get().components.get<0>().get<0>().getSize(); ++i) {
//                    print("--> ", val.get().components.get<0>().get<0>()(i).get().getHash(), val.get().components.get<0>().get<0>()(i).orientation);
//                }
//                for (int i = 0; i < val.get().components.get<1>().get<0>().getSize(); ++i) {
//                    print("----> ", val.get().components.get<1>().get<0>()(i).get().getHash(), val.get().components.get<1>().get<0>()(i).orientation);
//                }
//                for (auto const & ghj: val.get().neighbours.get<0>().get<0>().data) {
//                    print(" <-> ", ghj.ptr_neighbour().getHash());
//                }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
//            }
//            print("QUAS DONE");
//
//            print("TRIS");
//            for (auto const & [key, val] : elems.get<2>().get<0>().data) {
//                print(key, val.get().tag);
//                for (int i = 0; i < val.get().components.get<0>().get<0>().getSize(); ++i) {
//                    print("--> ", val.get().components.get<0>().get<0>()(i).get().getHash(), val.get().components.get<0>().get<0>()(i).orientation);
//                }
//                for (int i = 0; i < val.get().components.get<1>().get<0>().getSize(); ++i) {
//                    print("----> ", val.get().components.get<1>().get<0>()(i).get().getHash(), val.get().components.get<1>().get<0>()(i).orientation);
//                }
//                for (auto const & ghj: val.get().neighbours.get<0>().get<0>().data) {
//                    print(" <-> ", ghj.ptr_neighbour().getHash());
//                }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
//            }
//            print("TRIS DONE");
//
//            print("END ELEMS");
//
//            print("START ELSETS");
//            for (auto const & [key, val] : elsets.data) {
//                print(key);
//                for (auto const & [elhash, elptr] : val.get<0>().get<0>().data) {
//                    print(elhash);
//                }
//            }
//            for (auto const & [key, val] : elsets.data) {
//                print(key);
//                for (auto const & [elhash, elptr] : val.get<1>().get<0>().data) {
//                    print(elhash);
//                }
//            }
//            for (auto const & [key, val] : elsets.data) {
//                print(key);
//                for (auto const & [elhash, elptr] : val.get<2>().get<0>().data) {
//                    print(elhash);
//                }
//            }
//            print("END ELSETS");
//        }
    };

//    PrintIt{}.operator()(ggg);
    auto const & elem = ggg.element_collection.get<2>().get<0>().get("345").get();
    print("getCurrentCoordinates :");
    print(elem.getCurrentCoordinates());
    print("getShapeMappingEvaluation :");
    print(elem.getShapeMappingEvaluation(elem.getCurrentCoordinates().row(0), Matrix<Real, 2>(0, 0)));
    print("getComponentReferencePoint3 :");
    print(elem.getComponentReferenceQuadraturePoint<QuadratureDescription(Quadrature::Gauss, 1), 0, 0>(0, 0));
    print("getBasisEvaluation :");
    print(elem.getBasisEvaluation<BasisDescription(Basis::Monomial, 1)>(Matrix<Real, 2>::Zero()));
    print("isBound() :");
    print(elem.components.get<0>().get<0>().get(0).get().get<0>().degrees_of_freedom.unknowns.getComponentValues(0, 0));
    print("getDimBasis() :");
    print(elem.getDimBasis<BasisDescription(Basis::Monomial, 1)>());
    elem.hey<0>();
    print("stabi :");
    print(elem.template get<0>().stabilization);
    print("map grad :");
    print(elem.getMappingMatrix<0, Mapping::SmallStrain>());
    elem.mapss.get<0>().get<0>();
    print("cell unknowns :");
    print(elem.get<0>().degrees_of_freedom.unknowns.values);
//    print("all unknowns :");
//    print(elem.getUnknowns<0>());
//    print("SmallStrain vals :");
//    print(elem.getMappedField<0, Mapping::SmallStrain>());

    auto const constexpr collection_u_test = Collection(Collection(1));
//    auto valll = collection_u_test.get<0>().get<0>();

//    print(elem.unknowns.template get<0>().get(0, 0).isBound());
//    print("elem.getNumUnknowns<0>()", elem.getNumUnknowns<0>());
//    elem.components.get<0>().get<0>().get(0).get().unknowns.get<0>().get(0, 0).isBound();

//    ElemenT2<tri_3, mem_ud> elemn;
//    mem_ud.getIndex(fem_u);
    static_assert(fem_u == fem_u);
    static_assert(me_ud == me_ud);
//    static_assert(std::array<Indx, 1>{1} == std::array<Indx, 1>{1});
//    static_assert(mem_ud.getIndex(fem_u) == 1);
//    static_assert(aggregate.has('3'));
//    static_assert(aggregate.has('5'));
//    static_assert(aggregate.has(2.0));
//    print("getDimQuadrature :");
//    print(elem.getDimQuadrature<QuadratureDescription(Quadrature::Gauss, 1)>());
    auto ru = Matrix<Real, 3, 3>().setZero();
    auto du = Real(0);
    auto pu = Real();
    print("du :", du);
    print("ru :", ru);
    print("pu :", pu);

    auto const constexpr binom_res = getBinomial(6, 4);
    print(binom_res);

    for (auto i = 0; i < 5; ++i) {
        static_assert(std::is_same_v<decltype(i), int>);
        Indx a = i;
    }

//    print("sizeof(ggg) :", sizeof(ggg));

//    print("--- SEGS :");
//    for (auto const & [key, val] : ggg.elements.get<1>().get<0>().data) {
////        print(key, " - " ,val.tag);
//        print(val.get().tag);
//    }
//    print("--- QUAS :");
//    for (auto const & [key, val] : ggg.elements.get<2>().get<1>().data) {
////        print(key, " - " ,val.tag);
//        print(val.get().tag);
//
//    print("--- TRIS :");
//    for (auto const & [key, val] : ggg.elements.get<2>().get<0>().data) {
////        print(key, " - " ,val.tag);
//        print(val.get().tag);
//    }

//    print(typeid(ggg).name());
//    printStructName(ggg);

    feta::print("time (s) :", duration.count() * 1.e-3);

    ASSERT_EQ(1., 1.);

}
