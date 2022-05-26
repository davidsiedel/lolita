//
// Created by dsiedel on 03/05/22.
//

#include <span>
#include <ranges>
#include "gtest/gtest.h"

#include "lolita/lolita_user.hxx"
#include "lolita/lolita_element.hxx"
#include "lolita/lolita_core_finite_element.hxx"
#include "lolita/lolita_core_mesh.hxx"
#include "lolita/lolita_core_mesh2.hxx"

TEST(test_lolita_2, test_lolita_2)
{

//    using namespace lolita;
//    using namespace lolita::core;

    auto constexpr domain = lolita::geometry::Domain("Middle", 2, lolita::geometry::Frame::Cartesian);
//    auto constexpr displacement = lolita::field::DegreeOfFreedom("Displacement", 4, lolita::field::Mapping::Gradient);
//    auto constexpr displacement_field = lolita::field::Tensor("Displacement", 1);
    auto constexpr u_unknown = lolita::field::Unknown("Displacement", 1, lolita::field::Mapping::Gradient);
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

    std::cout << "sizeof(Middle)" << sizeof(std::basic_string_view<lolita::character>("MiddleMiddleMiddleMiddleMiddleMiddleMiddleMiddleMiddleMiddle")) << std::endl;
    std::cout << "sizeof(Middle)" << sizeof(std::basic_string<lolita::character>("MiddleMiddleMiddleMiddleMiddleMiddleMiddleMiddleMiddleMiddle")) << std::endl;
//    opt.stress_measure = mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
//    opt.tangent_operator = mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator::DPK1_DF;

//    auto bhv = mgis::behaviour::load(opt, path, name, hyp);
    auto mgipara = lolita::behaviour::MgisParameter("Temperature", [] (lolita::geometry::Point const & p, lolita::real const & t) {return 1.0;});
//    auto bhvs = lolita::behaviour::MgisBehaviours{
//            {"Middle", lolita::behaviour::MgisBehaviour("Displacement", path, name, hyp, {mgipara})}
//    };
    auto bhv = lolita::behaviour::MgisBehaviour("Displacement", "SQUARE", path, name, hyp, {mgipara});

    auto load = lolita::finite_element::Load("Displacement", "SQUARE", 0, 0, [] (lolita::geometry::Point const & p, lolita::real const & t) {return 1.0;}, lolita::finite_element::Loading::Natural);

//    auto loads = lolita::finite_element::Loads{
//            {"Middle", lolita::finite_element::Load2("Displacement", 0, 0, [] (lolita::geometry::Point const & p, lolita::real const & t) {return 1.0;}, lolita::finite_element::Loading::Natural)}
//    };

//    std::map<lolita::mesh::Domain, lolita::real> mymap = {
//            {lolita::mesh::Domain(1, "SQUARE"), 1.0}
//    };

    auto bhvvv = lolita::behaviour::Behaviour<u_voce>();

    auto bhvholder = lolita::behaviour::BehaviourHolder(bhvvv);

    lolita::matrix::Vector<lolita::real> vector(200);
    vector.setZero();

    auto t0 = std::chrono::high_resolution_clock::now();

    vector.segment<20>(40) * vector.segment<20>(50).transpose();

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration<lolita::real>(t1 - t0);

    std::cout << "duration : " << dt.count() << " s" << std::endl;

    lolita::matrix::Vector<lolita::real, 20> vector2;
    vector2.setZero();
    t0 = std::chrono::high_resolution_clock::now();

    vector2 * vector2.transpose();

    t1 = std::chrono::high_resolution_clock::now();
    dt = std::chrono::duration<lolita::real>(t1 - t0);

    std::cout << "duration : " << dt.count() << " s" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();

    vector.segment(20, 40) * vector.segment(20, 40).transpose();

    t1 = std::chrono::high_resolution_clock::now();
    dt = std::chrono::duration<lolita::real>(t1 - t0);

    std::cout << "duration : " << dt.count() << " s" << std::endl;

//    lolita::finite_element::LoadFunction f;
//    std::cout << "fun res : " << f(lolita::geometry::Point::Ones(), 1.0) << std::endl;


//    std::cout << lolita::core::base::neighbourPosition<lolita::core::element::pnt_00, lolita::core::element::seg_02>()[2] << std::endl;
//    std::cout << lolita::core::base::neighbourPosition<lolita::core::element::pnt_00, lolita::core::element::tri_03>()[2] << std::endl;
//    std::cout << lolita::core::base::neighbourPosition<lolita::core::element::pnt_00, lolita::core::element::tet_04>()[2] << std::endl;

    auto file_path = "";
//    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/perforated_strip_huge.msh";
    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";
    auto msh2 = lolita::core::mesh2::MeshParser<lolita::mesh::Format::Gmsh, domain>("/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh");
    auto msh_file = lolita::utility::File(file_path);
    auto msh = lolita::core::mesh::MeshBase<lolita::mesh::Format::Gmsh, domain, hho_u0, hho_d>(msh_file, {load}, {bhv});
    auto & elem = * msh.elements_.getElements<2, 0>()["345"];
//    auto & elem = * std::get<0>(std::get<2>(msh.elements_))["345"];
//    std::cout << elem.getCurrentCoordinates() << std::endl;
//    std::cout << elem.getReferenceCoordinates() << std::endl;
////
//    for (auto const & [n_hash, n] : msh.elements_.getElements<0, 0>()) {
////    for (auto const & [n_hash, n] : std::get<0>(std::get<0>(msh.elements_))) {
//        std::cout << "--" << n_hash << std::endl;
//        std::cout << n->getCurrentCoordinates() << std::endl;
//    }
//
//    std::cout << "std::tuple_size_v<std::tuple<>> : " << std::tuple_size_v<std::tuple<>> << std::endl;
//
//    if constexpr(lolita::index(0) < -1) {
//        std::cout << "hello" << std::endl;
//    }

    std::cout << " msh.unknown_indices[0].unknown_index_ : " << msh.unknown_indices[0].unknown_index_ << std::endl;
    std::cout << " msh.unknown_indices[0].binding_index_ : " << msh.unknown_indices[0].binding_index_ << std::endl;

    std::cout << elem.getElement<0>()->getCurrentCoordinates() << std::endl;
    std::cout << elem.getCurrentCoordinates() << std::endl;
    std::cout << elem.getCurrentCentroid() << std::endl;
    std::cout << elem.getCurrentDiameters() << std::endl;
    std::cout << elem.getCurrentQuadraturePoint<lolita::finite_element::Quadrature::Gauss, 2>(0) << std::endl;

//    std::get<0>(std::tuple<>());
//    static_assert(lolita::index(0) < -1);

//    std::cout << lolita::core::element::numNeighbours<lolita::core::element::pnt_00, 1, 2>() << std::endl;

    std::cout << msh << std::endl;

//    std::cout << u_unknown.tensor_.cardinality(domain).rows_ << u_unknown.tensor_.cardinality(domain).cols_ << std::endl;
//
//    std::cout << "index : " << lolita::behaviour::MgisBehaviour2(u_unknown).getUnknownIndex<u_unknown>() << std::endl;
//
//    auto constexpr newlabel = lolita::utility::label<lolita::character>(
//            std::forward<std::basic_string_view<lolita::character>>(lolita::utility::readLabel(u_unknown.tensor_.tag_)),
//            std::forward<std::basic_string_view<lolita::character>>(lolita::utility::readLabel(u_unknown.tensor_.tag_))
//    );
//
//    std::cout << lolita::utility::readLabel(newlabel) << std::endl;
//
//    std::cout << std::string(lolita::utility::readLabel(u_unknown.tensor_.tag_)).append("HJU") << std::endl;

//    auto constexpr agg = lolita::utility::Aggregate<long, double, int>{1, 2, 3};
//
////    std::cout << lolita::utility::get<0u>(agg) << std::endl;
//
////    lolita::utility::Mult<1, 2, 3>();
//
//    auto constexpr inpt = std::tuple<int, long, unsigned, long , char>(1, 2, 3, 4, 'A');
//
////    auto constexpr outp = lolita::utility::subtuple<1>(inpt);
//
////    auto constexpr outp = lolita::utility::subtuple2(inpt);
////    auto constexpr outp2 = lolita::utility::expandD(inpt);
////    auto constexpr outp2 = lolita::utility::expandTuple<1>(inpt);
//
//    auto constexpr enfin = lolita::utility::tupleSlice<1, 2>(inpt);
//
////    std::apply(([&](auto x) {std::cout << x <<std::endl;}; ...), outp);
//    std::apply([&](auto&&... args) {((std::cout << args << std::endl), ...);}, enfin);

//    lolita::geometry::Frame::Cartesian;
//    lolita::geometry::Domain::Cartesian;

//    auto constexpr domain = Domain("Middle", 3, Domain::AxiSymmetric);
//    auto constexpr displacement = DegreeOfFreedom("Displacement", 0, Mapping::Gradient);
//    auto constexpr damage = DegreeOfFreedom("Damage", 0, Mapping::Gradient, Mapping::Identity);
//    auto constexpr young = MaterialProperty("YoungModulus");
//    auto constexpr plastic_strain = InternalVariable("PlasticStrain", 2);
//    auto constexpr material = BehaviourData("Voce", displacement, damage, plastic_strain, young);
//    auto constexpr dis = Discretization::HybridHighOrder(1, 1);
//
//    auto integrator = [] ()
//    constexpr
//    {
//
//    };
//
//    auto a = lolita::index(1);
//    auto b = lolita::index(5);
////    std::cout << 3 / 2 << std::endl;
////    std::cout << 3 % 2 << std::endl;
//    std::cout << displacement.cardinality(domain).rows_ << displacement.cardinality(domain).cols_ << std::endl;
//
//    auto constexpr bhv = BHV<displacement, damage, plastic_strain, young>("Voce");

//    static_assert(young == MaterialProperty("YoungModulus"));

////    auto static constexpr ttt = lta::agg<char, int, double>{};
//    auto constexpr static ttt = lolita::utility::Aggregate2<char, std::size_t , std::size_t>('A', 1, 2);
//    auto constexpr static ttt2 = lolita::utility::Aggregate<char, std::size_t, std::size_t>{'A', 1, 2};
//
//    std::cout << ttt.get<0>() << std::endl;
//    std::cout << ttt.get<1>() << std::endl;
//    std::cout << ttt.get<2>() << std::endl;
//
//    auto constexpr label = lolita::utl::Label("hello");

//    Unknown('U', 1, Mapping::Type::G, Mapping::Type::H);
//    Field('E', 1, Field::Type::Internal);
//    static_assert(displacement == DegreeOfFreedom("Displacement", 1, Mapping::Gradient));
//    static_assert(Mapping(Mapping::Gradient, 1, 1) == Mapping(Mapping::Gradient, 1, 1));
//    static_assert(std::derived_from<decltype(damage2), detail::FieldBase>);
//    static_assert(std::derived_from<decltype(damage), detail::FieldBase>);



//    lolita::real a0 = 3;
//    lolita::real const & a = a0;
//    lolita::integer i0 = 3;
//    lolita::integer const i = i0;
//
////    static_assert(std::is_floating_point_v<decltype(a)>);
//    static_assert(std::is_integral_v<decltype(i)>);
//
//    lolita::numerics::pow(1, 2);
//
//    std::cout << lolita::numerics::abs(lolita::index(2)) << std::endl;
//
//    auto constexpr tpl = std::tuple(1, 2);
//    auto constexpr s = S<std::get<0>(tpl)>();
//    auto constexpr exp = expand<S, damage, displacement>();

//    struct Se {
//
//        lolita::real
//        operator()(
//                geometry::Point const & p
//        )
//        const
//        {
//            return p(0);
//        }
//
//    };
//
//    auto constexpr see = Se();
//
//    auto fun = [] (lolita::real, lolita::real, lolita::real) constexpr {return lolita::real (1); };
//
//    typedef double feE (double,double,double) const;
//
//    static_assert(std::is_same_v<decltype(fun), feE>);
//
//    std::function<lolita::real(geometry::Point)> fun2;
//
//    static_assert(ValueSetter<Se>);
//    static_assert(ValueSetter<decltype(fun2)>);
//    auto name = utility::setName("helooo");
//    std::cout << utility::getName(name) << std::endl;
//    std::cout << Domain::Cartesian << std::endl;
//    std::cout << Domain::AxiSymmetric << std::endl;

//    auto constexpr ffld = core::Field::Structure('A', 1, Mapping::Type::Divg, Mapping::Type::Grad);

//    lolita::Test<ttt> h;

}

