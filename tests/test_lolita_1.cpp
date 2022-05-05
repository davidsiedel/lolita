//
// Created by dsiedel on 03/05/22.
//

#include <span>
#include <ranges>
#include "gtest/gtest.h"

#include "lolita/lolita.hxx"

TEST(test_lolita_2, test_lolita_2)
{

//    using namespace lolita;
    using namespace lolita::core;

    auto constexpr domain = lolita::geometry::Domain("Middle", 3, lolita::geometry::Frame::Cartesian);
    auto constexpr displacement = lolita::field::DegreeOfFreedom("Displacement", 4, lolita::field::Mapping::Gradient);

    std::cout << displacement.cardinality(domain).rows_ << displacement.cardinality(domain).cols_ << std::endl;

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

