//
// Created by dsiedel on 03/05/22.
//

#include <span>
#include <ranges>
#include "gtest/gtest.h"

#include "lolita/lolita.hxx"

TEST(test_lolita_2, test_lolita_2)
{

    using namespace lolita;
    using namespace core;

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
    Domain("Middle", 1, Domain::AxiSymmetric);
//    TensorField("Displacement", 1);
    auto constexpr displacement = DegreeOfFreedom("Displacement", 1, Mapping::Gradient);
    auto constexpr damage = DegreeOfFreedom("Damage", 0, Mapping::Gradient, Mapping::Identity);
    auto constexpr damage2 = detail::FieldBase("Damage", 0);
    static_assert(FieldBaseDerivedType<decltype(damage2)>);
    static_assert(std::derived_from<decltype(damage2), detail::FieldBase>);
    static_assert(std::derived_from<decltype(damage), detail::FieldBase>);

    auto constexpr tpl = std::tuple(1, 2);
    auto constexpr s = S<std::get<0>(tpl)>();
    auto constexpr exp = expand<S, damage, displacement>();

    auto fun = [] () constexpr {};
//    auto name = utility::setName("helooo");
//    std::cout << utility::getName(name) << std::endl;
//    std::cout << Domain::Cartesian << std::endl;
//    std::cout << Domain::AxiSymmetric << std::endl;

//    auto constexpr ffld = core::Field::Structure('A', 1, Mapping::Type::Divg, Mapping::Type::Grad);

//    lolita::Test<ttt> h;

}

