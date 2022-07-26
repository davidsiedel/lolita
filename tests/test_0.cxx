#include "gtest/gtest.h"
#include "lolita/lolita_core_n_111.hxx"

TEST(t0, t0)
{
    
    auto constexpr domain = lolita2::Domain::cartesian(2);
    auto constexpr field = lolita2::Field::vector();
    auto constexpr displacement_field = lolita2::Field::vector();
    auto constexpr damage_field = lolita2::Field::scalar();
    auto constexpr basis = lolita2::Basis::monomial(1);

    // auto lmp = lolita2::Discretization::HJKL(lolita2::Quadrature(lolita2::Quadrature::Gauss, 2), 2);


    auto constexpr u1 = lolita2::Unknown(displacement_field, lolita2::Mapping::smallStrain());
    auto constexpr u2 = lolita2::Unknown(damage_field, lolita2::Mapping::gradient(), lolita2::Mapping::identity());

    auto constexpr bhv = lolita2::Behavior(u1, u2);
    auto constexpr strainSize = lolita2::geometry::BehaviorTraits<bhv>::template strainSize<domain>();
    std::cout << "strainSize : " << strainSize << std::endl;
    // auto constexpr llp = lolita2::S<bhv>();
    // auto constexpr agg = lolita::utility::Aggregate(1, '2', 3.0);
    // auto constexpr tpl = agg.asTuple();
    // std::cout << "tuple : " << std::endl;
    // std::apply([&](auto&&... args) {((std::cout << args << " "), ...);}, tpl);
    // std::cout << std::endl;
    // auto constexpr bhv = lolita2::UnknownCollection(u1, u2);

    // auto constexpr fe1 =  lolita2::HybridDiscontinuousGalerkin(u1, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2), basis, basis);
    // auto constexpr fe2 =  lolita2::HybridDiscontinuousGalerkin(u2, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2), basis, basis);

    auto constexpr hdg = lolita2::HybridDiscontinuousGalerkin(basis, basis, lolita2::HybridDiscontinuousGalerkin::Hdg);

    auto constexpr fe1 =  lolita2::FiniteElementMethod(u1, bhv, hdg, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2));
    auto constexpr fe2 =  lolita2::FiniteElementMethod(u2, bhv, hdg, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2));
    
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    auto elements = lolita2::geometry::MeshFileParser(file_path).template makeFiniteElementSet<domain, fe1, fe2>();
    // std::cout << elements << std::endl;
    
    for (auto const & elem : elements->getElements<1, 0>())
    {
        std::cout << elem.second->getHash() << std::endl;
        // std::cout << elem.second->template getElement<0>()->template getOuterNeighbors<0, 0>().size() << std::endl;
        for (auto const & domain : elem.second->template getFiniteElement<0>()->domains_)
        {
            std::cout << domain->tag_ << std::endl;
        }
    }

    auto constexpr t_element = lolita2::geometry::Element::triangle(1);

    std::cout << "ici : " << lolita2::geometry::Traits<fe2>::template getNumElementUnknowns<t_element, domain>() << std::endl;
    std::cout << "ici : " << lolita2::geometry::Traits<fe2>::template getNumCellUnknowns<t_element, domain>() << std::endl;
    std::cout << "ici fld size : " << lolita2::geometry::FieldTraits<field>::template size<domain>() << std::endl;
    std::cout << "ici basis size : " << lolita2::geometry::FiniteElementBasisTraits<basis>::template size<t_element>() << std::endl;
    // std::cout << "ici bhv size : " << lolita2::geometry::BehaviorTraits<bhv>::template strainSize<domain>() << std::endl;
    

    // std::cout << elements.getElements<2, 0>().at("304")->getElement<0>()->getOuterNeighbors<0, 0>().size() << std::endl;
    // elements.getElements<2, 0>()["304"]->domains_;
    
}