#include "gtest/gtest.h"
#include "lolita/lolita_core_n_111.hxx"

TEST(t0, t0)
{
    
    auto constexpr domain = lolita2::Domain(2, lolita2::Domain::Cartesian);
    auto constexpr field = lolita2::Field(2);
    auto constexpr basis = lolita2::Basis(lolita2::Basis::Monomial, 2);

    // auto lmp = lolita2::Discretization::HJKL(lolita2::Quadrature(lolita2::Quadrature::Gauss, 2), 2);


    auto constexpr u1 = lolita2::Unknown(lolita2::Field(1), lolita2::Mapping(lolita2::Mapping::Gradient));
    auto constexpr u2 = lolita2::Unknown(lolita2::Field(2), lolita2::Mapping(lolita2::Mapping::Gradient), lolita2::Mapping(lolita2::Mapping::Identity));

    // auto constexpr fe1 =  lolita2::HybridDiscontinuousGalerkin(u1, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2), basis, basis);
    // auto constexpr fe2 =  lolita2::HybridDiscontinuousGalerkin(u2, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2), basis, basis);

    auto constexpr hdg = lolita2::discretization::HybridDiscontinuousGalerkin(basis, basis, lolita2::discretization::HybridDiscontinuousGalerkin::Hdg);

    auto constexpr fe1 =  lolita2::FiniteElementMethod(u1, hdg, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2));
    auto constexpr fe2 =  lolita2::FiniteElementMethod(u2, hdg, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2));
    
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

    // std::cout << elements.getElements<2, 0>().at("304")->getElement<0>()->getOuterNeighbors<0, 0>().size() << std::endl;
    // elements.getElements<2, 0>()["304"]->domains_;
    
}