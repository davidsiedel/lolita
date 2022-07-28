#include "gtest/gtest.h"
#include "lolita/lolita_core_n_111.hxx"

TEST(t0, t0)
{
    
    auto constexpr domain = lolita2::Domain::cartesian(2);
    auto constexpr displacement_field = lolita2::Field::vector("Displacement");
    auto constexpr damage_field = lolita2::Field::scalar("Damage");
    auto constexpr basis = lolita2::Basis::monomial(1);

    // auto lmp = lolita2::Discretization::HJKL(lolita2::Quadrature(lolita2::Quadrature::Gauss, 2), 2);


    // auto constexpr u0 = lolita2::GeneralizedStrain(lolita2::Field::vector("A"), lolita2::Mapping::smallStrain());
    // auto constexpr u00 = lolita2::GeneralizedStrain(lolita2::Field::vector("B"), lolita2::Mapping::smallStrain());
    auto constexpr u1 = lolita2::GeneralizedStrain(displacement_field, lolita2::Mapping::smallStrain());
    auto constexpr u2 = lolita2::GeneralizedStrain(damage_field, lolita2::Mapping::gradient(), lolita2::Mapping::identity());

    auto constexpr bhv = lolita2::Behavior(u1, u2);

    auto constexpr hdg = lolita2::HybridDiscontinuousGalerkin(basis, basis, lolita2::HybridDiscontinuousGalerkin::Hdg);

    auto constexpr fe1 =  lolita2::FiniteElementMethod(u1, bhv, hdg, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2));
    auto constexpr fe2 =  lolita2::FiniteElementMethod(u2, bhv, hdg, lolita2::Quadrature(lolita2::Quadrature::Gauss, 2));

    // auto constexpr strainSize = lolita2::geometry::BehaviorTraits<bhv>::template getGeneralizedStrainSize<domain>();
    auto constexpr strainSize = lolita2::geometry::FiniteElementMethodTraits<fe1>::template getGeneralizedStrainSize<domain>();
    std::cout << "strainSize : " << strainSize << std::endl;
    auto constexpr strainOffset = lolita2::geometry::FiniteElementMethodTraits<fe1>::template getGeneralizedStrainOffset<domain>();
    std::cout << "strainOffset : " << strainOffset << std::endl;
    
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    auto elements = lolita2::geometry::MeshFileParser(file_path).template makeFiniteElementSet<domain, fe1, fe2>();
    auto degree_of_freedom = std::make_shared<lolita2::geometry::DegreeOfFreedom>("FaceDisplacement");
    elements->activate<2, 0, 1>("SQUARE");
    elements->activate<1, 0, 1>("SQUARE");
    elements->setDegreeOfFreedom<displacement_field, basis, 1, 0>("SQUARE", degree_of_freedom);
    // elements->setDegreeOfFreedom<displacement_field, basis, 0, 0>("SQUARE", degree_of_freedom);
    // std::cout << elements << std::endl;
    
    // for (auto const & elem : elements->getElements<1, 0>())
    // {
    //     std::cout << elem.second->getHash() << std::endl;
    //     // std::cout << elem.second->template getElement<0>()->template getOuterNeighbors<0, 0>().size() << std::endl;
    //     for (auto const & domain : elem.second->template getFiniteElement<0>()->domains_)
    //     {
    //         std::cout << domain->tag_ << std::endl;
    //     }
    // }

    // auto degree_of_freedom = std::make_shared<lolita2::geometry::DegreeOfFreedom>(lolita2::geometry::DegreeOfFreedom("FaceDisplacement"));

    std::cout << degree_of_freedom->coefficients_.size() << std::endl;

    auto constexpr t_element = lolita2::geometry::Element::triangle(1);

    // std::cout << "ici : " << lolita2::geometry::FiniteElementMethodTraits<fe2>::template getNumElementUnknowns<t_element, domain>() << std::endl;
    // std::cout << "ici : " << lolita2::geometry::FiniteElementMethodTraits<fe2>::template getNumCellUnknowns<t_element, domain>() << std::endl;
    std::cout << "ici fld size : " << lolita2::geometry::FieldTraits<displacement_field>::template size<domain>() << std::endl;
    std::cout << "ici basis size : " << lolita2::geometry::FiniteElementBasisTraits<basis>::template size<t_element>() << std::endl;
    // std::cout << "ici bhv size : " << lolita2::geometry::BehaviorTraits<bhv>::template strainSize<domain>() << std::endl;
    

    // std::cout << elements.getElements<2, 0>().at("304")->getElement<0>()->getOuterNeighbors<0, 0>().size() << std::endl;
    // elements.getElements<2, 0>()["304"]->domains_;

    auto constexpr size = 500;

    auto mat = lolita::matrix::Matrix<lolita::real>(size, size);
    auto vec = lolita::matrix::Vector<lolita::real>(size);
    auto mat_view = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, size, size> const>(mat.data());
    auto vec_view = lolita::matrix::Span<lolita::matrix::Vector<lolita::real, size> const>(vec.data());
    auto tick = std::chrono::high_resolution_clock::now();
    auto res = mat_view * vec_view;
    auto tock = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration<double>(tock - tick);
    std::cout << "time : " << time.count() << std::endl;

    // auto mat2 = std::make_unique<lolita::matrix::Matrix<lolita::real, size, size>>(lolita::matrix::Matrix<lolita::real, size, size>());
    // auto vec2 = std::make_unique<lolita::matrix::Vector<lolita::real, size>>(lolita::matrix::Vector<lolita::real, size>());
    // auto tick2 = std::chrono::high_resolution_clock::now();
    // auto res2 = (* mat2) * (* vec2);
    // auto tock2 = std::chrono::high_resolution_clock::now();
    // auto time2 = std::chrono::duration<double>(tock2 - tick2);
    // std::cout << "time : " << time2.count() << std::endl;

    auto mat3 = lolita::matrix::Matrix<lolita::real>(size, size);
    auto vec3 = lolita::matrix::Vector<lolita::real>(size);
    auto tick3 = std::chrono::high_resolution_clock::now();
    auto res3 = mat3 * vec3;
    auto tock3 = std::chrono::high_resolution_clock::now();
    auto time3 = std::chrono::duration<double>(tock3 - tick3);
    std::cout << "time : " << time3.count() << std::endl;
    
}