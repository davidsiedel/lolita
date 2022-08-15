#include "gtest/gtest.h"

TEST(tm, tm)
{
    // constants
    // auto constexpr domain = lolita::Domain::cartesian(2);
    // auto constexpr cell_basis = lolita::Basis::monomial(1);
    // auto constexpr face_basis = lolita::Basis::monomial(1);
    // auto constexpr quadrature = lolita::Quadrature(lolita::Quadrature::Gauss, 2);
    // // fields
    // auto constexpr displacement_field = lolita::Field::vector("Displacement");
    // auto constexpr damage_field = lolita::Field::scalar("Damage");
    // // generalized strains
    // auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(displacement_field, lolita::Mapping::smallStrain());
    // auto constexpr damage_generalized_strain = lolita::GeneralizedStrain(damage_field, lolita::Mapping::gradient(), lolita::Mapping::identity());
    // // behaviors
    // auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    // auto constexpr damage_behavior = lolita::Behavior(damage_generalized_strain);
    // // discretization
    // auto constexpr hdg = lolita::HybridDiscontinuousGalerkin(cell_basis, face_basis, lolita::HybridDiscontinuousGalerkin::Stabilization::Hdg);
    // // finite elements
    // auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    // auto constexpr damage_element =  lolita::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);

    // auto constexpr t_element = lolita::core::Element::triangle(1);

    // // std::cout << "ici : " << lolita::core::FiniteElementMethodTraits<damage_element>::template getNumElementUnknowns<t_element, domain>() << std::endl;
    // // std::cout << "ici : " << lolita::core::FiniteElementMethodTraits<damage_element>::template getNumCellUnknowns<t_element, domain>() << std::endl;
    // std::cout << "ici fld size : " << lolita::core::FieldTraits<displacement_field>::template size<domain>() << std::endl;
    // std::cout << "ici basis size : " << lolita::core::FiniteElementBasisTraits<cell_basis>::template size<t_element>() << std::endl;
    // // std::cout << "ici displacement_behavior size : " << lolita::core::BehaviorTraits<displacement_behavior>::template strainSize<domain>() << std::endl;

    // std::cout << "yo : " << std::endl;
    // std::cout << lolita::core::FiniteElementHolder<t_element, domain, displacement_element, damage_element>::template getArgIndex<displacement_element>();
    // std::cout << " , ";
    // std::cout << lolita::core::FiniteElementHolder<t_element, domain, displacement_element, damage_element>::template getArgIndex<damage_element>();
    // std::cout << std::endl;
    
    // std::cout << "yo : " << std::endl;
    // std::cout << lolita2::geometry::FiniteElementHolder<t_element, domain, displacement_element, damage_element>::template getArgIndex<displacement_element>();
    // std::cout << " , ";
    // std::cout << lolita2::geometry::FiniteElementHolder<t_element, domain, displacement_element, damage_element>::template getArgIndex<damage_element>();
    // // std::cout << "ici : " << lolita::core::FiniteElementMethodTraits<damage_element>::template getNumElementUnknowns<t_element, domain>() << std::endl;
    // // std::cout << "ici : " << lolita::core::FiniteElementMethodTraits<damage_element>::template getNumCellUnknowns<t_element, domain>() << std::endl;
    // std::cout << "ici fld size : " << lolita::core::FieldTraits<displacement_field>::template size<domain>() << std::endl;
    // std::cout << "ici basis size : " << lolita::core::FiniteElementBasisTraits<cell_basis>::template size<t_element>() << std::endl;
    // // std::cout << "ici displacement_behavior size : " << lolita::core::BehaviorTraits<displacement_behavior>::template strainSize<domain>() << std::endl;

    // // std::cout << elements.getElements<2, 0>().at("304")->getElement<0>()->getOuterNeighbors<0, 0>().size() << std::endl;
    // // elements.getElements<2, 0>()["304"]->domains_;

    // auto constexpr size = 20;

    // auto mat = lolita::matrix::Matrix<Real>(size, size);
    // auto vec = lolita::matrix::Vector<Real>(size);
    // auto mat_view = lolita::matrix::Span<lolita::matrix::Matrix<Real, size, size> const>(mat.data());
    // auto vec_view = lolita::matrix::Span<lolita::matrix::Vector<Real, size> const>(vec.data());
    // auto tick = std::chrono::high_resolution_clock::now();
    // auto res = mat_view * vec_view;
    // auto tock = std::chrono::high_resolution_clock::now();
    // auto time = std::chrono::duration<double>(tock - tick);
    // std::cout << "time : " << time.count() << std::endl;

    // // auto mat2 = std::make_unique<lolita::matrix::Matrix<Real, size, size>>(lolita::matrix::Matrix<Real, size, size>());
    // // auto vec2 = std::make_unique<lolita::matrix::Vector<Real, size>>(lolita::matrix::Vector<Real, size>());
    // // auto tick2 = std::chrono::high_resolution_clock::now();
    // // auto res2 = (* mat2) * (* vec2);
    // // auto tock2 = std::chrono::high_resolution_clock::now();
    // // auto time2 = std::chrono::duration<double>(tock2 - tick2);
    // // std::cout << "time : " << time2.count() << std::endl;

    // auto mat3 = lolita::matrix::Matrix<Real>(size, size);
    // auto vec3 = lolita::matrix::Vector<Real>(size);
    // auto tick3 = std::chrono::high_resolution_clock::now();
    // auto res3 = mat3 * vec3;
    // auto tock3 = std::chrono::high_resolution_clock::now();
    // auto time3 = std::chrono::duration<double>(tock3 - tick3);
    // std::cout << "time : " << time3.count() << std::endl;

    // std::cout << "size of : " << std::endl;
    // std::cout << "std::vector<std::shared_ptr<Loading>> " << std::endl;
    // std::cout << sizeof(std::vector<std::shared_ptr<lolita::Loading>>) << std::endl;
    // std::cout << "size of : " << std::endl;
    // std::cout << "std::vector<Real> " << std::endl;
    // std::cout << sizeof(std::vector<Real>) << std::endl;
    // std::cout << "std::map<std::string, std::shared_ptr<lolita::Loading>> " << std::endl;
    // std::cout << sizeof(std::map<lolita::matrix::Coordinates, std::shared_ptr<lolita::Loading>>) << std::endl;
    // std::cout << "std::array<std::shared_ptr<lolita::Loading>, 3, 1> " << std::endl;
    // std::cout << sizeof(std::array<std::array<std::shared_ptr<lolita::Loading>, 3>, 1>) << std::endl;

    // auto mapp = std::map<std::pair<Integer, Integer>, std::string>();
    // mapp[std::make_pair(0, 0)] = "HJK";

    // auto stdvec = std::vector<Integer>{1, 2, 3};
    // auto stdmap = std::map<std::string, Integer>{{"1", 1}, {"2", 2}, {"3", 3}};
    // // auto vec3 = lolita::matrix::Vector<Real>(size);
    // auto tickv = std::chrono::high_resolution_clock::now();
    // auto resv = stdvec[2];
    // auto tockv = std::chrono::high_resolution_clock::now();
    // auto timev = std::chrono::duration<double>(tockv - tickv);
    // std::cout << "time vector : " << timev.count() << std::endl;
    // auto tickm = std::chrono::high_resolution_clock::now();
    // auto resm = stdmap["3"];
    // auto tockm = std::chrono::high_resolution_clock::now();
    // auto timem = std::chrono::duration<double>(tockm - tickm);
    // std::cout << "time map : " << timem.count() << std::endl;
    // // std::cout << * stdvec.end() << std::endl;
    
}