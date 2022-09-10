#include "gtest/gtest.h"
#include "core/008_mesh.hxx"

TEST(test_quadrangle, test_quadrangle)
{
    std::cout << std::fixed << std::setprecision(3);
    //
    auto constexpr domain = lolita::Domain::cartesian(2);
    auto constexpr node = lolita::Element::node();
    auto constexpr segment = lolita::Element::segment(1);
    auto constexpr quadrangle = lolita::Element::quadrangle(1);
    //
    auto constexpr hdg = lolita::HybridDiscontinuousGalerkin::hybridDiscontinuousGalerkin(1, 1);
    auto constexpr quadrature = lolita::Quadrature::gauss(4);
    auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(lolita::Field::vector(), lolita::Mapping::smallStrain());
    auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    //
    auto node_0 = std::make_shared<lolita::FiniteElementHolder<node, domain>>();
    auto node_1 = std::make_shared<lolita::FiniteElementHolder<node, domain>>();
    auto node_2 = std::make_shared<lolita::FiniteElementHolder<node, domain>>();
    auto node_3 = std::make_shared<lolita::FiniteElementHolder<node, domain>>();
    //
    auto segment_0 = std::make_shared<lolita::FiniteElementHolder<segment, domain>>();
    auto segment_1 = std::make_shared<lolita::FiniteElementHolder<segment, domain>>();
    auto segment_2 = std::make_shared<lolita::FiniteElementHolder<segment, domain>>();
    auto segment_3 = std::make_shared<lolita::FiniteElementHolder<segment, domain>>();
    //
    auto quadrangle_0 = std::make_shared<lolita::FiniteElementHolder<quadrangle, domain>>();
    // node 0
    node_0->coordinates_ = std::make_shared<lolita::Point>(lolita::Point({0.0, 0.0, 0.0}));
    // node 1
    node_1->coordinates_ = std::make_shared<lolita::Point>(lolita::Point({1.0, 0.0, 0.0}));
    // node 2
    node_2->coordinates_ = std::make_shared<lolita::Point>(lolita::Point({1.0, 1.0, 0.0}));
    // node 3
    node_3->coordinates_ = std::make_shared<lolita::Point>(lolita::Point({0.0, 1.0, 0.0}));
    // segment 0 / nodes
    segment_0->getInnerNeighbors<0, 0>()[0] = node_0;
    node_0->getOuterNeighbors<0, 0>().push_back(segment_0);
    segment_0->getInnerNeighbors<0, 0>()[1] = node_1;
    node_1->getOuterNeighbors<0, 0>().push_back(segment_0);
    // segment 1 / nodes
    segment_1->getInnerNeighbors<0, 0>()[0] = node_1;
    node_1->getOuterNeighbors<0, 0>().push_back(segment_1);
    segment_1->getInnerNeighbors<0, 0>()[1] = node_2;
    node_2->getOuterNeighbors<0, 0>().push_back(segment_1);
    // segment 2 / nodes
    segment_2->getInnerNeighbors<0, 0>()[0] = node_2;
    node_2->getOuterNeighbors<0, 0>().push_back(segment_2);
    segment_2->getInnerNeighbors<0, 0>()[1] = node_3;
    node_3->getOuterNeighbors<0, 0>().push_back(segment_2);
    // segment 3 / nodes
    segment_3->getInnerNeighbors<0, 0>()[0] = node_3;
    node_3->getOuterNeighbors<0, 0>().push_back(segment_3);
    segment_3->getInnerNeighbors<0, 0>()[1] = node_0;
    node_0->getOuterNeighbors<0, 0>().push_back(segment_3);
    // quadrangle 0 / segments
    quadrangle_0->getInnerNeighbors<0, 0>()[0] = segment_0;
    segment_0->getOuterNeighbors<1, 1>().push_back(quadrangle_0);
    quadrangle_0->getInnerNeighbors<0, 0>()[1] = segment_1;
    segment_1->getOuterNeighbors<1, 1>().push_back(quadrangle_0);
    quadrangle_0->getInnerNeighbors<0, 0>()[2] = segment_2;
    segment_2->getOuterNeighbors<1, 1>().push_back(quadrangle_0);
    quadrangle_0->getInnerNeighbors<0, 0>()[3] = segment_3;
    segment_3->getOuterNeighbors<1, 1>().push_back(quadrangle_0);
    // quadrangle 0 / nodes
    quadrangle_0->getInnerNeighbors<1, 0>()[0] = node_0;
    node_0->getOuterNeighbors<1, 1>().push_back(quadrangle_0);
    quadrangle_0->getInnerNeighbors<1, 0>()[1] = node_1;
    node_1->getOuterNeighbors<1, 1>().push_back(quadrangle_0);
    quadrangle_0->getInnerNeighbors<1, 0>()[2] = node_2;
    node_2->getOuterNeighbors<1, 1>().push_back(quadrangle_0);
    quadrangle_0->getInnerNeighbors<1, 0>()[3] = node_3;
    node_3->getOuterNeighbors<1, 1>().push_back(quadrangle_0);

    // quadrangle_0->template getSymmetricGradientRhs<lolita::Field::vector(), hdg>(0, 1);

    // auto count_seg = 0;
    // for (auto const & seg : {segment_0, segment_1, segment_2})
    // {
    //     auto mass = lolita::Matrix<lolita::Real, 2, 2>();
    //     mass.setZero();
    //     std::cout << "--- seg " << count_seg << "\n\n";
    //     for (auto q = 0; q < lolita::ElementQuadratureRuleTraits<segment, quadrature>::getSize(); ++q)
    //     {
    //         auto point = seg->template getReferenceQuadraturePoint<quadrature>(q);
    //         auto c_point = seg->template getCurrentQuadraturePoint<quadrature>(q);
    //         auto weight = seg->template getCurrentQuadratureWeight<quadrature>(q);
    //         auto dist = seg->getLocalFrameDistance(seg->getReferenceCentroid(), point, 0);
    //         auto diams = seg->getLocalFrameDiameters();
    //         auto n = seg->getNormalVector(point);
    //         auto rot = seg->getRotationMatrix(seg->getReferenceCentroid());
    //         auto s_q_f = rot * c_point;
    //         auto s_f = rot * seg->getCurrentCentroid();
    //         std::cout << "weight " << q <<" is :\n" << weight << "\n";
    //         std::cout << "point " << q <<" is :\n" << c_point << "\n";
    //         std::cout << "n " << q <<" is :\n" << n << "\n";
    //         std::cout << "s_q_f " << q <<" is :\n" << s_q_f << "\n";
    //         std::cout << "s_f " << q <<" is :\n" << s_f << "\n";
    //         std::cout << "dist " << q <<" is :\n" << dist << "\n";
    //         std::cout << "diams " << q <<" is :\n" << diams << "\n";
    //         std::cout << "rot " << q <<" is :\n" << rot << "\n";
    //         auto basis_vector = seg->template getBasisEvaluation<lolita::Basis::monomial(1)>(point);
    //         mass += weight * basis_vector * basis_vector.transpose();
    //     }
    //     std::cout << "mass " << "\n";
    //     std::cout << mass << "\n\n";
    //     count_seg ++;
    // }

    for (auto q = 0; q < lolita::ElementQuadratureRuleTraits<quadrangle, lolita::Quadrature::gauss(2)>::getSize(); ++q)
    {
        auto point = quadrangle_0->template getCurrentQuadraturePoint<lolita::Quadrature::gauss(2)>(q);
        auto r_point = quadrangle_0->template getReferenceQuadraturePoint<lolita::Quadrature::gauss(2)>(q);
        auto weight = quadrangle_0->template getCurrentQuadratureWeight<lolita::Quadrature::gauss(2)>(q);
        auto basis_vector = quadrangle_0->template getBasisEvaluation<lolita::Basis::monomial(1)>(point);
        // std::cout << "basis_vector" << "\n";
        // std::cout << basis_vector << "\n";
        auto grad = quadrangle_0->template getMapping<lolita::Field::vector(), lolita::Mapping::smallStrain(), hdg>(r_point);
        std::cout << "grad ------------------ " << q << "\n";
        std::cout << lolita::mat2str(grad, 3) << "\n";
        std::cout << "------------------ " << q << "\n";
    }
    auto stab = quadrangle_0->template getStabilization<lolita::Field::vector(), hdg>();
    std::cout << "stab ------------------ " << "\n";
    std::cout << lolita::mat2str(stab, 3) << "\n";
    std::cout << "------------------ " << "\n";

    std::cout << lolita::mat2str(stab.col(1), 3) << "\n";

    // auto count_seg = 0;
    // for (auto const & seg : {segment_0, segment_1, segment_2})
    // {
    //     auto mass = lolita::Matrix<lolita::Real, 2, 2>();
    //     mass.setZero();
    //     std::cout << "--- seg " << count_seg << "\n\n";
    //     for (auto q = 0; q < lolita::ElementQuadratureRuleTraits<segment, lolita::Quadrature::gauss(4)>::getSize(); ++q)
    //     {
    //         auto point = seg->template getReferenceQuadraturePoint<lolita::Quadrature::gauss(4)>(q);
    //         auto i_point = quadrangle_0->getInnerNeighborReferenceQuadraturePoint<lolita::Quadrature::gauss(4), 0, 0>(count_seg, q);
    //         auto c_point = seg->template getCurrentQuadraturePoint<lolita::Quadrature::gauss(4)>(q);
    //         auto weight = seg->template getCurrentQuadratureWeight<lolita::Quadrature::gauss(4)>(q);
    //         auto basis_vector = quadrangle_0->template getBasisEvaluation<lolita::Basis::monomial(1)>(i_point);
    //         std::cout << "basis_vector " << q <<" is :\n" << basis_vector << "\n";
    //         // auto dist = seg->getLocalFrameDistance(seg->getReferenceCentroid(), point, 0);
    //         // auto diams = seg->getLocalFrameDiameters();
    //         // auto n = seg->getNormalVector(point);
    //         // auto rot = seg->getRotationMatrix(seg->getReferenceCentroid());
    //         // auto s_q_f = rot * c_point;
    //         // auto s_f = rot * seg->getCurrentCentroid();
    //         // std::cout << "weight " << q <<" is :\n" << weight << "\n";
    //         // std::cout << "point " << q <<" is :\n" << c_point << "\n";
    //         // std::cout << "n " << q <<" is :\n" << n << "\n";
    //         // std::cout << "s_q_f " << q <<" is :\n" << s_q_f << "\n";
    //         // std::cout << "s_f " << q <<" is :\n" << s_f << "\n";
    //         // std::cout << "dist " << q <<" is :\n" << dist << "\n";
    //         // std::cout << "diams " << q <<" is :\n" << diams << "\n";
    //         // std::cout << "rot " << q <<" is :\n" << rot << "\n";
    //         auto basis_vector2 = seg->template getBasisEvaluation<lolita::Basis::monomial(1)>(point);
    //         std::cout << "basis_vector2 " << q <<" is :\n" << basis_vector2 << "\n";
    //         mass += weight * basis_vector2 * basis_vector2.transpose();
    //     }
    //     // std::cout << "mass " << "\n";
    //     // std::cout << mass << "\n\n";
    //     count_seg ++;
    // }

}