#include "gtest/gtest.h"
#include "lolita_lolita/lolita_core/lolita_core_n_4004.hxx"

TEST(t_deb, t_deb)
{
    std::cout << std::fixed << std::setprecision(2);
    //
    auto tick = std::chrono::high_resolution_clock::now();
    auto tock = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration<double>(tock - tick);

    auto mat_ = lolita::Matrix<lolita::Real, 2, 2>();
    auto vec_ = lolita::Vector<lolita::Real, 2>();

    static_assert(lolita::RealVectorConcept<decltype(mat_ * lolita::Vector<lolita::Real, 2>::Zero())>);
    static_assert(lolita::RealMatrixConcept<decltype(mat_ * lolita::Vector<lolita::Real, 2>::Zero())>);

    // declaring behavior
    auto lib_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/bhv_elasticity/src/libBehaviour.so";
    auto lib_name = "Elasticity";
    auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
        mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
        mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    };
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    //
    //
    //
    // constants
    auto constexpr domain = lolita::Domain::cartesian(2);
    auto constexpr quadrature = lolita::Quadrature::gauss(2);
    auto constexpr cells = lolita::ElementType::cells(domain);
    auto constexpr faces = lolita::ElementType::faces(domain);
    // fields
    auto constexpr displacement_field = lolita::Field::vector();
    // generalized strains
    auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(lolita::Field::vector(), lolita::Mapping::smallStrain());
    // behaviors
    auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    // discretization
    auto constexpr hdg = lolita::HybridDiscontinuousGalerkin::hybridDiscontinuousGalerkin(1, 1);
    // finite elements
    auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    // mesh    
    // auto file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/quadrangle.msh";
    // auto file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/quadrangle004.msh";
    // auto file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/quadrangle050.msh";
    // auto file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/quadrangle_unstructured_00050.msh";
    // auto file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/triangle_structured_00050.msh";
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/square2.msh";
    // auto file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";
    // mesh build
    auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    // adding behavior
    auto micromorphic_damage = elements->setBehavior<cells, quadrature>("SQUARE", lib_path, lib_name, hyp);
    //making operators
    elements->setStrainOperators<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement");
    elements->setElementOperators<cells, displacement_element, hdg>("SQUARE", "Stabilization");
    // -------------------------------------------------------------------------------------------------------------------------------------------------------->
    auto constexpr t_quadrature = lolita::Quadrature::gauss(4);
    auto constexpr t_cell_basis = lolita::Basis::monomial(1);
    auto constexpr t_face_basis = lolita::Basis::monomial(1);
    auto constexpr t_grad_basis = lolita::Basis::monomial(1);
    auto constexpr t_face = lolita::Element::segment(1);
    // for (auto const & finite_element : {elements->template getElements<2, 0>()[0]})
    //
    //
    //
    auto printelem2 = [&] <int t_i = 0> (auto & self)
    {
        for (auto const & finite_element : elements->template getElements<2, t_i>())
        {
            finite_element->template getSymmetricGradientRhsDEBUG<displacement_field, hdg>(0, 0);
        }
        if constexpr (t_i < lolita::DomainTraits<domain>::template getNumElements<2>() - 1)
        {
            self.template operator()<t_i + 1>(self);
        }
    };
    printelem2(printelem2);
    //
    //
    //
    auto printelem = [&] <int t_i = 0> (auto & self)
    {
        // std::cout << "t_i ======================================== " << t_i << std::endl;
        // std::cout << "d ======================================== " << lolita::DomainTraits<domain>::template getNumElements<2>() << std::endl;
        auto row = 0;
        auto col = 0;
        auto celeem = 0;
        auto constexpr t_element = lolita::DomainTraits<domain>::template getElement<2, t_i>();
        for (auto const & finite_element : elements->template getElements<2, t_i>())
        {
            auto constexpr cell_cols_ = lolita::ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize();
            auto constexpr cell_rows_ = lolita::FiniteElementBasisTraits<t_cell_basis>::template getSize<t_element>();
            auto constexpr grad_rows_ = lolita::FiniteElementBasisTraits<t_grad_basis>::template getSize<t_element>();
            // auto mat = lolita::Matrix<lolita::Real, rows_, cols_>();
            auto cell_row_vectors = lolita::Matrix<lolita::Real, grad_rows_, cell_cols_>();
            auto cell_colI_vectors = lolita::Matrix<lolita::Real, cell_rows_, cell_cols_>();
            auto cell_colJ_vectors = lolita::Matrix<lolita::Real, cell_rows_, cell_cols_>();
            // auto matips2 = lolita::Matrix<lolita::Real, 3, cols_>();
            auto wts = lolita::Matrix<lolita::Real, 1, cell_cols_>();
            for (auto ipc = 0; ipc < lolita::ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); ipc++)
            {
                auto weight = finite_element->template getCurrentQuadratureWeight<t_quadrature>(ipc);
                auto point = finite_element->template getReferenceQuadraturePoint<t_quadrature>(ipc);
                auto row_cell_vector = finite_element->template getBasisEvaluation<t_grad_basis>(point);
                auto col_cell_vector_i = finite_element->template getBasisDerivative<t_cell_basis>(point, col);
                auto col_cell_vector_j = finite_element->template getBasisDerivative<t_cell_basis>(point, row);
                cell_row_vectors.col(ipc) = row_cell_vector;
                cell_colI_vectors.col(ipc) = col_cell_vector_i;
                cell_colJ_vectors.col(ipc) = col_cell_vector_j;
                wts(ipc) = weight;
            }
            auto current_coordinates = finite_element->getCurrentCoordinates();
            std::cout << "*********************** finite_element :" << std::endl;
            std::cout << current_coordinates << std::endl;
            // std::cout << "* cell_row_vectors :" << std::endl;
            // std::cout << cell_row_vectors << std::endl;
            // std::cout << "* cell_colI_vectors :" << std::endl;
            // std::cout << cell_colI_vectors << std::endl;
            // std::cout << "* cell_colJ_vectors :" << std::endl;
            // std::cout << cell_colJ_vectors << std::endl;
            // std::cout << "* wts :" << std::endl;
            // std::cout << wts << std::endl;
            // if (celeem == 1)
            // {
            //     finite_element->template getSymmetricGradientRhs<displacement_field, hdg>(0, 0);
            // }
            std::cout << "* op (0, 0) :" << std::endl;
            std::cout << finite_element->template getSymmetricGradientRhs<displacement_field, hdg>(0, 0) << std::endl;
            std::cout << "* op (0, 1) :" << std::endl;
            std::cout << finite_element->template getSymmetricGradientRhs<displacement_field, hdg>(0, 1) << std::endl;
            std::cout << "* op (1, 0) :" << std::endl;
            std::cout << finite_element->template getSymmetricGradientRhs<displacement_field, hdg>(1, 0) << std::endl;
            std::cout << "* op (1, 1) :" << std::endl;
            std::cout << finite_element->template getSymmetricGradientRhs<displacement_field, hdg>(1, 1) << std::endl;
            auto ipc = 0;
            for (auto & ip : finite_element->quadrature_.at("Elasticity").ips_)
            {
                std::cout << "* _x_q_c " << ipc << " :" << std::endl;
                std::cout << ip.coordinates_ << std::endl;
                std::cout << "* B " << ipc << " :" << std::endl;
                std::cout << ip.ops_.at("Displacement") << std::endl;
                ipc ++;
            }
            auto if_ = 0;
            for (auto const & face : finite_element->template getInnerNeighbors<0, 0>())
            {
                std::cout << "*** face :" << std::endl;
                auto face_current_coordinates = face->getCurrentCoordinates();
                std::cout << face_current_coordinates << std::endl;
                // auto unknowns = finite_element->degrees_of_freedom_.at("Displacement").template getCoefficients<lolita::Field::vector(), hdg.getFaceBasis()>();
                // std::cout << "* unknowns :" << std::endl;
                // std::cout << unknowns << std::endl;
                std::cout << "* sign :" << std::endl;
                std::cout << finite_element->template getInnerNeighborOrientation<0, 0>(if_) << std::endl;
                std::cout << "* normal_vector :" << std::endl;
                std::cout << face->getNormalVector(face->getReferenceCentroid()) << std::endl;
                auto constexpr cols_ = lolita::ElementQuadratureRuleTraits<t_face, t_quadrature>::getSize();
                auto constexpr rows_ = lolita::FiniteElementBasisTraits<t_face_basis>::template getSize<t_face>();
                auto mat = lolita::Matrix<lolita::Real, rows_, cols_>();
                auto matips = lolita::Matrix<lolita::Real, 3, cols_>();
                auto matips2 = lolita::Matrix<lolita::Real, 3, cols_>();
                auto wts = lolita::Matrix<lolita::Real, 1, cols_>();
                for (auto ipf = 0; ipf < lolita::ElementQuadratureRuleTraits<t_face, t_quadrature>::getSize(); ipf++)
                {
                    auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(ipf);
                    auto point = face->template getReferenceQuadraturePoint<t_quadrature>(ipf);
                    auto i_p = finite_element->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, 0>(if_, ipf);
                    // auto mp_0 = lolita::Point();
                    auto mp_1 = lolita::Point();
                    // mp_0.setZero();
                    mp_1.setZero();
                    for (auto ihelo = 0; ihelo < domain.getDim(); ++ihelo)
                    {
                        // mp_0(ihelo) = finite_element->getShapeMappingEvaluation(current_coordinates.row(ihelo), finite_element->getReferenceCentroid());
                        mp_1(ihelo) = finite_element->getShapeMappingEvaluation(current_coordinates.row(ihelo), i_p);
                    }
                    matips2.col(ipf) = mp_1;
                    auto row_cell_vector = finite_element->template getBasisEvaluation<t_grad_basis>(i_p);
                    auto col_face_vector = face->template getBasisEvaluation<t_face_basis>(point);
                    auto col_cell_vector = finite_element->template getBasisEvaluation<t_cell_basis>(i_p);
                    mat.col(ipf) = col_face_vector;
                    matips.col(ipf) = i_p;
                    wts(ipf) = weight;
                }
                std::cout << "* col_face_vector :" << std::endl;
                std::cout << mat << std::endl;
                std::cout << "* wts :" << std::endl;
                std::cout << wts << std::endl;
                std::cout << "* ips :" << std::endl;
                std::cout << matips2 << std::endl;
                if_ ++;
            }
            celeem ++;
        }
        if constexpr (t_i < lolita::DomainTraits<domain>::template getNumElements<2>() - 1)
        {
            self.template operator()<t_i + 1>(self);
        }
    };
    // printelem(printelem);
    //
    //
    //
    // <--------------------------------------------------------------------------------------------------------------------------------------------------------
    // dofs
    auto face_displacement = elements->setDegreeOfFreedom<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
    auto cell_displacement = elements->setDegreeOfFreedom<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
    //
    auto top_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
    auto left_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
    auto bottom_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
    // auto cell_rhs = lolita::Vector<lolita::Real>(cell_displacement->getCoefficients().size());
    // systems
    auto displacement_system = lolita::System::make();
    displacement_system->setUnknown("Displacement", face_displacement->size());
    displacement_system->setBinding("TopForce", top_force->size());
    displacement_system->setBinding("LeftForce", left_force->size());
    displacement_system->setBinding("BottomForce", bottom_force->size());
    displacement_system->initialize();
    std::cout << "displacement_system u size : " << displacement_system->getUnknownsSize() << std::endl;
    std::cout << "displacement_system b size : " << displacement_system->getBindingsSize() << std::endl;
    // load
    auto load0 = elements->setConstraint<faces>("TOP", "Pull", [](lolita::Point const & p, lolita::Real const & t) { return 1.0; }, 1, 0);
    auto load1 = elements->setConstraint<faces>("LEFT", "FixedL", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
    auto load2 = elements->setConstraint<faces>("BOTTOM", "FixedB", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 1, 0);
    // setting variable
    // elements->setMaterialProperty<cells>("SQUARE", "Elasticity", "YoungModulus", [](lolita::Point const & p) { return 1.0; });
    // elements->setMaterialProperty<cells>("SQUARE", "Elasticity", "PoissonRatio", [](lolita::Point const & p) { return 0.0; });
    elements->setExternalVariable<cells>("SQUARE", "Elasticity", "Temperature", [](lolita::Point const & p) { return 293.15; });
    // elements->setExternalVariable<cells>("SQUARE", "Elasticity", "Damage", [](lolita::Point const & p) { return 0.0; });
    // setting parameter
    elements->setParameter<faces>("TOP", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    elements->setParameter<faces>("LEFT", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    elements->setParameter<faces>("BOTTOM", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    // stab
    elements->setParameter<cells>("SQUARE", "Stabilization", [](lolita::Point const & p) { return 1.0; });
    //
    //
    elements->setStrainValues<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement");
    //
    elements->integrate<cells>("SQUARE", "Elasticity");
    //
    elements->assembleUnknownBlock<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement", displacement_system);
    //
    elements->assembleBindingBlock<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system);
    elements->assembleBindingBlock<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system);
    elements->assembleBindingBlock<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system);
    //    
    displacement_system->setCorrection();
    //
    elements->updateUnknown<cells, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
    elements->updateUnknown<faces, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce", displacement_system);
    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce", displacement_system);
    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce", displacement_system);

    // for (auto const & finite_element : elements->template getElements<1, 0>())
    // {
    //     std::cout << "*** finite_element :" << std::endl;
    //     std::cout << finite_element->getCurrentCoordinates() << std::endl;
    //     // auto unknowns = finite_element->degrees_of_freedom_.at("Displacement").template getCoefficients<lolita::Field::vector(), hdg.getFaceBasis()>();
    //     // std::cout << "* unknowns :" << std::endl;
    //     // std::cout << unknowns << std::endl;
    //     std::cout << "* sign :" << std::endl;
    //     std::cout << finite_element->getNormalVector(finite_element->getReferenceCentroid()) << std::endl;
    //     std::cout << "* normal_vector :" << std::endl;
    //     std::cout << finite_element->getNormalVector(finite_element->getReferenceCentroid()) << std::endl;
    // }

    //
    auto out_file = "/home/dsiedel/projetcs/lolita/lolita/tests/out1.msh";
    lolita::GmshFileParser::setOutput<domain>(out_file, elements, "Elasticity");
    // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, 0, 0.0, "Elasticity", 0);
    // lolita::GmshFileParser::addNodalDofOutput<2, domain, displacement_element, hdg>(out_file, elements, 0, 0.0, "Displacement", 0, 0);
    // lolita::GmshFileParser::addNodalDofOutput<2, domain, displacement_element, hdg>(out_file, elements, 0, 0.0, "Displacement", 1, 0);
    lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, 0, 0.0, "Displacement", "Elasticity", 0, 0);
    lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, 0, 0.0, "Displacement", "Elasticity", 1, 0);

    
    
}