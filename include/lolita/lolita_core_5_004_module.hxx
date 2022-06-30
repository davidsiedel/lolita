//
// Created by dsiedel on 14/06/22.
//

#ifndef LOLITA_LOLITA_CORE_5_004_MODULE_HXX
#define LOLITA_LOLITA_CORE_5_004_MODULE_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_4.hxx"
#include "lolita/lolita_core_5_000_connectivity.hxx"
#include "lolita/lolita_core_5_001_base.hxx"
#include "lolita/lolita_core_5_002_basis.hxx"
#include "lolita/lolita_core_5_003_unknown.hxx"

namespace lolita::core2::finite_element
{

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementBase2 : FiniteElementBasis<t_element, t_domain, t_finite_element>
    {};

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementCell_n : FiniteElementBase2<t_element, t_domain, t_finite_element>
    {

    private:

        /**
         * @brief
         */
        using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         */
//        using t_FiniteElementCellTraits = cell::FiniteElementCellTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         */
        using t_Quadrature = typename FiniteElementTraits<t_element, t_domain, t_finite_element>::Quadrature;

    };

//    namespace face
//    {
//
//        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//        struct FiniteElementFaceTraits;
//
//        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//        requires(t_finite_element.discretization_.isHHO())
//        struct FiniteElementFaceTraits<t_element, t_domain, t_finite_element>
//        {
//
//            struct Face
//            {};
//
//            struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
//            {
//
//                template<auto t_element_group>
//                void
//                makeDirichlet(
//                        lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
//                )
//                {
//
//                }
//
//            };
//
//        };
//
//    }
//
//    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//    struct FiniteElementFace
//    {
//
//    private:
//
//        /**
//         * @brief
//         */
//        using t_FiniteElementTraits = lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>;
//
//    public:
//
//        /**
//         * @brief
//         */
//        lolita::integer static constexpr dim_structural_unknowns = t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Structural()>();
//
//        /**
//         * @brief
//         */
//        lolita::integer static constexpr dim_subsidiary_unknowns = t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Subsidiary()>();
//
//        /**
//         * @brief
//         */
//        lolita::integer static constexpr num_structural_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Structural()>();
//
//        /**
//         * @brief
//         */
//        lolita::integer static constexpr num_subsidiary_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
//
//        /**
//         * @brief
//         */
//        struct IntegrationPoint
//        {
//
//            /**
//             * @brief
//             */
//            lolita::domain::Point coordinates_;
//
//            /**
//             * @brief
//             */
//            lolita::real weight_;
//
//            /**
//             * @brief
//             */
//            std::unique_ptr<mgis::behaviour::BehaviourData> material_point_;
//
//            /**
//             * @brief
//             */
//            lolita::matrix::Vector<lolita::real, dim_structural_unknowns> structural_load_operator_;
//
//            /**
//             * @brief
//             */
//            lolita::matrix::Vector<lolita::real, dim_subsidiary_unknowns> subsidiary_load_operator_;
//
//        };
//
//        /**
//         * @brief
//         */
//        std::array<IntegrationPoint, t_FiniteElementTraits::Quadrature::dim_> integration_points_;
//
//    };
//
//    /**
//     * @brief
//     * @tparam t_element
//     * @tparam t_domain
//     * @tparam t_finite_element
//     */
//    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//    struct FiniteElementCell;
//
//    namespace cell
//    {
//
//        /**
//         * @brief
//         * @tparam t_element
//         * @tparam t_domain
//         * @tparam t_finite_element
//         */
//        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//        struct FiniteElementCellTraits;
//
//        /**
//         * @brief
//         * @tparam t_element
//         * @tparam t_domain
//         * @tparam t_finite_element
//         */
//        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//        requires(t_finite_element.discretization_.isHHO())
//        struct FiniteElementCellTraits<t_element, t_domain, t_finite_element>
//        {
//
//        private:
//
//            /**
//             * @brief
//             */
//            using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element>;
//
//            /**
//             * @brief
//             */
//            using t_FiniteElementCell = FiniteElementCell<t_element, t_domain, t_finite_element>;
//
//            /**
//             * @brief
//             * @return
//             */
//            static constexpr
//            lolita::integer
//            getNumUnknowns()
//            {
//                return t_FiniteElementTraits::getNumUnknowns();
//            }
//
//            /**
//             * @brief
//             * @return
//             */
//            static constexpr
//            lolita::integer
//            getNumCellUnknowns()
//            {
//                return t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
//            }
//
//            /**
//             * @brief
//             * @return
//             */
//            static constexpr
//            lolita::integer
//            getNumFaceUnknowns()
//            {
//                return t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Structural()>();
//            }
//
//        public:
//
//            /**
//             * @brief
//             * @return
//             */
//            static constexpr
//            lolita::integer
//            getLoadVectorSize()
//            {
//                return t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Subsidiary()>();
//            }
//
//            /**
//             * @brief
//             * @return
//             */
//            static constexpr
//            lolita::integer
//            getUnknownVectorSize()
//            {
//                auto num_cell_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
//                auto num_face_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Structural()>();
//                return t_FiniteElementTraits::getNumUnknowns();
//            }
//
//            struct Cell
//            {
//
//                using Stabilization = lolita::matrix::Matrix<lolita::real, getUnknownVectorSize(), getUnknownVectorSize()>;
//
//                using KTTinv = lolita::matrix::Matrix<lolita::real, getNumCellUnknowns(), getNumCellUnknowns()>;
//
//                using KFT = lolita::matrix::Matrix<lolita::real, getNumFaceUnknowns(), getNumCellUnknowns()>;
//
//                Stabilization stabilization_;
//
//                KTTinv ktt_inv_;
//
//                KFT kft_;
//
//            };
//
//            struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
//            {
//
////                lolita::matrix::Vector<lolita::real, getUnknownVectorSize()>
////                getUnknowns()
////                const
////                {
////                    using t_UnknownVector = lolita::matrix::Vector<lolita::real, getUnknownVectorSize()>;
////                    auto unknowns = lolita::matrix::Zero<t_UnknownVector>();
////                    auto count = 0;
////                    auto set_cell_block = [&] () {
////                        using t_CellUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>::FieldUnknowns ;
////                        using t_CellUnknown = std::tuple_element_t<0, typename t_CellUnknowns::Unknowns>;
////                        for (int i = 0; i < t_FiniteElementTraits::Field::shape_.rows_; ++i) {
////                            for (int j = 0; j < t_FiniteElementTraits::Field::shape_.cols_; ++j) {
////                                auto const & rhs = this->unknowns_.template getUnknown<0>().getDegreeOfFreedom(i, j).getUnknownsCoefficients();
////                                unknowns.template segment<t_CellUnknown::dim_>(count) = rhs;
////                                count += t_CellUnknown::dim_;
////                            }
////                        }
////                    };
////                    auto set_faces_block = [&] <lolita::integer t_i = 0> (auto & self) mutable {
////                        auto const constexpr _face = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<0, t_i>();
////                        using t_FaceUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<_face, t_domain, t_finite_element>::FieldUnknowns;
////                        using t_FaceUnknown = std::tuple_element_t<0, typename t_FaceUnknowns::Unknowns>;
////                        auto const & faces = this->template getComponents<0, t_i>();
////                        for (auto const & face : faces) {
////                            for (int i = 0; i < t_FiniteElementTraits::Field::shape_.rows_; ++i) {
////                                for (int j = 0; j < t_FiniteElementTraits::Field::shape_.cols_; ++j) {
////                                    auto const & rhs = face->unknowns_.template getUnknown<0>().getDegreeOfFreedom(i, j).getUnknownsCoefficients();
////                                    unknowns.template segment<t_FaceUnknown::dim_>(count) = rhs;
////                                    count += t_FaceUnknown::dim_;
////                                }
////                            }
////                        }
////                        if constexpr (t_i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<0>() - 1) {
////                            self.template operator()<t_i + 1>(self);
////                        }
////                    };
////                    set_cell_block();
////                    set_faces_block(set_faces_block);
////                    return unknowns;
////                }
//
//                void
//                setGeneralizedGradients()
//                {
//                    using t_Strain = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows()>;
////                    auto unknowns = this->template getUnknowns<unknown::Unknown::Structural()>();
//                    auto structural_unknowns = this->template getUnknowns<unknown::Unknown::Structural()>();
//                    auto subsidiary_unknowns = this->template getUnknowns<unknown::Unknown::Subsidiary()>();
//                    for (int i = 0; i < t_FiniteElementTraits::Quadrature::dim_; ++i) {
//                        /*
//                         * update the generalized gradient vector for the current integration point
//                         */
////                        auto generalized_gradient_values = t_Strain(this->integration_points_[i].generalized_strain_operator_ * unknowns);
//                        auto structural_strain_values = t_Strain(this->integration_points_[i].structural_strain_operator_ * structural_unknowns);
//                        auto subsidiary_strain_values = t_Strain(this->integration_points_[i].subsidiary_strain_operator_ * subsidiary_unknowns);
//                        auto strain = t_Strain(structural_strain_values + subsidiary_strain_values);
////                        auto generalized_gradient_values = _Gradients::Zero();
////                        generalized_gradient_values.setZero();
//                        auto set_generalized_gradient = [&] <lolita::integer t_i = 0> (auto & self) mutable {
//                            auto constexpr t_mapping = t_finite_element.unknown_.mappings_[t_i];
//                            using t_MappingPolicy = lolita::core2::field::MappingPolicy<t_finite_element.unknown_.tensor_, t_domain, t_mapping>;
//                            auto constexpr t_mapping_vector_block = t_FiniteElementTraits::template getMappingBlock<t_mapping>();
//                            auto constexpr t_mapping_block_size = t_mapping_vector_block.j_ - t_mapping_vector_block.i_;
//                            auto strain_block = strain.template segment<t_mapping_block_size>(t_mapping_vector_block.i_);
//                            t_MappingPolicy::non_linear(strain_block);
//                            if constexpr (t_i < t_finite_element.unknown_.mappings_.size() - 1) {
//                                self.template operator()<t_i + 1>(self);
//                            }
//                        };
//                        set_generalized_gradient(set_generalized_gradient);
//                        auto generalized_gradients = lolita::matrix::Span<t_Strain>(this->integration_points_[i].material_point_->s1.gradients.data());
//                        generalized_gradients = strain;
//                        /*
//                         * update the generalized flux, tangent operator, and internal state variables for the current integration point
//                         */
//                        auto const & behaviour = this->behaviour_->behaviour_;
//                        auto behaviour_view = mgis::behaviour::make_view(* this->integration_points_[i].material_point_);
//                        auto result = mgis::behaviour::integrate(behaviour_view, behaviour);
//                    }
//                }
//
//                lolita::matrix::Vector<lolita::real, getUnknownVectorSize()>
//                getInternalForces()
//                const
//                {
//                    using t_InternalForcesVector = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getNumUnknowns()>;
//                    using t_Stress = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows()>;
//                    auto internal_forces = lolita::matrix::Zero<t_InternalForcesVector>();
//                    for (int i = 0; i < t_FiniteElementTraits::Quadrature::dim_; ++i) {
//                        auto generalized_flux = lolita::matrix::Span<t_Stress>(this->integration_points_[i].material_point_->s1.thermodynamic_forces.data());
//                        internal_forces += this->integration_points_[i].weight_ * this->integration_points_[i].operator_.transpose() * generalized_flux;
//                    }
//                    internal_forces += this->cell_.stabilization_ * this->getUnknowns();
//                    return internal_forces;
//                }
//
//                lolita::matrix::Vector<lolita::real, getUnknownVectorSize()>
//                getExternalForces()
//                const
//                {
//                    using t_ExternalForcesVector = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getNumUnknowns()>;
//                    auto external_forces_vector = lolita::matrix::Zero<t_ExternalForcesVector>();
//                    auto count = 0;
//                    auto TIME___ = 0.0;
//                    auto set_cell_block = [&] () {
//                        using t_CellFieldUnknownsTraits = unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>;
//                        using t_CellUnknowns = typename t_CellFieldUnknownsTraits::FieldUnknowns;
//                        using t_CellUnknown = std::tuple_element_t<0, typename t_CellUnknowns::Unknowns>;
//                        for (int i = 0; i < t_FiniteElementTraits::Field::shape_.rows_; ++i) {
//                            for (int j = 0; j < t_FiniteElementTraits::Field::shape_.cols_; ++j) {
//                                for (int k = 0; k < t_FiniteElementTraits::Quadrature::dim_; ++k) {
//                                    auto external_cell_forces_vector = external_forces_vector.template segment<t_CellUnknown::dim_>(count);
//                                    auto cell_load_value = this->getLoad(i, j)->getImposedValue(this->integration_points_[k].coordinates_, TIME___);
//                                    auto const & generalized_load_operator = this->integration_points_[k].generalized_load_operator_;
//                                    external_cell_forces_vector += this->integration_points_[k].weight_ * cell_load_value * generalized_load_operator;
//                                }
//                                count += t_CellUnknown::dim_;
//                            }
//                        }
//                    };
//                    auto set_faces_block = [&] <lolita::integer _i = 0> (auto & self) mutable {
//                        auto const constexpr t_face = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<0, _i>();
//                        using t_FaceFieldUnknownsTraits = unknown::FiniteElementFieldUnknownsTraits<t_face, t_domain, t_finite_element>;
//                        using t_FaceUnknowns = typename t_FaceFieldUnknownsTraits::FieldUnknowns;
//                        using t_FaceUnknown = std::tuple_element_t<0, typename t_FaceUnknowns::Unknowns>;
//                        using t_FaceElementTraits = lolita::core2::finite_element::FiniteElementTraits<t_face, t_domain, t_finite_element>;
//                        auto const & faces = this->template getComponents<0, _i>();
//                        for (auto const & face : faces) {
//                            for (int i = 0; i < t_FaceElementTraits::Field::shape_.rows_; ++i) {
//                                for (int j = 0; j < t_FaceElementTraits::Field::shape_.cols_; ++j) {
//                                    for (int k = 0; k < t_FaceElementTraits::Quadrature::dim_; ++k) {
//                                        auto external_face_forces_vector = external_forces_vector.template segment<t_FaceUnknown::dim_>(count);
//                                        auto face_load_value = face->getLoad(i, j)->getImposedValue(face->integration_points_[k].coordinates_, TIME___);
//                                        auto const & generalized_load_operator = face->integration_points_[k].generalized_load_operator_;
//                                        external_face_forces_vector += face->integration_points_[k].weight_ * face_load_value * generalized_load_operator;
//                                    }
//                                    count += t_FaceUnknown::dim_;
//                                }
//                            }
//                        }
//                        if constexpr (_i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<0>() - 1) {
//                            self.template operator()<_i + 1>(self);
//                        }
//                    };
//                    set_cell_block();
//                    set_faces_block(set_faces_block);
//                    return external_forces_vector;
//                }
////
////                lolita::matrix::Matrix<lolita::real, getUnknownVectorSize(), getUnknownVectorSize()>
////                getTangentOperator()
////                const
////                {
////                    using _TangentOperatorMatrix = lolita::matrix::Matrix<lolita::real, t_FiniteElementTraits::getNumUnknowns(), t_FiniteElementTraits::getNumUnknowns()>;
////                    using _TangentOperatorTensor = lolita::matrix::Matrix<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows(), t_FiniteElementTraits::getGeneralizedStrainNumRows()>;
////                    auto tangent_operator_matrix = lolita::matrix::Zero<_TangentOperatorMatrix>();
////                    for (int i = 0; i < t_FiniteElementTraits::Quadrature::dim_; ++i) {
////                        auto count = 0;
////                        auto tangent_operator_tensor = lolita::matrix::Zero<_TangentOperatorTensor>();
////                        auto set_tangent_operator_tensor = [&] <lolita::integer _i = 0> (auto & self) mutable {
////                            auto const constexpr _mapping = t_finite_element.unknown_.mappings_[_i];
////                            auto const constexpr _mapping_vector_block = t_FiniteElementTraits::template getMappingBlock<_mapping>();
////                            auto const constexpr _mapping_block_size = _mapping_vector_block.j_ - _mapping_vector_block.i_;
////                            using _MappingPolicy = lolita::core2::field::MappingPolicy<t_finite_element.unknown_.tensor_, t_domain, _mapping>;
////                            using _MappingMatrix = lolita::matrix::Matrix<lolita::real, _mapping_block_size, _mapping_block_size>;
////                            auto tangent_operator_block = tangent_operator_tensor.template block<_mapping_block_size, _mapping_block_size>(
////                                    _mapping_vector_block.i_,
////                                    _mapping_vector_block.i_
////                            );
////                            auto tangent_operator_values = lolita::matrix::Span<_MappingMatrix const>(this->integration_points_[i].material_point_->K.data() + count);
////                            tangent_operator_block += tangent_operator_values;
////                            count += _mapping_block_size * _mapping_block_size;
////                            if constexpr (_i < t_finite_element.unknown_.mappings_.size() - 1) {
////                                self.template operator()<_i + 1>(self);
////                            }
////                        };
////                        set_tangent_operator_tensor(set_tangent_operator_tensor);
////                        tangent_operator_matrix += this->integration_points_[i].weight_ * this->integration_points_[i].generalized_strain_operator_.transpose() * tangent_operator_tensor * this->integration_points_[i].generalized_strain_operator_;
////                    }
////                    tangent_operator_matrix += this->cell_.stabilization_;
////                    return tangent_operator_matrix;
////                }
////
////                template<auto _arg>
////                void
////                assemble(
////                        lolita::core2::mesh::Mesh<t_domain, _arg> & mesh
////                )
////                {
////                    auto const constexpr _cs = getNumCellUnknowns();
////                    auto const constexpr _fs = getNumFaceUnknowns();
////                    //
////                    using FullMatrix = lolita::matrix::Matrix<lolita::real, getNumUnknowns(), getNumUnknowns()>;
////                    using FullVector = lolita::matrix::Vector<lolita::real, getNumUnknowns()>;
////                    //
////                    using CondMatrix = lolita::matrix::Matrix<lolita::real, getNumFaceUnknowns(), getNumFaceUnknowns()>;
////                    using CondVector = lolita::matrix::Vector<lolita::real, getNumFaceUnknowns()>;
////                    //
////                    using CellCellMatrixBlock = lolita::matrix::Matrix<lolita::real, getNumCellUnknowns(), getNumCellUnknowns()>;
////                    using CellFaceMatrixBlock = lolita::matrix::Matrix<lolita::real, getNumCellUnknowns(), getNumFaceUnknowns()>;
////                    using FaceCellMatrixBlock = lolita::matrix::Matrix<lolita::real, getNumFaceUnknowns(), getNumCellUnknowns()>;
////                    using FaceFaceMatrixBlock = lolita::matrix::Matrix<lolita::real, getNumFaceUnknowns(), getNumFaceUnknowns()>;
////                    //
////                    using CellVectorBlock = lolita::matrix::Vector<lolita::real, getNumCellUnknowns()>;
////                    using FaceVectorBlock = lolita::matrix::Vector<lolita::real, getNumFaceUnknowns()>;
////                    //
////                    auto residual_vector = getInternalForces() - getExternalForces();
////                    auto tangent_matrix = getTangentOperator();
////                    auto Ktt = tangent_matrix.template block<_cs, _cs>(0, 0);
////                    auto Ktf = tangent_matrix.template block<_cs, _fs>(0, _cs);
////                    auto Kft = tangent_matrix.template block<_fs, _cs>(_cs, 0);
////                    auto Kff = tangent_matrix.template block<_fs, _fs>(_cs, _cs);
////                    auto Rt = residual_vector.template segment<_cs>(0);
////                    auto Rf = residual_vector.template segment<_fs>(_cs);
////                    auto Ktt_inv = Ktt.llt().solve(decltype(Ktt)::Identity());
////                    auto Kc = Kff - Kft * Ktt_inv * Ktf;
////                    auto Rc = Rf - Kft * Ktt_inv * Rt;
////                    this->kft_ = Kft;
////                    this->ktt_inv_ = Ktt_inv;
////                    auto constexpr _finite_element_index = _arg.template getFiniteElementIndex<t_finite_element>();
////                    auto & system = mesh.systems_[_finite_element_index];
////                    auto row_offset = 0;
////                    auto col_offset = 0;
////                    auto set_rows = [&] <lolita::integer _i = 0> (auto & t_set_rows) mutable {
////                        auto const & faces_row = this->template getComponents<0, _i>();
////                        auto const constexpr _face_r = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<0, _i>();
////                        using _RowFaceUnknowns = typename lolita::core2::finite_element::unknown::FiniteElementFieldUnknownsTraits<_face_r, t_domain, t_finite_element>::FieldUnknowns;
////                        using _RowFaceUnknown = std::tuple_element_t<0, typename _RowFaceUnknowns::Unknowns>;
//////                        using _RowFiniteElementFace = lolita::core2::finite_element::face::FiniteElementFace<_face_r, t_domain, t_finite_element>;
////                        using _RowFiniteElementFace = lolita::core2::finite_element::FiniteElementTraits<_face_r, t_domain, t_finite_element>;
////                        auto set_cols = [&] <lolita::integer _j = 0> (auto & t_set_cols) mutable {
////                            auto const & faces_col = this->template getComponents<0, _j>();
////                            auto const constexpr _face_c = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<0, _j>();
////                            using _ColFaceUnknowns = typename lolita::core2::finite_element::unknown::FiniteElementFieldUnknownsTraits<_face_c, t_domain, t_finite_element>::FieldUnknowns;
////                            using _ColFaceUnknown = std::tuple_element_t<0, typename _ColFaceUnknowns::Unknowns>;
//////                            using _ColFiniteElementFace = lolita::core2::finite_element::face::FiniteElementFace<_face_c, t_domain, t_finite_element>;
////                            using _ColFiniteElementFace = lolita::core2::finite_element::FiniteElementTraits<_face_c, t_domain, t_finite_element>;
////                            for (auto const & face_row : faces_row) {
////                                for (auto const & face_col : faces_col) {
////                                    for (int i = 0; i < t_FiniteElementTraits::Field::shape_.rows_; ++i) {
////                                        for (int j = 0; j < t_FiniteElementTraits::Field::shape_.cols_; ++j) {
////                                            auto const & i_rows = face_row->unknwons_.template getUnknown<0>().getUnknownComponent(i, j)->unknowns_coordinates_;
////                                            auto const & i_cols = face_col->unknwons_.template getUnknown<0>().getUnknownComponent(i, j)->unknowns_coordinates_;
////                                            for (auto const & idx_row : i_rows) {
////                                                for (auto const & idx_col : i_cols) {
////                                                    auto triplet = Eigen::Triplet<lolita::real>(row_offset, col_offset, Kc(row_offset, col_offset));
////                                                    system.lhs_values.push_back(Eigen::Triplet<lolita::real>(idx_row, idx_col, Kc(row_offset, col_offset)));
////                                                    system.lhs_values.push_back(lolita::matrix::SparseMatrixInput(idx_row, idx_col, Kc(row_offset, col_offset)));
////                                                    Kc(row_offset, col_offset);
////                                                    col_offset += 1;
////                                                }
////                                                system.rhs_values.push_back(lolita::matrix::SparseVectorInput(idx_row, Rc(row_offset)));
////                                                Rc(row_offset);
////                                                row_offset += 1;
////                                            }
////                                        }
////                                    }
////                                }
////                            }
////                            if constexpr (_j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<0>() - 1) {
////                                t_set_cols.template operator()<_j + 1>(t_set_cols);
////                            }
////                        };
////                        if constexpr (_i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<0>() - 1) {
////                            t_set_rows.template operator()<_i + 1>(t_set_rows);
////                        }
////                    };
////                }
////
////                void
////                setElement()
////                {
////                    for (int i = 0; i < t_FiniteElementTraits::Quadrature::dim_; ++i) {
////                        auto const & bhv = this->behaviour_->behaviour_;
////                        this->integration_points_[i].coordinates_ = this->template getCurrentQuadraturePoint<t_FiniteElementTraits::Quadrature::quadrature_, t_FiniteElementTraits::Quadrature::ord_>(i);
////                        std::cout << "my coordinates are : " << i << std::endl;
////                        std::cout << this->integration_points_[i].coordinates_ << std::endl;
////                        this->integration_points_[i].material_point_ = std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(bhv));
////                        this->integration_points_[i].material_point_->K[0] = 4;
////                        auto behaviour_view = mgis::behaviour::make_view(* this->integration_points_[i].material_point_);
////                        auto result = mgis::behaviour::integrate(behaviour_view, bhv);
////                        getUnknowns();
////                        auto res = this->integration_points_[i].operator_ * getUnknowns();
////                        setGeneralizedGradients();
////                        getInternalForces();
////                        getExternalForces();
////                        getTangentOperator();
////                    }
////                }
//
//            };
//
//        };
//
//    }
//
//    /**
//     * @brief
//     * @tparam t_element
//     * @tparam t_domain
//     * @tparam t_finite_element
//     */
//    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//    struct FiniteElementCell
//    {
//
//    private:
//
//        /**
//         * @brief
//         */
//        using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element>;
//
//        /**
//         * @brief
//         */
//        using t_FiniteElementCellTraits = cell::FiniteElementCellTraits<t_element, t_domain, t_finite_element>;
//
//        /**
//         * @brief
//         */
//        using t_Quadrature = typename FiniteElementTraits<t_element, t_domain, t_finite_element>::Quadrature;
//
//    public:
//
//        /**
//         * @brief
//         */
//        lolita::integer static constexpr dim_structural_unknowns = t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Structural()>();
//
//        /**
//         * @brief
//         */
//        lolita::integer static constexpr dim_subsidiary_unknowns = t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Subsidiary()>();
//
//        /**
//         * @brief
//         */
//        lolita::integer static constexpr num_structural_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Structural()>();
//
//        /**
//         * @brief
//         */
//        lolita::integer static constexpr num_subsidiary_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
//
//        /**
//         * @brief
//         */
//        struct IntegrationPoint
//        {
//
//            /**
//             * @brief
//             */
//            lolita::domain::Point coordinates_;
//
//            /**
//             * @brief
//             */
//            lolita::real weight_;
//
//            /**
//             * @brief
//             */
//            std::unique_ptr<mgis::behaviour::BehaviourData> material_point_;
//
//            /**
//             * @brief
//             */
//            lolita::matrix::Matrix<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows(), num_structural_unknowns> structural_strain_operator_;
//
//            /**
//             * @brief
//             */
//            lolita::matrix::Matrix<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows(), num_subsidiary_unknowns> subsidiary_strain_operator_;
//
//            /**
//             * @brief
//             */
//            lolita::matrix::Vector<lolita::real, dim_structural_unknowns> structural_load_operator_;
//
//            /**
//             * @brief
//             */
//            lolita::matrix::Vector<lolita::real, dim_subsidiary_unknowns> subsidiary_load_operator_;
//
//        };
//
//        /**
//         * @brief
//         */
//        struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
//        {
//
//            void
//            setIntegrationPoint()
//            {
//                auto const & behaviour = this->behaviour_->behaviour_;
//                for (auto & integration_point : this->integration_points_) {
//                    integration_point.material_point_ = std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(behaviour));
//                    integration_point.material_point_->K[0] = 4;
//                }
//            }
//
//        };
//
//        void
//        setIntegrationPoint()
//        {
//            static_cast<Implementation *>(this)->setIntegrationPoint();
//        }
//
//        void
//        setGeneralizedGradients()
//        {
//            setIntegrationPoint();
//            using t_Implementation = typename t_FiniteElementCellTraits::Implementation;
//            static_cast<t_Implementation *>(this)->setGeneralizedGradients();
//        }
//
//        /**
//         * @brief
//         */
//        std::array<IntegrationPoint, t_Quadrature::dim_> integration_points_;
//
//    };

}

#endif //LOLITA_LOLITA_CORE_5_004_MODULE_HXX
