#ifndef FEB10822_D09A_4468_9437_2ADB758EE6F4
#define FEB10822_D09A_4468_9437_2ADB758EE6F4

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/100_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/200_quadrature.hxx"
#include "core/201_finite_element_dof.hxx"
#include "core/202_finite_element_frm.hxx"
#include "core/300_finite_element.hxx"
#include "core/400_finite_element_basis.hxx"

namespace lolita
{

    template<DiscretizationConcept auto t_discretization>
    requires(t_discretization.is("HybridDiscontinuousGalerkin"))
    struct DiscretizationTraits3<t_discretization>
    {

    private:

        static constexpr
        Basis
        getCellBasis()
        {
            return t_discretization.getCellBasis();
        }

        static constexpr
        Basis
        getFaceBasis()
        {
            return t_discretization.getFaceBasis();
        }

        static constexpr
        Basis
        getGradBasis()
        {
            return t_discretization.getGradBasis();
        }

        static constexpr
        Basis
        getPotentialBasis()
        {
            return t_discretization.getFaceBasis().toOrder(t_discretization.getFaceBasis().getOrd() + 1);
        }

        template<Element t_element>
        static constexpr
        Integer
        getCellBasisSize()
        {
            return BasisTraits<getCellBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getFaceBasisSize()
        {
            return BasisTraits<getFaceBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getGradBasisSize()
        {
            return BasisTraits<getGradBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getPotentialBasisSize()
        {
            return BasisTraits<getPotentialBasis()>::template getSize<t_element>();
        }

    public:
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumElementCoefficients()
        {
            return 0;
        }
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumElementCoefficients()
        requires(t_element.isSub(t_domain, 0))
        {
            return BasisTraits<t_discretization.getCellBasis()>::template getSize<t_element>();
        }
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumElementCoefficients()
        requires(t_element.isSub(t_domain, 1))
        {
            return BasisTraits<t_discretization.getFaceBasis()>::template getSize<t_element>();
        }

        template<Element t_element, Domain t_domain, FieldConcept auto t_field>
        struct Implementation : FiniteElement<t_element, t_domain>
        {

        private:

            using Base = FiniteElement<t_element, t_domain>;

            using Traits = FiniteElementTraits<t_element, t_domain>;

        public:
            
            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, Traits::template getNumUnknownCoefficients<t_field>(), Traits::template getNumUnknownCoefficients<t_field>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "Stabilization")
            {
                return getStabilization();
            }
            
            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), Traits::template getNumUnknownCoefficients<t_field>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "Gradient" || t_mapping.getTransformation() == "LargeStrain")
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto rhs = getGradientRhs(mapping_value.row(), mapping_value.col());
                    auto line = mapping.template block<1, Traits::template getNumUnknownCoefficients<t_field>()>(mapping_value.rank(), 0);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return mapping;
            }

            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), Traits::template getNumUnknownCoefficients<t_field>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "SmallStrain")
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto rhs = getSymmetricGradientRhs(mapping_value.row(), mapping_value.col());
                    auto line = mapping.template block<1, Traits::template getNumUnknownCoefficients<t_field>()>(mapping_value.rank(), 0);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return mapping;
            }

            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), Traits::template getNumUnknownCoefficients<t_field>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "Identity")
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                mapping.setZero();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto left_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                    auto col_offset = mapping_value.row() * getCellBasisSize<t_element>();
                    auto line = mapping.template block<1, getCellBasisSize<t_element>()>(mapping_value.rank(), col_offset);
                    line = mapping_value.value() * this->template getBasisEvaluation<getCellBasis()>(point);
                }
                return mapping;
            }

            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, 1, Traits::template getNumUnknownCoefficients<t_field>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "Identity" && (t_mapping.getRow() >= 0 || t_mapping.getCol() >= 0))
            {
                auto mapping = DenseMatrix<Real, 1, Traits::template getNumUnknownCoefficients<t_field>()>();
                mapping.setZero();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    if (mapping_value.row() == t_mapping.getRow() && mapping_value.col() == t_mapping.getCol())
                    {
                        auto left_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                        auto col_offset = mapping_value.row() * getCellBasisSize<t_element>();
                        auto line = mapping.template block<1, getCellBasisSize<t_element>()>(0, col_offset);
                        line = mapping_value.value() * this->template getBasisEvaluation<getCellBasis()>(point);
                    }
                }
                return mapping;
            }
            
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            {}
            
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr num_cell_unknowns = FiniteElementTraits<t_element, t_domain>::template getNumElementCoefficients<t_field>();
                auto constexpr strain_operator_num_cols = Traits::template getNumUnknownCoefficients<t_field>();
                auto constexpr num_face_unknowns = Traits::template getNumUnknownCoefficients<t_field>() - num_cell_unknowns;
                auto faces_correction = DenseVector<Real, num_face_unknowns>();
                faces_correction.setZero();
                auto faces_correction_offset = 0;
                auto set_faces_increment = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_dof = face->template getDiscreteField<t_field>().getDegreeOfFreedom();
                        auto lhs = faces_correction.template segment<face_dof.template getSize<t_field, getFaceBasis()>()>(faces_correction_offset);
                        auto rhs = increment.template segment<face_dof.template getSize<t_field, getFaceBasis()>()>(face_dof.getOffset());
                        lhs = rhs;
                        faces_correction_offset += face_dof.template getSize<t_field, getFaceBasis()>();
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_increment(set_faces_increment);
                // auto k_tt_inv = this->template getDiscreteField<t_field>().getMatrix("KTT");
                // auto k_tf = this->template getDiscreteField<t_field>().getMatrix("KTF");
                // auto r_t = this->template getDiscreteField<t_field>().getVector("RT");
                // auto cell_corr = k_tt_inv * (r_t - k_tf * faces_correction);
                // this->template getDiscreteField<t_field>().getDegreeOfFreedom().addCoefficients(cell_corr);
            }
            
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            requires(t_element.isSub(t_domain, 1))
            {
                static_cast<Base *>(this)->template upgradeDiscreteFieldDegreeOfFreedom<t_field, getFaceBasis()>(increment);
            }
            
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            {
                return 0.0;
            }
            
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                return static_cast<Base const *>(this)->template getDiscreteFieldDegreeOfFreedomValue<t_field, getCellBasis()>(point, row, col);
            }
            
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                Point const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                return static_cast<Base const *>(this)->template getDiscreteFieldDegreeOfFreedomValue<t_field, getFaceBasis()>(point, row, col);
            }

            template<Quadrature t_quadrature>
            Real
            getDiscreteFieldDegreeOfFreedomIntegrationValue(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                // this->template getDiscreteFieldDegreeOfFreedomValue<t_field>(row, col);
                auto value = Real(0);
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    value += weight * this->template getDiscreteFieldDegreeOfFreedomValue<t_field>(point, row, col);
                }
                return value;
            }

            DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>
            getFieldDualVector(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto external_forces = DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>();
                auto vector = this->template getBasisEvaluation<getCellBasis()>(point);
                auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::template getCols<t_domain>() * row + col);
                external_forces.setZero();
                external_forces.template segment<getCellBasisSize<t_element>()>(offset) = vector;
                return external_forces;
            }

            DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>
            getFieldDualVector(
                PointConcept auto const & point
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto external_forces = DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>();
                auto vector = this->template getBasisEvaluation<getCellBasis()>(point);
                external_forces.setZero();
                for (auto row = 0; row < FieldTraits<t_field>::template getRows<t_domain>(); row++)
                {
                    for (auto col = 0; col < FieldTraits<t_field>::template getCols<t_domain>(); col++)
                    {
                        auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::template getCols<t_domain>() * row + col);
                        external_forces.template segment<getCellBasisSize<t_element>()>(offset) = vector;
                    }
                }
                return external_forces;
            }

            DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>
            getFieldDualVector(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto external_forces = DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>();
                auto vector = this->template getBasisEvaluation<getFaceBasis()>(point);
                auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::template getCols<t_domain>() * row + col);
                external_forces.setZero();
                external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = vector;
                return external_forces;
            }

            // DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>
            // getFieldDualVector(
            //     PointConcept auto const & point
            // )
            // const
            // requires(t_element.isSub(t_domain, 1))
            // {
            //     auto external_forces = DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>();
            //     auto vector = this->template getBasisEvaluation<getFaceBasis()>(point);
            //     external_forces.setZero();
            //     for (auto row = 0; row < FieldTraits<t_field>::getRows(); row++)
            //     {
            //         for (auto col = 0; col < FieldTraits<t_field>::getCols(); col++)
            //         {
            //             auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
            //             external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = vector;
            //         }
            //     }
            //     return external_forces;
            // }

            // DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>
            // getFieldPrimalVector(
            //     PointConcept auto const & point,
            //     Integer row,
            //     Integer col
            // )
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto external_forces = DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>();
            //     auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
            //     auto field_coefficients = cell_dof.template getCoefficients<t_element, t_field, getCellBasis()>(row, col);
            //     auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
            //     external_forces.setZero();
            //     external_forces.template segment<getCellBasisSize<t_element>()>(offset) = field_coefficients;
            //     return external_forces;
            // }

            // DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>
            // getFieldPrimalVector(
            //     PointConcept auto const & point
            // )
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto external_forces = DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>();
            //     auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
            //     external_forces.setZero();
            //     for (auto row = 0; row < FieldTraits<t_field>::getRows(); row++)
            //     {
            //         for (auto col = 0; col < FieldTraits<t_field>::getCols(); col++)
            //         {
            //             auto field_coefficients = cell_dof.template getCoefficients<t_element, t_field, getCellBasis()>(row, col);
            //             auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
            //             external_forces.template segment<getCellBasisSize<t_element>()>(offset) = field_coefficients;
            //         }
            //     }
            //     return external_forces;
            // }

            // DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>
            // getFieldPrimalVector(
            //     PointConcept auto const & point,
            //     Integer row,
            //     Integer col
            // )
            // const
            // requires(t_element.isSub(t_domain, 1))
            // {
            //     auto external_forces = DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>();
            //     auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
            //     auto field_coefficients = cell_dof.template getCoefficients<t_element, t_field, getFaceBasis()>(row, col);
            //     auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
            //     external_forces.setZero();
            //     external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = field_coefficients;
            //     return external_forces;
            // }

            // DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>
            // getFieldPrimalVector(
            //     PointConcept auto const & point
            // )
            // const
            // requires(t_element.isSub(t_domain, 1))
            // {
            //     auto external_forces = DenseVector<Real, Traits::template getNumUnknownCoefficients<t_field>()>();
            //     auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
            //     external_forces.setZero();
            //     for (auto row = 0; row < FieldTraits<t_field>::getRows(); row++)
            //     {
            //         for (auto col = 0; col < FieldTraits<t_field>::getCols(); col++)
            //         {
            //             auto field_coefficients = cell_dof.template getCoefficients<t_element, t_field, getFaceBasis()>(row, col);
            //             auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
            //             external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = field_coefficients;
            //         }
            //     }
            //     return external_forces;
            // }

            /**
             * Implementation
             * *************************************************************************************************************************************************
             */

            static constexpr
            Integer
            getGradientConstructionQuadratureOrder()
            {
                if constexpr (t_domain.isAxiSymmetric())
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1, getGradBasis().getOrd()) + 1;
                }
                else
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1, getGradBasis().getOrd());
                }
            }

            static constexpr
            Integer
            getPotentialConstructionQuadratureOrder()
            {
                if constexpr (t_domain.isAxiSymmetric())
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1) + 1;
                }
                else
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1);
                }
            }
            
            DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>
            getPotentialLhs()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto lhs = DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>();
                lhs.setZero();
                for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
                {
                    for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                        auto vector = this->template getBasisDerivative<getPotentialBasis()>(point, i_component);
                        auto vector_j = vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        lhs += weight * vector_j * vector_j.transpose();
                    }
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }
            
            DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, Traits::template getNumUnknownCoefficients<t_field>()>
            getPotentialRhs(
                Integer row
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto rhs = DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, Traits::template getNumUnknownCoefficients<t_field>()>();
                rhs.setZero();
                for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
                {
                    auto face_offset = t_field_size * getCellBasisSize<t_element>();
                    for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                        auto row_cell_vector = this->template getBasisDerivative<getPotentialBasis()>(point, i_component);
                        auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        auto col_cell_vector = this->template getBasisDerivative<getCellBasis()>(point, i_component);
                        auto col_offset = row * getCellBasisSize<t_element>();
                        auto block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getCellBasisSize<t_element>()>(0, col_offset);
                        block += weight * row_cell_vector_j * col_cell_vector.transpose();
                    }
                    auto set_faces_blocks = [&] <Integer t_i = 0> (
                        auto & self
                    )
                    constexpr mutable
                    {
                        auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                        auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                        auto i_face = 0;
                        for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                        {
                            auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                            auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                            for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                            {
                                auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                                auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                auto argf = i;
                                //
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                                auto normal_vector = face->getNormalVector(point);
                                //
                                auto row_cell_vector = this->template getBasisDerivative<getPotentialBasis()>(i_p, i_component);
                                auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                                auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                                auto normal_vector_component = face_orientation * normal_vector(i_component);
                                auto face_col_offset = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                                auto face_block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                                face_block += normal_vector_component * weight * row_cell_vector_j * col_face_vector.transpose();
                                //
                                auto cell_col_offset = row * getCellBasisSize<t_element>();
                                auto cell_block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getCellBasisSize<t_element>()>(0, cell_col_offset);
                                cell_block -= normal_vector_component * weight * row_cell_vector_j * col_cell_vector.transpose();
                            }
                            face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                            i_face ++;
                        }
                        if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                        {
                            self.template operator()<t_i + 1>(self);
                        }
                    };
                    set_faces_blocks(set_faces_blocks);
                }
                return rhs;
            }
            
            DenseMatrix<Real, Traits::template getNumUnknownCoefficients<t_field>(), Traits::template getNumUnknownCoefficients<t_field>()>
            getStabilization()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                //
                auto S_T_A_B_I_L_I_Z_A_T_I_O_N = DenseMatrix<Real, Traits::template getNumUnknownCoefficients<t_field>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                auto potential_operator = DenseMatrix<Real, getPotentialBasisSize<t_element>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                auto potential_lhs = getPotentialLhs();
                auto denom = 1.0 / QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize();
                S_T_A_B_I_L_I_Z_A_T_I_O_N.setZero();
                for (auto i_component = 0; i_component < FieldTraits<t_field>::template getSize<t_domain>(); i_component++)
                {
                    auto face_offset = t_field_size * getCellBasisSize<t_element>();
                    potential_operator.setZero();
                    auto oppp = potential_operator.template block<getPotentialBasisSize<t_element>() - 1, Traits::template getNumUnknownCoefficients<t_field>()>(1, 0);
                    oppp = potential_lhs * getPotentialRhs(i_component);
                    for (auto i_quadrature = 0; i_quadrature < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i_quadrature++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                        auto row_cell_vector = this->template getBasisEvaluation<getPotentialBasis()>(point);
                        auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        potential_operator.template block<1, Traits::template getNumUnknownCoefficients<t_field>()>(0, 0) -= denom * weight * row_cell_vector_j.transpose() * oppp;
                        potential_operator.template block<1, getCellBasisSize<t_element>()>(0, i_component * getCellBasisSize<t_element>()) +=
                        denom * weight * this->template getBasisEvaluation<getCellBasis()>(point).transpose();
                    }
                    //
                    auto cell_projection_lhs = DenseMatrix<Real, getCellBasisSize<t_element>(), getCellBasisSize<t_element>()>();
                    auto cell_projection_rhs = DenseMatrix<Real, getCellBasisSize<t_element>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                    cell_projection_lhs.setZero();
                    cell_projection_rhs.setZero();
                    for (auto i_quadrature = 0; i_quadrature < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i_quadrature++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                        auto row_cell_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                        cell_projection_lhs += weight * row_cell_vector * row_cell_vector.transpose();
                        cell_projection_rhs += weight * row_cell_vector * this->template getBasisEvaluation<getPotentialBasis()>(point).transpose() * potential_operator;
                    }
                    auto proj = cell_projection_lhs.llt().solve(decltype(cell_projection_lhs)::Identity()) * cell_projection_rhs;
                    //
                    auto set_faces_blocks = [&] <Integer t_i = 0> (
                        auto & self
                    )
                    constexpr mutable
                    {
                        auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                        auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                        auto i_face = 0;
                        for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                        {
                            auto constexpr quadrature_size = QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize();
                            auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                            auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                            auto face_projection_lhs = DenseMatrix<Real, getFaceBasisSize<t_inner_neighbor>(), getFaceBasisSize<t_inner_neighbor>()>();
                            auto face_projection_rhs = DenseMatrix<Real, getFaceBasisSize<t_inner_neighbor>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                            face_projection_lhs.setZero();
                            face_projection_rhs.setZero();
                            for (auto i_quadrature = 0; i_quadrature < quadrature_size; i_quadrature++)
                            {
                                auto argh = face_sign == 1 ? i_quadrature : quadrature_size - (i_quadrature + 1);
                                auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                auto argf = i_quadrature;
                                //
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                                auto normal_vector = face->getNormalVector(point);
                                auto row_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto col_potential_vector = this->template getBasisEvaluation<getPotentialBasis()>(i_p);
                                auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                                //
                                face_projection_lhs += weight * row_face_vector * row_face_vector.transpose();
                                face_projection_rhs -= weight * row_face_vector * col_potential_vector.transpose() * potential_operator;
                                face_projection_rhs += weight * row_face_vector * col_cell_vector.transpose() * proj;
                                //
                                auto cell_col_offset = i_component * getCellBasisSize<t_element>();
                                auto cell_block = face_projection_rhs.template block<getFaceBasisSize<t_inner_neighbor>(), getCellBasisSize<t_element>()>(0, cell_col_offset);
                                cell_block -= weight * row_face_vector * col_cell_vector.transpose();
                                //
                                auto face_col_offset = face_offset + i_component * getFaceBasisSize<t_inner_neighbor>();
                                auto face_block = face_projection_rhs.template block<getFaceBasisSize<t_inner_neighbor>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                                face_block += weight * row_face_vector * row_face_vector.transpose();
                            }
                            auto fa_ce_proj = face_projection_lhs.llt().solve(decltype(face_projection_lhs)::Identity()) * face_projection_rhs;
                            for (auto i_quadrature = 0; i_quadrature < quadrature_size; i_quadrature++)
                            {
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                                auto row_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto V_E_C_T_O_R = row_face_vector.transpose() * fa_ce_proj;
                                auto factor = (1.0 / face->getLocalFrameDiameters().norm());
                                S_T_A_B_I_L_I_Z_A_T_I_O_N += factor * weight * V_E_C_T_O_R.transpose() * V_E_C_T_O_R;
                            }
                            face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                            i_face ++;
                        }
                        if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                        {
                            self.template operator()<t_i + 1>(self);
                        }
                    };
                    set_faces_blocks(set_faces_blocks);
                }
                return S_T_A_B_I_L_I_Z_A_T_I_O_N;
            }
            
            DenseMatrix<Real, getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>
            getGradientLhs()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto lhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>();
                lhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    lhs += weight * vector * vector.transpose();
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }
            
            DenseMatrix<Real, getGradBasisSize<t_element>(), Traits::template getNumUnknownCoefficients<t_field>()>
            getGradientRhs(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = DenseMatrix<Real, getGradBasisSize<t_element>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                rhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    auto col_cell_vector = this->template getBasisDerivative<getCellBasis()>(point, col);
                    auto col_offset = row * getCellBasisSize<t_element>();
                    auto block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset);
                    block += weight * row_cell_vector * col_cell_vector.transpose();
                }
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                            auto argf = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            //
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            auto normal_vector = face->getNormalVector(point);
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            auto normal_vector_component = face_orientation * normal_vector(col);
                            auto face_col_offset = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                            face_block += normal_vector_component * weight * row_cell_vector * col_face_vector.transpose();
                            auto cell_col_offset = row * getCellBasisSize<t_element>();
                            auto cell_block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset);
                            cell_block -= normal_vector_component * weight * row_cell_vector * col_cell_vector.transpose();
                        }
                        face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                        i_face ++;
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator()<t_i + 1>(self);
                    }
                };
                set_faces_blocks(set_faces_blocks);
                return rhs;
            }
            
            DenseMatrix<Real, getGradBasisSize<t_element>(), Traits::template getNumUnknownCoefficients<t_field>()>
            getSymmetricGradientRhs(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = DenseMatrix<Real, getGradBasisSize<t_element>(), Traits::template getNumUnknownCoefficients<t_field>()>();
                rhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    auto col_cell_vector_i = this->template getBasisDerivative<getCellBasis()>(point, col);
                    auto col_cell_vector_j = this->template getBasisDerivative<getCellBasis()>(point, row);
                    auto col_offset_i = row * getCellBasisSize<t_element>();
                    auto block_i = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset_i);
                    block_i += (1.0 / 2.0) * weight * row_cell_vector * col_cell_vector_i.transpose();
                    auto col_offset_j = col * getCellBasisSize<t_element>();
                    auto block_j = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset_j);
                    block_j += (1.0 / 2.0) * weight * row_cell_vector * col_cell_vector_j.transpose();
                }
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                            auto argf = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            auto normal_vector = face->getNormalVector(point);
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            auto normal_vector_component_i = face_orientation * normal_vector(col);
                            auto normal_vector_component_j = face_orientation * normal_vector(row);
                            auto face_col_offset_i = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block_i = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset_i);
                            face_block_i += (1.0 / 2.0) * normal_vector_component_i * weight * row_cell_vector * col_face_vector.transpose();
                            //
                            auto face_col_offset_j = face_offset + col * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block_j = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset_j);
                            face_block_j += (1.0 / 2.0) * normal_vector_component_j * weight * row_cell_vector * col_face_vector.transpose();
                            //
                            auto cell_col_offset_i = row * getCellBasisSize<t_element>();
                            auto cell_block_i = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset_i);
                            cell_block_i -= (1.0 / 2.0) * normal_vector_component_i * weight * row_cell_vector * col_cell_vector.transpose();
                            //
                            auto cell_col_offset_j = col * getCellBasisSize<t_element>();
                            auto cell_block_j = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset_j);
                            cell_block_j -= (1.0 / 2.0) * normal_vector_component_j * weight * row_cell_vector * col_cell_vector.transpose();
                        }
                        face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                        i_face ++;
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator()<t_i + 1>(self);
                    }
                };
                set_faces_blocks(set_faces_blocks);
                return rhs;
            }

        };
        
    };

    template<DiscreteFieldConcept auto t_field>
    requires(t_field.getDiscretization().is("HybridDiscontinuousGalerkin"))
    struct DiscretizationTraits2<t_field>
    {

    private:

        static constexpr
        Basis
        getCellBasis()
        {
            return t_field.getDiscretization().getCellBasis();
        }

        static constexpr
        Basis
        getFaceBasis()
        {
            return t_field.getDiscretization().getFaceBasis();
        }

        static constexpr
        Basis
        getGradBasis()
        {
            return t_field.getDiscretization().getGradBasis();
        }

        static constexpr
        Basis
        getPotentialBasis()
        {
            return t_field.getDiscretization().getFaceBasis().toOrder(t_field.getDiscretization().getFaceBasis().getOrd() + 1);
        }

        template<Element t_element>
        static constexpr
        Integer
        getCellBasisSize()
        {
            return BasisTraits<getCellBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getFaceBasisSize()
        {
            return BasisTraits<getFaceBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getGradBasisSize()
        {
            return BasisTraits<getGradBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getPotentialBasisSize()
        {
            return BasisTraits<getPotentialBasis()>::template getSize<t_element>();
        }

        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumCellUnknowns()
        requires(t_element.isSub(t_domain, 0))
        {
            return FieldTraits<t_field>::template getSize<t_domain>() * BasisTraits<t_field.getDiscretization().getCellBasis()>::template getSize<t_element>();
        }

        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumFaceUnknowns()
        requires(t_element.isSub(t_domain, 1))
        {
            return FieldTraits<t_field>::template getSize<t_domain>() * BasisTraits<t_field.getDiscretization().getFaceBasis()>::template getSize<t_element>();
        }

    public:
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getStaticSize()
        {
            return 0;
        }

        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getStaticSize()
        requires(t_element.isSub(t_domain, 0))
        {
            auto num_element_unknowns = getNumCellUnknowns<t_element, t_domain>();
            auto set_num_faces_unknowns = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                auto constexpr num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                num_element_unknowns += getNumFaceUnknowns<inner_neighbor, t_domain>() * num_inner_neighbors;
                if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_num_faces_unknowns(set_num_faces_unknowns);
            return num_element_unknowns;
        }
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getStaticSize()
        requires(t_element.isSub(t_domain, 1))
        {
            return getNumFaceUnknowns<t_element, t_domain>();
        }

        /**
         * Scalar num coef
         * Unknown
         * Coefficient
         * UnknownCoefficients
         * ElementCoefficients
         * Component
         * getNumElementFieldCoefficients
         * getNumElementScalarCoefficients
         * 
         * *****************************************************************************************************************************************************
         */
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumElementDegreesOfFreedom()
        {
            return 0;
        }
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumElementDegreesOfFreedom()
        requires(t_element.isSub(t_domain, 0))
        {
            return BasisTraits<t_field.getDiscretization().getCellBasis()>::template getSize<t_element>();
            // return getNumCellUnknowns<t_element, t_domain>();
        }
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumElementDegreesOfFreedom()
        requires(t_element.isSub(t_domain, 1))
        {
            return BasisTraits<t_field.getDiscretization().getFaceBasis()>::template getSize<t_element>();
            // return getNumFaceUnknowns<t_element, t_domain>();
        }

        /**
         * Field num coef
         * *****************************************************************************************************************************************************
         */
        
        template<Element t_element, Domain t_domain>
        static constexpr
        Integer
        getNumElementDegreesOfFreedom2()
        {
            return getNumElementDegreesOfFreedom<t_element, t_domain>() * FieldTraits<t_field>::template getSize<t_domain>();
        }

        template<Element t_element, Domain t_domain>
        struct Implementation : FiniteElement<t_element, t_domain>
        {

        private:

            using Base = FiniteElement<t_element, t_domain>;

        public:

            Integer
            getDynamicSize()
            const
            {
                return DiscretizationTraits2::template getStaticSize<t_element, t_domain>();
            }
            
            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, getStaticSize<t_element, t_domain>(), getStaticSize<t_element, t_domain>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "Stabilization")
            {
                return getStabilization();
            }
            
            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getStaticSize<t_element, t_domain>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "Gradient" || t_mapping.getTransformation() == "LargeStrain")
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getStaticSize<t_element, t_domain>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto rhs = getGradientRhs(mapping_value.row(), mapping_value.col());
                    auto line = mapping.template block<1, getStaticSize<t_element, t_domain>()>(mapping_value.rank(), 0);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return mapping;
            }

            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getStaticSize<t_element, t_domain>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "SmallStrain")
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getStaticSize<t_element, t_domain>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto rhs = getSymmetricGradientRhs(mapping_value.row(), mapping_value.col());
                    auto line = mapping.template block<1, getStaticSize<t_element, t_domain>()>(mapping_value.rank(), 0);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return mapping;
            }

            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getStaticSize<t_element, t_domain>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "Identity")
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getStaticSize<t_element, t_domain>()>();
                mapping.setZero();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto left_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                    auto col_offset = mapping_value.row() * getCellBasisSize<t_element>();
                    auto line = mapping.template block<1, getCellBasisSize<t_element>()>(mapping_value.rank(), col_offset);
                    line = mapping_value.value() * this->template getBasisEvaluation<getCellBasis()>(point);
                }
                return mapping;
            }

            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, 1, getStaticSize<t_element, t_domain>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.getTransformation() == "Identity" && (t_mapping.getRow() >= 0 || t_mapping.getCol() >= 0))
            {
                auto mapping = DenseMatrix<Real, 1, getStaticSize<t_element, t_domain>()>();
                mapping.setZero();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    if (mapping_value.row() == t_mapping.getRow() && mapping_value.col() == t_mapping.getCol())
                    {
                        auto left_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                        auto col_offset = mapping_value.row() * getCellBasisSize<t_element>();
                        auto line = mapping.template block<1, getCellBasisSize<t_element>()>(0, col_offset);
                        line = mapping_value.value() * this->template getBasisEvaluation<getCellBasis()>(point);
                    }
                }
                return mapping;
            }

            
            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getDiscreteFieldDegreeOfFreedomCoefficients()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto offset = 0;
                auto unknown = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
                auto cell_block = unknown.template segment<getNumElementDegreesOfFreedom2<t_element, t_domain>()>(offset);
                cell_block = cell_dof.template getCoefficients<t_field>();
                offset += getNumElementDegreesOfFreedom2<t_element, t_domain>();
                auto set_faces_unknowns = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_dof = face->template getDiscreteField<t_field>().getDegreeOfFreedom();
                        auto face_block = unknown.template segment<getNumElementDegreesOfFreedom2<t_inner_neighbor, t_domain>()>(offset);
                        face_block = face_dof.template getCoefficients<t_field>();
                        offset += getNumElementDegreesOfFreedom2<t_inner_neighbor, t_domain>();
                        // auto face_block = unknown.template segment<face_dof.template getSize<t_inner_neighbor, t_field, getFaceBasis()>()>(offset);
                        // face_block = face_dof.template getCoefficients<t_inner_neighbor, t_field, getFaceBasis()>();
                        // offset += face_dof.template getSize<t_inner_neighbor, t_field, getFaceBasis()>();
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
                return unknown;
            }
            
            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getDiscreteFieldDegreeOfFreedomCoefficients()
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto const & face_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
                return face_dof.template getCoefficients<t_field>();
            }
            
            // void
            // addDiscreteFieldDegreeOfFreedom()
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     static_cast<Base *>(this)->template addDiscreteFieldDegreeOfFreedom<t_field, getCellBasis()>();
            // }

            // // template<Strategy t_s>
            // void
            // addDiscreteFieldDegreeOfFreedom(
            //     // std::unique_ptr<LinearSystem<t_s>> const & linear_system
            // )
            // requires(t_element.isSub(t_domain, 1))
            // {
            //     // static_cast<Base *>(this)->template addDiscreteFieldDegreeOfFreedom<t_field, getFaceBasis()>(linear_system);
            //     static_cast<Base *>(this)->template addDiscreteFieldDegreeOfFreedom<t_field, getFaceBasis()>();
            // }

            // template<Label t_label>
            // void
            // addDiscreteFieldOperator()
            // requires(t_label == "Stabilization")
            // {
            //     this->template getDiscreteField<t_field>().addMatrix(t_label, getStabilization());
            // }
            
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            {}
            
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr num_cell_unknowns = getNumCellUnknowns<t_element, t_domain>();
                auto constexpr strain_operator_num_cols = getStaticSize<t_element, t_domain>();
                auto constexpr num_face_unknowns = getStaticSize<t_element, t_domain>() - num_cell_unknowns;
                auto faces_correction = DenseVector<Real, num_face_unknowns>();
                faces_correction.setZero();
                auto faces_correction_offset = 0;
                auto set_faces_increment = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_dof = face->template getDiscreteField<t_field>().getDegreeOfFreedom();
                        auto lhs = faces_correction.template segment<face_dof.template getSize<t_field, getFaceBasis()>()>(faces_correction_offset);
                        auto rhs = increment.template segment<face_dof.template getSize<t_field, getFaceBasis()>()>(face_dof.getOffset());
                        lhs = rhs;
                        faces_correction_offset += face_dof.template getSize<t_field, getFaceBasis()>();
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_increment(set_faces_increment);
                auto k_tt_inv = this->template getDiscreteField<t_field>().getMatrix("KTT");
                auto k_tf = this->template getDiscreteField<t_field>().getMatrix("KTF");
                auto r_t = this->template getDiscreteField<t_field>().getVector("RT");
                auto cell_corr = k_tt_inv * (r_t - k_tf * faces_correction);
                this->template getDiscreteField<t_field>().getDegreeOfFreedom().addCoefficients(cell_corr);
            }
            
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            requires(t_element.isSub(t_domain, 1))
            {
                static_cast<Base *>(this)->template upgradeDiscreteFieldDegreeOfFreedom<t_field, getFaceBasis()>(increment);
            }
            
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            {
                return 0.0;
            }
            
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                return static_cast<Base const *>(this)->template getDiscreteFieldDegreeOfFreedomValue<t_field, getCellBasis()>(point, row, col);
            }
            
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                Point const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                return static_cast<Base const *>(this)->template getDiscreteFieldDegreeOfFreedomValue<t_field, getFaceBasis()>(point, row, col);
            }

            template<Quadrature t_quadrature>
            Real
            getDiscreteFieldDegreeOfFreedomIntegrationValue(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                this->template getDiscreteFieldDegreeOfFreedomValue<t_field>(row, col);
                auto value = Real(0);
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    value += weight * this->template getDiscreteFieldDegreeOfFreedomValue<t_field>(point, row, col);
                }
                return value;
            }

            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getFieldDualVector(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto vector = this->template getBasisEvaluation<getCellBasis()>(point);
                auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
                external_forces.setZero();
                external_forces.template segment<getCellBasisSize<t_element>()>(offset) = vector;
                return external_forces;
            }

            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getFieldDualVector(
                PointConcept auto const & point
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto vector = this->template getBasisEvaluation<getCellBasis()>(point);
                external_forces.setZero();
                for (auto row = 0; row < FieldTraits<t_field>::getRows(); row++)
                {
                    for (auto col = 0; col < FieldTraits<t_field>::getCols(); col++)
                    {
                        auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
                        external_forces.template segment<getCellBasisSize<t_element>()>(offset) = vector;
                    }
                }
                return external_forces;
            }

            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getFieldDualVector(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto vector = this->template getBasisEvaluation<getFaceBasis()>(point);
                auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
                external_forces.setZero();
                external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = vector;
                return external_forces;
            }

            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getFieldDualVector(
                PointConcept auto const & point
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto vector = this->template getBasisEvaluation<getFaceBasis()>(point);
                external_forces.setZero();
                for (auto row = 0; row < FieldTraits<t_field>::getRows(); row++)
                {
                    for (auto col = 0; col < FieldTraits<t_field>::getCols(); col++)
                    {
                        auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
                        external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = vector;
                    }
                }
                return external_forces;
            }

            //
            //
            //

            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getFieldPrimalVector(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
                auto field_coefficients = cell_dof.template getCoefficients<t_element, t_field, getCellBasis()>(row, col);
                auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
                external_forces.setZero();
                external_forces.template segment<getCellBasisSize<t_element>()>(offset) = field_coefficients;
                return external_forces;
            }

            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getFieldPrimalVector(
                PointConcept auto const & point
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
                external_forces.setZero();
                for (auto row = 0; row < FieldTraits<t_field>::getRows(); row++)
                {
                    for (auto col = 0; col < FieldTraits<t_field>::getCols(); col++)
                    {
                        auto field_coefficients = cell_dof.template getCoefficients<t_element, t_field, getCellBasis()>(row, col);
                        auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
                        external_forces.template segment<getCellBasisSize<t_element>()>(offset) = field_coefficients;
                    }
                }
                return external_forces;
            }

            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getFieldPrimalVector(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
                auto field_coefficients = cell_dof.template getCoefficients<t_element, t_field, getFaceBasis()>(row, col);
                auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
                external_forces.setZero();
                external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = field_coefficients;
                return external_forces;
            }

            DenseVector<Real, getStaticSize<t_element, t_domain>()>
            getFieldPrimalVector(
                PointConcept auto const & point
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
                auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
                external_forces.setZero();
                for (auto row = 0; row < FieldTraits<t_field>::getRows(); row++)
                {
                    for (auto col = 0; col < FieldTraits<t_field>::getCols(); col++)
                    {
                        auto field_coefficients = cell_dof.template getCoefficients<t_element, t_field, getFaceBasis()>(row, col);
                        auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::getCols() * row + col);
                        external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = field_coefficients;
                    }
                }
                return external_forces;
            }

            // auto offset = 0;
            // auto unknown = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
            // auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
            // auto cell_block = unknown.template segment<cell_dof.template getSize<t_element, t_field, getCellBasis()>()>(offset);
            // cell_block = cell_dof.template getCoefficients<t_element, t_field, getCellBasis()>();
            // offset += cell_dof.template getSize<t_element, t_field, getCellBasis()>();
            // auto set_faces_unknowns = [&] <Integer t_i = 0> (
            //     auto & self
            // )
            // constexpr mutable
            // {
            //     auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
            //     for (auto const & face : this->template getInnerNeighbors<0, t_i>())
            //     {
            //         auto const & face_dof = face->template getDiscreteField<t_field>().getDegreeOfFreedom();
            //         auto face_block = unknown.template segment<face_dof.template getSize<t_inner_neighbor, t_field, getFaceBasis()>()>(offset);
            //         face_block = face_dof.template getCoefficients<t_inner_neighbor, t_field, getFaceBasis()>();
            //         offset += face_dof.template getSize<t_inner_neighbor, t_field, getFaceBasis()>();
            //     }
            //     if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
            //     {
            //         self.template operator ()<t_i + 1>(self);
            //     }
            // };
            // set_faces_unknowns(set_faces_unknowns);
            // return unknown;

            //
            //
            //

            // DenseVector<Real, getStaticSize<t_element, t_domain>()>
            // getElementExternalForces(
            //     Real const & time
            // )
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto constexpr quadrature = Quadrature("Gauss", 4);
            //     auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
            //     external_forces.setZero();
            //     for (auto const & domain : this->getDomains())
            //     {
            //         if (domain->hasLoads())
            //         {
            //             for (auto const & load : domain->getLoads())
            //             {
            //                 for (auto i = 0; i < QuadratureTraits<quadrature>::template Rule<t_element>::getSize(); i++)
            //                 {
            //                     auto weight = this->template getCurrentQuadratureWeight<quadrature>(i);
            //                     auto reference_point = this->template getReferenceQuadraturePoint<quadrature>(i);
            //                     auto current_point = this->template getCurrentQuadraturePoint<quadrature>(i);
            //                     auto vector = this->template getBasisEvaluation<getCellBasis()>(reference_point);
            //                     auto offset = getCellBasisSize<t_element>() * FieldTraits<t_field>::getCols() * load.getRow() + load.getCol();
            //                     auto lhs = external_forces.template segment<getCellBasisSize<t_element>()>(offset);
            //                     lhs += weight * load.getValue(current_point, time) * vector;
            //                 }
            //             }
            //         }
            //     }
            //     return external_forces;
            // }

            // DenseVector<Real, getStaticSize<t_element, t_domain>()>
            // getElementExternalForces(
            //     Real const & time
            // )
            // const
            // requires(t_element.isSub(t_domain, 1))
            // {
            //     auto constexpr quadrature = Quadrature("Gauss", 4);
            //     auto external_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
            //     external_forces.setZero();
            //     for (auto const & domain : this->getDomains())
            //     {
            //         if (domain->hasLoads())
            //         {
            //             for (auto const & load : domain->getLoads())
            //             {
            //                 for (auto i = 0; i < QuadratureTraits<quadrature>::template Rule<t_element>::getSize(); i++)
            //                 {
            //                     auto weight = this->template getCurrentQuadratureWeight<quadrature>(i);
            //                     auto reference_point = this->template getReferenceQuadraturePoint<quadrature>(i);
            //                     auto current_point = this->template getCurrentQuadraturePoint<quadrature>(i);
            //                     auto vector = this->template getBasisEvaluation<getFaceBasis()>(reference_point);
            //                     auto offset = getFaceBasisSize<t_element>() * FieldTraits<t_field>::getCols() * load.getRow() + load.getCol();
            //                     auto lhs = external_forces.template segment<getFaceBasisSize<t_element>()>(offset);
            //                     lhs += weight * load.getValue(current_point, time) * vector;
            //                 }
            //             }
            //         }
            //     }
            //     return external_forces;
            // }

            // template<PotentialConcept auto... t_behaviors>
            // DenseVector<Real, getStaticSize<t_element, t_domain>()>
            // getInternalForces()
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto internal_forces = DenseVector<Real, getStaticSize<t_element, t_domain>()>();
            //     internal_forces.setZero();
            //     auto bb = [&] <PotentialConcept auto t_behavior> ()
            //     {
            //         auto aaa = [&] <Integer t_i = 0> (
            //             auto & self
            //         )
            //         constexpr mutable
            //         {
            //             auto constexpr strain = t_behavior.template getStrain<t_i>();
            //             if (strain.getField() == t_field)
            //             {
            //                 auto constexpr strain_operator_num_rows = MappingTraits<strain>::template getSize<t_domain>();
            //                 using StrainView = algebra::View<DenseVector<Real, strain_operator_num_rows>>;
            //                 for (auto const & ip : this->template getFormulation<t_behavior>().getIntegrationPoints())
            //                 {
            //                     auto stress_view = StrainView(ip.behavior_data_->s1.thermodynamic_forces.data());
            //                     internal_forces += ip.getCurrentWeight() * ip.template getStrainOperator<strain>().transpose() * stress_view;
            //                 }
            //                 if constexpr (t_i < t_behavior.getNumMappings() - 1)
            //                 {
            //                     self.template operator()<t_i + 1>(self);
            //                 }
            //             }
            //         };
            //         aaa(aaa);
            //     };
            //     bb.template operator()<t_behaviors...>();
            //     return internal_forces;
            // }

            // DenseVector<Real, getStaticSize<t_element, t_domain>()>
            // getRhs(
            //     Real const & time
            // )
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     return this->getExternalForces(time) - this->getInternalForces();
            // }

            // template<PotentialConcept auto... t_behaviors>
            // DenseMatrix<Real, getStaticSize<t_element, t_domain>(), getStaticSize<t_element, t_domain>()>
            // getLhs()
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto jacobian_matrix = DenseMatrix<Real, getStaticSize<t_element, t_domain>(), getStaticSize<t_element, t_domain>()>();
            //     jacobian_matrix.setZero();
            //     // auto bb = [&] <PotentialConcept auto t_behavior> ()
            //     // {
            //     //     auto aaa = [&] <Integer t_i = 0> (
            //     //         auto & self
            //     //     )
            //     //     constexpr mutable
            //     //     {
            //     //         auto constexpr strain = t_behavior.template getStrain<t_i>();
            //     //         auto constexpr strain_operator_num_rows = MappingTraits<strain>::template getSize<t_domain>();
            //     //         using StrainView = algebra::View<DenseVector<Real, strain_operator_num_rows>>;
            //     //         for (auto const & ip : this->template getFormulation<t_behavior>().getIntegrationPoints())
            //     //         {
            //     //             auto stress_view = StrainView(ip.behavior_data_->s1.thermodynamic_forces.data());
            //     //             internal_forces += ip.getCurrentWeight() * ip.template getStrainOperator<strain>().transpose() * stress_view;
            //     //         }
            //     //         if constexpr (t_i < t_behavior.getNumMappings() - 1)
            //     //         {
            //     //             self.template operator()<t_i + 1>(self);
            //     //         }
            //     //     };
            //     //     aaa(aaa);
            //     // };
            //     // bb.template operator()<t_behaviors...>();
            //     return jacobian_matrix;
            // }

            // DenseVector<Real, getStaticSize<t_element, t_domain>()>
            // getRhs(
            //     Real const & time
            // )
            // const
            // requires(t_element.isSub(t_domain, 1))
            // {
            //     auto set_faces_unknowns = [&] <Integer t_i = 0> (
            //         auto & self
            //     )
            //     constexpr mutable
            //     {
            //         for (auto const & cell : this->template getOuterNeighbors<1, t_i>())
            //         {
            //             auto cell_rhs = cell->getExternalForces(time) - cell->getInternalForces();
            //             auto constexpr t_outer_neighbor = ElementTraits<t_element>::template getOuterNeighbor<t_domain, 1, t_i>();
            //             auto set_faces_unknowns2 = [&] <Integer t_j = 0> (
            //                 auto & self2
            //             )
            //             constexpr mutable
            //             {
            //                 for (auto const & face : cell->template getInnerNeighbors<0, t_j>())
            //                 {
            //                     // auto face_unknown_size = face->degrees_of_freedom_.at(std::string(unknown_label)).template getSize<field, getFaceBasis()>();
            //                     // band_width += face_unknown_size;
            //                 }
            //                 if constexpr (t_j < ElementTraits<t_outer_neighbor>::template getNumInnerNeighbors<0>() - 1)
            //                 {
            //                     self2.template operator ()<t_j + 1>(self2);
            //                 }
            //             };
            //             set_faces_unknowns2(set_faces_unknowns2);
            //         }
            //         if constexpr (t_i < ElementTraits<t_element>::template getNumOuterNeighbors<t_domain, 1>() - 1)
            //         {
            //             self.template operator ()<t_i + 1>(self);
            //         }
            //     };
            //     set_faces_unknowns(set_faces_unknowns);
            //     return this->getExternalForces(time) - this->getInternalForces();
            // }

            /**
             * Implementation
             * *************************************************************************************************************************************************
             */

            static constexpr
            Integer
            getGradientConstructionQuadratureOrder()
            {
                if constexpr (t_domain.isAxiSymmetric())
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1, getGradBasis().getOrd()) + 1;
                }
                else
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1, getGradBasis().getOrd());
                }
            }

            static constexpr
            Integer
            getPotentialConstructionQuadratureOrder()
            {
                if constexpr (t_domain.isAxiSymmetric())
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1) + 1;
                }
                else
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1);
                }
            }
            
            DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>
            getPotentialLhs()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto lhs = DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>();
                lhs.setZero();
                for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
                {
                    for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                        auto vector = this->template getBasisDerivative<getPotentialBasis()>(point, i_component);
                        auto vector_j = vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        lhs += weight * vector_j * vector_j.transpose();
                    }
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }
            
            DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getStaticSize<t_element, t_domain>()>
            getPotentialRhs(
                Integer row
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto rhs = DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getStaticSize<t_element, t_domain>()>();
                rhs.setZero();
                for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
                {
                    auto face_offset = t_field_size * getCellBasisSize<t_element>();
                    for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                        auto row_cell_vector = this->template getBasisDerivative<getPotentialBasis()>(point, i_component);
                        auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        auto col_cell_vector = this->template getBasisDerivative<getCellBasis()>(point, i_component);
                        auto col_offset = row * getCellBasisSize<t_element>();
                        auto block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getCellBasisSize<t_element>()>(0, col_offset);
                        block += weight * row_cell_vector_j * col_cell_vector.transpose();
                    }
                    auto set_faces_blocks = [&] <Integer t_i = 0> (
                        auto & self
                    )
                    constexpr mutable
                    {
                        auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                        auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                        auto i_face = 0;
                        for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                        {
                            auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                            auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                            for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                            {
                                auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                                auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                auto argf = i;
                                //
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                                auto normal_vector = face->getNormalVector(point);
                                //
                                auto row_cell_vector = this->template getBasisDerivative<getPotentialBasis()>(i_p, i_component);
                                auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                                auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                                auto normal_vector_component = face_orientation * normal_vector(i_component);
                                auto face_col_offset = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                                auto face_block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                                face_block += normal_vector_component * weight * row_cell_vector_j * col_face_vector.transpose();
                                //
                                auto cell_col_offset = row * getCellBasisSize<t_element>();
                                auto cell_block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getCellBasisSize<t_element>()>(0, cell_col_offset);
                                cell_block -= normal_vector_component * weight * row_cell_vector_j * col_cell_vector.transpose();
                            }
                            face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                            i_face ++;
                        }
                        if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                        {
                            self.template operator()<t_i + 1>(self);
                        }
                    };
                    set_faces_blocks(set_faces_blocks);
                }
                return rhs;
            }
            
            DenseMatrix<Real, getStaticSize<t_element, t_domain>(), getStaticSize<t_element, t_domain>()>
            getStabilization()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                //
                auto S_T_A_B_I_L_I_Z_A_T_I_O_N = DenseMatrix<Real, getStaticSize<t_element, t_domain>(), getStaticSize<t_element, t_domain>()>();
                auto potential_operator = DenseMatrix<Real, getPotentialBasisSize<t_element>(), getStaticSize<t_element, t_domain>()>();
                auto potential_lhs = getPotentialLhs();
                auto denom = 1.0 / QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize();
                S_T_A_B_I_L_I_Z_A_T_I_O_N.setZero();
                for (auto i_component = 0; i_component < FieldTraits<t_field>::template getSize<t_domain>(); i_component++)
                {
                    auto face_offset = t_field_size * getCellBasisSize<t_element>();
                    potential_operator.setZero();
                    auto oppp = potential_operator.template block<getPotentialBasisSize<t_element>() - 1, getStaticSize<t_element, t_domain>()>(1, 0);
                    oppp = potential_lhs * getPotentialRhs(i_component);
                    for (auto i_quadrature = 0; i_quadrature < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i_quadrature++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                        auto row_cell_vector = this->template getBasisEvaluation<getPotentialBasis()>(point);
                        auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        potential_operator.template block<1, getStaticSize<t_element, t_domain>()>(0, 0) -= denom * weight * row_cell_vector_j.transpose() * oppp;
                        potential_operator.template block<1, getCellBasisSize<t_element>()>(0, i_component * getCellBasisSize<t_element>()) +=
                        denom * weight * this->template getBasisEvaluation<getCellBasis()>(point).transpose();
                    }
                    //
                    auto cell_projection_lhs = DenseMatrix<Real, getCellBasisSize<t_element>(), getCellBasisSize<t_element>()>();
                    auto cell_projection_rhs = DenseMatrix<Real, getCellBasisSize<t_element>(), getStaticSize<t_element, t_domain>()>();
                    cell_projection_lhs.setZero();
                    cell_projection_rhs.setZero();
                    for (auto i_quadrature = 0; i_quadrature < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i_quadrature++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                        auto row_cell_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                        cell_projection_lhs += weight * row_cell_vector * row_cell_vector.transpose();
                        cell_projection_rhs += weight * row_cell_vector * this->template getBasisEvaluation<getPotentialBasis()>(point).transpose() * potential_operator;
                    }
                    auto proj = cell_projection_lhs.llt().solve(decltype(cell_projection_lhs)::Identity()) * cell_projection_rhs;
                    //
                    auto set_faces_blocks = [&] <Integer t_i = 0> (
                        auto & self
                    )
                    constexpr mutable
                    {
                        auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                        auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                        auto i_face = 0;
                        for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                        {
                            auto constexpr quadrature_size = QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize();
                            auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                            auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                            auto face_projection_lhs = DenseMatrix<Real, getFaceBasisSize<t_inner_neighbor>(), getFaceBasisSize<t_inner_neighbor>()>();
                            auto face_projection_rhs = DenseMatrix<Real, getFaceBasisSize<t_inner_neighbor>(), getStaticSize<t_element, t_domain>()>();
                            face_projection_lhs.setZero();
                            face_projection_rhs.setZero();
                            for (auto i_quadrature = 0; i_quadrature < quadrature_size; i_quadrature++)
                            {
                                auto argh = face_sign == 1 ? i_quadrature : quadrature_size - (i_quadrature + 1);
                                auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                auto argf = i_quadrature;
                                //
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                                auto normal_vector = face->getNormalVector(point);
                                auto row_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto col_potential_vector = this->template getBasisEvaluation<getPotentialBasis()>(i_p);
                                auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                                //
                                face_projection_lhs += weight * row_face_vector * row_face_vector.transpose();
                                face_projection_rhs -= weight * row_face_vector * col_potential_vector.transpose() * potential_operator;
                                face_projection_rhs += weight * row_face_vector * col_cell_vector.transpose() * proj;
                                //
                                auto cell_col_offset = i_component * getCellBasisSize<t_element>();
                                auto cell_block = face_projection_rhs.template block<getFaceBasisSize<t_inner_neighbor>(), getCellBasisSize<t_element>()>(0, cell_col_offset);
                                cell_block -= weight * row_face_vector * col_cell_vector.transpose();
                                //
                                auto face_col_offset = face_offset + i_component * getFaceBasisSize<t_inner_neighbor>();
                                auto face_block = face_projection_rhs.template block<getFaceBasisSize<t_inner_neighbor>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                                face_block += weight * row_face_vector * row_face_vector.transpose();
                            }
                            auto fa_ce_proj = face_projection_lhs.llt().solve(decltype(face_projection_lhs)::Identity()) * face_projection_rhs;
                            for (auto i_quadrature = 0; i_quadrature < quadrature_size; i_quadrature++)
                            {
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                                auto row_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto V_E_C_T_O_R = row_face_vector.transpose() * fa_ce_proj;
                                auto factor = (1.0 / face->getLocalFrameDiameters().norm());
                                S_T_A_B_I_L_I_Z_A_T_I_O_N += factor * weight * V_E_C_T_O_R.transpose() * V_E_C_T_O_R;
                            }
                            face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                            i_face ++;
                        }
                        if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                        {
                            self.template operator()<t_i + 1>(self);
                        }
                    };
                    set_faces_blocks(set_faces_blocks);
                }
                return S_T_A_B_I_L_I_Z_A_T_I_O_N;
            }
            
            DenseMatrix<Real, getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>
            getGradientLhs()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto lhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>();
                lhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    lhs += weight * vector * vector.transpose();
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }
            
            DenseMatrix<Real, getGradBasisSize<t_element>(), getStaticSize<t_element, t_domain>()>
            getGradientRhs(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getStaticSize<t_element, t_domain>()>();
                rhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    auto col_cell_vector = this->template getBasisDerivative<getCellBasis()>(point, col);
                    auto col_offset = row * getCellBasisSize<t_element>();
                    auto block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset);
                    block += weight * row_cell_vector * col_cell_vector.transpose();
                }
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                            auto argf = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            //
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            auto normal_vector = face->getNormalVector(point);
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            auto normal_vector_component = face_orientation * normal_vector(col);
                            auto face_col_offset = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                            face_block += normal_vector_component * weight * row_cell_vector * col_face_vector.transpose();
                            auto cell_col_offset = row * getCellBasisSize<t_element>();
                            auto cell_block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset);
                            cell_block -= normal_vector_component * weight * row_cell_vector * col_cell_vector.transpose();
                        }
                        face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                        i_face ++;
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator()<t_i + 1>(self);
                    }
                };
                set_faces_blocks(set_faces_blocks);
                return rhs;
            }
            
            DenseMatrix<Real, getGradBasisSize<t_element>(), getStaticSize<t_element, t_domain>()>
            getSymmetricGradientRhs(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getStaticSize<t_element, t_domain>()>();
                rhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    auto col_cell_vector_i = this->template getBasisDerivative<getCellBasis()>(point, col);
                    auto col_cell_vector_j = this->template getBasisDerivative<getCellBasis()>(point, row);
                    auto col_offset_i = row * getCellBasisSize<t_element>();
                    auto block_i = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset_i);
                    block_i += (1.0 / 2.0) * weight * row_cell_vector * col_cell_vector_i.transpose();
                    auto col_offset_j = col * getCellBasisSize<t_element>();
                    auto block_j = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset_j);
                    block_j += (1.0 / 2.0) * weight * row_cell_vector * col_cell_vector_j.transpose();
                }
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                            auto argf = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            auto normal_vector = face->getNormalVector(point);
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            auto normal_vector_component_i = face_orientation * normal_vector(col);
                            auto normal_vector_component_j = face_orientation * normal_vector(row);
                            auto face_col_offset_i = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block_i = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset_i);
                            face_block_i += (1.0 / 2.0) * normal_vector_component_i * weight * row_cell_vector * col_face_vector.transpose();
                            //
                            auto face_col_offset_j = face_offset + col * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block_j = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset_j);
                            face_block_j += (1.0 / 2.0) * normal_vector_component_j * weight * row_cell_vector * col_face_vector.transpose();
                            //
                            auto cell_col_offset_i = row * getCellBasisSize<t_element>();
                            auto cell_block_i = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset_i);
                            cell_block_i -= (1.0 / 2.0) * normal_vector_component_i * weight * row_cell_vector * col_cell_vector.transpose();
                            //
                            auto cell_col_offset_j = col * getCellBasisSize<t_element>();
                            auto cell_block_j = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset_j);
                            cell_block_j -= (1.0 / 2.0) * normal_vector_component_j * weight * row_cell_vector * col_cell_vector.transpose();
                        }
                        face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                        i_face ++;
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator()<t_i + 1>(self);
                    }
                };
                set_faces_blocks(set_faces_blocks);
                return rhs;
            }

            template<PotentialConcept auto t_behavior>
            void
            hhh()
            const
            {

            }

            // template<PotentialConcept auto t_behavior, MappingConcept auto t_strain>
            // DenseVector<Real, getSize()>
            // getInternalForces()
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto constexpr strain_operator_num_rows = MappingTraits<t_strain>::template getSize<t_domain>();
            //     auto constexpr strain_operator_num_cols = getSize();
            //     // using StrainOperatorView = algebra::View<DenseMatrix<Real, strain_operator_num_rows, strain_operator_num_cols> const>;
            //     using StrainView = algebra::View<DenseVector<Real, strain_operator_num_rows>>;
            //     auto const unknown = this->template getDiscreteFieldDegreeOfFreedomCoefficients<t_field>();
            //     auto internal_forces = DenseVector<Real, strain_operator_num_cols>();
            //     internal_forces.setZero();
            //     for (auto const & ip : this->template getFormulation<t_behavior>().getIntegrationPoints())
            //     {
            //         // auto strain_operator_view = StrainOperatorView(ip.template getStrainOperator<t_strain>().data());
            //         auto stress_view = StrainView(ip.behavior_data_->s1.thermodynamic_forces.data());
            //         internal_forces += ip.getCurrentWeight() * ip.template getStrainOperator<t_strain>().transpose() * stress_view;
            //     }
            //     auto const & s_o = this->template getDiscreteField<t_field>().getMatrix("Stabilization");
            //     auto const & s_p = this->template getDiscreteField<t_field>().getScalar("Stabilization");
            //     internal_forces += s_p * s_o * unknown;
            //     return internal_forces;
            // }

            // template<PotentialConcept auto t_behavior, MappingConcept auto t_strain>
            // void
            // setInternalForces(
            //     DenseVector<Real> & internal_forces,
            //     Integer & offset
            // )
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto constexpr t_field = t_strain.getField();
            //     auto constexpr strain_operator_num_rows = MappingTraits<t_strain>::template getSize<t_domain>();
            //     auto constexpr strain_operator_num_cols = getSize();
            //     using StrainView = algebra::View<DenseVector<Real, strain_operator_num_rows>>;
            //     auto internal_forces_block = internal_forces.template segment<getSize()>(offset);
            //     auto const unknown = this->template getDiscreteFieldDegreeOfFreedomCoefficients<t_field>();
            //     for (auto const & ip : this->template getFormulation<t_behavior>().getIntegrationPoints())
            //     {
            //         auto stress_view = StrainView(ip.behavior_data_->s1.thermodynamic_forces.data());
            //         internal_forces_block += ip.getCurrentWeight() * ip.template getStrainOperator<t_strain>().transpose() * stress_view;
            //     }
            //     offset += getSize();
            // }

            // template<PotentialConcept auto t_behavior, MappingConcept auto t_strain>
            // void
            // setStabilizationForces(
            //     DenseVector<Real> & internal_forces,
            //     Integer & offset
            // )
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto constexpr t_field = t_strain.getField();
            //     auto constexpr strain_operator_num_cols = getSize();
            //     auto internal_forces_block = internal_forces.template segment<strain_operator_num_cols>(offset);
            //     auto const unknown = this->template getDiscreteFieldDegreeOfFreedomCoefficients<t_field>();
            //     auto const & s_o = this->template getDiscreteField<t_field>().getMatrix("Stabilization");
            //     auto const & s_p = this->template getDiscreteField<t_field>().getScalar("Stabilization");
            //     internal_forces_block += s_p * s_o * unknown;
            // }

        };
        
    };

} // namespace lolita

#endif /* FEB10822_D09A_4468_9437_2ADB758EE6F4 */
