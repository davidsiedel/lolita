#ifndef A50FB6B0_0CE6_4A1A_BC52_C8CF6B7F8812
#define A50FB6B0_0CE6_4A1A_BC52_C8CF6B7F8812

#include "2/core/_include.hxx"
#include "2/core/region.hxx"
// #include "2/core/frm.hxx"
#include "2/core/dof.hxx"

#include "2/core/element_integration_point.hxx"
#include "2/core/element_potential.hxx"
#include "2/core/element_lagrangian_implementation.hxx"
#include "2/core/element_lagrangian_interface.hxx"

#include "2/core/region_potential.hxx"
#include "2/core/region_lagrangian_implementation.hxx"
#include "2/core/region_lagrangian_interface.hxx"

namespace lolita::core
{
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    struct FiniteElement;

    template<LagrangeShapeConcept auto t_element, MeshConcept auto t_domain, FieldConcept auto t_field>
    struct DegreesOfFreedomImplementation;

    template<LagrangeShapeConcept auto t_element, MeshConcept auto t_domain>
    struct DegreesOfFreedomInterface;

    template<ShapeConcept auto t_element, MeshConcept auto t_domain, FieldConcept auto t_field>
    struct DiscretizationImplementation;

    template<ShapeConcept auto t_element, MeshConcept auto t_domain, FieldConcept auto t_field>
    requires(HybridDiscontinuousGalerkinDiscretizationConcept<decltype(t_field.getDiscretization())>)
    struct DiscretizationImplementation<t_element, t_domain, t_field> : DegreesOfFreedomInterface<t_element, t_domain>
    {

    private:

        static constexpr
        BasisConcept auto const &
        getCellBasis()
        {
            return t_field.getDiscretization().getCellBasis();
        }

        static constexpr
        BasisConcept auto const &
        getFaceBasis()
        {
            return t_field.getDiscretization().getFaceBasis();
        }

        static constexpr
        BasisConcept auto const &
        getGradBasis()
        {
            return t_field.getDiscretization().getFaceBasis();
        }

        static constexpr
        BasisConcept auto
        getPotentialBasis()
        {
            return t_field.getDiscretization().getFaceBasis().toOrder(t_field.getDiscretization().getFaceBasis().getOrd() + 1);
        }

        template<LagrangeShapeConcept auto shape_>
        static constexpr
        Integer
        getCellBasisSize()
        {
            return BasisTraits<getCellBasis()>::template getSize<shape_>();
        }

        template<LagrangeShapeConcept auto shape_>
        static constexpr
        Integer
        getFaceBasisSize()
        {
            return BasisTraits<getFaceBasis()>::template getSize<shape_>();
        }

        template<LagrangeShapeConcept auto shape_>
        static constexpr
        Integer
        getGradBasisSize()
        {
            return BasisTraits<getGradBasis()>::template getSize<shape_>();
        }

        template<LagrangeShapeConcept auto shape_>
        static constexpr
        Integer
        getPotentialBasisSize()
        {
            return BasisTraits<getPotentialBasis()>::template getSize<shape_>();
        }


        static constexpr
        Integer
        getTensorSize()
        {
            return Traits_::template getTensorSpaceSize<t_element, t_domain>();
        }

        using Traits_ = DiscretizationTraits<t_field>;

        using Base_ = FiniteElement<t_element, t_domain>;

    public:
            
        // template<StabilizationOperatorConcept auto t_mapping>
        // DenseMatrix<Real, getTensorSize(), getTensorSize()>
        // getMapping(
        //     PointConcept auto const & point
        // )
        // const
        // // requires(t_mapping.getTransformation() == "Stabilization")
        // {
        //     return getStabilization();
        // }
        
        template<LinearOperatorConcept auto t_mapping>
        DenseMatrix<Real, LinearOperatorTraits<t_mapping>::template getSize<t_domain, t_field>(), getTensorSize()>
        letLinearOperator(
            PointConcept auto const & point
        )
        const
        requires(GradientOperatorConcept<decltype(t_mapping)>)
        {
            auto mapping = DenseMatrix<Real, LinearOperatorTraits<t_mapping>::template getSize<t_domain, t_field>(), getTensorSize()>();
            mapping.setZero();
            auto lhs = this->getGradientLhs();
            for (auto const & mapping_value : LinearOperatorTraits<t_mapping>::template getValues<t_domain, t_field>())
            {
                auto rhs = this->getGradientRhs(mapping_value.row(), mapping_value.col());
                auto line = mapping.template block<1, getTensorSize()>(mapping_value.rank(), 0);
                line = mapping_value.value() * this->getFiniteElement().template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
            }
            return mapping;
        }
        
        template<LinearOperatorConcept auto t_mapping>
        DenseMatrix<Real, LinearOperatorTraits<t_mapping>::template getSize<t_domain, t_field>(), getTensorSize()>
        letLinearOperator(
            PointConcept auto const & point
        )
        const
        requires(LargeStrainOperatorConcept<decltype(t_mapping)>)
        {
            auto mapping = DenseMatrix<Real, LinearOperatorTraits<t_mapping>::template getSize<t_domain, t_field>(), getTensorSize()>();
            mapping.setZero();
            auto lhs = this->getGradientLhs();
            for (auto const & mapping_value : LinearOperatorTraits<t_mapping>::template getValues<t_domain, t_field>())
            {
                auto rhs = this->getGradientRhs(mapping_value.row(), mapping_value.col());
                auto line = mapping.template block<1, getTensorSize()>(mapping_value.rank(), 0);
                line = mapping_value.value() * this->getFiniteElement().template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
            }
            return mapping;
        }

        template<LinearOperatorConcept auto t_mapping>
        DenseMatrix<Real, LinearOperatorTraits<t_mapping>::template getSize<t_domain, t_field>(), getTensorSize()>
        letLinearOperator(
            PointConcept auto const & point
        )
        const
        requires(SmallStrainOperatorConcept<decltype(t_mapping)>)
        {
            auto mapping = DenseMatrix<Real, LinearOperatorTraits<t_mapping>::template getSize<t_domain, t_field>(), getTensorSize()>();
            mapping.setZero();
            auto lhs = this->getGradientLhs();
            for (auto const & mapping_value : LinearOperatorTraits<t_mapping>::template getValues<t_domain, t_field>())
            {
                auto rhs = this->getSymmetricGradientRhs(mapping_value.row(), mapping_value.col());
                auto line = mapping.template block<1, getTensorSize()>(mapping_value.rank(), 0);
                line = mapping_value.value() * this->getFiniteElement().template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
            }
            return mapping;
        }

        template<LinearOperatorConcept auto t_mapping>
        DenseMatrix<Real, LinearOperatorTraits<t_mapping>::template getSize<t_domain, t_field>(), getTensorSize()>
        letLinearOperator(
            PointConcept auto const & point
        )
        const
        requires(TraceOperatorConcept<decltype(t_mapping)>)
        {
            auto mapping = DenseMatrix<Real, LinearOperatorTraits<t_mapping>::template getSize<t_domain, t_field>(), getTensorSize()>();
            mapping.setZero();
            for (auto const & mapping_value : LinearOperatorTraits<t_mapping>::template getValues<t_domain, t_field>())
            {
                auto left_vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(point);
                auto col_offset = mapping_value.row() * getCellBasisSize<t_element>();
                auto line = mapping.template block<1, getCellBasisSize<t_element>()>(mapping_value.rank(), col_offset);
                line = mapping_value.value() * this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(point);
            }
            return mapping;
        }

        template<LinearOperatorConcept auto t_mapping>
        DenseMatrix<Real, 1, getTensorSize()>
        letLinearOperator(
            PointConcept auto const & point
        )
        const
        requires(t_mapping.getTransformation() == "Identity" && (t_mapping.getRow() >= 0 || t_mapping.getCol() >= 0))
        {
            auto mapping = DenseMatrix<Real, 1, getTensorSize()>();
            mapping.setZero();
            for (auto const & mapping_value : LinearOperatorTraits<t_mapping>::template getValues<t_domain, t_field>())
            {
                if (mapping_value.row() == t_mapping.getRow() && mapping_value.col() == t_mapping.getCol())
                {
                    auto left_vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(point);
                    auto col_offset = mapping_value.row() * getCellBasisSize<t_element>();
                    auto line = mapping.template block<1, getCellBasisSize<t_element>()>(0, col_offset);
                    line = mapping_value.value() * this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(point);
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
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            auto constexpr num_cell_unknowns = Traits_::template getTensorBasisSize<t_element, t_domain>();
            auto constexpr strain_operator_num_cols = getTensorSize();
            auto constexpr num_face_unknowns = getTensorSize() - num_cell_unknowns;
            auto faces_correction = DenseVector<Real, num_face_unknowns>();
            faces_correction.setZero();
            auto faces_correction_offset = 0;
            auto set_faces_increment = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                for (auto const & face : this->getFiniteElement().template getInnerNeighbors<0, t_i>())
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
        requires(t_element.getDim() == t_field.getDimDomain() - 1)
        {
            static_cast<Base_ *>(this)->template upgradeDiscreteFieldDegreeOfFreedom<t_field, getFaceBasis()>(increment);
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
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            return static_cast<Base_ const *>(this)->template getDiscreteFieldDegreeOfFreedomValue<t_field, getCellBasis()>(point, row, col);
        }
        
        Real
        getDiscreteFieldDegreeOfFreedomValue(
            Point const & point,
            Integer row,
            Integer col
        )
        const
        requires(t_element.getDim() == t_field.getDimDomain() - 1)
        {
            return static_cast<Base_ const *>(this)->template getDiscreteFieldDegreeOfFreedomValue<t_field, getFaceBasis()>(point, row, col);
        }

        template<QuadratureConcept auto t_quadrature>
        Real
        getDiscreteFieldDegreeOfFreedomIntegrationValue(
            Integer row,
            Integer col
        )
        const
        requires(t_element.getDim() == t_field.getDimDomain() - 1)
        {
            // this->template getDiscreteFieldDegreeOfFreedomValue<t_field>(row, col);
            auto value = Real(0);
            for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
            {
                auto point = this->getFiniteElement().template getReferenceQuadraturePoint<t_quadrature>(i);
                auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<t_quadrature>(i);
                value += weight * this->template getDiscreteFieldDegreeOfFreedomValue<t_field>(point, row, col);
            }
            return value;
        }

        DenseVector<Real, getTensorSize()>
        getFieldDualVector(
            PointConcept auto const & point,
            Integer row,
            Integer col
        )
        const
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            auto external_forces = DenseVector<Real, getTensorSize()>();
            auto vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(point);
            auto offset = getCellBasisSize<t_element>() * (FieldTraits<t_field>::template getCols<t_domain>() * row + col);
            external_forces.setZero();
            external_forces.template segment<getCellBasisSize<t_element>()>(offset) = vector;
            return external_forces;
        }

        DenseVector<Real, getTensorSize()>
        getFieldDualVector(
            PointConcept auto const & point
        )
        const
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            auto external_forces = DenseVector<Real, getTensorSize()>();
            auto vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(point);
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

        DenseVector<Real, getTensorSize()>
        getFieldDualVector(
            PointConcept auto const & point,
            Integer row,
            Integer col
        )
        const
        requires(t_element.getDim() == t_field.getDimDomain() - 1)
        {
            auto external_forces = DenseVector<Real, getTensorSize()>();
            auto vector = this->getFiniteElement().template getBasisEvaluation<getFaceBasis()>(point);
            auto offset = getFaceBasisSize<t_element>() * (FieldTraits<t_field>::template getCols<t_domain>() * row + col);
            external_forces.setZero();
            external_forces.template segment<getFaceBasisSize<t_element>()>(offset) = vector;
            return external_forces;
        }

        /**
         * Implementation
         * *************************************************************************************************************************************************
         */

        static constexpr
        Integer
        getGradientConstructionQuadratureOrder()
        {
            return MeshTraits<t_domain>::letQuadratureOrder(numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1, getGradBasis().getOrd()));
        }

        static constexpr
        Integer
        getPotentialConstructionQuadratureOrder()
        {
            return MeshTraits<t_domain>::letQuadratureOrder(numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1));
        }
        
        DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>
        getPotentialLhs()
        const
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
            auto constexpr t_quadrature = GaussQuadrature(getPotentialConstructionQuadratureOrder());
            auto lhs = DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>();
            lhs.setZero();
            for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
            {
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->getFiniteElement().template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vector = this->getFiniteElement().template getBasisDerivative<getPotentialBasis()>(point, i_component);
                    auto vector_j = vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                    lhs += weight * vector_j * vector_j.transpose();
                }
            }
            return lhs.llt().solve(decltype(lhs)::Identity());
        }
        
        DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getTensorSize()>
        getPotentialRhs(
            Integer row
        )
        const
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
            auto constexpr t_quadrature = GaussQuadrature(getPotentialConstructionQuadratureOrder());
            auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
            auto rhs = DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getTensorSize()>();
            rhs.setZero();
            for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
            {
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->getFiniteElement().template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto row_cell_vector = this->getFiniteElement().template getBasisDerivative<getPotentialBasis()>(point, i_component);
                    auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                    auto col_cell_vector = this->getFiniteElement().template getBasisDerivative<getCellBasis()>(point, i_component);
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
                    for (auto const & face : this->getFiniteElement().template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->getFiniteElement().template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->getFiniteElement().template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                            auto i_p = this->getFiniteElement().template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            auto argf = i;
                            //
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            auto normal_vector = face->getNormalVector(point);
                            //
                            auto row_cell_vector = this->getFiniteElement().template getBasisDerivative<getPotentialBasis()>(i_p, i_component);
                            auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(i_p);
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
        
        DenseMatrix<Real, getTensorSize(), getTensorSize()>
        getStabilization()
        const
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
            auto constexpr t_quadrature = GaussQuadrature(getPotentialConstructionQuadratureOrder());
            auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
            //
            auto S_T_A_B_I_L_I_Z_A_T_I_O_N = DenseMatrix<Real, getTensorSize(), getTensorSize()>();
            auto potential_operator = DenseMatrix<Real, getPotentialBasisSize<t_element>(), getTensorSize()>();
            auto potential_lhs = getPotentialLhs();
            auto denom = 1.0 / QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize();
            S_T_A_B_I_L_I_Z_A_T_I_O_N.setZero();
            for (auto i_component = 0; i_component < FieldTraits<t_field>::template getSize<t_domain>(); i_component++)
            {
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                potential_operator.setZero();
                auto oppp = potential_operator.template block<getPotentialBasisSize<t_element>() - 1, getTensorSize()>(1, 0);
                oppp = potential_lhs * getPotentialRhs(i_component);
                for (auto i_quadrature = 0; i_quadrature < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i_quadrature++)
                {
                    auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                    auto point = this->getFiniteElement().template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                    auto row_cell_vector = this->getFiniteElement().template getBasisEvaluation<getPotentialBasis()>(point);
                    auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                    potential_operator.template block<1, getTensorSize()>(0, 0) -= denom * weight * row_cell_vector_j.transpose() * oppp;
                    potential_operator.template block<1, getCellBasisSize<t_element>()>(0, i_component * getCellBasisSize<t_element>()) +=
                    denom * weight * this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(point).transpose();
                }
                //
                auto cell_projection_lhs = DenseMatrix<Real, getCellBasisSize<t_element>(), getCellBasisSize<t_element>()>();
                auto cell_projection_rhs = DenseMatrix<Real, getCellBasisSize<t_element>(), getTensorSize()>();
                cell_projection_lhs.setZero();
                cell_projection_rhs.setZero();
                for (auto i_quadrature = 0; i_quadrature < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i_quadrature++)
                {
                    auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                    auto point = this->getFiniteElement().template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                    auto row_cell_vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(point);
                    cell_projection_lhs += weight * row_cell_vector * row_cell_vector.transpose();
                    cell_projection_rhs += weight * row_cell_vector * this->getFiniteElement().template getBasisEvaluation<getPotentialBasis()>(point).transpose() * potential_operator;
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
                    for (auto const & face : this->getFiniteElement().template getInnerNeighbors<0, t_i>())
                    {
                        auto constexpr quadrature_size = QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize();
                        auto face_orientation = this->getFiniteElement().template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->getFiniteElement().template getInnerNeighborSign<0, t_i>(i_face);
                        auto face_projection_lhs = DenseMatrix<Real, getFaceBasisSize<t_inner_neighbor>(), getFaceBasisSize<t_inner_neighbor>()>();
                        auto face_projection_rhs = DenseMatrix<Real, getFaceBasisSize<t_inner_neighbor>(), getTensorSize()>();
                        face_projection_lhs.setZero();
                        face_projection_rhs.setZero();
                        for (auto i_quadrature = 0; i_quadrature < quadrature_size; i_quadrature++)
                        {
                            auto argh = face_sign == 1 ? i_quadrature : quadrature_size - (i_quadrature + 1);
                            auto i_p = this->getFiniteElement().template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            auto argf = i_quadrature;
                            //
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            auto normal_vector = face->getNormalVector(point);
                            auto row_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_potential_vector = this->getFiniteElement().template getBasisEvaluation<getPotentialBasis()>(i_p);
                            auto col_cell_vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(i_p);
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
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
            auto constexpr t_quadrature = GaussQuadrature(getGradientConstructionQuadratureOrder());
            auto lhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>();
            lhs.setZero();
            for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
            {
                auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<t_quadrature>(i);
                auto point = this->getFiniteElement().template getReferenceQuadraturePoint<t_quadrature>(i);
                auto vector = this->getFiniteElement().template getBasisEvaluation<getGradBasis()>(point);
                lhs += weight * vector * vector.transpose();
            }
            return lhs.llt().solve(decltype(lhs)::Identity());
        }
        
        DenseMatrix<Real, getGradBasisSize<t_element>(), getTensorSize()>
        getGradientRhs(
            Integer row,
            Integer col
        )
        const
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
            auto constexpr t_quadrature = GaussQuadrature(getGradientConstructionQuadratureOrder());
            auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
            auto face_offset = t_field_size * getCellBasisSize<t_element>();
            auto rhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getTensorSize()>();
            rhs.setZero();
            for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
            {
                auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<t_quadrature>(i);
                auto point = this->getFiniteElement().template getReferenceQuadraturePoint<t_quadrature>(i);
                auto row_cell_vector = this->getFiniteElement().template getBasisEvaluation<getGradBasis()>(point);
                auto col_cell_vector = this->getFiniteElement().template getBasisDerivative<getCellBasis()>(point, col);
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
                for (auto const & face : this->getFiniteElement().template getInnerNeighbors<0, t_i>())
                {
                    auto face_orientation = this->getFiniteElement().template getInnerNeighborOrientation<0, t_i>(i_face);
                    auto face_sign = this->getFiniteElement().template getInnerNeighborSign<0, t_i>(i_face);
                    for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                    {
                        auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                        auto argf = i;
                        auto i_p = this->getFiniteElement().template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                        //
                        auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                        auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                        auto normal_vector = face->getNormalVector(point);
                        auto row_cell_vector = this->getFiniteElement().template getBasisEvaluation<getGradBasis()>(i_p);
                        auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                        auto col_cell_vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(i_p);
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
        
        DenseMatrix<Real, getGradBasisSize<t_element>(), getTensorSize()>
        getSymmetricGradientRhs(
            Integer row,
            Integer col
        )
        const
        requires(t_element.getDim() == t_field.getDimDomain())
        {
            // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
            auto constexpr t_quadrature = GaussQuadrature(getGradientConstructionQuadratureOrder());
            auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
            auto face_offset = t_field_size * getCellBasisSize<t_element>();
            auto rhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getTensorSize()>();
            rhs.setZero();
            for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
            {
                auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<t_quadrature>(i);
                auto point = this->getFiniteElement().template getReferenceQuadraturePoint<t_quadrature>(i);
                auto row_cell_vector = this->getFiniteElement().template getBasisEvaluation<getGradBasis()>(point);
                auto col_cell_vector_i = this->getFiniteElement().template getBasisDerivative<getCellBasis()>(point, col);
                auto col_cell_vector_j = this->getFiniteElement().template getBasisDerivative<getCellBasis()>(point, row);
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
                for (auto const & face : this->getFiniteElement().template getInnerNeighbors<0, t_i>())
                {
                    auto face_orientation = this->getFiniteElement().template getInnerNeighborOrientation<0, t_i>(i_face);
                    auto face_sign = this->getFiniteElement().template getInnerNeighborSign<0, t_i>(i_face);
                    for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                    {
                        auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                        auto argf = i;
                        auto i_p = this->getFiniteElement().template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                        auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                        auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                        auto normal_vector = face->getNormalVector(point);
                        auto row_cell_vector = this->getFiniteElement().template getBasisEvaluation<getGradBasis()>(i_p);
                        auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                        auto col_cell_vector = this->getFiniteElement().template getBasisEvaluation<getCellBasis()>(i_p);
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

} // namespace lolita::core


#endif /* A50FB6B0_0CE6_4A1A_BC52_C8CF6B7F8812 */
