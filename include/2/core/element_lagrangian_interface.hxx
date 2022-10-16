#ifndef C58FD1C0_9C9C_4A68_A949_6F853DFBC423
#define C58FD1C0_9C9C_4A68_A949_6F853DFBC423

#include "2/core/_include.hxx"
#include "2/core/element_integration_point.hxx"
#include "2/core/element_potential.hxx"
#include "2/core/element_lagrangian_implementation.hxx"

namespace lolita::core
{
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    struct AbstractElementLagrangian
    {

    private:

        template<LagrangianConcept auto t_lag>
        static constexpr
        Integer
        getSystemSize()
        {
            return LagTraits<t_lag>::template getSize<t_element, t_domain>();
        }

        using FiniteElement_ = FiniteElement<t_element, t_domain>;

        template<LagrangianConcept auto t_lag>
        using ElementLagrangian_ = ElementLagrangian<t_element, t_domain, t_lag>;

    protected:

        AbstractElementLagrangian(
            LagrangianConcept auto const & lag,
            FiniteElement_ const & finite_element
        )
        :
        label_(lag.getLabel()),
        finite_element_(finite_element)
        {}

    public:

        Label const &
        getLabel()
        const
        {
            return label_;
        }

        FiniteElement_ const &
        getFiniteElement()
        const
        {
            return finite_element_;
        }

        template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        void
        setPotential()
        {
            static_cast<ElementLagrangian_<t_lag> *>(this)->template setPotential<t_potential>();
        }

        template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        auto const &
        getPotential()
        const
        {
            return static_cast<ElementLagrangian_<t_lag> const *>(this)->template getPotential<t_potential>();
        }

        template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        auto &
        getPotential()
        {
            return static_cast<ElementLagrangian_<t_lag> *>(this)->template getPotential<t_potential>();
        }

        template<LagrangianConcept auto t_lag>
        void
        setJacobianMatrix(
            DenseMatrixConcept<Real> auto && matrix
        )
        {
            static_cast<ElementLagrangian_<t_lag> *>(this)->setJacobianMatrix(std::forward<decltype(matrix)>(matrix));
        }

        template<LagrangianConcept auto t_lag>
        void
        setJacobianMatrix()
        {
            static_cast<ElementLagrangian_<t_lag> *>(this)->setJacobianMatrix(getElementJacobianMatrix<t_lag>());
        }

        template<LagrangianConcept auto t_lag>
        auto const &
        getJacobianMatrix()
        const
        {
            return static_cast<ElementLagrangian_<t_lag> const *>(this)->getJacobianMatrix();
        }

        template<LagrangianConcept auto t_lag>
        void
        setResidualVector(
            DenseVectorConcept<Real> auto && vector
        )
        {
            static_cast<ElementLagrangian_<t_lag> *>(this)->setResidualVector(std::forward<decltype(vector)>(vector));
        }

        template<LagrangianConcept auto t_lag>
        void
        setResidualVector()
        {
            static_cast<ElementLagrangian_<t_lag> *>(this)->setResidualVector(getElementInternalForces<t_lag>());
        }

        template<LagrangianConcept auto t_lag>
        auto const &
        getResidualVector()
        const
        {
            return static_cast<ElementLagrangian_<t_lag> const *>(this)->getResidualVector();
        }

        template<LagrangianConcept auto t_lag>
        DenseVector<Real, getSystemSize<t_lag>()>
        getElementInternalForces()
        const
        {
            auto internal_forces = DenseVector<Real, getSystemSize<t_lag>()>();
            internal_forces.setZero();
            auto set_i = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = LagTraits<t_lag>::template getPotential<t_i>();
                auto constexpr j_mapping = PotentialTraits2<potential>::template getLinearOperator<t_j>();
                auto constexpr j_field = j_mapping.getField();
                auto constexpr size_j = DiscretizationTraits<j_field>::template getTensorSpaceSize<t_element, t_domain>();
                auto constexpr offset_j = LagTraits<t_lag>::template getFieldOffset<t_element, t_domain, j_field>();
                for (auto const & ip : this->getPotential<t_lag, potential>().getIntegrationPoints())
                {
                    auto const & mat0 = ip->template getStrainOperator<j_mapping>();
                    auto const mat = ip->template getStressVectorBlock<j_mapping>();
                    internal_forces.template segment<size_j>(offset_j) += ip->getCurrentWeight() * mat0.transpose() * mat;
                }
                if constexpr (t_j < PotentialTraits2<potential>::getNumLinearOperators() - 1)
                {
                    t_set_i.template operator()<t_i, t_j + 1>(t_set_i);
                }
                else if constexpr (t_i < LagTraits<t_lag>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1, 0>(t_set_i);
                }
            };
            set_i(set_i);
            return internal_forces;
        }

        template<LagrangianConcept auto t_lag>
        DenseVector<Real, getSystemSize<t_lag>()>
        getElementExternalForces(
            Real const & time
        )
        const
        {
            auto constexpr quadrature = GaussQuadrature(4);
            auto external_forces = DenseVector<Real, getSystemSize<t_lag>()>();
            external_forces.setZero();
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr field = LagTraits<t_lag>::template getField<t_i>();
                auto constexpr size_j = DiscretizationTraits<field>::template getTensorSpaceSize<t_element, t_domain>();
                auto constexpr offset_j = LagTraits<t_lag>::template getFieldOffset<t_element, t_domain, field>();
                for (auto const & domain : this->getFiniteElement().getDomains())
                {
                    if (domain->template hasDiscreteField<field>())
                    {
                        if (domain->template getDiscreteField<field>().hasLoads())
                        {
                            for (auto const & load : domain->template getDiscreteField<field>().getLoads())
                            {
                                for (auto i = 0; i < QuadratureTraits<quadrature>::template Rule<t_element>::getSize(); i++)
                                {
                                    auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<quadrature>(i);
                                    auto reference_point = this->getFiniteElement().template getReferenceQuadraturePoint<quadrature>(i);
                                    auto current_point = this->getFiniteElement().template getCurrentQuadraturePoint<quadrature>(i);
                                    auto vector = this->getFiniteElement().template getFieldDualVector<field>(reference_point, load.getRow(), load.getCol());
                                    external_forces.template segment<size_j>(offset_j) += weight * load.getValue(current_point, time) * vector;
                                }
                            }
                        }
                    }
                }
                if constexpr (t_i < LagTraits<t_lag>::getNumFields() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
            return external_forces;
        }
        
        template<LagrangianConcept auto t_lag>
        DenseMatrix<Real, getSystemSize<t_lag>(), getSystemSize<t_lag>()>
        getElementJacobianMatrix()
        const
        {
            auto jacobian_matrix = DenseMatrix<Real, getSystemSize<t_lag>(), getSystemSize<t_lag>()>();
            jacobian_matrix.setZero();
            auto set_i = [&] <Integer t_i = 0, Integer t_j = 0, Integer t_k = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = LagTraits<t_lag>::template getPotential<t_i>();
                auto constexpr j_mapping = PotentialTraits2<potential>::template getLinearOperator<t_j>();
                auto constexpr j_field = j_mapping.getField();
                auto constexpr k_mapping = PotentialTraits2<potential>::template getLinearOperator<t_k>();
                auto constexpr k_field = k_mapping.getField();
                auto constexpr size_j = DiscretizationTraits<j_field>::template getTensorSpaceSize<t_element, t_domain>();
                auto constexpr size_k = DiscretizationTraits<k_field>::template getTensorSpaceSize<t_element, t_domain>();
                auto constexpr offset_j = LagTraits<t_lag>::template getFieldOffset<t_element, t_domain, j_field>();
                auto constexpr offset_k = LagTraits<t_lag>::template getFieldOffset<t_element, t_domain, k_field>();
                auto const & frm = this->template getPotential<t_lag, potential>();
                for (auto const & ip : frm.getIntegrationPoints())
                {
                    auto const & mat0 = ip->template getStrainOperator<j_mapping>();
                    auto const & mat1 = ip->template getStrainOperator<k_mapping>();
                    auto const mat = ip->template getJacobianMatrixBlock<j_mapping, k_mapping>();
                    jacobian_matrix.template block<size_j, size_k>(offset_j, offset_k) += ip->getCurrentWeight() * mat0.transpose() * mat * mat1;
                }
                if constexpr (t_k < PotentialTraits2<potential>::getNumLinearOperators() - 1)
                {
                    t_set_i.template operator()<t_i, t_j, t_k + 1>(t_set_i);
                }
                if constexpr (t_j < PotentialTraits2<potential>::getNumLinearOperators() - 1)
                {
                    t_set_i.template operator()<t_i, t_j + 1, 0>(t_set_i);
                }
                else if constexpr (t_i < LagTraits<t_lag>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1, 0, 0>(t_set_i);
                }
            };
            set_i(set_i);
            // std::cout << "jacobian_matrix : " << std::endl;
            // std::cout << jacobian_matrix << std::endl;
            return jacobian_matrix;
        }

    protected:

        Label const & label_;

        FiniteElement_ const & finite_element_;

    };
    
} // namespace lolita::core

#endif /* C58FD1C0_9C9C_4A68_A949_6F853DFBC423 */
