#ifndef A07BE16D_BAF0_4082_9FFC_4DA19EA578F3
#define A07BE16D_BAF0_4082_9FFC_4DA19EA578F3

#include "2/core/_include.hxx"
#include "2/core/element_integration_point.hxx"

namespace lolita::core
{
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    struct FiniteElement;

    template<DomainConcept auto t_dim, MeshConcept auto t_domain>
    struct FiniteDomain;

    template<DomainConcept auto t_dim, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct DomainPotential;
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct ElementPotential;
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain, LagrangianConcept auto t_lag>
    struct ElementLagrangian;
    
    template<DomainConcept auto t_dim, MeshConcept auto t_domain, LagrangianConcept auto t_lag>
    struct DomainLagrangian;
    
    template<DomainConcept auto t_dim, MeshConcept auto t_domain>
    struct AbstractDomainLagrangian;

    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    struct AbstractElementLagrangian;
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct ElementPotential
    {

        static constexpr
        Integer
        getNumQuadraturePoints()
        {
            return QuadratureTraits<t_potential.getQuadrature()>::template Rule<t_element>::getSize();
        }

        using ElementLagrangian_ = ElementLagrangian<t_element, t_domain, t_lag>;

        using FiniteElement_ = FiniteElement<t_element, t_domain>;

        using IntegrationPoint_ = IntegrationPoint2<t_element, t_domain, t_lag, t_potential>;

    private:

        void
        setIntegrationPoints()
        {
            auto quadrature_index = Integer(0);
            for (auto & ip : integration_points_)
            {
                ip = std::make_unique<IntegrationPoint_>(* this, quadrature_index);
                quadrature_index ++;
            }
        }

    public:
        
        explicit
        ElementPotential(
            ElementLagrangian_ const & lag
        )
        :
        lag_(lag)
        {
            setIntegrationPoints();
        }

        FiniteElement_ const &
        getFiniteElement()
        const
        {
            return lag_.getFiniteElement();
        }

        ElementLagrangian_ const &
        getLag()
        const
        {
            return lag_;
        }

        std::array<std::unique_ptr<IntegrationPoint_>, getNumQuadraturePoints()> const &
        getIntegrationPoints()
        const
        {
            return integration_points_;
        }

        std::array<std::unique_ptr<IntegrationPoint_>, getNumQuadraturePoints()> &
        getIntegrationPoints()
        {
            return integration_points_;
        }

        void
        setStrainOperators()
        {
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr mapping = PotentialTraits2<t_potential>::template getLinearOperator<t_i>();
                for (auto & ip : integration_points_)
                {
                    ip->template setStrainOperator<mapping>();
                }
                if constexpr (t_i < PotentialTraits2<t_potential>::getNumLinearOperators() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
        }

        void
        setStrains()
        {
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr mapping = PotentialTraits2<t_potential>::template getLinearOperator<t_i>();
                for (auto & ip : integration_points_)
                {
                    ip->template setStrainVectorBlock<mapping>();
                }
                if constexpr (t_i < PotentialTraits2<t_potential>::getNumLinearOperators() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
        }

        void
        setMaterialProperty(
            std::basic_string<Character> && label,
            std::function<Real(Point const &)> && function
        )
        {
            for (auto & ip : integration_points_)
            {
                ip->setMaterialProperty(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
            }
        }

        void
        setExternalVariable(
            std::basic_string<Character> && label,
            std::function<Real(Point const &)> && function
        )
        {
            for (auto & ip : integration_points_)
            {
                ip->setExternalVariable(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
            }
        }

        void
        integrateConstitutiveEquation(
            std::atomic<Boolean> & output_handler
        )
        {
            for (auto & ip : integration_points_)
            {
                ip->integrate(output_handler);
            }
        }

        void
        reserve()
        {
            for (auto & ip : integration_points_)
            {
                ip->reserve();
            }            
        }

        void
        recover()
        {
            for (auto & ip : integration_points_)
            {
                ip->recover();
            }            
        }

        ElementLagrangian_ const & lag_;

        std::array<std::unique_ptr<IntegrationPoint_>, getNumQuadraturePoints()> integration_points_;

    };

} // namespace lolita::core

#endif /* A07BE16D_BAF0_4082_9FFC_4DA19EA578F3 */
