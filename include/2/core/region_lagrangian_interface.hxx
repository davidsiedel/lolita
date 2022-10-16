#ifndef B1F5A3F6_7967_49B1_B46D_832D8837E4A8
#define B1F5A3F6_7967_49B1_B46D_832D8837E4A8

#include "2/core/_include.hxx"
#include "2/core/element_integration_point.hxx"
#include "2/core/element_potential.hxx"
#include "2/core/element_lagrangian_implementation.hxx"
#include "2/core/element_lagrangian_interface.hxx"
#include "2/core/region_potential.hxx"
#include "2/core/region_lagrangian_implementation.hxx"

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

    template<DomainConcept auto t_dim, MeshConcept auto t_domain>
    struct AbstractDomainLagrangian
    {

    private:

        template<LagrangianConcept auto t_lag>
        using DomainLagrangian_ = DomainLagrangian<t_dim, t_domain, t_lag>;

    protected:

        explicit
        AbstractDomainLagrangian(
            LagrangianConcept auto const & lag
        )
        :
        label_(lag.getLabel())
        {}

    public:

        Label const &
        getLabel()
        const
        {
            return label_;
        }

        template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        void
        setPotential(
            auto const &... args
        )
        {
            static_cast<DomainLagrangian_<t_lag> *>(this)->template setPotential<t_potential>(args...);
        }

        template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        auto const &
        getPotential()
        const
        {
            return static_cast<DomainLagrangian_<t_lag> const *>(this)->template getPotential<t_potential>();
        }

        template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        auto &
        getPotential()
        {
            return static_cast<DomainLagrangian_<t_lag> *>(this)->template getPotential<t_potential>();
        }

    protected:

        Label const & label_;

    };

} // namespace lolita::core

#endif /* B1F5A3F6_7967_49B1_B46D_832D8837E4A8 */
