#ifndef B04B8EFF_C85E_4FB3_BF25_2A2C6AE61045
#define B04B8EFF_C85E_4FB3_BF25_2A2C6AE61045

#include "2/core/_include.hxx"
#include "2/core/element_integration_point.hxx"
#include "2/core/element_potential.hxx"
#include "2/core/element_lagrangian_implementation.hxx"
#include "2/core/element_lagrangian_interface.hxx"

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

    template<DomainConcept auto t_dim, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct DomainPotential
    {

    private:

        using DomainLagrangian_ = DomainLagrangian<t_dim, t_domain, t_lag>;

    public:

        DomainPotential(
            DomainLagrangian_ const & lag,
            auto const &... args
        )
        :
        lag_(lag),
        behavior_(mgis::behaviour::load(args...))
        {}

        DomainLagrangian_ const &
        getLag()
        const
        {
            return lag_;
        }

        mgis::behaviour::Behaviour const &
        getMgisBehavior()
        const
        {
            return behavior_;
        }

    private:

        DomainLagrangian_ const & lag_;
        
        mgis::behaviour::Behaviour behavior_;

    };

} // namespace lolita::core


#endif /* B04B8EFF_C85E_4FB3_BF25_2A2C6AE61045 */
