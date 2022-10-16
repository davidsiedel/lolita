#ifndef D7E6A22A_F9B0_41A7_8145_DAE411A2152F
#define D7E6A22A_F9B0_41A7_8145_DAE411A2152F

#include "2/core/_include.hxx"
#include "2/core/element_integration_point.hxx"
#include "2/core/element_potential.hxx"
#include "2/core/element_lagrangian_implementation.hxx"
#include "2/core/element_lagrangian_interface.hxx"
#include "2/core/region_potential.hxx"

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
    
    template<DomainConcept auto t_dim, MeshConcept auto t_domain, LagrangianConcept auto t_lag>
    struct DomainLagrangian : AbstractDomainLagrangian<t_dim, t_domain>
    {

    private:

        template<PotentialConcept auto t_potential>
        static constexpr
        Integer
        getPotentialIndex()
        {
            return LagTraits<t_lag>::template getPotentialIndex<t_potential>();
        }

        using Base_ = AbstractDomainLagrangian<t_dim, t_domain>;

        using Domain_ = FiniteDomain<t_dim, t_domain>;

        template<PotentialConcept auto t_potential>
        using DomainPotential_ = DomainPotential<t_dim, t_domain, t_lag, t_potential>;

        template<PotentialConcept auto t_potential>
        using DomainPotentialPtr_ = std::unique_ptr<DomainPotential_<t_potential>>;

        using DomainPotentials_ = utility::tuple_expansion_t<std::tuple, DomainPotentialPtr_, t_lag.getPotentials()>;

    public:

        explicit
        DomainLagrangian(
            Domain_ const & domain
        )
        :
        Base_(t_lag),
        domain_(domain)
        {}

        Domain_ const &
        getFiniteDomain()
        const
        {
            return domain_;
        }

        template<PotentialConcept auto t_potential>
        void
        setPotential(
            auto const &... args
        )
        {
            std::get<getPotentialIndex<t_potential>()>(domain_potentials_) = std::make_unique<DomainPotential_<t_potential>>(* this, args...);
        }

        template<PotentialConcept auto t_potential>
        DomainPotential_<t_potential> const &
        getPotential()
        const
        {
            return * std::get<getPotentialIndex<t_potential>()>(domain_potentials_);
        }

        template<PotentialConcept auto t_potential>
        DomainPotential_<t_potential> &
        getPotential()
        {
            return * std::get<getPotentialIndex<t_potential>()>(domain_potentials_);
        }

    private:

        Domain_ const & domain_;

        DomainPotentials_ domain_potentials_;

    };

} // namespace lolita::core


#endif /* D7E6A22A_F9B0_41A7_8145_DAE411A2152F */
