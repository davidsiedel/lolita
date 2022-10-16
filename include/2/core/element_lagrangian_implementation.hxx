#ifndef A747EFCE_DCB2_4622_B340_8E41844643B9
#define A747EFCE_DCB2_4622_B340_8E41844643B9

#include "2/core/_include.hxx"
#include "2/core/element_integration_point.hxx"
#include "2/core/element_potential.hxx"

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
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain, LagrangianConcept auto t_lag>
    struct ElementLagrangian : AbstractElementLagrangian<t_element, t_domain>
    {

    private:

        static constexpr
        Integer
        getSystemSize()
        {
            return LagTraits<t_lag>::template getSize<t_element, t_domain>();
        }

        template<PotentialConcept auto t_potential>
        static constexpr
        Integer
        getPotentialIndex()
        {
            return LagTraits<t_lag>::template getPotentialIndex<t_potential>();
        }

        using Base_ = AbstractElementLagrangian<t_element, t_domain>;

        using FiniteElement_ = FiniteElement<t_element, t_domain>;

        template<PotentialConcept auto t_potential>
        using ElementPotential_ = ElementPotential<t_element, t_domain, t_lag, t_potential>;

        template<PotentialConcept auto t_potential>
        using ElementPotentialPtr_ = std::unique_ptr<ElementPotential_<t_potential>>;

        using ElementPotentials_ = utility::tuple_expansion_t<std::tuple, ElementPotentialPtr_, t_lag.getPotentials()>;

    public:

        explicit
        ElementLagrangian(
            FiniteElement_ const & finite_element
        )
        :
        Base_(t_lag, finite_element)
        {}

        template<PotentialConcept auto t_potential>
        void
        setPotential()
        {
            std::get<getPotentialIndex<t_potential>()>(element_potentials_) = std::make_unique<ElementPotential_<t_potential>>(* this);
        }

        template<PotentialConcept auto t_potential>
        ElementPotential_<t_potential> const &
        getPotential()
        const
        {
            return * std::get<getPotentialIndex<t_potential>()>(element_potentials_);
        }

        template<PotentialConcept auto t_potential>
        ElementPotential_<t_potential> &
        getPotential()
        {
            return * std::get<getPotentialIndex<t_potential>()>(element_potentials_);
        }

        void
        setResidualVector(
            Real const & time
        )
        {
            this->residual_vector_ = this->getElementExternalForces(time) - this->getElementInternalForces();
        }
        
        void
        setResidualVector(
            DenseVectorConcept<Real> auto && vector
        )
        {
            this->residual_vector_ = std::make_unique<DenseVector<Real, getSystemSize()>>(std::forward<decltype(vector)>(vector));
        }

        DenseVector<Real, getSystemSize()> const &
        getResidualVector()
        const
        {
            return * this->residual_vector_;
        }

        void
        setJacobianMatrix(
            DenseMatrixConcept<Real> auto && matrix
        )
        {
            this->jacobian_matrix_ = std::make_unique<DenseMatrix<Real, getSystemSize(), getSystemSize()>>(std::forward<decltype(matrix)>(matrix));
        }

        DenseMatrix<Real, getSystemSize(), getSystemSize()> const &
        getJacobianMatrix()
        const
        {
            return * this->jacobian_matrix_;
        }

        ElementPotentials_ element_potentials_;

        std::unique_ptr<DenseVector<Real, getSystemSize()>> residual_vector_;

        std::unique_ptr<DenseMatrix<Real, getSystemSize(), getSystemSize()>> jacobian_matrix_;

    };

} // namespace lolita::core

#endif /* A747EFCE_DCB2_4622_B340_8E41844643B9 */
