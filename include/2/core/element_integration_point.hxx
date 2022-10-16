#ifndef D233764A_D5AA_4760_ABD8_E3702E8533ED
#define D233764A_D5AA_4760_ABD8_E3702E8533ED

#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/BehaviourData.h>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

#include "2/core/_include.hxx"

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

    /**
     * *********************************************************************************************************************************************************
     */

    template<DomainConcept auto t_dim, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct BehaviourData
    {

        using DomainPotential_ = DomainPotential<t_dim, t_domain, t_lag, t_potential>;
        
        explicit
        BehaviourData(
            DomainPotential_ const & potential
        )
        :
        potential_(potential),
        behavior_data_(mgis::behaviour::BehaviourData(potential.getMgisBehavior()))
        {}

        mgis::behaviour::Behaviour const &
        getMgisBehavior()
        const
        {
            return potential_.getMgisBehavior();
        }

        mgis::behaviour::BehaviourData const &
        getMgisBehaviorData()
        const
        {
            return behavior_data_;
        }

        mgis::behaviour::BehaviourData &
        getMgisBehaviorData()
        {
            return behavior_data_;
        }

        void
        integrate(
            std::atomic<Boolean> & output
        )
        {
            this->behavior_data_.K[0] = 4;
            auto behaviour_data_view = mgis::behaviour::make_view(this->behavior_data_);
            auto res = mgis::behaviour::integrate(behaviour_data_view, this->getMgisBehavior());
            if (res < 1)
            {
                output = false;
            }
        }

    private:

        mgis::behaviour::BehaviourData behavior_data_;

        DomainPotential_ const & potential_;
        
    };

    template<ShapeConcept auto t_element, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct IntegrationPoint2
    {

    private:

        template<LinearOperatorConcept auto t_strain>
        static constexpr
        Integer
        getStrainOperatorNumRows()
        {
            return LinearOperatorTraits<t_strain>::template getSize<t_domain>();
        }

        template<LinearOperatorConcept auto t_strain>
        static constexpr
        Integer
        getStrainOperatorNumCols()
        {
            return DiscretizationTraits<t_strain.getField()>::template getTensorSpaceSize<t_element, t_domain>();
        }

        template<LinearOperatorConcept auto t_strain>
        static constexpr
        Integer
        getStrainIndex()
        {
            return PotentialTraits2<t_potential>::template getLinearOperatorIndex<t_strain>();
        }

        using FiniteElement_ = FiniteElement<t_element, t_domain>;

        using ElementPotential_ = ElementPotential<t_element, t_domain, t_lag, t_potential>;

        template<DomainConcept auto t_dim, MeshConcept auto t__domain>
        using BehaviorData_1 = BehaviourData<t_dim, t__domain, t_lag, t_potential>;

        using BehaviorData_ = ElementTraits<t_element>::template DomainConnectivity<BehaviorData_1, t_domain>;

        template<LinearOperatorConcept auto t_strain>
        using StrainOperator_ = DenseMatrix<Real, getStrainOperatorNumRows<t_strain>(), getStrainOperatorNumCols<t_strain>()>;

        template<LinearOperatorConcept auto t_strain>
        using StrainOperatorPtr_ = std::unique_ptr<StrainOperator_<t_strain>>;

        using StrainOperators_ = utility::tuple_expansion_t<std::tuple, StrainOperatorPtr_, t_potential.getLinearOperators()>;

    public:
        
        IntegrationPoint2(
            ElementPotential_ const & potential,
            Integer index
        )
        :
        potential_(potential),
        index_(index),
        behavior_data_(potential.getFiniteElement().getDomain()->template getLagrangian<t_lag>().template getPotential<t_lag, t_potential>())
        {}

        FiniteElement_ const &
        getFiniteElement()
        const
        {
            return potential_.getFiniteElement();
        }

        Integer
        getIndex()
        const
        {
            return index_;
        }

        auto const &
        getDomainPotential()
        {
            return this->getFiniteElement().getDomain()->template getLagrangian<t_lag>().template getPotential<t_lag, t_potential>();
        }

        BehaviorData_ const &
        getBhv()
        const
        {
            return behavior_data_;
        }

        BehaviorData_ &
        getBhv()
        {
            return behavior_data_;
        }

        auto const &
        getMgisBhv()
        const
        {
            return this->getBhv().getMgisBehaviorData();
        }

        auto &
        getMgisBhv()
        {
            return this->getBhv().getMgisBehaviorData();
        }

        void
        integrate(
            std::atomic<Boolean> & output
        )
        {
            this->getBhv().integrate(output);
        }

        void
        reserve()
        {
            mgis::behaviour::update(this->getMgisBhv());
        }

        void
        recover()
        {
            mgis::behaviour::revert(this->getMgisBhv());
        }

        void
        setMaterialProperty(
            std::basic_string<Character> && material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto value = std::forward<std::function<Real(Point const &)>>(function)(this->getCurrentCoordinates());
            mgis::behaviour::setMaterialProperty(this->getMgisBhv().s0, std::forward<std::basic_string<Character>>(material_property_label), value);
            mgis::behaviour::setMaterialProperty(this->getMgisBhv().s1, std::forward<std::basic_string<Character>>(material_property_label), value);
        }

        void
        setExternalVariable(
            std::basic_string<Character> && material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto value = std::forward<std::function<Real(Point const &)>>(function)(this->getCurrentCoordinates());
            mgis::behaviour::setExternalStateVariable(this->getMgisBhv().s0, std::forward<std::basic_string<Character>>(material_property_label), value);
            mgis::behaviour::setExternalStateVariable(this->getMgisBhv().s1, std::forward<std::basic_string<Character>>(material_property_label), value);
        }

        template<LinearOperatorConcept auto t_strain1, LinearOperatorConcept auto t_strain2>
        algebra::View<DenseMatrix<Real, LinearOperatorTraits<t_strain1>::template getSize<t_domain>(), LinearOperatorTraits<t_strain2>::template getSize<t_domain>()> const>
        getJacobianMatrixBlock()
        const
        {
            auto constexpr jacobian_size = PotentialTraits2<t_potential>::template getSize<t_domain>();
            auto constexpr jacobian_block_num_rows = LinearOperatorTraits<t_strain1>::template getSize<t_domain>();
            auto constexpr jacobian_block_num_cols = LinearOperatorTraits<t_strain2>::template getSize<t_domain>();
            auto jacobian_block_row_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain1>();
            auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain2>();
            auto block_offset = jacobian_size * jacobian_block_row_offset + jacobian_block_col_offset;
            auto const * mgis_data = this->getMgisBhv().K.data();
            return algebra::View<DenseMatrix<Real, jacobian_block_num_rows, jacobian_block_num_cols> const>(mgis_data + block_offset);
        }

        template<LinearOperatorConcept auto t_strain>
        algebra::View<DenseVector<Real, LinearOperatorTraits<t_strain>::template getSize<t_domain>()> const>
        getStressVectorBlock()
        const
        {
            auto constexpr jacobian_block_num_cols = LinearOperatorTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain>();
            auto const * mgis_data = this->getMgisBhv().s1.thermodynamic_forces.data();
            return algebra::View<DenseVector<Real, jacobian_block_num_cols> const>(mgis_data + jacobian_block_col_offset);
        }

        template<LinearOperatorConcept auto t_strain>
        algebra::View<DenseVector<Real, LinearOperatorTraits<t_strain>::template getSize<t_domain>()> const>
        getStrainVectorBlock()
        const
        {
            auto constexpr jacobian_block_num_cols = LinearOperatorTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain>();
            auto const * mgis_data = this->getMgisBhv().s1.gradients.data();
            return algebra::View<DenseVector<Real, jacobian_block_num_cols> const>(mgis_data + jacobian_block_col_offset);
        }

        template<LinearOperatorConcept auto t_strain>
        void
        setStrainVectorBlock()
        {
            auto constexpr jacobian_block_num_cols = LinearOperatorTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain>();
            auto * mgis_data = this->getMgisBhv().s1.gradients.data();
            auto strain_block = algebra::View<DenseVector<Real, jacobian_block_num_cols>>(mgis_data + jacobian_block_col_offset);
            auto ukns = this->getFiniteElement().template getDiscreteField<t_strain.getField()>().template getUnknownCoefficients<t_strain.getField()>();
            if (this->template hasStrainOperator<t_strain>())
            {
                strain_block = this->template getStrainOperator<t_strain>() * ukns;
            }
            else
            {
                strain_block = this->template letStrainOperator<t_strain>() * ukns;
            }
        }
        
        Point
        getCurrentCoordinates()
        const
        {
            return getFiniteElement().template getCurrentQuadraturePoint<t_potential.getQuadrature()>(this->getIndex());
        }
        
        algebra::View<Point const>
        getReferenceCoordinates()
        const
        {
            return getFiniteElement().template getReferenceQuadraturePoint<t_potential.getQuadrature()>(this->getIndex());
        }
        
        Real
        getCurrentWeight()
        const
        {
            return getFiniteElement().template getCurrentQuadratureWeight<t_potential.getQuadrature()>(this->getIndex());
        }
        
        Real
        getReferenceWeight()
        const
        {
            return getFiniteElement().template getReferenceQuadratureWeight<t_potential.getQuadrature()>(this->getIndex());
        }

        template<LinearOperatorConcept auto t_strain>
        Boolean
        hasStrainOperator()
        const
        {
            return std::get<getStrainIndex<t_strain>()>(strain_operators_) != nullptr;
        }
        
        template<LinearOperatorConcept auto t_strain>
        StrainOperator_<t_strain> const &
        getStrainOperator()
        const
        {
            return * std::get<getStrainIndex<t_strain>()>(strain_operators_);
        }

        template<LinearOperatorConcept auto t_strain>
        StrainOperator_<t_strain>
        letStrainOperator()
        const
        {
            return getFiniteElement().template getDiscreteField<t_strain.getField()>().template letLinearOperator<t_strain>(this->getReferenceCoordinates());
        }

        template<LinearOperatorConcept auto t_strain>
        void
        setStrainOperator()
        {
            std::get<getStrainIndex<t_strain>()>(strain_operators_) = std::make_unique<StrainOperator_<t_strain>>(this->template letStrainOperator<t_strain>());
        }

        Integer index_;

        StrainOperators_ strain_operators_;

        BehaviorData_ behavior_data_;

        ElementPotential_ const & potential_;
        
    };

} // namespace lolita::core

#endif /* D233764A_D5AA_4760_ABD8_E3702E8533ED */
