#ifndef A0EFE384_4F56_4FAA_84B9_65C5B2053C0D
#define A0EFE384_4F56_4FAA_84B9_65C5B2053C0D

#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/BehaviourData.h>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

#include "2/core/_include.hxx"
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

    // template<DomainConcept auto t_dim, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    // struct BehaviourData
    // {

    //     using DomainPotential_ = DomainPotential<t_dim, t_domain, t_lag, t_potential>;
        
    //     explicit
    //     BehaviourData(
    //         DomainPotential_ const & potential
    //     )
    //     :
    //     potential_(potential),
    //     behavior_data_(mgis::behaviour::BehaviourData(potential.getMgisBehavior()))
    //     {}

    //     mgis::behaviour::Behaviour const &
    //     getMgisBehavior()
    //     const
    //     {
    //         return potential_.getMgisBehavior();
    //     }

    //     mgis::behaviour::BehaviourData const &
    //     getMgisBehaviorData()
    //     const
    //     {
    //         return behavior_data_;
    //     }

    //     mgis::behaviour::BehaviourData &
    //     getMgisBehaviorData()
    //     {
    //         return behavior_data_;
    //     }

    //     void
    //     integrate(
    //         std::atomic<Boolean> & output
    //     )
    //     {
    //         this->behavior_data_.K[0] = 4;
    //         auto behaviour_data_view = mgis::behaviour::make_view(this->behavior_data_);
    //         auto res = mgis::behaviour::integrate(behaviour_data_view, this->getMgisBehavior());
    //         if (res < 1)
    //         {
    //             output = false;
    //         }
    //     }

    // private:

    //     mgis::behaviour::BehaviourData behavior_data_;

    //     DomainPotential_ const & potential_;
        
    // };

    // template<ShapeConcept auto t_element, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    // struct IntegrationPoint2
    // {

    // private:

    //     template<LinearOperatorConcept auto t_strain>
    //     static constexpr
    //     Integer
    //     getStrainOperatorNumRows()
    //     {
    //         return LinearOperatorTraits<t_strain>::template getSize<t_domain>();
    //     }

    //     template<LinearOperatorConcept auto t_strain>
    //     static constexpr
    //     Integer
    //     getStrainOperatorNumCols()
    //     {
    //         return DiscretizationTraits<t_strain.getField()>::template getTensorSpaceSize<t_element, t_domain>();
    //     }

    //     template<LinearOperatorConcept auto t_strain>
    //     static constexpr
    //     Integer
    //     getStrainIndex()
    //     {
    //         return PotentialTraits2<t_potential>::template getLinearOperatorIndex<t_strain>();
    //     }

    //     using FiniteElement_ = FiniteElement<t_element, t_domain>;

    //     using ElementPotential_ = ElementPotential<t_element, t_domain, t_lag, t_potential>;

    //     template<DomainConcept auto t_dim, MeshConcept auto t__domain>
    //     using BehaviorData_1 = BehaviourData<t_dim, t__domain, t_lag, t_potential>;

    //     using BehaviorData_ = ElementTraits<t_element>::template DomainConnectivity<BehaviorData_1, t_domain>;

    //     template<LinearOperatorConcept auto t_strain>
    //     using StrainOperator_ = DenseMatrix<Real, getStrainOperatorNumRows<t_strain>(), getStrainOperatorNumCols<t_strain>()>;

    //     template<LinearOperatorConcept auto t_strain>
    //     using StrainOperatorPtr_ = std::unique_ptr<StrainOperator_<t_strain>>;

    //     using StrainOperators_ = utility::tuple_expansion_t<std::tuple, StrainOperatorPtr_, t_potential.getLinearOperators()>;

    // public:
        
    //     IntegrationPoint2(
    //         ElementPotential_ const & potential,
    //         Integer index
    //     )
    //     :
    //     potential_(potential),
    //     index_(index),
    //     behavior_data_(potential.getFiniteElement().getDomain()->template getLagrangian<t_lag>().template getPotential<t_lag, t_potential>())
    //     {}

    //     FiniteElement_ const &
    //     getFiniteElement()
    //     const
    //     {
    //         return potential_.getFiniteElement();
    //     }

    //     Integer
    //     getIndex()
    //     const
    //     {
    //         return index_;
    //     }

    //     auto const &
    //     getDomainPotential()
    //     {
    //         return this->getFiniteElement().getDomain()->template getLagrangian<t_lag>().template getPotential<t_lag, t_potential>();
    //     }

    //     BehaviorData_ const &
    //     getBhv()
    //     const
    //     {
    //         return behavior_data_;
    //     }

    //     BehaviorData_ &
    //     getBhv()
    //     {
    //         return behavior_data_;
    //     }

    //     auto const &
    //     getMgisBhv()
    //     const
    //     {
    //         return this->getBhv().getMgisBehaviorData();
    //     }

    //     auto &
    //     getMgisBhv()
    //     {
    //         return this->getBhv().getMgisBehaviorData();
    //     }

    //     void
    //     integrate(
    //         std::atomic<Boolean> & output
    //     )
    //     {
    //         this->getBhv().integrate(output);
    //     }

    //     void
    //     reserve()
    //     {
    //         mgis::behaviour::update(this->getMgisBhv());
    //     }

    //     void
    //     recover()
    //     {
    //         mgis::behaviour::revert(this->getMgisBhv());
    //     }

    //     void
    //     setMaterialProperty(
    //         std::basic_string<Character> && material_property_label,
    //         std::function<Real(Point const &)> && function
    //     )
    //     {
    //         auto value = std::forward<std::function<Real(Point const &)>>(function)(this->getCurrentCoordinates());
    //         std::cout << material_property_label << " : " << value << std::endl;
    //         mgis::behaviour::setMaterialProperty(this->getMgisBhv().s0, std::forward<std::basic_string<Character>>(material_property_label), value);
    //         mgis::behaviour::setMaterialProperty(this->getMgisBhv().s1, std::forward<std::basic_string<Character>>(material_property_label), value);
    //     }

    //     void
    //     setExternalVariable(
    //         std::basic_string<Character> && material_property_label,
    //         std::function<Real(Point const &)> && function
    //     )
    //     {
    //         auto value = std::forward<std::function<Real(Point const &)>>(function)(this->getCurrentCoordinates());
    //         std::cout << material_property_label << " : " << value << std::endl;
    //         mgis::behaviour::setExternalStateVariable(this->getMgisBhv().s0, std::forward<std::basic_string<Character>>(material_property_label), value);
    //         mgis::behaviour::setExternalStateVariable(this->getMgisBhv().s1, std::forward<std::basic_string<Character>>(material_property_label), value);
    //     }

    //     template<LinearOperatorConcept auto t_strain1, LinearOperatorConcept auto t_strain2>
    //     algebra::View<DenseMatrix<Real, LinearOperatorTraits<t_strain1>::template getSize<t_domain>(), LinearOperatorTraits<t_strain2>::template getSize<t_domain>()> const>
    //     getJacobianMatrixBlock()
    //     const
    //     {
    //         auto constexpr jacobian_size = PotentialTraits2<t_potential>::template getSize<t_domain>();
    //         auto constexpr jacobian_block_num_rows = LinearOperatorTraits<t_strain1>::template getSize<t_domain>();
    //         auto constexpr jacobian_block_num_cols = LinearOperatorTraits<t_strain2>::template getSize<t_domain>();
    //         auto jacobian_block_row_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain1>();
    //         auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain2>();
    //         auto block_offset = jacobian_size * jacobian_block_row_offset + jacobian_block_col_offset;
    //         auto const * mgis_data = this->getMgisBhv().K.data();
    //         return algebra::View<DenseMatrix<Real, jacobian_block_num_rows, jacobian_block_num_cols> const>(mgis_data + block_offset);
    //     }

    //     template<LinearOperatorConcept auto t_strain>
    //     algebra::View<DenseVector<Real, LinearOperatorTraits<t_strain>::template getSize<t_domain>()> const>
    //     getStressVectorBlock()
    //     const
    //     {
    //         auto constexpr jacobian_block_num_cols = LinearOperatorTraits<t_strain>::template getSize<t_domain>();
    //         auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain>();
    //         auto const * mgis_data = this->getMgisBhv().s1.thermodynamic_forces.data();
    //         return algebra::View<DenseVector<Real, jacobian_block_num_cols> const>(mgis_data + jacobian_block_col_offset);
    //     }

    //     template<LinearOperatorConcept auto t_strain>
    //     algebra::View<DenseVector<Real, LinearOperatorTraits<t_strain>::template getSize<t_domain>()> const>
    //     getStrainVectorBlock()
    //     const
    //     {
    //         auto constexpr jacobian_block_num_cols = LinearOperatorTraits<t_strain>::template getSize<t_domain>();
    //         auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain>();
    //         auto const * mgis_data = this->getMgisBhv().s1.gradients.data();
    //         return algebra::View<DenseVector<Real, jacobian_block_num_cols> const>(mgis_data + jacobian_block_col_offset);
    //     }

    //     template<LinearOperatorConcept auto t_strain>
    //     void
    //     setStrainVectorBlock()
    //     {
    //         auto constexpr jacobian_block_num_cols = LinearOperatorTraits<t_strain>::template getSize<t_domain>();
    //         auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getLinearOperatorOffset<t_domain, t_strain>();
    //         auto * mgis_data = this->getMgisBhv().s1.gradients.data();
    //         auto strain_block = algebra::View<DenseVector<Real, jacobian_block_num_cols>>(mgis_data + jacobian_block_col_offset);
    //         auto ukns = this->getFiniteElement().template getDiscreteField<t_strain.getField()>().template getUnknownCoefficients<t_strain.getField()>();
    //         if (this->template hasStrainOperator<t_strain>())
    //         {
    //             strain_block = this->template getStrainOperator<t_strain>() * ukns;
    //         }
    //         else
    //         {
    //             strain_block = this->template letStrainOperator<t_strain>() * ukns;
    //         }
    //     }
        
    //     Point
    //     getCurrentCoordinates()
    //     const
    //     {
    //         return getFiniteElement().template getCurrentQuadraturePoint<t_potential.getQuadrature()>(this->getIndex());
    //     }
        
    //     algebra::View<Point const>
    //     getReferenceCoordinates()
    //     const
    //     {
    //         return getFiniteElement().template getReferenceQuadraturePoint<t_potential.getQuadrature()>(this->getIndex());
    //     }
        
    //     Real
    //     getCurrentWeight()
    //     const
    //     {
    //         return getFiniteElement().template getCurrentQuadratureWeight<t_potential.getQuadrature()>(this->getIndex());
    //     }
        
    //     Real
    //     getReferenceWeight()
    //     const
    //     {
    //         return getFiniteElement().template getReferenceQuadratureWeight<t_potential.getQuadrature()>(this->getIndex());
    //     }

    //     template<LinearOperatorConcept auto t_strain>
    //     Boolean
    //     hasStrainOperator()
    //     const
    //     {
    //         return std::get<getStrainIndex<t_strain>()>(strain_operators_) != nullptr;
    //     }
        
    //     template<LinearOperatorConcept auto t_strain>
    //     StrainOperator_<t_strain> const &
    //     getStrainOperator()
    //     const
    //     {
    //         return * std::get<getStrainIndex<t_strain>()>(strain_operators_);
    //     }

    //     template<LinearOperatorConcept auto t_strain>
    //     StrainOperator_<t_strain>
    //     letStrainOperator()
    //     const
    //     {
    //         return getFiniteElement().template getDiscreteField<t_strain.getField()>().template letLinearOperator<t_strain>(this->getReferenceCoordinates());
    //     }

    //     template<LinearOperatorConcept auto t_strain>
    //     void
    //     setStrainOperator()
    //     {
    //         std::get<getStrainIndex<t_strain>()>(strain_operators_) = std::make_unique<StrainOperator_<t_strain>>(this->template letStrainOperator<t_strain>());
    //     }

    //     Integer index_;

    //     StrainOperators_ strain_operators_;

    //     BehaviorData_ behavior_data_;

    //     ElementPotential_ const & potential_;
        
    // };
    
    // template<ShapeConcept auto t_element, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    // struct ElementPotential
    // {

    //     static constexpr
    //     Integer
    //     getNumQuadraturePoints()
    //     {
    //         return QuadratureTraits<t_potential.getQuadrature()>::template Rule<t_element>::getSize();
    //     }

    //     using ElementLagrangian_ = ElementLagrangian<t_element, t_domain, t_lag>;

    //     using FiniteElement_ = FiniteElement<t_element, t_domain>;

    //     using IntegrationPoint_ = IntegrationPoint2<t_element, t_domain, t_lag, t_potential>;

    // private:

    //     void
    //     setIntegrationPoints()
    //     {
    //         auto quadrature_index = Integer(0);
    //         for (auto & ip : integration_points_)
    //         {
    //             ip = std::make_unique<IntegrationPoint_>(* this, quadrature_index);
    //             quadrature_index ++;
    //         }
    //     }

    // public:
        
    //     explicit
    //     ElementPotential(
    //         ElementLagrangian_ const & lag
    //     )
    //     :
    //     lag_(lag)
    //     {
    //         setIntegrationPoints();
    //     }

    //     FiniteElement_ const &
    //     getFiniteElement()
    //     const
    //     {
    //         return lag_.getFiniteElement();
    //     }

    //     ElementLagrangian_ const &
    //     getLag()
    //     const
    //     {
    //         return lag_;
    //     }

    //     std::array<std::unique_ptr<IntegrationPoint_>, getNumQuadraturePoints()> const &
    //     getIntegrationPoints()
    //     const
    //     {
    //         return integration_points_;
    //     }

    //     std::array<std::unique_ptr<IntegrationPoint_>, getNumQuadraturePoints()> &
    //     getIntegrationPoints()
    //     {
    //         return integration_points_;
    //     }

    //     void
    //     setStrainOperators()
    //     {
    //         auto set_i = [&] <Integer t_i = 0> (
    //             auto & t_set_i
    //         )
    //         constexpr mutable
    //         {
    //             auto constexpr mapping = PotentialTraits2<t_potential>::template getLinearOperator<t_i>();
    //             for (auto & ip : integration_points_)
    //             {
    //                 ip->template setStrainOperator<mapping>();
    //             }
    //             if constexpr (t_i < PotentialTraits2<t_potential>::getNumLinearOperators() - 1)
    //             {
    //                 t_set_i.template operator()<t_i + 1>(t_set_i);
    //             }
    //         };
    //         set_i(set_i);
    //     }

    //     void
    //     setStrains()
    //     {
    //         auto set_i = [&] <Integer t_i = 0> (
    //             auto & t_set_i
    //         )
    //         constexpr mutable
    //         {
    //             auto constexpr mapping = PotentialTraits2<t_potential>::template getLinearOperator<t_i>();
    //             for (auto & ip : integration_points_)
    //             {
    //                 ip->template setStrainVectorBlock<mapping>();
    //             }
    //             if constexpr (t_i < PotentialTraits2<t_potential>::getNumLinearOperators() - 1)
    //             {
    //                 t_set_i.template operator()<t_i + 1>(t_set_i);
    //             }
    //         };
    //         set_i(set_i);
    //     }

    //     void
    //     setMaterialProperty(
    //         std::basic_string<Character> && label,
    //         std::function<Real(Point const &)> && function
    //     )
    //     {
    //         for (auto & ip : integration_points_)
    //         {
    //             ip->setMaterialProperty(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
    //         }
    //     }

    //     void
    //     setExternalVariable(
    //         std::basic_string<Character> && label,
    //         std::function<Real(Point const &)> && function
    //     )
    //     {
    //         for (auto & ip : integration_points_)
    //         {
    //             ip->setExternalVariable(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
    //         }
    //     }

    //     void
    //     integrateConstitutiveEquation(
    //         std::atomic<Boolean> & output_handler
    //     )
    //     {
    //         for (auto & ip : integration_points_)
    //         {
    //             ip->integrate(output_handler);
    //         }
    //     }

    //     void
    //     reserve()
    //     {
    //         for (auto & ip : integration_points_)
    //         {
    //             ip->reserve();
    //         }            
    //     }

    //     void
    //     recover()
    //     {
    //         for (auto & ip : integration_points_)
    //         {
    //             ip->recover();
    //         }            
    //     }

    //     ElementLagrangian_ const & lag_;

    //     std::array<std::unique_ptr<IntegrationPoint_>, getNumQuadraturePoints()> integration_points_;

    // };
    
    // template<ShapeConcept auto t_element, MeshConcept auto t_domain, LagrangianConcept auto t_lag>
    // struct ElementLagrangian : AbstractElementLagrangian<t_element, t_domain>
    // {

    // private:

    //     static constexpr
    //     Integer
    //     getSystemSize()
    //     {
    //         return LagTraits<t_lag>::template getSize<t_element, t_domain>();
    //     }

    //     template<PotentialConcept auto t_potential>
    //     static constexpr
    //     Integer
    //     getPotentialIndex()
    //     {
    //         return LagTraits<t_lag>::template getPotentialIndex<t_potential>();
    //     }

    //     using Base_ = AbstractElementLagrangian<t_element, t_domain>;

    //     using FiniteElement_ = FiniteElement<t_element, t_domain>;

    //     template<PotentialConcept auto t_potential>
    //     using ElementPotential_ = ElementPotential<t_element, t_domain, t_lag, t_potential>;

    //     template<PotentialConcept auto t_potential>
    //     using ElementPotentialPtr_ = std::unique_ptr<ElementPotential_<t_potential>>;

    //     using ElementPotentials_ = utility::tuple_expansion_t<std::tuple, ElementPotentialPtr_, t_lag.getPotentials()>;

    // public:

    //     explicit
    //     ElementLagrangian(
    //         FiniteElement_ const & finite_element
    //     )
    //     :
    //     Base_(t_lag, finite_element)
    //     {}

    //     template<PotentialConcept auto t_potential>
    //     void
    //     setPotential()
    //     {
    //         std::get<getPotentialIndex<t_potential>()>(element_potentials_) = std::make_unique<ElementPotential_<t_potential>>(* this);
    //     }

    //     template<PotentialConcept auto t_potential>
    //     ElementPotential_<t_potential> const &
    //     getPotential()
    //     const
    //     {
    //         return * std::get<getPotentialIndex<t_potential>()>(element_potentials_);
    //     }

    //     template<PotentialConcept auto t_potential>
    //     ElementPotential_<t_potential> &
    //     getPotential()
    //     {
    //         return * std::get<getPotentialIndex<t_potential>()>(element_potentials_);
    //     }

    //     void
    //     setResidualVector(
    //         Real const & time
    //     )
    //     {
    //         this->residual_vector_ = this->getElementExternalForces(time) - this->getElementInternalForces();
    //     }
        
    //     void
    //     setResidualVector(
    //         DenseVectorConcept<Real> auto && vector
    //     )
    //     {
    //         this->residual_vector_ = std::make_unique<DenseVector<Real, getSystemSize()>>(std::forward<decltype(vector)>(vector));
    //     }

    //     DenseVector<Real, getSystemSize()> const &
    //     getResidualVector()
    //     const
    //     {
    //         return * this->residual_vector_;
    //     }

    //     void
    //     setJacobianMatrix(
    //         DenseMatrixConcept<Real> auto && matrix
    //     )
    //     {
    //         this->jacobian_matrix_ = std::make_unique<DenseMatrix<Real, getSystemSize(), getSystemSize()>>(std::forward<decltype(matrix)>(matrix));
    //     }

    //     DenseMatrix<Real, getSystemSize(), getSystemSize()> const &
    //     getJacobianMatrix()
    //     const
    //     {
    //         return * this->jacobian_matrix_;
    //     }

    //     ElementPotentials_ element_potentials_;

    //     std::unique_ptr<DenseVector<Real, getSystemSize()>> residual_vector_;

    //     std::unique_ptr<DenseMatrix<Real, getSystemSize(), getSystemSize()>> jacobian_matrix_;

    // };

    // template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    // struct AbstractElementLagrangian
    // {

    // private:

    //     template<LagrangianConcept auto t_lag>
    //     static constexpr
    //     Integer
    //     getSystemSize()
    //     {
    //         return LagTraits<t_lag>::template getSize<t_element, t_domain>();
    //     }

    //     using FiniteElement_ = FiniteElement<t_element, t_domain>;

    //     template<LagrangianConcept auto t_lag>
    //     using ElementLagrangian_ = ElementLagrangian<t_element, t_domain, t_lag>;

    // protected:

    //     AbstractElementLagrangian(
    //         LagrangianConcept auto const & lag,
    //         FiniteElement_ const & finite_element
    //     )
    //     :
    //     label_(lag.getLabel()),
    //     finite_element_(finite_element)
    //     {}

    // public:

    //     Label const &
    //     getLabel()
    //     const
    //     {
    //         return label_;
    //     }

    //     FiniteElement_ const &
    //     getFiniteElement()
    //     const
    //     {
    //         return finite_element_;
    //     }

    //     template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    //     void
    //     setPotential()
    //     {
    //         static_cast<ElementLagrangian_<t_lag> *>(this)->template setPotential<t_potential>();
    //     }

    //     template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    //     auto const &
    //     getPotential()
    //     const
    //     {
    //         return static_cast<ElementLagrangian_<t_lag> const *>(this)->template getPotential<t_potential>();
    //     }

    //     template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    //     auto &
    //     getPotential()
    //     {
    //         return static_cast<ElementLagrangian_<t_lag> *>(this)->template getPotential<t_potential>();
    //     }

    //     template<LagrangianConcept auto t_lag>
    //     void
    //     setJacobianMatrix(
    //         DenseMatrixConcept<Real> auto && matrix
    //     )
    //     {
    //         static_cast<ElementLagrangian_<t_lag> *>(this)->setJacobianMatrix(std::forward<decltype(matrix)>(matrix));
    //     }

    //     template<LagrangianConcept auto t_lag>
    //     void
    //     setJacobianMatrix()
    //     {
    //         static_cast<ElementLagrangian_<t_lag> *>(this)->setJacobianMatrix(getElementJacobianMatrix<t_lag>());
    //     }

    //     template<LagrangianConcept auto t_lag>
    //     auto const &
    //     getJacobianMatrix()
    //     const
    //     {
    //         return static_cast<ElementLagrangian_<t_lag> const *>(this)->getJacobianMatrix();
    //     }

    //     template<LagrangianConcept auto t_lag>
    //     void
    //     setResidualVector(
    //         DenseVectorConcept<Real> auto && vector
    //     )
    //     {
    //         static_cast<ElementLagrangian_<t_lag> *>(this)->setResidualVector(std::forward<decltype(vector)>(vector));
    //     }

    //     template<LagrangianConcept auto t_lag>
    //     void
    //     setResidualVector()
    //     {
    //         static_cast<ElementLagrangian_<t_lag> *>(this)->setResidualVector(getElementInternalForces<t_lag>());
    //     }

    //     template<LagrangianConcept auto t_lag>
    //     auto const &
    //     getResidualVector()
    //     const
    //     {
    //         return static_cast<ElementLagrangian_<t_lag> const *>(this)->getResidualVector();
    //     }

    //     template<LagrangianConcept auto t_lag>
    //     DenseVector<Real, getSystemSize<t_lag>()>
    //     getElementInternalForces()
    //     const
    //     {
    //         auto internal_forces = DenseVector<Real, getSystemSize<t_lag>()>();
    //         internal_forces.setZero();
    //         auto set_i = [&] <Integer t_i = 0, Integer t_j = 0> (
    //             auto & t_set_i
    //         )
    //         constexpr mutable
    //         {
    //             auto constexpr potential = LagTraits<t_lag>::template getPotential<t_i>();
    //             auto constexpr j_mapping = PotentialTraits2<potential>::template getLinearOperator<t_j>();
    //             auto constexpr j_field = j_mapping.getField();
    //             auto constexpr size_j = DiscretizationTraits<j_field>::template getTensorSpaceSize<t_element, t_domain>();
    //             auto constexpr offset_j = LagTraits<t_lag>::template getFieldOffset<t_element, t_domain, j_field>();
    //             for (auto const & ip : this->getPotential<t_lag, potential>().getIntegrationPoints())
    //             {
    //                 auto const & mat0 = ip->template getStrainOperator<j_mapping>();
    //                 auto const mat = ip->template getStressVectorBlock<j_mapping>();
    //                 internal_forces.template segment<size_j>(offset_j) += ip->getCurrentWeight() * mat0.transpose() * mat;
    //             }
    //             if constexpr (t_j < PotentialTraits2<potential>::getNumLinearOperators() - 1)
    //             {
    //                 t_set_i.template operator()<t_i, t_j + 1>(t_set_i);
    //             }
    //             else if constexpr (t_i < LagTraits<t_lag>::getNumPotentials() - 1)
    //             {
    //                 t_set_i.template operator()<t_i + 1, 0>(t_set_i);
    //             }
    //         };
    //         set_i(set_i);
    //         return internal_forces;
    //     }

    //     template<LagrangianConcept auto t_lag>
    //     DenseVector<Real, getSystemSize<t_lag>()>
    //     getElementExternalForces(
    //         Real const & time
    //     )
    //     const
    //     {
    //         auto constexpr quadrature = GaussQuadrature(4);
    //         auto external_forces = DenseVector<Real, getSystemSize<t_lag>()>();
    //         external_forces.setZero();
    //         auto set_i = [&] <Integer t_i = 0> (
    //             auto & t_set_i
    //         )
    //         constexpr mutable
    //         {
    //             auto constexpr field = LagTraits<t_lag>::template getField<t_i>();
    //             auto constexpr size_j = DiscretizationTraits<field>::template getTensorSpaceSize<t_element, t_domain>();
    //             auto constexpr offset_j = LagTraits<t_lag>::template getFieldOffset<t_element, t_domain, field>();
    //             for (auto const & domain : this->getFiniteElement().getDomains())
    //             {
    //                 if (domain->template hasDiscreteField<field>())
    //                 {
    //                     if (domain->template getDiscreteField<field>().hasLoads())
    //                     {
    //                         for (auto const & load : domain->template getDiscreteField<field>().getLoads())
    //                         {
    //                             for (auto i = 0; i < QuadratureTraits<quadrature>::template Rule<t_element>::getSize(); i++)
    //                             {
    //                                 auto weight = this->getFiniteElement().template getCurrentQuadratureWeight<quadrature>(i);
    //                                 auto reference_point = this->getFiniteElement().template getReferenceQuadraturePoint<quadrature>(i);
    //                                 auto current_point = this->getFiniteElement().template getCurrentQuadraturePoint<quadrature>(i);
    //                                 auto vector = this->getFiniteElement().template getFieldDualVector<field>(reference_point, load.getRow(), load.getCol());
    //                                 external_forces.template segment<size_j>(offset_j) += weight * load.getValue(current_point, time) * vector;
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //             if constexpr (t_i < LagTraits<t_lag>::getNumFields() - 1)
    //             {
    //                 t_set_i.template operator()<t_i + 1>(t_set_i);
    //             }
    //         };
    //         set_i(set_i);
    //         return external_forces;
    //     }
        
    //     template<LagrangianConcept auto t_lag>
    //     DenseMatrix<Real, getSystemSize<t_lag>(), getSystemSize<t_lag>()>
    //     getElementJacobianMatrix()
    //     const
    //     {
    //         auto jacobian_matrix = DenseMatrix<Real, getSystemSize<t_lag>(), getSystemSize<t_lag>()>();
    //         jacobian_matrix.setZero();
    //         auto set_i = [&] <Integer t_i = 0, Integer t_j = 0, Integer t_k = 0> (
    //             auto & t_set_i
    //         )
    //         constexpr mutable
    //         {
    //             auto constexpr potential = LagTraits<t_lag>::template getPotential<t_i>();
    //             auto constexpr j_mapping = PotentialTraits2<potential>::template getLinearOperator<t_j>();
    //             auto constexpr j_field = j_mapping.getField();
    //             auto constexpr k_mapping = PotentialTraits2<potential>::template getLinearOperator<t_k>();
    //             auto constexpr k_field = k_mapping.getField();
    //             auto constexpr size_j = DiscretizationTraits<j_field>::template getTensorSpaceSize<t_element, t_domain>();
    //             auto constexpr size_k = DiscretizationTraits<k_field>::template getTensorSpaceSize<t_element, t_domain>();
    //             auto constexpr offset_j = LagTraits<t_lag>::template getFieldOffset<t_element, t_domain, j_field>();
    //             auto constexpr offset_k = LagTraits<t_lag>::template getFieldOffset<t_element, t_domain, k_field>();
    //             auto const & frm = this->template getPotential<t_lag, potential>();
    //             for (auto const & ip : frm.getIntegrationPoints())
    //             {
    //                 auto const & mat0 = ip->template getStrainOperator<j_mapping>();
    //                 auto const & mat1 = ip->template getStrainOperator<k_mapping>();
    //                 auto const mat = ip->template getJacobianMatrixBlock<j_mapping, k_mapping>();
    //                 jacobian_matrix.template block<size_j, size_k>(offset_j, offset_k) += ip->getCurrentWeight() * mat0.transpose() * mat * mat1;
    //             }
    //             if constexpr (t_k < PotentialTraits2<potential>::getNumLinearOperators() - 1)
    //             {
    //                 t_set_i.template operator()<t_i, t_j, t_k + 1>(t_set_i);
    //             }
    //             if constexpr (t_j < PotentialTraits2<potential>::getNumLinearOperators() - 1)
    //             {
    //                 t_set_i.template operator()<t_i, t_j + 1, 0>(t_set_i);
    //             }
    //             else if constexpr (t_i < LagTraits<t_lag>::getNumPotentials() - 1)
    //             {
    //                 t_set_i.template operator()<t_i + 1, 0, 0>(t_set_i);
    //             }
    //         };
    //         set_i(set_i);
    //         // std::cout << "jacobian_matrix : " << std::endl;
    //         // std::cout << jacobian_matrix << std::endl;
    //         return jacobian_matrix;
    //     }

    // protected:

    //     Label const & label_;

    //     FiniteElement_ const & finite_element_;

    // };
    
    /**
     * *********************************************************************************************************************************************************
     */

    // template<DomainConcept auto t_dim, MeshConcept auto t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    // struct DomainPotential
    // {

    // private:

    //     using DomainLagrangian_ = DomainLagrangian<t_dim, t_domain, t_lag>;

    // public:

    //     DomainPotential(
    //         DomainLagrangian_ const & lag,
    //         auto const &... args
    //     )
    //     :
    //     lag_(lag),
    //     behavior_(mgis::behaviour::load(args...))
    //     {}

    //     DomainLagrangian_ const &
    //     getLag()
    //     const
    //     {
    //         return lag_;
    //     }

    //     mgis::behaviour::Behaviour const &
    //     getMgisBehavior()
    //     const
    //     {
    //         return behavior_;
    //     }

    // private:

    //     DomainLagrangian_ const & lag_;
        
    //     mgis::behaviour::Behaviour behavior_;

    // };
    
    // template<DomainConcept auto t_dim, MeshConcept auto t_domain, LagrangianConcept auto t_lag>
    // struct DomainLagrangian : AbstractDomainLagrangian<t_dim, t_domain>
    // {

    // private:

    //     template<PotentialConcept auto t_potential>
    //     static constexpr
    //     Integer
    //     getPotentialIndex()
    //     {
    //         return LagTraits<t_lag>::template getPotentialIndex<t_potential>();
    //     }

    //     using Base_ = AbstractDomainLagrangian<t_dim, t_domain>;

    //     using Domain_ = FiniteDomain<t_dim, t_domain>;

    //     template<PotentialConcept auto t_potential>
    //     using DomainPotential_ = DomainPotential<t_dim, t_domain, t_lag, t_potential>;

    //     template<PotentialConcept auto t_potential>
    //     using DomainPotentialPtr_ = std::unique_ptr<DomainPotential_<t_potential>>;

    //     using DomainPotentials_ = utility::tuple_expansion_t<std::tuple, DomainPotentialPtr_, t_lag.getPotentials()>;

    // public:

    //     explicit
    //     DomainLagrangian(
    //         Domain_ const & domain
    //     )
    //     :
    //     Base_(t_lag),
    //     domain_(domain)
    //     {}

    //     Domain_ const &
    //     getFiniteDomain()
    //     const
    //     {
    //         return domain_;
    //     }

    //     template<PotentialConcept auto t_potential>
    //     void
    //     setPotential(
    //         auto const &... args
    //     )
    //     {
    //         std::get<getPotentialIndex<t_potential>()>(domain_potentials_) = std::make_unique<DomainPotential_<t_potential>>(* this, args...);
    //     }

    //     template<PotentialConcept auto t_potential>
    //     DomainPotential_<t_potential> const &
    //     getPotential()
    //     const
    //     {
    //         return * std::get<getPotentialIndex<t_potential>()>(domain_potentials_);
    //     }

    //     template<PotentialConcept auto t_potential>
    //     DomainPotential_<t_potential> &
    //     getPotential()
    //     {
    //         return * std::get<getPotentialIndex<t_potential>()>(domain_potentials_);
    //     }

    // private:

    //     Domain_ const & domain_;

    //     DomainPotentials_ domain_potentials_;

    // };

    // template<DomainConcept auto t_dim, MeshConcept auto t_domain>
    // struct AbstractDomainLagrangian
    // {

    // private:

    //     template<LagrangianConcept auto t_lag>
    //     using DomainLagrangian_ = DomainLagrangian<t_dim, t_domain, t_lag>;

    // protected:

    //     explicit
    //     AbstractDomainLagrangian(
    //         LagrangianConcept auto const & lag
    //     )
    //     :
    //     label_(lag.getLabel())
    //     {}

    // public:

    //     Label const &
    //     getLabel()
    //     const
    //     {
    //         return label_;
    //     }

    //     template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    //     void
    //     setPotential(
    //         auto const &... args
    //     )
    //     {
    //         static_cast<DomainLagrangian_<t_lag> *>(this)->template setPotential<t_potential>(args...);
    //     }

    //     template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    //     auto const &
    //     getPotential()
    //     const
    //     {
    //         return static_cast<DomainLagrangian_<t_lag> const *>(this)->template getPotential<t_potential>();
    //     }

    //     template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    //     auto &
    //     getPotential()
    //     {
    //         return static_cast<DomainLagrangian_<t_lag> *>(this)->template getPotential<t_potential>();
    //     }

    // protected:

    //     Label const & label_;

    // };

} // namespace lolita::core

#endif /* A0EFE384_4F56_4FAA_84B9_65C5B2053C0D */
