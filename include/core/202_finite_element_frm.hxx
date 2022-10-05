#ifndef B53096ED_E28A_48CC_B15B_29B622D775A3
#define B53096ED_E28A_48CC_B15B_29B622D775A3

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/100_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/200_quadrature.hxx"
#include "core/201_finite_element_dof.hxx"

namespace lolita
{

    template<Element t_element, Domain t_domain>
    struct FiniteElement;

    template<Integer t_dim, Domain t_domain>
    struct FiniteDomain;

    template<Integer t_dim, Domain t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct DomainPotential;
    
    template<Element t_element, Domain t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct ElementPotential;

    template<Integer t_dim, Domain t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct DomainPotential;
    
    template<Integer t_dim, Domain t_domain>
    struct AbstractDomainLagrangian;

    template<Element t_element, Domain t_domain>
    struct AbstractElementLagrangian;
    
    template<Element t_element, Domain t_domain, LagrangianConcept auto t_lag>
    struct ElementLagrangian;
    
    template<Integer t_dim, Domain t_domain, LagrangianConcept auto t_lag>
    struct DomainLagrangian;

    /**
     * *********************************************************************************************************************************************************
     */

    template<Domain t_domain>
    struct IntegrationPoint
    {

        IntegrationPoint(
            PointConcept auto const & reference_point,
            PointConcept auto const & current_point,
            Real weight,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        :
        reference_coordinates_(reference_point),
        coordinates_(current_point),
        weight_(weight),
        behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behavior)))
        {}

        IntegrationPoint(
            PointConcept auto && reference_point,
            PointConcept auto && current_point,
            Real weight,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        :
        reference_coordinates_(std::move(reference_point)),
        coordinates_(std::move(current_point)),
        weight_(weight),
        behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behavior)))
        {}
    
        Boolean
        operator==(
            IntegrationPoint const & other
        )
        const = default;
        
        Boolean
        operator!=(
            IntegrationPoint const & other
        )
        const = default;

        void
        integrate(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior,
            std::atomic<Boolean> & output_handler
        )
        {
            behavior_data_->K[0] = 4;
            auto behaviour_data_view = mgis::behaviour::make_view(* behavior_data_);
            auto res = mgis::behaviour::integrate(behaviour_data_view, * behavior);
            if (res < 1)
            {
                output_handler = false;
            }
        }

        void
        reserve()
        {
            mgis::behaviour::update(* behavior_data_);
        }

        void
        recover()
        {
            mgis::behaviour::revert(* behavior_data_);
        }

        void
        setMaterialProperty(
            std::basic_string<Character> && material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto value = std::forward<std::function<Real(Point const &)>>(function)(coordinates_);
            mgis::behaviour::setMaterialProperty(behavior_data_->s0, std::forward<std::basic_string<Character>>(material_property_label), value);
            mgis::behaviour::setMaterialProperty(behavior_data_->s1, std::forward<std::basic_string<Character>>(material_property_label), value);
        }

        void
        setExternalVariable(
            std::basic_string<Character> && material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto value = std::forward<std::function<Real(Point const &)>>(function)(coordinates_);
            mgis::behaviour::setExternalStateVariable(behavior_data_->s0, std::forward<std::basic_string<Character>>(material_property_label), value);
            mgis::behaviour::setExternalStateVariable(behavior_data_->s1, std::forward<std::basic_string<Character>>(material_property_label), value);
        }

        // template<PotentialConcept auto t_potential>
        // algebra::View<typename PotentialTraits<t_potential>::template JacobianMatrix<t_domain> const>
        // getJacobianMatrixView()
        // const
        // {
        //     return algebra::View<typename PotentialTraits<t_potential>::template JacobianMatrix<t_domain> const>(behavior_data_->K.data());
        // }

        template<PotentialConcept auto t_potential, MappingConcept auto t_strain1, MappingConcept auto t_strain2>
        algebra::View<DenseMatrix<Real, MappingTraits<t_strain1>::template getSize<t_domain>(), MappingTraits<t_strain2>::template getSize<t_domain>()> const>
        getJacobianMatrixBlock()
        const
        {
            auto constexpr jacobian_size = PotentialTraits<t_potential>::template getSize<t_domain>();
            auto constexpr jacobian_block_num_rows = MappingTraits<t_strain1>::template getSize<t_domain>();
            auto constexpr jacobian_block_num_cols = MappingTraits<t_strain2>::template getSize<t_domain>();
            auto jacobian_block_row_offset = PotentialTraits<t_potential>::template getOffset<t_domain, t_strain1>();
            auto jacobian_block_col_offset = PotentialTraits<t_potential>::template getOffset<t_domain, t_strain2>();
            auto block_offset = jacobian_size * jacobian_block_row_offset + jacobian_block_col_offset;
            return algebra::View<DenseMatrix<Real, jacobian_block_num_rows, jacobian_block_num_cols> const>(behavior_data_->K.data() + block_offset);
        }

        template<PotentialConcept auto t_potential, MappingConcept auto t_strain>
        algebra::View<DenseVector<Real, MappingTraits<t_strain>::template getSize<t_domain>()> const>
        getStressVectorBlock()
        const
        {
            auto constexpr jacobian_block_num_cols = MappingTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits<t_potential>::template getOffset<t_domain, t_strain>();
            return algebra::View<DenseVector<Real, jacobian_block_num_cols> const>(behavior_data_->s1.thermodynamic_forces.data() + jacobian_block_col_offset);
        }

        template<PotentialConcept auto t_potential, MappingConcept auto t_strain>
        algebra::View<DenseVector<Real, MappingTraits<t_strain>::template getSize<t_domain>()> const>
        getStrainVectorBlock()
        const
        {
            auto constexpr jacobian_block_num_cols = MappingTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits<t_potential>::template getOffset<t_domain, t_strain>();
            return algebra::View<DenseVector<Real, jacobian_block_num_cols> const>(behavior_data_->s1.gradients.data() + jacobian_block_col_offset);
        }

        template<PotentialConcept auto t_potential, MappingConcept auto t_strain>
        void
        setStrainVectorBlock(
            DenseVectorConcept<Real> auto && strain
        )
        {
            auto constexpr jacobian_block_num_cols = MappingTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits<t_potential>::template getOffset<t_domain, t_strain>();
            auto strain_block = algebra::View<DenseVector<Real, jacobian_block_num_cols>>(behavior_data_->s1.gradients.data() + jacobian_block_col_offset);
            strain_block = std::forward<decltype(strain)>(strain);
        }

        template<PotentialConcept auto t_potential, MappingConcept auto t_strain>
        algebra::View<DenseVector<Real, MappingTraits<t_strain>::template getSize<t_domain>()> const>
        getResidualVectorBlock()
        const
        {
            return this->getStressVectorBlock<t_potential, t_strain>();
        }
        
        Point const &
        getCurrentCoordinates()
        const
        {
            return coordinates_;
        }
        
        algebra::View<Point const>
        getReferenceCoordinates()
        const
        {
            return reference_coordinates_;
        }
        
        Real const &
        getCurrentWeight()
        const
        {
            return weight_;
        }
        
        Real const &
        getReferenceWeight()
        const
        {
            return weight_;
        }

        template<MappingConcept auto t_strain>
        void
        addStrainOperator(
            DenseMatrixConcept<Real> auto && input
        )
        {
            if (strain_matrix_list_ == nullptr)
            {
                strain_matrix_list_ = std::make_unique<std::vector<ElementaryOperator<Label, DenseMatrix<Real>>>>();
            }
            for (auto const & m : * strain_matrix_list_)
            {
                if (m.getTag() == t_strain.getLabel())
                {
                    return;
                }
            }
            strain_matrix_list_->push_back(ElementaryOperator<Label, DenseMatrix<Real>>(t_strain.getLabel(), std::forward<decltype(input)>(input)));
        }

        template<MappingConcept auto t_strain>
        DenseMatrix<Real> const &
        getStrainOperator()
        const
        {
            for (auto const & m : * strain_matrix_list_)
            {
                if (m.getTag() == t_strain.getLabel())
                {
                    return m.getOperator();
                }
            }
            throw std::runtime_error("No usch thing");
        }

        template<MappingConcept auto t_strain>
        Boolean
        hasStrainOperator()
        const
        {
            for (auto const & m : * strain_matrix_list_)
            {
                if (m.getTag() == t_strain.getLabel())
                {
                    return true;
                }
            }
            return false;
        }

        // template<FieldConcept auto t_field>
        // void
        // addDiscreteField()
        // {
        //     if (ptr_data_ == nullptr)
        //     {
        //         ptr_data_ = std::make_unique<std::vector<ElementDiscreteField<t_domain>>>();
        //     }
        //     for (auto const & item : * ptr_data_)
        //     {
        //         if (item.getLabel() == t_field.getLabel())
        //         {
        //             return;
        //         }
        //     }
        //     ptr_data_->push_back(ElementDiscreteField<t_domain>(t_field));
        // }

        // template<FieldConcept auto t_field>
        // ElementDiscreteField<t_domain> const &
        // getDiscreteField()
        // const
        // {
        //     if (ptr_data_ == nullptr)
        //     {
        //         throw std::runtime_error("Empty");
        //     }
        //     else
        //     {
        //         for (auto const & item : * ptr_data_)
        //         {
        //             if (item.getLabel() == t_field.getLabel())
        //             {
        //                 return item;
        //             }
        //         }
        //         throw std::runtime_error("No such field data");
        //     }
        // }

        // template<FieldConcept auto t_field>
        // ElementDiscreteField<t_domain> &
        // getDiscreteField()
        // {
        //     if (ptr_data_ == nullptr)
        //     {
        //         throw std::runtime_error("Empty");
        //     }
        //     else
        //     {
        //         for (auto & item : * ptr_data_)
        //         {
        //             if (item.getLabel() == t_field.getLabel())
        //             {
        //                 return item;
        //             }
        //         }
        //         throw std::runtime_error("No such field data");
        //     }
        // }

        // template<FieldConcept auto t_field, Basis t_basis>
        // void
        // addDiscreteFieldDegreeOfFreedom(
        //     auto const &... args
        // )
        // {
        //     this->template getDiscreteField<t_field>().template addDegreeOfFreedom<t_field, t_basis>(args...);
        // }

        algebra::View<Point const> reference_coordinates_;
    
        Point coordinates_;
        
        Real weight_;

        std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

        std::unique_ptr<std::vector<ElementaryOperator<Label, DenseMatrix<Real>>>> strain_matrix_list_;

        // std::unique_ptr<std::vector<ElementDiscreteField<t_domain>>> ptr_data_;

    };
    
    template<Domain t_domain>
    struct ElementFormulation
    {

    private:

        using t_IntegrationPoint = IntegrationPoint<t_domain>;

    public:

        explicit
        ElementFormulation(
            PotentialConcept auto const & potential
        )
        :
        label_(potential.getLabel()),
        integration_points_(),
        behavior_()
        {}
        
        ElementFormulation(
            PotentialConcept auto const & potential,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        :
        label_(potential.getLabel()),
        integration_points_(),
        behavior_(behavior)
        {}

        void
        addIntegrationPoint(
            PointConcept auto && reference_point,
            PointConcept auto && current_point,
            Real weight,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            if (integration_points_ == nullptr)
            {
                integration_points_ = std::make_unique<std::vector<t_IntegrationPoint>>();
            }
            integration_points_->push_back(t_IntegrationPoint(std::forward<decltype(reference_point)>(reference_point), std::forward<decltype(current_point)>(current_point), weight, behavior));
        }

        std::vector<t_IntegrationPoint> const &
        getIntegrationPoints()
        const
        {
            return * integration_points_;
        }

        std::vector<t_IntegrationPoint> &
        getIntegrationPoints()
        {
            return * integration_points_;
        }
        
        Boolean
        operator==(
            ElementFormulation const & other
        )
        const = default;
        
        Boolean
        operator!=(
            ElementFormulation const & other
        )
        const = default;

        Label const &
        getLabel()
        const
        {
            return label_;
        }

        void
        setMaterialProperty(
            std::basic_string<Character> && label,
            std::function<Real(Point const &)> && function
        )
        {
            for (auto & ip : * integration_points_)
            {
                ip.setMaterialProperty(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
            }
        }

        void
        setExternalVariable(
            std::basic_string<Character> && label,
            std::function<Real(Point const &)> && function
        )
        {
            for (auto & ip : * integration_points_)
            {
                ip.setExternalVariable(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
            }
        }

        void
        integrateConstitutiveEquation(
            std::atomic<Boolean> & output_handler
        )
        {
            for (auto & ip : * integration_points_)
            {
                ip.integrate(behavior_, output_handler);
            }
        }

        void
        reserve()
        {
            for (auto & ip : * integration_points_)
            {
                ip.reserve();
            }            
        }

        void
        recover()
        {
            for (auto & ip : * integration_points_)
            {
                ip.recover();
            }            
        }

        // template<Integer t_size>
        // void
        // setSystem(
        //     DenseMatrixConcept<Real, t_size, t_size> auto && matrix,
        //     DenseVectorConcept<Real, t_size> auto && vector
        // )
        // {
        //     system_ = std::make_unique<DenseSystem<Real, t_size>>(std::forward<decltype(matrix)>(matrix), std::forward<decltype(vector)>(vector));
        // }

        Label const & label_;

        std::unique_ptr<std::vector<t_IntegrationPoint>> integration_points_;

        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        std::unique_ptr<StaticCondensation> condensation_;

        // std::unique_ptr<DenseSystemBase> system_;

        // DenseMatrix<Real> jacobian_matrix_;

        // DenseVector<Real> residual_vector_;

    };

    /**
     * *********************************************************************************************************************************************************
     */

    template<Integer t_dim, Domain t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
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

    template<Element t_element, Domain t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
    struct IntegrationPoint2
    {

    private:

        template<MappingConcept auto t_strain>
        static constexpr
        Integer
        getStrainOperatorNumRows()
        {
            return MappingTraits<t_strain>::template getSize<t_domain>();
        }

        template<MappingConcept auto t_strain>
        static constexpr
        Integer
        getStrainOperatorNumCols()
        {
            return FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<t_strain.getField()>();
        }

        template<MappingConcept auto t_strain>
        static constexpr
        Integer
        getStrainIndex()
        {
            return PotentialTraits2<t_potential>::template getIndex<t_strain>();
        }

        using FiniteElement_ = FiniteElement<t_element, t_domain>;

        using ElementPotential_ = ElementPotential<t_element, t_domain, t_lag, t_potential>;

        using BehaviorData_ = BehaviourData<t_element.getDim(), t_domain, t_lag, t_potential>;

        template<MappingConcept auto t_strain>
        using StrainOperator_ = DenseMatrix<Real, getStrainOperatorNumRows<t_strain>(), getStrainOperatorNumCols<t_strain>()>;

        template<MappingConcept auto t_strain>
        using StrainOperatorPtr_ = std::unique_ptr<StrainOperator_<t_strain>>;

        using StrainOperators_ = utility::tuple_expansion_t<std::tuple, StrainOperatorPtr_, t_potential.getStrains()>;

    public:
        
        IntegrationPoint2(
            ElementPotential_ const & potential,
            Integer index
        )
        :
        potential_(potential),
        index_(index),
        behavior_data_(potential.getFiniteElement().getDomain().template getLagrangian<t_lag>().template getPotential<t_lag, t_potential>())
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
            return this->getFiniteElement().getDomain().template getLagrangian<t_lag>().template getPotential<t_lag, t_potential>();
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

        template<MappingConcept auto t_strain1, MappingConcept auto t_strain2>
        algebra::View<DenseMatrix<Real, MappingTraits<t_strain1>::template getSize<t_domain>(), MappingTraits<t_strain2>::template getSize<t_domain>()> const>
        getJacobianMatrixBlock()
        const
        {
            auto constexpr jacobian_size = PotentialTraits2<t_potential>::template getSize<t_domain>();
            auto constexpr jacobian_block_num_rows = MappingTraits<t_strain1>::template getSize<t_domain>();
            auto constexpr jacobian_block_num_cols = MappingTraits<t_strain2>::template getSize<t_domain>();
            auto jacobian_block_row_offset = PotentialTraits2<t_potential>::template getOffset<t_domain, t_strain1>();
            auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getOffset<t_domain, t_strain2>();
            auto block_offset = jacobian_size * jacobian_block_row_offset + jacobian_block_col_offset;
            auto const * mgis_data = this->getMgisBhv().K.data();
            return algebra::View<DenseMatrix<Real, jacobian_block_num_rows, jacobian_block_num_cols> const>(mgis_data + block_offset);
        }

        template<MappingConcept auto t_strain>
        algebra::View<DenseVector<Real, MappingTraits<t_strain>::template getSize<t_domain>()> const>
        getStressVectorBlock()
        const
        {
            auto constexpr jacobian_block_num_cols = MappingTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getOffset<t_domain, t_strain>();
            auto const * mgis_data = this->getMgisBhv().s1.thermodynamic_forces.data();
            return algebra::View<DenseVector<Real, jacobian_block_num_cols> const>(mgis_data + jacobian_block_col_offset);
        }

        template<MappingConcept auto t_strain>
        algebra::View<DenseVector<Real, MappingTraits<t_strain>::template getSize<t_domain>()> const>
        getStrainVectorBlock()
        const
        {
            auto constexpr jacobian_block_num_cols = MappingTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getOffset<t_domain, t_strain>();
            auto const * mgis_data = this->getMgisBhv().s1.gradients.data();
            return algebra::View<DenseVector<Real, jacobian_block_num_cols> const>(mgis_data + jacobian_block_col_offset);
        }

        template<MappingConcept auto t_strain>
        void
        setStrainVectorBlock()
        {
            auto constexpr jacobian_block_num_cols = MappingTraits<t_strain>::template getSize<t_domain>();
            auto jacobian_block_col_offset = PotentialTraits2<t_potential>::template getOffset<t_domain, t_strain>();
            auto * mgis_data = this->getMgisBhv().s1.gradients.data();
            auto strain_block = algebra::View<DenseVector<Real, jacobian_block_num_cols>>(mgis_data + jacobian_block_col_offset);
            if (this->template hasStrainOperator<t_strain>())
            {
                strain_block = this->template getStrainOperator<t_strain>() * this->getFiniteElement().template getUnknownCoefficients<t_strain.getField()>();
            }
            else
            {
                strain_block = this->template letStrainOperator<t_strain>() * this->getFiniteElement().template getUnknownCoefficients<t_strain.getField()>();
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

        template<MappingConcept auto t_strain>
        Boolean
        hasStrainOperator()
        const
        {
            return std::get<getStrainIndex<t_strain>()>(strain_operators_) != nullptr;
        }
        
        template<MappingConcept auto t_strain>
        StrainOperator_<t_strain> const &
        getStrainOperator()
        const
        {
            return * std::get<getStrainIndex<t_strain>()>(strain_operators_);
        }

        template<MappingConcept auto t_strain>
        StrainOperator_<t_strain>
        letStrainOperator()
        const
        {
            return getFiniteElement().template getMapping<t_strain>(this->getReferenceCoordinates());
        }

        template<MappingConcept auto t_strain>
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
    
    template<Element t_element, Domain t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
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
                auto constexpr mapping = PotentialTraits2<t_potential>::template getStrain<t_i>();
                for (auto & ip : integration_points_)
                {
                    ip->template setStrainOperator<mapping>();
                }
                if constexpr (t_i < PotentialTraits2<t_potential>::getNumMappings() - 1)
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
                auto constexpr mapping = PotentialTraits2<t_potential>::template getStrain<t_i>();
                for (auto & ip : integration_points_)
                {
                    ip->template setStrainVectorBlock<mapping>();
                }
                if constexpr (t_i < PotentialTraits2<t_potential>::getNumMappings() - 1)
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
    
    template<Element t_element, Domain t_domain, LagrangianConcept auto t_lag>
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
            return LagTraits<t_lag>::template getIndex<t_potential>();
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
        Base_(t_lag),
        finite_element_(finite_element)
        {}

        FiniteElement_ const &
        getFiniteElement()
        const
        {
            return finite_element_;
        }

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

        DenseVector<Real, getSystemSize()>
        getElementInternalForces()
        const
        {
            auto internal_forces = DenseVector<Real, getSystemSize()>();
            internal_forces.setZero();
            auto set_i = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = LagTraits<t_lag>::template getPotential<t_i>();
                auto constexpr j_mapping = PotentialTraits2<potential>::template getStrain<t_j>();
                auto constexpr j_field = j_mapping.getField();
                // --> making block
                auto constexpr size_j = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<j_field>();
                auto constexpr offset_j = LagTraits<t_lag>::template getOffset<t_element, t_domain, j_field>();
                for (auto const & ip : this->getPotential<potential>().getIntegrationPoints())
                {
                    auto const & mat0 = ip->template getStrainOperator<j_mapping>();
                    auto const mat = ip->template getStressVectorBlock<j_mapping>();
                    internal_forces.template segment<size_j>(offset_j) += ip->getCurrentWeight() * mat0.transpose() * mat;
                }
                if constexpr (t_j < PotentialTraits2<potential>::getNumMappings() - 1)
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

        DenseVector<Real, getSystemSize()>
        getElementExternalForces(
            Real const & time
        )
        const
        {
            auto constexpr quadrature = Quadrature("Gauss", 4);
            auto external_forces = DenseVector<Real, getSystemSize()>();
            external_forces.setZero();
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr field = LagTraits<t_lag>::template getField<t_i>();
                auto constexpr size_j = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<field>();
                auto constexpr offset_j = LagTraits<t_lag>::template getOffset<t_element, t_domain, field>();
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
        
        DenseMatrix<Real, getSystemSize(), getSystemSize()>
        getElementJacobianMatrix()
        const
        {
            auto jacobian_matrix = DenseMatrix<Real, getSystemSize(), getSystemSize()>();
            jacobian_matrix.setZero();
            auto set_i = [&] <Integer t_i = 0, Integer t_j = 0, Integer t_k = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = LagTraits<t_lag>::template getPotential<t_i>();
                auto constexpr j_mapping = PotentialTraits2<potential>::template getStrain<t_j>();
                auto constexpr j_field = j_mapping.getField();
                auto constexpr k_mapping = PotentialTraits2<potential>::template getStrain<t_k>();
                auto constexpr k_field = k_mapping.getField();
                auto constexpr size_j = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<j_field>();
                auto constexpr size_k = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<k_field>();
                auto constexpr offset_j = LagTraits<t_lag>::template getOffset<t_element, t_domain, j_field>();
                auto constexpr offset_k = LagTraits<t_lag>::template getOffset<t_element, t_domain, k_field>();
                auto const & frm = this->template getPotential<potential>();
                for (auto const & ip : frm.getIntegrationPoints())
                {
                    auto const & mat0 = ip->template getStrainOperator<j_mapping>();
                    auto const & mat1 = ip->template getStrainOperator<k_mapping>();
                    auto const mat = ip->template getJacobianMatrixBlock<potential, j_mapping, k_mapping>();
                    jacobian_matrix.template block<size_j, size_k>(offset_j, offset_k) += ip->getCurrentWeight() * mat0.transpose() * mat * mat1;
                }
                if constexpr (t_k < PotentialTraits2<potential>::getNumMappings() - 1)
                {
                    t_set_i.template operator()<t_i, t_j, t_k + 1>(t_set_i);
                }
                if constexpr (t_j < PotentialTraits2<potential>::getNumMappings() - 1)
                {
                    t_set_i.template operator()<t_i, t_j + 1, 0>(t_set_i);
                }
                else if constexpr (t_i < LagTraits<t_lag>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1, 0, 0>(t_set_i);
                }
            };
            set_i(set_i);
            return jacobian_matrix;
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

        FiniteElement_ const & finite_element_;

        ElementPotentials_ element_potentials_;

        std::unique_ptr<DenseVector<Real, getSystemSize()>> residual_vector_;

        std::unique_ptr<DenseMatrix<Real, getSystemSize(), getSystemSize()>> jacobian_matrix_;

    };

    template<Element t_element, Domain t_domain>
    struct AbstractElementLagrangian
    {

    private:

        template<LagrangianConcept auto t_lag>
        using ElementLagrangian_ = ElementLagrangian<t_element, t_domain, t_lag>;

    protected:

        explicit
        AbstractElementLagrangian(
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
            return static_cast<ElementLagrangian_<t_lag> *>(this)->setJacobianMatrix(std::forward<decltype(matrix)>(matrix));
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
            return static_cast<ElementLagrangian_<t_lag> *>(this)->setResidualVector(std::forward<decltype(vector)>(vector));
        }

        template<LagrangianConcept auto t_lag>
        auto const &
        getResidualVector()
        const
        {
            return static_cast<ElementLagrangian_<t_lag> const *>(this)->getResidualVector();
        }

    protected:

        Label const & label_;

    };
    
    /**
     * *********************************************************************************************************************************************************
     */

    template<Integer t_dim, Domain t_domain, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
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
    
    template<Integer t_dim, Domain t_domain, LagrangianConcept auto t_lag>
    struct DomainLagrangian : AbstractDomainLagrangian<t_dim, t_domain>
    {

    private:

        template<PotentialConcept auto t_potential>
        static constexpr
        Integer
        getPotentialIndex()
        {
            return LagTraits<t_lag>::template getIndex<t_potential>();
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

    template<Integer t_dim, Domain t_domain>
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
    
} // namespace lolita


#endif /* B53096ED_E28A_48CC_B15B_29B622D775A3 */
