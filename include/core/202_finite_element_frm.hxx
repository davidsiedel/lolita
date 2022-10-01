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

        // template<PotentialConcept auto t_potential>
        // DenseMatrix<Real>
        // getJacobian()
        // const
        // {
        //     // auto size = PotentialTraits<t_potential>::template getSize<t_domain>();
        //     // auto jacobian_matrix = DenseMatrix<Real>(size, size);
        //     // jacobian_matrix.setZero();
        //     // auto set_mapping_block = [&] <Integer t_i = 0> (
        //     //     auto & self
        //     // )
        //     // constexpr mutable
        //     // {
        //     //     auto constexpr mapping = PotentialTraits<t_potential>::template getStrain<t_i>();
        //     //     auto constexpr mapping_size = MappingTraits<mapping>::template getSize<t_domain>();
        //     //     auto constexpr offset = PotentialTraits<t_potential>::template getOffset<t_domain, mapping>();
        //     //     auto set_mapping_block = [&] <Integer t_i = 0> (
        //     //         auto & self
        //     //     )
        //     //     constexpr mutable
        //     //     {

        //     //     };
        //     //     auto rhs = algebra::View<DenseMatrix<Real, mapping_size, mapping_size> const>(behavior_data_->K.data() + offset * offset);
        //     //     auto lhs = jacobian_matrix.template block<mapping_size, mapping_size>(offset, offset);
        //     //     lhs = rhs;
        //     //     if constexpr (t_i < t_potential.getNumMappings() - 1)
        //     //     {
        //     //         self.template operator ()<t_i + 1>(self);
        //     //     }
        //     // };
        //     // set_mapping_block(set_mapping_block);
        //     // return jacobian_matrix;
        // }

        // template<PotentialConcept auto t_potential>
        // DenseMatrix<Real, PotentialTraits<t_potential>::template getSize<t_domain>(), PotentialTraits<t_potential>::template getSize<t_domain>()>
        // getJacobian()
        // const
        // {
        //     auto constexpr size = PotentialTraits<t_potential>::template getSize<t_domain>();
        //     auto row_offset = 0;
        //     auto col_offset = 0;
        //     auto offset = 0;
        //     auto jacobian_matrix = DenseMatrix<Real, size, size>();
        //     jacobian_matrix.setZero();
        //     auto set_mapping_rows = [&] <Integer t_row = 0> (
        //         auto & t_set_mapping_rows
        //     )
        //     constexpr mutable
        //     {
        //         auto constexpr row_mapping_size = MappingTraits<t_potential.template getStrain<t_row>()>::template getSize<t_domain>();
        //         auto set_jacobian_cols = [&] <Integer t_col = 0> (
        //             auto & t_set_jacobian_cols
        //         )
        //         constexpr mutable
        //         {
        //             auto constexpr col_mapping_size = MappingTraits<t_potential.template getStrain<t_col>()>::template getSize<t_domain>();
        //             auto rhs = algebra::View<DenseMatrix<Real, row_mapping_size, col_mapping_size> const>(behavior_data_->K.data() + offset);
        //             auto lhs = jacobian_matrix.template block<row_mapping_size, col_mapping_size>(row_offset, col_offset);
        //             lhs = rhs;
        //             offset += row_mapping_size * col_mapping_size;
        //             col_offset += col_mapping_size;
        //             if constexpr (t_col < t_potential.getNumMappings() - 1)
        //             {
        //                 t_set_jacobian_cols.template operator ()<t_col + 1>(t_set_jacobian_cols);
        //             }
        //         };
        //         set_jacobian_cols(set_jacobian_cols);
        //         row_offset += row_mapping_size;
        //         if constexpr (t_row < t_potential.getNumMappings() - 1)
        //         {
        //             t_set_mapping_rows.template operator ()<t_row + 1>(t_set_mapping_rows);
        //         }
        //     };
        //     set_mapping_rows(set_mapping_rows);
        //     return jacobian_matrix;
        // }

        // Point &
        // getCurrentCoordinates()
        // {
        //     return coordinates_;
        // }
        
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

        template<FieldConcept auto t_field>
        void
        addDiscreteField()
        {
            if (ptr_data_ == nullptr)
            {
                ptr_data_ = std::make_unique<std::vector<ElementDiscreteField<t_domain>>>();
            }
            for (auto const & item : * ptr_data_)
            {
                if (item.getLabel() == t_field.getLabel())
                {
                    return;
                }
            }
            ptr_data_->push_back(ElementDiscreteField<t_domain>(t_field));
        }

        template<FieldConcept auto t_field>
        ElementDiscreteField<t_domain> const &
        getDiscreteField()
        const
        {
            if (ptr_data_ == nullptr)
            {
                throw std::runtime_error("Empty");
            }
            else
            {
                for (auto const & item : * ptr_data_)
                {
                    if (item.getLabel() == t_field.getLabel())
                    {
                        return item;
                    }
                }
                throw std::runtime_error("No such field data");
            }
        }

        template<FieldConcept auto t_field>
        ElementDiscreteField<t_domain> &
        getDiscreteField()
        {
            if (ptr_data_ == nullptr)
            {
                throw std::runtime_error("Empty");
            }
            else
            {
                for (auto & item : * ptr_data_)
                {
                    if (item.getLabel() == t_field.getLabel())
                    {
                        return item;
                    }
                }
                throw std::runtime_error("No such field data");
            }
        }

        template<FieldConcept auto t_field, Basis t_basis>
        void
        addDiscreteFieldDegreeOfFreedom(
            auto const &... args
        )
        {
            this->template getDiscreteField<t_field>().template addDegreeOfFreedom<t_field, t_basis>(args...);
        }

        algebra::View<Point const> reference_coordinates_;
    
        Point coordinates_;
        
        Real weight_;

        std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

        std::unique_ptr<std::vector<ElementaryOperator<Label, DenseMatrix<Real>>>> strain_matrix_list_;

        std::unique_ptr<std::vector<ElementDiscreteField<t_domain>>> ptr_data_;

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

        Label const & label_;

        std::unique_ptr<std::vector<t_IntegrationPoint>> integration_points_;

        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        DenseMatrix<Real> jacobian_matrix_;

        DenseVector<Real> residual_vector_;

    };
    
} // namespace lolita


#endif /* B53096ED_E28A_48CC_B15B_29B622D775A3 */
