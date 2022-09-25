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

        template<PotentialConcept auto t_potential>
        Matrix<Real>
        getJacobian()
        const
        {
            auto size = PotentialTraits<t_potential>::template getSize<t_domain>();
            auto jacobian_matrix = Matrix<Real>(size, size);
            jacobian_matrix.setZero();
            auto set_mapping_block = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr mapping = t_potential.template getStrain<t_i>();
                auto constexpr mapping_size = PotentialTraits<t_potential>::template getMappingSize<t_domain, mapping>();
                auto constexpr offset = PotentialTraits<t_potential>::template getMappingOffset<t_domain, mapping>();
                auto rhs = algebra::View<Matrix<Real, mapping_size, mapping_size> const>(behavior_data_->K.data() + offset * offset);
                auto lhs = jacobian_matrix.template block<mapping_size, mapping_size>(offset, offset);
                lhs = rhs;
                if constexpr (t_i < t_potential.getNumMappings() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_mapping_block(set_mapping_block);
            return jacobian_matrix;
        }

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

        template<Mapping t_strain>
        void
        addStrainOperator(
            MatrixConcept<Real> auto && input
        )
        {
            if (strain_matrix_list_ == nullptr)
            {
                strain_matrix_list_ = std::make_unique<std::vector<ElementaryOperator<Mapping, Matrix<Real>>>>();
            }
            for (auto const & m : * strain_matrix_list_)
            {
                if (m.getTag() == t_strain)
                {
                    return;
                }
            }
            strain_matrix_list_->push_back(ElementaryOperator<Mapping, Matrix<Real>>(t_strain, std::forward<decltype(input)>(input)));
        }

        template<Mapping t_strain>
        Matrix<Real> const &
        getStrainOperator()
        const
        {
            for (auto const & m : * strain_matrix_list_)
            {
                if (m.getTag() == t_strain)
                {
                    return m.getOperator();
                }
            }
            throw std::runtime_error("No usch thing");
        }

        template<Mapping t_strain>
        Boolean
        hasStrainOperator()
        const
        {
            for (auto const & m : * strain_matrix_list_)
            {
                if (m.getTag() == t_strain)
                {
                    return true;
                }
            }
            return false;
        }

        template<Field t_field>
        void
        addDiscreteField()
        {
            if (ptr_data_ == nullptr)
            {
                ptr_data_ = std::make_unique<std::vector<ElementDiscreteField<t_element, t_domain>>>();
            }
            for (auto const & item : * ptr_data_)
            {
                if (item.getLabel() == t_field.getLabel())
                {
                    return;
                }
            }
            ptr_data_->push_back(ElementDiscreteField<t_element, t_domain>(t_field));
        }

        template<Field t_field>
        ElementDiscreteField<t_element, t_domain> const &
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

        template<Field t_field>
        ElementDiscreteField<t_element, t_domain> &
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

        template<Field t_field, Basis t_basis>
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

        std::unique_ptr<std::vector<ElementaryOperator<Mapping, Matrix<Real>>>> strain_matrix_list_;

        std::unique_ptr<std::vector<ElementDiscreteField<t_element, t_domain>>> ptr_data_;

    };
    
    template<Element t_element, Domain t_domain>
    struct ElementFormulation
    {

    private:

        using t_IntegrationPoint = IntegrationPoint<t_element, t_domain>;

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

        utility::Label const &
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

        utility::Label const & label_;

        std::unique_ptr<std::vector<t_IntegrationPoint>> integration_points_;

        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        Matrix<Real> jacobian_matrix_;

        Vector<Real> residual_vector_;

    };
    
} // namespace lolita


#endif /* B53096ED_E28A_48CC_B15B_29B622D775A3 */
