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
    struct ElementFormulation
    {

        struct IntegrationPoint
        {

        private:

            struct StrainOperator
            {

                // template<GeneralizedStrainConcept auto t_strain>
                // static
                // StrainOperator
                // make(
                //     MatrixConcept<Real> auto const & input
                // )
                // {
                //     return StrainOperator(t_strain.getTag(), input);
                // }

                template<GeneralizedStrainConcept auto t_strain>
                static
                StrainOperator
                make(
                    MatrixConcept<Real> auto && input
                )
                {
                    return StrainOperator(t_strain.getTag(), std::forward<std::decay_t<decltype(input)>>(input));
                }

            private:

                StrainOperator(
                    Integer tag,
                    MatrixConcept<Real> auto const & input
                )
                :
                tag_(tag),
                matrix_(input)
                {}

                // StrainOperator(
                //     Integer tag,
                //     MatrixConcept<Real> auto && input
                // )
                // :
                // tag_(tag),
                // matrix_(std::forward<std::decay_t<decltype(input)>>(input))
                // {}

                StrainOperator(
                    Integer tag,
                    MatrixConcept<Real> auto && input
                )
                :
                tag_(tag),
                matrix_(std::move(input))
                {}

            public:

                void
                setMatrix(
                    MatrixConcept<Real> auto const & input
                )
                {
                    matrix_ = input;
                }

                void
                setMatrix(
                    MatrixConcept<Real> auto && input
                )
                {
                    matrix_ = std::forward<std::decay_t<decltype(input)>>(input);
                }

                Matrix<Real> const &
                getMatrix()
                const
                {
                    return matrix_;
                }

                Matrix<Real> &
                getMatrix()
                {
                    return matrix_;
                }

                Integer
                getTag()
                const
                {
                    return tag_;
                }

            private:

                Integer tag_;

                Matrix<Real> matrix_;

            };

        public:

            template<FiniteElementMethodConcept auto t_finite_element_method>
            static constexpr
            Integer
            getGeneralizedStrainSize()
            {
                return FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            }

            template<BehaviorConcept auto t_behavior>
            static constexpr
            Integer
            getGeneralizedStrainSize()
            {
                return BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
            }

            static
            IntegrationPoint
            make(
                PointConcept auto && ref_pt,
                PointConcept auto && coordinates,
                Real weight,
                std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
            )
            {
                return IntegrationPoint(std::forward<std::decay_t<decltype(ref_pt)>>(ref_pt), std::forward<std::decay_t<decltype(coordinates)>>(coordinates), weight, behavior);
                // return IntegrationPoint(std::move(ref_pt), std::forward<std::decay_t<decltype(coordinates)>>(coordinates), weight, behavior);
            }

        private:

            IntegrationPoint(
                PointConcept auto const & ref_pt,
                PointConcept auto const & coordinates,
                Real weight,
                std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
            )
            :
            reference_coordinates_(ref_pt),
            coordinates_(coordinates),
            weight_(weight),
            behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behavior)))
            {}

            // IntegrationPoint(
            //     PointConcept auto && ref_pt,
            //     PointConcept auto && coordinates,
            //     Real weight,
            //     std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
            // )
            // :
            // reference_coordinates_(std::forward<std::decay_t<decltype(ref_pt)>>(ref_pt)),
            // coordinates_(std::forward<std::decay_t<decltype(coordinates)>>(coordinates)),
            // weight_(weight),
            // behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behavior)))
            // {}

            IntegrationPoint(
                PointConcept auto && ref_pt,
                PointConcept auto && coordinates,
                Real weight,
                std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
            )
            :
            reference_coordinates_(std::move(ref_pt)),
            coordinates_(std::move(coordinates)),
            weight_(weight),
            behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behavior)))
            {}

        public:
        
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
                std::basic_string_view<Character> material_property_label,
                std::function<Real(Point const &)> && function
            )
            {
                auto value = std::forward<std::function<Real(Point const &)>>(function)(coordinates_);
                mgis::behaviour::setMaterialProperty(behavior_data_->s0, std::string(material_property_label), value);
                mgis::behaviour::setMaterialProperty(behavior_data_->s1, std::string(material_property_label), value);
            }

            void
            setExternalVariable(
                std::basic_string_view<Character> material_property_label,
                std::function<Real(Point const &)> && function
            )
            {
                auto value = std::forward<std::function<Real(Point const &)>>(function)(coordinates_);
                mgis::behaviour::setExternalStateVariable(behavior_data_->s0, std::string(material_property_label), value);
                mgis::behaviour::setExternalStateVariable(behavior_data_->s1, std::string(material_property_label), value);
            }

            template<GeneralizedStrainConcept auto t_strain>
            Matrix<Real>
            getJacobian()
            const
            {
                auto size = GeneralizedStrainTraits<t_strain>::template getSize<t_domain>();
                auto jacobian_matrix = Matrix<Real>(size, size);
                jacobian_matrix.setZero();
                auto set_mapping_block = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr mapping = t_strain.template getMapping<t_i>();
                    auto constexpr mapping_size = GeneralizedStrainTraits<t_strain>::template getMappingSize<t_domain, mapping>();
                    auto constexpr offset = GeneralizedStrainTraits<t_strain>::template getMappingOffset<t_domain, mapping>();
                    auto rhs = algebra::View<Matrix<Real, mapping_size, mapping_size> const>(behavior_data_->K.data() + offset * offset);
                    auto lhs = jacobian_matrix.template block<mapping_size, mapping_size>(offset, offset);
                    lhs = rhs;
                    if constexpr (t_i < t_strain.getNumMappings() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_mapping_block(set_mapping_block);
                return jacobian_matrix;
            }

            Point &
            getCurrentCoordinates()
            {
                return coordinates_;
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

            template<GeneralizedStrainConcept auto t_strain>
            Matrix<Real> const &
            getStrainOperator()
            const
            {
                auto find_it = [&] (auto const & item)
                {
                    return item.getTag() == t_strain.getTag();
                };
                return std::find_if(strain_matrix_list_.begin(), strain_matrix_list_.end(), find_it)->getMatrix();
            }

            template<GeneralizedStrainConcept auto t_strain>
            void
            addStrainOperator(
                MatrixConcept<Real> auto && input
            )
            {
                strain_matrix_list_.push_back(StrainOperator::template make<t_strain>(std::forward<decltype(input)>(input)));
            }

            algebra::View<Point const> reference_coordinates_;
        
            Point coordinates_;
            
            Real weight_;
        
            // std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

            std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

            // std::map<std::basic_string<Character>, Matrix<Real>> ops_;

            std::vector<StrainOperator> strain_matrix_list_;

        };

        static
        ElementFormulation
        make(
            Integer tag,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            return ElementFormulation(tag, behavior);
        }

        static
        ElementFormulation
        make(
            Integer tag
        )
        {
            return ElementFormulation(tag);
        }

    private:

        ElementFormulation(
            Integer tag
        )
        :
        tag_(tag),
        integration_points_(),
        behavior_()
        {}

        explicit
        ElementFormulation(
            Integer tag,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        :
        tag_(tag),
        integration_points_(),
        behavior_(behavior)
        {}

    public:

        template<PointConcept t_ReferencePoint, PointConcept t_CurrentPoint>
        void
        addIntegrationPoint(
            t_ReferencePoint && ref_pt,
            t_CurrentPoint && coordinates,
            Real weight,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            auto ip = IntegrationPoint::make(std::forward<t_ReferencePoint>(ref_pt), std::forward<t_CurrentPoint>(coordinates), weight, behavior);
            integration_points_.push_back(std::move(ip));
        }

        std::vector<IntegrationPoint> const &
        getIntegrationPoints()
        const
        {
            return integration_points_;
        }

        std::vector<IntegrationPoint> &
        getIntegrationPoints()
        {
            return integration_points_;
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

        Integer
        getTag()
        const
        {
            return tag_;
        }

        void
        integrate(
            std::atomic<Boolean> & output_handler
        )
        {
            for (auto & ip : integration_points_)
            {
                ip.integrate(output_handler, behavior_);
            }            
        }

        void
        reserve()
        {
            for (auto & ip : integration_points_)
            {
                ip.reserve();
            }            
        }

        void
        recover()
        {
            for (auto & ip : integration_points_)
            {
                ip.recover();
            }            
        }

        Integer tag_;

        std::vector<IntegrationPoint> integration_points_;

        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        Matrix<Real> jacobian_matrix_;

        Vector<Real> residual_vector_;

    };
    
} // namespace lolita


#endif /* B53096ED_E28A_48CC_B15B_29B622D775A3 */
