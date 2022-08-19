#ifndef B9B48BEA_09F5_41DB_84A1_45A6E708901C
#define B9B48BEA_09F5_41DB_84A1_45A6E708901C

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_3000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4001.hxx"

namespace lolita
{

    struct ElementDegreeOfFreedom
    {

        template<Element t_element, Domain t_domain, Field t_field, Basis t_basis>
        static constexpr
        Integer
        getSize()
        {
            auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            return t_field_size * t_basis_size;
        }

        template<Element t_element, Domain t_domain, Field t_field, Basis t_basis>
        static inline
        ElementDegreeOfFreedom
        make(
            std::shared_ptr<DegreeOfFreedom> const & degree_of_freedom
        )
        {
            auto element_degree_of_freedom = ElementDegreeOfFreedom{degree_of_freedom, Natural(degree_of_freedom->coefficients_.size())};
            degree_of_freedom->coefficients_.resize(degree_of_freedom->coefficients_.size() + getSize<t_element, t_domain, t_field, t_basis>());
            return element_degree_of_freedom;
        }
        
        inline
        Boolean
        operator==(
            ElementDegreeOfFreedom const & other
        )
        const = default;
        
        inline
        Boolean
        operator!=(
            ElementDegreeOfFreedom const & other
        )
        const = default;

        inline
        std::shared_ptr<DegreeOfFreedom> const &
        getDegreeOfFreedom()
        const
        {
            return degree_of_freedom_;
        }

        inline
        std::shared_ptr<DegreeOfFreedom> &
        getDegreeOfFreedom()
        {
            return degree_of_freedom_;
        }
        
        inline
        std::basic_string_view<Character>
        getLabel()
        const
        {
            return degree_of_freedom_->getLabel();
        }

        inline
        Natural
        getTag()
        const
        {
            return tag_;
        }

        template<Element t_element, Domain t_domain, Field t_field, Basis t_basis>
        lolita::algebra::Span<lolita::algebra::Vector<Real, getSize<t_element, t_domain, t_field, t_basis>()> const>
        getCoefficients()
        const
        {
            auto const & data = degree_of_freedom_->coefficients_.data() + tag_;
            return lolita::algebra::Span<lolita::algebra::Vector<Real, getSize<t_element, t_domain, t_field, t_basis>()> const>(data);
        }

        template<Element t_element, Domain t_domain, Field t_field, Basis t_basis>
        lolita::algebra::Span<lolita::algebra::Vector<Real, getSize<t_element, t_domain, t_field, t_basis>()>>
        getCoefficients()
        {
            auto const & data = degree_of_freedom_->coefficients_.data() + tag_;
            return lolita::algebra::Span<lolita::algebra::Vector<Real, getSize<t_element, t_domain, t_field, t_basis>()>>(data);
        }

        template<Element t_element, Domain t_domain, Field t_field, Basis t_basis>
        lolita::algebra::Span<lolita::algebra::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>
        getCoefficients(
            Integer row,
            Integer col
        )
        {
            auto constexpr t_field_shape = FieldTraits<t_field>::template shape<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            auto const & data = degree_of_freedom_->coefficients_.data() + tag_ + (t_field_shape.cols() * row  + col) * t_basis_size;
            return lolita::algebra::Span<lolita::algebra::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
        }

        template<Element t_element, Domain t_domain, Field t_field, Basis t_basis>
        lolita::algebra::Span<lolita::algebra::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()> const>
        getCoefficients(
            Integer row,
            Integer col
        )
        const
        {
            auto constexpr t_field_shape = FieldTraits<t_field>::template shape<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            auto const & data = degree_of_freedom_->coefficients_.data() + tag_ + (t_field_shape.cols() * row  + col) * t_basis_size;
            return lolita::algebra::Span<lolita::algebra::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
        }

        std::shared_ptr<DegreeOfFreedom> degree_of_freedom_;

        Natural tag_;

    };

    struct ElementLoad
    {

        static inline
        ElementLoad
        make(
            std::shared_ptr<Load> const & load
        )
        {
            return ElementLoad{load};
        }
        
        inline
        Boolean
        operator==(
            ElementLoad const & other
        )
        const = default;
        
        inline
        Boolean
        operator!=(
            ElementLoad const & other
        )
        const = default;

        std::shared_ptr<Load> const &
        getLoad()
        const
        {
            return load_;
        }

        std::shared_ptr<Load> &
        getLoad()
        {
            return load_;
        }
        
        inline
        Real
        getImposedValue(
            Point const & point,
            Real const & time
        )
        const
        {
            return load_->loading_(point, time);
        }

        std::shared_ptr<Load> load_;

    };

    struct QuadratureElement
    {

        struct IntegrationPoint
        {

            // static inline
            // IntegrationPoint
            // make(
            //     Point const & coordinates,
            //     Real const & weight,
            //     std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
            // )
            // {
            //     auto behavior_data = mgis::behaviour::BehaviourData(* behavior);
            //     return IntegrationPoint{
            //         .coordinates_ = coordinates,
            //         .weight_ = weight,
            //         .behavior_ = behavior,
            //         .behavior_data_ = behavior_data,
            //         .behavior_data_view_ = mgis::behaviour::make_view(behavior_data)
            //     };
            // }

            IntegrationPoint(
                Point const & coordinates,
                Real const & weight,
                std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
            )
            :
            coordinates_(coordinates),
            weight_(weight),
            behavior_(behavior),
            behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behavior))),
            behavior_data_view_(std::make_unique<mgis::behaviour::BehaviourDataView>(mgis::behaviour::make_view(* behavior_data_)))
            {}
        
            inline
            Boolean
            operator==(
                IntegrationPoint const & other
            )
            const = default;
            
            inline
            Boolean
            operator!=(
                IntegrationPoint const & other
            )
            const = default;
            
            template<Domain t_domain, BehaviorConcept auto t_behavior>
            lolita::algebra::Span<lolita::algebra::Vector<Real, BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>()> const>
            getGeneralizedStrain()
            const
            {
                auto constexpr size = BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
                return lolita::algebra::Span<lolita::algebra::Vector<Real, size> const>(behavior_data_->s1.gradients.data());
            }
            
            template<Domain t_domain, auto t_finite_element_method>
            lolita::algebra::Span<RealVector<FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()> const>
            getGeneralizedStrain()
            const
            {
                auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
                auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
                return lolita::algebra::Span<lolita::algebra::Vector<Real, size> const>(behavior_data_->s1.gradients.data() + offset);
            }
            
            template<Domain t_domain, auto t_finite_element_method>
            lolita::algebra::Span<RealVector<FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()>>
            getGeneralizedStrain()
            {
                auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
                auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
                return lolita::algebra::Span<lolita::algebra::Vector<Real, size>>(behavior_data_->s1.gradients.data() + offset);
            }
        
            Point coordinates_;
            
            Real weight_;
        
            std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

            std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

            std::unique_ptr<mgis::behaviour::BehaviourDataView> behavior_data_view_;

            std::map<std::basic_string<Character>, RealMatrix<>> ops_;

        };

        static inline
        QuadratureElement
        make(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            return QuadratureElement{.ips_ = {}, .behavior_ = behavior};
        }
        
        inline
        Boolean
        operator==(
            QuadratureElement const & other
        )
        const = default;
        
        inline
        Boolean
        operator!=(
            QuadratureElement const & other
        )
        const = default;

        // QuadratureElement()
        // {}

        // explicit
        // QuadratureElement(
        //     std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        // )
        // :
        // behavior_(behavior)
        // {}

        // explicit
        // QuadratureElement(
        //     std::shared_ptr<mgis::behaviour::Behaviour> && behavior
        // )
        // :
        // behavior_(std::forward<std::shared_ptr<mgis::behaviour::Behaviour>>(behavior))
        // {}

        // void inline
        // integrate(
        //     Integer & res
        // )
        // {
        //     res = mgis::behaviour::integrate(behavior_data_view_, * behavior_);
        // }

        std::vector<IntegrationPoint> ips_;
        
        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        // QuadratureElement(
        //     Point const & coordinates,
        //     Real const & weight,
        //     std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        // )
        // :
        // coordinates_(coordinates),
        // weight_(weight),
        // behavior_(behavior),
        // behavior_data_(mgis::behaviour::BehaviourData(* behavior)),
        // behavior_data_view_(mgis::behaviour::make_view(behavior_data_))
        // {}

        
        
        // Point coordinates_;
        
        // Real weight_;
        
        // std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        // mgis::behaviour::BehaviourData behavior_data_;

        // mgis::behaviour::BehaviourDataView behavior_data_view_;

        // std::vector<RealMatrix<>> operators_;

    };
    
} // namespace lolita

#endif /* B9B48BEA_09F5_41DB_84A1_45A6E708901C */
