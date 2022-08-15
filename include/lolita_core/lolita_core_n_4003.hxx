#ifndef AB07387F_F0B7_4F5A_875F_1CC5F0206148
#define AB07387F_F0B7_4F5A_875F_1CC5F0206148

#include "lolita_core/lolita.hxx"
#include "lolita_core/lolita_core_n_0000.hxx"
#include "lolita_core/lolita_core_n_1000.hxx"
#include "lolita_core/lolita_core_n_2000.hxx"
#include "lolita_core/lolita_core_n_3000.hxx"
#include "lolita_core/lolita_core_n_4001.hxx"
#include "lolita_core/lolita_core_n_4002.hxx"

namespace lolita
{

    struct QuadratureOperator
    {

        inline
        std::basic_string_view<Character>
        getLabel()
        const
        {
            return label_;
        }

        inline
        lolita::matrix::Matrix<Real>
        getOperator(
            Integer i
        )
        const
        {
            return operators_[i];
        }

        std::basic_string_view<Character> label_;

        std::vector<lolita::matrix::Matrix<Real>> operators_;

    };
    
    template<Domain t_domain>
    struct IntegrationPoint
    {

        IntegrationPoint(
            Point const & coordinates,
            Real const & weight,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behaviour
        )
        :
        coordinates_(coordinates),
        weight_(weight),
        behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behaviour)))
        {}

        std::unique_ptr<mgis::behaviour::BehaviourData> const &
        getMaterialPoint()
        const
        {
            return behavior_data_;
        }

        std::unique_ptr<mgis::behaviour::BehaviourData> &
        getMaterialPoint()
        {
            return behavior_data_;
        }
        
        template<BehaviorConcept auto t_behavior>
        lolita::matrix::Span<lolita::matrix::Vector<Real, BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<Real, size> const>(behavior_data_->s1.gradients.data());
        }
        
        template<auto t_finite_element_method>
        lolita::matrix::Span<RealVector<FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<Real, size> const>(behavior_data_->s1.gradients.data() + offset);
        }
        
        template<auto t_finite_element_method>
        lolita::matrix::Span<RealVector<FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()>>
        getGeneralizedStrain()
        {
            auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<Real, size>>(behavior_data_->s1.gradients.data() + offset);
        }
        
        Point coordinates_;
        
        Real weight_;

        std::vector<QuadratureOperator> quadrature_operators_;
        
        std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

    };

    struct QuadraturePoint
    {

        void
        integrate(
            Integer & res
        )
        {
            res = mgis::behaviour::integrate(* std::make_unique<mgis::behaviour::BehaviourDataView>(mgis::behaviour::make_view(* behavior_data_)), * behavior_);
        }
        
        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;
        
        std::shared_ptr<mgis::behaviour::BehaviourData> behavior_data_;

    };

    struct QuadraturePoint2 : QuadraturePoint
    {

        std::vector<lolita::matrix::Matrix<Real>> element_operators_;

    };

    template<Domain t_domain>
    struct ElementIntegrationPoints
    {

        ElementIntegrationPoints(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behaviour
        )
        :
        behavior_(behaviour)
        {}

        Integer
        getSize()
        const
        {
            return integration_points_.size();
        }

        std::basic_string_view<Character>
        getLabel()
        const
        {
            return behavior_->behaviour;
        }

        void
        integrate()
        {
            for (auto const & integration_point : integration_points_)
            {
                auto behaviour_view = mgis::behaviour::make_view(integration_point->behavior_data_);
                auto result = mgis::behaviour::integrate(behaviour_view, * behavior_);
            }
        }
        
        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        std::vector<IntegrationPoint<t_domain>> integration_points_;

    };
    
}

#endif /* AB07387F_F0B7_4F5A_875F_1CC5F0206148 */
