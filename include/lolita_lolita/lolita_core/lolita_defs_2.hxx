#ifndef AB960F4D_38AB_4BE8_9B74_CB77F60CB560
#define AB960F4D_38AB_4BE8_9B74_CB77F60CB560

#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/BehaviourData.h>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"

namespace lolita2
{

    using Integer = lolita::integer;

    using Real = lolita::real;

    using Boolean = lolita::boolean;
    
    using Point = Eigen::Matrix<Real, 3, 1>;

    using Loading = std::function<lolita::real(lolita2::Point const &, lolita::real const &)>;

    enum struct DomainFrame
    {

        Cartesian,
        AxiSymmetric,

    };

    template<Frame t_frame, Integer t_dim>
    struct Domain
    {

        static constexpr
        Domain<t_frame, t_dim - 1>
        getSubDomain()
        const
        {
            return Domain<t_frame, t_dim - 1>();
        }

        static constexpr
        Boolean
        hasDim(
            Integer dim
        )
        const
        {
            return t_dim == dim;
        }
        
        static constexpr
        Boolean
        isAxiSymmetric()
        const
        {
            return t_frame == DomainFrame::AxiSymmetric;
        }
        
        static constexpr
        Boolean
        isCartesian()
        const
        {
            return t_frame == DomainFrame::Cartesian;
        }

        static constexpr
        Integer
        getDim()
        {
            return t_dim;
        }

        static constexpr
        DomainFrame
        getFrame()
        {
            return t_frame;
        }

    };

    namespace detail
    {

        template<typename t_T>
        struct DomainConceptTraits : std::false_type {};
        
        template<auto... t_args>
        struct DomainConceptTraits<Domain<t_args...>> : std::true_type {};

    }

    template<typename t_T>
    concept DomainConcept = detail::DomainConceptTraits<t_T>::value;

    enum struct QuadratureRule
    {

        Gauss,

    };

    template<QuadratureRule t_rule, Integer t_ord>
    struct Quadrature
    {
        
        static constexpr
        Boolean
        isGauss()
        {
            return t_rule == QuadratureRule::Gauss;
        }

        static constexpr
        Integer
        getOrd()
        {
            return t_ord;
        }

        static constexpr
        QuadratureRule
        getRule()
        {
            return t_rule;
        }

    };

    namespace detail
    {

        template<typename t_T>
        struct QuadratureConceptTraits : std::false_type {};
        
        template<auto... t_args>
        struct QuadratureConceptTraits<Quadrature<t_args...>> : std::true_type {};

    }

    template<typename t_T>
    concept QuadratureConcept = detail::QuadratureConceptTraits<t_T>::value;

    enum struct Polynomial
    {

        Monomial,

    };

    template<Polynomial t_polynomial, Integer t_ord>
    struct Basis
    {
        
        static constexpr
        Boolean
        isMonomial()
        {
            return t_rule == Polynomial::Monomial;
        }

        static constexpr
        Integer
        getOrd()
        {
            return t_ord;
        }

        static constexpr
        Polynomial
        getPolynomial()
        {
            return t_rule;
        }

    };

    namespace detail
    {

        template<typename t_T>
        struct BasisConceptTraits : std::false_type {};
        
        template<auto... t_args>
        struct BasisConceptTraits<Basis<t_args...>> : std::true_type {};

    }

    template<typename t_T>
    concept BasisConcept = detail::BasisConceptTraits<t_T>::value;

    enum struct Mapping
    {

        Gradient,
        Identity,
        Divergence,
        SmallStrain,
        LargeStrain,

    };

    template<Mapping t_mapping>
    struct DifferentialOperator
    {

        static constexpr
        Mapping
        getMapping()
        {
            return t_mapping;
        }
        

        static constexpr
        Boolean
        isGradient()
        {
            return t_mapping == Mapping::Gradient;
        }

        static constexpr
        Boolean
        isIdentity()
        {
            return t_mapping == Mapping::Identity;
        }

        static constexpr
        Boolean
        isDivergence()
        {
            return t_mapping == Mapping::Divergence;
        }

        static constexpr
        Boolean
        isSmallStrain()
        {
            return t_mapping == Mapping::SmallStrain;
        }

        static constexpr
        Boolean
        isLargeStrain()
        {
            return t_mapping == Mapping::LargeStrain;
        }

    };

    namespace detail
    {

        template<typename t_T>
        struct DifferentialOperatorConceptTraits : std::false_type {};
        
        template<auto... t_args>
        struct DifferentialOperatorConceptTraits<DifferentialOperator<t_args...>> : std::true_type {};

    }

    template<typename t_T>
    concept DifferentialOperatorConcept = detail::DifferentialOperatorConceptTraits<t_T>::value;

    template<Integer t_dim>
    struct Field
    {
        
        static constexpr
        Boolean
        isTensor(
            Integer dim
        )
        {
            return t_dim == dim;
        }

    };

    namespace detail
    {

        template<typename t_T>
        struct FieldConceptTraits : std::false_type {};
        
        template<auto... t_args>
        struct FieldConceptTraits<Field<t_args...>> : std::true_type {};

    }

    template<typename t_T>
    concept FieldConcept = detail::FieldConceptTraits<t_T>::value;

    template<FieldConcept auto t_field, DifferentialOperatorConcept auto... t_differential_operators>
    struct GeneralizedStrain
    {
        
        static constexpr
        Integer
        getNumMappings()
        {
            return sizeof...(t_differential_operators);
        }

        static constexpr
        FieldConcept auto
        getField()
        {
            return t_field;
        }

        template<Integer t_i>
        static constexpr
        DifferentialOperatorConcept auto
        getDifferentialOperator()
        {
            return std::get<t_i>(std::tuple<decltype(t_differential_operators)...>(t_differential_operators...));
        }

    };

    namespace detail
    {

        template<typename t_T>
        struct GeneralizedStrainConceptTraits : std::false_type {};
        
        template<auto... t_args>
        struct GeneralizedStrainConceptTraits<GeneralizedStrain<t_args...>> : std::true_type {};

    }

    template<typename t_T>
    concept GeneralizedStrainConcept = detail::GeneralizedStrainConceptTraits<t_T>::value;
    
    template<GeneralizedStrainConcept auto... t_generalized_strains>
    struct Behavior
    {
        
        static constexpr
        Integer
        getNumGeneralizedStrains()
        {
            return sizeof...(t_generalized_strains);
        }

    };
    
    template<BasisConcept auto t_cell_basis, BasisConcept auto t_face_basis, auto... t_args>
    struct HybridDiscontinuousGalerkin
    {
        
        static constexpr
        BasisConcept auto
        getCellBasis()
        {
            return t_cell_basis;
        }
        
        static constexpr
        BasisConcept auto
        getFaceBasis()
        {
            return t_face_basis;
        }

        template<Integer t_i>
        static constexpr
        auto
        getArg()
        {
            return return std::get<t_i>(std::tuple<decltype(t_args)...>(t_args...));
        }

    };

    // template<typename T, typename U>
    // struct Exp;

    // template<template<auto... t_args> typename T, template<auto> typename U>
    // struct Exp<T, U>
    // {
    //     using type = std::tuple<U<t_args>...>;
    // };

}

#endif /* AB960F4D_38AB_4BE8_9B74_CB77F60CB560 */
