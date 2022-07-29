#ifndef AD31B868_1D56_48DB_A2CC_F898A668702D
#define AD31B868_1D56_48DB_A2CC_F898A668702D

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

    // struct MeshFileFormat
    // {

    //     enum Format
    //     {

    //         Gmsh,

    //     };

    //     constexpr
    //     MeshFileFormat(
    //         Format format
    //     )
    //     :
    //     format_(format)
    //     {}

    //     constexpr
    //     lolita::boolean
    //     operator==(
    //         MeshFileFormat const & other
    //     )
    //     const = default;

    //     constexpr
    //     lolita::boolean
    //     operator!=(
    //         MeshFileFormat const & other
    //     )
    //     const = default;

    //     constexpr
    //     lolita::boolean
    //     isGmsh()
    //     const
    //     {
    //         return format_ == Format::Gmsh;
    //     }

    //     Format format_;

    // };
    
    using Point = lolita::matrix::Vector<lolita::real, 3>;

    using Loading = std::function<lolita::real(lolita2::Point const &, lolita::real const &)>;

    struct Domain
    {

        enum Frame
        {

            Cartesian,
            AxiSymmetric

        };

        constexpr
        Domain(
            lolita::integer dim,
            Frame frame
        )
        :
        dim_(dim),
        frame_(frame)
        {}

        constexpr
        lolita::boolean
        operator==(
            Domain const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Domain const & other
        )
        const = default;

        constexpr
        Domain
        getSubDomain()
        const
        {
            return Domain(dim_ - 1, frame_);
        }

        constexpr
        lolita::boolean
        hasDim(
            lolita::integer dim
        )
        const
        {
            return dim_ == dim;
        }

        static constexpr
        Domain
        axiSymmetric(
            lolita::integer dim
        )
        {
            return Domain(dim, Frame::AxiSymmetric);
        }
        
        constexpr
        lolita::boolean
        isAxiSymmetric()
        const
        {
            return frame_ == Frame::AxiSymmetric;
        }

        static constexpr
        Domain
        cartesian(
            lolita::integer dim
        )
        {
            return Domain(dim, Frame::Cartesian);
        }
        
        constexpr
        lolita::boolean
        isCartesian()
        const
        {
            return frame_ == Frame::Cartesian;
        }

        lolita::integer dim_;

        Frame frame_;

    };

    struct Quadrature
    {
        
        enum Rule
        {

            Gauss,

        };
        
        constexpr
        Quadrature(
            Quadrature::Rule rule,
            lolita::integer ord
        )
        :
        rule_(rule),
        ord_(ord)
        {}

        constexpr
        lolita::boolean
        operator==(
            Quadrature const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Quadrature const & other
        )
        const = default;

        static constexpr
        Quadrature
        gauss(
            lolita::integer ord
        )
        {
            return Quadrature(Rule::Gauss, ord);
        }
        
        constexpr
        lolita::boolean
        isGauss()
        const
        {
            return rule_ == Rule::Gauss;
        }
        
        Quadrature::Rule rule_;
        
        lolita::integer ord_;

    };

    struct Basis
    {

        enum Polynomial
        {

            Lagrange,
            Monomial,

        };
        
        constexpr
        Basis(
            Polynomial polynomial,
            lolita::integer ord
        )
        :
        polynomial_(polynomial),
        ord_(ord)
        {}

        constexpr
        lolita::boolean
        operator==(
            Basis const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Basis const & other
        )
        const = default;

        static constexpr
        Basis
        monomial(
            lolita::integer ord
        )
        {
            return Basis(Polynomial::Monomial, ord);
        }
        
        constexpr
        lolita::boolean
        isMonomial()
        const
        {
            return polynomial_ == Polynomial::Monomial;
        }

        constexpr
        lolita::integer
        getOrder()
        const
        {
            return ord_;
        }
        
        Polynomial polynomial_;
        
        lolita::integer ord_;

    };

    struct Mapping
    {

        enum Type
        {

            Gradient,
            Identity,
            Divergence,
            SmallStrain,
            LargeStrain,

        };

        constexpr
        Mapping(
            Type type
        )
        :
        type_(type)
        {}

        constexpr
        lolita::boolean
        operator==(
            Mapping const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Mapping const & other
        )
        const = default;

        static constexpr
        Mapping
        gradient()
        {
            return Mapping(Type::Gradient);
        }

        constexpr
        lolita::boolean
        isGradient()
        const
        {
            return type_ == Type::Gradient;
        }

        static constexpr
        Mapping
        identity()
        {
            return Mapping(Type::Identity);
        }

        constexpr
        lolita::boolean
        isIdentity()
        const
        {
            return type_ == Type::Identity;
        }

        static constexpr
        Mapping
        divergence()
        {
            return Mapping(Type::Divergence);
        }

        constexpr
        lolita::boolean
        isDivergence()
        const
        {
            return type_ == Type::Divergence;
        }

        static constexpr
        Mapping
        smallStrain()
        {
            return Mapping(Type::SmallStrain);
        }

        constexpr
        lolita::boolean
        isSmallStrain()
        const
        {
            return type_ == Type::SmallStrain;
        }

        static constexpr
        Mapping
        largeStrain()
        {
            return Mapping(Type::LargeStrain);
        }

        constexpr
        lolita::boolean
        isLargeStrain()
        const
        {
            return type_ == Type::LargeStrain;
        }

        Type type_;

    };

    struct Field
    {

        static constexpr
        Field
        scalar(
            std::basic_string_view<lolita::character> label
        )
        {
            return Field(label, 0);
        }

        static constexpr
        Field
        vector(
            std::basic_string_view<lolita::character> label
        )
        {
            return Field(label, 1);
        }
        
        constexpr
        Field(
            lolita::utility::Label const & label,
            lolita::integer dim
        )
        :
        label_(label),
        dim_(dim)
        {}
        
        constexpr
        Field(
            std::basic_string_view<lolita::character> label,
            lolita::integer dim
        )
        :
        label_(label),
        dim_(dim)
        {}

        constexpr
        lolita::boolean
        operator==(
            Field const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Field const & other
        )
        const = default;

        constexpr
        lolita::boolean
        isTensor(
            lolita::integer dim
        )
        const
        {
            return dim_ == dim;
        }

        lolita::utility::Label label_;

        lolita::integer dim_;

    };

    struct Unknown
    {

        enum Type
        {

            Cell,
            Face,
            Edge,
            Node,

        };
        
        constexpr
        Unknown(
            Field field,
            Basis basis
        )
        :
        is_active_(true),
        field_(field),
        basis_(basis)
        {}

        lolita::boolean is_active_;

        Field field_;

        Basis basis_;

    };

    struct HybridDiscontinuousGalerkin
    {

        enum Stabilization
        {

            Hdg,
            Hho,

        };

        constexpr
        HybridDiscontinuousGalerkin(
            Basis cell_basis,
            Basis face_basis,
            Stabilization stabilization
        )
        :
        cell_basis_(cell_basis),
        face_basis_(face_basis),
        grad_basis_(face_basis),
        stabilization_(stabilization)
        {}

        constexpr
        lolita::boolean
        operator==(
            HybridDiscontinuousGalerkin const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            HybridDiscontinuousGalerkin const & other
        )
        const = default;

        constexpr
        lolita::boolean
        isHybridDiscontinuousGalerkin()
        const
        {
            return true;
        }

        constexpr
        lolita::boolean
        isHdg()
        const
        {
            return stabilization_ == Stabilization::Hdg;
        }

        constexpr
        lolita::boolean
        isHho()
        const
        {
            return stabilization_ == Stabilization::Hho;
        }

        Basis cell_basis_;

        Basis face_basis_;

        Basis grad_basis_;

        Stabilization stabilization_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsHybridDiscontinuousGalerkin : std::false_type {};
        
        template<>
        struct IsHybridDiscontinuousGalerkin<HybridDiscontinuousGalerkin> : std::true_type {};

    }

    template<typename t_T>
    concept HybridDiscontinuousGalerkinConcept = detail::IsHybridDiscontinuousGalerkin<std::decay_t<t_T>>::value;

    template<typename... t_Mappings>
    struct GeneralizedStrain
    {

        using Mappings = lolita::utility::Aggregate<t_Mappings...>;
        
        static constexpr
        lolita::integer
        getNumMappings()
        {
            return sizeof...(t_Mappings);
        }

        constexpr
        GeneralizedStrain(
            Field field,
            t_Mappings... mappings
        )
        :
        field_(field),
        mappings_(mappings...)
        {}

        constexpr
        lolita::boolean
        operator==(
            GeneralizedStrain const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            GeneralizedStrain const & other
        )
        const = default;

        constexpr
        Field
        getField()
        const
        {
            return field_;
        }

        template<lolita::integer t_i>
        constexpr
        Mapping
        getMapping()
        const
        {
            return mappings_.template get<t_i>();
        }

        Field field_;

        Mappings mappings_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsGeneralizedStrain : std::false_type {};
        
        template<typename... t_T>
        struct IsGeneralizedStrain<GeneralizedStrain<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept GeneralizedStrainConcept = detail::IsGeneralizedStrain<t_T>::value;

    template<GeneralizedStrainConcept... t_GeneralizedStrains>
    struct Behavior
    {

        using GeneralizedStrains = lolita::utility::Aggregate<t_GeneralizedStrains...>;

        static constexpr
        lolita::integer
        getNumGeneralizedStrains()
        {
            return sizeof...(t_GeneralizedStrains);
        }

        constexpr
        Behavior(
            t_GeneralizedStrains... generalized_strains
        )
        :
        generalized_strains_(generalized_strains...)
        {}

        constexpr
        lolita::boolean
        operator==(
            Behavior const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Behavior const & other
        )
        const = default;

        template<lolita::integer t_i>
        constexpr
        std::tuple_element_t<t_i, std::tuple<t_GeneralizedStrains...>> const &
        getGeneralizedStrain()
        const
        {
            return generalized_strains_.template get<t_i>();
        }

        GeneralizedStrains generalized_strains_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsBehavior : std::false_type {};
        
        template<typename... t_T>
        struct IsBehavior<Behavior<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept BehaviorConcept = detail::IsBehavior<t_T>::value;

    template<GeneralizedStrainConcept t_GeneralizedStrain, BehaviorConcept t_Behavior, typename t_Discretization>
    struct FiniteElementMethod
    {

        constexpr
        FiniteElementMethod(
            t_GeneralizedStrain generalized_strain,
            t_Behavior behavior,
            t_Discretization discretization,
            Quadrature quadrature
        )
        :
        generalized_strain_(generalized_strain),
        behavior_(behavior),
        discretization_(discretization),
        quadrature_(quadrature)
        {}

        constexpr
        lolita::boolean
        operator==(
            FiniteElementMethod const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            FiniteElementMethod const & other
        )
        const = default;

        constexpr
        lolita::boolean
        isHdg()
        const
        {
            return HybridDiscontinuousGalerkinConcept<t_Discretization>;
        }

        constexpr
        Field
        getField()
        const
        {
            return generalized_strain_.field_;
        }

        template<lolita::integer t_i>
        constexpr
        Mapping
        getMapping()
        const
        {
            return generalized_strain_.mappings_.template get<t_i>();
        }

        constexpr
        t_GeneralizedStrain const &
        getGeneralizedStrain()
        const
        {
            return generalized_strain_;
        }

        constexpr
        t_Behavior const &
        getBehavior()
        const
        {
            return behavior_;
        }

        constexpr
        t_Discretization const &
        getDiscretization()
        const
        {
            return discretization_;
        }

        constexpr
        Quadrature const &
        getQuadrature()
        const
        {
            return quadrature_;
        }

        // Unknown cell_unknown_;

        // Unknown face_unknown_;

        // Unknown edge_unknown_;

        // Unknown node_unknown_;

        t_GeneralizedStrain generalized_strain_;

        t_Behavior behavior_;

        t_Discretization discretization_;

        Quadrature quadrature_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsFiniteElementMethod : std::false_type {};
        
        template<typename... t_T>
        struct IsFiniteElementMethod<FiniteElementMethod<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept FiniteElementMethodConcept = detail::IsFiniteElementMethod<t_T>::value;

}


#endif /* AD31B868_1D56_48DB_A2CC_F898A668702D */
