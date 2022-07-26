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

    struct Load
    {

        auto static constexpr zero = [] (auto const &, auto const &) constexpr { return lolita::real(0); };

        Load()
        :
        loading_(std::make_shared<Loading>(zero))
        {}

        Load(
            Loading const & loading
        )
        :
        loading_(std::make_shared<Loading>(loading))
        {}

        Load(
            Loading && loading
        )
        :
        loading_(std::make_shared<Loading>(std::forward<Loading>(loading)))
        {}
        
        lolita::real
        getImposedValue(
            lolita2::Point const & point,
            lolita::real const & time
        )
        const
        {
            return loading_->operator ()(point, time);
        }

        std::shared_ptr<Loading> loading_;

    };

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
        scalar()
        {
            return Field(0);
        }

        static constexpr
        Field
        vector()
        {
            return Field(1);
        }
        
        constexpr
        Field(
            lolita::integer dim
        )
        :
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

        lolita::integer dim_;

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
    concept HybridDiscontinuousGalerkinConcept = detail::IsHybridDiscontinuousGalerkin<t_T>::value;

    template<typename... t_Mappings>
    struct Unknown
    {

        using Mappings = lolita::utility::Aggregate<t_Mappings...>;
        
        static constexpr
        lolita::integer
        getNumMappings()
        {
            return sizeof...(t_Mappings);
        }

        constexpr
        Unknown(
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
            Unknown const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Unknown const & other
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
        struct IsUnknown : std::false_type {};
        
        template<typename... t_T>
        struct IsUnknown<Unknown<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept UnknownConcept = detail::IsUnknown<t_T>::value;

    template<UnknownConcept... t_Unknowns>
    struct Behavior
    {

        using Unknowns = lolita::utility::Aggregate<t_Unknowns...>;

        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            return sizeof...(t_Unknowns);
        }

        constexpr
        Behavior(
            t_Unknowns... unknowns
        )
        :
        unknowns_(unknowns...)
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
        std::tuple_element_t<t_i, std::tuple<t_Unknowns...>> const &
        getUnknown()
        const
        {
            return unknowns_.template get<t_i>();
        }

        Unknowns unknowns_;

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

    template<UnknownConcept t_Unknown, BehaviorConcept t_Behavior, typename t_Discretization>
    struct FiniteElementMethod
    {

        constexpr
        FiniteElementMethod(
            t_Unknown unknown,
            t_Behavior behavior,
            t_Discretization discretization,
            Quadrature quadrature
        )
        :
        unknown_(unknown),
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

        t_Unknown unknown_;

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
