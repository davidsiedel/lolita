#ifndef A4BCD9B5_985A_4D19_B3E9_7C559F45A353
#define A4BCD9B5_985A_4D19_B3E9_7C559F45A353

#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/BehaviourData.h>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

#include "lolita_lolita/utility.hxx"
#include "lolita_lolita/config.hxx"
#include "lolita_lolita/numerics.hxx"
#include "lolita_lolita/algebra.hxx"

namespace lolita
{
    
    using Point = RealVector<3>;

    using Loading = std::function<Real(Point const &, Real const &)>;

    struct Domain
    {

        enum Frame
        {

            Cartesian,
            AxiSymmetric

        };

        constexpr
        Domain(
            Integer dim,
            Frame frame
        )
        :
        dim_(dim),
        frame_(frame)
        {}

        constexpr
        Boolean
        operator==(
            Domain const & other
        )
        const = default;

        constexpr
        Boolean
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
        Boolean
        hasDim(
            Integer dim
        )
        const
        {
            return dim_ == dim;
        }

        static constexpr
        Domain
        axiSymmetric(
            Integer dim
        )
        {
            return Domain(dim, Frame::AxiSymmetric);
        }
        
        constexpr
        Boolean
        isAxiSymmetric()
        const
        {
            return frame_ == Frame::AxiSymmetric;
        }

        static constexpr
        Domain
        cartesian(
            Integer dim
        )
        {
            return Domain(dim, Frame::Cartesian);
        }
        
        constexpr
        Boolean
        isCartesian()
        const
        {
            return frame_ == Frame::Cartesian;
        }

        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }

        constexpr
        Frame
        getFrame()
        const
        {
            return frame_;
        }

        Integer dim_;

        Frame frame_;

    };

    struct ElementType
    {

        enum Type
        {

            // Points,
            // Curves,
            // Facets,
            // Solids,
            Cells,
            Faces,
            Edges,
            Nodes,

        };

        static constexpr
        ElementType
        cells(
            Domain domain
        )
        {
            return ElementType(domain, Type::Cells);
        }

        static constexpr
        ElementType
        faces(
            Domain domain
        )
        {
            return ElementType(domain, Type::Faces);
        }

        constexpr
        ElementType(
            Domain domain,
            Type type
        )
        :
        domain_(domain),
        type_(type)
        {}
        
        constexpr
        Integer
        getDim()
        const
        {
            switch (type_)
            {
                case Type::Cells: return domain_.getDim() - 0;
                case Type::Faces: return domain_.getDim() - 1;
                case Type::Edges: return 1;
                case Type::Nodes: return 0;
                default : return -1;
            }
        }

        Domain domain_;

        Type type_;

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
            Integer ord
        )
        :
        rule_(rule),
        ord_(ord)
        {}

        constexpr
        Boolean
        operator==(
            Quadrature const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Quadrature const & other
        )
        const = default;

        static constexpr
        Quadrature
        gauss(
            Integer ord
        )
        {
            return Quadrature(Rule::Gauss, ord);
        }
        
        constexpr
        Boolean
        isGauss()
        const
        {
            return rule_ == Rule::Gauss;
        }
        
        Quadrature::Rule rule_;
        
        Integer ord_;

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
            Integer ord
        )
        :
        polynomial_(polynomial),
        ord_(ord)
        {}

        constexpr
        Boolean
        operator==(
            Basis const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Basis const & other
        )
        const = default;

        static constexpr
        Basis
        monomial(
            Integer ord
        )
        {
            return Basis(Polynomial::Monomial, ord);
        }
        
        constexpr
        Boolean
        isMonomial()
        const
        {
            return polynomial_ == Polynomial::Monomial;
        }

        constexpr
        Integer
        getOrder()
        const
        {
            return ord_;
        }
        
        Polynomial polynomial_;
        
        Integer ord_;

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
        Boolean
        operator==(
            Mapping const & other
        )
        const = default;

        constexpr
        Boolean
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
        Boolean
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
        Boolean
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
        Boolean
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
        Boolean
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
        Boolean
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
            Integer dim
        )
        :
        dim_(dim)
        {}

        constexpr
        Boolean
        operator==(
            Field const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Field const & other
        )
        const = default;

        constexpr
        Boolean
        isTensor(
            Integer dim
        )
        const
        {
            return dim_ == dim;
        }

        Integer dim_;

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
        Boolean
        operator==(
            HybridDiscontinuousGalerkin const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            HybridDiscontinuousGalerkin const & other
        )
        const = default;

        constexpr
        Boolean
        isHybridDiscontinuousGalerkin()
        const
        {
            return true;
        }

        constexpr
        Boolean
        isHdg()
        const
        {
            return stabilization_ == Stabilization::Hdg;
        }

        constexpr
        Boolean
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
        Integer
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
        Boolean
        operator==(
            GeneralizedStrain const & other
        )
        const = default;

        constexpr
        Boolean
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

        template<Integer t_i>
        constexpr
        Mapping
        getMapping()
        const
        {
            return utility::get<t_i>(mappings_);
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
        Integer
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
        Boolean
        operator==(
            Behavior const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Behavior const & other
        )
        const = default;

        template<Integer t_i>
        constexpr
        std::tuple_element_t<t_i, std::tuple<t_GeneralizedStrains...>> const &
        getGeneralizedStrain()
        const
        {
            return utility::get<t_i>(generalized_strains_);
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
        Boolean
        operator==(
            FiniteElementMethod const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            FiniteElementMethod const & other
        )
        const = default;

        constexpr
        Boolean
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

        template<Integer t_i>
        constexpr
        Mapping
        getMapping()
        const
        {
            return utility::get<t_i>(generalized_strain_.mappings_);
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
    concept FiniteElementMethodConcept = detail::IsFiniteElementMethod<std::decay_t<t_T>>::value;

} // namespace lolita

#endif /* A4BCD9B5_985A_4D19_B3E9_7C559F45A353 */
