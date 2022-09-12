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

#include "utility.hxx"
#include "config.hxx"
#include "numerics.hxx"
#include "algebra.hxx"

namespace lolita
{
    
    using Point = Vector<Real, 3>;

    template<typename t_T>
    concept PointConcept = VectorConcept<t_T, Real, 3>;

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
        
        constexpr
        Boolean
        isGauss(
            Integer ord
        )
        const
        {
            return rule_ == Rule::Gauss && ord_ == ord;
        }
        
        constexpr
        Boolean
        hasOrd(
            Integer ord
        )
        const
        {
            return ord_ == ord;
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
        Polynomial
        getPolynomial()
        const
        {
            return polynomial_;
        }
        
        constexpr
        Boolean
        isMonomial()
        const
        {
            return polynomial_ == Polynomial::Monomial;
        }
        
        constexpr
        Boolean
        isMonomial(
            Integer order
        )
        const
        {
            return polynomial_ == Polynomial::Monomial && ord_ == order;
        }

        constexpr
        Integer
        getOrd()
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

        static constexpr
        Mapping
        gradient()
        {
            return Mapping(Type::Gradient);
        }

        static constexpr
        Mapping
        identity()
        {
            return Mapping(Type::Identity);
        }

        static constexpr
        Mapping
        divergence()
        {
            return Mapping(Type::Divergence);
        }

        static constexpr
        Mapping
        smallStrain()
        {
            return Mapping(Type::SmallStrain);
        }

        static constexpr
        Mapping
        largeStrain()
        {
            return Mapping(Type::LargeStrain);
        }

    private:

        constexpr
        Mapping(
            Type type
        )
        :
        type_(type),
        tag_(-1)
        {}

        constexpr
        Mapping(
            Type type,
            Integer tag
        )
        :
        type_(type),
        tag_(tag)
        {}

    public:

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

        constexpr
        Integer
        getTag()
        const
        {
            return tag_;
        }

        constexpr
        Boolean
        isGradient()
        const
        {
            return type_ == Type::Gradient;
        }

        constexpr
        Boolean
        isIdentity()
        const
        {
            return type_ == Type::Identity;
        }

        constexpr
        Boolean
        isDivergence()
        const
        {
            return type_ == Type::Divergence;
        }

        constexpr
        Boolean
        isSmallStrain()
        const
        {
            return type_ == Type::SmallStrain;
        }

        constexpr
        Boolean
        isLargeStrain()
        const
        {
            return type_ == Type::LargeStrain;
        }

        Type type_;

        Integer tag_;

    };

    struct Field
    {

        static constexpr
        Field
        tensor(
            Integer dim
        )
        {
            return Field(dim);
        }

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

    private:
        
        constexpr
        Field(
            Integer dim
        )
        :
        dim_(dim),
        tag_(-1)
        {}
        
        constexpr
        Field(
            Integer dim,
            Integer tag
        )
        :
        dim_(dim),
        tag_(tag)
        {}

    public:

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
        Integer
        getDim()
        const
        {
            return dim_;
        }

        constexpr
        Integer
        getTag()
        const
        {
            return tag_;
        }

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

        Integer tag_;

    };

    struct HybridDiscontinuousGalerkin
    {

        enum Stabilization
        {

            Hdg,
            HybridHighOrder,

        };

        static constexpr
        HybridDiscontinuousGalerkin
        hybridDiscontinuousGalerkin(
            Integer ord_cell,
            Integer ord_face
        )
        {
            return HybridDiscontinuousGalerkin(Basis::monomial(ord_cell), Basis::monomial(ord_face), Basis::monomial(ord_face), Hdg);
        }

        static constexpr
        HybridDiscontinuousGalerkin
        hybridHighOrder(
            Integer ord_cell,
            Integer ord_face
        )
        {
            return HybridDiscontinuousGalerkin(Basis::monomial(ord_cell), Basis::monomial(ord_face), Basis::monomial(ord_face), HybridHighOrder);
        }

        constexpr
        HybridDiscontinuousGalerkin(
            Basis cell_basis,
            Basis face_basis,
            Basis grad_basis,
            Stabilization stabilization
        )
        :
        cell_basis_(cell_basis),
        face_basis_(face_basis),
        grad_basis_(grad_basis),
        stabilization_(stabilization)
        {}

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
        Basis
        getCellBasis()
        const
        {
            return cell_basis_;
        }

        constexpr
        Basis
        getFaceBasis()
        const
        {
            return face_basis_;
        }

        constexpr
        Basis
        getGradBasis()
        const
        {
            return grad_basis_;
        }

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
            return stabilization_ == Stabilization::HybridHighOrder;
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
            Integer tag,
            Field field,
            t_Mappings... mappings
        )
        :
        tag_(tag),
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
        Integer
        getTag()
        const
        {
            return tag_;
        }

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

        Integer tag_;

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
            Integer tag,
            t_GeneralizedStrains... generalized_strains
        )
        :
        tag_(tag),
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

        constexpr
        Integer
        getTag()
        const
        {
            return tag_;
        }

        template<Integer t_i>
        constexpr
        std::tuple_element_t<t_i, std::tuple<t_GeneralizedStrains...>> const &
        getGeneralizedStrain()
        const
        {
            return utility::get<t_i>(generalized_strains_);
        }

        Integer tag_;

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

    struct Strategy
    {

        enum Type
        {
            EigenLU,
            EigenLDLT,
        };

        static constexpr
        Strategy
        eigenLU()
        {
            return Strategy(EigenLU);
        }

        static constexpr
        Strategy
        eigenLU(
            Integer tag
        )
        {
            return Strategy(EigenLU, tag);
        }

    private:

        constexpr
        Strategy(
            Type type
        )
        :
        type_(type),
        tag_(-1)
        {}

        constexpr
        Strategy(
            Type type,
            Integer tag
        )
        :
        type_(type),
        tag_(tag)
        {}

    public:

        Type type_;

        Integer tag_;

    };

    template<Strategy t_s>
    struct LinearSystem
    {

        static
        std::unique_ptr<LinearSystem>
        make_unique()
        {
            return std::make_unique<LinearSystem>();
        }

        LinearSystem()
        :
        size_(0)
        {}

        std::atomic<Natural> const &
        getSize()
        const
        {
            return size_;
        }

        std::atomic<Natural> &
        getSize()
        {
            return size_;
        }

    private:

        std::atomic<Natural> size_;

    };

} // namespace lolita

#endif /* A4BCD9B5_985A_4D19_B3E9_7C559F45A353 */
