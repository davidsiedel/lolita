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

        constexpr
        lolita::integer
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

        template<GeneralizedStrainConcept auto t_generalized_strain>
        constexpr
        lolita::integer
        getGeneralizedStrainIndex()
        {
            auto index = lolita::integer(0);
            auto found = lolita::boolean(false);
            auto set_index = [&] <lolita::integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                if (lolita::utility::areEqual(t_generalized_strain, lolita::utility::get<t_i>(generalized_strains_)))
                {
                    found = true;
                }
                if (!found)
                {
                    index += 1;
                }
                if constexpr (t_i < getNumGeneralizedStrains() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_index(set_index);
            return index;
        }

        constexpr
        GeneralizedStrains const &
        getGeneralizedStrains()
        const
        {
            return generalized_strains_;
        }

        template<lolita::integer t_i>
        constexpr
        std::tuple_element_t<t_i, std::tuple<t_GeneralizedStrains...>> const &
        getGeneralizedStrain()
        const
        {
            return generalized_strains_.template get<t_i>();
        }

        template<GeneralizedStrainConcept auto t_generalized_strain>
        constexpr
        std::tuple_element_t<getGeneralizedStrainIndex<t_generalized_strain>(), std::tuple<t_GeneralizedStrains...>> const &
        getGeneralizedStrain()
        const
        {
            return generalized_strains_.template get<getGeneralizedStrainIndex<t_generalized_strain>()>();
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
    concept FiniteElementMethodConcept = detail::IsFiniteElementMethod<std::decay_t<t_T>>::value;

    using ParameterFunction = std::function<lolita::real(Point const &, lolita::real const &)>;

    struct MgisParameter
    {

        MgisParameter(
            std::basic_string_view<lolita::character> parameter_tag,
            ParameterFunction && function
        )
        :
        parameter_tag_(parameter_tag),
        function_(std::forward<ParameterFunction>(function))
        {}

        constexpr
        lolita::boolean
        operator==(
            MgisParameter const & other
        )
        const
        {
            return other.parameter_tag_ == this->parameter_tag_;
        }

        constexpr
        lolita::boolean
        operator!=(
            MgisParameter const & other
        )
        const
        {
            return !(* this == other);
        }

        std::basic_string_view<lolita::character> parameter_tag_;

        ParameterFunction function_;

    };

    struct MgisBehaviourData
    {

        mgis::behaviour::Behaviour behaviour_;

        std::vector<MgisParameter> parameters_;

    };

    struct MgisBehaviour
    {

        MgisBehaviour(
                std::basic_string_view<lolita::character> unknown_tag,
                std::basic_string_view<lolita::character> domain_tag,
                std::basic_string<lolita::character> const & path,
                std::basic_string<lolita::character> const & name,
                mgis::behaviour::Hypothesis hypothesis,
                std::vector<MgisParameter> && parameters
        )
        :
        unknown_tag_(unknown_tag),
        domain_tag_(domain_tag),
        behaviour_data_(std::make_shared<MgisBehaviourData>(MgisBehaviourData{mgis::behaviour::load(path, name, hypothesis), std::forward<std::vector<MgisParameter>>(parameters)}))
        {}

        // MgisBehaviour(
        //         std::basic_string<lolita::character> && unknown_tag,
        //         std::basic_string<lolita::character> && domain_tag,
        //         std::string const & path,
        //         std::string const & name,
        //         mgis::behaviour::Hypothesis hypothesis,
        //         mgis::behaviour::FiniteStrainBehaviourOptions finite_strain_behaviour_options,
        //         std::vector<lolita::behaviour::MgisParameter> && parameters
        // )
        // :
        // unknown_tag_(std::forward<std::basic_string<lolita::character>>(unknown_tag)),
        // domain_tag_(std::forward<std::basic_string<lolita::character>>(domain_tag)),
        // behaviour_data_(std::make_shared<lolita::behaviour::MgisBehaviourData>(lolita::behaviour::MgisBehaviourData{
        //     mgis::behaviour::load(finite_strain_behaviour_options, path, name, hypothesis),
        //     std::forward<std::vector<lolita::behaviour::MgisParameter>>(parameters)
        // }))
        // {}

        std::basic_string<lolita::character> unknown_tag_;

        std::basic_string<lolita::character> domain_tag_;

        std::shared_ptr<MgisBehaviourData> behaviour_data_;

    };

}


#endif /* AD31B868_1D56_48DB_A2CC_F898A668702D */
