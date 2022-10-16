#ifndef CA28089D_B58E_41D5_86C7_5B3F507D8B1A
#define CA28089D_B58E_41D5_86C7_5B3F507D8B1A

#include "config.hxx"
#include "utility.hxx"

#include "2/label.hxx"
#include "2/quadrature.hxx"
#include "2/field.hxx"

namespace lolita
{

    template<typename Implementation_, QuadratureConcept Quadrature_>
    struct PotentialBase;

    namespace detail
    {
        
        struct PotentialType
        {

        private:

            constexpr
            PotentialType()
            {}

            constexpr
            Boolean
            operator==(
                PotentialType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                PotentialType const & other
            )
            const
            = default;

            template<typename Implementation_, QuadratureConcept Quadrature_>
            friend class lolita::PotentialBase;

        };
        
    } // namespace detail

    template<typename T>
    concept PotentialConcept = std::derived_from<std::decay_t<T>, detail::PotentialType>;

    template<typename Implementation_, QuadratureConcept Quadrature_>
    struct PotentialBase : detail::PotentialType
    {

        using Quadrature = Quadrature_;

    protected:

        constexpr
        PotentialBase(
            Integer dim_domain,
            Quadrature quadrature
        )
        :
        detail::PotentialType(),
        label_(),
        dim_domain_(dim_domain),
        quadrature_(quadrature)
        {}

        constexpr
        PotentialBase(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Quadrature quadrature
        )
        :
        detail::PotentialType(),
        label_(std::forward<std::basic_string_view<Character>>(label)),
        dim_domain_(dim_domain),
        quadrature_(quadrature)
        {}

    public:

        constexpr
        Boolean
        operator==(
            PotentialBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            PotentialBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            PotentialBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            PotentialBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

        constexpr
        Label const &
        getLabel()
        const
        {
            return label_;
        }

        constexpr
        Integer
        getDimDomain()
        const
        {
            return dim_domain_;
        }

        constexpr
        Quadrature_ const &
        getQuadrature()
        const
        {
            return quadrature_;
        }

        Label label_;

        Integer dim_domain_;

        Quadrature quadrature_;

    };

    template<typename Quadrature_, typename... LinearOperator_>
    struct InternalPotential : PotentialBase<InternalPotential<Quadrature_, LinearOperator_...>, Quadrature_>
    {

    private:

        using Base_ = PotentialBase<InternalPotential<Quadrature_, LinearOperator_...>, Quadrature_>;

    public:

        using Quadrature = typename Base_::Quadrature;
        
        static constexpr
        Integer
        getNumMappings()
        {
            return sizeof...(LinearOperator_);
        }

    private:

        using LinearOperators_ = utility::Aggregate<LinearOperator_...>;

    public:

        constexpr
        InternalPotential(
            Integer dim_domain,
            Quadrature_ quadrature,
            LinearOperator_ const &... linear_operators
        )
        :
        Base_(dim_domain, quadrature),
        linear_operators_(linear_operators...)
        {}

        constexpr
        InternalPotential(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Quadrature_ quadrature,
            LinearOperator_ const &... linear_operators
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, quadrature),
        linear_operators_(linear_operators...)
        {}

        constexpr
        LinearOperators_ const &
        getLinearOperators()
        const
        {
            return linear_operators_;
        }

        template<Integer i>
        constexpr
        std::tuple_element_t<i, std::tuple<LinearOperator_...>> const &
        getLinearOperator()
        const
        {
            return utility::get<i>(linear_operators_);
        }

        LinearOperators_ linear_operators_;

    };

    namespace detail
    {

        template<typename... T>
        struct IsInternalPotentialTraits : std::false_type
        {};

        template<typename... T>
        struct IsInternalPotentialTraits<InternalPotential<T...>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept InternalPotentialConcept = detail::IsInternalPotentialTraits<T>::value;
    
    template<typename Quadrature_, typename UnknownField_>
    struct ExternalPotential : PotentialBase<ExternalPotential<Quadrature_, UnknownField_>, Quadrature_>
    {

        using UnknownField = UnknownField_;

    private:

        using Base_ = PotentialBase<ExternalPotential<Quadrature_, UnknownField_>, Quadrature_>;

    public:

        using Quadrature = typename Base_::Quadrature;

        constexpr
        ExternalPotential(
            Integer dim_domain,
            Quadrature quadrature,
            UnknownField const & unknown_field,
            ImposedField const & imposed_field
        )
        :
        Base_(dim_domain, quadrature),
        unknown_field_(unknown_field),
        imposed_field_(imposed_field)
        {}

        constexpr
        ExternalPotential(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Quadrature quadrature,
            UnknownField const & unknown_field,
            ImposedField const & imposed_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, quadrature),
        unknown_field_(unknown_field),
        imposed_field_(imposed_field)
        {}

        constexpr
        UnknownField const &
        getUnknownField()
        const
        {
            return unknown_field_;
        }

        constexpr
        ImposedField const &
        getImposedField()
        const
        {
            return imposed_field_;
        }

        UnknownField unknown_field_;

        ImposedField imposed_field_;

    };

    namespace detail
    {

        template<typename... T>
        struct IsExternalPotentialTraits : std::false_type
        {};

        template<typename... T>
        struct IsExternalPotentialTraits<ExternalPotential<T...>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept ExternalPotentialConcept = detail::IsExternalPotentialTraits<T>::value;

} // namespace lolita

#endif /* CA28089D_B58E_41D5_86C7_5B3F507D8B1A */
