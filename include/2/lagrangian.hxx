#ifndef B8D8CED7_F054_41EE_8A5E_A9C0C5887C9C
#define B8D8CED7_F054_41EE_8A5E_A9C0C5887C9C

#include "config.hxx"
#include "utility.hxx"

#include "2/label.hxx"
#include "2/quadrature.hxx"
#include "2/field.hxx"
#include "2/potential.hxx"

namespace lolita
{

    template<typename Implementation_, PotentialConcept... Potential_>
    struct LagrangianBase;

    namespace detail
    {
        
        struct LagrangianType
        {
            
        private:

            constexpr
            LagrangianType()
            {}

            constexpr
            Boolean
            operator==(
                LagrangianType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                LagrangianType const & other
            )
            const
            = default;

            template<typename Implementation_, PotentialConcept... Potential_>
            friend class lolita::LagrangianBase;

        };
        
    } // namespace detail

    template<typename T>
    concept LagrangianConcept = std::derived_from<std::decay_t<T>, detail::LagrangianType>;

    template<typename Implementation_, PotentialConcept... Potential_>
    struct LagrangianBase : detail::LagrangianType
    {
        
        static constexpr
        Integer
        getNumPotentials()
        {
            return sizeof...(Potential_);
        }

        using Potentials = utility::Aggregate<Potential_...>;

        constexpr
        LagrangianBase(
            std::basic_string_view<Character> && label,
            Potential_ const &... potentials
        )
        :
        detail::LagrangianType(),
        label_(std::forward<std::basic_string_view<Character>>(label)),
        potentials_(potentials...)
        {}

        constexpr
        Boolean
        operator==(
            LagrangianBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            LagrangianBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            LagrangianBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            LagrangianBase<T...> const & other
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

        template<Integer t_i>
        constexpr
        utility::aggregate_element_t<t_i, Potentials> const &
        getPotential()
        const
        {
            return utility::get<t_i>(potentials_);
        }

        constexpr
        Potentials const &
        getPotentials()
        const
        {
            return potentials_;
        }

        Label label_;

        Potentials potentials_;

    };

    template<PotentialConcept... Potential_>
    struct Lagrangian : LagrangianBase<Lagrangian<Potential_...>, Potential_...>
    {

    private:

        using Base_ = LagrangianBase<Lagrangian<Potential_...>, Potential_...>;

    public:

        constexpr
        Lagrangian(
            std::basic_string_view<Character> && label,
            Potential_ const &... potentials
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), potentials...)
        {}

    };
    
} // namespace lolita


#endif /* B8D8CED7_F054_41EE_8A5E_A9C0C5887C9C */
