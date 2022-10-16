#ifndef C5BF0FEB_78CF_4CCE_9C79_32F158BA48F9
#define C5BF0FEB_78CF_4CCE_9C79_32F158BA48F9

#include "config.hxx"

namespace lolita
{

    template<typename Implementation_>
    struct BasisBase;

    namespace detail
    {
        
        struct BasisType
        {

        private:

            constexpr
            BasisType()
            {}

            constexpr
            Boolean
            operator==(
                BasisType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                BasisType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::BasisBase;

        };
        
    } // namespace detail

    template<typename T>
    concept BasisConcept = std::derived_from<std::decay_t<T>, detail::BasisType>;

    template<typename Implementation_>
    struct BasisBase : detail::BasisType
    {

    protected:
        
        constexpr
        BasisBase(
            Integer ord
        )
        :
        detail::BasisType(),
        ord_(ord)
        {}

    public:

        constexpr
        Boolean
        operator==(
            BasisBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            BasisBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            BasisBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            BasisBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

        constexpr
        Implementation_
        toOrder(
            Integer order,
            auto const &... args
        )
        const
        {
            return Implementation_(order, args...);
        }

        constexpr
        Integer
        getOrd()
        const
        {
            return ord_;
        }
        
        Integer ord_;

    };

    struct MonomialBasis : BasisBase<MonomialBasis>
    {

    private:

        using Base_ = BasisBase<MonomialBasis>;

    public:

        explicit constexpr
        MonomialBasis(
            Integer dim
        )
        :
        Base_(dim)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsMonomialBasisTraits : std::false_type
        {};

        template<>
        struct IsMonomialBasisTraits<MonomialBasis> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept MonomialBasisConcept = detail::IsMonomialBasisTraits<std::decay_t<T>>::value;

    struct LagrangeBasis : BasisBase<LagrangeBasis>
    {

    private:

        using Base_ = BasisBase<LagrangeBasis>;

    public:

        explicit constexpr
        LagrangeBasis(
            Integer dim
        )
        :
        Base_(dim)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsLagrangeBasisTraits : std::false_type
        {};

        template<>
        struct IsLagrangeBasisTraits<LagrangeBasis> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept LagrangeBasisConcept = detail::IsLagrangeBasisTraits<std::decay_t<T>>::value;
    
} // namespace lolita


#endif /* C5BF0FEB_78CF_4CCE_9C79_32F158BA48F9 */
