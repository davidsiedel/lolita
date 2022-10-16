#ifndef A48E7CDD_8385_4FD9_B6CD_1EE89B816AA8
#define A48E7CDD_8385_4FD9_B6CD_1EE89B816AA8

#include "config.hxx"

namespace lolita
{

    template<typename Implementation_>
    struct QuadratureBase;

    namespace detail
    {
        
        struct QuadratureType
        {

        private:

            constexpr
            QuadratureType()
            {}

            constexpr
            Boolean
            operator==(
                QuadratureType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                QuadratureType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::QuadratureBase;

        };
        
    } // namespace detail

    template<typename T>
    concept QuadratureConcept = std::derived_from<std::decay_t<T>, detail::QuadratureType>;
    
    template<typename Implementation_>
    struct QuadratureBase : detail::QuadratureType
    {

    protected:
        
        constexpr explicit
        QuadratureBase(
            Integer ord
        )
        :
        detail::QuadratureType(),
        ord_(ord)
        {}

    public:

        constexpr
        Boolean
        operator==(
            QuadratureBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            QuadratureBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            QuadratureBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            QuadratureBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

        constexpr
        Integer
        getOrder()
        const
        {
            return ord_;
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
        
        Integer ord_;

    };

    struct GaussQuadrature : QuadratureBase<GaussQuadrature>
    {

    private:

        using Base_ = QuadratureBase<GaussQuadrature>;

    public:

        constexpr explicit
        GaussQuadrature(
            Integer order
        )
        :
        Base_(order)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsGaussQuadratureTraits : std::false_type
        {};

        template<>
        struct IsGaussQuadratureTraits<GaussQuadrature> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept GaussQuadratureConcept = detail::IsGaussQuadratureTraits<std::decay_t<T>>::value;

} // namespace lolita

#endif /* A48E7CDD_8385_4FD9_B6CD_1EE89B816AA8 */
