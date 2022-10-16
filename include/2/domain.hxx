#ifndef C85E0DF8_E510_4C08_859A_1EA578E62679
#define C85E0DF8_E510_4C08_859A_1EA578E62679

#include "config.hxx"

namespace lolita
{

    template<typename Implementation_>
    struct DomainBase;

    namespace detail
    {
        
        struct DomainType
        {

        private:

            constexpr
            DomainType()
            {}

            constexpr
            Boolean
            operator==(
                DomainType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                DomainType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::DomainBase;

        };
        
    } // namespace detail

    template<typename T>
    concept DomainConcept = std::derived_from<std::decay_t<T>, detail::DomainType>;
    
    template<typename Implementation_>
    struct DomainBase : detail::DomainType
    {

    protected:

        constexpr
        DomainBase(
            Integer dim
        )
        :
        detail::DomainType(),
        dim_(dim)
        {}

    public:

        constexpr
        Boolean
        operator==(
            DomainBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            DomainBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            DomainBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            DomainBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }

        Integer dim_;

    };

    struct PointDomain : DomainBase<PointDomain>
    {

    private:

        using Base_ = DomainBase<PointDomain>;

    public:

        constexpr
        PointDomain()
        :
        Base_(0)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsPointDomainTraits : std::false_type
        {};

        template<>
        struct IsPointDomainTraits<PointDomain> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept PointDomainConcept = detail::IsPointDomainTraits<T>::value;

    struct CurveDomain : DomainBase<CurveDomain>
    {

    private:

        using Base_ = DomainBase<CurveDomain>;

    public:

        constexpr
        CurveDomain()
        :
        Base_(1)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsCurveDomainTraits : std::false_type
        {};

        template<>
        struct IsCurveDomainTraits<CurveDomain> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept CurveDomainConcept = detail::IsCurveDomainTraits<T>::value;

    struct FacetDomain : DomainBase<FacetDomain>
    {

    private:

        using Base_ = DomainBase<FacetDomain>;

    public:

        constexpr
        FacetDomain()
        :
        Base_(1)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsFacetDomainTraits : std::false_type
        {};

        template<>
        struct IsFacetDomainTraits<FacetDomain> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept FacetDomainConcept = detail::IsFacetDomainTraits<T>::value;

    struct SolidDomain : DomainBase<SolidDomain>
    {

    private:

        using Base_ = DomainBase<SolidDomain>;

    public:

        constexpr
        SolidDomain()
        :
        Base_(1)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsSolidDomainTraits : std::false_type
        {};

        template<>
        struct IsSolidDomainTraits<SolidDomain> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept SolidDomainConcept = detail::IsSolidDomainTraits<T>::value;

    struct DomainLibrary
    {

        template<template<DomainConcept auto, auto...> typename T, auto... t_arg>
        using Domains = std::tuple<
            T<PointDomain{}, t_arg...>,
            T<CurveDomain{}, t_arg...>,
            T<FacetDomain{}, t_arg...>,
            T<SolidDomain{}, t_arg...>
        >;

    };

} // namespace lolita


#endif /* C85E0DF8_E510_4C08_859A_1EA578E62679 */
