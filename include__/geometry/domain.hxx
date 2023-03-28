#ifndef F065DF76_9691_41CB_B3C5_B48C5C36066D
#define F065DF76_9691_41CB_B3C5_B48C5C36066D

#include "config.hxx"
#include "geometry/frame.hxx"

namespace lolita::geometry
{

    template<typename Domain_>
    struct DomainTraits
    {

        static constexpr
        Integer
        getDimension()
        {
            return Domain_::dimension_;
        }

    };
    
    template<Integer dim_>
    struct Domain
    {

        static constexpr
        Integer dimension_ = dim_;

    };

    template<template<typename, typename...> typename T_, typename... U_>
    using Domains = std::tuple<
        T_<Domain<0>, U_...>,
        T_<Domain<1>, U_...>,
        T_<Domain<2>, U_...>,
        T_<Domain<3>, U_...>
    >;

    template<typename Domain_, FrameConcept Frame_>
    struct FiniteDomain
    {

        explicit
        FiniteDomain(
            std::basic_string<Character> const & tag
        )
        :
        tag_(tag)
        {}

        explicit
        FiniteDomain(
            std::basic_string<Character> && tag
        )
        :
        tag_(std::move(tag))
        {}

        std::basic_string<Character> const &
        getLabel()
        const
        {
            return tag_;
        }

    private:

        std::basic_string<Character> tag_;

    };

} // namespace lolita::geometry


#endif /* F065DF76_9691_41CB_B3C5_B48C5C36066D */
