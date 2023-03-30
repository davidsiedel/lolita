#ifndef F065DF76_9691_41CB_B3C5_B48C5C36066D
#define F065DF76_9691_41CB_B3C5_B48C5C36066D

#include "config.hxx"
#include "geometry/frame.hxx"

namespace lolita::geometry
{

    /**
     * @brief 
     * 
     * @tparam T_ 
     */
    template<typename T_>
    concept DomainConcept = requires
    {

        { T_::getDimDomain() } -> std::same_as<lolita::Integer>;
        
    };
    
    template<Integer dim_>
    struct Domain
    {

        static constexpr
        Integer
        getDimDomain()
        {
            return dim_;
        }

    };

    using PointDomain = Domain<0>;
    using CurveDomain = Domain<1>;
    using FacetDomain = Domain<2>;
    using SolidDomain = Domain<3>;

    template<template<DomainConcept, typename...> typename T_, typename... U_>
    using Domains = std::tuple<
        T_<Domain<0>, U_...>,
        T_<Domain<1>, U_...>,
        T_<Domain<2>, U_...>,
        T_<Domain<3>, U_...>
    >;

    template<FrameConcept Frame_>
    struct DomainLibrary
    {

        template<template<typename, typename...> typename T_, typename... U_>
        using MyDomains = tuple_slice_t<Domains<T_, U_...>, 0, Frame_::getDimEuclidean() + 1>;
        
        template<Integer i_>
        using Domain = std::tuple_element_t<i_, MyDomains<TypeView>>::type;
        
        template<DomainConcept Domain_>
        static constexpr
        Boolean
        hasDomain()
        {
            return Domain_::getDimDomain() <= Frame_::getDimEuclidean();
        }
        
        // template<Integer... t_i>
        static constexpr
        Integer
        getNumDomains()
        // requires(sizeof...(t_i) == 0)
        {
            return std::tuple_size_v<MyDomains<TypeView>>;
        }

    };

} // namespace lolita::geometry


#endif /* F065DF76_9691_41CB_B3C5_B48C5C36066D */
