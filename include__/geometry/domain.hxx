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
        T_<PointDomain, U_...>,
        T_<CurveDomain, U_...>,
        T_<FacetDomain, U_...>,
        T_<SolidDomain, U_...>
    >;
    
    template<FrameConcept Frame_>
    struct DomainCollectionTraits
    {

        template<template<typename, typename...> typename T_ = TypeView, typename... U_>
        using Collection = tuple_slice_t<Domains<T_, U_...>, 0, Frame_::getDimEuclidean() + 1>;
        
        template<Integer i_>
        using Domain = std::tuple_element_t<i_, Collection<>>;

    private:
    
        template<DomainConcept Domain_>
        struct DomainTraits
        {
            
            static constexpr
            Boolean
            hasCoordinate()
            {
                return Domain_::getDimDomain() < Frame_::getDimEuclidean() + 1;
            }
            
            static constexpr
            Integer coordinate_ = Domain_::getDimDomain();

        };

    public:
    
        template<DomainConcept Domain_>
        static constexpr
        Integer
        getDomainIndex()
        requires(DomainTraits<Domain_>::hasCoordinate())
        {
            return DomainTraits<Domain_>::coordinate_;
        }

        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<Collection<>>;
        }

        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer i_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<Domain<i_>>();
                if constexpr (i_ < getNumComponents() - 1)
                {
                    apply_.template operator()<i_ + 1>(apply_);
                }
            };
            apply(apply);
        }

        template<Integer i_>
        static
        void
        apply(
            auto const & fun
        )
        {
            fun.template operator()<Domain<i_>>();
        }

    };

} // namespace lolita::geometry


#endif /* F065DF76_9691_41CB_B3C5_B48C5C36066D */
