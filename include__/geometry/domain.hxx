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

    /**
     * @brief 
     * 
     * @tparam Frame_ 
     * @tparam T_ 
     * @tparam U_ 
     */
    template<FrameConcept Frame_, template<typename, typename...> typename T_ = TypeView2, typename... U_>
    struct DomainCollection
    {
        
        /**
         * @brief 
         * 
         */
        using Components = tuple_slice_t<Domains<T_, U_...>, 0, Frame_::getDimEuclidean() + 1>;
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         */
        template<Integer i_>
        using Component = std::tuple_element_t<i_, Components>;

    private:

        /**
         * @brief 
         * 
         * @tparam Shape_ 
         */
        template<DomainConcept Domain_>
        struct DomainTraits
        {
            
            /**
             * @brief 
             * 
             * @return constexpr Boolean 
             */
            static constexpr
            Boolean
            hasCoordinate()
            {
                return Domain_::getDimDomain() < Frame_::getDimEuclidean() + 1;
            }

            /**
             * @brief 
             * 
             */
            static constexpr
            Integer coordinate_ = Domain_::getDimDomain();

        };

    public:

        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<Components>;
        }
    
        template<Integer i_>
        std::tuple_element_t<i_, Components> const &
        getComponent()
        const
        {
            return std::get<i_>(components_);
        }
    
        template<DomainConcept Domain_>
        std::tuple_element_t<DomainTraits<Domain_>::coordinate_, Components> const &
        getComponent()
        const
        requires(DomainTraits<Domain_>::hasCoordinate())
        {
            return std::get<DomainTraits<Domain_>::coordinate_>(components_);
        }
        
        template<Integer i_>
        std::tuple_element_t<i_, Components> &
        getComponent()
        {
            return std::get<i_>(components_);
        }
    
        template<DomainConcept Domain_>
        std::tuple_element_t<DomainTraits<Domain_>::coordinate_, Components> &
        getComponent()
        requires(DomainTraits<Domain_>::hasCoordinate())
        {
            return std::get<DomainTraits<Domain_>::coordinate_>(components_);
        }

    private:

        Components components_;

    };

} // namespace lolita::geometry


#endif /* F065DF76_9691_41CB_B3C5_B48C5C36066D */
