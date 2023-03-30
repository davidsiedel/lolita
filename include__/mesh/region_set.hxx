/**
 * @file region_map.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-03-29
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef A0FD5970_E642_4680_BA8F_71D94646EABA
#define A0FD5970_E642_4680_BA8F_71D94646EABA

#include "geometry/shape.hxx"
#include "geometry/point.hxx"
#include "geometry/domain.hxx"
#include "geometry/frame.hxx"
#include "mesh/region.hxx"


namespace lolita::mesh
{

    template<template<geometry::DomainConcept, geometry::FrameConcept, typename...> typename T_, geometry::FrameConcept Frame_, typename... V_>
    struct RegionMap
    {

    private:

        template<typename... U_>
        // using RegionMap_ = std::unordered_map<std::basic_string<Character>, std::shared_ptr<T_<U_...>>>;
        using RegionMap_ = std::unordered_map<std::basic_string<Character>, T_<U_...>>;

        using Regions_ = typename geometry::DomainLibrary<Frame_>::template MyDomains<RegionMap_, Frame_, V_...>;

    public:

        RegionMap()
        {}
    
        template<Integer i_>
        std::tuple_element_t<i_, Regions_> const &
        getDomains()
        const
        {
            return std::get<i_>(domains_);
        }
        
        template<Integer i_>
        std::tuple_element_t<i_, Regions_> &
        getDomains()
        {
            return std::get<i_>(domains_);
        }
        
        Regions_ domains_;
        
    };

    template<template<geometry::DomainConcept, geometry::FrameConcept, typename...> typename T_, geometry::FrameConcept Frame_, typename... V_>
    struct RegionSet
    {

    private:

        template<typename... U_>
        // using RegionSet_ = std::vector<std::shared_ptr<T_<U_...>>>;
        using RegionSet_ = std::vector<T_<U_...>>;

        using Regions_ = typename geometry::DomainLibrary<Frame_>::template MyDomains<RegionSet_, Frame_, V_...>;

    public:

        RegionSet()
        {}
    
        template<Integer i_>
        std::tuple_element_t<i_, Regions_> const &
        getDomains()
        const
        {
            return std::get<i_>(domains_);
        }
        
        template<Integer i_>
        std::tuple_element_t<i_, Regions_> &
        getDomains()
        {
            return std::get<i_>(domains_);
        }
        
        Regions_ domains_;
        
    };
    
} // namespace lolita::mesh


#endif /* A0FD5970_E642_4680_BA8F_71D94646EABA */
