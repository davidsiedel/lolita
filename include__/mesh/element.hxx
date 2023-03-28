#ifndef AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A
#define AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A

#include "geometry/frame.hxx"
#include "geometry/shape.hxx"
#include "geometry/domain.hxx"

namespace lolita::mesh
{

    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam t_domain 
     */
    template<template<typename, typename...> typename T_, FrameConcept Frame_, typename... U_>
    struct DomainMap
    {

    private:

        struct DomainsTraits
        {

            template<typename Domain_, typename... V_>
            using Neighbors = std::unordered_map<std::basic_string<Character>, std::shared_ptr<T_<Domain_, V_...>>>;

            using Type = tuple_slice_t<Domains<Neighbors, Frame_, U_...>, 0, Frame_::getDimEuclidean() + 1>;

        };

    public:

        using Domains = typename DomainsTraits::Type;

        DomainMap()
        {}
    
        template<Integer t_i>
        std::tuple_element_t<t_i, Domains> const &
        getDomains()
        const
        {
            return std::get<t_i>(domains_);
        }
        
        template<Integer t_i>
        std::tuple_element_t<t_i, Domains> &
        getDomains()
        {
            return std::get<t_i>(domains_);
        }
        
        Domains domains_;
        
    };

    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam t_domain 
     */
    template<template<typename, typename> typename T, typename t_domain>
    struct DomainSet
    {

    private:

        template<typename t__dim, typename t__domain>
        using t__RegionSet = std::vector<std::shared_ptr<T<t__dim, t__domain>>>;

        using t_RegionSet = lolita::utility::tuple_slice_t<typename DomainLibrary::Domains<t__RegionSet, t_domain>, 0, t_domain.getDim() + 1>;

    public:

        DomainSet()
        {}
    
        template<Integer t_i>
        std::tuple_element_t<t_i, t_RegionSet> const &
        getDomains()
        const
        {
            return std::get<t_i>(domains_);
        }
        
        template<Integer t_i>
        std::tuple_element_t<t_i, t_RegionSet> &
        getDomains()
        {
            return std::get<t_i>(domains_);
        }
        
        t_RegionSet domains_;
        
    };
    
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam t_domain 
     */
    template<template<ShapeConcept auto, typename> typename T, typename t_domain>
    struct ElementMap
    {

    private:

        template<ShapeConcept auto t_element, typename t__domain, auto... t__args>
        using t_ElementMap = std::unordered_map<std::basic_string<Character>, std::shared_ptr<T<t_element, t__domain, t__args...>>>;

        using t_Elements = lolita::utility::tuple_slice_t<typename ShapeLibrary::Elements<t_ElementMap, t_domain>, 0, t_domain.getDim() + 1>;

    public:

        ElementMap()
        {}
    
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> const &
        getElements()
        const
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> &
        getElements()
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        t_Elements elements_;

    };

    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam t_domain 
     */
    template<template<ShapeConcept auto, typename> typename T, typename t_domain>
    struct ElementSet
    {

    private:

        template<ShapeConcept auto t_element, typename t__domain>
        using t_ElementSet = std::vector<std::shared_ptr<T<t_element, t__domain>>>;

        using t_Elements = lolita::utility::tuple_slice_t<typename ShapeLibrary::Elements<t_ElementSet, t_domain>, 0, t_domain.getDim() + 1>;

    public:

        ElementSet()
        {}
    
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> const &
        getElements()
        const
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> &
        getElements()
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        t_Elements elements_;

    };
    
} // namespace lolita::mesh


#endif /* AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A */
