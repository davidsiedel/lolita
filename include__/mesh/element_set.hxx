#ifndef F4B914C6_936A_4791_8796_770B48F19366
#define F4B914C6_936A_4791_8796_770B48F19366

#include "geometry/shape.hxx"
#include "geometry/frame.hxx"
#include "geometry/domain.hxx"
#include "geometry/point.hxx"
#include "mesh/element.hxx"

namespace lolita::mesh
{

    template<template<geometry::ShapeConcept, geometry::FrameConcept, typename...> typename T_, geometry::FrameConcept Frame_, typename... V_>
    struct ElementMap
    {

    private:

        template<typename... U_>
        // using ElementMap_ = std::unordered_map<std::basic_string<Character>, std::shared_ptr<T_<U_...>>>;
        using ElementMap_ = std::unordered_map<std::basic_string<Character>, T_<U_...>>;

        using Elements_ = typename geometry::ShapeLibraryTraits<Frame_>::template MyShapes<ElementMap_, Frame_, V_...>;

    public:

        ElementMap()
        {}
    
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Elements_>> const &
        getElements()
        const
        {
            return std::get<j_>(std::get<i_>(domains_));
        }
        
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Elements_>> &
        getElements()
        {
            return std::get<j_>(std::get<i_>(domains_));
        }
        
        Elements_ domains_;
        
    };

    template<template<geometry::ShapeConcept, geometry::FrameConcept, typename...> typename T_, geometry::FrameConcept Frame_, typename... V_>
    struct ElementSet
    {

    private:

        template<typename... U_>
        // using ElementSet_ = std::vector<std::shared_ptr<T_<U_...>>>;
        using ElementSet_ = std::vector<T_<U_...>>;

        using Elements_ = typename geometry::ShapeLibraryTraits<Frame_>::template MyShapes<ElementSet_, Frame_, V_...>;

    public:

        ElementSet()
        {}
    
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Elements_>> const &
        getElements()
        const
        {
            return std::get<j_>(std::get<i_>(domains_));
        }
        
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Elements_>> &
        getElements()
        {
            return std::get<j_>(std::get<i_>(domains_));
        }
        
        Elements_ domains_;
        
    };

} // namespace lolita::mesh


#endif /* F4B914C6_936A_4791_8796_770B48F19366 */
