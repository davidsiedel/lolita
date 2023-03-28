#ifndef D58F56E2_01E1_4EFA_803C_662777D072C4
#define D58F56E2_01E1_4EFA_803C_662777D072C4

#include <iostream>
#include <fstream>
#include <execution>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <ostream>
#include <iomanip>
#include <filesystem>
#include <span>

namespace lolita
{

    namespace config
    {

        using Character =   char;

        using Integer =     int;

        using Index =       unsigned short;

        using Natural =     unsigned long long;

        using Real =        double;

        using Boolean =     bool;

    } // namespace config

    using namespace config;

    template<typename U_, typename...>
    struct TypeView
    {

        using type = U_;

    };

    /**
     * Hello
     * *********************************************************************************************************************************************************
     */

    template<typename... T_>
    struct tuple_concatenation_traits;

    template<typename T_, typename... U_>
    struct tuple_concatenation_traits<T_, U_...>
    {

        using type = typename tuple_concatenation_traits<T_, typename tuple_concatenation_traits<U_...>::type>::type;
        // using type = typename tuple_concatenation_traits<std::tuple<T_>, typename tuple_concatenation_traits<U_...>::type>::type;

    };

    template<typename... T_, typename... U_>
    struct tuple_concatenation_traits<std::tuple<T_...>, std::tuple<U_...>>
    {

        using type = std::tuple<T_..., U_...>;

    };

    template<typename... T_>
    struct tuple_concatenation_traits<std::tuple<T_...>>
    {

        using type = std::tuple<T_...>;

    };

    template<template<typename...> typename T_, typename... U_, typename... V_>
    struct tuple_concatenation_traits<T_<U_...>, T_<V_...>>
    {

        using type = T_<U_..., V_...>;

    };

    template<template<typename...> typename T_, typename... U_>
    struct tuple_concatenation_traits<T_<U_...>>
    {

        using type = T_<U_...>;

    };

    template<typename... T_>
    using tuple_cat_t = typename tuple_concatenation_traits<T_...>::type;

    template<typename... T_>
    struct tuple_merge_traits;

    template<typename T_, typename... U_>
    struct tuple_merge_traits<T_, U_...>
    {

        using type = typename tuple_merge_traits<std::tuple<T_>, typename tuple_merge_traits<U_...>::type>::type;
        
    };

    template<typename... T_, typename... U_>
    struct tuple_merge_traits<std::tuple<T_...>, U_...>
    {

        using type = typename tuple_merge_traits<std::tuple<>, typename tuple_merge_traits<T_..., U_...>::type>::type;
        
    };
    
    template<typename... T_, typename... U_>
    struct tuple_merge_traits<std::tuple<T_...>, std::tuple<U_...>>
    {

        using type = std::tuple<T_..., U_...>;
        
    };
    
    template<typename T_>
    struct tuple_merge_traits<std::tuple<T_>>
    {

        using type = T_;
        
    };

    template<typename T_>
    using tuple_merge_t = typename tuple_merge_traits<T_>::type;

    template<typename T_, Integer i_, Integer j_>
    struct tuple_slice_traits;

    template<typename... T_, Integer i_, Integer j_>
    struct tuple_slice_traits<std::tuple<T_...>, i_, j_>
    {

        using type = tuple_cat_t<std::tuple<std::tuple_element_t<i_, std::tuple<T_...>>>, typename tuple_slice_traits<std::tuple<T_...>, i_ + 1, j_>::type>;

    };
    
    template<typename... T_, Integer i_, Integer j_>
    requires(i_ == j_ - 1)
    struct tuple_slice_traits<std::tuple<T_...>, i_, j_>
    {

        using type = std::tuple<std::tuple_element_t<j_ - 1, std::tuple<T_...>>>;

    };
    
    template<typename... T_, Integer i_, Integer j_>
    requires(i_ == j_)
    struct tuple_slice_traits<std::tuple<T_...>, i_, j_>
    {

        using type = std::tuple<>;

    };

    template<typename T_, Integer i_, Integer j_>
    using tuple_slice_t = typename tuple_slice_traits<T_, i_, j_>::type;

} // namespace lolita

#endif /* D58F56E2_01E1_4EFA_803C_662777D072C4 */
