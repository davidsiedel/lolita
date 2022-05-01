//
// Created by dsiedel on 13/04/2022.
//

#ifndef FETA_FETA_FETA_TYPES_HXX
#define FETA_FETA_FETA_TYPES_HXX

#include <iostream>
#include <tuple>

namespace lolita
{

    /**
     * @brief
     */
    struct Void{};

    /**
     * @brief
     */
    using Char = char;

    /**
     * @brief
     */
    using Intg = int;

    /**
     * @brief
     */
    using Indx = std::size_t;

    /**
     * @brief
     */
    using Long = unsigned long;

    /**
     * @brief
     */
    using Real = double;

    /**
     * @brief
     */
    using Strg = std::string;

    /**
     * @brief
     */
    using Bool = bool;

    /**
     * @brief
     */
    using StrgStream = std::stringstream;

    /**
     * @brief
     */
    using OuptStream = std::ostream;

    /**
     * @brief
     */
    using InptStream = std::istream;

    /**
     * @brief
     */
    using FileStream = std::ifstream;

    template<typename T>
    using Raww = std::remove_cvref_t<T>;

    template<typename T, typename... U>
    concept IsSameAs = std::conjunction_v<std::is_same<T, U>...>;

    template<typename T, typename... U>
    concept IsContainedIn = std::disjunction_v<std::is_same<T, U>...>;

    template<typename T>
    concept NaturalType = IsSameAs<Raww<T>, Indx>;

    template<typename T>
    concept IntegerType = IsSameAs<Raww<T>, Intg, Indx>;

    template<typename T>
    concept RealType = IsSameAs<Raww<T>, Real, Intg, Indx>;

    template<typename T, typename... U>
    static constexpr
    auto
    areSame(
            T &&
            val,
            U &&...
            values
    )
    {
        auto res = Bool(true);
        auto set_res = [& val, & res] (auto && x) constexpr mutable {
            if (val != static_cast<T>(x)) res = false;
        };
        (set_res(std::forward<U>(values)), ...);
        return res;
    }

    template<typename T, typename... U>
    static constexpr
    auto
    isIn(
            T &&
            val,
            U &&...
            values
    )
    {
        auto res = Bool(false);
        auto set_res = [&] (auto && x) constexpr mutable {
            if (val == static_cast<T>(x)) res = true;
        };
        (set_res(std::forward<U>(values)), ...);
        return res;
    }

    template<typename T, typename... U>
    inline static
    void
    print(
            std::ostream &
            out,
            T &&
            arg,
            U &&...
            args
    )
    {
        out << std::forward<T>(arg);
        ((out << " " << std::forward<U>(args)), ...);
    }

    template<typename T, typename... U>
    inline static
    void
    print(
            T &&
            arg,
            U &&...
            args
    )
    {
        std::cout << std::forward<T>(arg);
        ((std::cout << " " << std::forward<U>(args)), ...);
        std::cout << std::endl;
    }

}

#endif //FETA_FETA_FETA_TYPES_HXX
