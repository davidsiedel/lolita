//
// Created by dsiedel on 13/04/2022.
//

#ifndef LOLITA_LOLITA_LOLITA_HXX
#define LOLITA_LOLITA_LOLITA_HXX

#include "lolita/lolita.hxx"

namespace lolita
{

    template<typename T>
    static constexpr
    auto
    isNumeric()
    {
        return std::is_arithmetic_v<T>;
    }

    template<typename T>
    using Raww = std::remove_cvref_t<T>;

    template<typename T, typename... U>
    concept IsSameAs = std::conjunction_v<std::is_same<T, U>...>;

    template<typename T, typename... U>
    concept IsContainedIn = std::disjunction_v<std::is_same<T, U>...>;

    template<typename T, typename... U>
    static constexpr
    auto
    areSame()
    {
        return std::conjunction_v<std::is_same<T, U>...>;
    }

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
    isIn()
    {
        return std::disjunction_v<std::is_same<T, U>...>;
    }

    template<typename T, typename... U>
    static constexpr
    auto
    isIn(
            T const &
            val,
            U &&...
            values
    )
    {
        auto res = Bool(false);
        auto set_res = [& val, & res] (auto && x) constexpr mutable {
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

#endif //LOLITA_LOLITA_LOLITA_HXX
