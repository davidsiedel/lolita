#ifndef E5351E45_F2AD_48E4_8AB5_ED530F0DBCDA
#define E5351E45_F2AD_48E4_8AB5_ED530F0DBCDA

#include <iostream>
#include <fstream>
#include <execution>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include "lolita/lolita.hxx"
#include "lolita/lolita_algebra.hxx"

namespace lolita::utility
{

    template<auto t_value>
    struct Holder
    {

        auto static constexpr value_ = t_value;

    };

    template<typename t_T, typename t_U>
    static constexpr
    lolita::boolean
    areEqual(
        t_T const & x,
        t_U const & y
    )
    {
        if constexpr (std::is_same_v<t_T, t_U>)
        {
            if (x == y)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, typename t_U>
    struct tuple_concatenation_traits;

    template<typename... t_T, typename... t_U>
    struct tuple_concatenation_traits<std::tuple<t_T...>, std::tuple<t_U...>>
    {

        using type = std::tuple<t_T..., t_U...>;

    };

    template<typename t_T, typename t_U>
    using tuple_cat_t = tuple_concatenation_traits<t_T, t_U>::type;

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, lolita::integer t_a, lolita::integer t_b>
    struct TupleSliceTraits;

    template<typename... t_T, lolita::integer t_a, lolita::integer t_b>
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = tuple_cat_t<std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>, typename TupleSliceTraits<std::tuple<t_T...>, t_a + 1, t_b>::type>;

    };
    
    template<typename... t_T, lolita::integer t_a, lolita::integer t_b>
    requires(t_a == t_b - 1)
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = std::tuple<std::tuple_element_t<t_b - 1, std::tuple<t_T...>>>;

    };
    
    template<typename... t_T, lolita::integer t_a, lolita::integer t_b>
    requires(t_a == t_b)
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = std::tuple<>;

    };

    template<typename t_T, lolita::integer t_a, lolita::integer t_b>
    using tuple_slice_t = TupleSliceTraits<t_T, t_a, t_b>::type;

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, typename t_U>
    struct TupleHasTraits;

    template<typename... t_T, typename t_U>
    struct TupleHasTraits<std::tuple<t_T...>, t_U>
    {

        lolita::boolean static constexpr value = (std::is_same_v<t_U, t_T> || ...);

    };

    template<typename t_T, typename t_U>
    static constexpr
    lolita::boolean tuple_has_v = TupleHasTraits<t_T, t_U>::value;

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, lolita::integer t_a>
    struct TupleUniqueTraits;

    template<typename... t_T, lolita::integer t_a>
    requires(t_a < sizeof...(t_T) - 1 && !tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = tuple_cat_t<std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>, typename TupleUniqueTraits<std::tuple<t_T...>, t_a + 1>::type>;

    };

    template<typename... t_T, lolita::integer t_a>
    requires(t_a < sizeof...(t_T) - 1 && tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = tuple_cat_t<std::tuple<>, typename TupleUniqueTraits<std::tuple<t_T...>, t_a + 1>::type>;

    };

    template<typename... t_T, lolita::integer t_a>
    requires(t_a == sizeof...(t_T) - 1 && !tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>;

    };

    template<typename... t_T, lolita::integer t_a>
    requires(t_a == sizeof...(t_T) - 1 && tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = std::tuple<>;

    };

    template<typename t_T>
    using tuple_unique_t = TupleUniqueTraits<t_T, 0>::type;

    // ------------------------------------------------------------------------------------------------------

    template<typename>
    void
    TD()
    {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    // ------------------------------------------------------------------------------------------------------

    namespace detail
    {

        template<lolita::integer t_i, typename t_T>
        struct AggregateTail
        {

            constexpr
            lolita::boolean
            operator==(
                AggregateTail const &
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                AggregateTail const &
            )
            const = default;

            t_T value_;

        };

        template<lolita::integer t_i, typename... t_T>
        struct AggregateHead;

        template<lolita::integer t_i>
        struct AggregateHead<t_i>
        {

            constexpr
            lolita::boolean
            operator==(
                AggregateHead const &
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                AggregateHead const &
            )
            const = default;

        };

        template<lolita::integer t_i, typename t_T, typename... t_U>
        struct AggregateHead<t_i, t_T, t_U...> : AggregateTail<t_i, t_T>, AggregateHead<t_i + 1, t_U...>
        {

            constexpr
            lolita::boolean
            operator==(
                AggregateHead const &
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                AggregateHead const &
            )
            const = default;

        };

        template<lolita::integer t_i, typename t_T, typename... t_U>
        static constexpr
        t_T &
        get(
            AggregateHead<t_i, t_T, t_U...> & aggregate
        )
        {
            return aggregate.AggregateTail<t_i, t_T>::value_;
        }

        template<lolita::integer t_i, typename t_T, typename... t_U>
        static constexpr
        t_T const &
        get(
            AggregateHead<t_i, t_T, t_U...> const & aggregate
        )
        {
            return aggregate.AggregateTail<t_i, t_T>::value_;
        }

    }

    // ------------------------------------------------------------------------------------------------------

    template<typename... T>
    struct Aggregate : detail::AggregateHead<0, T...>
    {

        static constexpr
        lolita::integer
        getSize()
        {
            return sizeof...(T);
        }

        constexpr
        Aggregate()
        :
        detail::AggregateHead<0, T...>()
        {}

        constexpr
        Aggregate(
            T &&... args
        )
        :
        detail::AggregateHead<0, T...>{std::forward<T>(args)...}
        {}

        constexpr
        Aggregate(
            T const &... args
        )
        :
        detail::AggregateHead<0, T...>{args...}
        {}

        constexpr
        lolita::boolean
        operator==(
            Aggregate const &
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Aggregate const &
        )
        const = default;

        template<lolita::integer t_i>
        constexpr
        std::tuple_element_t<t_i, std::tuple<T...>> &
        get()
        {
            return detail::get<t_i>(* this);
        }

        template<lolita::integer t_i>
        constexpr
        std::tuple_element_t<t_i, std::tuple<T...>> const &
        get()
        const
        {
            return detail::get<t_i>(* this);
        }

    };

    template<lolita::integer t_i, typename... t_U>
    static constexpr
    std::tuple_element_t<t_i, std::tuple<t_U...>> &
    get(
        Aggregate<t_U...> & aggregate
    )
    {
        return detail::get<t_i>(aggregate);
    }

    template<lolita::integer t_i, typename... t_U>
    static constexpr
    std::tuple_element_t<t_i, std::tuple<t_U...>> const &
    get(
        Aggregate<t_U...> const & aggregate
    )
    {
        return detail::get<t_i>(aggregate);
    }

    template<typename T>
    struct AggregateSizeTraits;

    template<typename... T>
    struct AggregateSizeTraits<Aggregate<T...>>
    {

        lolita::integer static constexpr value = sizeof...(T);

    };

    template<typename T>
    static constexpr
    lolita::index aggregate_size_v = AggregateSizeTraits<T>::value;

    template<lolita::index I, typename T>
    struct AggregateElementTraits;

    template<lolita::index I, typename... T>
    struct AggregateElementTraits<I, Aggregate<T...>>
    {

        using type = std::tuple_element_t<I, std::tuple<T...>>;

    };

    template<lolita::index I, typename... T>
    using aggregate_element_t = typename AggregateElementTraits<I, T...>::type;

    template<template<auto> typename t_T, auto t_aggregate, lolita::integer t_i>
    struct aggregate_expansion_traits
    {

        using type = tuple_cat_t<std::tuple<t_T<get<t_i>(t_aggregate)>>, typename aggregate_expansion_traits<t_T, t_aggregate, t_i + 1>::type>;

    };
    
    template<template<auto> typename t_T, auto t_aggregate, lolita::integer t_i>
    requires(t_i == aggregate_size_v<std::decay_t<decltype(t_aggregate)>> - 1)
    struct aggregate_expansion_traits<t_T, t_aggregate, t_i>
    {

        using type = std::tuple<t_T<get<t_i>(t_aggregate)>>;

    };

    template<template<auto> typename t_T, auto t_aggregate>
    using aggregate_expansion_t = aggregate_expansion_traits<t_T, t_aggregate, 0>::type;

    // ------------------------------------------------------------------------------------------------------

    struct Label
    {

        using Tag = std::array<lolita::character, 50>;

    private:

        static constexpr
        Label::Tag
        setTag(
            std::basic_string_view<lolita::character> str
        )
        {
            auto tag = Label::Tag();
            auto count = lolita::integer(0);
            for (auto c : str) {
                tag[count] = c;
                count ++;
            }
            return tag;
        }

    public:

        constexpr
        Label(
            std::basic_string_view<lolita::character> str
        )
        :
        tag_(setTag(str))
        {}

        constexpr
        lolita::boolean
        operator==(
            Label const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Label const & other
        )
        const = default;

        constexpr inline
        lolita::boolean
        operator==(
            std::basic_string_view<lolita::character> str
        )
        const
        {
            return str == this->view();
        }

        constexpr inline
        lolita::boolean
        operator!=(
            std::basic_string_view<lolita::character> str
        )
        const
        {
            return !(* this == str);
        }

        constexpr inline
        std::basic_string_view<lolita::character>
        view()
        const
        {
            return std::basic_string_view<lolita::character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), lolita::character())));
        }

        friend inline
        std::ostream &
        operator<<(
            std::ostream & os,
            Label const & label
        )
        {
            os << label.view();
            return os;
        }

        Label::Tag tag_;

    };

    // ------------------------------------------------------------------------------------------------------

    // static void inline
    // removeCharacter(
    //         std::basic_string<lolita::character> & line,
    //         lolita::character character
    // )
    // {
    //     line.erase(std::remove(line.begin(), line.end(), character), line.end());
    // }

    struct File
    {

        File() = default;

        explicit
        File(
            std::basic_string_view<lolita::character> file_path
        )
        :
        file_path_(file_path)
        {
            readLines();
        }

        lolita::boolean
        operator==(
            File const &
            other
        )
        const = default;

        lolita::boolean
        operator!=(
            File const &
            other
        )
        const = default;

    private:

        inline
        void
        readLines()
        {
            std::basic_ifstream<lolita::character> file(file_path_);
            if (!file)
            {
                throw std::runtime_error("Could not open file");
            }
            for (std::basic_string<lolita::character> line; std::getline(file, line); )
            {
                lines_.push_back(line);
            }
        }

    public:

        std::vector<std::basic_string<lolita::character>> lines_;

        std::basic_string<lolita::character> file_path_;

    };

}


#endif /* E5351E45_F2AD_48E4_8AB5_ED530F0DBCDA */