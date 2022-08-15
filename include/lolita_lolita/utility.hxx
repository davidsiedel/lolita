#ifndef BD1D6E93_67D6_418F_A38E_F9C11C70DB4B
#define BD1D6E93_67D6_418F_A38E_F9C11C70DB4B

#include <iostream>
#include <fstream>
#include <execution>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <ostream>
#include <iomanip>

#include "lolita_lolita/config.hxx"

namespace lolita::utility
{

    template<auto t_value>
    struct Holder
    {

        auto static constexpr value_ = t_value;

    };

    template<typename t_T, typename t_U>
    static constexpr
    Boolean
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

    template<typename t_T, Integer t_a, Integer t_b>
    struct TupleSliceTraits;

    template<typename... t_T, Integer t_a, Integer t_b>
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = tuple_cat_t<std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>, typename TupleSliceTraits<std::tuple<t_T...>, t_a + 1, t_b>::type>;

    };
    
    template<typename... t_T, Integer t_a, Integer t_b>
    requires(t_a == t_b - 1)
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = std::tuple<std::tuple_element_t<t_b - 1, std::tuple<t_T...>>>;

    };
    
    template<typename... t_T, Integer t_a, Integer t_b>
    requires(t_a == t_b)
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = std::tuple<>;

    };

    template<typename t_T, Integer t_a, Integer t_b>
    using tuple_slice_t = TupleSliceTraits<t_T, t_a, t_b>::type;

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, typename t_U>
    struct TupleHasTraits;

    template<typename... t_T, typename t_U>
    struct TupleHasTraits<std::tuple<t_T...>, t_U>
    {

        Boolean static constexpr value = (std::is_same_v<t_U, t_T> || ...);

    };

    template<typename t_T, typename t_U>
    static constexpr
    Boolean tuple_has_v = TupleHasTraits<t_T, t_U>::value;

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, Integer t_a>
    struct TupleUniqueTraits;

    template<typename... t_T, Integer t_a>
    requires(t_a < sizeof...(t_T) - 1 && !tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = tuple_cat_t<std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>, typename TupleUniqueTraits<std::tuple<t_T...>, t_a + 1>::type>;

    };

    template<typename... t_T, Integer t_a>
    requires(t_a < sizeof...(t_T) - 1 && tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = tuple_cat_t<std::tuple<>, typename TupleUniqueTraits<std::tuple<t_T...>, t_a + 1>::type>;

    };

    template<typename... t_T, Integer t_a>
    requires(t_a == sizeof...(t_T) - 1 && !tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>;

    };

    template<typename... t_T, Integer t_a>
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

    template<Integer t_i, typename t_T>
    struct AggregateTail
    {

        constexpr
        AggregateTail()
        :
        value_(t_T())
        {}

        constexpr explicit
        AggregateTail(
            t_T const & value
        )
        :
        value_(value)
        {}

        constexpr explicit
        AggregateTail(
            t_T && value
        )
        :
        value_(std::forward<t_T>(value))
        {}

        constexpr
        Boolean
        operator==(
            AggregateTail const &
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            AggregateTail const &
        )
        const = default;

        t_T value_;

    };

    template<Integer t_i, typename... t_T>
    struct AggregateHead;

    template<Integer t_i>
    struct AggregateHead<t_i>
    {

        constexpr
        AggregateHead()
        {}

        constexpr
        Boolean
        operator==(
            AggregateHead const &
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            AggregateHead const &
        )
        const = default;

    };

    template<Integer t_i, typename t_T, typename... t_U>
    struct AggregateHead<t_i, t_T, t_U...> : AggregateTail<t_i, t_T>, AggregateHead<t_i + 1, t_U...>
    {

        constexpr
        AggregateHead(
            t_T const & value,
            t_U const &... values
        )
        :
        AggregateTail<t_i, t_T>(value),
        AggregateHead<t_i + 1, t_U...>(values...)
        {}

        constexpr
        AggregateHead(
            t_T && value,
            t_U &&... values
        )
        :
        AggregateTail<t_i, t_T>(std::forward<t_T>(value)),
        AggregateHead<t_i + 1, t_U...>(std::forward<t_U>(values)...)
        {}

        constexpr
        Boolean
        operator==(
            AggregateHead const &
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            AggregateHead const &
        )
        const = default;

    };

    template<Integer t_i, typename t_T, typename... t_U>
    static constexpr
    t_T &
    get(
        AggregateHead<t_i, t_T, t_U...> & aggregate
    )
    {
        return aggregate.AggregateTail<t_i, t_T>::value_;
    }

    template<Integer t_i, typename t_T, typename... t_U>
    static constexpr
    t_T const &
    get(
        AggregateHead<t_i, t_T, t_U...> const & aggregate
    )
    {
        return aggregate.AggregateTail<t_i, t_T>::value_;
    }

    template<typename... T>
    using Aggregate = AggregateHead<0, T...>;

    // ------------------------------------------------------------------------------------------------------

    template<typename T>
    struct AggregateSizeTraits;

    template<typename... T>
    struct AggregateSizeTraits<Aggregate<T...>>
    {

        Integer static constexpr value = sizeof...(T);

    };

    template<typename T>
    static constexpr
    Integer aggregate_size_v = AggregateSizeTraits<T>::value;

    template<Integer I, typename T>
    struct AggregateElementTraits;

    template<Integer I, typename... T>
    struct AggregateElementTraits<I, Aggregate<T...>>
    {

        using type = std::tuple_element_t<I, std::tuple<T...>>;

    };

    template<Integer I, typename... T>
    using aggregate_element_t = typename AggregateElementTraits<I, T...>::type;

    template<template<auto> typename t_T, auto t_aggregate, Integer t_i>
    struct aggregate_expansion_traits
    {

        using type = tuple_cat_t<std::tuple<t_T<get<t_i>(t_aggregate)>>, typename aggregate_expansion_traits<t_T, t_aggregate, t_i + 1>::type>;

    };
    
    template<template<auto> typename t_T, auto t_aggregate, Integer t_i>
    requires(t_i == aggregate_size_v<std::decay_t<decltype(t_aggregate)>> - 1)
    struct aggregate_expansion_traits<t_T, t_aggregate, t_i>
    {

        using type = std::tuple<t_T<get<t_i>(t_aggregate)>>;

    };

    template<template<auto> typename t_T, auto t_aggregate>
    using aggregate_expansion_t = aggregate_expansion_traits<t_T, t_aggregate, 0>::type;

    // ------------------------------------------------------------------------------------------------------

    // struct Label
    // {

    //     using Tag = std::array<Character, 50>;

    // private:

    //     static constexpr
    //     Label::Tag
    //     setTag(
    //         std::basic_string_view<Character> str
    //     )
    //     {
    //         auto tag = Label::Tag();
    //         auto count = Integer(0);
    //         for (auto c : str) {
    //             tag[count] = c;
    //             count ++;
    //         }
    //         return tag;
    //     }

    // public:

    //     constexpr
    //     Label(
    //         std::basic_string_view<Character> str
    //     )
    //     :
    //     tag_(setTag(str))
    //     {}

    //     constexpr
    //     Boolean
    //     operator==(
    //         Label const & other
    //     )
    //     const = default;

    //     constexpr
    //     Boolean
    //     operator!=(
    //         Label const & other
    //     )
    //     const = default;

    //     constexpr inline
    //     Boolean
    //     operator==(
    //         std::basic_string_view<Character> str
    //     )
    //     const
    //     {
    //         return str == this->view();
    //     }

    //     constexpr inline
    //     Boolean
    //     operator!=(
    //         std::basic_string_view<Character> str
    //     )
    //     const
    //     {
    //         return !(* this == str);
    //     }

    //     constexpr inline
    //     std::basic_string_view<Character>
    //     view()
    //     const
    //     {
    //         return std::basic_string_view<Character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), Character())));
    //     }

    //     friend inline
    //     std::ostream &
    //     operator<<(
    //         std::ostream & os,
    //         Label const & label
    //     )
    //     {
    //         os << label.view();
    //         return os;
    //     }

    //     Label::Tag tag_;

    // };

    // ------------------------------------------------------------------------------------------------------

    struct File
    {

        File() = default;

        explicit
        File(
            std::basic_string_view<Character> file_path
        )
        :
        file_path_(file_path)
        {
            readLines();
        }

        Boolean
        operator==(
            File const &
            other
        )
        const = default;

        Boolean
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
            std::basic_ifstream<Character> file(file_path_);
            if (!file)
            {
                throw std::runtime_error("Could not open file");
            }
            for (std::basic_string<Character> line; std::getline(file, line); )
            {
                lines_.push_back(line);
            }
        }

    public:

        std::vector<std::basic_string<Character>> lines_;

        std::basic_string<Character> file_path_;

    };

} // namespace lolita::utility


#endif /* BD1D6E93_67D6_418F_A38E_F9C11C70DB4B */
