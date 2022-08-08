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

    namespace detail
    {

        template<typename T>
        struct AggregateSizeTraits;

        template<typename... T>
        struct AggregateSizeTraits<Aggregate<T...>>
        {

            lolita::index value = sizeof...(T);

        };

        template<lolita::index I, typename T>
        struct AggregateElementTraits;

        template<lolita::index I, typename... T>
        struct AggregateElementTraits<I, Aggregate<T...>>
        {

            using type = std::tuple_element_t<I, std::tuple<T...>>;

        };

    }

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
    static constexpr
    lolita::index aggregate_size_v = detail::AggregateSizeTraits<T>::value;

    template<lolita::index I, typename... T>
    using aggregate_element_t = typename detail::AggregateElementTraits<I, T...>::type;

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

    namespace detail
    {

        template<typename... t_T, typename... t_U, lolita::index... t_i, lolita::index... t_j>
        static constexpr
        std::tuple<t_T..., t_U...>
        getMergedTuple(
                std::tuple<t_T...> const & first_tuple,
                std::tuple<t_U...> const & second_tuple,
                std::integer_sequence<lolita::index, t_i...>,
                std::integer_sequence<lolita::index, t_j...>
        )
        {
            return std::make_tuple(std::get<t_i>(first_tuple)..., std::get<t_j>(second_tuple)...);
        }

        template<lolita::index t_offset, typename... t_T, lolita::index... t_i>
        constexpr
        auto
        getSlicedTuple(
                std::tuple<t_T...> const & tuple,
                std::integer_sequence<lolita::index, t_i...>
        )
        {
            return std::make_tuple(std::get<t_i + t_offset>(tuple)...);
        }

    }

    template<typename... t_T, typename... t_U>
    static constexpr
    std::tuple<t_T..., t_U...>
    getMergedTuple(
            std::tuple<t_T...> const & first_tuple,
            std::tuple<t_U...> const & second_tuple
    )
    {
        auto const constexpr fi = std::make_integer_sequence<lolita::index, std::tuple_size_v<std::tuple<t_T...>>>{};
        auto const constexpr si = std::make_integer_sequence<lolita::index, std::tuple_size_v<std::tuple<t_U...>>>{};
        return lolita::utility::detail::getMergedTuple(first_tuple, second_tuple, fi, si);
    }

    template<lolita::index _begin, lolita::index _end, typename... t_T>
    constexpr
    auto
    getSlicedTuple(
            std::tuple<t_T...> const & tuple
    )
    requires(_end >= _begin && sizeof...(t_T) >= _end)
    {
        return detail::getSlicedTuple<_begin>(tuple, std::make_integer_sequence<lolita::index, _end - _begin>{});
    }

    // template<lolita::index _begin, lolita::index _end>
    // constexpr
    // auto
    // getSlicedTuple(
    //         std::tuple<> const & tuple
    // )
    // {
    //     return tuple;
    // }

    template<typename t_T, lolita::index _begin, lolita::index _end>
    requires(_end >= _begin && std::tuple_size_v<t_T> >= _end)
    using tuple_slice_t = decltype(lolita::utility::getSlicedTuple<_begin, _end>(std::declval<t_T>()));

    template<typename t_T, typename t_U>
    using tuple_merge_t = decltype(lolita::utility::getMergedTuple(std::declval<t_T>(), std::declval<t_U>()));

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

    template<typename t_T>
    struct tuple_singleton_traits;

    template<typename t_T, typename... t_U>
    struct tuple_singleton_traits<std::tuple<t_T, t_U...>>
    {

        using type = tuple_cat_t<std::tuple<t_T>, typename tuple_singleton_traits<std::tuple<t_U...>>::type>;

    };

    template<typename t_T, typename... t_U>
    requires((std::is_same_v<t_T, t_U> || ...))
    struct tuple_singleton_traits<std::tuple<t_T, t_U...>>
    {

        using type = typename tuple_singleton_traits<std::tuple<t_U...>>::type;

    };

    template<typename t_T>
    struct tuple_singleton_traits<std::tuple<t_T>>
    {

        using type = std::tuple<t_T>;

    };
    
    template<typename t_T>
    using tuple_sort_t2 = tuple_singleton_traits<t_T>::type;

    // ------------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------

    // template<lolita::integer t_i, lolita::integer t_a, lolita::integer t_b, typename t_T>
    // struct Filt3;

    // template<lolita::integer t_i, lolita::integer t_a, lolita::integer t_b, typename t_T, typename... t_U>
    // struct Filt3<t_i, t_a, t_b, std::tuple<t_T, t_U...>>
    // {

    //     using type = typename Filt3<t_i - 1, t_a, t_b, std::tuple<t_U...>>::type;

    // };

    // template<lolita::integer t_i, lolita::integer t_a, lolita::integer t_b, typename t_T, typename... t_U>
    // requires(t_a <= t_i <= t_b)
    // struct Filt3<t_i, t_a, t_b, std::tuple<t_T, t_U...>>
    // {

    //     using type = tuple_cat_t<std::tuple<t_T>, typename Filt3<t_i - 1, t_a, t_b, std::tuple<t_U...>>::type>;

    // };

    // template<lolita::integer t_i, lolita::integer t_a, lolita::integer t_b, typename t_T>
    // struct Filt3<t_i, t_a, t_b, std::tuple<t_T>>
    // {

    //     using type = std::tuple<t_T>;

    // };
    
    // template<lolita::integer t_a, lolita::integer t_b, typename t_T>
    // using tuple_slice_t2 = Filt3<std::tuple_size_v<t_T>, t_a, t_b, t_T>::type;

    // ------------------------------------------------------------------------------------------------------

    template<lolita::integer N, typename Tuple_Type>
    struct tuple_trunc_head;
    // {};

    template<lolita::integer N, typename Head, typename... Tail>
    struct tuple_trunc_head<N, std::tuple<Head, Tail...>>
    {
        using type = typename tuple_trunc_head<N - 1, std::tuple<Tail...>>::type;
    };

    template<typename Head, typename... Tail>
    struct tuple_trunc_head<0, std::tuple<Head, Tail...>>
    {
        using type = std::tuple<Head, Tail...>;
    };

    // ------------------------------------------------------------------------------------------------------

    static void inline
    removeCharacter(
            std::basic_string<lolita::character> & line,
            lolita::character character
    )
    {
        line.erase(std::remove(line.begin(), line.end(), character), line.end());
    }

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
            if (!file) {
                throw std::runtime_error("Could not open file");
            }
            for (std::basic_string<lolita::character> line; std::getline(file, line); ) {
                lines_.push_back(line);
            }
        }

    public:

        std::vector<std::basic_string<lolita::character>> lines_;

        std::basic_string<lolita::character> file_path_;

    };

}


#endif /* E5351E45_F2AD_48E4_8AB5_ED530F0DBCDA */
