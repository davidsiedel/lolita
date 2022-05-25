//
// Created by dsiedel on 06/05/22.
//

#ifndef LOLITA_LOLITA_UTILITY_HXX
#define LOLITA_LOLITA_UTILITY_HXX

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include "lolita/lolita.hxx"
#include "lolita/lolita_algebra.hxx"

namespace lolita::utility
{

    enum struct Output
    {

        Success,
        Failure

    };

    template<typename _T>
    struct OutputA
    {

        using ReturnType = _T;

        lolita::utility::Output result_;

        _T value_;

    };

    using Label = std::array<lolita::character, 150>;

    template<typename... _U>
    static constexpr
    lolita::utility::Label
    label(
            std::basic_string_view<_U> &&... str
    )
    requires(std::same_as<lolita::character, _U> && ...)
    {
        auto label = lolita::utility::Label();
        auto count = lolita::index(0);
        auto make = [&] (auto && s) constexpr mutable {
            for (auto i = 0; i < s.size(); ++i) {
                label[i + count] = s[i];
            }
            count += s.size();
        };
        (make(std::forward<std::basic_string_view<_U>>(str)), ...);
        for (auto i = count; i < label.size(); ++i) {
            label[i] = '#';
        }
        return label;
    }

    static constexpr
    std::basic_string_view<lolita::character>
    readLabel(
            lolita::utility::Label const & label
    )
    {
        return std::basic_string_view<lolita::character>(label.data(), std::distance(label.begin(), std::find(label.begin(), label.end(), '#')));
    }

    namespace detail
    {

        template<lolita::index _offset, typename... _T, lolita::index... _i>
        constexpr
        auto
        tupleSlice(
                std::tuple<_T...> const & tuple,
                std::integer_sequence<lolita::index, _i...>
        )
        {
            return std::make_tuple(std::get<_i + _offset>(tuple)...);
        }

    }

    template<lolita::index _begin, lolita::index _end, typename... _T>
    constexpr
    auto
    tupleSlice(
            std::tuple<_T...> const & tuple
    )
    requires(_end >= _begin && sizeof...(_T) >= _end)
    {
        return detail::tupleSlice<_begin>(tuple, std::make_integer_sequence<lolita::index, _end - _begin>{});
    }

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
                std::basic_string_view<lolita::character> const & file_path_arg
        )
        :
        file_path_(file_path_arg)
        {
            readLines();
        }

        explicit
        File(
                std::basic_string_view<lolita::character> && file_path_arg
        )
        :
        file_path_(std::forward<std::basic_string_view<lolita::character>>(file_path_arg))
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

    /*
     *
     */

    namespace detail
    {

        template<lolita::index I, typename T>
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

            T value_;

        };

        template<lolita::index I, typename... T>
        struct AggregateHead;

        template<lolita::index I>
        struct AggregateHead<I>
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

        template<lolita::index I, typename T, typename... U>
        struct AggregateHead<I, T, U...> : public AggregateTail<I, T>, public AggregateHead<I + 1, U...>
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

        template<lolita::index I, typename T, typename... U>
        static constexpr
        T &
        get(
                AggregateHead<I, T, U...> &
                aggregate_arg
        )
        {
            return aggregate_arg.AggregateTail<I, T>::value_;
        }

        template<lolita::index I, typename T, typename... U>
        static constexpr
        T const &
        get(
                AggregateHead<I, T, U...> const &
                aggregate_arg
        )
        {
            return aggregate_arg.AggregateTail<I, T>::value_;
        }

    }

    template<typename... T>
    using Aggregate = detail::AggregateHead<0u, T...>;

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

    template<typename T>
    static constexpr
    lolita::index aggregate_size_v = detail::AggregateSizeTraits<T>::value;

    template<lolita::index I, typename... T>
    using aggregate_element_t = typename detail::AggregateElementTraits<I, T...>::type;

    template<typename... T>
    static constexpr
    lolita::index
    size(
            Aggregate<T...> const &
            aggregate
    )
    {
        return sizeof...(T);
    }

    template<lolita::index I, typename... T>
    static constexpr
    std::tuple_element_t<I, std::tuple<T...>> const &
    get(
            Aggregate<T...> const &
            aggregate
    )
    {
        return detail::get<I>(aggregate);
    }

    template<lolita::index I, typename... T>
    static constexpr
    std::tuple_element_t<I, std::tuple<T...>> &
    get(
            Aggregate<T...> &
            aggregate
    )
    {
        return detail::get<I>(aggregate);
    }

    template<auto A>
    struct S{};



    template<typename T, lolita::index... I>
    using TupleParts = std::tuple<std::tuple_element_t<I, T>...>;

    template<auto... A>
    struct Expand
    {

        using E = std::tuple<S<A>...>;

    };



    template<lolita::index... s>
    struct seq
    {

        using type = seq<s...>;

    };

    template<lolita::index max, lolita::index... s>
    struct make_seq : make_seq<max - 1, max - 1, s...> {};

    template<lolita::index...s>
    struct make_seq<0, s...> : seq<s...> {};

    template<lolita::index... s, typename Tuple>
    auto
    extract_tuple(
            seq<s...>,
            Tuple &
            tup
    )
    {
        return std::make_tuple(std::get<s>(tup)...);
    }

//        template<typename... T>
//        struct TupleExpansion
//        {
//
//            TupleExpansion(
//                    std::tuple<T...> const & t
//            )
//            :
//            tpl(t)
//            {}
//
//            std::tuple<T...> const & tpl;
//
////            static constexpr
////            auto subtuple_(const std::tuple<T...>& t, std::index_sequence<I...>) {
////                return std::make_tuple(std::get<I>(t)...);
////            }
//            template<lolita::index... K>
//            constexpr
//            auto
//            subtuple_2(std::index_sequence<K...>)
//            {
//                return std::make_tuple(S<std::get<K>(tpl)>()...);
//            }
//
////            template <int Trim = 0>
////            static constexpr
////            auto subtuple(const std::tuple<T...>& t) {
////                return subtuple_(t, std::make_index_sequence<sizeof...(T) - Trim>());
////            }
//
//            constexpr
//            auto
//            subtuple_n()
//            {
//                return subtuple_2(std::make_index_sequence<sizeof...(T)>());
//            }
//
//        };
//
//        template<typename... T, lolita::index... K>
//        constexpr
//        auto
//        subtuple_2(std::tuple<T...> const & tpl, std::index_sequence<K...>)
//        {
//            return std::make_tuple(S<std::get<K>(tpl)>()...);
//        }
//
//        constexpr
//        auto
//        subtuple_n()
//        {
//            return subtuple_2(std::make_index_sequence<sizeof...(T)>());
//        }

//        template<lolita::index I, typename... T>
//        static constexpr
//        auto
//        expandTuple(
//                std::tuple<T...> const &
//                tuple
//        )
//        {
//
//            auto doit = [&] <lolita::index... K> (std::index_sequence<K...> seq) constexpr mutable {
//                return std::make_tuple(S<std::get<K>(tuple)>()...);
//            };
////            auto doit2 = [&] (auto kkk) constexpr mutable {
////                return std::make_tuple(S<kkk>()...);
////            };
//            return doit.template(std::make_index_sequence<I>());
//        }
//
//        template<typename... T, lolita::index... I>
//        static constexpr
//        auto
//        expandD(
//                std::tuple<T...> const &
//                tuple,
//                std::index_sequence<I...>
//                sequence
//        )
//        {
//            return std::make_tuple(S<std::get<I>(tuple)>()...);
//        }

//        template <typename... T, lolita::index... I>
//        static constexpr
//        expandD(
//                std::tuple<T...> const &
//                tuple,
//                std::index_sequence<I...>
//                sequence
//        )
//        {
//            return std::make_tuple(S<std::get<I>(t)>()...);
//        }

//        template <typename... T, std::size_t... I>
//        static constexpr
//        auto subtuple_(const std::tuple<T...>& t, std::index_sequence<I...>) {
//            return std::make_tuple(std::get<I>(t)...);
//        }
//
//        template <typename... T, std::size_t... I>
//        static constexpr
//        auto subtuple__(const std::tuple<T...>& t, std::index_sequence<I...>) {
//            return std::make_tuple(S<std::get<I>(t)>()...);
//        }
//
//        template <int Trim = 0, typename... T>
//        static constexpr
//        auto subtuple(const std::tuple<T...>& t) {
//            return subtuple_(t, std::make_index_sequence<sizeof...(T) - Trim>());
//        }
//
//        template<typename... T>
//        static constexpr
//        auto subtuple2(const std::tuple<T...>& t) {
//            return subtuple__(t, std::make_index_sequence<sizeof...(T)>());
//        }

//        template<auto... I>
//        struct Mult {
//                std::array<lolita::index, sizeof...(I)> const static constexpr arr = {static_cast<lolita::index>(I)...};
//        };

//        template<typename Tuple, lolita::index... I>
//        using SubTuple = std::tuple<std::tuple_element_t<I, Tuple>...>;
//
//        template<typename Tuple, lolita::index... I>
//        struct Strip
//        {
//
//            constexpr
//            Strip(Tuple t, std::index_sequence<I...> s)
//            :
//            h(std::make_tuple(std::get<I>(t)...))
//            {}
//
//
//            using HH = std::tuple<std::tuple_element_t<I, Tuple>...>;
//
//            HH h;
//
//        };


//        template<lolita::index I, typename... T>
//        using TupleStrip = SubTuple<std::tuple<T...>, std::index_sequence<I>>;

//        template<typename... T>
//        auto
//        extract_tuple(
//                std::tuple<T...> &
//                tup
//        )
//        {
//            return std::make_tuple(std::get<std::make_index_sequence<2>>(tup)...);
//        }

//        template<typename T>
//        static constexpr
//        auto
//        extract_tuple2(
//                T &
//                tup
//        )
//        {
//            return std::make_tuple(std::get<s>(tup)...);
//        }

//        template<typename... T>
//        struct Expand
//        {
//
//            constexpr
//            Expand(T const &... a) : values({a...}) {}
//
//
//
//            std::tuple<T...> values;
//
//            std::tuple<S<T{}>...> output;
//
//        };


//        template<typename... T>
//        struct Aggregate
//        {
//
//            using Data = detail::AggregateHead<0, T...>;
//
//            template<lolita::numerics::NaturalType auto I>
//            using Type = typename std::tuple_element<I, std::tuple<T...>>::type;
//
//            static constexpr
//            lolita::index
//            size()
//            {
//                return sizeof...(T);
//            }
//
//            constexpr
//            lolita::boolean
//            operator==(
//                    Aggregate const &
//                    other
//            )
//            const
//            {
//                return data == other.data;
//            }
////            const = default;
//
//            constexpr
//            lolita::boolean
//            operator!=(
//                    Aggregate const &
//                    other
//            )
//            const
//            {
//                return !(* this == other);
//            }
////            const = default;
//
//            template<lolita::numerics::NaturalType auto I>
//            constexpr
//            auto const &
//            get()
//            const
//            {
//                return detail::get<I>(data);
//            }
//
//            template<lolita::numerics::NaturalType auto I>
//            constexpr
//            auto &
//            get()
//            {
//                return detail::get<I>(data);
//            }
//
//            Data data;
//
//        };


}

#endif //LOLITA_LOLITA_UTILITY_HXX
