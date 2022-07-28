//
// Created by dsiedel on 06/05/22.
//

#ifndef LOLITA_LOLITA_UTILITY_HXX
#define LOLITA_LOLITA_UTILITY_HXX

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

    template<typename... T>
    struct Aggregate : detail::AggregateHead<0, T...>
    {

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

        // constexpr
        // std::tuple<T...>
        // asTuple()
        // const
        // {
        //     auto tuple = std::tuple<T...>();
        //     auto fill = [&] <lolita::integer t_i = 0> (
        //         auto & self
        //     )
        //     constexpr mutable
        //     {
        //         std::get<t_i>(tuple) = detail::get<t_i>(* this);
        //         if constexpr (t_i < sizeof...(T) - 1)
        //         {
        //             self.template operator ()<t_i + 1>(self);
        //         }
        //     };
        //     fill(fill);
        //     return tuple;
        // }

    };

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

//     namespace detail
//     {

//         struct _Void {};

//     }

//     enum struct Output
//     {

//         Success,
//         Failure

//     };

//     template<typename _T>
//     struct OutputA
//     {

//         using ReturnType = _T;

//         lolita::utility::Output result_;

//         _T value_;

//     };

//     struct EnumA
//     {

//     private:

//         using Label = std::array<lolita::character, 50>;

//         template<typename... _U>
//         static constexpr
//         Label
//         makeLabel(
//                 std::basic_string_view<_U> const &... str
//         )
//         requires(std::same_as<lolita::character, _U> && ...)
//         {
//             auto label = Label();
//             auto count = lolita::index(0);
//             auto make = [&] (auto const & s) constexpr mutable {
//                 for (auto i = 0; i < s.size(); ++i) {
//                     label[i + count] = s[i];
//                 }
//                 count += s.size();
//             };
//             (make(str), ...);
//             return label;
//         }

//         template<typename... _U>
//         static constexpr
//         Label
//         makeLabel(
//                 std::basic_string_view<_U> &&... str
//         )
//         requires(std::same_as<lolita::character, _U> && ...)
//         {
//             auto label = Label();
//             auto count = lolita::index(0);
//             auto make = [&] (auto && s) constexpr mutable {
//                 for (auto i = 0; i < s.size(); ++i) {
//                     label[i + count] = s[i];
//                 }
//                 count += s.size();
//             };
//             (make(std::forward<std::basic_string_view<_U>>(str)), ...);
//             return label;
//         }

//     public:

//         constexpr
//         EnumA(
//                 std::basic_string_view<lolita::character> const & str
//         )
//         :
//         tag_(makeLabel(str))
//         {}

//         constexpr
//         EnumA(
//                 std::basic_string_view<lolita::character> && str
//         )
//         :
//         tag_(makeLabel(std::forward<std::basic_string_view<lolita::character>>(str)))
//         {}

// //        template<typename... _U>
// //        constexpr
// //        EnumA(
// //                std::basic_string_view<_U> const &... str
// //        )
// //        :
// //        tag_(makeLabel(str...))
// //        {}
// //
// //
// //        template<typename... _U>
// //        constexpr
// //        EnumA(
// //                std::basic_string_view<_U> &&... str
// //        )
// //        :
// //        tag_(makeLabel(std::forward<std::basic_string_view<lolita::character>>(str)...))
// //        {}

//         constexpr
//         lolita::boolean
//         operator==(
//                 EnumA const & other
//         )
//         const = default;

//         constexpr
//         lolita::boolean
//         operator!=(
//                 EnumA const & other
//         )
//         const = default;

//         constexpr
//         lolita::boolean
//         operator==(
//                 std::basic_string_view<lolita::character> const & other
//         )
//         const
//         {
//             auto view = std::basic_string_view<lolita::character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), lolita::character())));
//             return other == view;
//         }

//         constexpr
//         lolita::boolean
//         operator!=(
//                 std::basic_string_view<lolita::character> const & other
//         )
//         const
//         {
//             return !(* this == other);
//         }

//         friend
//         std::ostream &
//         operator<<(
//                 std::ostream & os,
//                 EnumA const & enuma
//         )
//         {
//             auto const & tag = enuma.tag_;
//             os << std::basic_string_view<lolita::character>(tag.data(), std::distance(tag.begin(), std::find(tag.begin(), tag.end(), lolita::character())));
//             return os;
//         }

//         Label tag_;

//     };

//     using Tag = std::array<lolita::character, 50>;

//     struct Labell
//     {

//     private:

//         static constexpr
//         Tag
//         setTag(
//                 std::basic_string_view<lolita::character> str
//         )
//         {
//             auto tag = Tag();
//             std::copy_n(std::make_move_iterator(str.begin()), str.size(), tag.begin());
// //            auto count = lolita::index(0);
// //            for (auto c : std::forward<std::basic_string_view<lolita::character>>(str)) {
// //                tag[count] = c;
// //                count ++;
// //            }
//             return tag;
//         }

//     public:

//         constexpr
//         Labell(
//                 std::basic_string_view<lolita::character> str
//         )
//         :
//         tag_(setTag(str))
//         {}

//         constexpr
//         lolita::boolean
//         operator==(
//                 Labell const & other
//         )
//         const = default;

//         constexpr
//         lolita::boolean
//         operator!=(
//                 Labell const & other
//         )
//         const = default;

//         constexpr
//         lolita::boolean
//         operator==(
//                 std::basic_string_view<lolita::character> const & tag
//         )
//         const
//         {
//             return tag == this->view();
//         }

//         constexpr
//         lolita::boolean
//         operator!=(
//                 std::basic_string_view<lolita::character> const & other
//         )
//         const
//         {
//             return !(* this == other);
//         }

//         constexpr
//         std::basic_string_view<lolita::character>
//         view()
//         const
//         {
//             return std::basic_string_view<lolita::character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), lolita::character())));
//         }

//         friend
//         std::ostream &
//         operator<<(
//                 std::ostream & os,
//                 Labell const & label
//         )
//         {
//             os << label.view();
//             return os;
//         }

//         Tag tag_;

//     };

//     template<typename t_Base>
//     struct Enumeration
//     {

//     private:

//         template<typename... _U>
//         static constexpr
//         Tag
//         setTag(
//                 std::basic_string_view<_U> const &... str
//         )
//         requires(std::same_as<lolita::character, _U> && ...)
//         {
//             auto label = Tag();
//             auto count = lolita::index(0);
//             auto make = [&] (auto const & s) constexpr mutable {
//                 for (auto i = 0; i < s.size(); ++i) {
//                     label[i + count] = s[i];
//                 }
//                 count += s.size();
//             };
//             (make(str), ...);
//             return label;
//         }

//         template<typename... _U>
//         static constexpr
//         Tag
//         setTag(
//                 std::basic_string_view<_U> &&... str
//         )
//         requires(std::same_as<lolita::character, _U> && ...)
//         {
//             auto label = Tag();
//             auto count = lolita::index(0);
//             auto make = [&] (auto && s) constexpr mutable {
//                 for (auto i = 0; i < s.size(); ++i) {
//                     label[i + count] = s[i];
//                 }
//                 count += s.size();
//             };
//             (make(std::forward<std::basic_string_view<_U>>(str)), ...);
//             return label;
//         }

//     public:

//         constexpr
//         Enumeration(
//                 std::basic_string_view<lolita::character> const & str
//         )
//         :
//         tag_(setTag(str))
//         {}

//         constexpr
//         Enumeration(
//                 std::basic_string_view<lolita::character> && str
//         )
//         :
//         tag_(setTag(std::forward<std::basic_string_view<lolita::character>>(str)))
//         {}

//         constexpr
//         std::basic_string_view<lolita::character>
//         view()
//         const
//         {
//             return std::basic_string_view<lolita::character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), lolita::character())));
//         }

// //        template<typename... _U>
// //        constexpr
// //        EnumA(
// //                std::basic_string_view<_U> const &... str
// //        )
// //        :
// //        tag_(makeLabel(str...))
// //        {}
// //
// //
// //        template<typename... _U>
// //        constexpr
// //        EnumA(
// //                std::basic_string_view<_U> &&... str
// //        )
// //        :
// //        tag_(makeLabel(std::forward<std::basic_string_view<lolita::character>>(str)...))
// //        {}

//         constexpr
//         lolita::boolean
//         operator==(
//                 Enumeration const & other
//         )
//         const = default;

//         constexpr
//         lolita::boolean
//         operator!=(
//                 Enumeration const & other
//         )
//         const = default;

//         constexpr
//         lolita::boolean
//         operator==(
//                 std::basic_string_view<lolita::character> const & other
//         )
//         const
//         {
//             return other == this->view();
//         }

//         constexpr
//         lolita::boolean
//         operator!=(
//                 std::basic_string_view<lolita::character> const & other
//         )
//         const
//         {
//             return !(* this == other);
//         }

//         constexpr
//         lolita::boolean
//         operator==(
//                 lolita::utility::Tag const & tag
//         )
//         const
//         {
//             auto constexpr t_null = lolita::character();
//             auto view = std::basic_string_view<lolita::character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), t_null)));
//             auto view1 = std::basic_string_view<lolita::character>(tag.data(), std::distance(tag.begin(), std::find(tag.begin(), tag.end(), t_null)));
//             return view1 == view;
//         }

//         constexpr
//         lolita::boolean
//         operator!=(
//                 lolita::utility::Tag const & tag
//         )
//         const
//         {
//             return !(* this == tag);
//         }

//         friend
//         std::ostream &
//         operator<<(
//                 std::ostream & os,
//                 Enumeration const & enuma
//         )
//         {
//             auto constexpr t_null = lolita::character();
//             auto const & tag = enuma.tag_;
//             os << std::basic_string_view<lolita::character>(tag.data(), std::distance(tag.begin(), std::find(tag.begin(), tag.end(), t_null)));
//             return os;
//         }

//         Tag tag_;

//     };

//     struct TupleSlice
//     {

//         struct detail
//         {

//             template<lolita::index _offset, typename... _T, lolita::index... _i>
//             static constexpr
//             auto
//             tupleSlice(
//                     std::tuple<_T...> const & tuple,
//                     std::integer_sequence<lolita::index, _i...>
//             )
//             {
//                 return std::make_tuple(std::get<_i + _offset>(tuple)...);
//             }

//         };

//         template<lolita::index _begin, lolita::index _end, typename... _T>
//         static constexpr
//         auto
//         tupleSlice(
//                 std::tuple<_T...> const & tuple
//         )
//         requires(_end >= _begin && sizeof...(_T) >= _end)
//         {
//             return detail::template tupleSlice<_begin, _T...>(tuple, std::make_integer_sequence<lolita::index, _end - _begin>{});
//         }

//     };

//     using Label = std::array<lolita::character, 20>;

//     template<typename... _U>
//     static constexpr
//     lolita::utility::Label
//     makeLabel(
//             std::basic_string_view<_U> const &... str
//     )
//     requires(std::same_as<lolita::character, _U> && ...)
//     {
//         auto label = lolita::utility::Label();
//         auto count = lolita::index(0);
//         auto make = [&] (auto const & s) constexpr mutable {
//             for (auto i = 0; i < s.size(); ++i) {
//                 label[i + count] = s[i];
//             }
//             count += s.size();
//         };
//         (make(str), ...);
//         return label;
//     }

//     template<typename... _U>
//     static constexpr
//     lolita::utility::Label
//     makeLabel(
//             std::basic_string_view<_U> &&... str
//     )
//     requires(std::same_as<lolita::character, _U> && ...)
//     {
//         auto label = lolita::utility::Label();
//         auto count = lolita::index(0);
//         auto make = [&] (auto && s) constexpr mutable {
//             for (auto i = 0; i < s.size(); ++i) {
//                 label[i + count] = s[i];
//             }
//             count += s.size();
//         };
//         (make(std::forward<std::basic_string_view<_U>>(str)), ...);
//         return label;
//     }

//     static constexpr
//     std::basic_string_view<lolita::character>
//     readLabel(
//             lolita::utility::Label const & label
//     )
//     {
//         return std::basic_string_view<lolita::character>(label.data(), std::distance(label.begin(), std::find(label.begin(), label.end(), lolita::character())));
//     }

//    namespace detail
//    {
//
//        template<lolita::index _offset, typename... _T, lolita::index... _i>
//        constexpr
//        auto
//        tupleSlice(
//                std::tuple<_T...> const & tuple,
//                std::integer_sequence<lolita::index, _i...>
//        )
//        {
//            return std::make_tuple(std::get<_i + _offset>(tuple)...);
//        }
//
//    }
//
//    template<lolita::index _begin, lolita::index _end, typename... _T>
//    constexpr
//    auto
//    tupleSlice(
//            std::tuple<_T...> const & tuple
//    )
//    requires(_end >= _begin && sizeof...(_T) >= _end)
//    {
//        return detail::tupleSlice<_begin>(tuple, std::make_integer_sequence<lolita::index, _end - _begin>{});
//    }

    namespace detail
    {

        template<typename... t_T, typename... t_U, lolita::index... t_i, lolita::index... t_j>
        static constexpr
        std::tuple<t_T..., t_U...>
        tupleMerge(
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
        tupleSlice(
                std::tuple<t_T...> const & tuple,
                std::integer_sequence<lolita::index, t_i...>
        )
        {
            return std::make_tuple(std::get<t_i + t_offset>(tuple)...);
        }

    }

    template<typename... _T, typename... _U>
    static constexpr
    std::tuple<_T..., _U...>
    tupleMerge(
            std::tuple<_T...> const & first_tuple,
            std::tuple<_U...> const & second_tuple
    )
    {
        auto const constexpr fi = std::make_integer_sequence<lolita::index, std::tuple_size_v<std::tuple<_T...>>>{};
        auto const constexpr si = std::make_integer_sequence<lolita::index, std::tuple_size_v<std::tuple<_U...>>>{};
        return lolita::utility::detail::tupleMerge(first_tuple, second_tuple, fi, si);
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

    template<typename _T, lolita::index _begin, lolita::index _end>
            requires(_end >= _begin && std::tuple_size_v<_T> >= _end)
    using tuple_slice_t = decltype(lolita::utility::tupleSlice<_begin, _end>(std::declval<_T>()));

    template<typename _T, typename _U>
    using tuple_merge_t = decltype(lolita::utility::tupleMerge(std::declval<_T>(), std::declval<_U>()));

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

    /*
     *
     */

    

    // template<typename... T>
    // using Aggregate = detail::AggregateHead<0u, T...>;

    // namespace detail
    // {

    //     template<typename T>
    //     struct AggregateSizeTraits;

    //     template<typename... T>
    //     struct AggregateSizeTraits<Aggregate<T...>>
    //     {

    //         lolita::index value = sizeof...(T);

    //     };

    //     template<lolita::index I, typename T>
    //     struct AggregateElementTraits;

    //     template<lolita::index I, typename... T>
    //     struct AggregateElementTraits<I, Aggregate<T...>>
    //     {

    //         using type = std::tuple_element_t<I, std::tuple<T...>>;

    //     };

    // }

    // template<typename T>
    // static constexpr
    // lolita::index aggregate_size_v = detail::AggregateSizeTraits<T>::value;

    // template<lolita::index I, typename... T>
    // using aggregate_element_t = typename detail::AggregateElementTraits<I, T...>::type;

    // template<typename... T>
    // static constexpr
    // lolita::index
    // size(
    //         Aggregate<T...> const &
    //         aggregate
    // )
    // {
    //     return sizeof...(T);
    // }

    // template<lolita::index I, typename... T>
    // static constexpr
    // std::tuple_element_t<I, std::tuple<T...>> const &
    // get(
    //         Aggregate<T...> const &
    //         aggregate
    // )
    // {
    //     return detail::get<I>(aggregate);
    // }

    // template<lolita::index I, typename... T>
    // static constexpr
    // std::tuple_element_t<I, std::tuple<T...>> &
    // get(
    //         Aggregate<T...> &
    //         aggregate
    // )
    // {
    //     return detail::get<I>(aggregate);
    // }

    // template<auto A>
    // struct S{};



    // template<typename T, lolita::index... I>
    // using TupleParts = std::tuple<std::tuple_element_t<I, T>...>;

    // template<auto... A>
    // struct Expand
    // {

    //     using E = std::tuple<S<A>...>;

    // };



    // template<lolita::index... s>
    // struct seq
    // {

    //     using type = seq<s...>;

    // };

    // template<lolita::index max, lolita::index... s>
    // struct make_seq : make_seq<max - 1, max - 1, s...> {};

    // template<lolita::index...s>
    // struct make_seq<0, s...> : seq<s...> {};

    // template<lolita::index... s, typename Tuple>
    // auto
    // extract_tuple(
    //         seq<s...>,
    //         Tuple &
    //         tup
    // )
    // {
    //     return std::make_tuple(std::get<s>(tup)...);
    // }

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
