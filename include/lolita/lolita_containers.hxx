//
// Created by dsiedel on 18/04/22.
//

#ifndef LOLITA_LOLITA_CONTAINERS_HXX
#define LOLITA_LOLITA_CONTAINERS_HXX

#include <cassert>
#include <functional>
#include <array>
#include <vector>
#include <tuple>
#include <unordered_set>
#include <unordered_map>

#include "lolita/lolita.hxx"
#include "lolita/lolita_numerics.hxx"

namespace lolita
{

    namespace pair
    {

        template<typename T>
        struct Pair
        {

            constexpr
            Pair()
            :
            i(T(0)),
            j(T(0))
            {}

            constexpr
            Pair(
                    T
                    i_arg,
                    T
                    j_arg
            )
            :
            i(i_arg),
            j(j_arg)
            {}

            constexpr
            Bool
            operator==(
                    Pair const &
                    other
            )
            const = default;
            constexpr
            Bool
            operator!=(
                    Pair const &
                    other
            )
            const = default;

            T i;

            T j;

        };

    }

    namespace set
    {

        template<typename T, typename... U>
        struct Set
        {

            constexpr
            Set(
                    T
                    arg,
                    U...
                    args
            )
                    :
                    data({arg, static_cast<T>(args)...})
            {}

            constexpr
            Bool
            operator==(
                    Set const &
                    other
            )
            const
            {
                return data == other.data;
            }

            constexpr
            Bool
            operator!=(
                    Set const &
                    other
            )
            const
            {
                return !(* this == other);
            }

            static constexpr
            auto
            size()
            {
                return Indx(1 + sizeof...(U));
            }

            constexpr
            auto
            get(
                    auto
                    index_arg
            )
            const
            {
                return data[index_arg];
            }

            constexpr
            auto
            get(
                    auto
                    index_arg
            )
            {
                return data[index_arg];
            }

            std::array<T, size()> data;

        };

    }

    namespace array
    {

        template<typename T, auto... Is>
        struct Array
        {

            static_assert(0 < sizeof...(Is) <= 2);

            using Type = T;

            using Data = std::array<T, numerics::prod(Is...)>;

            static constexpr
            auto
            size()
            {
                return numerics::prod(static_cast<Indx>(Is)...);
            }

            static constexpr
            auto
            order()
            {
                return Indx(sizeof...(Is));
            }

            static constexpr
            auto
            dim(
                    auto
                    i
            )
            {
                return std::array<Indx, order()>{static_cast<Indx>(Is)...}[i];
            }

            constexpr
            Bool
            operator==(
                    Array const &
                    other
            )
            const
            {
                return data == other.data;
            }
//            = default;

            constexpr
            Bool
            operator!=(
                    Array const &
                    other
            )
            const
            {
                return !(* this == other);
            }
//            = default;

            constexpr
            T const &
            get(
                    auto ...
                    is
            )
            const
            {
                static_assert(sizeof...(is) == order());
                auto const indices = std::array<Indx, order()>{static_cast<Indx>(is)...};
                if constexpr(order() == 1) {
                    return data[indices[0]];
                }
                else if constexpr(order() == 2) {
                    return data[indices[0] * dim(1) + indices[1]];
                }
                else if constexpr(order() == 3) {
                    return data[indices[0] * dim(1) * dim(2) + indices[1] * dim(2) + indices[2]];
                }
                else {

                }
            }

            constexpr
            T &
            get(
                    auto ...
                    is
            )
            {
                static_assert(sizeof...(is) == order());
                auto const indices = std::array<Indx, order()>{static_cast<Indx>(is)...};
                if constexpr(order() == 1) {
                    return data[indices[0]];
                }
                else if constexpr(order() == 2) {
                    return data[indices[0] * dim(1) + indices[1]];
                }
                else if constexpr(order() == 3) {
                    return data[indices[0] * dim(1) * dim(2) + indices[1] * dim(2) + indices[2]];
                }
                else {

                }
            }

            Data data;

        };

        template<typename T>
        struct Array<T>
        {

            using Type = T;

            using Data = std::vector<T>;

            constexpr
            Array()
                    :
                    data({})
            {}

            constexpr
            Array(
                    std::initializer_list<T> const &
                    list
            )
                    :
                    data(list)
            {}

            constexpr
            Array(
                    std::initializer_list<T> &&
                    list
            )
                    :
                    data(list)
            {}

            constexpr
            Array(
                    Indx
                    size_arg,
                    T const &
                    value_arg
            )
                    :
                    data(size_arg, value_arg)
            {}

            Bool
            operator==(
                    Array const &
                    other
            )
            const
            {
                return data == other.data;
            }
//            = default;

            Bool
            operator!=(
                    Array const &
                    other
            )
            const
            {
                return !(* this == other);
            }
//            = default;

            constexpr
            auto
            size()
            const
            {
                return data.size();
            }

            constexpr
            T const &
            get(
                    Indx
                    i
            )
            const
            {
                return data.operator[](i);
            }

            constexpr
            T &
            get(
                    Indx
                    i
            )
            {
                return data.operator[](i);
            }

            Data data;

        };

        template<typename T>
        struct IsArray
        {
            auto const static constexpr value = false;
        };

        template<typename T, auto N>
        struct IsArray<std::array<T, N>>
        {
            auto const static constexpr value = true;
        };

        template<typename T, auto N>
        struct IsArray<Array<T, N>>
        {
            auto const static constexpr value = true;
        };

        template<typename T, typename U>
        static constexpr
        Indx
        index(
                T const &
                array_arg,
                U const &
                arg
        )
        {
            auto itr = std::find(array_arg.data.begin(), array_arg.data.end(), arg);
            if (itr == array_arg.data.end()) {
                return array_arg.size();
            }
            else {
                return Indx(std::distance(array_arg.data.begin(), itr));
            }
        }

        template<auto I, typename T>
        static constexpr
        auto
        range(
                T
                start_arg
        )
        {
            auto a = Array<T, I>();
            for (T i = 0; i < I; ++i) {
                a.get(i) = i + start_arg;
            }
            return a;
        }

        template<typename T>
        static constexpr
        auto
        range(
                T
                start_arg,
                T
                last_arg
        )
        {
            auto a = Array<T>();
            for (T i = 0; i < last_arg - start_arg; ++i) {
                a.data.push_back(i + start_arg);
            }
            return a;
        }

        template<typename T>
        static constexpr
        auto
        max(
                T const &
                array_arg
        )
        {
            auto value = typename T::Data(0);
            for (auto const & i: array_arg.data) {
                if (i > value) {
                    value = i;
                }
            }
            return value;
        }

        template<typename T>
        static constexpr
        auto
        sum(
                T const &
                array_arg
        )
        {
            auto value = typename T::Data(0);
            for (auto const & i: array_arg.data) {
                if (i > value) {
                    value += i;
                }
            }
            return value;
        }

    }

    namespace collection
    {

        template<typename... T>
        struct Collection
        {

            using Data = std::tuple<T...>;

            template<Indx I>
            using Type = typename std::tuple_element<I, Data>::type;

            static constexpr
            auto
            size()
            {
                return Indx(sizeof...(T));
            }

            constexpr
            Collection()
            :
            data()
            {}

            constexpr
            Collection(
                    T &&...
                    args
            )
            :
            data(std::forward<T>(args)...)
            {}

            constexpr
            Collection(
                    T const &...
                    args
            )
            :
            data(args...)
            {}

            constexpr
            Bool
            operator==(
                    Collection const &
                    other
            )
            const
            {
                return data == other.data;
            }

            constexpr
            Bool
            operator!=(
                    Collection const &
                    other
            )
            const
            {
                return !(* this == other);
            }

            template<Indx I>
            constexpr
            auto &
            get()
            {
                return std::get<I>(data);
            }

            template<Indx I>
            constexpr
            auto const &
            get()
            const
            {
                return std::get<I>(data);
            }

            Data data;

        };

        template<>
        struct Collection<>
        {

            static constexpr
            auto
            size()
            {
                return 0;
            }

            constexpr
            Collection()
            {}

            constexpr
            Bool
            operator==(
                    Collection const &
                    other
            )
            const
            {
                return true;
            }

            constexpr
            Bool
            operator!=(
                    Collection const &
                    other
            )
            const
            {
                return !(* this == other);
            }

        };

        template<typename... T>
        struct IsAggregate;

        template<typename... T>
        struct IsCollection
        {
            auto const static constexpr value = false;
        };

        template<typename... T>
        struct IsCollection<Collection<T...>>
        {
            auto const static constexpr value = true;
        };

        static constexpr
        void
        apply(
                auto const &
                collection_arg,
                auto &
                fun
        )
        {
            std::apply([&](auto const &... x){(..., fun(x));}, collection_arg.data);
        }

        static constexpr
        void
        apply(
                auto &
                collection_arg,
                auto &
                fun
        )
        {
            std::apply([&](auto &... x){(..., fun(x));}, collection_arg.data);
        }

        namespace detail
        {

            template<typename T, Indx I>
            static constexpr
            void
            apply(
                    auto &
                    fun
            )
            {
                using CollectionT = std::remove_cvref_t<T>;
                fun.template operator()<typename CollectionT::template Type<I>>();
                if constexpr (I < CollectionT::size() - 1) {
                    apply<T, I + 1>(fun);
                }
            }

            template<auto K, auto N>
            static constexpr
            void
            apply(
                    auto &
                    fun,
                    auto &&...
                    args
            )
            {
                fun.template operator()<K>(args...);
                if constexpr (K < N - 1) {
                    apply<K + 1, N>(fun, args...);
                }
            }

            template<auto K, auto N>
            static constexpr
            void
            apply(
                    auto &
                    fun,
                    auto const & col
            )
            {
                fun.template operator()<K>(col.template get<K>());
                if constexpr (K < N - 1) {
                    apply<K + 1, N>(fun);
                }
            }

            template<auto K, auto N>
            static constexpr
            void
            apply(
                    auto &
                    fun,
                    auto & col
            )
            {
                fun.template operator()<K>(col.template get<K>());
                if constexpr (K < N - 1) {
                    apply<K + 1, N>(fun);
                }
            }

        }

        template<typename T>
        static constexpr
        void
        apply(
                auto &
                fun
        )
        {
            detail::apply<T, 0>(fun);
        }

        template<auto N, auto K = 0>
        static constexpr
        void
        apply(
                auto &
                fun,
                auto &&...
                args
        )
        {
            detail::apply<K, N>(fun, args...);
        }

//        template<auto N>
//        static constexpr
//        void
//        apply(
//                auto &
//                fun,
//                auto & col
//        )
//        {
//            detail::apply<0, N>(fun, col);
//        }
//
//        template<auto N>
//        static constexpr
//        void
//        apply(
//                auto &
//                fun,
//                auto const & col
//        )
//        {
//            detail::apply<0, N>(fun, col);
//        }

        template<typename T, typename U>
        static constexpr
        Indx
        index()
        {
            using CollectionT = std::remove_cvref_t<T>;
            auto index = CollectionT::size();
            auto count = Indx(0);
            auto set_index = [& index, & count] <typename V>() constexpr mutable {
                if constexpr (std::is_same_v<std::remove_cvref_t<U>, std::remove_cvref_t<V>>) {
                    index = count;
                }
                count ++;
            };
            apply<CollectionT>(set_index);
            return index;
        }

        template<typename T, typename U>
        static constexpr
        Indx
        index(
                T const &
                C,
                U const &
                E
        )
        {
            auto index = C.size();
            auto count = Indx(0);
            auto find_index = [& index, & count, & E](auto const & x) constexpr mutable {
                if constexpr (std::is_same_v<U, std::remove_cvref_t<decltype(x)>>) {
                    if (x == E) {
                        index = count;
                    }
                }
                count ++;
            };
            apply(C, find_index);
            return index;
        }

        template<typename T, typename U>
        static constexpr
        Indx
        has()
        {
            using CollectionT = std::remove_cvref_t<T>;
            auto value = false;
            auto find_type = [& value] <typename V>() constexpr mutable {
                if constexpr (std::is_same_v<std::remove_cvref_t<U>, std::remove_cvref_t<V>>) {
                    value = true;
                }
            };
            apply<CollectionT>(find_type);
            return value;
        }

        template<typename T, typename U>
        static constexpr
        Indx
        has(
                T const &
                C,
                U const &
                E
        )
        {
            auto value = false;
            auto find_index = [& value, & E](auto const & x) constexpr mutable {
                if constexpr (std::is_same_v<U, std::remove_cvref_t<decltype(x)>>) {
                    if (x == E) {
                        value = true;
                    }
                }
            };
            apply(C, find_index);
            return value;
        }

        template<auto StaticArrayArg>
        struct ArrayCollectionWrapper
        {

        private:

            template<typename Indices = std::make_index_sequence<StaticArrayArg.size()>>
            struct CollectionGenerator;

            template<Indx... I>
            struct CollectionGenerator<std::index_sequence<I...>>
            {

                template<auto... s>
                struct CollectionWrapper
                {

                    template<template<auto>typename T>
                    using Wrapper = Collection<T<s>...>;

                };

                using CollectionWrapperType = CollectionWrapper<StaticArrayArg.data[I]...>;

            };

        public:

            template<template<auto>typename T>
            using Wrapper = typename CollectionGenerator<>::CollectionWrapperType::template Wrapper<T>;

        };

        template<auto A, template<auto>typename T>
        using ArrayWrapper = typename collection::ArrayCollectionWrapper<A>::template Wrapper<T>;

    }

    namespace aggregate
    {

        namespace detail
        {

            template<Indx I, typename T>
            struct Item
            {

                constexpr
                Bool
                operator==(
                        Item const &
                )
                const = default;

                constexpr
                Bool
                operator!=(
                        Item const &
                )
                const = default;

                T value;

            };

            template<Indx I, typename... T>
            struct Implementation;

            template<Indx I>
            struct Implementation<I>
            {


                constexpr
                Bool
                operator==(
                        Implementation const &
                )
                const = default;

                constexpr
                Bool
                operator!=(
                        Implementation const &
                )
                const = default;

            };

            template<Indx I, typename T, typename... U>
            struct Implementation<I, T, U...> : public Item<I, T>, public Implementation<I + 1, U...>
            {


                constexpr
                Bool
                operator==(
                        Implementation const &
                )
                const = default;

                constexpr
                Bool
                operator!=(
                        Implementation const &
                )
                const = default;

            };

            template<Indx I, typename T, typename... U>
            constexpr
            auto &
            get(
                    Implementation<I, T, U...> &
                    aggregate_arg
            )
            {
                return aggregate_arg.Item<I, T>::value;
            }

            template<Indx I, typename T, typename... U>
            constexpr
            auto const &
            get(
                    Implementation<I, T, U...> const &
                    aggregate_arg
            )
            {
                return aggregate_arg.Item<I, T>::value;
            }

            template<typename F, typename Tuple, std::size_t... I>
            constexpr
            decltype(auto) apply_impl(
                    F &&
                    f,
                    Tuple &&
                    t,
                    std::index_sequence<I...>
            )
            {
                return std::invoke(std::forward<F>(f), t.template get<I>()...);
            }

        }

        template<typename... T>
        struct Aggregate
        {

            using Data = detail::Implementation<0, T...>;

            template<Indx I>
            using Type = typename std::tuple_element<I, std::tuple<T...>>::type;

            static constexpr
            Indx
            size()
            {
                return sizeof...(T);
            }

            constexpr
            Bool
            operator==(
                    Aggregate const &
                    other
            )
            const
            {
                return data == other.data;
            }
//            const = default;

            constexpr
            Bool
            operator!=(
                    Aggregate const &
                    other
            )
            const
            {
                return !(* this == other);
            }
//            const = default;

            template<Indx I>
            constexpr
            auto const &
            get()
            const
            {
                return aggregate::detail::get<I>(data);
            }

            template<Indx I>
            constexpr
            auto &
            get()
            {
                return aggregate::detail::get<I>(data);
            }

            Data data;

        };

        template<typename... T>
        struct IsAggregate
        {
            auto const static constexpr value = false;
        };

        template<typename... T>
        struct IsAggregate<Aggregate<T...>>
        {
            auto const static constexpr value = true;
        };

        namespace detail
        {

            template<typename F, typename Tuple>
            constexpr
            decltype(auto)
            apply(
                    F &&
                    f,
                    Tuple &&
                    t
            )
            {
                auto const constexpr size = std::remove_cvref_t<Tuple>::size();
                return detail::apply_impl(std::forward<F>(f), std::forward<Tuple>(t), std::make_index_sequence<size>{});
            }

        }

        template<typename F, typename Tuple>
        constexpr
        decltype(auto)
        apply(
                F &&
                f,
                Tuple &&
                t
        )
        {
            detail::apply([&](auto const &... x){(..., f(x));}, t);
        }

        template<typename T, typename U>
        static constexpr
        Indx
        index()
        {
            auto index = std::remove_cvref_t<T>::size();
            auto count = Indx(0);
            auto find_index = [& index, & count](auto const & x) constexpr mutable {
                if constexpr (std::is_same_v<std::remove_cvref_t<U>, std::remove_cvref_t<decltype(x)>>) {
                    index = count;
                }
                count ++;
            };
            //apply([&](auto const &... x){(..., find_index(x));}, T());
            apply(find_index, T());
            return index;
        }

        template<typename T, typename U>
        static constexpr
        Indx
        index(
                T const &
                C,
                U const &
                E
        )
        {
            auto index = std::remove_cvref_t<T>::size();
            auto count = Indx(0);
            auto find_index = [& index, & count, & E](auto const & x) constexpr mutable {
                if constexpr (std::is_same_v<U, std::remove_cvref_t<decltype(x)>>) {
                    if (x == E) {
                        index = count;
                    }
                }
                count ++;
            };
            //apply([&](auto const &... x){(..., find_index(x));}, C);
            apply(find_index, C);
            return index;
        }

        template<typename T, typename U>
        static constexpr
        Indx
        has()
        {
            auto value = false;
            auto find_index = [& value](auto const & x) constexpr mutable {
                if constexpr (std::is_same_v<std::remove_cvref_t<U>, std::remove_cvref_t<decltype(x)>>) {
                    value = true;
                }
            };
            //apply([&](auto const &... x){(..., find_index(x));}, T());
            apply(find_index, T());
            return value;
        }

        template<typename T, typename U>
        static constexpr
        Indx
        has(
                T const &
                C,
                U const &
                E
        )
        {
            auto value = false;
            auto find_index = [& value, & E](auto const & x) constexpr mutable {
                if constexpr (std::is_same_v<U, std::remove_cvref_t<decltype(x)>>) {
                    if (x == E) {
                        value = true;
                    }
                }
            };
            //apply([&](auto const &... x){(..., find_index(x));}, C);
            apply(find_index, C);
            return value;
        }

    }

    namespace unordered_map
    {

        template<typename K, typename T>
        struct UnorderedMap
        {

            using Type = T;

            using Data = std::unordered_map<K, T>;

            UnorderedMap() = default;

            Bool
            operator==(
                    UnorderedMap const &
                    other
            )
//            const = default;
            const
            {
                return data == other.data;
            }

            Bool
            operator!=(
                    UnorderedMap const &
                    other
            )
//            const = default;
            const
            {
                return !(* this == other);
            }

            T &
            get(
                    K const &
                    i
            )
            {
                return data.operator[](i);
            }

            T &
            get(
                    K &&
                    i
            )
            {
                return data.operator[](i);
            }

            T const &
            get(
                    K const &
                    i
            )
            const
            {
                return data.at(i);
            }

            T const &
            get(
                    K &&
                    i
            )
            const
            {
                return data.at(i);
            }

            Indx
            getSize()
            const
            {
                return data.size();
            }

            Indx
            size()
            const
            {
                return data.size();
            }

            Data data;

        };

    }

    namespace unordered_set
    {

        template<typename Type, typename HashType>
        struct UnorderedSetHashHelper
        {

            using UnorderedSetDataT = std::unordered_set<Type, HashType>;

        };

        template<typename Type>
        struct UnorderedSetHashHelper<Type, void>
        {

            using UnorderedSetDataT = std::unordered_set<Type>;

        };

        template<typename T, typename HashType = void>
        struct UnorderedSet
        {

            using Type = T;

            using Data = typename UnorderedSetHashHelper<T, HashType>::UnorderedSetDataT;

            UnorderedSet() = default;

            UnorderedSet(
                    std::initializer_list<T> &&
                    init
            )
                    :
                    data(init)
            {}

            Bool
            operator==(
                    UnorderedSet const &
                    other
            )
//            const = default;
            const
            {
                return data == other.data;
            }

            Bool
            operator!=(
                    UnorderedSet const &
                    other
            )
//            const = default;
            const
            {
                return !(* this == other);
            }

            T &
            get(
                    Indx
                    i
            )
            {
                auto it = data.begin();
                std::advance(it, i);
                return data.extract(it).value();
            }

            T const &
            get(
                    Indx
                    i
            )
            const
            {
                auto it = data.begin();
                std::advance(it, i);
                return data.extract(it).value();
            }

            constexpr
            Indx
            getSize()
            const
            {
                return data.size();
            }

            Data data;

        };

        template<typename Type>
        using PointerUnorederedSet = UnorderedSet<Type, typename Type::Hash>;

    }

    template<typename T>
    using Pair = pair::Pair<T>;

    template<typename T, typename... U>
    using Set = set::Set<T, U...>;

    template<typename T, auto... I>
    using Array = array::Array<T, I...>;

    template<typename... T>
    using Collection = collection::Collection<T...>;

    template<typename... T>
    using Aggregate = aggregate::Aggregate<T...>;

    template<typename K, typename T>
    using UnorderedMap = unordered_map::UnorderedMap<K, T>;

    template<typename T, typename HashType = void>
    using UnorderedSet = unordered_set::UnorderedSet<T, HashType>;

}

#endif //LOLITA_LOLITA_CONTAINERS_HXX
