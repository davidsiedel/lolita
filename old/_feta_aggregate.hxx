//
// Created by dsiedel on 27/02/2022.
//

#ifndef FETA__FETA_AGGREGATE_HXX
#define FETA__FETA_AGGREGATE_HXX

#include "new/_feta.hxx"
#include "new/_feta_collection.hxx"

namespace feta
{

    namespace aggregate::detail
    {

        template<Indx I, typename T>
        struct Item
        {

            T value;

        };

        template<Indx I, typename... T>
        struct Implementation;

        template<Indx I>
        struct Implementation<I>
        {};

        template<Indx I, typename T, typename... U>
        struct Implementation<I, T, U...> : public Item<I, T>, public Implementation<I + 1, U...>
        {};

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

    }

    template<typename... T>
    struct Aggregate
    {

    private:

        using DataT = aggregate::detail::Implementation<0, T...>;

    public:

        template<Indx I>
        using Type = typename std::tuple_element<I, std::tuple<T...>>::type;

        constexpr
        Bool
        operator==(
                Aggregate const &
                other
        )
        const
        {
            auto eq_0 = data == other.data;
            return eq_0;
        }

        constexpr
        Bool
        operator!=(
                Aggregate const &
                other
        )
        const
        {
            return !(other == * this);
        }

//        constexpr
//        bool operator==(
//                Aggregate const &
//        )
//        const = default;
//
//        constexpr
//        bool operator!=(
//                Aggregate const &
//        )
//        const = default;

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

    private:

        template<Indx I = 0>
        constexpr
        void
        setIndex(
                auto const &
                E,
                auto &
                i
        )
        const
        {
            if constexpr (std::is_same_v<std::remove_cvref_t<decltype(E)>, Type<I>>) {
                if (get<I>() == E) {
                    i = I;
                }
            }
            if constexpr (I < size() - 1) {
                setIndex<I + 1>(E, i);
            }
        }

    public:

        constexpr
        auto
        index(
                auto const &
                E
        )
        const
        {
            auto i = size();
            setIndex(E, i);
            return i;
        }

        constexpr
        auto
        has(
                auto const &
                E
        )
        const
        {
            if (index(E) < size()) {
                return true;
            } else {
                return false;
            }
        }

        static constexpr
        Indx
        getSize()
        {
            return sizeof...(T);
        }

        static constexpr
        Indx
        size()
        {
            return sizeof...(T);
        }

    private:

        template<typename U, Indx I = 0>
        static constexpr
        void
        hasType(
                auto &
                value_arg
        )
        {
            if constexpr (std::is_same_v<U, Type<I>>) {
                value_arg = true;
            }
            if constexpr (I < getSize() - 1) {
                hasType<I + 1>(value_arg);
            }
        }

    public:

        template<typename U>
        static constexpr
        auto
        has()
        {
            auto value_arg = false;
            hasType<U>(value_arg);
            return value_arg;
        }

        DataT data;

    };





    template<auto StaticArrayArg>
    struct TupleCollectionWrapper
    {

    private:

        template<typename Indices = std::make_index_sequence<StaticArrayArg.getSize()>>
        struct CollectionGenerator;

        template<Indx... I>
        struct CollectionGenerator<std::index_sequence<I...>>
        {

            template<auto... A>
            struct CollectionWrapper
            {

                template<template<auto> typename T>
                using Wrapper = Collection<T<A>...>;

            };

//            using CollectionWrapperType = CollectionWrapper<StaticArrayArg.template get<I>()...>;

            using CollectionWrapperType = CollectionWrapper<aggregate::detail::get<I>(StaticArrayArg.data)...>;

        };

    public:

        template<template<auto> typename T>
        using Wrapper = typename CollectionGenerator<>::CollectionWrapperType::template Wrapper<T>;

    };

}

#endif //FETA__FETA_AGGREGATE_HXX
