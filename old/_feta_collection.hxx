//
// Created by dsiedel on 09/02/2022.
//

#ifndef FETA__FETA_COLLECTION_HXX
#define FETA__FETA_COLLECTION_HXX

//#include <unordered_set>
//#include <variant>
//#include <ostream>
#include "new/_feta.hxx"
#include "new/_feta_array.hxx"

namespace feta
{

    template<typename... T>
    struct Collection
    {

        using Data = std::tuple<T...>;

        template<Indx I>
        using Type = typename std::tuple_element<I, Data>::type;

        static constexpr
        auto
        getSize()
        {
            return Indx(sizeof...(T));
        }

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
            auto eq_0 = data == other.data;
            return eq_0;
        }

        constexpr
        Bool
        operator!=(
                Collection const &
                other
        )
        const
        {
            return !(other == * this);
        }

//        template<Indx I>
//        constexpr
//        auto &
//        operator() ()
//        {
//            return std::get<I>(data);
//        }
//
//        template<Indx I>
//        constexpr
//        auto const &
//        operator() ()
//        const
//        {
//            return std::get<I>(data);
//        }

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

        template<Indx I = 0>
        void
        loop(
                auto
                fun,
                auto
                ...
                args
        )
        const
        {
            fun(get<I>(), args...);
            if constexpr (I < size() - 1) {
                loop<I + 1>(fun , args...);
            }
        }

        /**
         * @brief
         */
        Data data;

    };

    template<>
    struct Collection<>
    {

        static constexpr
        auto
        getSize()
        {
            return 0;
        }

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
            return !(other == * this);
        }

    };

//    template<typename T>
//    struct Collection<Collection<T>>
//    {
//
//        using Data = std::tuple<Collection<T>>;
//
//        template<Indx I>
//        using Type = typename std::tuple_element<I, Data>::type;
//
//        static constexpr
//        auto
//        getSize()
//        {
//            return 1;
//        }
//
//        constexpr
//        Collection()
//        :
//        data()
//        {}
//
//        constexpr
//        Collection(
//                Collection<T> &&
//                arg
//        )
//        :
//        data(arg)
//        {
//            print("------------------------------------------- here move");
//        }
//
//        constexpr
//        Collection(
//                Collection<T> const &
//                arg
//        )
//        :
//        data(arg)
//        {
//            print("------------------------------------------- here copy");
//        }
//
//        template<Indx I>
//        constexpr
//        auto &
//        get()
//        {
//            static_assert(I == 0);
//            return std::get<I>(data);
//        }
//
//        template<Indx I>
//        constexpr
//        auto const &
//        get()
//        const
//        {
//            static_assert(I == 0);
//            return std::get<I>(data);
//        }
//
//        Data data;
//
//    };

//    template<typename...>
//    struct CollectionUnionPolicy;
//
//    template<typename ...T, typename ...U>
//    struct CollectionUnionPolicy<Collection<T...>, Collection<U...>>
//    {
//
//        using Union = Collection<T..., U...>;
//
//    };
//
//    template<typename ...T, typename ...U, typename ...V>
//    struct CollectionUnionPolicy<Collection<T...>, Collection<U...>, Collection<V...>>
//    {
//
//        using Union = Collection<T..., U..., V...>;
//
//    };
//
////    template<typename T, typename U>
////    using CollectionUnion = typename CollectionUnionPolicy<T, U>::Union;
//
//    template<typename ...T>
//    using CollectionUnion = typename CollectionUnionPolicy<T...>::Union;








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

    template<Indx I>
    struct TESTT
    {
        Indx const static constexpr value = I;
    };














//    template<
//            auto StaticArrayArg
//    >
//    struct ArrayCollectionWrapper
//    {
//
//    private:
//
//        using StructuralTypeF = typename StructuralType<StaticArrayArg>::Type;
//
//        template<
//                typename Indices = std::make_index_sequence<
//                        StaticArrayArg.getSize()
//                >
//        >
//        struct CollectionGenerator;
//
//        template<
//                Indx... I
//        >
//        struct CollectionGenerator<
//                std::index_sequence<
//                        I...
//                >
//        >
//        {
//
//            template<
//                    StructuralTypeF... s
//            >
//            struct CollectionWrapper
//            {
//
//                template<
//                        template<
//                                StructuralTypeF
//                        >
//                        typename T
//                >
//                using Wrapper = Collection<
//                        T<
//                                s
//                        >...
//                >;
//
//            };
//
//            using CollectionWrapperType = CollectionWrapper<
//                    StaticArrayArg.data[I]...
//            >;
//
//        };
//
//    public:
//
//        template<
//                template<
//                        StructuralTypeF
//                >
//                typename T
//        >
//        using Wrapper = typename CollectionGenerator<>::CollectionWrapperType::template Wrapper<
//                T
//        >;
//
//    };


}

#endif //FETA__FETA_COLLECTION_HXX
