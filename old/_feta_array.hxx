//
// Created by dsiedel on 24/03/2022.
//

#ifndef FETA__FETA_ARRAY_HXX
#define FETA__FETA_ARRAY_HXX

#include <array>
#include <vector>

#include "new/_feta.hxx"

namespace feta
{

//    template<typename T, auto... Is>
//    struct Array;

    namespace detail
    {

        template<typename T, auto... I>
        struct ArrayImpl
        {

            static_assert(0 < sizeof...(I) <= 2);

            /**
             * @brief
             */
            using Type = T;

            /**
             * @brief
             */
            using Data = std::array<T, getProduct(static_cast<Indx>(I)...)>;

            /**
             * @brief
             */
            using Shape = std::array<Indx, sizeof...(I)>;

            static constexpr
            Indx
            getOrder()
            {
                return sizeof...(I);
            }

            /***
             * @brief
             * @return
             */
            static constexpr
            Shape
            getShape()
            {
                return {I...};
            }

            /***
             * @brief
             * @return
             */
            static constexpr
            Indx
            getSize()
            {
                return getProduct(static_cast<Indx>(I)...);
            }

            constexpr
            T const &
            operator() (
                    auto ...
                    is
            )
            const
            {
                return fetch(is...);
            }

            constexpr
            T &
            operator() (
                    auto ...
                    is
            )
            {
                return fetch(is...);
            }

            constexpr
            T const &
            get(
                    auto ...
                    is
            )
            const
            {
                return fetch(is...);
            }

            constexpr
            T &
            get(
                    auto ...
                    is
            )
            {
                return fetch(is...);
            }

        private:

            constexpr
            T &
            fetch(
                    auto ...
                    is
            )
            {
                static_assert(sizeof...(is) == getOrder());
                auto const indices = Shape{static_cast<Indx>(is)...};
                if constexpr(getOrder() == 1) {
                    return data[indices[0]];
                }
                if constexpr(getOrder() == 2) {
                    return data[indices[0] * getShape()[1] + indices[1]];
                }
                if constexpr(getOrder() == 3) {
                    return data[indices[0] * getShape()[1] * getShape()[2] + indices[1] * getShape()[2] + indices[2]];
                }
            }

            constexpr
            T const &
            fetch(
                    auto ...
                    is
            )
            const
            {
                static_assert(sizeof...(is) == getOrder());
                auto const indices = Shape{static_cast<Indx>(is)...};
                if constexpr(getOrder() == 1) {
                    return data[indices[0]];
                }
                if constexpr(getOrder() == 2) {
                    return data[indices[0] * getShape()[1] + indices[1]];
                }
                if constexpr(getOrder() == 3) {
                    return data[indices[0] * getShape()[1] * getShape()[2] + indices[1] * getShape()[2] + indices[2]];
                }
            }

        public:

            Data data;

        };

    }

//    template<typename T, auto... Is>
//    struct Array : public detail::ArrayImpl<T, Is...>
//    {};
//
//    template<auto... Is>
//    struct Array<Real, Is...> : public detail::ArrayImpl<Real, Is...>
//    {
//
//        constexpr
//        Array()
//        :
//        data(Data{})
//        {}
//
//        constexpr
//        Array(
//                auto &&...
//                values
//        )
//        :
//        data(Data{static_cast<Type>(values)...})
////        data(Data{values...})
//        {}
//
//        constexpr
//        Array(
//                auto const &...
//                values
//        )
//        :
//        data(Data{static_cast<Type>(values)...})
////                data(Data{values...})
//        {}
//
//    };

    template<typename T, auto... I>
    struct Array2Pol;

    template<typename T, auto... I>
    using Array2 = typename Array2Pol<T, I...>::Type;

    template<typename T, auto... I>
    struct Array2Pol
    {

        using Type = std::array<T, getProduct(static_cast<Indx>(I)...)>;

    };

    template<typename T>
    struct Array2Pol<T>
    {

        using Type = std::vector<T>;

    };

    template<typename T, auto... I>
    struct ArrayMap
    {

        static constexpr
        auto
        order()
        {
            return Indx(sizeof...(I));
        }

        static constexpr
        auto
        shape()
        {
            return Array2<T, sizeof...(I)>{I...};
        }

        static constexpr
        auto
        size()
        {
            return getProduct(static_cast<Indx>(I)...);
        }

        constexpr
        ArrayMap(
                Array2<T, I...> const &
                array_arg
        )
        :
        data(* array_arg)
        {}

        constexpr
        T const &
        operator() (
                auto ...
                is
        )
        const
        {
            return fetch(is...);
        }

        constexpr
        T &
        operator() (
                auto ...
                is
        )
        {
            return fetch(is...);
        }

        constexpr
        T const &
        get(
                auto ...
                is
        )
        const
        {
            return fetch(is...);
        }

        constexpr
        T &
        get(
                auto ...
                is
        )
        {
            return fetch(is...);
        }

    private:

        constexpr
        T &
        fetch(
                auto ...
                is
        )
        {
            static_assert(sizeof...(is) == order());
            auto const indices = Array2<Indx, order()>{static_cast<Indx>(is)...};
            if constexpr(order() == 1) {
                return data[indices[0]];
            }
            if constexpr(order() == 2) {
                return data[indices[0] * shape()[1] + indices[1]];
            }
            if constexpr(order() == 3) {
                return data[indices[0] * shape()[1] * shape()[2] + indices[1] * shape()[2] + indices[2]];
            }
        }

        constexpr
        T const &
        fetch(
                auto ...
                is
        )
        const
        {
            static_assert(sizeof...(is) == order());
            auto const indices = Array2<Indx, order()>{static_cast<Indx>(is)...};
            if constexpr(order() == 1) {
                return data[indices[0]];
            }
            if constexpr(order() == 2) {
                return data[indices[0] * shape()[1] + indices[1]];
            }
            if constexpr(order() == 3) {
                return data[indices[0] * shape()[1] * shape()[2] + indices[1] * shape()[2] + indices[2]];
            }
        }

        Array2<T, I...> & data;

    };

    template<auto data, auto... I>
    struct ArrayConstMap
    {

        static constexpr
        auto
        order()
        {
            return Indx(sizeof...(I));
        }

        static constexpr
        auto
        shape()
        {
            return Array2<typename decltype(data)::type, sizeof...(I)>{I...};
        }

        static constexpr
        auto
        size()
        {
            return getProduct(static_cast<Indx>(I)...);
        }

        constexpr
        auto const &
        operator() (
                auto ...
                is
        )
        const
        {
            return fetch(is...);
        }

        constexpr
        auto &
        operator() (
                auto ...
                is
        )
        {
            return fetch(is...);
        }

        constexpr
        auto const &
        get(
                auto ...
                is
        )
        const
        {
            return fetch(is...);
        }

        constexpr
        auto &
        get(
                auto ...
                is
        )
        {
            return fetch(is...);
        }

    private:

        constexpr
        auto &
        fetch(
                auto ...
                is
        )
        {
            static_assert(sizeof...(is) == order());
            auto const indices = Array2<Indx, order()>{static_cast<Indx>(is)...};
            if constexpr(order() == 1) {
                return data[indices[0]];
            }
            if constexpr(order() == 2) {
                return data[indices[0] * shape()[1] + indices[1]];
            }
            if constexpr(order() == 3) {
                return data[indices[0] * shape()[1] * shape()[2] + indices[1] * shape()[2] + indices[2]];
            }
        }

        constexpr
        auto const &
        fetch(
                auto ...
                is
        )
        const
        {
            static_assert(sizeof...(is) == order());
            auto const indices = Array2<Indx, order()>{static_cast<Indx>(is)...};
            if constexpr(order() == 1) {
                return data[indices[0]];
            }
            if constexpr(order() == 2) {
                return data[indices[0] * shape()[1] + indices[1]];
            }
            if constexpr(order() == 3) {
                return data[indices[0] * shape()[1] * shape()[2] + indices[1] * shape()[2] + indices[2]];
            }
        }

    };

    template<typename T, auto... Is>
    struct Array
    {

        static_assert(0 < sizeof...(Is) <= 2);

        /**
         * @brief
         */
        using Type = T;

        /**
         * @brief
         */
        using Data = std::array<T, getProduct(static_cast<Indx>(Is)...)>;

        /***
         * @brief
         * @return
         */
        static constexpr
        auto
        getSize()
        {
            return getProduct(static_cast<Indx>(Is)...);
        }

        /***
         * @brief
         * @return
         */
        static constexpr
        auto
        size()
        {
            return getProduct(static_cast<Indx>(Is)...);
        }

        static constexpr
        auto
        order()
        {
            return Indx(sizeof...(Is));
        }

        /***
         * @brief
         * @return
         */
        static constexpr
        auto
        dim(
                auto
                i
        )
        {
            return std::array<Indx, order()>{static_cast<Indx>(Is)...}[i];
        }

//        /***
//         * @brief
//         * @return
//         */
//        static constexpr
//        Shape
//        getShape()
//        {
//            return {Is...};
//        }

//        constexpr
//        Array()
//        :
//        data(Data{})
//        {}
//
//        constexpr
//        Array(
//                auto &&...
//                values
//        )
//        :
//        data(Data{static_cast<Type>(values)...})
////        data(Data{values...})
//        {}
//
//        constexpr
//        Array(
//                auto const &...
//                values
//        )
//                :
//        data(Data{static_cast<Type>(values)...})
////                data(Data{values...})
//        {}

        constexpr
        T const &
        operator() (
                auto ...
                is
        )
        const
        {
            return fetch(is...);
        }

        constexpr
        T &
        operator() (
                auto ...
                is
        )
        {
            return fetch(is...);
        }

        constexpr
        T const &
        get(
                auto ...
                is
        )
        const
        {
            return fetch(is...);
        }

        constexpr
        T &
        get(
                auto ...
                is
        )
        {
            return fetch(is...);
        }

//        constexpr
//        Bool
//        operator==(
//                Array const &
//        )
//        const = default;
//
//        constexpr
//        Bool
//        operator!=(
//                Array const &
//        )
//        const = default;

        constexpr
        Bool
        operator==(
                Array const &
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
                Array const &
                other
        )
        const
        {
            return !(other == * this);
        }

    private:

        constexpr
        T &
        fetch(
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

        constexpr
        T const &
        fetch(
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

    public:

        Data data;

    };

    template<typename T>
    struct Array<T>
    {

        using Type = T;

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

//        constexpr
//        bool operator==(
//                Array const &
//        )
//        const = default;
//
//        constexpr
//        bool operator!=(
//                Array const &
//        )
//        const = default;

//        Bool
//        operator==(
//                Array const &
//        )
//        const = default;
//
//        Bool
//        operator!=(
//                Array const &
//        )
//        const = default;

        Bool
        operator==(
                Array const &
                other
        )
        const
        {
            auto eq_0 = data == other.data;
            return eq_0;
        }

        Bool
        operator!=(
                Array const &
                other
        )
        const
        {
            return !(other == * this);
        }

        constexpr
        Indx
        getSize()
        const
        {
            return data.size();
        }

        constexpr
        T const &
        operator() (
                Indx
                i
        )
        const
        {
            return fetch(i);
        }

        constexpr
        T &
        operator() (
                Indx
                i
        )
        {
            return fetch(i);
        }

        constexpr
        T const &
        get(
                Indx
                i
        )
        const
        {
            return fetch(i);
        }

        constexpr
        T &
        get(
                Indx
                i
        )
        {
            return fetch(i);
        }

        constexpr
        Indx
        index(
                T const &
                t
        )
        const
        {
            auto itr = std::find(data.begin(), data.end(), t);
            if (itr == data.end()) {
                return getSize();
            } else {
                return Indx(std::distance(data.begin(), itr));
            }
        }

    private:

        T &
        fetch(
                Indx
                i
        )
        {
            return data.operator[](i);
        }

        T const &
        fetch(
                Indx
                i
        )
        const
        {
            return data.operator[](i);
        }

    public:

        using Data = std::vector<
                T
        >;

        Data data;

//        Shape shape;

    };

}

#endif //FETA__FETA_ARRAY_HXX
