//
// Created by dsiedel on 04/05/22.
//

#ifndef LOLITA_LOLITA_HXX
#define LOLITA_LOLITA_HXX

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>
#include <ostream>



//#define EIGEN_USE_MKL_ALL
//#define EIGEN_USE_LAPACKE
//#define EIGEN_USE_MKL_VML
#define EIGEN_USE_BLAS
#define EIGEN_USE_MKL

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/PardisoSupport>

//#include "lolita_matrix.hxx"
//#include "lolita_pointers.hxx"

namespace lolita
{

    namespace config
    {

        using character = char;

        using integer = int;

        using index = unsigned short;

        using large = unsigned long long;

        using real = double;

        using boolean = bool;

    }

    using namespace config;

    namespace utility
    {

        template<typename T>
        using Label = std::array<T, 150>;

        template<typename T>
        static constexpr
        Label<T>
        makeLabel(
                std::basic_string_view<T> &&
                str
        )
        {
            auto label = Label<T>();
            for (auto i = 0; i < label.size(); ++i) {
                i < str.size() ? label[i] = str[i] : label[i] = '#';
            }
            return label;
        }

        namespace detail
        {

            template<std::unsigned_integral auto I, typename T>
            struct AggregateHead
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

                T value_;

            };

            template<std::unsigned_integral auto I, typename... T>
            struct AggregateTail;

            template<std::unsigned_integral auto I>
            struct AggregateTail<I>
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

            };

            template<std::unsigned_integral auto I, typename T, typename... U>
            struct AggregateTail<I, T, U...> : public AggregateHead<I, T>, public AggregateTail<I + 1, U...>
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

            };

            template<std::unsigned_integral auto I, typename T, typename... U>
            constexpr
            auto &
            get(
                    AggregateTail<I, T, U...> &
                    aggregate_tail
            )
            {
                return aggregate_tail.AggregateHead<I, T>::value_;
            }

            template<std::unsigned_integral auto I, typename T, typename... U>
            constexpr
            auto const &
            get(
                    AggregateTail<I, T, U...> const &
                    aggregate_tail
            )
            {
                return aggregate_tail.AggregateHead<I, T>::value_;
            }

        }

        template<typename... T>
        struct Aggregate
        {

            using Data = detail::AggregateTail<0, T...>;

            template<lolita::index I>
            using Type = typename std::tuple_element<I, std::tuple<T...>>::type;

            static constexpr
            lolita::index
            size()
            {
                return sizeof...(T);
            }

            constexpr
            lolita::boolean
            operator==(
                    Aggregate const &
                    other
            )
            const
            {
                return data_ == other.data;
            }
//            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    Aggregate const &
                    other
            )
            const
            {
                return !(* this == other);
            }
//            const = default;

            template<std::unsigned_integral I>
            constexpr
            auto const &
            get()
            const
            {
                return detail::get<I>(data_);
            }

            template<std::unsigned_integral I>
            constexpr
            auto &
            get()
            {
                return detail::get<I>(data_);
            }

            Data data_;

        };

    }

    namespace matrix
    {



    }

}

namespace lolita::core
{

    struct Quadrature
    {

        enum Type
        {

            Gauss,

        };

        constexpr
        Quadrature(
                lolita::index
                dim,
                Quadrature::Type
                typ
        )
        :
        dim_(dim),
        typ_(typ)
        {}

        constexpr
        lolita::boolean
        operator==(
                Quadrature const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                Quadrature const &
                other
        )
        const = default;

        lolita::index dim_;

        Quadrature::Type typ_;

    };

    struct Domain
    {

        enum Type
        {

            AxiSymmetric,
            Cartesian,

        };

        constexpr
        Domain(
                std::basic_string_view<lolita::character> &&
                tag,
                lolita::index
                dim,
                Type
                typ
        )
        :
        tag_(utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag))),
        dim_(dim),
        typ_(typ)
        {}

        constexpr
        lolita::boolean
        operator==(
                Domain const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                Domain const &
                other
        )
        const = default;

        constexpr
        lolita::index
        ordIntegration(
                lolita::index
                ord
        )
        const
        {
            return typ_ == Domain::AxiSymmetric ? 2 * ord + 1 : 2 * ord;
        }

        utility::Label<lolita::character> tag_;

        lolita::index dim_;

        Domain::Type typ_;

    };

    struct Mapping
    {

        enum Type
        {

            Gradient,
            Identity,

        };

        constexpr
        Mapping(
                Mapping::Type
                typ,
                lolita::index
                ord_field,
                lolita::index
                dim_field
        )
        :
        typ_(typ),
        ord_field_(ord_field),
        dim_field_(dim_field)
        {}

        Mapping::Type typ_;

        lolita::index ord_field_;

        lolita::index dim_field_;

    };

    namespace detail
    {

        struct FieldBase
        {

            enum Type
            {

                Material,
                External,
                Internal,

            };

            constexpr
            FieldBase(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim
            )
            :
            tag_(utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag))),
            dim_(dim)
            {}

            constexpr
            lolita::boolean
            operator==(
                    FieldBase const &
                    other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    FieldBase const &
                    other
            )
            const = default;

            utility::Label<lolita::character> tag_;

            lolita::index dim_;

        };

    }

    template<typename T>
    concept FieldBaseDerivedType = std::is_base_of_v<detail::FieldBase, T>;

    template<typename... MappingT>
    requires((std::same_as<MappingT, Mapping::Type> && ...) && sizeof...(MappingT) > 0)
    struct DegreeOfFreedom : public detail::FieldBase
    {

        constexpr
        DegreeOfFreedom(
                std::basic_string_view<lolita::character> &&
                tag,
                lolita::index
                dim,
                MappingT &&...
                mappings
        )
        :
        detail::FieldBase(std::forward<std::basic_string_view<lolita::character>>(tag), dim),
        mappings_({std::forward<Mapping::Type>(mappings)...})
        {}

        constexpr
        lolita::boolean
        operator==(
                DegreeOfFreedom const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                DegreeOfFreedom const &
                other
        )
        const = default;

        std::array<Mapping::Type, sizeof...(MappingT)> mappings_;

    };

    struct MaterialProperty : public detail::FieldBase
    {

        constexpr
        MaterialProperty(
                std::basic_string_view<lolita::character> &&
                tag
        )
        :
        detail::FieldBase(std::forward<std::basic_string_view<lolita::character>>(tag), 0)
        {}

        constexpr
        lolita::boolean
        operator==(
                MaterialProperty const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                MaterialProperty const &
                other
        )
        const = default;

    };

    struct SomeField : public detail::FieldBase
    {

        constexpr
        SomeField(
                std::basic_string_view<lolita::character> &&
                tag,
                lolita::index
                dim
        )
        :
        detail::FieldBase(std::forward<std::basic_string_view<lolita::character>>(tag), dim)
        {}

        constexpr
        lolita::boolean
        operator==(
                SomeField const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                SomeField const &
                other
        )
        const = default;

    };

    template<typename... FieldT>
    requires((std::derived_from<FieldT, detail::FieldBase> && ...) && sizeof...(FieldT) > 0)
    struct BehaviourData
    {

        constexpr
        BehaviourData(
                std::basic_string_view<lolita::character> &&
                tag,
                lolita::index
                dim
        )
        :
        tag_(utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag)))
        {}

        utility::Label<lolita::character> tag_;

    };

    template<auto A>
    struct S {};

    template<template<auto> typename T, auto... H>
    static constexpr
    std::tuple<T<H>...>
    expand(
//            auto...
//            items
    )
    {
        return std::tuple<T<H>...>();
    }

    template<typename... T>
    struct Holder;

    template<typename Ta, typename Tb>
    struct Holder<Ta, Tb>
    {

        Ta a;
        Tb b;

        template<std::unsigned_integral auto I>
        constexpr
        auto
        get()
        {
            if constexpr (I == 0) {
                return a;
            }
//            else if constexpr (I == 1) {
//                return b;
//            }
            else {
                return b;
            }
        }


    };

    template<template<auto> typename T, auto H>
    struct Expander
    {



    };

}

#endif //LOLITA_LOLITA_HXX
