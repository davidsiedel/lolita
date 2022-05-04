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

        using Char = char;

        using Intg = int;

        using Indx = std::size_t;

        using Long = unsigned long long;

        using Real = double;

        using Bool = bool;

        using Strg = std::basic_string<Char>;

        using Strv = std::basic_string_view<Char>;

        namespace detail
        {

            Indx const static constexpr name_size = 150;

        }

        using Name = std::array<Char, detail::name_size>;

    }

    using namespace config;

    namespace utility
    {

        static constexpr
        Name
        setName(
                Strv &&
                strv
        )
        {
            auto name = Name();
            for (auto i = 0; i < config::detail::name_size; ++i) {
                i < strv.size() ? name[i] = strv[i] : name[i] = '#';
            }
            return name;
        }

        static constexpr
        Strv
        getName(
                Name const &
                name
        )
        {
            return Strv(name.data(), std::distance(std::begin(name), std::find(std::begin(name), std::end(name), '#')));
        }

        namespace detail
        {

            template<typename T, auto... A>
            struct ArrayPolicy;

            template<typename T>
            struct ArrayPolicy<T>
            {

                using type = std::vector<T>;

            };

            template<typename T, std::integral auto M>
            struct ArrayPolicy<T, M>
            {

                using type = std::array<T, M>;

            };

            template<typename T, std::integral auto M, std::integral auto N>
            struct ArrayPolicy<T, M, N>
            {

                using type = std::array<std::array<T, M>, N>;

            };

        }

        template<typename T, auto... A>
        using Array = typename detail::ArrayPolicy<T, A...>::type;

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
                std::integral auto
                dim,
                Type
                typ
        )
        :
        dim_(dim),
        typ_(typ)
        {}

        constexpr
        Bool
        operator==(
                Quadrature const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Quadrature const &
                other
        )
        const = default;

        Indx dim_;

        Type typ_;

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
                Strv &&
                tag,
                std::integral auto
                dim,
                Type
                typ
        )
        :
        tag_(utility::setName(std::forward<Strv>(tag))),
        dim_(dim),
        typ_(typ)
        {}

        constexpr
        Bool
        operator==(
                Domain const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Domain const &
                other
        )
        const = default;

        constexpr
        Indx
        ordIntegration(
                std::integral auto
                ord
        )
        const
        {
            return typ_ == Type::AxiSymmetric ? 2 * ord + 1 : 2 * ord;
        }

        Name tag_;

        Indx dim_;

        Type typ_;

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
                Type
                typ,
                std::integral auto
                ord_field,
                std::integral auto
                dim_field
        )
        :
        typ_(typ),
        ord_field_(ord_field),
        dim_field_(dim_field)
        {}

        Type typ_;

        Indx ord_field_;

        Indx dim_field_;

    };

    template<typename... MappingT>
    requires(std::same_as<MappingT, Mapping::Type> && ...)
    struct TensorField
    {

        enum Type
        {

            Material,
            External,
            Internal,

        };

        constexpr
        TensorField(
                Strv &&
                tag,
                Indx
                dim,
                MappingT &&...
                mappings
        )
        :
        tag_(utility::setName(std::forward<Strv>(tag))),
        dim_(dim),
        mappings_({std::forward<Mapping::Type>(mappings)...})
        {}

        constexpr
        Bool
        operator==(
                TensorField const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                TensorField const &
                other
        )
        const = default;

        Name tag_;

        Indx dim_;

        std::array<Mapping::Type, sizeof...(MappingT)> mappings_;

    };

    struct FieldBase
    {

        constexpr
        FieldBase(
                Char
                tag,
                Indx
                dim
        )
        :
        tag_(tag),
        dim_(dim)
        {}

        constexpr
        Bool
        operator==(
                FieldBase const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                FieldBase const &
                other
        )
        const = default;

        Char tag_;

        Indx dim_;

    };

    struct Field : public FieldBase
    {

        enum struct Type
        {

            Material,
            External,
            Internal,

        };

        constexpr
        Field(
                Char
                tag,
                Indx
                dim,
                Type
                typ
        )
        :
        FieldBase(tag, dim),
        typ_(typ)
        {}

        constexpr
        Bool
        operator==(
                Field const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Field const &
                other
        )
        const = default;

        Type typ_;

    };

    template<typename... MappingT>
    requires(std::same_as<MappingT, Mapping::Type> && ...)
    struct Unknown : public FieldBase
    {

        constexpr
        Unknown(
                Char
                tag,
                Indx
                dim,
                MappingT &&...
                mappings
        )
        :
        FieldBase(tag, dim),
        mappings_({std::forward<Mapping::Type>(mappings)...})
        {}

        constexpr
        Bool
        operator==(
                Unknown const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Unknown const &
                other
        )
        const = default;

        std::array<Mapping::Type, sizeof...(MappingT)> mappings_;

    };

}

#endif //LOLITA_LOLITA_HXX
