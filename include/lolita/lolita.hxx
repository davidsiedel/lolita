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

        using Indx = unsigned short;

        using Long = unsigned long long;

        using Real = double;

        using Bool = bool;

    }

    using namespace config;

    namespace utility
    {

        namespace detail
        {

            Indx const static constexpr name_size = 150;

        }

        template<typename T>
        struct Labell2
        {

            constexpr
            Labell2(
                    std::basic_string_view<T> &&
                    str
            )
            noexcept
            {
                for (auto i = 0; i < detail::name_size; ++i) {
                    i < str.size() ? this->operator ()[i] = str[i] : this->operator ()[i] = '#';
                }
            }

            constexpr
            Bool
            operator==(
                    Labell2 const &
                    other
            )
            const = default;

            constexpr
            Bool
            operator!=(
                    Labell2 const &
                    other
            )
            const = default;

            constexpr
            std::basic_string_view<T>
            view()
            const
            {
                return Strv(this->data(), std::distance(std::begin(* this), std::find(std::begin(* this), std::end(* this), '#')));
            }

            Char name_ [detail::name_size];

        };

        template<typename T>
        struct Labell : public std::array<T, detail::name_size>
        {

            constexpr
            Labell(
                    std::basic_string_view<T> &&
                    str
            )
            noexcept
            {
                for (auto i = 0; i < detail::name_size; ++i) {
                    i < str.size() ? this->operator ()[i] = str[i] : this->operator ()[i] = '#';
                }
            }

            constexpr
            Bool
            operator==(
                    Labell const &
                    other
            )
            const = default;

            constexpr
            Bool
            operator!=(
                    Labell const &
                    other
            )
            const = default;

            constexpr
            std::basic_string_view<T>
            view()
            const
            {
                return Strv(this->data(), std::distance(std::begin(* this), std::find(std::begin(* this), std::end(* this), '#')));
            }

        };

        using Label = Labell<Char>;

//        static constexpr
//        Name
//        setName(
//                Strv &&
//                strv
//        )
//        {
//            auto name = Name();
//            for (auto i = 0; i < config::detail::name_size; ++i) {
//                i < strv.size() ? name[i] = strv[i] : name[i] = '#';
//            }
//            return name;
//        }
//
//        static constexpr
//        Strv
//        getName(
//                Name const &
//                name
//        )
//        {
//            return Strv(name.data(), std::distance(std::begin(name), std::find(std::begin(name), std::end(name), '#')));
//        }

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
                Indx
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
                Indx
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
                Indx
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

}

#endif //LOLITA_LOLITA_HXX
