#ifndef DA9C5D5A_13CE_4129_A31C_3550B18DAB24
#define DA9C5D5A_13CE_4129_A31C_3550B18DAB24

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"

namespace lolita2::geometry
{

    template<Field t_field, Domain t_domain>
    struct FieldTraits;

    template<Field t_field, Domain t_domain>
    requires(t_field.isTensor(0))
    struct FieldTraits<t_field, t_domain>
    {

        lolita::matrix::Shape static constexpr shape_ = {
                lolita::numerics::pow(t_domain.dim_, 0),
                lolita::numerics::pow(t_domain.dim_, 0)
        };

        lolita::index static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {1, 1};
        }

        static constexpr
        lolita::integer
        size()
        {
            return shape().size();
        }

    };

    template<Field t_field, Domain t_domain>
    requires(t_field.isTensor(1))
    struct FieldTraits<t_field, t_domain>
    {

        lolita::matrix::Shape static constexpr shape_ = {
                lolita::numerics::pow(t_domain.dim_, 0),
                lolita::numerics::pow(t_domain.dim_, 1)
        };

        lolita::index static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {1, t_domain.dim_};
        }

        static constexpr
        lolita::integer
        size()
        {
            return shape().size();
        }

    };

    template<Field t_field, Domain t_domain>
    requires(t_field.isTensor(2))
    struct FieldTraits<t_field, t_domain>
    {

        lolita::matrix::Shape static constexpr shape_ = {
                lolita::numerics::pow(t_domain.dim_, 1),
                lolita::numerics::pow(t_domain.dim_, 1)
        };

        lolita::index static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {t_domain.dim_, t_domain.dim_};
        }

        static constexpr
        lolita::integer
        size()
        {
            return shape().size();
        }

    };

    template<Field t_field, Domain t_domain>
    requires(t_field.isTensor(3))
    struct FieldTraits<t_field, t_domain>
    {

        lolita::matrix::Shape static constexpr shape_ = {
                lolita::numerics::pow(t_domain.dim_, 1),
                lolita::numerics::pow(t_domain.dim_, 2)
        };

        lolita::index static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {t_domain.dim_, t_domain.dim_ * t_domain.dim_};
        }

        static constexpr
        lolita::integer
        size()
        {
            return shape().size();
        }

    };

    template<Field t_field, Domain t_domain>
    requires(t_field.isTensor(4))
    struct FieldTraits<t_field, t_domain>
    {

        lolita::matrix::Shape static constexpr shape_ = {
                lolita::numerics::pow(t_domain.dim_, 2),
                lolita::numerics::pow(t_domain.dim_, 2)
        };

        lolita::index static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {t_domain.dim_ * t_domain.dim_, t_domain.dim_ * t_domain.dim_};
        }

        static constexpr
        lolita::integer
        size()
        {
            return shape().size();
        }

    };

    struct MappingValues
    {

        constexpr
        lolita::integer
        row()
        const
        {
            return row_;
        }

        constexpr
        lolita::integer
        col()
        const
        {
            return col_;
        }

        constexpr
        lolita::integer
        rank()
        const
        {
            return position_;
        }

        constexpr
        lolita::real
        value()
        const
        {
            return coefficient_;
        }

        lolita::index row_;

        lolita::index col_;

        lolita::index position_;

        lolita::real coefficient_;

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    struct MappingPolicy;

    /*
     * SCALAR
     */

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(0) && t_mapping.isIdentity())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {1, 1};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{1, 1};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(0) && t_domain.hasDim(2) && t_mapping.isGradient())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {2, 1};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{2, 1};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(0) && t_domain.hasDim(3) && t_mapping.isGradient())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {3, 1};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
                MappingValues{0, 2, 2, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{3, 1};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    /*
     * VECTOR
     */

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(1) && t_domain.hasDim(2) && t_mapping.isIdentity())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {2, 1};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{2, 1};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(1) && t_domain.hasDim(3) && t_mapping.isIdentity())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {3, 1};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
                MappingValues{0, 2, 2, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{3, 1};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(1) && t_domain.hasDim(2) && t_mapping.isGradient())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {2, 2};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
                MappingValues{1, 0, 2, 1},
                MappingValues{1, 1, 3, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{2, 2};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(1) && t_domain.hasDim(3) && t_mapping.isGradient())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {3, 3};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
                MappingValues{0, 2, 2, 1},
                MappingValues{1, 0, 3, 1},
                MappingValues{1, 1, 4, 1},
                MappingValues{1, 2, 5, 1},
                MappingValues{2, 0, 6, 1},
                MappingValues{2, 1, 7, 1},
                MappingValues{2, 2, 8, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{3, 3};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(1) && t_domain.hasDim(2) && t_mapping.isSmallStrain())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {1, 4};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{0, 0, 2, 0},
                MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{1, 4};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(1) && t_domain.hasDim(3) && t_mapping.isSmallStrain())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {1, 6};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{2, 2, 2, 0},
                MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
                MappingValues{0, 2, 4, lolita::numerics::sqrt_2},
                MappingValues{1, 2, 5, lolita::numerics::sqrt_2},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{1, 6};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            return lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {}

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(1) && t_domain.hasDim(2) && t_mapping.isLargeStrain())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {1, 5};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{0, 0, 2, 0},
                MappingValues{0, 1, 3, 1},
                MappingValues{1, 0, 4, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{1, 5};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            auto res = lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
            res(0) += 1.0;
            res(1) += 1.0;
            res(2) += 1.0;
            return res;
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {
            gradient(0) += 1.0;
            gradient(1) += 1.0;
            gradient(2) += 1.0;
        }

    };

    template<Field t_field, Domain t_domain, Mapping t_mapping>
    requires(t_field.isTensor(1) && t_domain.hasDim(3) && t_mapping.isLargeStrain())
    struct MappingPolicy<t_field, t_domain, t_mapping>
    {

        lolita::matrix::Shape static constexpr shape_ = {1, 9};

        std::array<MappingValues, shape_.size_> static constexpr values_ = {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{2, 2, 2, 0},
                MappingValues{0, 1, 3, 1},
                MappingValues{0, 2, 4, 1},
                MappingValues{1, 2, 5, 1},
                MappingValues{1, 0, 6, 1},
                MappingValues{2, 0, 7, 1},
                MappingValues{2, 1, 8, 1},
        };

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return lolita::matrix::Shape{1, 9};
        }

        static
        lolita::matrix::Vector<lolita::real, shape_.size_>
        non_linear()
        {
            auto res = lolita::matrix::Zero<lolita::matrix::Vector<lolita::real, shape_.size_>>();
            res(0) += 1.0;
            res(1) += 1.0;
            res(2) += 1.0;
            return res;
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {
            gradient(0) += 1.0;
            gradient(1) += 1.0;
            gradient(2) += 1.0;
        }

    };

}

#endif /* DA9C5D5A_13CE_4129_A31C_3550B18DAB24 */
