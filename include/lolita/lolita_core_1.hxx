//
// Created by dsiedel on 04/06/22.
//

#ifndef LOLITA_LOLITA_CORE_1_HXX
#define LOLITA_LOLITA_CORE_1_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"

namespace lolita::core::field
{

    namespace core_fld = lolita::core::field;

    template<lolita::field::Field _tensor, lolita::index _dim_euclidean>
    struct TensorPolicy;

    template<lolita::field::Field _tensor, lolita::index _dim_euclidean>
    requires(_tensor.isTensor(0))
    struct TensorPolicy<_tensor, _dim_euclidean>
    {



        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 0),
                lolita::numerics::pow(_dim_euclidean, 0)
        };

        lolita::index const static constexpr size_ = shape_.size_;

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

    template<lolita::field::Field _tensor, lolita::index _dim_euclidean>
    requires(_tensor.isTensor(1))
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 0),
                lolita::numerics::pow(_dim_euclidean, 1)
        };

        lolita::index const static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {1, _dim_euclidean};
        }

        static constexpr
        lolita::integer
        size()
        {
            return shape().size();
        }

    };

    template<lolita::field::Field _tensor, lolita::index _dim_euclidean>
    requires(_tensor.isTensor(2))
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 1),
                lolita::numerics::pow(_dim_euclidean, 1)
        };

        lolita::index const static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {_dim_euclidean, _dim_euclidean};
        }

        static constexpr
        lolita::integer
        size()
        {
            return shape().size();
        }

    };

    template<lolita::field::Field _tensor, lolita::index _dim_euclidean>
    requires(_tensor.isTensor(3))
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 1),
                lolita::numerics::pow(_dim_euclidean, 2)
        };

        lolita::index const static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {_dim_euclidean, _dim_euclidean * _dim_euclidean};
        }

        static constexpr
        lolita::integer
        size()
        {
            return shape().size();
        }

    };

    template<lolita::field::Field _tensor, lolita::index _dim_euclidean>
    requires(_tensor.isTensor(4))
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 2),
                lolita::numerics::pow(_dim_euclidean, 2)
        };

        lolita::index const static constexpr size_ = shape_.size_;

        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return {_dim_euclidean * _dim_euclidean, _dim_euclidean * _dim_euclidean};
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    struct MappingPolicy;

    /*
     * SCALAR
     */

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(0) && _mapping.isIdentity())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(0) && _domain.hasDim(2) && _mapping.isGradient())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(0) && _domain.hasDim(3) && _mapping.isGradient())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(1) && _domain.hasDim(2) && _mapping.isIdentity())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(1) && _domain.hasDim(3) && _mapping.isIdentity())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(1) && _domain.hasDim(2) && _mapping.isGradient())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 2};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{1, 0, 2, 1},
                lolita::core::field::MappingValues{1, 1, 3, 1},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(1) && _domain.hasDim(3) && _mapping.isGradient())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 3};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
                lolita::core::field::MappingValues{1, 0, 3, 1},
                lolita::core::field::MappingValues{1, 1, 4, 1},
                lolita::core::field::MappingValues{1, 2, 5, 1},
                lolita::core::field::MappingValues{2, 0, 6, 1},
                lolita::core::field::MappingValues{2, 1, 7, 1},
                lolita::core::field::MappingValues{2, 2, 8, 1},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(1) && _domain.hasDim(2) && _mapping.isSmallStrain())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 4};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{0, 0, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(1) && _domain.hasDim(3) && _mapping.isSmallStrain())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 6};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{2, 2, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
                lolita::core::field::MappingValues{0, 2, 4, lolita::numerics::sqrt_2},
                lolita::core::field::MappingValues{1, 2, 5, lolita::numerics::sqrt_2},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(1) && _domain.hasDim(2) && _mapping.isLargeStrain())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 5};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{0, 0, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, 1},
                lolita::core::field::MappingValues{1, 0, 4, 1},
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

    template<lolita::field::Field _tensor, lolita::domain::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.isTensor(1) && _domain.hasDim(3) && _mapping.isLargeStrain())
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 9};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{2, 2, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, 1},
                lolita::core::field::MappingValues{0, 2, 4, 1},
                lolita::core::field::MappingValues{1, 2, 5, 1},
                lolita::core::field::MappingValues{1, 0, 6, 1},
                lolita::core::field::MappingValues{2, 0, 7, 1},
                lolita::core::field::MappingValues{2, 1, 8, 1},
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

#endif //LOLITA_LOLITA_CORE_1_HXX
