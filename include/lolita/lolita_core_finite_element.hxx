//
// Created by dsiedel on 09/05/22.
//

#ifndef LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
#define LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_element.hxx"
#include "lolita/lolita_core_mesh2.hxx"

namespace lolita::core::field
{

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    struct TensorPolicy;

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 0)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 0),
                lolita::numerics::pow(_dim_euclidean, 0)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 1)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 0),
                lolita::numerics::pow(_dim_euclidean, 1)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 2)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 1),
                lolita::numerics::pow(_dim_euclidean, 1)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 3)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 1),
                lolita::numerics::pow(_dim_euclidean, 2)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 4)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 2),
                lolita::numerics::pow(_dim_euclidean, 2)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    struct MappingValues
    {

        lolita::index row_;

        lolita::index col_;

        lolita::index position_;

        lolita::real coefficient_;

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    struct MappingPolicy;

    /*
     * SCALAR
     */

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 0 && _mapping == lolita::field::Mapping::Identity)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 0 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::Gradient)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 0 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::Gradient)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    /*
     * VECTOR
     */

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::Identity)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::Identity)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::Gradient)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 2};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{1, 0, 2, 1},
                lolita::core::field::MappingValues{1, 1, 3, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::Gradient)
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

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _mapping == lolita::field::Mapping::SmallStrain)
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

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::SmallStrainPlane)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 4};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{0, 0, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::SmallStrainSolid)
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

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::LargeStrainPlane)
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

        static
        void
        non_linear(
//                lolita::matrix::VectorConcept2<shape_.size_> auto & gradient
                auto gradient
        )
        {
            gradient(0) += 1.0;
            gradient(1) += 1.0;
            gradient(2) += 1.0;
        }

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::LargeStrainSolid)
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

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {
            gradient(0) += 1.0;
            gradient(1) += 1.0;
            gradient(2) += 1.0;
        }

    };
    
}

namespace lolita::core::element
{

    /**
     * @brief A simple structure holding constants for the element, and defining the base element type
     * @tparam _element
     * @tparam _domain
     * @tparam _finite_element
     */
    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FiniteElementBase
    {

        /**
         * @brief The finite element object
         */
        lolita::finite_element::FiniteElementConcept auto const static constexpr finite_element_ = _finite_element;

        /**
         * @brief The field object of the finite element argument
         */
        lolita::field::Tensor const static constexpr field_ = _finite_element.unknown_.tensor_;

        /**
         * @brief The array of mappings to which the unknown field is subjected
         */
        std::array<lolita::field::Mapping, _finite_element.unknown_.mappings_.size()> const static constexpr mappings_ = _finite_element.unknown_.mappings_;

        /**
         * @brief The discretization method
         */
        lolita::finite_element::Method const static constexpr method_ = _finite_element.discretization_.method_;

        /**
         * @brief The discretization object
         */
        lolita::finite_element::FiniteElementMethodConcept auto const static constexpr discretization_ = _finite_element.discretization_;

        /**
         * @brief The type of quadrature for the current element
         */
        using Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _finite_element.ord_quadrature_>;

        /**
         * @brief the FieldPolicy given the current domain
         */
        using Field = lolita::core::field::TensorPolicy<field_, _domain.dim_>;

    };

    /**
     * @brief A structure providing a scalar-valued vector of degrees of freedom, and an index to locate the unknown vector in a global framework (e.g.
     * if the unknown is involved at the structural level)
     * @tparam _dim
     */
    template<lolita::integer _dim>
    struct DegreeOfFreedom
    {

    private:

        /**
         * @brief The unknown coefficient vector type
         */
        using _Coefficicents = lolita::matrix::Vector<lolita::real, _dim>;

    public:

        /**
         * @brief initialize the compile-time sized unknown coefficient vector to zero
         */
        void
        setCoefficients()
        requires(_dim > 0)
        {
            coefficients_ = _Coefficicents::Zero();
        }

        /**
         * @brief initialize the run-time sized unknown coefficient vector to zero
         * @param dim_unknown
         */
        void
        setCoefficients(
                lolita::index dim_unknown
        )
        requires(_dim == -1)
        {
            coefficients_ = _Coefficicents::Zero(dim_unknown);
        }

        /**
         * @brief
         * @param mesh
         */
        void
        setCoefficients(
                lolita::numerics::RealConcept auto & mesh
        )
        requires(_dim > 0)
        {
            setCoefficients();
            index_ = mesh;
            mesh += _dim;
        }

        /**
         * @brief
         * @param dim_unknown
         * @param mesh
         */
        void
        setCoefficients(
                lolita::index dim_unknown,
                lolita::numerics::RealConcept auto & mesh
        )
        requires(_dim == -1)
        {
            setCoefficients(dim_unknown);
            index_ = mesh;
            mesh += dim_unknown;
        }

        /**
         * @brief
         */
        _Coefficicents coefficients_;

        /**
         * @brief
         */
        lolita::integer index_;

    };

    /**
     * @brief
     * @tparam _domain
     * @tparam _field
     * @tparam _dim
     */
    template<lolita::geometry::Domain _domain, lolita::field::Tensor _field, lolita::integer _dim>
    struct Unknowns
    {

    private:

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        _dim_binding()
        {
            return _dim == 0 ? 0 : -1;
        }

        /**
         * @brief
         */
        using _Field = lolita::core::field::TensorPolicy<_field, _domain.dim_>;

    public:

        /**
         * @brief
         */
        struct UnknownComponent
        {

            /**
             * @brief
             */
            DegreeOfFreedom<_dim> unknowns_;

            /**
             * @brief
             */
            DegreeOfFreedom<_dim_binding()> bindings_;

        };

        /**
         * @brief
         */
        std::array<std::array<UnknownComponent, _Field::shape_.rows_>, _Field::shape_.cols_> a;

    };

    template<lolita::core::element::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto>
    struct FiniteElementPolicy;

    template<lolita::core::element::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto>
    struct FiniteElementModule;

    template<lolita::core::element::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto>
    struct FiniteElementDiscreteField;

//    template<
//            template<lolita::core::element::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto...> typename,
//            lolita::core::element::Element,
//            lolita::geometry::Domain,
//            lolita::finite_element::FiniteElementConcept auto...
//    >
//    struct FiniteElementConnectivity2;
//
//    template<
//            template<lolita::core::element::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto...> typename _T,
//            lolita::core::element::Element _element,
//            lolita::geometry::Domain _domain,
//            lolita::finite_element::FiniteElementConcept auto... _finite_element
//    >
//    requires(lolita::core::element::PointConcept<_element>)
//    struct FiniteElementConnectivity2<_T, _element, _domain, _finite_element...>
//    {
//
//        lolita::index arg_;
//
//    };

    /**
     * @brief
     * @tparam ...
     */
    template<
            template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename,
            lolita::core::element::Element,
            lolita::geometry::Domain,
            auto...
    >
    struct FiniteElementGeometry;

    /**
     * @brief
     * @tparam _element
     * @tparam _basis
     * @tparam _ord
     */
    template<lolita::core::element::Element _element, lolita::finite_element::Basis _basis, lolita::index _ord>
    struct FiniteElementBasis;

    /**
     * @brief
     * @tparam _element
     * @tparam _basis
     * @tparam _ord
     */
    template<lolita::core::element::Element _element, lolita::finite_element::Basis _basis, lolita::index _ord>
    requires(_basis == lolita::finite_element::Basis::Monomial)
    struct FiniteElementBasis<_element, _basis, _ord>
    {

        /**
         * @brief
         */
        lolita::finite_element::Basis const static constexpr basis_ = _basis;

        /**
         * @brief
         */
        lolita::index const static constexpr ord_ = _ord;

        /**
         * @brief
         */
        lolita::index const static constexpr dim_ = numerics::binomial(_element.dim_ + _ord, _element.dim_);

        /**
         * @brief
         * @tparam _T
         * @tparam _domain
         * @tparam _arg
         */
        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        struct Implementation : public FiniteElementGeometry<_T, _element, _domain, _arg...>
        {

        private:

            /**
             * @brief
             * @return
             */
            static constexpr
            std::array<std::array<lolita::index, 3>, dim_>
            set_exponents()
            {
                auto exponents = std::array<std::array<lolita::index, 3>, dim_>();
                auto row = lolita::index(0);
                if constexpr (_element.dim_ == 0) {
                    exponents[row][0] = 0;
                    exponents[row][1] = 0;
                    exponents[row][2] = 0;
                }
                else if constexpr (_element.dim_ == 1) {
                    for (lolita::index i = 0; i < _ord + 1; ++i) {
                        exponents[row][0] = i;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                        row += 1;
                    }
                }
                else if constexpr (_element.dim_ == 2) {
                    for (lolita::index i = 0; i < _ord + 1; ++i) {
                        for (lolita::index j = 0; j < i + 1; ++j) {
                            exponents[row][0] = i - j;
                            exponents[row][1] = j;
                            exponents[row][2] = 0;
                            row += 1;
                        }
                    }
                }
                else if constexpr (_element.dim_ == 3) {
                    for (lolita::index i = 0; i < _ord + 1; ++i) {
                        for (lolita::index j = 0; j < i + 1; ++j) {
                            for (lolita::index k = 0; k < i + 1; ++k) {
                                if (j + k < i + 1) {
                                    exponents[row][0] = i - (j + k);
                                    exponents[row][1] = k;
                                    exponents[row][2] = j;
                                    row += 1;
                                }
                            }
                        }
                    }
                }
                return exponents;
            }

            /**
             * @brief
             */
            std::array<std::array<lolita::index, 3>, dim_> const static constexpr exponents_ = set_exponents();

        public:

            /**
             * @brief
             * @param point
             * @return
             */
            lolita::matrix::Vector<lolita::real, dim_>
            evaluate(
                    lolita::geometry::Point const & point
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::index i = 0; i < dim_; ++i) {
                    auto value = lolita::real(1);
                    for (lolita::index j = 0; j < _element.dim_; ++j) {
                        auto dist = this->getRiemannianDistance(centroid, point, j);
                        //value *= numerics::pow(2.0 * dist / diameters(j), exponents.get(i, j));
                        value *= std::pow(2.0 * dist / diameters(j), exponents_[i][j]);
                    }
                    basis_vector_values(i) = value;
                }
                return basis_vector_values;
            }

            /**
             * @brief
             * @param point
             * @param derivative_direction
             * @return
             */
            lolita::matrix::Vector<lolita::real, dim_>
            evaluate(
                    lolita::geometry::Point const & point,
                    lolita::index derivative_direction
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::index i = 0; i < dim_; ++i) {
                    auto value = lolita::real(1);
                    for (lolita::index j = 0; j < _element.dim_; ++j) {
                        if (j != derivative_direction) {
                            auto dist = this->getRiemannianDistance(centroid, point, j);
                            value *= std::pow(2.0 * (dist) / diameters(j), exponents_[i][j]);
                        }
                        else {
                            if (exponents_[i][j] > 0) {
                                auto c = 2.0 * exponents_[i][j] / diameters(j);
                                auto dist = this->getRiemannianDistance(centroid, point, j);
                                value *= c * std::pow(2.0 * (dist) / diameters(j), exponents_[i][j] - 1);
//                                value *= c * numerics::pow(2.0 * (dist) / diameters(j), exponents.get(i, j) - 1);
                            }
                            else {
                                value *= 0.0;
                            }
                        }
                    }
                    basis_vector_values(i) = value;
                }
                return basis_vector_values;
            }

        };

    };

    template<
            template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T,
            lolita::core::element::Element _element,
            lolita::geometry::Domain _domain,
            auto... _finite_element
    >
    struct FiniteElementGeometry
    {

    private:

        template<lolita::core::element::Element ___element, lolita::geometry::Domain __domain, auto... __arg>
        using _ElementPointer = std::shared_ptr<_T<___element, __domain, __arg...>>;

        template<lolita::core::element::Element>
        struct _ConnexionPolicy;

        template<lolita::core::element::Element __element>
        requires(!lolita::core::element::PointConcept<__element>)
        struct _ConnexionPolicy<__element>
        {

            using Components = typename lolita::core::element::ElementGeometry<__element>::template Components<_ElementPointer, _domain, _finite_element...>;

            using Neighbours = typename lolita::core::element::ElementGeometry<__element>::template Neighbours<_ElementPointer, _domain, _finite_element...>;

        };

        template<lolita::core::element::Element __element>
        requires(lolita::core::element::PointConcept<__element>)
        struct _ConnexionPolicy<__element>
        {

            using Components = std::shared_ptr<typename lolita::core::element::ElementGeometry<__element>::Components>;

            using Neighbours = typename lolita::core::element::ElementGeometry<__element>::template Neighbours<_ElementPointer, _domain, _finite_element...>;

        };

    public:

        lolita::core::element::Element const static constexpr element_ = _element;

        using Components = typename _ConnexionPolicy<_element>::Components;

        using Neighbours = typename _ConnexionPolicy<_element>::Neighbours;

        lolita::geometry::Point &
        getCurrentCoordinates()
        requires(lolita::core::element::PointConcept<_element>)
        {
            return * components_;
        }

        lolita::geometry::Point const &
        getCurrentCoordinates()
        const
        requires(lolita::core::element::PointConcept<_element>)
        {
            return * components_;
        }

        lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_>
        getCurrentCoordinates()
        const
        requires(!lolita::core::element::PointConcept<_element>)
        {
            auto current_nodes_coordinates = lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_>();
            auto count = lolita::index(0);
            for (auto const & node : getComponents<_element.dim_ - 1, 0>()) {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }

        static
        lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_, matrix::col_major> const>
        getReferenceCoordinates()
        {
            using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_, matrix::col_major> const>;
            return _ReferenceCoordinates(lolita::core::element::ElementGeometry<_element>::reference_nodes_.begin()->begin());
        }

        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::index
        getComponentNodeConnection(
                lolita::index i,
                lolita::index j
        )
        requires(!lolita::core::element::PointConcept<_element>)
        {
            return std::get<_j>(std::get<_i>(lolita::core::element::ElementGeometry<_element>::node_connectivity_))[i][j];
        }

        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Components>> &
        getComponents()
        requires(!lolita::core::element::PointConcept<_element>)
        {
            return std::get<_j>(std::get<_i>(components_));
        }

        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Components>> const &
        getComponents()
        const
        requires(!lolita::core::element::PointConcept<_element>)
        {
            return std::get<_j>(std::get<_i>(components_));
        }

        template<lolita::index _i, lolita::index _j>
        lolita::index
        getComponentIndex(
                lolita::index i
        )
        const
        requires(!lolita::core::element::PointConcept<_element>)
        {
            auto const constexpr _component = lolita::core::element::component<_element, _domain, _i, _j>();
            auto const constexpr _position = lolita::core::element::neighbourPosition<_component, _domain, _element>();
            auto const & items = getComponents<_i, _j>()[i]->template getNeighbours<_position[0], _position[1]>();
            auto is_equal = [&] (
                    std::shared_ptr<FiniteElementGeometry> const & ptr_element
            )
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }

        /**
         * @brief
         * @tparam _basis
         * @tparam _ord
         * @param point
         * @return
         */
        template<lolita::finite_element::Basis _basis, lolita::index _ord>
        lolita::matrix::Vector<lolita::real, FiniteElementBasis<_element, _basis, _ord>::dim_>
        getBasisEvaluation(
                lolita::geometry::Point const & point
        )
        const
        {
            using FiniteElementBasis = typename FiniteElementBasis<_element, _basis, _ord>::template Implementation<_T, _domain, _finite_element...>;
            return static_cast<FiniteElementBasis const *>(this)->evaluate(point);
        }

        /**
         * @brief
         * @tparam _basis
         * @tparam _ord
         * @param point
         * @param derivative_direction
         * @return
         */
        template<lolita::finite_element::Basis _basis, lolita::index _ord>
        lolita::matrix::Vector<lolita::real, FiniteElementBasis<_element, _basis, _ord>::dim_>
        getBasisDerivative(
                lolita::geometry::Point const & point,
                lolita::index derivative_direction
        )
        const
        {
            using FiniteElementBasis = typename FiniteElementBasis<_element, _basis, _ord>::template Implementation<_T, _domain, _finite_element...>;
            return static_cast<FiniteElementBasis const *>(this)->evaluate(point, derivative_direction);
        }

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point
        )
        {
            return lolita::core::element::ElementGeometry<_element>::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @param derivative_direction
         * @return
         */
        static
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point,
                lolita::index derivative_direction
        )
        {
            return lolita::core::element::ElementGeometry<_element>::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }

        /*
         *
         */
        lolita::real
        getShapeMappingDifferential(
                lolita::geometry::Point const & point
        )
        const
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>();
            auto du = lolita::real(0);
            ru.setZero();
            for (lolita::index i = 0; i < _domain.dim_; ++i) {
                for (lolita::index j = 0; j < _element.dim_; ++j) {
                    ru(i, j) = FiniteElementGeometry::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            if constexpr (_element.dim_ == 3) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            else if constexpr (_element.dim_ == 2) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else {
                du = numerics::abs(ru.col(0).norm());
            }
            if constexpr (_domain.frame_ == lolita::geometry::Frame::AxiSymmetric) {
                lolita::real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10) {
                    r0 = 1.e-10;
                }
                du *= 2.0 * numerics::pi * r0;
            }
            return du;
        }

        lolita::real
        getRiemannianDistance(
                lolita::geometry::Point const & first_point_arg,
                lolita::geometry::Point const & second_point_arg,
                lolita::integer direction_arg = -1
        )
        const
        {
            if constexpr (lolita::core::element::CellConcept<_element, _domain>) {
                auto const & current_coordinates = this->getCurrentCoordinates();
                auto distance = lolita::real();
                auto mp0 = lolita::geometry::Point();
                auto mp1 = lolita::geometry::Point();
                for (lolita::index i = 0; i < _element.dim_; ++i) {
                    mp0(i) = FiniteElementGeometry::getShapeMappingEvaluation(current_coordinates.row(i), first_point_arg);
                    mp1(i) = FiniteElementGeometry::getShapeMappingEvaluation(current_coordinates.row(i), second_point_arg);
                }
                direction_arg == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction_arg);
                return distance;
            }
            else {
                using SegmentQuadrature = lolita::core::element::ElementQuadrature<lolita::core::element::seg_02, lolita::finite_element::Quadrature::Gauss, 4>;
                auto distance = lolita::real(0);
                auto dt = lolita::real();
                auto const current_nodes_coordinates = this->getCurrentCoordinates();
                for (lolita::index q = 0; q < SegmentQuadrature::dim_; ++q) {
                    auto pq = SegmentQuadrature::reference_points_[q][0];
                    auto wq = SegmentQuadrature::reference_weights_[q];
                    auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>().setZero();
                    auto difference = second_point_arg - first_point_arg;
                    auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                    for (lolita::index i = 0; i < _domain.dim_; ++i) {
                        for (lolita::index j = 0; j < _element.dim_; ++j) {
                            if (direction_arg == -1 || i == static_cast<lolita::index>(direction_arg)) {
                                auto du = (1.0 / 2.0) * (second_point_arg(j) - first_point_arg(j));
                                auto dx = FiniteElementGeometry::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                                ru(i, j) = dx * du;
                            }
                        }
                    }
                    if constexpr (lolita::core::element::SurfaceConcept<_element>) {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        auto Fff = ru.col(0).template dot(ru.col(1));
                        auto Gff = ru.col(1).template dot(ru.col(1));
                        dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                    }
                    else if constexpr (lolita::core::element::SegmentConcept<_element>) {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        dt = std::sqrt(Eff);
                    }
                    else {
                        dt = 0;
                    }
                    distance += wq * dt;
                }
                return distance;
            }
        }

        /*
         * GEOMETRY
         */

        lolita::geometry::Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElementGeometry::getReferenceCoordinates();
            auto current_diameters = lolita::geometry::Point().setZero();
            for (lolita::index i = 0; i < _element.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < _element.num_nodes_; ++j) {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (lolita::index k = 0; k < 3; ++k) {
                        auto new_value = numerics::abs(getRiemannianDistance(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }

        lolita::geometry::Point
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            return lolita::geometry::getBarycenter(current_nodes_coordinates);
        }

        static
        lolita::geometry::Point
        getReferenceCentroid()
        {
            auto nds = FiniteElementGeometry::getReferenceCoordinates();
            return lolita::geometry::getBarycenter(nds);
        }

        static
        lolita::geometry::Point
        getReferenceDiameters()
        {
            auto dts = lolita::geometry::Point().setZero();//Vector<lolita::real, E.dim>().setZero();
            auto nds = FiniteElementGeometry::getReferenceCoordinates();
            for (lolita::index i = 0; i < _element.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < _element.num_nodes_; ++j) {
                    for (lolita::index k = 0; k < 3; ++k) {
//                        auto & a = dts(k);
//                        auto b = numerics::abs(nds(k, i) - nds(k, j));
//                        if (b > a) {
////                            reference_diameters(k) = b;
//                            a = b;
//                        }
                        auto new_value = numerics::abs(nds(k, i) - nds(k, j));
                        auto & current_value = dts(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return dts;
        }

        /*
         * QUADRATURE
         */

        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        static
        lolita::real
        getReferenceQuadratureWeight(
                lolita::index index_arg
        )
        {
            return ElementQuadrature<_element, _quadrature, _ord>::reference_weights_[index_arg];
        }

        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        static
        lolita::matrix::Span<lolita::geometry::Point const>
        getReferenceQuadraturePoint(
                lolita::index index_arg
        )
        {
//            auto offset = E.dimPoint() * index_arg;
//            //lolita::real const constexpr * const d = ShapeQuadrature<E.shape_description, Q>::reference_points.data.data();
//            auto const constexpr * const d = ElementQuadrature<E, Q>::reference_points.data.data();
//            return matrix::MatMap<Vector<lolita::real, E.dimPoint()> const>(d + offset);
            return lolita::matrix::Span<lolita::geometry::Point const>(ElementQuadrature<_element, _quadrature, _ord>::reference_points_[index_arg].begin());
        }

        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        lolita::real
        getCurrentQuadratureWeight(
                lolita::index index_arg
        )
        const
        {
            auto w = getReferenceQuadratureWeight<_quadrature, _ord>(index_arg);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<_quadrature, _ord>(index_arg));
        }

        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        lolita::geometry::Point
        getCurrentQuadraturePoint(
                lolita::index index_arg
        )
        const
        {
            auto p = lolita::geometry::Point();
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<_quadrature, _ord>(index_arg));
            }
            return p;
        }

        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord, lolita::index _i, lolita::index _j>
        static
        lolita::real
        getComponentReferenceQuadratureWeight(
                lolita::index index_arg
        )
        requires(!lolita::core::element::PointConcept<_element>)
        {
            using ComponentGeometry = FiniteElementGeometry<_T, element::component<_element, _domain, _i, _j>(), _domain, _finite_element...>;
            return ComponentGeometry::template getReferenceQuadratureWeight<_quadrature, _ord>(index_arg);
        }

        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord, lolita::index _i, lolita::index _j>
        static
        lolita::geometry::Point
        getComponentReferenceQuadraturePoint(
                lolita::index component_index_arg,
                lolita::index index_arg
        )
        requires(!lolita::core::element::PointConcept<_element>)
        {
            auto const constexpr _component = lolita::core::element::component<_element, _domain, _i, _j>();
            auto p = lolita::geometry::Point();
            using ComponentGeometry = FiniteElementGeometry<_T, _component, _domain, _finite_element...>;
            auto const & elt_reference_nodes = ElementGeometry<_element>::reference_nodes_;
            for (lolita::index i = 0; i < 3; ++i) {
                auto cpt_coordinates = lolita::matrix::Vector<lolita::real, _component.num_nodes_>();
                for (lolita::index j = 0; j < _component.num_nodes_; ++j) {
                    auto const node_tag = FiniteElementGeometry::template getComponentNodeConnection<_i, _j>(component_index_arg, j);//.get(component_index_arg).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<_quadrature, _ord>(index_arg);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }

        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord, lolita::index _i, lolita::index _j>
        lolita::real
        getComponentCurrentQuadratureWeight(
                lolita::index component_index_arg,
                lolita::index index_arg
        )
        const
        requires(!lolita::core::element::PointConcept<_element>)
        {
            auto const & cmp =  this->template getComponents<_i, _j>()[component_index_arg];//.template get<I>().template get<J>().get(component_index_arg).get();
            return cmp->template getCurrentQuadratureWeight<_quadrature, _ord>(index_arg);
        }

        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord, lolita::index _i, lolita::index _j>
        lolita::geometry::Point
        getComponentCurrentQuadraturePoint(
                lolita::index component_index_arg,
                lolita::index index_arg
        )
        const
        requires(!lolita::core::element::PointConcept<_element>)
        {
            auto p = lolita::geometry::Point();
            auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<_quadrature, _ord, _i, _j>(component_index_arg, index_arg);
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

        lolita::geometry::Point
        getNormalVector(
                lolita::geometry::Point const & point_arg
        )
        const
        requires(lolita::core::element::FaceConcept<_element, _domain>)
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, _element.dim_>();
            ru.setZero();
            for (lolita::index i = 0; i < 3; ++i) {
                for (lolita::index j = 0; j < _element.dim_; ++j) {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, j);
                }
            }
            if constexpr (lolita::core::element::PointConcept<_element>) {
                return lolita::geometry::Point{0, 0, 0};
            }
            else if constexpr (lolita::core::element::CurveConcept<_element>) {
                return lolita::geometry::Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
            }
            else {
                return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            }
        }

        lolita::geometry::Point
        getTangentVector(
                lolita::geometry::Point const & point_arg,
                lolita::index direction_arg
        )
        const
        requires(!lolita::core::element::CellConcept<_element, _domain>)
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = lolita::geometry::Point();
            for (lolita::index i = 0; i < 3; ++i) {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, direction_arg);
            }
            return tangent_vector;
        }

        Neighbours neighbours_;

        Components components_;

        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

        lolita::natural tag_;

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct UnknownsPolicy
    {

        struct _Unknowns
        {

        };

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(_finite_element.discretization_.method_ == lolita::finite_element::Method::HHO && _domain.dim_ - _element.dim_ == 0)
    struct UnknownsPolicy<_element, _domain, _finite_element>
    {

        using _Basis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.discretization_.ord_cell_>;

        using _Unknowns = Unknowns<_domain, _finite_element.unknown_.tensor_, _Basis::dim_>;

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(_finite_element.discretization_.method_ == lolita::finite_element::Method::HHO && _domain.dim_ - _element.dim_ == 1)
    struct UnknownsPolicy<_element, _domain, _finite_element>
    {

    };

    std::shared_ptr<lolita::finite_element::LoadComponent> const static null_load_ptr = std::make_shared<lolita::finite_element::LoadComponent>(lolita::finite_element::LoadComponent());

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FiniteElementPolicy
    {

        static constexpr
        lolita::integer
        _dim_unknown()
        {
            return 0;
        }

        struct UnknownDescription
        {

            lolita::field::Tensor field_;

            lolita::integer dim_;

        };

        struct Module
        {

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

        };

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FiniteElementModule : public FiniteElementDiscreteField<_element, _domain, _finite_element>
    {

    private:

        using _Policy = FiniteElementPolicy<_element, _domain, _finite_element>;

    public:

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        initialize(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            this->setLoads(mesh);
            this->setBehaviour(mesh);
        }

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FiniteElementDiscre : public lolita::core::element::FiniteElementGeometry<lolita::core::element::FiniteElementModule, _element, _domain, _finite_element>
    {

        /**
         * @brief
         */
        using Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _finite_element.ord_quadrature_>;

        /**
         * @brief
         */
        using Field = lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>;

        /**
         * @brief
         * @tparam __finite_element
         * @param mesh
         */
        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setLoads(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            for (auto const & load : mesh.loads_) {
                auto is_equal = [&] (
                        std::shared_ptr<std::basic_string<lolita::character>> const & domain
                )
                {
                    return * domain == load.domain_tag_;
                };
                auto const has_domain = std::find_if(this->domains_.begin(), this->domains_.end(), is_equal) != this->domains_.end();
                auto const has_unknown = load.unknown_tag_ == lolita::utility::readLabel(_finite_element.unknown_.tensor_.tag_);
                if (has_domain && has_unknown) {
                    std::cout << "setting load : " << load.domain_tag_ << std::endl;
                    auto const i = load.components_.row_;
                    auto const j = load.components_.col_;
                    loads_[i][j] = load.load_;
                }
            }
        }

        /**
         * @brief
         * @tparam __finite_element
         * @param mesh
         */
        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setBehaviour(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            for (auto const & behaviour : mesh.behaviours_) {
                auto is_equal = [&] (
                        std::shared_ptr<std::basic_string<lolita::character>> const & domain
                )
                {
                    return * domain == behaviour.domain_tag_;
                };
                auto const has_domain = std::find_if(this->domains_.begin(), this->domains_.end(), is_equal) != this->domains_.end();
                auto const has_unknown = behaviour.unknown_tag_ == lolita::utility::readLabel(_finite_element.unknown_.tensor_.tag_);
                if (has_domain && has_unknown) {
                    std::cout << "setting behaviour : " << behaviour.domain_tag_ << std::endl;
                    behaviour_ = behaviour.behaviour_data_;
                }
            }
        }

        /**
         * @brief
         */
        std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, Field::shape_.rows_>, Field::shape_.cols_> loads_;

        /**
         * @brief
         */
        std::shared_ptr<lolita::behaviour::MgisBehaviourData> behaviour_;

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FiniteElementUU
    {



    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FiniteElementDiscreteField : public lolita::core::element::FiniteElementGeometry<lolita::core::element::FiniteElementModule, _element, _domain, _finite_element>
    {

    private:

        using _Policy = lolita::core::element::FiniteElementPolicy<_element, _domain, _finite_element>;

    public:

        lolita::finite_element::FiniteElementConcept auto const static constexpr finite_element_ = _finite_element;

        lolita::field::UnknownConcept auto const static constexpr unknown_ = _finite_element.unknown_;

        lolita::field::Tensor const static constexpr field_ = _finite_element.unknown_.tensor_;

        lolita::matrix::Shape const static constexpr shp_field_ = lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::shape_;

        struct DegreeOfFreedom
        {

            inline
            lolita::boolean
            isUnknownLocal()
            const
            {
                return unknown_index_ == -1 ? true : false;
            }

            inline
            lolita::boolean
            isUnknownActive()
            const
            {
                return unknown_coefficients_.size() > 0 ? true : false;
            }

            inline
            lolita::boolean
            isBindingLocal()
            const
            {
                return binding_index_ == -1 ? true : false;
            }

            inline
            lolita::boolean
            isBindingActive()
            const
            {
                return binding_coefficients_.size() > 0 ? true : false;
            }

            lolita::matrix::Vector<lolita::real, _Policy::dimUnknown()> unknown_coefficients_;

            lolita::matrix::Vector<lolita::real, _Policy::dimUnknown()> binding_coefficients_;

            lolita::integer unknown_index_ = -1;

            lolita::integer binding_index_ = -1;

            std::shared_ptr<lolita::finite_element::LoadComponent> load_ = null_load_ptr;

        };

        constexpr
        std::basic_string_view<lolita::character>
        getUnknownTag()
        {
            return lolita::utility::readLabel(_finite_element.unknown_.tensor_.tag_);
        }

        void
        setUnknown(
                lolita::index i,
                lolita::index j
        )
        requires(_Policy::dimUnknown() > -1)
        {
            degrees_of_freedom_[i][j].unknown_coefficients_ = lolita::matrix::Vector<lolita::real, _Policy::dimUnknown()>::Zero();
        }

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setUnknown(
                lolita::index i,
                lolita::index j,
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        requires(_Policy::dimUnknown() > -1)
        {
            auto const constexpr _finite_element_index = lolita::core::element::MixedElementPolicy<__finite_element...>::template index<_finite_element>();
            setUnknown(i, j);
            if constexpr(_Policy::dimUnknown() > 0) {
                degrees_of_freedom_[i][j].unknown_index_ = mesh.dof_indices_[_finite_element_index].num_unknowns_;
                mesh.dof_indices_[_finite_element_index].num_unknowns_ += _Policy::dimUnknown();
            }
        }

        void
        setUnknown(
                lolita::index i,
                lolita::index j,
                lolita::index dim_unknowns
        )
        requires(_Policy::dimUnknown() == -1)
        {
            degrees_of_freedom_[i][j].unknown_coefficients_ = lolita::matrix::Vector<lolita::real, _Policy::dimUnknown()>::Zero(dim_unknowns);
        }

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setUnknown(
                lolita::index i,
                lolita::index j,
                lolita::index dim_unknowns,
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        requires(_Policy::dimUnknown() == -1)
        {
            auto const constexpr _finite_element_index = lolita::core::element::MixedElementPolicy<__finite_element...>::template index<_finite_element>();
            setUnknown(i, j, dim_unknowns);
            if (dim_unknowns > 0) {
                degrees_of_freedom_[i][j].unknown_index_ = mesh.dof_indices_[_finite_element_index].num_unknowns_;
                mesh.dof_indices_[_finite_element_index].num_unknowns_ += dim_unknowns;
            }
        }

        void
        setBinding(
                lolita::index i,
                lolita::index j
        )
        requires(_Policy::dimUnknown() > -1)
        {
            degrees_of_freedom_[i][j].binding_coefficients_ = lolita::matrix::Vector<lolita::real, _Policy::dimUnknown()>::Zero();
        }

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setBinding(
                lolita::index i,
                lolita::index j,
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        requires(_Policy::dimUnknown() > -1)
        {
            auto const constexpr _finite_element_index = lolita::core::element::MixedElementPolicy<__finite_element...>::template index<_finite_element>();
            setBinding(i, j);
            if constexpr(_Policy::dimUnknown() > 0) {
                degrees_of_freedom_[i][j].binding_index_ = mesh.dof_indices_[_finite_element_index].num_bindings_;
                mesh.dof_indices_[_finite_element_index].num_bindings_ += _Policy::dimUnknown();
            }
        }

        void
        setBinding(
                lolita::index i,
                lolita::index j,
                lolita::index dim_unknowns
        )
        requires(_Policy::dimUnknown() == -1)
        {
            degrees_of_freedom_[i][j].binding_coefficients_ = lolita::matrix::Vector<lolita::real, _Policy::dimUnknown()>::Zero(dim_unknowns);
        }

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setBinding(
                lolita::index i,
                lolita::index j,
                lolita::index dim_unknowns,
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        requires(_Policy::dimUnknown() == -1)
        {
            auto const constexpr _finite_element_index = lolita::core::element::MixedElementPolicy<__finite_element...>::template index<_finite_element>();
            setBinding(i, j, dim_unknowns);
            if (dim_unknowns > 0) {
                degrees_of_freedom_[i][j].binding_index_ = mesh.dof_indices_[_finite_element_index].num_bindings_;
                mesh.dof_indices_[_finite_element_index].num_bindings_ += dim_unknowns;
            }
        }

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setLoads(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            for (auto const & load : mesh.loads_) {
                auto is_equal = [&] (
                        std::shared_ptr<std::basic_string<lolita::character>> const & domain
                )
                {
                    return * domain == load.domain_tag_;
                };
                auto const has_domain = std::find_if(this->domains_.begin(), this->domains_.end(), is_equal) != this->domains_.end();
                auto const has_unknown = load.unknown_tag_ == lolita::utility::readLabel(_finite_element.unknown_.tensor_.tag_);
                if (has_domain && has_unknown) {
                    std::cout << "setting load : " << load.domain_tag_ << std::endl;
                    auto const i = load.components_.row_;
                    auto const j = load.components_.col_;
                    degrees_of_freedom_[i][j].load_ = load.load_;
                    if (load.load_->loading_ == lolita::finite_element::Loading::Constraint) {
                        if constexpr (_Policy::dimUnknown() > -1) {
                            setBinding(i, j, mesh);
                        }
                        else {
                            setBinding(i, j, degrees_of_freedom_[i][j].unknown_coefficients_.size(), mesh);
                        }
                    }
                }
            }
        }

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setBehaviour(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            for (auto const & behaviour : mesh.behaviours_) {
                auto is_equal = [&] (
                        std::shared_ptr<std::basic_string<lolita::character>> const & domain
                )
                {
                    return * domain == behaviour.domain_tag_;
                };
                auto const has_domain = std::find_if(this->domains_.begin(), this->domains_.end(), is_equal) != this->domains_.end();
                auto const has_unknown = behaviour.unknown_tag_ == lolita::utility::readLabel(_finite_element.unknown_.tensor_.tag_);
                if (has_domain && has_unknown) {
                    std::cout << "setting behaviour : " << behaviour.domain_tag_ << std::endl;
                    behaviour_ = behaviour.behaviour_data_;
                }
            }
        }

        std::array<std::array<DegreeOfFreedom, shp_field_.rows_>, shp_field_.cols_> degrees_of_freedom_;

        std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, shp_field_.rows_>, shp_field_.cols_> loads__;

        std::shared_ptr<lolita::behaviour::MgisBehaviourData> behaviour_;

    };

    /*
     * CELL
     */

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(lolita::core::element::CellConcept<_element, _domain> && _finite_element.discretization_.method_ == lolita::finite_element::Method::HHO)
    struct FiniteElementPolicy<_element, _domain, _finite_element>
    {

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        dimUnknown()
        {
            return lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.discretization_.ord_cell_>::dim_;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        numUnknown()
        {
            return dimUnknown() * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        dimElementUnknown()
        {
            using _CellBasis = FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.discretization_.ord_cell_>;
            auto dim_cell_element_unknown = _CellBasis::dim_;
            auto set_dim_field_mapping = [&] <lolita::index _k = 0> (
                    auto & self
            )
            constexpr mutable
            {
                auto const constexpr _component = lolita::core::element::component<_element, _domain, 0, _k>();
                auto const constexpr _num_components = lolita::core::element::numComponents<_element, _domain, 0, _k>();
                using _FaceBasis = FiniteElementBasis<_component, lolita::finite_element::Basis::Monomial, _finite_element.discretization_.ord_face_>;
                dim_cell_element_unknown += _FaceBasis::dim_ * _num_components;
                if constexpr (_k < element::numComponents<_element, _domain, 0>() - 1) {
                    self.template operator ()<_k + 1>(self);
                }
            };
            set_dim_field_mapping(set_dim_field_mapping);
            return dim_cell_element_unknown;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        numElementUnknown()
        {
            return dimElementUnknown() * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_;
        }

        struct Module
        {

            lolita::matrix::Matrix<lolita::real, numElementUnknown(), numElementUnknown()> stabilization_;

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

            template<lolita::field::Mapping _mapping>
            struct GradientOperatorFactory
            {

                auto const static constexpr ord_cell_ = _finite_element.discretization_.ord_cell_;
                auto const static constexpr ord_face_ = _finite_element.discretization_.ord_face_;
                auto const static constexpr ord_grad_ = _finite_element.discretization_.ordMapping(_mapping);
                auto const static constexpr ord_max_ = lolita::numerics::max(ord_cell_, ord_face_, ord_grad_);
                auto const static constexpr ord_qad_ = _domain.ordIntegration(ord_max_);
                using _Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, ord_qad_>;
                using _CellBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_cell_>;
                using _GradBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_grad_>;

            };

            template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
            void
            setUnknowns(
                    lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
            )
            {
                for (int i = 0; i < this->shp_field_.rows_; ++i) {
                    for (int j = 0; j < this->shp_field_.cols_; ++j) {
                        this->setUnknown(i, j);
                    }
                }
            }

            lolita::matrix::Vector<lolita::real, numElementUnknown()>
            getElementUnknown()
            const
            {
                auto element_unknown = lolita::matrix::Vector<lolita::real, numElementUnknown()>();
                auto count = lolita::index(0);
                for (int i = 0; i < this->shp_field_.rows_; ++i) {
                    for (int j = 0; j < this->shp_field_.cols_; ++j) {
                        element_unknown.template segment<dimUnknown()>(count) = this->degrees_of_freedom_[i][j].unknown_coefficients_;
                        count += dimUnknown();
                    }
                }
                auto set_faces_unknowns = [&] <lolita::index _i = 0u> (
                        auto & self
                )
                mutable
                {
                    auto const constexpr _component = lolita::core::element::component<_element, _domain, 0, _i>();
                    using _FacePolicy = lolita::core::element::FiniteElementPolicy<_component, _domain, _finite_element>;
                    for (auto const & face : this->template getComponents<0, _i>()) {
                        for (int i = 0; i < this->shp_field_.rows_; ++i) {
                            for (int j = 0; j < this->shp_field_.cols_; ++j) {
                                element_unknown.template segment<_FacePolicy::dimUnknown()>(count) = face->degrees_of_freedom_[i][j].unknown_coefficients_;
                                count += _FacePolicy::dimUnknown();
                            }
                        }
                        if constexpr (_i < lolita::core::element::numComponents<_element, _domain, 0>() - 1) {
                            self.template operator()<_i + 1u>(self);
                        }
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
                return element_unknown;
            }

            void
            setModule()
            {
                this->module_.stabilization_.setOnes();
                this->module_.stabilization_(0, 0) = 3;
            }

            template<lolita::field::Mapping _mapping>
            void
            setMapping();

            template<lolita::field::Mapping _mapping>
            void
            setMapping()
            requires(_mapping == lolita::field::Mapping::Identity)
            {}

            template<lolita::field::Mapping _mapping>
            void
            setMapping()
            requires(_mapping == lolita::field::Mapping::Gradient)
            {
                /*
                 * Defining constants
                 */
                auto const constexpr ord_cell = _finite_element.discretization_.ord_cell_;
                auto const constexpr ord_face = _finite_element.discretization_.ord_face_;
                auto const constexpr ord_grad = _finite_element.discretization_.ordMapping(_mapping);
                auto const constexpr ord_max = lolita::numerics::max(ord_cell, ord_face, ord_grad);
                auto const constexpr ord_qad = _domain.ordIntegration(ord_max);
                using _Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, ord_qad>;
                using _Quadrature2 = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _finite_element.ord_quadrature_>;
                using _CellBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_cell>;
                using _GradBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_grad>;
                /*
                 *
                 */
                using _MapPol = lolita::core::field::MappingPolicy<_finite_element.unknown_.tensor_, _domain, _mapping>;
                /*
                 * LHS
                 */
                auto get_lhs = [&] ()
                mutable
                {
                    auto lhs = lolita::matrix::Matrix<lolita::real, _GradBasis::dim_, _GradBasis::dim_>().setZero();
                    for (int i = 0; i < _Quadrature::dim_; ++i) {
                        auto pt = this->template getReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                        auto wt = this->template getCurrentQuadratureWeight<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                        auto vr = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pt);
                        lhs += wt * vr * vr.transpose();
                    }
                    lhs = lhs.llt().solve(decltype(lhs)::Identity());
                    return lhs;
                };
                /*
                 * RHS
                 */
                auto get_rhs = [&] (
                        auto
                        i_arg,
                        auto
                        j_arg
                )
                mutable
                {
                    /*
                     * Initializing the RHS part of the gradient operator
                     */
                    auto rhs = lolita::matrix::Matrix<lolita::real, _GradBasis::dim_, numElementUnknown()>().setZero();
                    /*
                     * Setting the cell part of the gradient operator
                     */
                    auto set_rhs_cell = [&] ()
                    mutable
                    {
                        for (int i = 0; i < _Quadrature::dim_; ++i) {
                            auto pc = this->template getReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                            auto wc = this->template getCurrentQuadratureWeight<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                            auto vr = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pc);
                            auto vc = this->template getBasisDerivative<_CellBasis::basis_, _CellBasis::ord_>(pc, j_arg);
                            auto prt = rhs.template block<_GradBasis::dim_, _CellBasis::dim_>(0, _CellBasis::dim_ * i_arg);
                            prt += wc * vr * vc.transpose();
                        }
                    };
                    /*
                     * Defining the face offset
                     */
                    auto get_faces_offset2 = [&] <lolita::index K> ()
                            constexpr
                    {
                        auto offset = _CellBasis::dim_ * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_;
                        auto get_faces_offset22 = [&] <auto L = 0> (
                                auto & self
                        )
                        constexpr mutable
                        {
                            if constexpr (L > 0) {
                                auto const constexpr _component = lolita::core::element::component<_element, _domain, 0, L - 1>();
                                using _FaceBasis = lolita::core::element::FiniteElementBasis<_component, lolita::finite_element::Basis::Monomial, ord_face>;
                                auto nj = _FaceBasis::dim_ * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_;
                                auto num_c = lolita::core::element::numComponents<_element, _domain, 0, L - 1>();
                                offset += num_c * nj;
                            }
                            if constexpr (L < K) {
                                self.template operator()<K + 1>(self);
                            }
                        };
                        get_faces_offset22(get_faces_offset22);
                        return offset;
                    };
                    /*
                     * Setting the jump part of the gradient operator, at a given face
                     */
                    auto set_rhs_face = [&] <lolita::index K = 0> (
                            auto &
                            self
                    )
                    mutable
                    {
                        auto const constexpr _component = lolita::core::element::component<_element, _domain, 0, K>();
                        using _FaceBasis = lolita::core::element::FiniteElementBasis<_component, lolita::finite_element::Basis::Monomial, ord_face>;
                        auto count = lolita::index(0);
                        for (auto & face : this->template getComponents<0, K>()) {
                            auto n_dir = lolita::real(1);
                            for (int i = 0; i < _Quadrature::dim_; ++i) {
                                auto pf = face->template getReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                                auto pc = this->template getComponentReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_, 0, K>(count, i);
                                auto wf = face->template getCurrentQuadratureWeight<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                                auto n = face->getNormalVector(pf);
                                auto vr = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pc);
                                auto vc = this->template getBasisEvaluation<_CellBasis::basis_, _CellBasis::ord_>(pc);
                                auto vf = face->template getBasisEvaluation<_FaceBasis::basis_, _FaceBasis::ord_>(pf);
                                auto offset = _CellBasis::dim_ * i_arg;
                                auto prt_cell = rhs.template block<_GradBasis::dim_, _CellBasis::dim_>(0, offset);
                                prt_cell -= wf * vr * vc.transpose() * n(j_arg) * n_dir;
                                //offset = get_faces_offset.template operator()<K>();
                                offset = get_faces_offset2.template operator()<K>();
                                offset += _FaceBasis::dim_ * (count * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_ + i_arg);
                                auto prt_face = rhs.template block<_GradBasis::dim_, _FaceBasis::dim_>(0, offset);
                                prt_face += wf * vr * vf.transpose() * n(j_arg) * n_dir;
                            }
                            count ++;
                        }
                        if constexpr (K < numComponents<_element, _domain, 0>() - 1) {
                            self.template operator()<K + 1>(i_arg, j_arg, self);
                        }
                    };
                    set_rhs_cell();
                    set_rhs_face(set_rhs_face);
                    return rhs;
                };
                /*
                 *
                 */
                auto lhs = get_lhs();
                for (lolita::core::field::MappingValues const & item : _MapPol::values_) {
//                    auto rhs = get_rhs(item.col_, item.row_);
                    auto rhs = get_rhs(item.row_, item.col_);
                    for (int k = 0; k < _Quadrature::dim_; ++k) {
                        auto pnt = this->template getReferenceQuadraturePoint<_Quadrature2::quadrature_, _Quadrature2::ord_>(k);
                        auto vct = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pnt);
                        this->operators_[k].row(this->template rowMapping<_mapping>() + item.position_) = vct.transpose() * item.coefficient_ * lhs * rhs;
//                        std::cout << "operator : " << k << std::endl;
//                        std::cout << this->operators_[k] << std::endl;
                    }
                }
            }

            template<lolita::field::Mapping _mapping>
            void
            setMapping()
            requires(_mapping == lolita::field::Mapping::LargeStrainPlane)
            {
                /*
                 * Defining constants
                 */
                auto const constexpr ord_cell = _finite_element.discretization_.ord_cell_;
                auto const constexpr ord_face = _finite_element.discretization_.ord_face_;
                auto const constexpr ord_grad = _finite_element.discretization_.ordMapping(_mapping);
                auto const constexpr ord_max = lolita::numerics::max(ord_cell, ord_face, ord_grad);
                auto const constexpr ord_qad = _domain.ordIntegration(ord_max);
                using _Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, ord_qad>;
                using _Quadrature2 = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _finite_element.ord_quadrature_>;
                using _CellBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_cell>;
                using _GradBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_grad>;
                /*
                 *
                 */
                using _MapPol = lolita::core::field::MappingPolicy<_finite_element.unknown_.tensor_, _domain, _mapping>;
                /*
                 * LHS
                 */
                auto get_lhs = [&] ()
                        mutable
                {
                    auto lhs = lolita::matrix::Matrix<lolita::real, _GradBasis::dim_, _GradBasis::dim_>().setZero();
                    for (int i = 0; i < _Quadrature::dim_; ++i) {
                        auto pt = this->template getReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                        auto wt = this->template getCurrentQuadratureWeight<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                        auto vr = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pt);
                        lhs += wt * vr * vr.transpose();
                    }
                    lhs = lhs.llt().solve(decltype(lhs)::Identity());
                    return lhs;
                };
                /*
                 * RHS
                 */
                auto get_rhs = [&] (
                        auto
                        i_arg,
                        auto
                        j_arg
                )
                        mutable
                {
                    /*
                     * Initializing the RHS part of the gradient operator
                     */
                    auto rhs = lolita::matrix::Matrix<lolita::real, _GradBasis::dim_, numElementUnknown()>().setZero();
                    /*
                     * Setting the cell part of the gradient operator
                     */
                    auto set_rhs_cell = [&] ()
                            mutable
                    {
                        for (int i = 0; i < _Quadrature::dim_; ++i) {
                            auto pc = this->template getReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                            auto wc = this->template getCurrentQuadratureWeight<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                            auto vr = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pc);
                            auto vc = this->template getBasisDerivative<_CellBasis::basis_, _CellBasis::ord_>(pc, j_arg);
                            auto prt = rhs.template block<_GradBasis::dim_, _CellBasis::dim_>(0, _CellBasis::dim_ * i_arg);
                            prt += wc * vr * vc.transpose();
                        }
                    };
                    /*
                     * Defining the face offset
                     */
                    auto get_faces_offset2 = [&] <lolita::index K> ()
                            constexpr
                    {
                        auto offset = _CellBasis::dim_ * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_;
                        auto get_faces_offset22 = [&] <auto L = 0> (
                                auto & self
                        )
                                constexpr mutable
                        {
                            if constexpr (L > 0) {
                                auto const constexpr _component = lolita::core::element::component<_element, _domain, 0, L - 1>();
                                using _FaceBasis = lolita::core::element::FiniteElementBasis<_component, lolita::finite_element::Basis::Monomial, ord_face>;
                                auto nj = _FaceBasis::dim_ * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_;
                                auto num_c = lolita::core::element::numComponents<_element, _domain, 0, L - 1>();
                                offset += num_c * nj;
                            }
                            if constexpr (L < K) {
                                self.template operator()<K + 1>(self);
                            }
                        };
                        get_faces_offset22(get_faces_offset22);
                        return offset;
                    };
                    /*
                     * Setting the jump part of the gradient operator, at a given face
                     */
                    auto set_rhs_face = [&] <lolita::index K = 0> (
                            auto &
                            self
                    )
                            mutable
                    {
                        auto const constexpr _component = lolita::core::element::component<_element, _domain, 0, K>();
                        using _FaceBasis = lolita::core::element::FiniteElementBasis<_component, lolita::finite_element::Basis::Monomial, ord_face>;
                        auto count = lolita::index(0);
                        for (auto & face : this->template getComponents<0, K>()) {
                            auto n_dir = lolita::real(1);
                            for (int i = 0; i < _Quadrature::dim_; ++i) {
                                auto pf = face->template getReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                                auto pc = this->template getComponentReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_, 0, K>(count, i);
                                auto wf = face->template getCurrentQuadratureWeight<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                                auto n = face->getNormalVector(pf);
                                auto vr = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pc);
                                auto vc = this->template getBasisEvaluation<_CellBasis::basis_, _CellBasis::ord_>(pc);
                                auto vf = face->template getBasisEvaluation<_FaceBasis::basis_, _FaceBasis::ord_>(pf);
                                auto offset = _CellBasis::dim_ * i_arg;
                                auto prt_cell = rhs.template block<_GradBasis::dim_, _CellBasis::dim_>(0, offset);
                                prt_cell -= wf * vr * vc.transpose() * n(j_arg) * n_dir;
                                //offset = get_faces_offset.template operator()<K>();
                                offset = get_faces_offset2.template operator()<K>();
                                offset += _FaceBasis::dim_ * (count * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_ + i_arg);
                                auto prt_face = rhs.template block<_GradBasis::dim_, _FaceBasis::dim_>(0, offset);
                                prt_face += wf * vr * vf.transpose() * n(j_arg) * n_dir;
                            }
                            count ++;
                        }
                        if constexpr (K < numComponents<_element, _domain, 0>() - 1) {
                            self.template operator()<K + 1>(i_arg, j_arg, self);
                        }
                    };
                    set_rhs_cell();
                    set_rhs_face(set_rhs_face);
                    return rhs;
                };
                /*
                 *
                 */
                auto lhs = get_lhs();
                for (lolita::core::field::MappingValues const & item : _MapPol::values_) {
//                    auto rhs = get_rhs(item.col_, item.row_);
                    auto rhs = get_rhs(item.row_, item.col_);
                    for (int k = 0; k < _Quadrature::dim_; ++k) {
                        auto pnt = this->template getReferenceQuadraturePoint<_Quadrature2::quadrature_, _Quadrature2::ord_>(k);
                        auto vct = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pnt);
                        this->operators_[k].row(this->template rowMapping<_mapping>() + item.position_) = vct.transpose() * item.coefficient_ * lhs * rhs;
//                        std::cout << "operator : " << k << std::endl;
//                        std::cout << this->operators_[k] << std::endl;
                    }
                }
            }

        };

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(lolita::core::element::CellConcept<_element, _domain>)
    struct FiniteElementModule<_element, _domain, _finite_element> : public lolita::core::element::FiniteElementDiscreteField<_element, _domain, _finite_element>
    {

    private:

        using _Policy = FiniteElementPolicy<_element, _domain, _finite_element>;

        using _Module = typename _Policy::Module;

        using _Implementation = typename _Policy::Implementation;

    public:

        template<lolita::field::Mapping _mapping>
        static constexpr
        lolita::index
        dimMapping()
        {
            return lolita::core::field::MappingPolicy<_finite_element.unknown_.tensor_, _domain, _mapping>::shape_.size_;
        }

        template<lolita::field::Mapping _mapping>
        static constexpr
        lolita::index
        rowMapping()
        {
            auto value = lolita::index(0);
            for (auto m : _finite_element.unknown_.mappings_) {
                if (m == _mapping) {
                    return value;
                }
                value += dimMapping<_mapping>();
            }
            return value;
        }

        static constexpr
        lolita::index
        dimOperator()
        {
            auto dim_mapping = lolita::index(0);
            auto set_dim_mapping = [&] <lolita::index _i = 0> (auto & self)
            constexpr mutable
            {
                dim_mapping += dimMapping<_finite_element.unknown_.mappings_[_i]>();
                if constexpr (_i < _finite_element.unknown_.mappings_.size() - 1) {
                    self.template operator()<_i + 1>(self);
                }
            };
            set_dim_mapping(set_dim_mapping);
            return dim_mapping;
        }

        using Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _finite_element.ord_quadrature_>;

        using Operator = lolita::matrix::Matrix<lolita::real, dimOperator(), _Policy::numElementUnknown()>;

        using BilinearFormOperator = lolita::matrix::Matrix<lolita::real, dimOperator(), _Policy::numElementUnknown()>;

        using LinearFormOperator = lolita::matrix::Vector<lolita::real, _Policy::numElementUnknown()>;

        struct MaterialPoint
        {

            BilinearFormOperator bilinear_form_operator_;

            LinearFormOperator linear_form_operator_;

            lolita::geometry::Point point_;

            lolita::real weight_;

            std::unique_ptr<mgis::behaviour::BehaviourData> material_point_;

        };

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        initialize(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            this->setLoads(mesh);
            this->setBehaviour(mesh);
            static_cast<_Implementation *>(this)->setUnknowns(mesh);
            static_cast<_Implementation *>(this)->setModule();
            auto set_operator = [&] <lolita::index _i = 0> (
                    auto & self
            )
            mutable
            {
                static_cast<_Implementation *>(this)->template setMapping<_finite_element.unknown_.mappings_[_i]>();
                if constexpr (_i < _finite_element.unknown_.mappings_.size() - 1) {
                    self.template operator()<_i + 1>(self);
                }
            };
            set_operator(set_operator);
            for (int i = 0; i < Quadrature::dim_; ++i) {
                material_points_[i] = std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(this->behaviour_->behaviour_));
                material_points_[i]->K[0] = 4;
            }
        }

        lolita::matrix::Vector<lolita::real, _Policy::numElementUnknown()>
        getElementUnknown()
        const
        {
            return static_cast<_Implementation const *>(this)->getElementUnknown();
        }

        static
        void
        ttest(
                lolita::matrix::MatrixConceptE<lolita::real, 2, 2> auto matrix
        )
        {

        }

        void
        integrateBehaviour()
        {
            auto u_vec = getElementUnknown();
            for (int i = 0; i < Quadrature::dim_; ++i) {
                auto grad_vec = lolita::matrix::Vector<lolita::real, dimOperator()>();
                grad_vec = operators_[i] * u_vec;
                std::cout << "grad_vec : " << std::endl;
                std::cout << grad_vec << std::endl;
                auto comp_grad = [&] <lolita::index _i = 0u> (
                        auto & self
                )
                mutable
                {
                    auto const constexpr _mapping = _finite_element.unknown_.mappings_[_i];
                    using _MappingPolicy = lolita::core::field::MappingPolicy<_finite_element.unknown_.tensor_, _domain, _mapping>;
                    ttest(lolita::matrix::Matrix<lolita::real, 2, 2>::Zero());
//                    static_assert(lolita::matrix::MatrixConcept<decltype(grad_vec)>);
//                    static_assert(lolita::matrix::MatrixConcept<decltype(grad_vec.template segment<_MappingPolicy::shape_.size_>(rowMapping<_mapping>()))>);
//                    static_assert(lolita::matrix::MatrixConceptE<decltype(grad_vec.template segment<_MappingPolicy::shape_.size_>(rowMapping<_mapping>())), lolita::real, _MappingPolicy::shape_.size_, 1>);
//                    static_assert(lolita::matrix::size(grad_vec.template segment<_MappingPolicy::shape_.size_>(rowMapping<_mapping>())) == _MappingPolicy::shape_.size_);
                    _MappingPolicy::non_linear(grad_vec.template segment<_MappingPolicy::shape_.size_>(rowMapping<_mapping>()));
//                    _MappingPolicy::non_linear(lolita::matrix::Span<lolita::matrix::Vector<lolita::real, _MappingPolicy::shape_.size_>>(grad_vec.data() + rowMapping<_mapping>()));
                    std::cout << "grad_vec : " << std::endl;
                    std::cout << grad_vec << std::endl;
                    if constexpr (_i < _finite_element.unknown_.mappings_.size() - 1) {
                        self.template operator()<_i + 1>(self);
                    }
                };
                comp_grad(comp_grad);
                auto grad_view = lolita::matrix::Span<lolita::matrix::Vector<lolita::real, dimOperator()>>(material_points_[i]->s1.gradients.data());
                grad_view = grad_vec;
                auto behaviour_view = mgis::behaviour::make_view(* material_points_[i]);
                this->template getCurrentQuadraturePoint<Quadrature::quadrature_, Quadrature::dim_>(i);
                auto result = mgis::behaviour::integrate(behaviour_view, this->behaviour_->behaviour_);
            }
        }

        _Module module_;

        std::array<BilinearFormOperator, Quadrature::dim_> operators_;

        std::array<LinearFormOperator, Quadrature::dim_> linear_form_operators_;

        std::array<MaterialPoint, Quadrature::dim_> mats_;

        std::array<std::unique_ptr<mgis::behaviour::BehaviourData>, Quadrature::dim_> material_points_;

    };

    /*
     * FACE
     */

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(lolita::core::element::FaceConcept<_element, _domain> && _finite_element.discretization_.method_ == lolita::finite_element::Method::HHO)
    struct FiniteElementPolicy<_element, _domain, _finite_element>
    {

        /**
         * @brief the number of unknown coefficient per component. -1 is not known at compile time, 0 is empty, and > 0 is known at compile time
         * @return
         */
        static constexpr
        lolita::integer
        dimUnknown()
        {
            return FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.discretization_.ord_face_>::dim_;
        }

        /**
         * @brief the total number of unknowns from all components and neighbours contribution.
         * @return
         */
        static constexpr
        lolita::integer
        numElementUnknown()
        {
            return dimUnknown() * lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>::dim_;
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

            template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
            void
            setUnknowns(
                    lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
            )
            {
                for (int i = 0; i < this->shp_field_.rows_; ++i) {
                    for (int j = 0; j < this->shp_field_.cols_; ++j) {
                        this->setUnknown(i, j, mesh);
                    }
                }
            }
        };

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(lolita::core::element::FaceConcept<_element, _domain>)
    struct FiniteElementModule<_element, _domain, _finite_element> : public lolita::core::element::FiniteElementDiscreteField<_element, _domain, _finite_element>
    {

    private:

        using _Policy = FiniteElementPolicy<_element, _domain, _finite_element>;

        using _Module = typename _Policy::Module;

        using _Implementation = typename _Policy::Implementation;

    public:

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        initialize(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            this->setLoads(mesh);
            this->setBehaviour(mesh);
            static_cast<_Implementation *>(this)->setUnknowns(mesh);
        }

    };

    /*
     * EDGE
     */

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(lolita::core::element::EdgeConcept<_element, _domain> && _finite_element.discretization_.method_ == lolita::finite_element::Method::HHO)
    struct FiniteElementPolicy<_element, _domain, _finite_element>
    {

        static constexpr
        lolita::integer
        dimUnknown()
        {
            return 0;
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

            template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
            void
            setUnknowns(
                    lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
            )
            {

            }

        };

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(lolita::core::element::EdgeConcept<_element, _domain>)
    struct FiniteElementModule<_element, _domain, _finite_element> : public lolita::core::element::FiniteElementDiscreteField<_element, _domain, _finite_element>
    {

    private:

        using _Policy = FiniteElementPolicy<_element, _domain, _finite_element>;

        using _Module = typename _Policy::Module;

        using _Implementation = typename _Policy::Implementation;

    public:

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        initialize(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            this->setLoads(mesh);
            this->setBehaviour(mesh);
            static_cast<_Implementation *>(this)->setUnknowns(mesh);
        }

    };

    /*
     * NODE
     */

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(lolita::core::element::NodeConcept<_element, _domain> && _finite_element.discretization_.method_ == lolita::finite_element::Method::HHO)
    struct FiniteElementPolicy<_element, _domain, _finite_element>
    {

        static constexpr
        lolita::integer
        dimUnknown()
        {
            return 0;
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

            template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
            void
            setUnknowns(
                    lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
            )
            {

            }

        };

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(lolita::core::element::NodeConcept<_element, _domain>)
    struct FiniteElementModule<_element, _domain, _finite_element> : public lolita::core::element::FiniteElementDiscreteField<_element, _domain, _finite_element>
    {

    private:

        using _Policy = FiniteElementPolicy<_element, _domain, _finite_element>;

        using _Module = typename _Policy::Module;

        using _Implementation = typename _Policy::Implementation;

    public:

        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        initialize(
                lolita::core::mesh::MeshData2<_domain, __finite_element...> & mesh
        )
        {
            this->setLoads(mesh);
            this->setBehaviour(mesh);
            static_cast<_Implementation *>(this)->setUnknowns(mesh);
        }

    };

    /*
     *
     */

    template<lolita::core::element::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto...>
    struct FiniteElementFinal;

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    requires(!lolita::core::element::PointConcept<_element>)
    struct FiniteElementFinal<_element, _domain, _finite_element...> : public FiniteElementGeometry<FiniteElementFinal, _element, _domain, _finite_element...>
    {

    private:

        using _ElementsPtr = std::tuple<std::shared_ptr<lolita::core::element::FiniteElementModule<_element, _domain, _finite_element>>...>;

        using _Elements = std::tuple<lolita::core::element::FiniteElementModule<_element, _domain, _finite_element>...>;

    public:

        std::basic_string<lolita::character>
        hash()
        const
        {
            std::basic_stringstream<lolita::character> hash;
            auto const & nodes = this->template getComponents<_element.dim_ - 1, 0>();
            for (auto const & node : nodes) {
                hash << node->hash();
            }
            return hash.str();
        }

        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementsPtr> const &
        getElement()
        const
        {
            return std::get<_i>(elements_);
        }

        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementsPtr> &
        getElement()
        {
            return std::get<_i>(elements_);
        }

        static
        void
        makeElement(
                lolita::core::mesh::ElementInitializationData<_element> const & initialization_data,
                lolita::core::mesh::MeshData2<_domain, _finite_element...> & mesh_data_arg
        )
        {
            /*
             *
             */
            auto get_element_hash = [&] <lolita::core::element::Element __element> (
                    std::array<lolita::index, __element.num_nodes_> node_tags
            )
            {
                auto element_hash = std::basic_stringstream<lolita::character>();
                if constexpr (node_tags.size() > 0) {
                    std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
                }
                for (int i = 0; i < node_tags.size(); ++i) {
                    element_hash << node_tags[i];
                }
                return element_hash.str();
            };
            /*
             * make element
             */
            auto make_element = [&] <lolita::core::element::Element __element = _element, lolita::index _i = 0, lolita::index _j = 0> (
                    std::shared_ptr<lolita::core::element::FiniteElementFinal<__element, _domain, _finite_element...>> & ptr_element,
                    lolita::core::mesh::ElementInitializationData<__element> const & element_initialization_data,
                    auto & make_element_imp
            )
            mutable
            {
                if constexpr (!lolita::core::element::PointConcept<__element>) {
                    auto const constexpr _is_initialized = _i == 0 && _j == 0;
                    auto const constexpr _component = element::component<__element, _domain, _i, _j>();
                    auto const constexpr _component_coordinates = element::elementPosition<_domain, _component>();
                    auto const constexpr _neighbour_coordinates = element::neighbourPosition<_component, _domain, __element>();
                    auto const constexpr _element_coordinates = element::elementPosition<_domain, __element>();
                    auto const constexpr _node_coordinates = element::componentPosition<__element, _domain, element::pnt_00>();
                    using __Component = lolita::core::element::FiniteElementFinal<_component, _domain, _finite_element...>;
                    using __Self = lolita::core::element::FiniteElementFinal<__element, _domain, _finite_element...>;
                    using __Elements = std::tuple<lolita::core::element::FiniteElementModule<__element, _domain, _finite_element>...>;
                    auto & components = mesh_data_arg.elements_.template getElements<_component_coordinates[0], _component_coordinates[1]>();
                    auto & element_component_array = ptr_element->template getComponents<_i, _j>();
                    for (auto i = 0; i < element_component_array.size(); ++i) {
                        auto component_hash = std::basic_string<lolita::character>();
                        if constexpr(!lolita::core::element::PointConcept<_component>) {
                            auto component_initialization_data = lolita::core::mesh::ElementInitializationData<_component>();
                            for (int j = 0; j < _component.num_nodes_; ++j) {
                                auto const k = __Self::template getComponentNodeConnection<_i, _j>(i, j);
                                component_initialization_data.node_tags_[j] = element_initialization_data.node_tags_[k];
                            }
                            component_initialization_data.tag_ = components.size();
                            component_hash = get_element_hash.template operator ()<_component>(component_initialization_data.node_tags_);
                            if (!components.contains(component_hash)) {
                                auto ptr_component = std::make_shared<__Component>(__Component());
                                make_element_imp.template operator ()<_component, 0, 0>(ptr_component, component_initialization_data, make_element_imp);
                            }
                        }
                        else {
                            component_hash = std::to_string(element_initialization_data.node_tags_[__Self::template getComponentNodeConnection<_i, _j>(i, 0)]);
                        }
                        element_component_array[i] = components[component_hash];
                        components[component_hash]->template getNeighbours<_neighbour_coordinates[0], _neighbour_coordinates[1]>().push_back(ptr_element);
                    }
                    if constexpr (_j < element::numComponents<__element, _domain, _i>() - 1) {
                        make_element_imp.template operator()<__element, _i, _j + 1u>(ptr_element, element_initialization_data, make_element_imp);
                    }
                    else if constexpr (_i < element::numComponents<__element, _domain>() - 1) {
                        make_element_imp.template operator()<__element, _i + 1u, 0u>(ptr_element, element_initialization_data, make_element_imp);
                    }
                    if constexpr (_is_initialized) {
                        auto & elements = mesh_data_arg.elements_.template getElements<_element_coordinates[0], _element_coordinates[1]>();
                        auto const & nodes = ptr_element->template getComponents<_node_coordinates[0], _node_coordinates[1]>();
                        auto domains = std::unordered_set<std::shared_ptr<std::basic_string<lolita::character>>>();
                        for (auto const & domain : nodes[0]->domains_) {
                            auto has_domain = true;
                            for (int j = 1; j < __element.num_nodes_; ++j) {
                                has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                            }
                            if (has_domain) {
                                domains.insert(domain);
                            }
                        }
                        ptr_element->tag_ = element_initialization_data.tag_;
                        ptr_element->domains_.assign(domains.begin(), domains.end());
                        /*
                         *
                         */
                        auto make_elements = [&] <lolita::index _k = 0> (
                                auto & make_elements_imp
                        )
                        {
                            using __Element = std::tuple_element_t<_k, __Elements>;
                            ptr_element->template getElement<_k>() = std::make_shared<__Element>(__Element());
                            if constexpr (_k < sizeof...(_finite_element) - 1) {
                                make_elements_imp.template operator()<_k + 1u>(make_elements_imp);
                            }
                        };
                        /*
                         *
                         */
                        make_elements(make_elements);
                        auto element_hash = get_element_hash.template operator ()<__element>(element_initialization_data.node_tags_);
                        elements.insert({element_hash, ptr_element});
                    }
                }
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FiniteElementFinal>(FiniteElementFinal());
            make_element(ptr_element, initialization_data, make_element);
        }

        void
        init(
                lolita::core::mesh::MeshData2<_domain, _finite_element...> & mesh_data_arg
        )
        {
            /*
             *
             */
            auto initialize_element = [&] <lolita::index _k = 0u> (
                    auto & initialize_element_imp
            )
            mutable
            {
                /*
                 *
                 */
                auto initialize_components = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_components_imp
                )
                mutable
                {
                    for (int i = 0; i < lolita::core::element::numComponents<_element, _domain, __i, __j>(); ++i) {
                        auto & rhs = this->template getComponents<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getComponents<__i, __j>()[i];
                        lhs = rhs;
                    }
                    if constexpr (__j < lolita::core::element::numComponents<_element, _domain, __i>() - 1) {
                        initialize_components_imp.template operator()<__i, __j + 1u>(initialize_components_imp);
                    }
                    else if constexpr (__i < lolita::core::element::numComponents<_element, _domain>() - 1) {
                        initialize_components_imp.template operator()<__i + 1u, 0u>(initialize_components_imp);
                    }
                };
                /*
                 *
                 */
                auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_neighbours_imp
                )
                mutable
                {
                    for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                        auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (__j < lolita::core::element::numNeighbours<_element, _domain, __i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core::element::numNeighbours<_element, _domain>() - 1) {
                        initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                    }
                };
                /*
                 *
                 */
                initialize_components(initialize_components);
                initialize_neighbours(initialize_neighbours);
                this->template getElement<_k>()->tag_ = this->tag_;
                this->template getElement<_k>()->domains_ = this->domains_;
                this->template getElement<_k>()->initialize(mesh_data_arg);
                if constexpr (_k < sizeof...(_finite_element) - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        _ElementsPtr elements_;

    };

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    requires(lolita::core::element::PointConcept<_element>)
    struct FiniteElementFinal<_element, _domain, _finite_element...> : public FiniteElementGeometry<FiniteElementFinal, _element, _domain, _finite_element...>
    {

    private:

        using _ElementsPtr = std::tuple<std::shared_ptr<FiniteElementModule<_element, _domain, _finite_element>>...>;

        using _Elements = std::tuple<FiniteElementModule<_element, _domain, _finite_element>...>;

    public:

        std::basic_string<lolita::character>
        hash()
        const
        {
            return std::to_string(this->tag_);
        }

        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementsPtr> const &
        getElement()
        const
        {
            return std::get<_i>(elements_);
        }

        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementsPtr> &
        getElement()
        {
            return std::get<_i>(elements_);
        }

        static
        void
        makeElement(
                lolita::core::mesh::ElementInitializationData<_element> const & initialization_data,
                lolita::core::mesh::MeshData2<_domain, _finite_element...> & mesh_data_arg
        )
        {
            /*
             *
             */
            auto make_element = [&] (
                    std::shared_ptr<lolita::core::element::FiniteElementFinal<_element, _domain, _finite_element...>> & ptr_element
            )
            mutable
            {
                auto const constexpr _element_coordinates = element::elementPosition<_domain, _element>();
                ptr_element->tag_ = initialization_data.tag_;
                ptr_element->domains_ = initialization_data.domains_;
                ptr_element->components_ = std::make_shared<lolita::geometry::Point>(initialization_data.coordinates_);
                /*
                 *
                 */
                auto make_elements = [&] <lolita::index _k = 0> (
                        auto & make_elements_imp
                )
                {
                    ptr_element->template getElement<_k>() = std::make_shared<std::tuple_element_t<_k, _Elements>>(std::tuple_element_t<_k, _Elements>());
                    if constexpr (_k < sizeof...(_finite_element) - 1) {
                        make_elements_imp.template operator()<_k + 1u>(make_elements_imp);
                    }
                };
                /*
                 *
                 */
                make_elements(make_elements);
                auto elem_hash = std::to_string(initialization_data.tag_);
                mesh_data_arg.elements_.template getElements<_element_coordinates[0], _element_coordinates[1]>().insert({elem_hash, ptr_element});
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FiniteElementFinal>(FiniteElementFinal());
            make_element(ptr_element);
        }

        void
        init(
                lolita::core::mesh::MeshData2<_domain, _finite_element...> & mesh_data_arg
        )
        {
            /*
             *
             */
            auto initialize_element = [&] <lolita::index _k = 0u> (
                    auto & initialize_element_imp
            )
            constexpr mutable
            {
                /*
                 *
                 */
                auto initialize_coordinates = [&] ()
                {
                    this->template getElement<_k>()->components_ = this->components_;
                };
                /*
                 *
                 */
                auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_neighbours_imp
                )
                mutable
                {
                    for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                        auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (__j < lolita::core::element::numNeighbours<_element, _domain, __i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core::element::numNeighbours<_element, _domain>() - 1) {
                        initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                    }
                };
                /*
                 *
                 */
                initialize_coordinates();
                initialize_neighbours(initialize_neighbours);
                this->template getElement<_k>()->tag_ = this->tag_;
                this->template getElement<_k>()->domains_ = this->domains_;
                this->template getElement<_k>()->initialize(mesh_data_arg);
                if constexpr (_k < sizeof...(_finite_element) - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        _ElementsPtr elements_;

    };

}

#endif //LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
