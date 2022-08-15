//
// Created by dsiedel on 01/04/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_QUADRATURE_GAUSS_SEGMENT_HXX
#define FETA_FETA_CORE_ELEMENT_QUADRATURE_GAUSS_SEGMENT_HXX

#include "new/feta_core_element_quadrature_description.hxx"
#include "new/feta_core_element_shape_description.hxx"
#include "new/feta_core_element_element_description.hxx"
#include "new/_feta_array.hxx"

namespace feta::core::element
{

    namespace detail
    {

        template<ShapeDescription, Quadrature, Indx>
        struct ShapeQuadraturePolicy;

    }

    template<ShapeDescription S, QuadratureDescription Q>
    using ShapeQuadrature = detail::ShapeQuadraturePolicy<S, Q.quadrature, Q.order>;

    namespace detail
    {

        /**
         * @brief
         */
        template<>
        struct ShapeQuadraturePolicy<shape::seg, Quadrature::Gauss, 1> {

            /**
             * @brief
             */
            auto const static constexpr dim_quadrature = 1;

            /**
             * @brief
             */
            auto const static constexpr reference_points = Array<Real, dim_quadrature, shape::seg.ord_shape>{
                    +0.0000000000000000
            };

            /**
             * @brief
             */
            auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                    +2.0000000000000000
            };

        };

        /**
         * @brief
         */
        template<>
        struct ShapeQuadraturePolicy<shape::seg, Quadrature::Gauss, 2> {

            /**
             * @brief
             */
            auto const static constexpr dim_quadrature = 2;

            /**
             * @brief
             */
            auto const static constexpr reference_points = Array<Real, dim_quadrature, shape::seg.ord_shape>{
                    -0.5773502691896257,
                    +0.5773502691896257
            };

            /**
             * @brief
             */
            auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                    +1.0000000000000000,
                    +1.0000000000000000
            };

        };

        /**
         * @brief
         */
        template<>
        struct ShapeQuadraturePolicy<shape::seg, Quadrature::Gauss, 3> {

            /**
             * @brief
             */
            auto const static constexpr dim_quadrature = 3;

            /**
             * @brief
             */
            auto const static constexpr reference_points = Array<Real, dim_quadrature, shape::seg.ord_shape>{
                    -0.7745966692414834,
                    +0.0000000000000000,
                    +0.7745966692414834
            };

            /**
             * @brief
             */
            auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                    +0.5555555555555557,
                    +0.8888888888888888,
                    +0.5555555555555557
            };

        };

        /**
         * @brief
         */
        template<>
        struct ShapeQuadraturePolicy<shape::seg, Quadrature::Gauss, 4> {

            /**
             * @brief
             */
            auto const static constexpr dim_quadrature = 3;

            /**
             * @brief
             */
            auto const static constexpr reference_points = Array<Real, dim_quadrature, shape::seg.ord_shape>{
                    -0.7745966692414834,
                    +0.0000000000000000,
                    +0.7745966692414834
            };

            /**
             * @brief
             */
            auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                    +0.5555555555555557,
                    +0.8888888888888888,
                    +0.5555555555555557
            };

        };

    }

}

#endif //FETA_FETA_CORE_ELEMENT_QUADRATURE_GAUSS_SEGMENT_HXX
