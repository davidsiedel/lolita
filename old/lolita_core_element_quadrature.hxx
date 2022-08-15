//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_CORE_ELEMENT_QUADRATURE_HXX
#define LOLITA_LOLITA_CORE_ELEMENT_QUADRATURE_HXX

#include "lolita/lolita_core_element_geometry.hxx"
#include "lolita/lolita_core.hxx"

namespace lolita::core::element::quadrature
{

    template<Element E, Quadrature Q>
    struct ShapeQuadrature;

    template<Element E, Quadrature Q>
    requires(E.dim == 1 && Q.quadrature == QuadratureRule::Gauss && Q.ord == 1)
    struct ShapeQuadrature<E, Q>
    {

        auto const static constexpr dim_quadrature = 1;

        auto const static constexpr reference_points = Array<Real, dim_quadrature, 1>{
                +0.0000000000000000
        };

        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +2.0000000000000000
        };

    };

    template<Element E, Quadrature Q>
    requires(E.dim == 1 && Q.quadrature == QuadratureRule::Gauss && Q.ord == 2)
    struct ShapeQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 2;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, 1>{
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

    template<Element E, Quadrature Q>
    requires(E.dim == 1 && Q.quadrature == QuadratureRule::Gauss && Q.ord == 3)
    struct ShapeQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 3;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, 1>{
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

    template<Element E, Quadrature Q>
    requires(E.dim == 1 && Q.quadrature == QuadratureRule::Gauss && Q.ord == 4)
    struct ShapeQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 3;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, 1>{
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

    namespace detail
    {

        template<Indx, Indx, QuadratureRule, Indx>
        struct ShapeQuadraturePolicy;

        /**
         * @brief
         */
        template<>
        struct ShapeQuadraturePolicy<1, 1, QuadratureRule::Gauss, 1>
        {

            /**
             * @brief
             */
            auto const static constexpr dim_quadrature = 1;

            /**
             * @brief
             */
            auto const static constexpr reference_points = Array<Real, dim_quadrature, 1>{
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
        struct ShapeQuadraturePolicy<1, 1, QuadratureRule::Gauss, 2>
        {

            /**
             * @brief
             */
            auto const static constexpr dim_quadrature = 2;

            /**
             * @brief
             */
            auto const static constexpr reference_points = Array<Real, dim_quadrature, 1>{
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
        struct ShapeQuadraturePolicy<1, 1, QuadratureRule::Gauss, 3>
        {

            /**
             * @brief
             */
            auto const static constexpr dim_quadrature = 3;

            /**
             * @brief
             */
            auto const static constexpr reference_points = Array<Real, dim_quadrature, 1>{
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
        struct ShapeQuadraturePolicy<1, 1, QuadratureRule::Gauss, 4>
        {

            /**
             * @brief
             */
            auto const static constexpr dim_quadrature = 3;

            /**
             * @brief
             */
            auto const static constexpr reference_points = Array<Real, dim_quadrature, 1>{
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

#endif //LOLITA_LOLITA_CORE_ELEMENT_QUADRATURE_HXX
