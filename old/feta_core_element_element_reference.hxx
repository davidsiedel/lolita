//
// Created by dsiedel on 29/03/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_ELEMENT_REFERENCE_HXX
#define FETA_FETA_CORE_ELEMENT_ELEMENT_REFERENCE_HXX

#include "new/_feta_aggregate.hxx"
#include "new/_feta_collection.hxx"
#include "new/_feta_array.hxx"
//#include "new/_feta_matrix.hxx"
#include "new/feta_core_element_element_description.hxx"
#include "new/feta_core_elements.hxx"

namespace feta::core::element
{

    template<ElementDescription E, auto M>
    struct Element;

    template<ElementDescription E, auto M>
    struct MixedElement;

    namespace detail
    {

        template<Indx I, Indx D, template<ElementDescription> typename T>
        struct ElementNeighbourhood;

        /*
         * OD
         */

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<0, 3, T>
        {

            using Type = Collection<CurveElements<T>, PlaneElements<T>, SolidElements<T>>;

        };

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<0, 2, T>
        {

            using Type = Collection<CurveElements<T>, PlaneElements<T>>;

        };

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<0, 1, T>
        {

            using Type = Collection<CurveElements<T>>;

        };

        /*
         * 1D
         */

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<1, 3, T>
        {

            using Type = Collection<CurveElements<T>, PlaneElements<T>, SolidElements<T>>;

        };

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<1, 2, T>
        {

            using Type = Collection<CurveElements<T>, PlaneElements<T>>;

        };

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<1, 1, T>
        {

            using Type = Collection<CurveElements<T>>;

        };

        /*
         * 2D
         */

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<2, 3, T>
        {

            using Type = Collection<PlaneElements<T>, SolidElements<T>>;

        };

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<2, 2, T>
        {

            using Type = Collection<PlaneElements<T>>;

        };

        /*
         * 3D
         */

        template<template<ElementDescription> typename T>
        struct ElementNeighbourhood<3, 3, T>
        {

            using Type = Collection<SolidElements<T>>;

        };

        template<ElementDescription E, Indx D, template<ElementDescription> typename T>
        struct ElementNeighbourArrayPolicy
        {

        private:

            auto const static constexpr I = E.shape_description.ord_shape;

            template<ElementDescription EE>
            using NeighbourArray = Array<T<EE>>;

        public:

            using Type = typename ElementNeighbourhood<I, D, NeighbourArray>::Type;

        };

    }

    template<ElementDescription E, Indx D, template<ElementDescription> typename T>
    using ElementNeighbourArray = typename detail::ElementNeighbourArrayPolicy<E, D, T>::Type;

    namespace detail
    {

        template<ShapeDescription, Indx, Indx>
        struct ElementReferencePolicy;

        template<>
        struct ElementReferencePolicy<ShapeDescription(Shape::Point), 0, 0>
        {

            Indx const static constexpr reference_tag = 0;

            /**
             * @brief
             */
            Indx const static constexpr num_nodes = 1;

            /**
             * @brief
             */
            auto const static constexpr reference_nodes = Array<Real, 1, 1>{
                    +0.0000000000000000
            };

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            Real
            getShapeMappingEvaluation(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point
            )
            {
                return nodal_field_values(0);
            }

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @param derivative_direction
             * @return
             */
            static
            Real
            getShapeMappingDerivative(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point,
                    auto
                    derivative_direction
            )
            {
                return Real(0);
            }

        };

        template<>
        struct ElementReferencePolicy<ShapeDescription(Shape::Segment), 1, 0>
        {

            Indx const static constexpr reference_tag = 0;

            /**
             * @brief
             */
            Indx const static constexpr num_nodes = 2;

            /**
             * @brief
             */
            auto const static constexpr reference_nodes = Array<Real, 2, 1>{
                    -1.0000000000000000,
                    +1.0000000000000000
            };

            template<template<ElementDescription> typename T>
            using Components = Collection<
                    Collection<
                            Array<T<pnt_0>, 2>
                    >
            >;

        private:

            using Connectivity = Collection<
                    Collection<Array<Indx, 2, 1>>
            >;

        public:

            auto const static constexpr node_connectivity = Connectivity(
                    Collection(
                            Array<Indx, 2, 1>{
                                    0,
                                    1
                            }
                    )
            );

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            Real
            getShapeMappingEvaluation(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point
            )
            {
                auto value = Real(0);
                value += nodal_field_values(0) * (1.0 / 2.0) * (1.0 - reference_point(0));
                value += nodal_field_values(1) * (1.0 / 2.0) * (1.0 + reference_point(0));
                return value;
            }

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @param derivative_direction
             * @return
             */
            static
            Real
            getShapeMappingDerivative(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point,
                    auto
                    derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 0);
                auto value = Real(0);
                value += - nodal_field_values(0) * (1.0 / 2.0);
                value += + nodal_field_values(1) * (1.0 / 2.0);
                return value;
            }

        };

        template<>
        struct ElementReferencePolicy<ShapeDescription(Shape::Triangle), 1, 0>
        {

            Indx const static constexpr reference_tag = 0;

            /**
             * @brief
             */
            Indx const static constexpr num_nodes = 3;

            /**
             * @brief
             */
            auto const static constexpr reference_nodes = Array<Real, 3, 2>{
                +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000
            };

            template<template<ElementDescription> typename T>
            using Components = Collection<
                    Collection<
                            Array<T<seg_2>, 3>
                    >,
                    Collection<
                            Array<T<pnt_0>, 3>
                    >
            >;

        private:

            using Connectivity = Collection<
                    Collection<Array<Indx, 3, 2>>,
                    Collection<Array<Indx, 3, 1>>
            >;

        public:

            auto const static constexpr node_connectivity = Connectivity(
                    Collection(
                            Array<Indx, 3, 2>{
                                    0, 1,
                                    1, 2,
                                    2, 0
                            }
                    ),
                    Collection(
                            Array<Indx, 3, 1>{
                                    0,
                                    1,
                                    2
                            }
                    )
            );

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            Real
            getShapeMappingEvaluation(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point
            )
            {
                auto value = Real(0);
                value += nodal_field_values(0) * (1.0 - reference_point(0) - reference_point(1));
                value += nodal_field_values(1) * reference_point(0);
                value += nodal_field_values(2) * reference_point(1);
                return value;
            }

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @param derivative_direction
             * @return
             */
            static
            Real
            getShapeMappingDerivative(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point,
                    auto
                    derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 1);
                auto value = Real(0);
                if (derivative_direction == 0) {
                    value += - nodal_field_values(0);
                    value += + nodal_field_values(1);
                    value += + 0.0;
                } else {
                    value += - nodal_field_values(0);
                    value += + 0.0;
                    value += + nodal_field_values(2);
                }
                return value;
            }

        };

        template<>
        struct ElementReferencePolicy<ShapeDescription(Shape::Quadrangle), 1, 0>
        {

            Indx const static constexpr reference_tag = 1;

            /**
             * @brief
             */
            Indx const static constexpr num_nodes = 4;

            /**
             * @brief
             */
            auto const static constexpr reference_nodes = Array<Real, 4, 2>{
                -1.0000000000000000, -1.0000000000000000,
                +1.0000000000000000, -1.0000000000000000,
                +1.0000000000000000, +1.0000000000000000,
                -1.0000000000000000, +1.0000000000000000
            };

            template<template<ElementDescription> typename T>
            using Components = Collection<
                    Collection<
                            Array<T<seg_2>, 4>
                    >,
                    Collection<
                            Array<T<pnt_0>, 4>
                    >
            >;

        private:

            using Connectivity = Collection<
                    Collection<Array<Indx, 4, 2>>,
                    Collection<Array<Indx, 4, 1>>
            >;

        public:

            auto const static constexpr node_connectivity = Connectivity(
                    Collection(
                            Array<Indx, 4, 2>{
                                    0, 1,
                                    1, 2,
                                    2, 3,
                                    3, 0
                            }
                    ),
                    Collection(
                            Array<Indx, 4, 1>{
                                    0,
                                    1,
                                    2,
                                    3
                            }
                    )
            );

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            Real
            getShapeMappingEvaluation(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point
            )
            {
                auto value = Real(0);
                value += nodal_field_values(0) * (1.0 / 4.0) * (1.0 - reference_point(0)) * (1.0 - reference_point(1));
                value += nodal_field_values(1) * (1.0 / 4.0) * (1.0 + reference_point(0)) * (1.0 - reference_point(1));
                value += nodal_field_values(2) * (1.0 / 4.0) * (1.0 + reference_point(0)) * (1.0 + reference_point(1));
                value += nodal_field_values(3) * (1.0 / 4.0) * (1.0 - reference_point(0)) * (1.0 + reference_point(1));
                return value;
            }

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @param derivative_direction
             * @return
             */
            static
            Real
            getShapeMappingDerivative(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point,
                    auto
                    derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 1);
                auto value = Real(0);
                if (derivative_direction == 0) {
                    value += - nodal_field_values(0) * (1.0 / 4.0) * (1.0 - reference_point(1));
                    value += + nodal_field_values(1) * (1.0 / 4.0) * (1.0 - reference_point(1));
                    value += + nodal_field_values(2) * (1.0 / 4.0) * (1.0 + reference_point(1));
                    value += - nodal_field_values(3) * (1.0 / 4.0) * (1.0 + reference_point(1));
                } else {
                    value += - nodal_field_values(0) * (1.0 / 4.0) * (1.0 - reference_point(0));
                    value += - nodal_field_values(1) * (1.0 / 4.0) * (1.0 + reference_point(0));
                    value += + nodal_field_values(2) * (1.0 / 4.0) * (1.0 + reference_point(0));
                    value += + nodal_field_values(3) * (1.0 / 4.0) * (1.0 - reference_point(0));
                }
                return value;
            }

        };

        template<>
        struct ElementReferencePolicy<ShapeDescription(Shape::Tetrahedron), 1, 0>
        {

            Indx const static constexpr reference_tag = 0;

            /**
             * @brief
             */
            Indx const static constexpr num_nodes = 4;

            /**
             * @brief
             */
            auto const static constexpr reference_nodes = Array<Real, 4, 3>{
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000, +1.0000000000000000
            };

            template<template<ElementDescription> typename T>
            using Components = Collection<
                    Collection<
                            Array<T<tri_3>, 4>
                    >,
                    Collection<
                            Array<T<seg_2>, 6>
                    >,
                    Collection<
                            Array<T<pnt_0>, 4>
                    >
            >;

        private:

            using Connectivity = Collection<
                    Collection<Array<Indx, 4, 3>>,
                    Collection<Array<Indx, 6, 2>>,
                    Collection<Array<Indx, 4, 1>>
            >;

        public:

            auto const static constexpr node_connectivity = Connectivity(
                    Collection(
                            Array<Indx, 4, 3>{
                                    0, 1, 3,
                                    0, 3, 2,
                                    0, 2, 1,
                                    1, 2, 3
                            }
                    ),
                    Collection(
                            Array<Indx, 6, 2>{
                                    0, 1,
                                    1, 2,
                                    2, 0,
                                    0, 3,
                                    3, 2,
                                    1, 3
                            }
                    ),
                    Collection(
                            Array<Indx, 4, 1>{
                                    0,
                                    1,
                                    2,
                                    3
                            }
                    )
            );

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            Real
            getShapeMappingEvaluation(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point
            )
            {
                auto value = Real(0);
                return value;
            }

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @param derivative_direction
             * @return
             */
            static
            Real
            getShapeMappingDerivative(
                    auto const &
                    nodal_field_values,
                    auto const &
                    reference_point,
                    auto
                    derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 2);
                auto value = Real(0);
                if (derivative_direction == 0) {
                } else if (derivative_direction == 1) {
                }
                return value;
            }

        };

    }

    template<ElementDescription E>
    using ElementReference = detail::ElementReferencePolicy<E.shape_description, E.ord_element, E.tag>;

    template<ElementDescription E>
    constexpr inline
    auto
    tag()
    {
        return ElementReference<E>::reference_tag;
    }

    template<ElementDescription E>
    constexpr inline
    auto
    numNodes()
    {
        return ElementReference<E>::num_nodes;
    }

    template<ElementDescription E>
    constexpr inline
    auto
    ordShape()
    {
        return E.shape_description.ord_shape;
    }

    template<ElementDescription E>
    constexpr inline
    auto
    dimShape()
    {
        return E.shape_description.dim_shape;
    }

    template<ElementDescription E, Indx I, Indx J>
    constexpr inline
    auto
    numComponents()
    {
        using L = typename ElementReference<E>::template Components<ElementReference>;
        return L::template Type<I>::template Type<J>::size();
    }

    template<ElementDescription E, Indx I>
    constexpr inline
    auto
    numComponents()
    {
        using L = typename ElementReference<E>::template Components<ElementReference>;
        return L::template Type<I>::size();
    }

    template<ElementDescription E>
    constexpr inline
    auto
    numComponents()
    {
        using L = typename ElementReference<E>::template Components<ElementReference>;
        return L::size();
    }

}

#endif //FETA_FETA_CORE_ELEMENT_ELEMENT_REFERENCE_HXX
