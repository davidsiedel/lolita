//
// Created by dsiedel on 05/05/22.
//

#ifndef LOLITA_LOLITA_ELEMENT_HXX
#define LOLITA_LOLITA_ELEMENT_HXX

#include "lolita/lolita.hxx"

namespace lolita::core::element
{

    struct Element
    {

        constexpr
        lolita::boolean
        operator==(
                Element const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                Element const &
                other
        )
        const = default;

        constexpr
        lolita::index 
        dimPoint()
        const
        {
            return dim == 0 ? 1 : dim;
        }

        lolita::index tag;

        lolita::index dim;

        lolita::index ord;

        lolita::index num_nodes;

    };

    template<Element element, lolita::geometry::Domain domain, auto C>
    struct FiniteElement;

    template<Element element, lolita::geometry::Domain domain, auto C>
    struct FiniteElement{};

    Element const static constexpr pnt_00 = Element{0, 0, 0, 1};
    //
    Element const static constexpr seg_02 = Element{1, 1, 1, 2};
    Element const static constexpr seg_03 = Element{5, 1, 2, 3};
    //
    Element const static constexpr tri_03 = Element{2, 2, 1, 3};
    Element const static constexpr tri_06 = Element{6, 2, 2, 6};
    Element const static constexpr qua_04 = Element{3, 2, 1, 4};
    Element const static constexpr qua_08 = Element{7, 2, 2, 8};
    //
    Element const static constexpr tet_04 = Element{4, 3, 1, 4};
    Element const static constexpr tet_12 = Element{8, 3, 2, 12};
    //
    // auto const static constexpr pg5_5 = Element{5, 2, 1, 5};

    /*
     *
     */

    template<Element element, lolita::geometry::Domain domain>
    concept Cell = domain.dim_ - element.dim == 0;

    template<Element element, lolita::geometry::Domain domain>
    concept Face = domain.dim_ - element.dim == 1;

    template<Element element, lolita::geometry::Domain domain>
    concept Edge = domain.dim_ - element.dim == 2;

    template<Element element, lolita::geometry::Domain domain>
    concept Node = domain.dim_ - element.dim == 3;

    /*
     *
     */

    template<Element element>
    concept Point = element == pnt_00;

    template<Element element>
    concept Segment = element == seg_02 || element == seg_03;

    template<Element element>
    concept Triangle = element == tri_03 || element == tri_06;

    template<Element element>
    concept Quadrangle = element == qua_04 || element == qua_08;

    template<Element element>
    concept Tetrahedron = element == tet_04 || element == tet_12;

    /*
     *
     */

    template<Element element>
    concept Curve = Segment<element>;

    template<Element element>
    concept Surface = Triangle<element> || Quadrangle<element>;

    template<Element element>
    concept Volume = Tetrahedron<element>;

    /*
     *
     */

    template<template<Element> typename T>
    using Points = std::tuple<T<pnt_00>>;

    template<template<Element> typename T>
    using Curves = std::tuple<T<seg_02>>;

    template<template<Element> typename T>
    using Facets = std::tuple<T<tri_03>, T<qua_04>>;

    template<template<Element> typename T>
    using Solids = std::tuple<T<tet_04>>;

    template<typename T>
    struct test_struct : public std::false_type {};

    template<template<Element, lolita::geometry::Domain, auto...> typename T, Element E, lolita::geometry::Domain D, auto... A>
    struct test_struct<T<E, D, A...>> : public std::true_type {};

    template<typename T>
    concept FiniteElementType = test_struct<T>::value;

    /*
     *
     */

    template<lolita::geometry::Domain domain, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
    using Points2 = std::tuple<T<pnt_00, domain, A...>>;

    template<lolita::geometry::Domain domain, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
    using Curves2 = std::tuple<T<seg_02, domain, A...>>;

    template<lolita::geometry::Domain domain, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
    using Facets2 = std::tuple<T<tri_03, domain, A...>, T<qua_04, domain, A...>>;

    template<lolita::geometry::Domain domain, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
    using Solids2 = std::tuple<T<tet_04, domain, A...>>;

    namespace detail
    {

        template<Element E>
        using ElementNodeConnectivity = std::array<lolita::index, E.num_nodes>;

//        template<lolita::index D, template<Element> typename T>
        template<lolita::index D, template<Element, auto...> typename T>
        struct ElementsPolicy;

//        template<template<Element> typename T>
        template<template<Element, auto...> typename T>
        struct ElementsPolicy<1, T>
        {

            using Type = std::tuple<Points<T>, Curves<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementsPolicy<2, T>
        {

            using Type = std::tuple<Points<T>, Curves<T>, Facets<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementsPolicy<3, T>
        {

            using Type = std::tuple<Points<T>, Curves<T>, Facets<T>, Solids<T>>;

        };

        template<lolita::index I, lolita::index D, template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy;

        /*
         * OD
         */

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<0, 3, T>
        {

            using Type = std::tuple<Curves<T>, Facets<T>, Solids<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<0, 2, T>
        {

            using Type = std::tuple<Curves<T>, Facets<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<0, 1, T>
        {

            using Type = std::tuple<Curves<T>>;

        };

        /*
         * 1D
         */

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<1, 3, T>
        {

            using Type = std::tuple<Curves<T>, Facets<T>, Solids<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<1, 2, T>
        {

            using Type = std::tuple<Curves<T>, Facets<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<1, 1, T>
        {

            using Type = std::tuple<Curves<T>>;

        };

        /*
         * 2D
         */

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<2, 3, T>
        {

            using Type = std::tuple<Facets<T>, Solids<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<2, 2, T>
        {

            using Type = std::tuple<Facets<T>>;

        };

        /*
         * 3D
         */

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<3, 3, T>
        {

            using Type = std::tuple<Solids<T>>;

        };

        template<Element E, lolita::index D, template<Element, auto...> typename T, auto... A>
        struct ElementNeighbourArrayPolicy
        {

        private:

            lolita::index const static constexpr I = E.dim;

            template<Element EE>
            using NeighbourArray = std::vector<T<EE, A...>>;

        public:

            using Type = typename ElementNeighbourhoodPolicy<I, D, NeighbourArray>::Type;

        };

        template<Element E, lolita::geometry::Domain D, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
        struct ElementNeighbourArrayPolicy2
        {

        private:

            template<Element EE>
            using NeighbourArray = std::vector<T<EE, D, A...>>;

        public:

            using Type = typename ElementNeighbourhoodPolicy<E.dim, D.dim_, NeighbourArray>::Type;

        };

//        template<lolita::geometry::Domain Dmn, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
//        struct ElementsPolicy2;
//
//        template<lolita::geometry::Domain Dmn, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
//        requires(Dmn.dim_ == 1)
//        struct ElementsPolicy2<Dmn, T, A...>
//        {
//
//            using Type = std::tuple<Points2<Dmn, T, A...>, Curves2<Dmn, T, A...>>;
//
//        };
//
//        template<lolita::geometry::Domain Dmn, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
//        requires(Dmn.dim_ == 2)
//        struct ElementsPolicy2<Dmn, T, A...>
//        {
//
//            using Type = std::tuple<Points2<Dmn, T, A...>, Curves2<Dmn, T, A...>, Facets2<Dmn, T, A...>>;
//
//        };
//
//        template<lolita::geometry::Domain Dmn, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
//        requires(Dmn.dim_ == 3)
//        struct ElementsPolicy2<Dmn, T, A...>
//        {
//
//            using Type = std::tuple<Points2<Dmn, T, A...>, Curves2<Dmn, T, A...>, Facets2<Dmn, T, A...>, Solids2<Dmn, T, A...>>;
//
//        };
//
////        template<lolita::index I, lolita::index D, template<Element, auto...> typename T>
//        template<Element E, lolita::geometry::Domain Dmn, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
//        struct ElementNeighbourhoodPolicy2;
//
//        /*
//         * OD
//         */
//
//        template<Element E, lolita::geometry::Domain Dmn, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
//        requires(Point<E>)
//        struct ElementNeighbourhoodPolicy2<E, Dmn, T, A...>
//        {
//
//            using Type = std::tuple<Curves2<Dmn, T, A...>, Facets2<Dmn, T, A...>, Solids2<Dmn, T, A...>>;
//
//        };
//
//        template<template<Element, auto...> typename T>
//        struct ElementNeighbourhoodPolicy2<0, 2, T>
//        {
//
//            using Type = std::tuple<Curves<T>, Facets<T>>;
//
//        };
//
//        template<template<Element, auto...> typename T>
//        struct ElementNeighbourhoodPolicy2<0, 1, T>
//        {
//
//            using Type = std::tuple<Curves<T>>;
//
//        };
//
//        /*
//         * 1D
//         */
//
//        template<template<Element, auto...> typename T>
//        struct ElementNeighbourhoodPolicy2<1, 3, T>
//        {
//
//            using Type = std::tuple<Curves<T>, Facets<T>, Solids<T>>;
//
//        };
//
//        template<template<Element, auto...> typename T>
//        struct ElementNeighbourhoodPolicy2<1, 2, T>
//        {
//
//            using Type = std::tuple<Curves<T>, Facets<T>>;
//
//        };
//
//        template<template<Element, auto...> typename T>
//        struct ElementNeighbourhoodPolicy2<1, 1, T>
//        {
//
//            using Type = std::tuple<Curves<T>>;
//
//        };
//
//        /*
//         * 2D
//         */
//
//        template<template<Element, auto...> typename T>
//        struct ElementNeighbourhoodPolicy2<2, 3, T>
//        {
//
//            using Type = std::tuple<Facets<T>, Solids<T>>;
//
//        };
//
//        template<template<Element, auto...> typename T>
//        struct ElementNeighbourhoodPolicy2<2, 2, T>
//        {
//
//            using Type = std::tuple<Facets<T>>;
//
//        };
//
//        /*
//         * 3D
//         */
//
//        template<template<Element, auto...> typename T>
//        struct ElementNeighbourhoodPolicy2<3, 3, T>
//        {
//
//            using Type = std::tuple<Solids<T>>;
//
//        };

    }

    template<Element E, lolita::index D, template<Element, auto...> typename T, auto... A>
    using ElementNeighbourArray = typename detail::ElementNeighbourArrayPolicy<E, D, T, A...>::Type;

    template<Element E, lolita::geometry::Domain Dmn, template<Element, lolita::geometry::Domain, auto...> typename T, auto... A>
    using ElementNeighbourArray2 = typename detail::ElementNeighbourArrayPolicy2<E, Dmn, T, A...>::Type;

    template<Element E>
    struct ElementGeometry;

    template<Element E>
    requires(E == pnt_00)
    struct ElementGeometry<E>
    {

        Element const static constexpr element = E;

        constexpr
        lolita::boolean
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        /**
         * @brief
         */
        std::array<std::array<lolita::real, 3>, 1> const static constexpr reference_nodes = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
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
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point,
                lolita::index
                derivative_direction
        )
        {
            return lolita::real(0);
        }

    };

    template<Element E>
    requires(E == seg_02)
    struct ElementGeometry<E>
    {

        Element const static constexpr element = E;

        constexpr
        lolita::boolean
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        auto const static constexpr lowers = std::array<Element, 0>{};
        auto const static constexpr higher = std::array<Element, 1>{seg_03};

        std::array<std::array<lolita::real, 3>, 2> const static constexpr reference_nodes = {
                -1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        template<template<Element, auto...> typename T>
        using Components = std::tuple<
                std::tuple<
                        std::array<T<pnt_00>, 2>
                >
        >;

        template<template<Element, auto...> typename T, auto... A>
        requires(std::same_as<std::tuple_element_t<0, std::tuple<std::remove_cvref_t<decltype(A)>...>>, lolita::geometry::Domain>)
        using Components2 = std::tuple<
                std::tuple<
                        std::array<T<pnt_00, A...>, 2>
                >
        >;

        template<lolita::geometry::Domain Dmn, template<Element, auto...> typename T, auto... A>
        using Neighbours2 = std::tuple<Components2<T, Dmn, A...>, ElementNeighbourArray2<E, Dmn, T, A...>>;

        Components<detail::ElementNodeConnectivity> const static constexpr node_connectivity = {
                {
                        {
                                0,
                                1,
                        }
                }
        };

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point
        )
        {
            auto value = lolita::real(0);
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
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point,
                lolita::index
                derivative_direction
        )
        {
            assert(derivative_direction == 0);
            auto value = lolita::real(0);
            value += - nodal_field_values(0) * (1.0 / 2.0);
            value += + nodal_field_values(1) * (1.0 / 2.0);
            return value;
        }

    };

    template<Element E>
    requires(E == tri_03)
    struct ElementGeometry<E>
    {

        Element const static constexpr element = E;

        constexpr
        lolita::boolean
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        std::array<std::array<lolita::real, 3>, 3> const static constexpr reference_nodes = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
        };

        template<template<Element, auto...> typename T>
        using Components = std::tuple<
        std::tuple<
                std::array<T<seg_02>, 3>
                >,
                std::tuple<
                std::array<T<pnt_00>, 3>
                >
        >;

        Components<detail::ElementNodeConnectivity> const static constexpr node_connectivity = {
                {
                        {
                                0, 1,
                                1, 2,
                                2, 0
                        }
                },
                {
                        {
                                0,
                                   1,
                                2,
                        }
                }
        };

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point
        )
        {
            auto value = lolita::real(0);
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
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point,
                lolita::index
                derivative_direction
        )
        {
            assert(0 <= derivative_direction <= 1);
            auto value = lolita::real(0);
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

    template<Element E>
    requires(E == qua_04)
    struct ElementGeometry<E>
    {

        auto const static constexpr element = E;

        constexpr
        lolita::boolean
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes = {
                -1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                -1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
        };

        template<template<Element, auto...> typename T>
        using Components = std::tuple<
        std::tuple<
                std::array<T<seg_02>, 4>
        >,
        std::tuple<
        std::array<T<pnt_00>, 4>
        >
        >;

        auto const static constexpr node_connectivity = Components<detail::ElementNodeConnectivity>{
                {
                        {
                                0, 1,
                                1, 2,
                                2, 3,
                                3, 0,
                        }
                },
                {
                        {
                                0,
                                   1,
                                2,
                                   3,
                        }
                }
        };

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point
        )
        {
            auto value = lolita::real(0);
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
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point,
                lolita::index
                derivative_direction
        )
        {
            assert(0 <= derivative_direction <= 1);
            auto value = lolita::real(0);
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

    template<Element E>
    requires(E == tet_04)
    struct ElementGeometry<E>
    {

        auto const static constexpr element = E;

        constexpr
        lolita::boolean
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000, +1.0000000000000000,
        };

        template<template<Element, auto...> typename T>
        using Components = std::tuple<
        std::tuple<
                std::array<T<tri_03>, 4>
        >,
        std::tuple<
        std::array<T<seg_02>, 6>
        >,
        std::tuple<
        std::array<T<pnt_00>, 4>
        >
        >;

        Components<detail::ElementNodeConnectivity> const static constexpr node_connectivity = {
                {
                        {
                                0, 1, 3,
                                         0, 3, 2,
                                0, 2, 1,
                                         1, 2, 3,
                        }
                },
                {
                        {
                                0, 1,
                                      1, 2,
                                            2, 0,
                                0, 3,
                                      3, 2,
                                            1, 3,
                        }
                },
                {
                        {
                                0,
                                   1,
                                      2,
                                         3,
                        }
                }
        };

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point
        )
        {
            auto value = lolita::real(0);
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
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, element.num_nodes> const &
                nodal_field_values,
                lolita::geometry::Point const &
                reference_point,
                lolita::index
                derivative_direction
        )
        {
            assert(0 <= derivative_direction <= 2);
            auto value = lolita::real(0);
            if (derivative_direction == 0) {
            } else if (derivative_direction == 1) {
            }
            return value;
        }

    };

    template<lolita::index D, template<Element> typename T>
    using Elements = typename detail::ElementsPolicy<D, T>::Type;

    namespace detail
    {

        template<template<Element, auto...> typename T, Element E, auto... A>
        struct ElementComponentsPolicy
        {

        private:

            template<Element B>
            using ElementType = T<B, A...>;

        public:

            using Type = typename ElementGeometry<E>::template Components<ElementType>;

        };

        template<template<Element, auto...> typename T, lolita::index I, lolita::index J, Element E, auto... A>
        struct ElementComponentPolicy
        {

        private:

            using ElementComponents = typename ElementComponentsPolicy<T, E, A...>::Type;

        public:

            using Type = typename ElementComponents::template Type<I>::template Type<J>::Type;

        };

    }

    template<Element E, template<Element, auto...> typename T, auto... A>
    struct ElemCom
    {

    private:

        template<Element B>
        using ElementType = T<B, A...>;

        using Type = typename ElementGeometry<E>::template Components<ElementType>;

        template<lolita::index... K>
        auto const static constexpr layers = std::array<lolita::index, sizeof...(K)>{static_cast<lolita::index>(K)...};

        template<lolita::index... K>
        struct Pol;

        template<lolita::index... K>
        requires(sizeof...(K) == 0)
        struct Pol<K...>
        {

            using Type = ElemCom::Type;

        };

        template<lolita::index... K>
        requires(sizeof...(K) == 1)
        struct Pol<K...>
        {

            using Type = typename ElemCom::Type::template Type<layers<K...>[0]>;

        };

        template<lolita::index... K>
        requires(sizeof...(K) == 2)
        struct Pol<K...>
        {

            using Type = typename ElemCom::Type::template Type<layers<K...>[0]>::template Type<layers<K...>[1]>::Type;

        };

    public:

        template<lolita::index... K>
        using Lay = typename Pol<K...>::Type;

    };

    template<Element E, template<Element, auto...> typename T, auto... A>
    using ElementComponents = typename detail::ElementComponentsPolicy<T, E, A...>::Type;

    template<Element E, lolita::index I, lolita::index J, template<Element, auto...> typename T, auto... A>
    using ElementComponent = typename detail::ElementComponentPolicy<T, I, J, E, A...>::Type;

    using ElementPosition = std::array<lolita::index, 3>;

    template<Element E, Element B>
    static constexpr
    ElementPosition
    elementIndex()
    {
        auto get_position = [] <lolita::index I = 0, lolita::index J = 0, lolita::index K = 0> (auto & self)
                constexpr
        {
        };
        get_position(get_position);
        return ElementPosition{0, 0, 0};
    }

    template<Element E>
    static constexpr
    auto
    elementIndex()
    {
        return 1;
//        return collection::index<typename Elements<3, ElementGeometry>::template Type<E.dim>, ElementGeometry<E>>();
    }

    template<lolita::index I, lolita::index J>
    static constexpr
    Element
    element()
    {
        return lolita::core::element::seg_02;
//        return Elements<3, ElementGeometry>::template Type<I>::template Type<J>::element;
    }

    template<lolita::index I>
    constexpr inline
    lolita::index
    numElements()
    {
        return 1;
//        return Elements<3, ElementGeometry>::template Type<I>::size();
    }

    template<Element E, lolita::index I, lolita::index J>
    constexpr inline
    lolita::index
    component()
    {
        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
        return L::template Type<I>::template Type<J>::Type::element;
    }

    template<Element E, lolita::index I, lolita::index J>
    constexpr inline
    lolita::index
    numComponents()
    {
        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
        return L::template Type<I>::template Type<J>::size();
    }

    template<Element E, lolita::index I>
    constexpr inline
    lolita::index
    numComponents()
    {
        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
        return L::template Type<I>::size();
    }

    template<Element E>
    constexpr inline
    lolita::index
    numComponents()
    {
        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
        return L::size();
    }

    template<Element E, auto D>
    constexpr inline
    lolita::index
    numNeighbours()
    {
        using L = ElementNeighbourArray<E, D, ElementGeometry>;
        return L::size();
    }

    template<Element E, auto D, lolita::index I>
    constexpr inline
    lolita::index
    numNeighbours()
    {
        using L = ElementNeighbourArray<E, D, ElementGeometry>;
        return L::template Type<I>::size();
    }

//    template<Element E, auto D, auto I, auto J>
//    constexpr inline
//    auto
//    numNeighbours()
//    {
//        using L = ElementNeighbourstd::array<E, D, ElementGeometry>;
//        return L::template Type<I>::template Type<J>::size();
//    }

    template<Element element, lolita::finite_element::Quadrature quadrature, lolita::index ord_quadrature>
    struct ElementQuadrature;

    template<Element element, lolita::finite_element::Quadrature quadrature, lolita::index ord_quadrature>
    requires(Point<element>)
    struct ElementQuadrature<element, quadrature, ord_quadrature>
    {

        lolita::index const static constexpr dim_quadrature = 1;

        std::array<std::array<lolita::real, 3>, dim_quadrature> const static constexpr reference_points = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        std::array<lolita::real, dim_quadrature> const static constexpr reference_weights = {
            +1.0000000000000000
        };

    };

    template<Element element, lolita::finite_element::Quadrature quadrature, lolita::index ord_quadrature>
    constexpr inline
    lolita::index
    dimQuadrature()
    {
        return ElementQuadrature<element, quadrature, ord_quadrature>::dim_quadrature;
    }

}

#endif //LOLITA_LOLITA_ELEMENT_HXX
