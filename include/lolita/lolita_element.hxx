//
// Created by dsiedel on 05/05/22.
//

#ifndef LOLITA_LOLITA_ELEMENT_HXX
#define LOLITA_LOLITA_ELEMENT_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"

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
            return dim_ == 0 ? 1 : dim_;
        }

        lolita::index tag_;

        lolita::index dim_;

        lolita::index ord_;

        lolita::index num_nodes_;

    };

    template<Element _element, lolita::geometry::Domain _domain, auto C>
    struct FiniteElement;

    template<Element _element, lolita::geometry::Domain _domain, auto C>
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

    template<Element _element, lolita::geometry::Domain _domain>
    concept CellConcept = _domain.dim_ - _element.dim_ == 0;

    template<Element _element, lolita::geometry::Domain _domain>
    concept FaceConcept = _domain.dim_ - _element.dim_ == 1;

    template<Element _element, lolita::geometry::Domain _domain>
    concept EdgeConcept = _domain.dim_ - _element.dim_ == 2;

    template<Element _element, lolita::geometry::Domain _domain>
    concept NodeConcept = _domain.dim_ - _element.dim_ == 3;

    /*
     *
     */

    template<Element _element>
    concept PointConcept = _element == pnt_00;

    template<Element _element>
    concept SegmentConcept = _element == seg_02 || _element == seg_03;

    template<Element _element>
    concept TriangleConcept = _element == tri_03 || _element == tri_06;

    template<Element _element>
    concept QuadrangleConcept = _element == qua_04 || _element == qua_08;

    template<Element _element>
    concept TetrahedronConcept = _element == tet_04 || _element == tet_12;

    /*
     *
     */

    template<Element _element>
    concept CurveConcept = SegmentConcept<_element>;

    template<Element _element>
    concept SurfaceConcept = TriangleConcept<_element> || QuadrangleConcept<_element>;

    template<Element _element>
    concept VolumeConcept = TetrahedronConcept<_element>;

    /*
     *
     */

//    template<typename T>
//    struct test_struct : public std::false_type {};
//
//    template<template<Element, lolita::geometry::Domain, auto...> typename T, Element E, lolita::geometry::Domain D, auto... A>
//    struct test_struct<T<E, D, A...>> : public std::true_type {};
//
//    template<typename T>
//    concept FiniteElementType = test_struct<T>::value;

    /*
     *
     */

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    using Points = std::tuple<
            _T<pnt_00, _domain, _arg...>
    >;

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    using Curves = std::tuple<
            _T<seg_02, _domain, _arg...>
    >;

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    using Facets = std::tuple<
            _T<tri_03, _domain, _arg...>,
            _T<qua_04, _domain, _arg...>
    >;

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    using Solids = std::tuple<
            _T<tet_04, _domain, _arg...>
    >;

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    using Elements = std::tuple<
            Points<_T, _domain, _arg...>,
            Curves<_T, _domain, _arg...>,
            Facets<_T, _domain, _arg...>,
            Solids<_T, _domain, _arg...>
    >;

    namespace detail
    {

//        lolita::geometry::Domain const static constexpr domain__ = lolita::geometry::Domain("", 3, lolita::geometry::Frame::Cartesian);

//        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
//        using Elements = std::tuple<
//                Points<_T, _domain, _arg...>,
//                Curves<_T, _domain, _arg...>,
//                Facets<_T, _domain, _arg...>,
//                Solids<_T, _domain, _arg...>
//        >;

    }

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    struct ElementCollection
    {

    private:

        template<template<Element, lolita::geometry::Domain, auto...> typename __T, lolita::geometry::Domain __domain, auto... __arg>
        using _Elements = decltype(lolita::utility::tupleSlice<0, _domain.dim_ + 1>(std::declval<Elements<__T, __domain, __arg...>>()));

    public:

        ElementCollection()
        :
        elements_()
        {}

        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, _Elements<_T, _domain, _arg...>>> const &
        getElements()
        const
        {
            return std::get<_j>(std::get<_i>(elements_));
        }

        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, _Elements<_T, _domain, _arg...>>> &
        getElements()
        {
            return std::get<_j>(std::get<_i>(elements_));
        }

        _Elements<_T, _domain, _arg...> elements_;

    };

    namespace detail
    {

        template<Element _element, auto...>
        using ElementNodeConnectivity = std::array<lolita::index, _element.num_nodes_>;

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        struct ElementNeighboursTraits
        {

        private:

            template<Element __element, lolita::geometry::Domain __domain, auto... __arg>
            using _Neighbours = std::vector<_T<__element, __domain, __arg...>>;

            template<Element _element>
            struct _NeighboursTraits;

            template<Element _element>
            requires(PointConcept<_element>)
            struct _NeighboursTraits<_element>
            {

                using Type = decltype(lolita::utility::tupleSlice<1, _domain.dim_ + 1>(std::declval<Elements<_Neighbours, _domain, _arg...>>()));

            };

            template<Element _element>
            requires(!PointConcept<_element>)
            struct _NeighboursTraits<_element>
            {

                using Type = decltype(lolita::utility::tupleSlice<_element.dim_, _domain.dim_ + 1>(std::declval<Elements<_Neighbours, _domain, _arg...>>()));

            };

        public:

            template<Element _element>
            using Type = typename _NeighboursTraits<_element>::Type;

        };

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, Element _element, lolita::geometry::Domain _domain, auto... _arg>
        using ElementNeighbours = typename detail::ElementNeighboursTraits<_T, _domain, _arg...>::template Type<_element>;

    }

    template<Element, auto...>
    struct ElementGeometry;

    template<Element _element, auto... _arg>
    requires(_element == pnt_00)
    struct ElementGeometry<_element, _arg...>
    {

        Element const static constexpr element_ = _element;

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
        std::array<std::array<lolita::real, 3>, 1> const static constexpr reference_nodes_ = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

//        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
//        using Neighbours = std::tuple<
////                std::tuple<>,
//                        detail::ElementNeighbourArray<_T, element, _domain, __arg...>>;

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

//        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
//        using Neighbours = detail::ElementNeighbourArray<_T, element, _domain, __arg...>;

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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

    template<Element _element, auto... _arg>
    requires(_element == seg_02)
    struct ElementGeometry<_element, _arg...>
    {

        Element const static constexpr element_ = _element;

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

        std::array<std::array<lolita::real, 3>, 2> const static constexpr reference_nodes_ = {
                -1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<pnt_00, __arg...>, 2>
                >
        >;

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

//        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
//        using Neighbours = std::tuple<Components<_T, _domain, __arg...>, detail::ElementNeighbourArray<_T, element, _domain, __arg...>>;

        Components<detail::ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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

    template<Element _element, auto... _arg>
    requires(_element == tri_03)
    struct ElementGeometry<_element, _arg...>
    {

        Element const static constexpr element_ = _element;

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

        std::array<std::array<lolita::real, 3>, 3> const static constexpr reference_nodes_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
        };

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<seg_02, __arg...>, 3>
                >,
                std::tuple<
                        std::array<_T<pnt_00, __arg...>, 3>
                >
        >;

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

//        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
//        using Neighbours = std::tuple<Components<_T, _domain, __arg...>, detail::ElementNeighbourArray<_T, element, _domain, __arg...>>;

        Components<detail::ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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

    template<Element _element, auto... _arg>
    requires(_element == qua_04)
    struct ElementGeometry<_element, _arg...>
    {

        auto const static constexpr element_ = _element;

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

        std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                -1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                -1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
        };

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<seg_02, __arg...>, 4>
                >,
                std::tuple<
                        std::array<_T<pnt_00, __arg...>, 4>
                >
        >;

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

//        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
//        using Neighbours = std::tuple<Components<_T, _domain, __arg...>, detail::ElementNeighbourArray<_T, element, _domain, __arg...>>;

        Components<detail::ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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

    template<Element _element, auto... _arg>
    requires(_element == tet_04)
    struct ElementGeometry<_element, _arg...>
    {

        auto const static constexpr element_ = _element;

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

        std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000, +1.0000000000000000,
        };

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<tri_03, __arg...>, 4>
                >,
                std::tuple<
                        std::array<_T<seg_02, __arg...>, 6>
                >,
                std::tuple<
                        std::array<_T<pnt_00, __arg...>, 4>
                >
        >;

        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

//        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
//        using Neighbours = std::tuple<Components<_T, _domain, __arg...>, detail::ElementNeighbourArray<_T, element, _domain, __arg...>>;

        Components<detail::ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const &
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

    /*
     *
     */

    using ElementPosition = std::array<lolita::integer, 2>;

    struct ElementPosition2
    {

        lolita::integer i;

        lolita::integer j;

    };

    template<lolita::geometry::Domain _domain, Element _element>
    static constexpr
    ElementPosition
    elementPosition()
    {
        using _Elements = Elements<ElementGeometry, _domain>;
        auto position = ElementPosition{_element.dim_, -1};
        auto f0 = [&] <lolita::index _i = 0u> (auto & self) constexpr mutable {
            using _Element = std::tuple_element_t<_i, std::tuple_element_t<_element.dim_, _Elements>>;
            if (_Element::element_ == _element) {
                position[1] = _i;
            }
            if constexpr (_i < std::tuple_size_v<std::tuple_element_t<_element.dim_, _Elements>> - 1) {
                self.template operator()<_i + 1u>(self);
            }
        };
        f0(f0);
        return position;
    }

    template<lolita::geometry::Domain _domain, lolita::index I, lolita::index J>
    static constexpr
    Element
    element()
    {
        return std::tuple_element_t<J, std::tuple_element_t<I, Elements<ElementGeometry, _domain>>>::element_;
    }

    template<lolita::geometry::Domain _domain, lolita::index I>
    static constexpr
    lolita::index
    numElements()
    {
        return std::tuple_size_v<std::tuple_element_t<I, Elements<ElementGeometry, _domain>>>;
    }

    template<lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numElements()
    {
        return std::tuple_size_v<Elements<ElementGeometry, _domain>>;
    }

    /*
     *
     */

    template<Element _element, lolita::geometry::Domain _domain, Element _component>
    static constexpr
    ElementPosition
    componentPosition()
    {
        using _Components = typename ElementGeometry<_element>::template Components<ElementGeometry, _domain>;
        auto position = ElementPosition{-1, -1};
        auto f0 = [&] <lolita::index _I = 0u, lolita::index _J = 0u> (
                auto & self
        )
        constexpr mutable
        {
            using _Component = typename std::tuple_element_t<_J, std::tuple_element_t<_I, _Components>>::value_type;
            if (_Component::element_ == _component) {
                position[0] = _I;
                position[1] = _J;
            }
            if constexpr (_J < std::tuple_size_v<std::tuple_element_t<_I, _Components>> - 1) {
                self.template operator()<_I, _J + 1u>(self);
            }
            else if constexpr (_I < std::tuple_size_v<_Components> - 1) {
                self.template operator()<_I + 1u, 0u>(self);
            }
        };
        f0(f0);
        return position;
    }

    template<Element _element, lolita::geometry::Domain _domain, lolita::index _I, lolita::index _J>
    static constexpr
    lolita::core::element::Element
    component()
    {
        using _Components = typename ElementGeometry<_element>::template Components<ElementGeometry, _domain>;
        return std::tuple_element_t<_J, std::tuple_element_t<_I, _Components>>::value_type::element_;
    }

    template<Element _element, lolita::geometry::Domain _domain, lolita::index _I, lolita::index _J>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename ElementGeometry<_element>::template Components<ElementGeometry, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_J, std::tuple_element_t<_I, _Components>>>;
    }

    template<Element _element, lolita::geometry::Domain _domain, lolita::index _I>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename ElementGeometry<_element>::template Components<ElementGeometry, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_I, _Components>>;
    }

    template<Element _element, lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename ElementGeometry<_element>::template Components<ElementGeometry, _domain>;
        return std::tuple_size_v<_Components>;
    }

    /*
     *
     */

    template<Element _element, lolita::geometry::Domain _domain, Element _neighbour>
    static constexpr
    ElementPosition
    neighbourPosition()
    {
        using _Neighbours = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
        auto position = ElementPosition{-1, -1};
        auto f0 = [&] <lolita::index _I = 0u, lolita::index _J = 0u> (
                auto & self
        )
        constexpr mutable
        {
            using _Neighbour = typename std::tuple_element_t<_J, std::tuple_element_t<_I, _Neighbours>>::value_type;
            if (_Neighbour::element_ == _neighbour) {
                position[0] = _I;
                position[1] = _J;
            }
            if constexpr (_J < std::tuple_size_v<std::tuple_element_t<_I, _Neighbours>> - 1) {
                self.template operator()<_I, _J + 1u>(self);
            }
            else if constexpr (_I < std::tuple_size_v<_Neighbours> - 1) {
                self.template operator()<_I + 1u, 0u>(self);
            }
        };
        f0(f0);
        return position;
    }

    template<Element _element, lolita::geometry::Domain _domain, lolita::index _I, lolita::index _J>
    static constexpr
    lolita::core::element::Element
    neighbour()
    {
        using _Neighbours = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
        return std::tuple_element_t<_J, std::tuple_element_t<_I, _Neighbours>>::value_type::element_;
    }

    template<Element _element, lolita::geometry::Domain _domain, lolita::index _I>
    static constexpr
    lolita::index
    numNeighbours()
    {
        using _Neighbours = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_I, _Neighbours>>;
    }

    template<Element _element, lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numNeighbours()
    {
        using _Neighbours = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
        return std::tuple_size_v<_Neighbours>;
    }

    /*
     *
     */

//    using NeighbourPosition = std::array<lolita::integer , 3>;
//
//    template<Element _element, lolita::geometry::Domain _domain, Element _neighbour>
//    static constexpr
//    NeighbourPosition
//    neighbourPosition()
//    {
//        using NeighbourArray = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
//        auto constexpr layer = _element.dim - _neighbour.dim > 0 || PointConcept<_element> ? 0 : 1;
////        auto constexpr layer = _element.dim - _neighbour.dim > 0 ? 0 : 1;
//        auto position = NeighbourPosition{-1, -1, -1};
////        if constexpr(!PointConcept<_element>) {
//            auto f0 = [&] <lolita::index _I = 0u, lolita::index _J = 0u> (
//                    auto & self
//            )
//            constexpr mutable
//            {
//                using Neighbour = typename std::tuple_element_t<_J, std::tuple_element_t<_I, std::tuple_element_t<layer, NeighbourArray>>>::value_type;
//                if (Neighbour::element == _neighbour) {
//                    position[0] = layer;
//                    position[1] = _I;
//                    position[2] = _J;
//                }
//                if constexpr (_J < std::tuple_size_v<std::tuple_element_t<_I, std::tuple_element_t<layer, NeighbourArray>>> - 1) {
//                    self.template operator()<_I, _J + 1u>(self);
//                }
//                else if constexpr (_I < std::tuple_size_v<std::tuple_element_t<layer, NeighbourArray>> - 1) {
//                    self.template operator()<_I + 1u, 0u>(self);
//                }
//            };
//            f0(f0);
////        }
////        else {
////            auto f0 = [&] <lolita::index _I = 0u> (
////                    auto & self
////            )
////            constexpr mutable
////            {
////                using Neighbour = typename std::tuple_element_t<_I, std::tuple_element_t<layer, NeighbourArray>>::value_type;
////                if (Neighbour::element == _neighbour) {
////                    position[0] = layer;
////                    position[1] = _I;
////                }
////                else if constexpr (_I < std::tuple_size_v<std::tuple_element_t<layer, NeighbourArray>> - 1) {
////                    self.template operator()<_I + 1u>(self);
////                }
////            };
////            f0(f0);
////        }
//        return position;
//    }
//
//    template<Element _element, lolita::geometry::Domain _domain, lolita::index _I, lolita::index _J, lolita::index _K>
//    static constexpr
//    lolita::core::element::Element
//    neighbour()
//    {
//        using NeighbourArray = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
//        return std::tuple_element_t<_K, std::tuple_element_t<_J, std::tuple_element_t<_I, NeighbourArray>>>::value_type::element;
//    }
//
////    template<Element _element, lolita::index _I, lolita::index _J, lolita::index _K>
////    static constexpr
////    lolita::integer
////    numNeighbours()
////    {
////        using NeighbourArray = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, detail::domain__>;
////        return std::tuple_size_v<std::tuple_element_t<_K, std::tuple_element_t<_J, std::tuple_element_t<_I, NeighbourArray>>>>;
////    }
//
//    template<Element _element, lolita::geometry::Domain _domain, lolita::index _I, lolita::index _J>
//    static constexpr
//    lolita::integer
//    numNeighbours()
//    {
//        using NeighbourArray = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
//        return std::tuple_size_v<std::tuple_element_t<_J, std::tuple_element_t<_I, NeighbourArray>>>;
//    }
//
//    template<Element _element, lolita::geometry::Domain _domain, lolita::index _I>
//    static constexpr
//    lolita::integer
//    numNeighbours()
//    {
//        using NeighbourArray = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
//        return std::tuple_size_v<std::tuple_element_t<_I, NeighbourArray>>;
//    }
//
//    template<Element _element, lolita::geometry::Domain _domain>
//    static constexpr
//    lolita::integer
//    numNeighbours()
//    {
//        using NeighbourArray = typename ElementGeometry<_element>::template Neighbours<ElementGeometry, _domain>;
//        return std::tuple_size_v<NeighbourArray>;
//    }

    /*
     *
     */

    template<Element, lolita::finite_element::Quadrature, lolita::index>
    struct ElementQuadrature;

    template<Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord_quadrature>
    requires(PointConcept<_element>)
    struct ElementQuadrature<_element, _quadrature, _ord_quadrature>
    {

        lolita::finite_element::Quadrature const static constexpr quadrature_ = _quadrature;

        lolita::index const static constexpr ord_ = _ord_quadrature;

        lolita::index const static constexpr dim_ = 1;

        std::array<std::array<lolita::real, 3>, dim_> const static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        std::array<lolita::real, dim_> const static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };

    template<Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord_quadrature>
    requires(!PointConcept<_element>)
    struct ElementQuadrature<_element, _quadrature, _ord_quadrature>
    {

        lolita::finite_element::Quadrature const static constexpr quadrature_ = _quadrature;

        lolita::index const static constexpr ord_ = _ord_quadrature;

        lolita::index const static constexpr dim_ = 1;

        std::array<std::array<lolita::real, 3>, dim_> const static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        std::array<lolita::real, dim_> const static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };

//    template<Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord_quadrature>
//    constexpr inline
//    lolita::index
//    dimQuadrature()
//    {
//        return ElementQuadrature<_element, _quadrature, _ord_quadrature>::dim_quadrature;
//    }

//    namespace detail
//    {
//
//        template<template<Element, auto...> typename T, Element E, auto... A>
//        struct ElementComponentsPolicy
//        {
//
//        private:
//
//            template<Element B>
//            using ElementType = T<B, A...>;
//
//        public:
//
//            using Type = typename ElementGeometry<E>::template Components<ElementType>;
//
//        };
//
//        template<template<Element, auto...> typename T, lolita::index I, lolita::index J, Element E, auto... A>
//        struct ElementComponentPolicy
//        {
//
//        private:
//
//            using ElementComponents = typename ElementComponentsPolicy<T, E, A...>::Type;
//
//        public:
//
//            using Type = typename ElementComponents::template Type<I>::template Type<J>::Type;
//
//        };
//
//    }
//
//    template<Element E, template<Element, auto...> typename T, auto... A>
//    struct ElemCom
//    {
//
//    private:
//
//        template<Element B>
//        using ElementType = T<B, A...>;
//
//        using Type = typename ElementGeometry<E>::template Components<ElementType>;
//
//        template<lolita::index... K>
//        auto const static constexpr layers = std::array<lolita::index, sizeof...(K)>{static_cast<lolita::index>(K)...};
//
//        template<lolita::index... K>
//        struct Pol;
//
//        template<lolita::index... K>
//        requires(sizeof...(K) == 0)
//        struct Pol<K...>
//        {
//
//            using Type = ElemCom::Type;
//
//        };
//
//        template<lolita::index... K>
//        requires(sizeof...(K) == 1)
//        struct Pol<K...>
//        {
//
//            using Type = typename ElemCom::Type::template Type<layers<K...>[0]>;
//
//        };
//
//        template<lolita::index... K>
//        requires(sizeof...(K) == 2)
//        struct Pol<K...>
//        {
//
//            using Type = typename ElemCom::Type::template Type<layers<K...>[0]>::template Type<layers<K...>[1]>::Type;
//
//        };
//
//    public:
//
//        template<lolita::index... K>
//        using Lay = typename Pol<K...>::Type;
//
//    };
//
//    template<Element E, template<Element, auto...> typename T, auto... A>
//    using ElementComponents = typename detail::ElementComponentsPolicy<T, E, A...>::Type;
//
//    template<Element E, lolita::index I, lolita::index J, template<Element, auto...> typename T, auto... A>
//    using ElementComponent = typename detail::ElementComponentPolicy<T, I, J, E, A...>::Type;



//    template<Element E, Element B>
//    static constexpr
//    ElementPosition
//    elementIndex()
//    {
//        auto get_position = [] <lolita::index I = 0, lolita::index J = 0, lolita::index K = 0> (auto & self)
//                constexpr
//        {
//        };
//        get_position(get_position);
//        return ElementPosition{0, 0, 0};
//    }
//
//    template<Element E>
//    static constexpr
//    auto
//    elementIndex()
//    {
//        return 1;
////        return collection::index<typename Elements<3, ElementGeometry>::template Type<E.dim>, ElementGeometry<E>>();
//    }
//
//    template<lolita::index I, lolita::index J>
//    static constexpr
//    Element
//    element()
//    {
//        return lolita::core::element::seg_02;
////        return Elements<3, ElementGeometry>::template Type<I>::template Type<J>::element;
//    }
//
//    template<lolita::index I>
//    constexpr inline
//    lolita::index
//    numElements()
//    {
//        return 1;
////        return Elements<3, ElementGeometry>::template Type<I>::size();
//    }
//
//    template<Element E, lolita::index I, lolita::index J>
//    constexpr inline
//    lolita::index
//    component()
//    {
//        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
//        return L::template Type<I>::template Type<J>::Type::element;
//    }
//
//    template<Element E, lolita::index I, lolita::index J>
//    constexpr inline
//    lolita::index
//    numComponents()
//    {
//        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
//        return L::template Type<I>::template Type<J>::size();
//    }
//
//    template<Element E, lolita::index I>
//    constexpr inline
//    lolita::index
//    numComponents()
//    {
//        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
//        return L::template Type<I>::size();
//    }
//
//    template<Element E>
//    constexpr inline
//    lolita::index
//    numComponents()
//    {
//        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
//        return L::size();
//    }
//
//    template<Element E, auto D>
//    constexpr inline
//    lolita::index
//    numNeighbours()
//    {
//        using L = ElementNeighbourArray<E, D, ElementGeometry>;
//        return L::size();
//    }
//
//    template<Element E, auto D, lolita::index I>
//    constexpr inline
//    lolita::index
//    numNeighbours()
//    {
//        using L = ElementNeighbourArray<E, D, ElementGeometry>;
//        return L::template Type<I>::size();
//    }

}

#endif //LOLITA_LOLITA_ELEMENT_HXX
