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
                Element const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                Element const & other
        )
        const = default;

        lolita::index tag_;

        lolita::index dim_;

        lolita::index ord_;

        lolita::index num_nodes_;

    };

    template<lolita::core::element::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto...>
    struct FiniteElementFinal;

    template<lolita::finite_element::FiniteElementConcept auto... _finite_element>
    struct MixedElementPolicy
    {

        template<
                template<lolita::core::element::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto> typename _T,
                lolita::core::element::Element _element,
                lolita::geometry::Domain _domain
        >
        using MixedElement = std::tuple<_T<_element, _domain, _finite_element>...>;

        template<lolita::finite_element::FiniteElementConcept auto __finite_element>
        static constexpr
        lolita::index
        index()
        {
            auto index = sizeof...(_finite_element);
            using _Elements = std::tuple<std::remove_cvref_t<decltype(_finite_element)>...>;
            auto const constexpr elements = _Elements{_finite_element...};
            auto set_index = [&] <lolita::index _i = 0u> (auto & self)
            constexpr mutable
            {
                if constexpr (std::is_same_v<std::tuple_element_t<_i, _Elements>, std::remove_cvref_t<decltype(__finite_element)>>) {
                    if (__finite_element == std::get<_i>(elements)) {
                        index = _i;
                    }
                }
                if constexpr(_i < sizeof...(_finite_element) - 1) {
                    self.template operator()<_i + 1u>(self);
                }
            };
            set_index(set_index);
            return index;
        }

    };

    lolita::core::element::Element const static constexpr pnt_00 = lolita::core::element::Element{0, 0, 0, 1};
    //
    lolita::core::element::Element const static constexpr seg_02 = lolita::core::element::Element{1, 1, 1, 2};
    lolita::core::element::Element const static constexpr seg_03 = lolita::core::element::Element{5, 1, 2, 3};
    //
    lolita::core::element::Element const static constexpr tri_03 = lolita::core::element::Element{2, 2, 1, 3};
    lolita::core::element::Element const static constexpr tri_06 = lolita::core::element::Element{6, 2, 2, 6};
    lolita::core::element::Element const static constexpr qua_04 = lolita::core::element::Element{3, 2, 1, 4};
    lolita::core::element::Element const static constexpr qua_08 = lolita::core::element::Element{7, 2, 2, 8};
    //
    lolita::core::element::Element const static constexpr tet_04 = lolita::core::element::Element{4, 3, 1, 4};
    lolita::core::element::Element const static constexpr tet_12 = lolita::core::element::Element{8, 3, 2, 12};
    //
    // auto const static constexpr pg5_5 = Element{5, 2, 1, 5};

    /*
     *
     */

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain>
    concept CellConcept = _domain.dim_ - _element.dim_ == 0;

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain>
    concept FaceConcept = _domain.dim_ - _element.dim_ == 1;

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain>
    concept EdgeConcept = _domain.dim_ - _element.dim_ == 2;

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain>
    concept NodeConcept = _domain.dim_ - _element.dim_ == 3;

    /*
     *
     */

    template<lolita::core::element::Element _element>
    concept PointConcept = _element == pnt_00;

    template<lolita::core::element::Element _element>
    concept SegmentConcept = _element == seg_02 || _element == seg_03;

    template<lolita::core::element::Element _element>
    concept TriangleConcept = _element == tri_03 || _element == tri_06;

    template<lolita::core::element::Element _element>
    concept QuadrangleConcept = _element == qua_04 || _element == qua_08;

    template<lolita::core::element::Element _element>
    concept TetrahedronConcept = _element == tet_04 || _element == tet_12;

    /*
     *
     */

    template<lolita::core::element::Element _element>
    concept CurveConcept = SegmentConcept<_element>;

    template<lolita::core::element::Element _element>
    concept SurfaceConcept = TriangleConcept<_element> || QuadrangleConcept<_element>;

    template<lolita::core::element::Element _element>
    concept VolumeConcept = TetrahedronConcept<_element>;

    /*
     *
     */

    namespace detail
    {

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using Points = std::tuple<
                _T<pnt_00, _domain, _arg...>
        >;

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using Curves = std::tuple<
                _T<seg_02, _domain, _arg...>
        >;

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using Facets = std::tuple<
                _T<tri_03, _domain, _arg...>,
                _T<qua_04, _domain, _arg...>
        >;

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using Solids = std::tuple<
                _T<tet_04, _domain, _arg...>
        >;

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using Elements = std::tuple<
                Points<_T, _domain, _arg...>,
                Curves<_T, _domain, _arg...>,
                Facets<_T, _domain, _arg...>,
                Solids<_T, _domain, _arg...>
        >;

    }

    template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    using Elements = decltype(lolita::utility::tupleSlice<0, _domain.dim_ + 1>(std::declval<detail::Elements<_T, _domain, _arg...>>()));

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    struct ElementCollection
    {

    private:

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename __T, lolita::geometry::Domain __domain, auto... __arg>
        using _Elements = decltype(lolita::utility::tupleSlice<0, _domain.dim_ + 1>(std::declval<lolita::core::element::Elements<__T, __domain, __arg...>>()));

    public:

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

        template<lolita::core::element::Element _element, auto...>
        using ElementNodeConnectivity = std::array<lolita::index, _element.num_nodes_>;

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        struct ElementNeighboursTraits
        {

        private:

            template<lolita::core::element::Element __element, lolita::geometry::Domain __domain, auto... __arg>
            using _Neighbours = std::vector<_T<__element, __domain, __arg...>>;

            template<lolita::core::element::Element _element>
            struct _NeighboursTraits;

            template<lolita::core::element::Element _element>
            requires(PointConcept<_element>)
            struct _NeighboursTraits<_element>
            {

                using Type = decltype(lolita::utility::tupleSlice<1, _domain.dim_ + 1>(std::declval<lolita::core::element::Elements<_Neighbours, _domain, _arg...>>()));

            };

            template<lolita::core::element::Element _element>
            requires(!PointConcept<_element>)
            struct _NeighboursTraits<_element>
            {

                using Type = decltype(lolita::utility::tupleSlice<_element.dim_, _domain.dim_ + 1>(std::declval<lolita::core::element::Elements<_Neighbours, _domain, _arg...>>()));

            };

        public:

            template<lolita::core::element::Element _element>
            using Type = typename _NeighboursTraits<_element>::Type;

        };

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::core::element::Element _element, lolita::geometry::Domain _domain, auto... _arg>
        using ElementNeighbours = typename detail::ElementNeighboursTraits<_T, _domain, _arg...>::template Type<_element>;

    }

    template<lolita::core::element::Element, auto...>
    struct ElementGeometry;

    template<lolita::core::element::Element _element, auto... _arg>
    requires(_element == pnt_00)
    struct ElementGeometry<_element, _arg...>
    {

        lolita::core::element::Element const static constexpr element_ = _element;

//        constexpr
//        lolita::boolean
//        operator==(
//                ElementGeometry const & other
//        )
//        const = default;
//
//        constexpr
//        lolita::boolean
//        operator!=(
//                ElementGeometry const & other
//        )
//        const = default;

        /**
         * @brief
         */
        std::array<std::array<lolita::real, 3>, 1> const static constexpr reference_nodes_ = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

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

    template<lolita::core::element::Element _element, auto... _arg>
    requires(_element == seg_02)
    struct ElementGeometry<_element, _arg...>
    {

        lolita::core::element::Element const static constexpr element_ = _element;

//        constexpr
//        lolita::boolean
//        operator==(
//                ElementGeometry const & other
//        )
//        const = default;
//
//        constexpr
//        lolita::boolean
//        operator!=(
//                ElementGeometry const & other
//        )
//        const = default;

        std::array<std::array<lolita::real, 3>, 2> const static constexpr reference_nodes_ = {
                -1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<pnt_00, __arg...>, 2>
                >
        >;

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
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

    template<lolita::core::element::Element _element, auto... _arg>
    requires(_element == tri_03)
    struct ElementGeometry<_element, _arg...>
    {

        lolita::core::element::Element const static constexpr element_ = _element;

//        constexpr
//        lolita::boolean
//        operator==(
//                ElementGeometry const &
//                other
//        )
//        const = default;
//
//        constexpr
//        lolita::boolean
//        operator!=(
//                ElementGeometry const &
//                other
//        )
//        const = default;

        std::array<std::array<lolita::real, 3>, 3> const static constexpr reference_nodes_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
        };

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<seg_02, __arg...>, 3>
                >,
                std::tuple<
                        std::array<_T<pnt_00, __arg...>, 3>
                >
        >;

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

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

    template<lolita::core::element::Element _element, auto... _arg>
    requires(_element == qua_04)
    struct ElementGeometry<_element, _arg...>
    {

        lolita::core::element::Element const static constexpr element_ = _element;

//        constexpr
//        lolita::boolean
//        operator==(
//                ElementGeometry const &
//                other
//        )
//        const = default;
//
//        constexpr
//        lolita::boolean
//        operator!=(
//                ElementGeometry const &
//                other
//        )
//        const = default;

        std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                -1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                -1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
        };

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<seg_02, __arg...>, 4>
                >,
                std::tuple<
                        std::array<_T<pnt_00, __arg...>, 4>
                >
        >;

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

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

    template<lolita::core::element::Element _element, auto... _arg>
    requires(_element == tet_04)
    struct ElementGeometry<_element, _arg...>
    {

        lolita::core::element::Element const static constexpr element_ = _element;

//        constexpr
//        lolita::boolean
//        operator==(
//                ElementGeometry const &
//                other
//        )
//        const = default;
//
//        constexpr
//        lolita::boolean
//        operator!=(
//                ElementGeometry const &
//                other
//        )
//        const = default;

        std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000, +1.0000000000000000,
        };

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
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

        template<template<lolita::core::element::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = detail::ElementNeighbours<_T, _element, _domain, __arg...>;

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

    template<lolita::geometry::Domain _domain, lolita::core::element::Element _element>
    static constexpr
    ElementPosition
    elementPosition()
    {
        using _Elements = lolita::core::element::Elements<lolita::core::element::ElementGeometry, _domain>;
        auto position = lolita::core::element::ElementPosition{_element.dim_, -1};
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
        return std::tuple_element_t<J, std::tuple_element_t<I, lolita::core::element::Elements<lolita::core::element::ElementGeometry, _domain>>>::element_;
    }

    template<lolita::geometry::Domain _domain, lolita::index I>
    static constexpr
    lolita::index
    numElements()
    {
        return std::tuple_size_v<std::tuple_element_t<I, lolita::core::element::Elements<lolita::core::element::ElementGeometry, _domain>>>;
    }

    template<lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numElements()
    {
        return std::tuple_size_v<lolita::core::element::Elements<lolita::core::element::ElementGeometry, _domain>>;
    }

    /*
     *
     */

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, Element _component>
    static constexpr
    ElementPosition
    componentPosition()
    {
        using _Components = typename lolita::core::element::ElementGeometry<_element>::template Components<lolita::core::element::ElementGeometry, _domain>;
        auto position = lolita::core::element::ElementPosition{-1, -1};
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

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::index _I, lolita::index _J>
    static constexpr
    lolita::core::element::Element
    component()
    {
        using _Components = typename lolita::core::element::ElementGeometry<_element>::template Components<lolita::core::element::ElementGeometry, _domain>;
        return std::tuple_element_t<_J, std::tuple_element_t<_I, _Components>>::value_type::element_;
    }

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::index _I, lolita::index _J>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename lolita::core::element::ElementGeometry<_element>::template Components<lolita::core::element::ElementGeometry, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_J, std::tuple_element_t<_I, _Components>>>;
    }

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::index _I>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename lolita::core::element::ElementGeometry<_element>::template Components<lolita::core::element::ElementGeometry, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_I, _Components>>;
    }

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename lolita::core::element::ElementGeometry<_element>::template Components<lolita::core::element::ElementGeometry, _domain>;
        return std::tuple_size_v<_Components>;
    }

    /*
     *
     */

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, Element _neighbour>
    static constexpr
    ElementPosition
    neighbourPosition()
    {
        using _Neighbours = typename lolita::core::element::ElementGeometry<_element>::template Neighbours<lolita::core::element::ElementGeometry, _domain>;
        auto position = lolita::core::element::ElementPosition{-1, -1};
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

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::index _I, lolita::index _J>
    static constexpr
    lolita::core::element::Element
    neighbour()
    {
        using _Neighbours = typename lolita::core::element::ElementGeometry<_element>::template Neighbours<lolita::core::element::ElementGeometry, _domain>;
        return std::tuple_element_t<_J, std::tuple_element_t<_I, _Neighbours>>::value_type::element_;
    }

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain, lolita::index _I>
    static constexpr
    lolita::index
    numNeighbours()
    {
        using _Neighbours = typename lolita::core::element::ElementGeometry<_element>::template Neighbours<lolita::core::element::ElementGeometry, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_I, _Neighbours>>;
    }

    template<lolita::core::element::Element _element, lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numNeighbours()
    {
        using _Neighbours = typename lolita::core::element::ElementGeometry<_element>::template Neighbours<lolita::core::element::ElementGeometry, _domain>;
        return std::tuple_size_v<_Neighbours>;
    }

    /*
     *
     */

    template<lolita::core::element::Element, lolita::finite_element::Quadrature, lolita::index>
    struct ElementQuadrature;

    template<lolita::core::element::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord_quadrature>
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

    template<lolita::core::element::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord_quadrature>
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

}

#endif //LOLITA_LOLITA_ELEMENT_HXX
