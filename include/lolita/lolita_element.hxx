//
// Created by dsiedel on 05/05/22.
//

#ifndef LOLITA_LOLITA_ELEMENT_HXX
#define LOLITA_LOLITA_ELEMENT_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"

namespace lolita::core::base
{

    /**
     * @brief Basic structure to define an element
     */
    struct Element
    {

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        constexpr
        lolita::boolean
        operator==(
                Element const & other
        )
        const = default;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        constexpr
        lolita::boolean
        operator!=(
                Element const & other
        )
        const = default;

        /**
         * @brief The element tag, that fully defines it. Two elements cannot have the same tag
         */
        lolita::index tag_;

        /**
         * @brief The euclidean element dimension
         */
        lolita::index dim_;

        /**
         * @brief The element polynomial order
         */
        lolita::index ord_;

        /**
         * @brief The element number of nodes
         */
        lolita::index num_nodes_;

    };

    /**
     * @brief
     * @tparam _finite_element
     */
    template<lolita::finite_element::FiniteElementConcept auto... _finite_element>
    struct MixedElementPolicy
    {

        /**
         * @brief
         */
        template<
                template<lolita::core::base::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto> typename _T,
                lolita::core::base::Element _element,
                lolita::geometry::Domain _domain
        >
        using MixedElement = std::tuple<_T<_element, _domain, _finite_element>...>;

        /**
         * @brief Fetch the finite element index within the _finite_element list
         * @tparam __finite_element
         * @return the finite_element index if found, and the size of the _finite_element list otherwise
         */
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

    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr pnt_00 = lolita::core::base::Element{0, 0, 0, 1};
    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr seg_02 = lolita::core::base::Element{1, 1, 1, 2};
    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr seg_03 = lolita::core::base::Element{5, 1, 2, 3};
    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr tri_03 = lolita::core::base::Element{2, 2, 1, 3};
    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr tri_06 = lolita::core::base::Element{6, 2, 2, 6};
    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr qua_04 = lolita::core::base::Element{3, 2, 1, 4};
    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr qua_08 = lolita::core::base::Element{7, 2, 2, 8};
    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr tet_04 = lolita::core::base::Element{4, 3, 1, 4};
    /**
     * @brief
     */
    lolita::core::base::Element const static constexpr tet_12 = lolita::core::base::Element{8, 3, 2, 12};

    /*
     *
     */

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain>
    concept CellConcept = _domain.dim_ - _element.dim_ == 0;

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain>
    concept FaceConcept = _domain.dim_ - _element.dim_ == 1;

    /*
     *
     */

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept PointConcept = _element.dim_ == 0;

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept CurveConcept = _element.dim_ == 1;

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept FacetConcept = _element.dim_ == 2;

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept SolidConcept = _element.dim_ == 3;

    /*
     *
     */

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept NodeConcept = _element == pnt_00;

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept SegmentConcept = _element == seg_02 || _element == seg_03;

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept TriangleConcept = _element == tri_03 || _element == tri_06;

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept QuadrangleConcept = _element == qua_04 || _element == qua_08;

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core::base::Element _element>
    concept TetrahedronConcept = _element == tet_04 || _element == tet_12;

    /*
     *
     */

    namespace detail
    {

        template<lolita::core::base::Element _element, lolita::geometry::Domain _domain, auto... _args>
        struct _Span
        {

            lolita::core::base::Element const static constexpr element_ = _element;

        };

        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using _Points = std::tuple<
                _T<pnt_00, _domain, _arg...>
        >;

        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using _Curves = std::tuple<
                _T<seg_02, _domain, _arg...>
        >;

        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using _Facets = std::tuple<
                _T<tri_03, _domain, _arg...>,
                _T<qua_04, _domain, _arg...>
        >;

        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using _Solids = std::tuple<
                _T<tet_04, _domain, _arg...>
        >;

        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        using _Elements = std::tuple<
                _Points<_T, _domain, _arg...>,
                _Curves<_T, _domain, _arg...>,
                _Facets<_T, _domain, _arg...>,
                _Solids<_T, _domain, _arg...>
        >;

    }

    /**
     * @brief
     */
    template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
    using Elements = decltype(lolita::utility::tupleSlice<0, _domain.dim_ + 1>(std::declval<lolita::core::base::detail::_Elements<_T, _domain, _arg...>>()));

    /**
     * @brief Holder for the position of an element in a mesh, or with respect to some other element
     */
    struct ElementCoordinates
    {

        /**
         * @brief The relative dimension of the element in the mesh, or with respect to some other element
         */
        lolita::integer dim_;

        /**
         * @brief The relative position of the element in the mesh, or with respect to some other element
         */
        lolita::integer tag_;

    };

    /**
     *
     * @tparam _domain
     * @tparam _element
     * @return
     */
    template<lolita::geometry::Domain _domain, lolita::core::base::Element _element>
    static constexpr
    lolita::core::base::ElementCoordinates
    elementPosition()
    {
        using _Elements = lolita::core::base::Elements<lolita::core::base::detail::_Span, _domain>;
        auto position = lolita::core::base::ElementCoordinates{_element.dim_, -1};
        auto f0 = [&] <lolita::index _i = 0u> (auto & self) constexpr mutable {
            using _Element = std::tuple_element_t<_i, std::tuple_element_t<_element.dim_, _Elements>>;
            if (_Element::element_ == _element) {
                position.tag_ = _i;
            }
            if constexpr (_i < std::tuple_size_v<std::tuple_element_t<_element.dim_, _Elements>> - 1) {
                self.template operator()<_i + 1u>(self);
            }
        };
        f0(f0);
        return position;
    }

    /**
     * @brief
     * @tparam _domain
     * @tparam _i
     * @tparam _j
     * @return
     */
    template<lolita::geometry::Domain _domain, lolita::index _i, lolita::index _j>
    static constexpr
    lolita::core::base::Element
    element()
    {
        return std::tuple_element_t<_j, std::tuple_element_t<_i, lolita::core::base::Elements<lolita::core::base::detail::_Span, _domain>>>::element_;
    }

    /**
     * @brief
     * @tparam _domain
     * @tparam _i
     * @return
     */
    template<lolita::geometry::Domain _domain, lolita::index _i>
    static constexpr
    lolita::index
    numElements()
    {
        return std::tuple_size_v<std::tuple_element_t<_i, lolita::core::base::Elements<lolita::core::base::detail::_Span, _domain>>>;
    }

    /**
     * @brief
     * @tparam _domain
     * @return
     */
    template<lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numElements()
    {
        return std::tuple_size_v<lolita::core::base::Elements<lolita::core::base::detail::_Span, _domain>>;
    }

    namespace detail
    {

        template<lolita::core::base::Element _element, auto...>
        using _ElementNodeConnectivity = std::array<lolita::index, _element.num_nodes_>;

        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        struct _ElementNeighboursPolicy
        {

        private:

            template<lolita::core::base::Element __element, lolita::geometry::Domain __domain, auto... __arg>
            using __Neighbours = std::vector<_T<__element, __domain, __arg...>>;

            template<lolita::core::base::Element _element>
            struct _NeighboursTraits;

            template<lolita::core::base::Element _element>
            requires(lolita::core::base::PointConcept<_element>)
            struct _NeighboursTraits<_element>
            {

                using HHH = lolita::core::base::Elements<__Neighbours, _domain, _arg...>;

                using Neighbours = decltype(lolita::utility::tupleSlice<1, _domain.dim_ + 1>(std::declval<HHH>()));

            };

            template<lolita::core::base::Element _element>
            requires(!lolita::core::base::PointConcept<_element>)
            struct _NeighboursTraits<_element>
            {

                using HHH = lolita::core::base::Elements<__Neighbours, _domain, _arg...>;

                using Neighbours = decltype(lolita::utility::tupleSlice<_element.dim_, _domain.dim_ + 1>(std::declval<HHH>()));

            };

        public:

            template<lolita::core::base::Element _element>
            using _Neighbours = typename _NeighboursTraits<_element>::Neighbours;

        };

        template<
                template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T,
                lolita::core::base::Element _element,
                lolita::geometry::Domain _domain,
                auto... _arg
        >
        using _ElementNeighbours = typename detail::_ElementNeighboursPolicy<_T, _domain, _arg...>::template _Neighbours<_element>;

    }

    /**
     * @brief
     * @tparam ...
     */
    template<lolita::core::base::Element, auto...>
    struct ElementGeometry;

    /**
     * @brief Element description of the point
     * @tparam _element
     * @tparam _arg
     */
    template<lolita::core::base::Element _element, auto... _arg>
    requires(_element == lolita::core::base::pnt_00)
    struct ElementGeometry<_element, _arg...>
    {

        /**
         * @brief The element object
         */
//        lolita::core::base::Element const static constexpr element_ = _element;

        /**
         * @brief Reference nodes for the reference element
         */
        std::array<std::array<lolita::real, 3>, 1> const static constexpr reference_nodes_ = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        /**
         * @brief Sub-elements composition
         */
        using Components = lolita::geometry::Point;

        /**
         * @brief Sup-elements composition
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = lolita::core::base::detail::_ElementNeighbours<_T, _element, _domain, __arg...>;

        /**
         * @brief Shape function evaluation given some scalar nodal values
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
            return nodal_field_values(0);
        }

        /**
         * @brief Shape function derivative given some scalar nodal values
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
            return lolita::real(0);
        }

    };

    /**
     * @brief
     * @tparam _element
     * @tparam _arg
     */
    template<lolita::core::base::Element _element, auto... _arg>
    requires(_element == lolita::core::base::seg_02)
    struct ElementGeometry<_element, _arg...>
    {

        /**
         * @brief The element object
         */
//        lolita::core::base::Element const static constexpr element_ = _element;

        /**
         * @brief Reference nodes for the reference element
         */
        std::array<std::array<lolita::real, 3>, 2> const static constexpr reference_nodes_ = {
                -1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        /**
         * @brief Sub-elements composition
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<lolita::core::base::pnt_00, __arg...>, 2>
                >
        >;

        /**
         * @brief Sup-elements composition
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = lolita::core::base::detail::_ElementNeighbours<_T, _element, _domain, __arg...>;

        /**
         * @brief Node connexions with sub-elements
         */
        Components<detail::_ElementNodeConnectivity> const static constexpr node_connectivity_ = {
                {
                        {
                                0,
                                1,
                        }
                }
        };

        /**
         * @brief Shape function evaluation given some scalar nodal values
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
            auto value = lolita::real(0);
            value += nodal_field_values(0) * (1.0 / 2.0) * (1.0 - reference_point(0));
            value += nodal_field_values(1) * (1.0 / 2.0) * (1.0 + reference_point(0));
            return value;
        }

        /**
         * @brief Shape function derivative given some scalar nodal values
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
            assert(derivative_direction == 0);
            auto value = lolita::real(0);
            value += - nodal_field_values(0) * (1.0 / 2.0);
            value += + nodal_field_values(1) * (1.0 / 2.0);
            return value;
        }

    };

    /**
     * @brief
     * @tparam _element
     * @tparam _arg
     */
    template<lolita::core::base::Element _element, auto... _arg>
    requires(_element == lolita::core::base::tri_03)
    struct ElementGeometry<_element, _arg...>
    {

        /**
         * @brief The element object
         */
//        lolita::core::base::Element const static constexpr element_ = _element;

        /**
         * @brief Reference nodes for the reference element
         */
        std::array<std::array<lolita::real, 3>, 3> const static constexpr reference_nodes_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
        };

        /**
         * @brief Sub-elements composition
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<lolita::core::base::seg_02, __arg...>, 3>
                >,
                std::tuple<
                        std::array<_T<lolita::core::base::pnt_00, __arg...>, 3>
                >
        >;

        /**
         * @brief Sup-elements composition
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = lolita::core::base::detail::_ElementNeighbours<_T, _element, _domain, __arg...>;

        /**
         * @brief Node connexions with sub-elements
         */
        Components<detail::_ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
         * @brief Shape function evaluation given some scalar nodal values
         * @param nodal_field_values nodal scalar values as a vector
         * @param reference_point a point of the reference element as a vector
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point
        )
        {
            auto value = lolita::real(0);
            value += nodal_field_values(0) * (1.0 - reference_point(0) - reference_point(1));
            value += nodal_field_values(1) * reference_point(0);
            value += nodal_field_values(2) * reference_point(1);
            return value;
        }

        /**
         * @brief Shape function derivative given some scalar nodal values
         * @param nodal_field_values the nodal scalar values of some field
         * @param reference_point the point in the reference element where to evaluate the field derivative
         * @param derivative_direction the derivative direction, that must be less than or equal to the element dimension
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

    /**
     * @brief
     * @tparam _element
     * @tparam _arg
     */
    template<lolita::core::base::Element _element, auto... _arg>
    requires(_element == lolita::core::base::qua_04)
    struct ElementGeometry<_element, _arg...>
    {

        /**
         * @brief The element object
         */
//        lolita::core::base::Element const static constexpr element_ = _element;

        /**
         * @brief Reference nodes for the reference element
         */
        std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                -1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                -1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
        };

        /**
         * @brief Sub-elements composition
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<lolita::core::base::seg_02, __arg...>, 4>
                >,
                std::tuple<
                        std::array<_T<lolita::core::base::pnt_00, __arg...>, 4>
                >
        >;

        /**
         * @brief Sup-elements composition
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = lolita::core::base::detail::_ElementNeighbours<_T, _element, _domain, __arg...>;

        /**
         * @brief Node connexions with sub-elements
         */
        Components<detail::_ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point,
                lolita::index derivative_direction
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

    /**
     * @brief
     * @tparam _element
     * @tparam _arg
     */
    template<lolita::core::base::Element _element, auto... _arg>
    requires(_element == lolita::core::base::tet_04)
    struct ElementGeometry<_element, _arg...>
    {

        /**
         * @brief
         */
//        lolita::core::base::Element const static constexpr element_ = _element;

        /**
         * @brief
         */
        std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000, +1.0000000000000000,
        };

        /**
         * @brief
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, auto... __arg>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<lolita::core::base::tri_03, __arg...>, 4>
                >,
                std::tuple<
                        std::array<_T<lolita::core::base::seg_02, __arg...>, 6>
                >,
                std::tuple<
                        std::array<_T<lolita::core::base::pnt_00, __arg...>, 4>
                >
        >;

        /**
         * @brief
         */
        template<template<lolita::core::base::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... __arg>
        using Neighbours = lolita::core::base::detail::_ElementNeighbours<_T, _element, _domain, __arg...>;

        /**
         * @brief
         */
        Components<detail::_ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point
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
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point,
                lolita::index derivative_direction
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

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _component
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain, lolita::core::base::Element _component>
    static constexpr
    lolita::core::base::ElementCoordinates
    componentPosition()
    {
        using _Components = typename lolita::core::base::ElementGeometry<_element>::template Components<lolita::core::base::detail::_Span, _domain>;
        auto position = lolita::core::base::ElementCoordinates{-1, -1};
        auto f0 = [&] <lolita::index _i = 0u, lolita::index _j = 0u> (
                auto & self
        )
        constexpr mutable
        {
            using _Component = typename std::tuple_element_t<_j, std::tuple_element_t<_i, _Components>>::value_type;
            if (_Component::element_ == _component) {
                position.dim_ = _i;
                position.tag_ = _j;
            }
            if constexpr (_j < std::tuple_size_v<std::tuple_element_t<_i, _Components>> - 1) {
                self.template operator()<_i, _j + 1u>(self);
            }
            else if constexpr (_i < std::tuple_size_v<_Components> - 1) {
                self.template operator()<_i + 1u, 0u>(self);
            }
        };
        f0(f0);
        return position;
    }

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _i
     * @tparam _j
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain, lolita::index _i, lolita::index _j>
    static constexpr
    lolita::core::base::Element
    component()
    {
        using _Components = typename lolita::core::base::ElementGeometry<_element>::template Components<lolita::core::base::detail::_Span, _domain>;
        return std::tuple_element_t<_j, std::tuple_element_t<_i, _Components>>::value_type::element_;
    }

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _i
     * @tparam _j
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain, lolita::index _i, lolita::index _j>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename lolita::core::base::ElementGeometry<_element>::template Components<lolita::core::base::detail::_Span, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_j, std::tuple_element_t<_i, _Components>>>;
    }

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _i
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain, lolita::index _i>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename lolita::core::base::ElementGeometry<_element>::template Components<lolita::core::base::detail::_Span, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_i, _Components>>;
    }

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numComponents()
    {
        using _Components = typename lolita::core::base::ElementGeometry<_element>::template Components<lolita::core::base::detail::_Span, _domain>;
        return std::tuple_size_v<_Components>;
    }

    /*
     *
     */

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _neighbour
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain, lolita::core::base::Element _neighbour>
    static constexpr
    lolita::core::base::ElementCoordinates
    neighbourPosition()
    {
        using _Neighbours = typename lolita::core::base::ElementGeometry<_element>::template Neighbours<lolita::core::base::detail::_Span, _domain>;
        auto position = lolita::core::base::ElementCoordinates{-1, -1};
        auto f0 = [&] <lolita::index _i = 0u, lolita::index _j = 0u> (
                auto & self
        )
        constexpr mutable
        {
            using _Neighbour = typename std::tuple_element_t<_j, std::tuple_element_t<_i, _Neighbours>>::value_type;
            if (_Neighbour::element_ == _neighbour) {
                position.dim_ = _i;
                position.tag_ = _j;
            }
            if constexpr (_j < std::tuple_size_v<std::tuple_element_t<_i, _Neighbours>> - 1) {
                self.template operator()<_i, _j + 1u>(self);
            }
            else if constexpr (_i < std::tuple_size_v<_Neighbours> - 1) {
                self.template operator()<_i + 1u, 0u>(self);
            }
        };
        f0(f0);
        return position;
    }

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _i
     * @tparam _j
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain, lolita::index _i, lolita::index _j>
    static constexpr
    lolita::core::base::Element
    neighbour()
    {
        using _Neighbours = typename lolita::core::base::ElementGeometry<_element>::template Neighbours<lolita::core::base::detail::_Span, _domain>;
        return std::tuple_element_t<_j, std::tuple_element_t<_i, _Neighbours>>::value_type::element_;
    }

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _i
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain, lolita::index _i>
    static constexpr
    lolita::index
    numNeighbours()
    {
        using _Neighbours = typename lolita::core::base::ElementGeometry<_element>::template Neighbours<lolita::core::base::detail::_Span, _domain>;
        return std::tuple_size_v<std::tuple_element_t<_i, _Neighbours>>;
    }

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @return
     */
    template<lolita::core::base::Element _element, lolita::geometry::Domain _domain>
    static constexpr
    lolita::index
    numNeighbours()
    {
        using _Neighbours = typename lolita::core::base::ElementGeometry<_element>::template Neighbours<lolita::core::base::detail::_Span, _domain>;
        return std::tuple_size_v<_Neighbours>;
    }

    /*
     *
     */

    /**
     * @brief
     * @tparam _element
     * @tparam _quadrature
     * @tparam _ord
     */
    template<lolita::core::base::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
    struct ElementQuadrature;

    /**
     * @brief
     * @tparam _quadrature
     * @tparam _ord
     */
    template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
    struct ElementQuadratureBase
    {

        /**
         * @brief
         */
        lolita::finite_element::Quadrature const static constexpr quadrature_ = _quadrature;

        /**
         * @brief
         */
        lolita::index const static constexpr ord_ = _ord;

    };

    /**
     * @brief
     * @tparam _element
     * @tparam _quadrature
     * @tparam _ord_quadrature
     */
    template<lolita::core::base::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
    requires(lolita::core::base::PointConcept<_element>)
    struct ElementQuadrature<_element, _quadrature, _ord> : public lolita::core::base::ElementQuadratureBase<_quadrature, _ord>
    {

        /**
         * @brief
         */
        lolita::index const static constexpr dim_ = 1;

        /**
         * @brief
         */
        std::array<std::array<lolita::real, 3>, dim_> const static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        /**
         * @brief
         */
        std::array<lolita::real, dim_> const static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };

    /**
     * @brief
     * @tparam _element
     * @tparam _quadrature
     * @tparam _ord_quadrature
     */
    template<lolita::core::base::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
    requires(!lolita::core::base::PointConcept<_element>)
    struct ElementQuadrature<_element, _quadrature, _ord> : public lolita::core::base::ElementQuadratureBase<_quadrature, _ord>
    {

        /**
         * @brief
         */
        lolita::index const static constexpr dim_ = 1;

        /**
         * @brief
         */
        std::array<std::array<lolita::real, 3>, dim_> const static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        /**
         * @brief
         */
        std::array<lolita::real, dim_> const static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };

}

#endif //LOLITA_LOLITA_ELEMENT_HXX
