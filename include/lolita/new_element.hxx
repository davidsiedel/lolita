//
// Created by dsiedel on 28/05/22.
//

#ifndef LOLITA_NEW_ELEMENT_HXX
#define LOLITA_NEW_ELEMENT_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/new.hxx"

namespace lolita::core::element
{
    
    namespace detail
    {



        /**
         * @brief
         * @tparam _T
         * @tparam _element
         * @tparam _domain
         * @tparam _args
         */
        template<
                template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T,
                lolita::core::Element _element,
                lolita::geometry::Domain _domain,
                auto... _args
        >
        struct _ElementNeighboursPolicy
        {

        private:

            template<lolita::core::Element __element, lolita::geometry::Domain __domain, auto... __args>
            using __Neighbours = std::vector<_T<__element, __domain, __args...>>;

            template<lolita::core::Element __element>
            struct _NeighboursTraits;

            template<lolita::core::Element __element>
            requires(_element == lolita::core::Element::Node())
            struct _NeighboursTraits<__element>
            {

                using HHH = lolita::core::detail::_Elements<__Neighbours, _domain, _args...>;

                using Neighbours = decltype(lolita::utility::tupleSlice<1, _domain.dim_ + 1>(std::declval<HHH>()));

            };

            template<lolita::core::Element __element>
            requires(_element != lolita::core::Element::Node())
            struct _NeighboursTraits<__element>
            {

                using HHH = lolita::core::detail::_Elements<__Neighbours, _domain, _args...>;

                using Neighbours = decltype(lolita::utility::tupleSlice<__element.dim_, _domain.dim_ + 1>(std::declval<HHH>()));

            };

        public:

            using _Neighbours = typename _NeighboursTraits<_element>::Neighbours;

        };

        /**
         * @brief
         */
        template<
                template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T,
                lolita::core::Element _element,
                lolita::geometry::Domain _domain,
                auto... _args
        >
        using _ElementNeighbours = typename detail::_ElementNeighboursPolicy<_T, _element, _domain, _args...>::_Neighbours;
        
    }

    /**
     * @brief
     */
    template<lolita::core::Element _element, auto...>
    using ElementNodeConnectivity = std::array<lolita::index, _element.num_nodes_>;


    template<lolita::core::Element _element, lolita::geometry::Domain _domain>
    requires(_element == lolita::core::Element::Node())
    struct ElementGeometry<_element, _domain>
    {

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
        template<template<lolita::core::Element __element, lolita::geometry::Domain __domain, auto...> typename _T, auto... __args>
        using Neighbours = lolita::core::element::detail::_ElementNeighbours<_T, _element, _domain, __args...>;

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

    template<lolita::core::Element _element, lolita::geometry::Domain _domain>
    requires(_element == lolita::core::Element::LinearSegment())
    struct ElementGeometry<_element, _domain>
    {

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
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, auto... __args>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<lolita::core::Element::Node(), _domain, __args...>, 2>
                >
        >;

        /**
         * @brief Sup-elements composition
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, auto... __args>
        using Neighbours = lolita::core::element::detail::_ElementNeighbours<_T, _element, _domain, __args...>;

        /**
         * @brief Node connexions with sub-elements
         */
        Components<lolita::core::element::ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
     * @tparam _domain 
     * @tparam _args 
     */
    template<lolita::core::Element _element, lolita::geometry::Domain _domain>
    requires(_element == lolita::core::Element::LinearTriangle())
    struct ElementGeometry<_element, _domain>
    {

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
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, auto... __args>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<lolita::core::Element::LinearSegment(), _domain, __args...>, 3>
                >,
                std::tuple<
                        std::array<_T<lolita::core::Element::Node(), _domain, __args...>, 3>
                >
        >;

        /**
         * @brief Sup-elements composition
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, auto... __args>
        using Neighbours = lolita::core::element::detail::_ElementNeighbours<_T, _element, _domain, __args...>;

        /**
         * @brief Node connexions with sub-elements
         */
        Components<lolita::core::element::ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
     * @tparam _domain 
     * @tparam _args 
     */
    template<lolita::core::Element _element, lolita::geometry::Domain _domain>
    requires(_element == lolita::core::Element::LinearQuadrangle())
    struct ElementGeometry<_element, _domain>
    {

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
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, auto... __args>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<lolita::core::Element::LinearSegment(), _domain, __args...>, 4>
                >,
                std::tuple<
                        std::array<_T<lolita::core::Element::Node(), _domain, __args...>, 4>
                >
        >;

        /**
         * @brief Sup-elements composition
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, auto... __args>
        using Neighbours = lolita::core::element::detail::_ElementNeighbours<_T, _element, _domain, __args...>;

        /**
         * @brief Node connexions with sub-elements
         */
        Components<lolita::core::element::ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
     * @tparam _domain 
     * @tparam _args 
     */
    template<lolita::core::Element _element, lolita::geometry::Domain _domain>
    requires(_element == lolita::core::Element::LinearTetrahedron())
    struct ElementGeometry<_element, _domain>
    {

        /**
         * @brief
         */
//        lolita::core::Element const static constexpr element_ = _element;

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
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, auto... __args>
        using Components = std::tuple<
                std::tuple<
                        std::array<_T<lolita::core::Element::LinearTriangle(), _domain, __args...>, 4>
                >,
                std::tuple<
                        std::array<_T<lolita::core::Element::LinearSegment(), _domain, __args...>, 6>
                >,
                std::tuple<
                        std::array<_T<lolita::core::Element::Node(), _domain, __args...>, 4>
                >
        >;

        /**
         * @brief
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, auto... __args>
        using Neighbours = lolita::core::element::detail::_ElementNeighbours<_T, _element, _domain, __args...>;

        /**
         * @brief
         */
        Components<lolita::core::element::ElementNodeConnectivity> const static constexpr node_connectivity_ = {
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
     * @tparam _quadrature
     * @tparam _ord
     */
    template<lolita::core::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
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
    template<lolita::core::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
    requires(_element.isPoint())
    struct ElementQuadrature<_element, _quadrature, _ord> : public lolita::core::element::ElementQuadratureBase<_quadrature, _ord>
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
    template<lolita::core::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
    requires(!_element.isPoint())
    struct ElementQuadrature<_element, _quadrature, _ord> : public lolita::core::element::ElementQuadratureBase<_quadrature, _ord>
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

#endif //LOLITA_NEW_ELEMENT_HXX
