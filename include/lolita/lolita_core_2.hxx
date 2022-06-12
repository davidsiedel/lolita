//
// Created by dsiedel on 04/06/22.
//

#ifndef LOLITA_LOLITA_CORE_2_HXX
#define LOLITA_LOLITA_CORE_2_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"

namespace lolita::core2::geometry
{

    /**
     * @brief Holder for the position of an element in a mesh, or with respect to some other element
     */
    struct ElementCoordinates
    {

        /**
         * @brief The relative dimension of the element in the mesh, or with respect to some other element
         */
        lolita::index dim_;

        /**
         * @brief The relative position of the element in the mesh, or with respect to some other element
         */
        lolita::index tag_;

        /**
         * @brief Simple getter
         * @return The relative dimension of the element in the mesh, or with respect to some other element
         */
        constexpr
        lolita::index
        dim()
        const
        {
            return dim_;
        }

        /**
         * @brief Simple getter
         * @return The relative position of the element in the mesh, or with respect to some other element
         */
        constexpr
        lolita::index
        tag()
        const
        {
            return tag_;
        }

        friend
        std::ostream &
        operator<<(
                std::ostream & os,
                ElementCoordinates const & coordinates
        )
        {
            os << "c0 : " << coordinates.dim_ << ", c1 : " << coordinates.tag_;
            return os;
        }

    };

    /**
     * @brief Basic structure to define an element
     */
    struct Element : public lolita::utility::Enumeration<Element>
    {

        /**
         * @brief
         * @param tag
         * @param dim
         * @param ord
         * @param num_nodes
         */
        constexpr
        Element(
                std::basic_string_view<lolita::character> && tag,
                lolita::integer dim,
                lolita::integer ord,
                lolita::integer num_nodes
        )
        :
        lolita::utility::Enumeration<Element>(std::forward<std::basic_string_view<lolita::character>>(tag)),
        dim_(dim),
        ord_(ord),
        num_nodes_(num_nodes)
        {}

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core2::geometry::Element
        Node()
        {
            return Element("Node", 0, 0, 1);
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isNode()
        const
        {
            return * this == Node();
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core2::geometry::Element
        LinearSegment()
        {
            return Element("LinearSegment", 1, 1, 2);
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isLinearSegment()
        const
        {
            return * this == LinearSegment();
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core2::geometry::Element
        LinearTriangle()
        {
            return Element("LinearTriangle", 2, 1, 3);
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isLinearTriangle()
        const
        {
            return * this == LinearTriangle();
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core2::geometry::Element
        LinearQuadrangle()
        {
            return Element("LinearQuadrangle", 2, 1, 4);
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isLinearQuadrangle()
        const
        {
            return * this == LinearQuadrangle();
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core2::geometry::Element
        LinearTetrahedron()
        {
            return Element("LinearTetrahedron", 3, 1, 4);
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isLinearTetrahedron()
        const
        {
            return * this == LinearTetrahedron();
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isPoint()
        const
        {
            return dim_ == 0;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isCurve()
        const
        {
            return dim_ == 1;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isFacet()
        const
        {
            return dim_ == 2;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isSolid()
        const
        {
            return dim_ == 3;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isSegment()
        const
        {
            return isLinearSegment();
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isTriangle()
        const
        {
            return isLinearTriangle();
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isQuadrangle()
        const
        {
            return isLinearQuadrangle();
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isTetrahedron()
        const
        {
            return isLinearTetrahedron();
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::integer
        dim()
        const
        {
            return dim_;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::integer
        ord()
        const
        {
            return ord_;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::integer
        numNodes()
        const
        {
            return num_nodes_;
        }

        /**
         * @brief
         * @param ord
         * @return
         */
        constexpr
        lolita::boolean
        hasOrd(
                lolita::integer ord
        )
        const
        {
            return ord_ == ord;
        }

        /**
         * @brief
         * @param domain
         * @param level
         * @return
         */
        constexpr
        lolita::boolean
        isSub(
                lolita::domain::Domain domain,
                lolita::integer level
        )
        const
        {
            return domain.dim() - dim_ == level;
        }

        /**
         * @brief The euclidean element dimension
         */
        lolita::integer const dim_;

        /**
         * @brief The element polynomial order
         */
        lolita::integer const ord_;

        /**
         * @brief The element number of nodes
         */
        lolita::integer const num_nodes_;

    };

    /**
     * @brief
     */
    template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, lolita::domain::Domain t_domain, auto... t_args>
    using Points = std::tuple<
            t_T<lolita::core2::geometry::Element::Node(), t_domain, t_args...>
    >;

    /**
     * @brief
     */
    template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, lolita::domain::Domain t_domain, auto... t_args>
    using Curves = std::tuple<
            t_T<lolita::core2::geometry::Element::LinearSegment(), t_domain, t_args...>
    >;

    /**
     * @brief
     */
    template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, lolita::domain::Domain t_domain, auto... t_args>
    using Facets = std::tuple<
            t_T<lolita::core2::geometry::Element::LinearTriangle(), t_domain, t_args...>,
            t_T<lolita::core2::geometry::Element::LinearQuadrangle(), t_domain, t_args...>
    >;

    /**
     * @brief
     */
    template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, lolita::domain::Domain t_domain, auto... t_args>
    using Solids = std::tuple<
            t_T<lolita::core2::geometry::Element::LinearTetrahedron(), t_domain, t_args...>
    >;

    namespace detail
    {

        template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, lolita::domain::Domain t_domain, auto... t_args>
        using Elements = std::tuple<
                Points<t_T, t_domain, t_args...>,
                Curves<t_T, t_domain, t_args...>,
                Facets<t_T, t_domain, t_args...>,
                Solids<t_T, t_domain, t_args...>
        >;

    }

    /**
     * @brief
     */
    template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, lolita::domain::Domain t_domain, auto... t_args>
    using Elements = lolita::utility::tuple_slice_t<detail::Elements<t_T, t_domain, t_args...>, 0, t_domain.dim() + 1>;

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
    struct ElementTraits;

    namespace detail
    {

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        struct ElementView
        {

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::core2::geometry::Element
            getElement()
            {
                return t_element;
            }

        };

        /**
         * @brief Forward declaration
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        struct ElementOuterGeometryTraits;

        /**
         * @brief Forward declaration
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        struct ElementInnerGeometryTraits;

        /**
         * @brief
         */
        template<lolita::core2::geometry::Element t_element, auto...>
        using ElementNodeConnectivity = std::array<lolita::index, t_element.numNodes()>;

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        requires (t_element.isPoint())
        struct ElementOuterGeometryTraits<t_element, t_domain>
        {

        private:

            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            struct OuterConnectivityTraits {

            private:

                template<lolita::core2::geometry::Element t__element, lolita::domain::Domain t__domain, auto... t__args>
                using NeighbourVector = std::vector<t_T<t__element, t__domain, t__args...>>;

            public:

                using OuterConnectivity = lolita::utility::tuple_slice_t<Elements<NeighbourVector, t_domain, t_args...>, 1, t_domain.dim() + 1>;

            };

        public:

            /**
             * @brief
             */
            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            using OuterConnectivity = typename OuterConnectivityTraits<t_T, t_args...>::OuterConnectivity;

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        requires (!t_element.isPoint())
        struct ElementOuterGeometryTraits<t_element, t_domain>
        {

        private:

            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            struct OuterConnectivityTraits {

            private:

                template<lolita::core2::geometry::Element t__element, lolita::domain::Domain t__domain, auto... t__args>
                using NeighbourVector = std::vector<t_T<t__element, t__domain, t__args...>>;

            public:

                using OuterConnectivity = lolita::utility::tuple_slice_t<Elements<NeighbourVector, t_domain, t_args...>, t_element.dim(), t_domain.dim() + 1>;

            };

        public:

            /**
             * @brief
             */
            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            using OuterConnectivity = typename OuterConnectivityTraits<t_T, t_args...>::OuterConnectivity;

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        requires (t_element.isPoint())
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {

            /**
             * @brief
             */
            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<>;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:

            /**
             * @brief function evaluation given some scalar nodal values
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point
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
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point,
                    lolita::index derivative_direction
            )
            {
                return lolita::real(0);
            }

            /**
             * @brief
             */
            NodeConnectivity const static constexpr node_connectivity_ = NodeConnectivity{};

            /**
             * @brief Reference nodes for the reference element
             */
            std::array<std::array<lolita::real, 3>, 1> const static constexpr reference_nodes_ = {
                    +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
            };

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        requires (t_element.isLinearSegment())
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {

            /**
             * @brief
             */
            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<
                    std::tuple<
                            std::array<t_T<lolita::core2::geometry::Element::Node(), t_domain, t_args...>, 2>
                    >
            >;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:

            /**
             * @brief Shape function evaluation given some scalar nodal values
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point
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
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point,
                    lolita::index derivative_direction
            )
            {
                assert(derivative_direction == 0);
                auto value = lolita::real(0);
                value += -nodal_field_values(0) * (1.0 / 2.0);
                value += +nodal_field_values(1) * (1.0 / 2.0);
                return value;
            }

            /**
             * @brief
             */
            NodeConnectivity const static constexpr node_connectivity_ = NodeConnectivity{
                    {
                            {
                                    0,
                                    1,
                            }
                    }
            };

            /**
             * @brief Reference nodes for the reference element
             */
            std::array<std::array<lolita::real, 3>, 2> const static constexpr reference_nodes_ = {
                    -1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, +0.0000000000000000, +0.0000000000000000
            };

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        requires (t_element.isLinearTriangle())
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {

            /**
             * @brief
             */
            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<
                    std::tuple<
                            std::array<t_T<lolita::core2::geometry::Element::LinearSegment(), t_domain, t_args...>, 3>
                    >,
                    std::tuple<
                            std::array<t_T<lolita::core2::geometry::Element::Node(), t_domain, t_args...>, 3>
                    >
            >;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:

            /**
             * @brief Shape function evaluation given some scalar nodal values
             * @param nodal_field_values nodal scalar values as a vector
             * @param reference_point a point of the reference element as a vector
             * @return
             */
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point
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
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point,
                    lolita::index derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 1);
                auto value = lolita::real(0);
                if (derivative_direction == 0) {
                    value += -nodal_field_values(0);
                    value += +nodal_field_values(1);
                    value += +0.0;
                } else {
                    value += -nodal_field_values(0);
                    value += +0.0;
                    value += +nodal_field_values(2);
                }
                return value;
            }

            NodeConnectivity const static constexpr node_connectivity_ = NodeConnectivity{
                    {
                            {
                                    0, 1,
                                    1, 2,
                                    2, 0,
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
             * @brief Reference nodes for the reference element
             */
            std::array<std::array<lolita::real, 3>, 3> const static constexpr reference_nodes_ = {
                    +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
            };

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        requires (t_element.isLinearQuadrangle())
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {

            /**
             * @brief
             */
            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<
                    std::tuple<
                            std::array<t_T<lolita::core2::geometry::Element::LinearSegment(), t_domain, t_args...>, 4>
                    >,
                    std::tuple<
                            std::array<t_T<lolita::core2::geometry::Element::Node(), t_domain, t_args...>, 4>
                    >
            >;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point
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
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point,
                    lolita::index derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 1);
                auto value = lolita::real(0);
                if (derivative_direction == 0) {
                    value += -nodal_field_values(0) * (1.0 / 4.0) * (1.0 - reference_point(1));
                    value += +nodal_field_values(1) * (1.0 / 4.0) * (1.0 - reference_point(1));
                    value += +nodal_field_values(2) * (1.0 / 4.0) * (1.0 + reference_point(1));
                    value += -nodal_field_values(3) * (1.0 / 4.0) * (1.0 + reference_point(1));
                }
                else {
                    value += -nodal_field_values(0) * (1.0 / 4.0) * (1.0 - reference_point(0));
                    value += -nodal_field_values(1) * (1.0 / 4.0) * (1.0 + reference_point(0));
                    value += +nodal_field_values(2) * (1.0 / 4.0) * (1.0 + reference_point(0));
                    value += +nodal_field_values(3) * (1.0 / 4.0) * (1.0 - reference_point(0));
                }
                return value;
            }

            /**
             * @brief Node connexions with sub-elements
             */
            NodeConnectivity const static constexpr node_connectivity_ = {
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
             * @brief Reference nodes for the reference element
             */
            std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                    -1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                    -1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
            };

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
        requires (t_element.isLinearTetrahedron())
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {

            /**
             * @brief
             */
            template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<
                    std::tuple<
                            std::array<t_T<lolita::core2::geometry::Element::LinearTriangle(), t_domain, t_args...>, 4>
                    >,
                    std::tuple<
                            std::array<t_T<lolita::core2::geometry::Element::LinearSegment(), t_domain, t_args...>, 6>
                    >,
                    std::tuple<
                            std::array<t_T<lolita::core2::geometry::Element::Node(), t_domain, t_args...>, 4>
                    >
            >;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:

            /**
             * @brief
             * @param nodal_field_values
             * @param reference_point
             * @return
             */
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point
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
                    lolita::matrix::Vector<lolita::real, t_element.numNodes()> const &nodal_field_values,
                    lolita::domain::Point const &reference_point,
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

            /**
             * @brief
             */
            NodeConnectivity const static constexpr node_connectivity_ = {
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
             */
            std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                    +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                    +0.0000000000000000, +0.0000000000000000, +1.0000000000000000,
            };

        };

    }

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain>
    struct ElementTraits : detail::ElementOuterGeometryTraits<t_element, t_domain>, detail::ElementInnerGeometryTraits<t_element, t_domain>
    {

    private:

        using ElementOuterNeighbourhoodT = typename ElementTraits<t_element, t_domain>::template OuterConnectivity<detail::ElementView>;

        using ElementInnerNeighbourhoodT = typename ElementTraits<t_element, t_domain>::template InnerConnectivity<detail::ElementView>;

        using ElementsT = lolita::core2::geometry::Elements<detail::ElementView, t_domain>;

    public:

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core2::geometry::Element
        getElement()
        {
            return t_element;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::boolean
        hasDim(
                lolita::integer i
        )
        {
            return t_domain.dim_ - t_element.dim_ == i;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core2::geometry::ElementCoordinates
        getCoordinates()
        {
//            auto coordinates = lolita::core2::element::ElementCoordinates{-1, -1};
            auto coordinates = lolita::core2::geometry::ElementCoordinates{0, 0};
            auto set_coordinates = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & t_set_coordinates) constexpr mutable {
                using ElementT = typename std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementsT>>;
                if (ElementT::getElement() == t_element) {
                    coordinates.dim_ = t_i;
                    coordinates.tag_ = t_j;
                }
                if constexpr (t_j < std::tuple_size_v<std::tuple_element_t<t_i, ElementsT>> - 1) {
                    t_set_coordinates.template operator()<t_i, t_j + 1>(t_set_coordinates);
                }
                else if constexpr (t_i < std::tuple_size_v<ElementsT> - 1) {
                    t_set_coordinates.template operator()<t_i + 1, 0>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }

        /**
         * @brief
         * @tparam t_neighbour
         * @return
         */
        template<lolita::core2::geometry::Element t_neighbour>
        static constexpr
        lolita::core2::geometry::ElementCoordinates
        getComponentCoordinates()
        requires(!t_element.isPoint())
        {
//            auto coordinates = lolita::core2::element::ElementCoordinates{-1, -1};
            auto coordinates = lolita::core2::geometry::ElementCoordinates{0, 0};
            auto set_coordinates = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & t_set_coordinates) constexpr mutable {
                using NeighbourT = typename std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementInnerNeighbourhoodT>>::value_type;
                if (NeighbourT::getElement() == t_neighbour) {
                    coordinates.dim_ = t_i;
                    coordinates.tag_ = t_j;
                }
                if constexpr (t_j < std::tuple_size_v<std::tuple_element_t<t_i, ElementInnerNeighbourhoodT>> - 1) {
                    t_set_coordinates.template operator()<t_i, t_j + 1>(t_set_coordinates);
                }
                else if constexpr (t_i < std::tuple_size_v<ElementInnerNeighbourhoodT> - 1) {
                    t_set_coordinates.template operator()<t_i + 1, 0>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::core2::geometry::Element
        getComponent()
        requires(!t_element.isPoint())
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementInnerNeighbourhoodT>>::value_type::getElement();
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::core2::geometry::ElementTraits<getComponent<t_i, t_j>(), t_domain>
        getComponentDescription()
        requires(!t_element.isPoint())
        {
            return lolita::core2::geometry::ElementTraits<getComponent<t_i, t_j>(), t_domain>();
        }

        /**
         * @brief
         * @tparam t_i
         * @return
         */
        template<lolita::integer... t_i>
        static constexpr
        lolita::index
        getNumComponents()
        requires(!t_element.isPoint() && 0 <= sizeof...(t_i) <= 2)
        {
            if constexpr (sizeof...(t_i) == 2) {
                auto const constexpr t_coordinates = std::array<lolita::index, 2>{t_i...};
                return std::tuple_size_v<std::tuple_element_t<t_coordinates[1], std::tuple_element_t<t_coordinates[0], ElementInnerNeighbourhoodT>>>;
            }
            else if constexpr (sizeof...(t_i) == 1) {
                auto const constexpr t_coordinates = std::array<lolita::index, 1>{t_i...};
                return std::tuple_size_v<std::tuple_element_t<t_coordinates[0], ElementInnerNeighbourhoodT>>;
            }
            else {
                return std::tuple_size_v<ElementInnerNeighbourhoodT>;
            }
        }

        /**
         * @brief
         * @tparam _neighbour
         * @return
         */
        template<lolita::core2::geometry::Element t_neighbour>
        static constexpr
        lolita::core2::geometry::ElementCoordinates
        getNeighbourCoordinates()
        {
//            auto coordinates = lolita::core2::element::ElementCoordinates{-1, -1};
            auto coordinates = lolita::core2::geometry::ElementCoordinates{0, 0};
            auto set_coordinates = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & t_set_coordinates) constexpr mutable {
                using NeighbourT = typename std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementOuterNeighbourhoodT>>::value_type;
                if (NeighbourT::getElement() == t_neighbour) {
                    coordinates.dim_ = t_i;
                    coordinates.tag_ = t_j;
                }
                if constexpr (t_j < std::tuple_size_v<std::tuple_element_t<t_i, ElementOuterNeighbourhoodT>> - 1) {
                    t_set_coordinates.template operator()<t_i, t_j + 1>(t_set_coordinates);
                }
                else if constexpr (t_i < std::tuple_size_v<ElementOuterNeighbourhoodT> - 1) {
                    t_set_coordinates.template operator()<t_i + 1, 0>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::core2::geometry::Element
        getNeighbour()
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementOuterNeighbourhoodT>>::value_type::getElement();
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::core2::geometry::ElementTraits<getNeighbour<t_i, t_j>(), t_domain>
        getNeighbourDescription()
        {
            return lolita::core2::geometry::ElementTraits<getNeighbour<t_i, t_j>(), t_domain>();
        }

        /**
         * @brief
         * @tparam t_i
         * @return
         */
        template<lolita::integer... t_i>
        static constexpr
        lolita::index
        getNumNeighbours()
        requires(0 <= sizeof...(t_i) <= 1)
        {
            if constexpr (sizeof...(t_i) == 1) {
                auto const constexpr _coordinates = std::array<lolita::index, 1>{t_i...};
                return std::tuple_size_v<std::tuple_element_t<_coordinates[0], ElementOuterNeighbourhoodT>>;
            }
            else {
                return std::tuple_size_v<ElementOuterNeighbourhoodT>;
            }
        }

    };

    /**
     * @brief
     * @tparam t_domain
     */
    template<lolita::domain::Domain t_domain>
    struct DomainTraits
    {

    private:

        using ElementsT = lolita::core2::geometry::Elements<detail::ElementView, t_domain>;

    public:

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::domain::Domain
        getDomain()
        {
            return t_domain;
        }

        /**
         * @brief
         * @tparam t_element
         * @return
         */
        template<lolita::core2::geometry::Element t_element>
        static constexpr
        lolita::core2::geometry::ElementCoordinates
        getElementCoordinates()
        {
//            auto coordinates = lolita::core2::element::ElementCoordinates{t_element.dim(), -1};
            auto coordinates = lolita::core2::geometry::ElementCoordinates{t_element.dim(), 0};
            auto set_coordinates = [&] <lolita::integer t_i = 0> (auto & t_set_coordinates) constexpr mutable {
                using ElementT = std::tuple_element_t<t_i, std::tuple_element_t<t_element.dim(), ElementsT>>;
                if (ElementT::getElement() == t_element) {
                    coordinates.tag_ = t_i;
                }
                if constexpr (t_i < std::tuple_size_v<std::tuple_element_t<t_element.dim(), ElementsT>> - 1) {
                    t_set_coordinates.template operator()<t_i + 1>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::core2::geometry::Element
        getElement()
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementsT>>::getElement();
        }

//        /**
//         * @brief
//         * @tparam t_i
//         * @tparam t_j
//         * @return
//         */
//        template<lolita::integer t_i, lolita::integer t_j>
//        static constexpr
//        lolita::core2::element::ElementGeometryTraits<getElement<t_i, t_j>(), t_domain>
//        getElementDescription()
//        {
//            return lolita::core2::element::ElementGeometryTraits<getElement<t_i, t_j>(), t_domain>();
//        }

        /**
         * @brief
         * @tparam t_i
         * @return
         */
        template<lolita::index... t_i>
        static constexpr
        lolita::integer
        getNumElements()
        requires(0 <= sizeof...(t_i) <= 1)
        {
            if constexpr (sizeof...(t_i) == 1) {
                auto const constexpr _coordinates = std::array<lolita::index, 1>{t_i...};
                return std::tuple_size_v<std::tuple_element_t<_coordinates[0], ElementsT>>;
            }
            else {
                return std::tuple_size_v<ElementsT>;
            }
        }

    };

    /**
     * @brief
     * @tparam t_T
     * @tparam t_domain
     * @tparam t_arg
     */
    template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto> typename t_T, lolita::domain::Domain t_domain, auto t_args>
    struct ElementCollection
    {

    private:

        using ElementsT = lolita::core2::geometry::Elements<t_T, t_domain, t_args>;

    public:

        /**
         * @brief
         * @tparam t_i
         * @tparam _j
         * @return
         */
        template<lolita::index t_i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<t_i, ElementsT>> const &
        getElements()
        const
        {
            return std::get<_j>(std::get<t_i>(elements_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam _j
         * @return
         */
        template<lolita::index t_i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<t_i, ElementsT>> &
        getElements()
        {
            return std::get<_j>(std::get<t_i>(elements_));
        }

        /**
         * @brief The list of all available elements in a given domain
         */
        ElementsT elements_;

    };

}

#endif //LOLITA_LOLITA_CORE_2_HXX
