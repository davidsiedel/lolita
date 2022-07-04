#ifndef D95DBEAA_910F_4201_AFBF_FEDDCE22F29E
#define D95DBEAA_910F_4201_AFBF_FEDDCE22F29E

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"

namespace lolita2::geometry
{
    
    struct ElementCoordinates
    {
        
        lolita::index dim_;
        
        lolita::index tag_;
        
        constexpr
        lolita::index
        dim()
        const
        {
            return dim_;
        }
        
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

    struct Element
    {

        enum Shape
        {
            
            Point,
            Segment,
            Triangle,
            Quadrangle,
            Polygon,
            Tetrahedron,
            Polyhedron,

        };

        Shape shape_;

        lolita::integer dim_;

        lolita::integer ord_;

        lolita::integer num_nodes_;

        lolita::integer num_edges_;

        lolita::integer num_faces_;
        
        constexpr
        Element(
                Shape shape,
                lolita::integer dim,
                lolita::integer ord,
                lolita::integer num_nodes,
                lolita::integer num_edges,
                lolita::integer num_faces
        )
        :
        shape_(shape),
        dim_(dim),
        ord_(ord),
        num_nodes_(num_nodes),
        num_edges_(num_edges),
        num_faces_(num_faces)
        {}

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
        
        static constexpr
        Element
        node()
        {
            return Element(Shape::Point, 0, 0, 1, 0, 0);
        }
        
        constexpr
        lolita::boolean
        isNode()
        const
        {
            return * this == node();
        }

        static constexpr
        Element
        segment(
            lolita::integer ord
        )
        {
            return Element(Shape::Segment, 1, ord, ord + 1, 0, 0);
        }
        
        constexpr
        lolita::boolean
        isSegment(
                lolita::integer ord
        )
        const
        {
            return * this == segment(ord);
        }

        static constexpr
        Element
        triangle(
            lolita::integer ord
        )
        {
            return Element(Shape::Triangle, 2, ord, 3 * ord, 3, 0);
        }
        
        constexpr
        lolita::boolean
        isTriangle(
                lolita::integer ord
        )
        const
        {
            return * this == triangle(ord);
        }

        static constexpr
        Element
        quadrangle(
            lolita::integer ord
        )
        {
            return Element(Shape::Quadrangle, 2, ord, 4 * ord, 4, 0);
        }
        
        constexpr
        lolita::boolean
        isQuadrangle(
                lolita::integer ord
        )
        const
        {
            return * this == quadrangle(ord);
        }

        static constexpr
        Element
        polygon(
            lolita::integer ord,
            lolita::integer num_edges
        )
        {
            return Element(Shape::Polygon, 1, ord, num_edges * ord, num_edges, 0);
        }

        static constexpr
        Element
        tetrahedron(
            lolita::integer ord
        )
        {
            return Element(Shape::Tetrahedron, 3, ord, 4, 6, 4);
        }

        constexpr
        lolita::boolean
        isTetrahedron(
            lolita::integer ord
        )
        const
        {
            return * this == tetrahedron(ord);
        }

        static constexpr
        Element
        solid(
            lolita::integer ord,
            lolita::integer num_edges,
            lolita::integer num_faces,
            lolita::integer num_nodes
        )
        {
            return Element(Shape::Segment, 1, ord, num_edges * (ord - 1), num_edges, 0);
        }
        
        constexpr
        lolita::boolean
        hasShape(
                Shape shape
        )
        const
        {
            return shape_ == shape;
        }
        
        constexpr
        lolita::boolean
        hasOrd(
                lolita::integer ord
        )
        const
        {
            return ord_ == ord;
        }
        
        constexpr
        lolita::boolean
        isSub(
                Domain domain,
                lolita::integer level
        )
        const
        {
            return domain.dim_ - dim_ == level;
        }

        friend
        std::ostream &
        operator<<(
                std::ostream & os,
                Element const & element
        )
        {
            if (element.shape_ == Shape::Segment)
            {
                os << "Segment " << element.ord_;
            }
            else if (element.shape_ == Shape::Triangle)
            {
                os << "Triangle " << element.ord_;
            }
            else if (element.shape_ == Shape::Quadrangle)
            {
                os << "Quadrangle " << element.ord_;
            }
            else if (element.shape_ == Shape::Point)
            {
                os << "Node " << element.ord_;
            }
            return os;
        }

    };
    
    template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    using Points = std::tuple<
            t_T<Element::node(), t_domain, t_args...>
    >;
    
    template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    using Curves = std::tuple<
            t_T<Element::segment(1), t_domain, t_args...>
    >;
    
    template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    using Facets = std::tuple<
            t_T<Element::triangle(1), t_domain, t_args...>,
            t_T<Element::quadrangle(1), t_domain, t_args...>
    >;
    
    template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    using Solids = std::tuple<
            t_T<Element::tetrahedron(1), t_domain, t_args...>
    >;

    namespace detail
    {

        template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
        using Elements = std::tuple<
                Points<t_T, t_domain, t_args...>,
                Curves<t_T, t_domain, t_args...>,
                Facets<t_T, t_domain, t_args...>,
                Solids<t_T, t_domain, t_args...>
        >;

    }

    template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    using Elements = lolita::utility::tuple_slice_t<detail::Elements<t_T, t_domain, t_args...>, 0, t_domain.dim_ + 1>;
    
    template<Element t_element, Domain t_domain>
    struct ElementTraits;

    namespace detail
    {
        
        template<Element t_element, Domain t_domain>
        struct ElementView
        {

            Element static constexpr element_ = t_element;
            
            static constexpr
            Element
            getElement()
            {
                return t_element;
            }

        };

    }

    namespace element_geometry
    {
        
        template<Element t_element, Domain t_domain>
        struct ElementOuterGeometryTraits;
        
        template<Element t_element, Domain t_domain>
        requires (t_element.isNode())
        struct ElementOuterGeometryTraits<t_element, t_domain>
        {

        private:

            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            struct OuterConnectivityTraits {

            private:

                template<Element t__element, Domain t__domain, auto... t__args>
                using NeighbourVector = std::vector<t_T<t__element, t__domain, t__args...>>;

            public:

                using OuterConnectivity = lolita::utility::tuple_slice_t<Elements<NeighbourVector, t_domain, t_args...>, 1, t_domain.dim_ + 1>;

            };

        public:
        
            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            using OuterConnectivity = typename OuterConnectivityTraits<t_T, t_args...>::OuterConnectivity;

        };
        
        template<Element t_element, Domain t_domain>
        requires (!t_element.isNode())
        struct ElementOuterGeometryTraits<t_element, t_domain>
        {

        private:

            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            struct OuterConnectivityTraits {

            private:

                template<Element t__element, Domain t__domain, auto... t__args>
                using NeighbourVector = std::vector<t_T<t__element, t__domain, t__args...>>;

            public:

                using OuterConnectivity = lolita::utility::tuple_slice_t<Elements<NeighbourVector, t_domain, t_args...>, t_element.dim_, t_domain.dim_ + 1>;

            };

        public:
        
            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            using OuterConnectivity = typename OuterConnectivityTraits<t_T, t_args...>::OuterConnectivity;

        };
        
        template<Element t_element, Domain t_domain>
        struct ElementInnerGeometryTraits;
        
        template<Element t_element, auto...>
        using ElementNodeConnectivity = std::array<lolita::index, t_element.num_nodes_>;
        
        template<Element t_element, Domain t_domain>
        requires (t_element.isNode())
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {
            
            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<>;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:
        
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point
            )
            {
                return nodal_field_values(0);
            }
            
            static
            lolita::real
            getShapeMappingDerivative(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point,
                    lolita::index derivative_direction
            )
            {
                return lolita::real(0);
            }
            
            NodeConnectivity const static constexpr node_connectivity_ = NodeConnectivity{};
            
            std::array<std::array<lolita::real, 3>, 1> const static constexpr reference_nodes_ = {
                    +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
            };

        };
        
        template<Element t_element, Domain t_domain>
        requires (t_element.isSegment(1))
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {
            
            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<
                    std::tuple<
                            std::array<t_T<Element::node(), t_domain, t_args...>, 2>
                    >
            >;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:
        
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point
            )
            {
                auto value = lolita::real(0);
                value += nodal_field_values(0) * (1.0 / 2.0) * (1.0 - reference_point(0));
                value += nodal_field_values(1) * (1.0 / 2.0) * (1.0 + reference_point(0));
                return value;
            }
            
            static
            lolita::real
            getShapeMappingDerivative(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point,
                    lolita::index derivative_direction
            )
            {
                assert(derivative_direction == 0);
                auto value = lolita::real(0);
                value += -nodal_field_values(0) * (1.0 / 2.0);
                value += +nodal_field_values(1) * (1.0 / 2.0);
                return value;
            }
            
            NodeConnectivity const static constexpr node_connectivity_ = NodeConnectivity{
                    {
                            {
                                    0,
                                    1,
                            }
                    }
            };
            
            std::array<std::array<lolita::real, 3>, 2> const static constexpr reference_nodes_ = {
                    -1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, +0.0000000000000000, +0.0000000000000000
            };

        };
        
        template<Element t_element, Domain t_domain>
        requires (t_element.isTriangle(1))
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {
            
            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<
                    std::tuple<
                            std::array<t_T<Element::segment(1), t_domain, t_args...>, 3>
                    >,
                    std::tuple<
                            std::array<t_T<Element::node(), t_domain, t_args...>, 3>
                    >
            >;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:
        
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point
            )
            {
                auto value = lolita::real(0);
                value += nodal_field_values(0) * (1.0 - reference_point(0) - reference_point(1));
                value += nodal_field_values(1) * reference_point(0);
                value += nodal_field_values(2) * reference_point(1);
                return value;
            }
            
            static
            lolita::real
            getShapeMappingDerivative(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point,
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
            
            std::array<std::array<lolita::real, 3>, 3> const static constexpr reference_nodes_ = {
                    +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
            };

        };
        
        template<Element t_element, Domain t_domain>
        requires (t_element.isQuadrangle(1))
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {
            
            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<
                    std::tuple<
                            std::array<t_T<Element::segment(1), t_domain, t_args...>, 4>
                    >,
                    std::tuple<
                            std::array<t_T<Element::node(), t_domain, t_args...>, 4>
                    >
            >;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:
        
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point
            )
            {
                auto value = lolita::real(0);
                value += nodal_field_values(0) * (1.0 / 4.0) * (1.0 - reference_point(0)) * (1.0 - reference_point(1));
                value += nodal_field_values(1) * (1.0 / 4.0) * (1.0 + reference_point(0)) * (1.0 - reference_point(1));
                value += nodal_field_values(2) * (1.0 / 4.0) * (1.0 + reference_point(0)) * (1.0 + reference_point(1));
                value += nodal_field_values(3) * (1.0 / 4.0) * (1.0 - reference_point(0)) * (1.0 + reference_point(1));
                return value;
            }
            
            static
            lolita::real
            getShapeMappingDerivative(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point,
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
            
            std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                    -1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, -1.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                    -1.0000000000000000, +1.0000000000000000, +0.0000000000000000,
            };

        };
        
        template<Element t_element, Domain t_domain>
        requires (t_element.isTetrahedron(1))
        struct ElementInnerGeometryTraits<t_element, t_domain>
        {

            template<template<Element, Domain, auto...> typename t_T, auto... t_args>
            using InnerConnectivity = std::tuple<
                    std::tuple<
                            std::array<t_T<Element::triangle(1), t_domain, t_args...>, 4>
                    >,
                    std::tuple<
                            std::array<t_T<Element::segment(1), t_domain, t_args...>, 6>
                    >,
                    std::tuple<
                            std::array<t_T<Element::node(), t_domain, t_args...>, 4>
                    >
            >;

        private:

            using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

        public:
        
            static
            lolita::real
            getShapeMappingEvaluation(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point
            )
            {
                auto value = lolita::real(0);
                return value;
            }
            
            static
            lolita::real
            getShapeMappingDerivative(
                    lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const &nodal_field_values,
                    lolita2::Point const &reference_point,
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
            
            std::array<std::array<lolita::real, 3>, 4> const static constexpr reference_nodes_ = {
                    +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                    +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
                    +0.0000000000000000, +0.0000000000000000, +1.0000000000000000,
            };

        };

    }
    
    template<Element t_element, Domain t_domain>
    struct ElementTraits
    :
    element_geometry::ElementOuterGeometryTraits<t_element, t_domain>,
    element_geometry::ElementInnerGeometryTraits<t_element, t_domain>
    {

    private:

        using ElementOuterNeighbourhoodT = typename ElementTraits<t_element, t_domain>::template OuterConnectivity<detail::ElementView>;

        using ElementInnerNeighbourhoodT = typename ElementTraits<t_element, t_domain>::template InnerConnectivity<detail::ElementView>;

        using ElementsT = lolita2::geometry::Elements<detail::ElementView, t_domain>;

    public:
    
        static constexpr
        Element
        getElement()
        {
            return t_element;
        }
        
        static constexpr
        lolita::boolean
        hasDim(
                lolita::integer i
        )
        {
            return t_domain.dim_ - t_element.dim_ == i;
        }
        
        static constexpr
        lolita2::geometry::ElementCoordinates
        getCoordinates()
        {
            auto coordinates = lolita2::geometry::ElementCoordinates{0, 0};
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
        
        template<Element t_neighbour>
        static constexpr
        lolita2::geometry::ElementCoordinates
        getComponentCoordinates()
        requires(!t_element.isNode())
        {
            auto coordinates = lolita2::geometry::ElementCoordinates{0, 0};
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
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        Element
        getComponent()
        requires(!t_element.isNode())
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementInnerNeighbourhoodT>>::value_type::getElement();
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita2::geometry::ElementTraits<getComponent<t_i, t_j>(), t_domain>
        getComponentDescription()
        requires(!t_element.isNode())
        {
            return lolita2::geometry::ElementTraits<getComponent<t_i, t_j>(), t_domain>();
        }
        
        template<lolita::integer... t_i>
        static constexpr
        lolita::index
        getNumComponents()
        requires(!t_element.isNode() && 0 <= sizeof...(t_i) <= 2)
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
        
        template<lolita2::geometry::Element t_neighbour>
        static constexpr
        lolita2::geometry::ElementCoordinates
        getNeighbourCoordinates()
        {
            auto coordinates = lolita2::geometry::ElementCoordinates{0, 0};
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
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita2::geometry::Element
        getNeighbour()
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementOuterNeighbourhoodT>>::value_type::getElement();
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita2::geometry::ElementTraits<getNeighbour<t_i, t_j>(), t_domain>
        getNeighbourDescription()
        {
            return lolita2::geometry::ElementTraits<getNeighbour<t_i, t_j>(), t_domain>();
        }
        
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
    
    template<Domain t_domain>
    struct DomainTraits
    {

    private:

        using ElementsT = lolita2::geometry::Elements<detail::ElementView, t_domain>;

    public:
    
        static constexpr
        Domain
        getDomain()
        {
            return t_domain;
        }
        
        template<lolita2::geometry::Element t_element>
        static constexpr
        lolita2::geometry::ElementCoordinates
        getElementCoordinates()
        {
            auto coordinates = lolita2::geometry::ElementCoordinates{t_element.dim_, 0};
            auto set_coordinates = [&] <lolita::integer t_i = 0> (auto & t_set_coordinates) constexpr mutable {
                using ElementT = std::tuple_element_t<t_i, std::tuple_element_t<t_element.dim_, ElementsT>>;
                if (ElementT::getElement() == t_element) {
                    coordinates.tag_ = t_i;
                }
                if constexpr (t_i < std::tuple_size_v<std::tuple_element_t<t_element.dim_, ElementsT>> - 1) {
                    t_set_coordinates.template operator()<t_i + 1>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita2::geometry::Element
        getElement()
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, ElementsT>>::getElement();
        }
        
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
    
    template<template<lolita2::geometry::Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    struct ElementCollection
    {

    private:

        using ElementsT = lolita2::geometry::Elements<t_T, t_domain, t_args...>;

    public:
    
        template<lolita::index t_i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<t_i, ElementsT>> const &
        getElements()
        const
        {
            return std::get<_j>(std::get<t_i>(elements_));
        }
        
        template<lolita::index t_i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<t_i, ElementsT>> &
        getElements()
        {
            return std::get<_j>(std::get<t_i>(elements_));
        }
        
        ElementsT elements_;

    };

}

#endif /* D95DBEAA_910F_4201_AFBF_FEDDCE22F29E */
