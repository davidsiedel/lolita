#ifndef D95DBEAA_910F_4201_AFBF_FEDDCE22F29E
#define D95DBEAA_910F_4201_AFBF_FEDDCE22F29E

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"

namespace lolita
{
    
    struct ElementCoordinates
    {

        constexpr
        ElementCoordinates()
        :
        dim_(-1),
        tag_(-1)
        {}
        
        constexpr
        ElementCoordinates(
            Integer dim,
            Integer tag
        )
        :
        dim_(dim),
        tag_(tag)
        {}
        
        constexpr
        Integer
        dim()
        const
        {
            return dim_;
        }
        
        constexpr
        Integer
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
        
        Integer dim_;
        
        Integer tag_;

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

        static constexpr
        std::basic_string_view<Character>
        getShapeLabel(
            Shape shape
        )
        noexcept
        {
            switch (shape)
            {
                case Shape::Point:          return "Point       ";
                case Shape::Segment:        return "Segment     ";
                case Shape::Triangle:       return "Triangle    ";
                case Shape::Quadrangle:     return "Quadrangle  ";
                case Shape::Polygon:        return "Polygon     ";
                case Shape::Tetrahedron:    return "Tetrahedron ";
                case Shape::Polyhedron:     return "Polyhedron  ";
                default : return "";
            }
        }
        
        constexpr
        Element(
                Shape shape,
                Integer dim,
                Integer ord,
                Integer num_nodes,
                Integer num_edges,
                Integer num_faces
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
        Boolean
        operator==(
                Element const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
                Element const & other
        )
        const = default;

        constexpr
        std::basic_string_view<Character>
        getShapeLabel()
        const noexcept
        {
            return getShapeLabel(shape_);
        }
        
        static constexpr
        Element
        node()
        {
            return Element(Shape::Point, 0, 0, 1, 0, 0);
        }
        
        constexpr
        Boolean
        isNode()
        const
        {
            return * this == node();
        }
        
        constexpr
        Boolean
        isCurve()
        const
        {
            return dim_ == 1;
        }
        
        constexpr
        Boolean
        isFacet()
        const
        {
            return dim_ == 2;
        }
        
        constexpr
        Boolean
        isSolid()
        const
        {
            return dim_ == 3;
        }

        static constexpr
        Element
        segment(
            Integer ord
        )
        {
            return Element(Shape::Segment, 1, ord, ord + 1, 0, 0);
        }
        
        constexpr
        Boolean
        isSegment(
            Integer ord
        )
        const
        {
            return * this == segment(ord);
        }

        static constexpr
        Element
        triangle(
            Integer ord
        )
        {
            return Element(Shape::Triangle, 2, ord, 3 * ord, 3, 0);
        }
        
        constexpr
        Boolean
        isTriangle(
            Integer ord
        )
        const
        {
            return * this == triangle(ord);
        }

        static constexpr
        Element
        quadrangle(
            Integer ord
        )
        {
            return Element(Shape::Quadrangle, 2, ord, 4 * ord, 4, 0);
        }
        
        constexpr
        Boolean
        isQuadrangle(
            Integer ord
        )
        const
        {
            return * this == quadrangle(ord);
        }

        static constexpr
        Element
        polygon(
            Integer ord,
            Integer num_edges
        )
        {
            return Element(Shape::Polygon, 1, ord, num_edges * ord, num_edges, 0);
        }

        static constexpr
        Element
        tetrahedron(
            Integer ord
        )
        {
            return Element(Shape::Tetrahedron, 3, ord, 4, 6, 4);
        }

        constexpr
        Boolean
        isTetrahedron(
            Integer ord
        )
        const
        {
            return * this == tetrahedron(ord);
        }

        static constexpr
        Element
        solid(
            Integer ord,
            Integer num_edges,
            Integer num_faces,
            Integer num_nodes
        )
        {
            return Element(Shape::Segment, 1, ord, num_edges * (ord - 1), num_edges, 0);
        }
        
        constexpr
        Boolean
        hasShape(
            Shape shape
        )
        const
        {
            return shape_ == shape;
        }
        
        constexpr
        Boolean
        hasOrd(
            Integer ord
        )
        const
        {
            return ord_ == ord;
        }
        
        constexpr
        Boolean
        isSub(
            Domain domain,
            Integer level
        )
        const
        {
            return domain.getDim() - dim_ == level;
        }

        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            Element const & element
        )
        {
            os << "{ " << element.getShapeLabel() << " | ";
            os << element.num_nodes_ << " nodes | ";
            os << element.num_edges_ << " edges | ";
            os << element.num_faces_ << " faces }";
            return os;
        }
        
        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }
        
        constexpr
        Integer
        getOrd()
        const
        {
            return ord_;
        }
        
        constexpr
        Integer
        getNumNodes()
        const
        {
            return num_nodes_;
        }
        
        constexpr
        Integer
        getNumEdges()
        const
        {
            return num_edges_;
        }
        
        constexpr
        Integer
        getNumFaces()
        const
        {
            return num_faces_;
        }

        Shape shape_;

        Integer dim_;

        Integer ord_;

        Integer num_nodes_;

        Integer num_edges_;

        Integer num_faces_;

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
            Element const &
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
        using ElementNodeConnectivity = std::array<Integer, t_element.getNumNodes()>;
        
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
            Real
            getShapeMappingEvaluation(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point
            )
            {
                return nodal_field_values(0);
            }
            
            static
            Real
            getShapeMappingDerivative(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point,
                Integer derivative_direction
            )
            {
                return Real(0);
            }
            
            NodeConnectivity const static constexpr node_connectivity_ = NodeConnectivity{};
            
            std::array<std::array<Real, 3>, 1> const static constexpr reference_nodes_ = {
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
            Real
            getShapeMappingEvaluation(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point
            )
            {
                auto value = Real(0);
                value += nodal_field_values(0) * (1.0 / 2.0) * (1.0 - reference_point(0));
                value += nodal_field_values(1) * (1.0 / 2.0) * (1.0 + reference_point(0));
                return value;
            }
            
            static
            Real
            getShapeMappingDerivative(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point,
                Integer derivative_direction
            )
            {
                assert(derivative_direction == 0);
                auto value = Real(0);
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
            
            std::array<std::array<Real, 3>, 2> const static constexpr reference_nodes_ = {
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
            Real
            getShapeMappingEvaluation(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point
            )
            {
                auto value = Real(0);
                value += nodal_field_values(0) * (1.0 - reference_point(0) - reference_point(1));
                value += nodal_field_values(1) * reference_point(0);
                value += nodal_field_values(2) * reference_point(1);
                return value;
            }
            
            static
            Real
            getShapeMappingDerivative(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point,
                Integer derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 1);
                auto value = Real(0);
                if (derivative_direction == 0)
                {
                    value += -nodal_field_values(0);
                    value += +nodal_field_values(1);
                    value += +0.0;
                }
                else
                {
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
            
            std::array<std::array<Real, 3>, 3> const static constexpr reference_nodes_ = {
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
            Real
            getShapeMappingEvaluation(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point
            )
            {
                auto value = Real(0);
                value += nodal_field_values(0) * (1.0 / 4.0) * (1.0 - reference_point(0)) * (1.0 - reference_point(1));
                value += nodal_field_values(1) * (1.0 / 4.0) * (1.0 + reference_point(0)) * (1.0 - reference_point(1));
                value += nodal_field_values(2) * (1.0 / 4.0) * (1.0 + reference_point(0)) * (1.0 + reference_point(1));
                value += nodal_field_values(3) * (1.0 / 4.0) * (1.0 - reference_point(0)) * (1.0 + reference_point(1));
                return value;
            }
            
            static
            Real
            getShapeMappingDerivative(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point,
                Integer derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 1);
                auto value = Real(0);
                if (derivative_direction == 0)
                {
                    value += -nodal_field_values(0) * (1.0 / 4.0) * (1.0 - reference_point(1));
                    value += +nodal_field_values(1) * (1.0 / 4.0) * (1.0 - reference_point(1));
                    value += +nodal_field_values(2) * (1.0 / 4.0) * (1.0 + reference_point(1));
                    value += -nodal_field_values(3) * (1.0 / 4.0) * (1.0 + reference_point(1));
                }
                else
                {
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
            
            std::array<std::array<Real, 3>, 4> const static constexpr reference_nodes_ = {
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
            Real
            getShapeMappingEvaluation(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point
            )
            {
                auto value = Real(0);
                return value;
            }
            
            static
            Real
            getShapeMappingDerivative(
                lolita::algebra::Vector<Real, t_element.getNumNodes()> const &nodal_field_values,
                Point const &reference_point,
                Integer derivative_direction
            )
            {
                assert(0 <= derivative_direction <= 2);
                auto value = Real(0);
                if (derivative_direction == 0)
                {

                }
                else if (derivative_direction == 1)
                {

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
            
            std::array<std::array<Real, 3>, 4> const static constexpr reference_nodes_ = {
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

        using t_ElementOuterNeighborhood = typename ElementTraits<t_element, t_domain>::template OuterConnectivity<detail::ElementView>;

        using t_ElementInnerNeighborhood = typename ElementTraits<t_element, t_domain>::template InnerConnectivity<detail::ElementView>;

        using t_Elements = Elements<detail::ElementView, t_domain>;

    public:
    
        static constexpr
        Element
        getElement()
        {
            return t_element;
        }
        
        static constexpr
        Boolean
        hasDim(
                Integer i
        )
        {
            return t_domain.dim_ - t_element.dim_ == i;
        }
        
        static constexpr
        ElementCoordinates
        getCoordinates()
        {
            auto coordinates = ElementCoordinates();
            auto set_coordinates = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_coordinates
            )
            constexpr mutable
            {
                using ElementT = typename std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>>;
                if (ElementT::getElement() == t_element)
                {
                    coordinates.dim_ = t_i;
                    coordinates.tag_ = t_j;
                }
                if constexpr (t_j < std::tuple_size_v<std::tuple_element_t<t_i, t_Elements>> - 1)
                {
                    t_set_coordinates.template operator()<t_i, t_j + 1>(t_set_coordinates);
                }
                else if constexpr (t_i < std::tuple_size_v<t_Elements> - 1)
                {
                    t_set_coordinates.template operator()<t_i + 1, 0>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }
        
        template<Element t_neighbour>
        static constexpr
        ElementCoordinates
        getInnerNeighborCoordinates()
        requires(!t_element.isNode())
        {
            auto coordinates = ElementCoordinates();
            auto set_coordinates = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_coordinates
            )
            constexpr mutable
            {
                using NeighbourT = typename std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_ElementInnerNeighborhood>>::value_type;
                if (NeighbourT::getElement() == t_neighbour)
                {
                    coordinates.dim_ = t_i;
                    coordinates.tag_ = t_j;
                }
                if constexpr (t_j < std::tuple_size_v<std::tuple_element_t<t_i, t_ElementInnerNeighborhood>> - 1)
                {
                    t_set_coordinates.template operator()<t_i, t_j + 1>(t_set_coordinates);
                }
                else if constexpr (t_i < std::tuple_size_v<t_ElementInnerNeighborhood> - 1)
                {
                    t_set_coordinates.template operator()<t_i + 1, 0>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Element
        getInnerNeighbor()
        requires(!t_element.isNode())
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_ElementInnerNeighborhood>>::value_type::getElement();
        }
        
        // template<Integer t_i, Integer t_j>
        // static constexpr
        // ElementTraits<getInnerNeighbor<t_i, t_j>(), t_domain>
        // getInnerNeighborDescription()
        // requires(!t_element.isNode())
        // {
        //     return ElementTraits<getInnerNeighbor<t_i, t_j>(), t_domain>();
        // }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(!t_element.isNode() && sizeof...(t_i) == 2)
        {
            auto const constexpr t_coordinates = std::array<Integer, 2>{t_i...};
            return std::tuple_size_v<std::tuple_element_t<t_coordinates[1], std::tuple_element_t<t_coordinates[0], t_ElementInnerNeighborhood>>>;
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(!t_element.isNode() && sizeof...(t_i) == 1)
        {
            auto const constexpr t_coordinates = std::array<Integer, 1>{t_i...};
            return std::tuple_size_v<std::tuple_element_t<t_coordinates[0], t_ElementInnerNeighborhood>>;
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(!t_element.isNode() && sizeof...(t_i) == 0)
        {
            return std::tuple_size_v<t_ElementInnerNeighborhood>;
        }
        
        template<Element t_neighbour>
        static constexpr
        ElementCoordinates
        getOuterNeighborCoordinates()
        {
            auto coordinates = ElementCoordinates();
            auto set_coordinates = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_coordinates
            )
            constexpr mutable
            {
                using NeighbourT = typename std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_ElementOuterNeighborhood>>::value_type;
                if (NeighbourT::getElement() == t_neighbour)
                {
                    coordinates.dim_ = t_i;
                    coordinates.tag_ = t_j;
                }
                if constexpr (t_j < std::tuple_size_v<std::tuple_element_t<t_i, t_ElementOuterNeighborhood>> - 1)
                {
                    t_set_coordinates.template operator()<t_i, t_j + 1>(t_set_coordinates);
                }
                else if constexpr (t_i < std::tuple_size_v<t_ElementOuterNeighborhood> - 1)
                {
                    t_set_coordinates.template operator()<t_i + 1, 0>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Element
        getOuterNeighbor()
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_ElementOuterNeighborhood>>::value_type::getElement();
        }
        
        // template<Integer t_i, Integer t_j>
        // static constexpr
        // ElementTraits<getOuterNeighbor<t_i, t_j>(), t_domain>
        // getOuterNeighborDescription()
        // {
        //     return ElementTraits<getOuterNeighbor<t_i, t_j>(), t_domain>();
        // }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumOuterNeighbors()
        requires(sizeof...(t_i) == 1)
        {
            auto const constexpr _coordinates = std::array<Integer, 1>{t_i...};
            return std::tuple_size_v<std::tuple_element_t<_coordinates[0], t_ElementOuterNeighborhood>>;
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumOuterNeighbors()
        requires(sizeof...(t_i) == 0)
        {
            return std::tuple_size_v<t_ElementOuterNeighborhood>;
        }

    };
    
    template<Domain t_domain>
    struct DomainTraits
    {

    private:

        using t_Elements = Elements<detail::ElementView, t_domain>;

    public:
    
        static constexpr
        Domain
        getDomain()
        {
            return t_domain;
        }
        
        template<Element t_element>
        static constexpr
        ElementCoordinates
        getElementCoordinates()
        {
            auto coordinates = ElementCoordinates(t_element.dim_, -1);
            auto set_coordinates = [&] <Integer t_i = 0> (
                auto & t_set_coordinates
            )
            constexpr mutable
            {
                using ElementT = std::tuple_element_t<t_i, std::tuple_element_t<t_element.dim_, t_Elements>>;
                if (ElementT::getElement() == t_element)
                {
                    coordinates.tag_ = t_i;
                }
                if constexpr (t_i < std::tuple_size_v<std::tuple_element_t<t_element.dim_, t_Elements>> - 1)
                {
                    t_set_coordinates.template operator()<t_i + 1>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Element
        getElement()
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>>::getElement();
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumElements()
        requires(0 <= sizeof...(t_i) <= 1)
        {
            if constexpr (sizeof...(t_i) == 1)
            {
                auto const constexpr _coordinates = std::array<Integer, 1>{t_i...};
                return std::tuple_size_v<std::tuple_element_t<_coordinates[0], t_Elements>>;
            }
            else
            {
                return std::tuple_size_v<t_Elements>;
            }
        }

    };

    struct MeshDomain
    {

        MeshDomain(
            Integer dim,
            std::basic_string<Character> const & tag
        )
        :
        dim_(dim),
        tag_(tag)
        {}

        Integer dim_;

        std::basic_string<Character> tag_;

    };
    
    template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    struct ElementMap
    {

    private:

        template<Element t_element, Domain t__domain, auto... t__args>
        using t_ElementMap = std::map<std::basic_string<Character>, std::shared_ptr<t_T<t_element, t__domain, t__args...>>>;

        using t_Elements = Elements<t_ElementMap, t_domain, t_args...>;

    public:
    
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> const &
        getElements()
        const
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> &
        getElements()
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        t_Elements elements_;

    };
    
    template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    struct ElementSet
    {

    private:

        template<Element t_element, Domain t__domain, auto... t__args>
        using t__Elements = std::vector<std::shared_ptr<t_T<t_element, t__domain, t__args...>>>;

        using t_Elements = Elements<t__Elements, t_domain, t_args...>;

    public:
    
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> const &
        getElements()
        const
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> &
        getElements()
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        t_Elements elements_;

    };

} // namespace lolita

#endif /* D95DBEAA_910F_4201_AFBF_FEDDCE22F29E */
