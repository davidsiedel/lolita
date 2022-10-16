#ifndef B940B5CA_76B5_4BC8_AB1B_916938C08B2D
#define B940B5CA_76B5_4BC8_AB1B_916938C08B2D

#include "2/core/traits/_include.hxx"

namespace lolita::core
{

    template<LagrangeShapeConcept auto t_element, auto...>
    struct ElementView
    {
        
        static constexpr
        LagrangeShapeConcept auto const &
        getElement()
        {
            return t_element;
        }

    };
    
    template<LagrangeShapeConcept auto t_element>
    struct ElementTraits;
    
    template<LagrangeShapeConcept auto t_element, auto...>
    using ElementNodeConnectivity = std::array<Natural, t_element.getNumNodes()>;
    
    template<LagrangeShapeConcept auto t_element>
    struct ElementInnerGeometryTraits;
    
    template<>
    struct ElementInnerGeometryTraits<Node{}>
    {
        
        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using InnerConnectivity = std::tuple<>;

    private:

        using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

    public:
    
        static
        Real
        getShapeMappingEvaluation(
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point
        )
        {
            return nodal_field_values(0);
        }
        
        static
        Real
        getShapeMappingDerivative(
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point,
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
    
    template<>
    struct ElementInnerGeometryTraits<LinearSegment{}>
    {
        
        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using InnerConnectivity = std::tuple<
            std::tuple<
                std::array<T<Node{}, t_args...>, 2>
            >
        >;

    private:

        using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

    public:
    
        static
        Real
        getShapeMappingEvaluation(
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point
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
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point,
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
    
    template<>
    struct ElementInnerGeometryTraits<LinearTriangle{}>
    {
        
        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using InnerConnectivity = std::tuple<
            std::tuple<
                std::array<T<LinearSegment{}, t_args...>, 3>
            >,
            std::tuple<
                std::array<T<Node{}, t_args...>, 3>
            >
        >;

    private:

        using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

    public:
    
        static
        Real
        getShapeMappingEvaluation(
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point
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
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point,
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
    
    template<>
    struct ElementInnerGeometryTraits<LinearQuadrangle{}>
    {
        
        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using InnerConnectivity = std::tuple<
            std::tuple<
                std::array<T<LinearSegment{}, t_args...>, 4>
            >,
            std::tuple<
                std::array<T<Node{}, t_args...>, 4>
            >
        >;

    private:

        using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

    public:
    
        static
        Real
        getShapeMappingEvaluation(
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point
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
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point,
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
    
    template<>
    struct ElementInnerGeometryTraits<LinearTetrahedron{}>
    {

        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using InnerConnectivity = std::tuple<
            std::tuple<
                std::array<T<LinearTriangle{}, t_args...>, 4>
            >,
            std::tuple<
                std::array<T<LinearSegment{}, t_args...>, 6>
            >,
            std::tuple<
                std::array<T<Node{}, t_args...>, 4>
            >
        >;

    private:

        using NodeConnectivity = InnerConnectivity<ElementNodeConnectivity>;

    public:
    
        static
        Real
        getShapeMappingEvaluation(
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point
        )
        {
            auto value = Real(0);
            return value;
        }
        
        static
        Real
        getShapeMappingDerivative(
            DenseVectorConcept<Real> auto const & nodal_field_values,
            PointConcept auto const & reference_point,
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
        
    template<LagrangeShapeConcept auto t_element>
    struct ElementOuterGeometryTraits
    {

    private:

        template<template<LagrangeShapeConcept auto, MeshConcept auto> typename T>
        struct OuterConnectivityTraits
        {

        private:
        
            template<LagrangeShapeConcept auto t__element, MeshConcept auto t__domain>
            using NeighbourVector = std::vector<T<t__element, t__domain>>;

            template<MeshConcept auto t__domain>
            using Expansion = typename ShapeLibrary::template Elements<NeighbourVector, t__domain>;

        public:

            template<MeshConcept auto t_domain>
            using OuterConnectivity = lolita::utility::tuple_slice_t<Expansion<t_domain>, t_element.getDim(), t_domain.getDim() + 1>;

        };

    public:
    
        template<template<LagrangeShapeConcept auto, MeshConcept auto> typename T, MeshConcept auto t__domain>
        using OuterConnectivity = typename OuterConnectivityTraits<T>::template OuterConnectivity<t__domain>;

    };
    
    template<>
    struct ElementOuterGeometryTraits<Node{}>
    {

    private:

        template<template<LagrangeShapeConcept auto, MeshConcept auto> typename T>
        struct OuterConnectivityTraits
        {

        private:

            template<LagrangeShapeConcept auto t__element, MeshConcept auto t__domain>
            using NeighbourVector = std::vector<T<t__element, t__domain>>;

            template<MeshConcept auto t__domain>
            using Expansion = typename ShapeLibrary::template Elements<NeighbourVector, t__domain>;

        public:

            template<MeshConcept auto t_domain>
            using OuterConnectivity = lolita::utility::tuple_slice_t<Expansion<t_domain>, 1, t_domain.getDim() + 1>;

        };

    public:
    
        template<template<LagrangeShapeConcept auto, MeshConcept auto> typename T, MeshConcept auto t__domain>
        using OuterConnectivity = typename OuterConnectivityTraits<T>::template OuterConnectivity<t__domain>;

    };

    template<LagrangeShapeConcept auto t_element>
    struct ElementTraits : ElementOuterGeometryTraits<t_element>, ElementInnerGeometryTraits<t_element>
    {

    private:

        using t_ElementInnerNeighborhood = typename ElementTraits<t_element>::template InnerConnectivity<ElementView>;

        template<MeshConcept auto t_domain>
        using t_ElementOuterNeighborhood = typename ElementTraits<t_element>::template OuterConnectivity<ElementView, t_domain>;
        
        template<MeshConcept auto t_domain>
        using t_Elements = lolita::utility::tuple_slice_t<typename ShapeLibrary::template Elements<ElementView>, 0, t_domain.getDim() + 1>;

    public:

        template<template<DomainConcept auto, auto...> typename T, auto... args_>
        using DomainConnectivity = std::tuple_element_t<t_element.getDim(), typename DomainLibrary::template Domains<T, args_...>>;
    
        static constexpr
        LagrangeShapeConcept auto const &
        getElement()
        {
            return t_element;
        }
        
        template<MeshConcept auto t_domain>
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
                using ElementT = typename std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements<t_domain>>>;
                if (ElementT::getElement() == t_element)
                {
                    coordinates.dim_ = t_i;
                    coordinates.tag_ = t_j;
                }
                if constexpr (t_j < std::tuple_size_v<std::tuple_element_t<t_i, t_Elements<t_domain>>> - 1)
                {
                    t_set_coordinates.template operator()<t_i, t_j + 1>(t_set_coordinates);
                }
                else if constexpr (t_i < std::tuple_size_v<t_Elements<t_domain>> - 1)
                {
                    t_set_coordinates.template operator()<t_i + 1, 0>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }
        
        template<LagrangeShapeConcept auto t_neighbour>
        static constexpr
        ElementCoordinates
        getInnerNeighborCoordinates()
        requires(t_element != Node())
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
        LagrangeShapeConcept auto
        getInnerNeighbor()
        requires(t_element != Node())
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_ElementInnerNeighborhood>>::value_type::getElement();
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(t_element != Node() && sizeof...(t_i) == 2)
        {
            auto constexpr t_coordinates = std::array<Integer, 2>{t_i...};
            return std::tuple_size_v<std::tuple_element_t<t_coordinates[1], std::tuple_element_t<t_coordinates[0], t_ElementInnerNeighborhood>>>;
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(t_element != Node() && sizeof...(t_i) == 1)
        {
            auto constexpr t_coordinates = std::array<Integer, 1>{t_i...};
            return std::tuple_size_v<std::tuple_element_t<t_coordinates[0], t_ElementInnerNeighborhood>>;
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(t_element != Node() && sizeof...(t_i) == 0)
        {
            return std::tuple_size_v<t_ElementInnerNeighborhood>;
        }
        
        template<MeshConcept auto t_domain, LagrangeShapeConcept auto t_neighbour>
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
                using NeighbourT = typename std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_ElementOuterNeighborhood<t_domain>>>::value_type;
                if (NeighbourT::getElement() == t_neighbour)
                {
                    coordinates.dim_ = t_i;
                    coordinates.tag_ = t_j;
                }
                if constexpr (t_j < std::tuple_size_v<std::tuple_element_t<t_i, t_ElementOuterNeighborhood<t_domain>>> - 1)
                {
                    t_set_coordinates.template operator()<t_i, t_j + 1>(t_set_coordinates);
                }
                else if constexpr (t_i < std::tuple_size_v<t_ElementOuterNeighborhood<t_domain>> - 1)
                {
                    t_set_coordinates.template operator()<t_i + 1, 0>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }
        
        template<MeshConcept auto t_domain, Integer t_i, Integer t_j>
        static constexpr
        LagrangeShapeConcept auto
        getOuterNeighbor()
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_ElementOuterNeighborhood<t_domain>>>::value_type::getElement();
        }
        
        template<MeshConcept auto t_domain, Integer... t_i>
        static constexpr
        Integer
        getNumOuterNeighbors()
        requires(sizeof...(t_i) == 1)
        {
            auto const constexpr _coordinates = std::array<Integer, 1>{t_i...};
            return std::tuple_size_v<std::tuple_element_t<_coordinates[0], t_ElementOuterNeighborhood<t_domain>>>;
        }
        
        template<MeshConcept auto t_domain, Integer... t_i>
        static constexpr
        Integer
        getNumOuterNeighbors()
        requires(sizeof...(t_i) == 0)
        {
            return std::tuple_size_v<t_ElementOuterNeighborhood<t_domain>>;
        }

    };
    
} // namespace lolita::core


#endif /* B940B5CA_76B5_4BC8_AB1B_916938C08B2D */
