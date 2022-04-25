//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_CORE_ELEMENT_GEOMETRY_HXX
#define LOLITA_LOLITA_CORE_ELEMENT_GEOMETRY_HXX

#include "lolita/lolita_lolita.hxx"
#include "lolita/lolita_containers.hxx"
#include "lolita/lolita_core.hxx"
//#include "lolita/lolita_collection.hxx"

namespace lolita::core::element
{

    struct Element
    {

        constexpr
        Bool
        operator==(
                Element const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Element const &
                other
        )
        const = default;

        constexpr inline
        auto
        dimPoint()
        const
        {
            if (dim == 0) {
                return Indx(1);
            }
            else {
                return dim;
            }
        }

        Indx tag;

        Indx dim;

        Indx ord;

        Indx num_nodes;

    };

    template<Element E, Domain D, auto C>
    struct FiniteElement;

    auto const static constexpr pnt_00 = Element{0, 0, 0, 1};
    //
    auto const static constexpr seg_02 = Element{1, 1, 1, 2};
    auto const static constexpr seg_03 = Element{5, 1, 2, 3};
    //
    auto const static constexpr tri_03 = Element{2, 2, 1, 3};
    auto const static constexpr tri_06 = Element{6, 2, 2, 6};
    auto const static constexpr qua_04 = Element{3, 2, 1, 4};
    auto const static constexpr qua_08 = Element{7, 2, 2, 8};
    //
    auto const static constexpr tet_04 = Element{4, 3, 1, 4};
    auto const static constexpr tet_12 = Element{8, 3, 2, 12};
    //
    // auto const static constexpr pg5_5 = Element{5, 2, 1, 5};

    /*
     *
     */

    template<Element E, Domain D>
    concept Cell = D.dim - E.dim == 0;

    template<Element E, Domain D>
    concept Face = D.dim - E.dim == 1;

    template<Element E, Domain D>
    concept Edge = D.dim - E.dim == 2;

    template<Element E, Domain D>
    concept Node = D.dim - E.dim == 3;

    /*
     *
     */

    template<Element E>
    concept Point = E == pnt_00;

    template<Element E>
    concept Segment = E == seg_02 || E == seg_03;

    template<Element E>
    concept Triangle = E == tri_03 || E == tri_06;

    template<Element E>
    concept Quadrangle = E == qua_04 || E == qua_08;

    template<Element E>
    concept Tetrahedron = E == tet_04 || E == tet_12;

    /*
     *
     */

    template<Element E>
    concept Curve = Segment<E>;

    template<Element E>
    concept Surface = Triangle<E> || Quadrangle<E>;

    template<Element E>
    concept Volume = Tetrahedron<E>;

    /*
     *
     */

    template<template<Element> typename T>
    using Points = Collection<T<pnt_00>>;

    template<template<Element> typename T>
    using Curves = Collection<T<seg_02>>;

    template<template<Element> typename T>
    using Facets = Collection<T<tri_03>, T<qua_04>>;

    template<template<Element> typename T>
    using Solids = Collection<T<tet_04>>;

    namespace detail
    {

        template<Element E>
        using ElementNodeConnectivity = Array<Indx, E.num_nodes>;

//        template<Indx D, template<Element> typename T>
        template<Indx D, template<Element, auto...> typename T>
        struct ElementsPolicy;

//        template<template<Element> typename T>
        template<template<Element, auto...> typename T>
        struct ElementsPolicy<1, T>
        {

            using Type = Collection<Points<T>, Curves<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementsPolicy<2, T>
        {

            using Type = Collection<Points<T>, Curves<T>, Facets<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementsPolicy<3, T>
        {

            using Type = Collection<Points<T>, Curves<T>, Facets<T>, Solids<T>>;

        };

        template<Indx I, Indx D, template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy;

        /*
         * OD
         */

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<0, 3, T>
        {

            using Type = Collection<Curves<T>, Facets<T>, Solids<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<0, 2, T>
        {

            using Type = Collection<Curves<T>, Facets<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<0, 1, T>
        {

            using Type = Collection<Curves<T>>;

        };

        /*
         * 1D
         */

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<1, 3, T>
        {

            using Type = Collection<Curves<T>, Facets<T>, Solids<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<1, 2, T>
        {

            using Type = Collection<Curves<T>, Facets<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<1, 1, T>
        {

            using Type = Collection<Curves<T>>;

        };

        /*
         * 2D
         */

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<2, 3, T>
        {

            using Type = Collection<Facets<T>, Solids<T>>;

        };

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<2, 2, T>
        {

            using Type = Collection<Facets<T>>;

        };

        /*
         * 3D
         */

        template<template<Element, auto...> typename T>
        struct ElementNeighbourhoodPolicy<3, 3, T>
        {

            using Type = Collection<Solids<T>>;

        };

        template<Element E, Indx D, template<Element, auto...> typename T, auto... A>
        struct ElementNeighbourArrayPolicy
        {

        private:

            auto const static constexpr I = E.dim;

            template<Element EE>
            using NeighbourArray = Array<T<EE, A...>>;

        public:

            using Type = typename ElementNeighbourhoodPolicy<I, D, NeighbourArray>::Type;

        };

    }

    template<Element E>
    struct ElementGeometry;

    template<Element E>
    requires(E == pnt_00)
    struct ElementGeometry<E>
    {

        auto const static constexpr element = E;

        constexpr
        Bool
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

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

    template<Element E>
    requires(E == seg_02)
    struct ElementGeometry<E>
    {

        auto const static constexpr element = E;

        constexpr
        Bool
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        auto const static constexpr lowers = Array<Element, 0>{};
        auto const static constexpr higher = Array<Element, 1>{seg_03};

        /**
         * @brief
         */
        auto const static constexpr reference_nodes = Array<Real, 2, 1>{
            -1.0000000000000000,
            +1.0000000000000000
        };

        template<template<Element, auto...> typename T>
        using Components = Collection<
                Collection<
                        Array<T<pnt_00>, 2>
                >
        >;

//        template<template<Element, auto...> typename T>
//        using Components2 = Collection<
//                Collection<
//                        T<pnt_00>,
//                        T<pnt_00>
//                >
//        >;
//
//    private:
//
//        using Connectivity = Collection<
//                Collection<Array<Indx, 2, 1>>
//        >;
//
//    public:

//        auto const static constexpr node_connectivity = Connectivity(
//                Collection(
//                        Array<Indx, 2, 1>{
//                                0,
//                                1
//                        }
//                )
//        );

//        auto const static constexpr node_connectivity = Components<detail::ElementNodeConnectivity>(
//                Collection(
//                        Array<Indx, 2, 1>{
//                                0,
//                                1
//                        }
//                )
//        );

        auto const static constexpr node_connectivity = Components<detail::ElementNodeConnectivity>{
                {
                    {
                        0,
                        1,
                    }
                }
        };


//        auto const static constexpr node_con = Components2<detail::ElementNodeConnectivity>{
//                {
//                    {0},
//                    {1}
//                }
//        };

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

    template<Element E>
    requires(E == tri_03)
    struct ElementGeometry<E>
    {

        auto const static constexpr element = E;

        constexpr
        Bool
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        /**
         * @brief
         */
        auto const static constexpr reference_nodes = Array<Real, 3, 2>{
            +0.0000000000000000, +0.0000000000000000,
            +1.0000000000000000, +0.0000000000000000,
            +0.0000000000000000, +1.0000000000000000
        };

        template<template<Element, auto...> typename T>
        using Components = Collection<
                Collection<
                        Array<T<seg_02>, 3>
                >,
                Collection<
                        Array<T<pnt_00>, 3>
                >
        >;

//        template<template<Element, auto...> typename T>
//        using Components2 = Collection<
//                Collection<
//                        T<seg_02>,
//                        T<seg_02>,
//                        T<seg_02>
//                >,
//                Collection<
//                        T<pnt_00>,
//                        T<pnt_00>,
//                        T<pnt_00>
//                >
//        >;
//
//    private:
//
//        using Connectivity = Collection<
//                Collection<Array<Indx, 3, 2>>,
//                Collection<Array<Indx, 3, 1>>
//        >;
//
//    public:

//        auto const static constexpr node_connectivity = Connectivity(
//                Collection(
//                        Array<Indx, 3, 2>{
//                                0, 1,
//                                1, 2,
//                                2, 0
//                        }
//                ),
//                Collection(
//                        Array<Indx, 3, 1>{
//                                0,
//                                1,
//                                2
//                        }
//                )
//        );

        auto const static constexpr node_connectivity = Components<detail::ElementNodeConnectivity>{
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

//        auto const static constexpr node_con = Components2<detail::ElementNodeConnectivity>{
//                {
//                        {0, 1},
//                        {1, 2},
//                        {2, 0},
//                },
//                {
//                        {0},
//                        {1},
//                        {2},
//                }
//        };

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

    template<Element E>
    requires(E == qua_04)
    struct ElementGeometry<E>
    {

        auto const static constexpr element = E;

        constexpr
        Bool
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        /**
         * @brief
         */
        auto const static constexpr reference_nodes = Array<Real, 4, 2>{
            -1.0000000000000000, -1.0000000000000000,
            +1.0000000000000000, -1.0000000000000000,
            +1.0000000000000000, +1.0000000000000000,
            -1.0000000000000000, +1.0000000000000000
        };

        template<template<Element, auto...> typename T>
        using Components = Collection<
                Collection<
                        Array<T<seg_02>, 4>
                >,
                Collection<
                        Array<T<pnt_00>, 4>
                >
        >;

//        template<template<Element, auto...> typename T>
//        using Components2 = Collection<
//                Collection<
//                        T<seg_02>,
//                        T<seg_02>,
//                        T<seg_02>,
//                        T<seg_02>
//                >,
//                Collection<
//                        T<pnt_00>,
//                        T<pnt_00>,
//                        T<pnt_00>,
//                        T<pnt_00>
//                >
//        >;
//
//    private:
//
//        using Connectivity = Collection<
//                Collection<Array<Indx, 4, 2>>,
//                Collection<Array<Indx, 4, 1>>
//        >;
//
//    public:

//        auto const static constexpr node_connectivity = Connectivity(
//                Collection(
//                        Array<Indx, 4, 2>{
//                                0, 1,
//                                1, 2,
//                                2, 3,
//                                3, 0
//                        }
//                ),
//                Collection(
//                        Array<Indx, 4, 1>{
//                                0,
//                                1,
//                                2,
//                                3
//                        }
//                )
//        );

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

//        auto const static constexpr node_con = Components2<detail::ElementNodeConnectivity>{
//                {
//                        {0, 1},
//                        {1, 2},
//                        {2, 3},
//                        {3, 0},
//                },
//                {
//                        {0},
//                        {1},
//                        {2},
//                        {3},
//                }
//        };

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

    template<Element E>
    requires(E == tet_04)
    struct ElementGeometry<E>
    {

        auto const static constexpr element = E;

        constexpr
        Bool
        operator==(
                ElementGeometry const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                ElementGeometry const &
                other
        )
        const = default;

        /**
         * @brief
         */
        auto const static constexpr reference_nodes = Array<Real, 4, 3>{
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
            +1.0000000000000000, +0.0000000000000000, +0.0000000000000000,
            +0.0000000000000000, +1.0000000000000000, +0.0000000000000000,
            +0.0000000000000000, +0.0000000000000000, +1.0000000000000000
        };

        template<template<Element, auto...> typename T>
        using Components = Collection<
                Collection<
                        Array<T<tri_03>, 4>
                >,
                Collection<
                        Array<T<seg_02>, 6>
                >,
                Collection<
                        Array<T<pnt_00>, 4>
                >
        >;

//        template<template<Element, auto...> typename T>
//        using Components2 = Collection<
//                Collection<
//                        T<tri_03>,
//                        T<tri_03>,
//                        T<tri_03>,
//                        T<tri_03>
//                >,
//                Collection<
//                        T<seg_02>,
//                        T<seg_02>,
//                        T<seg_02>,
//                        T<seg_02>,
//                        T<seg_02>,
//                        T<seg_02>
//                >,
//                Collection<
//                        T<pnt_00>,
//                        T<pnt_00>,
//                        T<pnt_00>,
//                        T<pnt_00>
//                >
//        >;
//
//    private:
//
//        using Connectivity = Collection<
//                Collection<Array<Indx, 4, 3>>,
//                Collection<Array<Indx, 6, 2>>,
//                Collection<Array<Indx, 4, 1>>
//        >;
//
//    public:

//        auto const static constexpr node_connectivity = Connectivity(
//                Collection(
//                        Array<Indx, 4, 3>{
//                                0, 1, 3,
//                                0, 3, 2,
//                                0, 2, 1,
//                                1, 2, 3
//                        }
//                ),
//                Collection(
//                        Array<Indx, 6, 2>{
//                                0, 1,
//                                1, 2,
//                                2, 0,
//                                0, 3,
//                                3, 2,
//                                1, 3
//                        }
//                ),
//                Collection(
//                        Array<Indx, 4, 1>{
//                                0,
//                                1,
//                                2,
//                                3
//                        }
//                )
//        );

        auto const static constexpr node_connectivity = Components<detail::ElementNodeConnectivity>{
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

//        auto const static constexpr node_con = Components2<detail::ElementNodeConnectivity>{
//                {
//                        {0, 1, 3},
//                        {0, 3, 2},
//                        {0, 2, 1},
//                        {1, 2, 3},
//                },
//                {
//                        {0, 1},
//                        {1, 2},
//                        {2, 0},
//                        {0, 3},
//                        {3, 2},
//                        {1, 3},
//                },
//                {
//                        {0},
//                        {1},
//                        {2},
//                        {3},
//                }
//        };

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

    template<Indx D, template<Element> typename T>
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

        template<template<Element, auto...> typename T, auto I, auto J, Element E, auto... A>
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

        template<auto... K>
        auto const static constexpr layers = std::array<Indx, sizeof...(K)>{static_cast<Indx>(K)...};

        template<auto... K>
        struct Pol;

        template<auto... K>
        requires(sizeof...(K) == 0)
        struct Pol<K...>
        {

            using Type = ElemCom::Type;

        };

        template<auto... K>
        requires(sizeof...(K) == 1)
        struct Pol<K...>
        {

            using Type = typename ElemCom::Type::template Type<layers<K...>[0]>;

        };

        template<auto... K>
        requires(sizeof...(K) == 2)
        struct Pol<K...>
        {

            using Type = typename ElemCom::Type::template Type<layers<K...>[0]>::template Type<layers<K...>[1]>::Type;

        };

    public:

        template<auto... K>
        using Lay = typename Pol<K...>::Type;

    };

    template<Element E, template<Element, auto...> typename T, auto... A>
    using ElementComponents = typename detail::ElementComponentsPolicy<T, E, A...>::Type;

    template<Element E, Indx I, Indx J, template<Element, auto...> typename T, auto... A>
    using ElementComponent = typename detail::ElementComponentPolicy<T, I, J, E, A...>::Type;

    template<Element E, Indx D, template<Element, auto...> typename T, auto... A>
    using ElementNeighbourArray = typename detail::ElementNeighbourArrayPolicy<E, D, T, A...>::Type;

    template<Element E>
    static constexpr
    auto
    elementIndex()
    {
        return collection::index<typename Elements<3, ElementGeometry>::template Type<E.dim>, ElementGeometry<E>>();
    }

    template<Element E, Indx I, Indx J>
    constexpr inline
    auto
    component()
    {
        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
        return L::template Type<I>::template Type<J>::Type::element;
    }

    template<Element E, Indx I, Indx J>
    constexpr inline
    auto
    numComponents()
    {
        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
        return L::template Type<I>::template Type<J>::size();
    }

    template<Element E, Indx I>
    constexpr inline
    auto
    numComponents()
    {
        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
        return L::template Type<I>::size();
    }

    template<Element E>
    constexpr inline
    auto
    numComponents()
    {
        using L = typename ElementGeometry<E>::template Components<ElementGeometry>;
        return L::size();
    }

    template<Element E, Quadrature Q>
    struct ElementQuadrature;

    template<Element E, Quadrature Q>
    requires(Point<E>)
    struct ElementQuadrature<E, Q>
    {

        auto const static constexpr dim_quadrature = 1;

        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                +0.0000000000000000
        };

        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +1.0000000000000000
        };

    };

    template<Element E, Quadrature Q>
    requires(Segment<E> && Q.quadrature == QuadratureRule::Gauss && Q.ord == 1)
    struct ElementQuadrature<E, Q>
    {

        auto const static constexpr dim_quadrature = 1;

        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                +0.0000000000000000
        };

        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +2.0000000000000000
        };

    };

    template<Element E, Quadrature Q>
    requires(Segment<E> && Q.quadrature == QuadratureRule::Gauss && Q.ord == 2)
    struct ElementQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 2;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                -0.5773502691896257,
                +0.5773502691896257
        };

        /**
         * @brief
         */
        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +1.0000000000000000,
                +1.0000000000000000
        };

    };

    template<Element E, Quadrature Q>
    requires(Segment<E> && Q.quadrature == QuadratureRule::Gauss && Q.ord == 3)
    struct ElementQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 3;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                -0.7745966692414834,
                +0.0000000000000000,
                +0.7745966692414834
        };

        /**
         * @brief
         */
        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +0.5555555555555557,
                +0.8888888888888888,
                +0.5555555555555557
        };

    };

    template<Element E, Quadrature Q>
    requires(Segment<E> && Q.quadrature == QuadratureRule::Gauss && Q.ord == 4)
    struct ElementQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 3;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                -0.7745966692414834,
                +0.0000000000000000,
                +0.7745966692414834
        };

        /**
         * @brief
         */
        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +0.5555555555555557,
                +0.8888888888888888,
                +0.5555555555555557
        };

    };

    template<Element E, Quadrature Q>
    requires(Triangle<E> && Q.quadrature == QuadratureRule::Gauss && Q.ord == 1)
    struct ElementQuadrature<E, Q>
    {

        auto const static constexpr dim_quadrature = 1;

        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                +0.0000000000000000, +0.0000000000000000
        };

        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +2.0000000000000000
        };

    };

    template<Element E, Quadrature Q>
    requires(Triangle<E> && Q.quadrature == QuadratureRule::Gauss && Q.ord == 2)
    struct ElementQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 2;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000
        };

        /**
         * @brief
         */
        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +1.0000000000000000,
                +1.0000000000000000
        };

    };

    template<Element E, Quadrature Q>
    requires(Triangle<E> && Q.quadrature == QuadratureRule::Gauss && Q.ord == 3)
    struct ElementQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 3;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000,
        };

        /**
         * @brief
         */
        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +0.5555555555555557,
                +0.8888888888888888,
                +0.5555555555555557
        };

    };

    template<Element E, Quadrature Q>
    requires(Triangle<E> && Q.quadrature == QuadratureRule::Gauss && Q.ord == 4)
    struct ElementQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 3;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000,
        };

        /**
         * @brief
         */
        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +0.5555555555555557,
                +0.8888888888888888,
                +0.5555555555555557
        };

    };

    template<Element E, Quadrature Q>
    requires(Quadrangle<E>)
    struct ElementQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 3;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000,
        };

        /**
         * @brief
         */
        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +0.5555555555555557,
                +0.8888888888888888,
                +0.5555555555555557
        };

    };

    template<Element E, Quadrature Q>
    requires(Tetrahedron<E>)
    struct ElementQuadrature<E, Q>
    {

        /**
         * @brief
         */
        auto const static constexpr dim_quadrature = 3;

        /**
         * @brief
         */
        auto const static constexpr reference_points = Array<Real, dim_quadrature, E.dimPoint()>{
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
        };

        /**
         * @brief
         */
        auto const static constexpr reference_weights = Array<Real, dim_quadrature>{
                +0.5555555555555557,
                +0.8888888888888888,
                +0.5555555555555557
        };

    };

    template<Element E, Quadrature Q>
    constexpr inline
    auto
    dimQuadrature()
    {
        return ElementQuadrature<E, Q>::dim_quadrature;
    }

}

#endif //LOLITA_LOLITA_CORE_ELEMENT_GEOMETRY_HXX
