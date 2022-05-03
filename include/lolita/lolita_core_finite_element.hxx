//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
#define LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX

#include "lolita/lolita_pointers.hxx"
#include "lolita/lolita_matrix.hxx"

#include "lolita/lolita_core_element_geometry.hxx"

namespace lolita::core::element
{

    template<Element E, Domain D, auto F>
    struct FiniteElementUnknownMod;

    template<Element E, Domain D, template<Element, Domain, auto...> typename T, auto... A>
    struct FiniteElementNeighbour
    {

        FiniteElementNeighbour()
        :
        ptr(),
        index()
        {}

        FiniteElementNeighbour(
                SharedPointer<T<E, D, A...>> const &
                ptr_arg
        )
        :
        ptr(ptr_arg),
        index(0)
        {}

        FiniteElementNeighbour(
                SharedPointer<T<E, D, A...>> const &
                ptr_arg,
                Indx
                index_arg
        )
        :
        ptr(ptr_arg),
        index(index_arg)
        {}

        Bool
        operator==(
                FiniteElementNeighbour const &
                other
        )
        const = default;

        Bool
        operator!=(
                FiniteElementNeighbour const &
                other
        )
        const = default;

        SharedPointer<T<E, D, A...>>   ptr;

        Indx index;

        Intg
        orientation()
        const
        {
            return index == 0 ? 1 : -1;
        }

    };

    template<Element E, Domain D, template<Element, Domain, auto...> typename T, auto... A>
    struct FiniteElementConnectivity;

    template<Element E, Domain D, template<Element, Domain, auto...> typename T, auto... A>
    struct FiniteElementGeometry;

    template<Element E, Domain D, template<Element, Domain, auto...> typename T, auto... A>
    requires(!Point<E>)
    struct FiniteElementConnectivity<E, D, T, A...> : public ElementGeometry<E>
    {

    private:

        using Self = FiniteElementConnectivity;

        template<Element Eb>
        using FiniteElementPointer = SharedPointer<T<Eb, D, A...>>;

        template<Element Eb>
        using Neighbourggg = FiniteElementNeighbour<Eb, D, T, A...>;

    public:

        template<Element Eb>
        struct Component
        {

            Component()
                    :
                    ptr(),
                    orientation()
            {}

            Component(
                    SharedPointer<T<Eb, D, A...>> const &
                    ptr_arg,
                    Intg
                    orientation_arg
            )
                    :
                    ptr(ptr_arg),
                    orientation(orientation_arg)
            {}

            Bool
            operator==(
                    Component const &
                    other
            )
            const = default;

            Bool
            operator!=(
                    Component const &
                    other
            )
            const = default;

            SharedPointer<T<Eb, D, A...>> ptr;

            Intg orientation;

        };

        template<Element Eb>
        struct Neighbour
        {

            Neighbour()
            :
            ptr()
            {}

            Neighbour(
                    SharedPointer<T<Eb, D, A...>> const &
                    ptr_arg
            )
            :
            ptr(ptr_arg)
            {}

            Bool
            operator==(
                    Neighbour const &
                    other
            )
            const = default;

            Bool
            operator!=(
                    Neighbour const &
                    other
            )
            const = default;

            SharedPointer<T<Eb, D, A...>>   ptr;
            Indx                            index;

        };

        template<Element Eb>
        struct Neighbour2
        {

            Neighbour2()
            :
            ptr(),
            index()
            {}

            Neighbour2(
                    SharedPointer<T<Eb, D, A...>> const &
                    ptr_arg,
                    Indx
                    index_arg
            )
            :
            ptr(ptr_arg),
            index(index_arg)
            {}

            Bool
            operator==(
                    Neighbour2 const &
                    other
            )
            const = default;

            Bool
            operator!=(
                    Neighbour2 const &
                    other
            )
            const = default;

            SharedPointer<T<Eb, D, A...>>   ptr;

            Indx const index;

        };

//        using Components = ElementComponents<E, Component>;
//
//        using Neighbours = ElementNeighbourArray<E, D.dim, Neighbour>;

        using Components = ElementComponents<E, Neighbourggg>;

        using Neighbours = ElementNeighbourArray<E, D.dim, Neighbourggg>;

        //using Adjacency = Collection<ElementComponents<E, Neighbourggg>, ElementNeighbourArray<E, D.dim, Neighbourggg>>;

        FiniteElementConnectivity()
        :
        components(),
        neighbours(),
        count(0)
        {}

        template<auto I, auto J>
        auto
        getComponentIndex(
                auto
                index_arg
        )
        const
        {
            auto const & cts = components.template get<I>().template get<J>().get(index_arg).ptr.get();
            auto const & klm = cts.neighbours.template get<I>().template get<J>();
            auto cnt = Indx(0);
            for (auto const & item : klm) {
                if (item.ptr.get() == * this) {

                }
                cnt ++;
            }
            return cnt;
        }

        template<auto I, auto J>
        auto
        getComponentPointer(
                auto
                index_arg
        )
        {
            return components.template get<I>().template get<J>().get(index_arg);
        }

        template<auto I, auto J>
        auto
        getComponentPointer(
                auto
                index_arg
        )
        const
        {
            return components.template get<I>().template get<J>().get(index_arg);
        }

        template<auto I, auto J>
        auto
        getNeighbourPointer(
                auto
                index_arg
        )
        {
            return neighbours.template get<I>().template get<J>().get(index_arg);
        }

        template<auto I, auto J>
        auto
        getNeighbourPointer(
                auto
                index_arg
        )
        const
        {
            return neighbours.template get<I>().template get<J>().get(index_arg);
        }

        auto
        getCurrentCoordinates()
        const
        {
            auto current_nodes_coordinates = Matrix<Real, D.dim, E.num_nodes>();
            auto const & nodes = components.template get<E.dim - 1>().template get<0>();
            for (Indx i = 0; i < E.num_nodes; ++i) {
                current_nodes_coordinates.col(i) = nodes.get(i).ptr.get().getCurrentCoordinates();
            }
            return current_nodes_coordinates;
        }

        static
        auto
        getReferenceCoordinates()
        {
            using ReferenceCoordinates = matrix::MatMap<Matrix<Real, E.dimPoint(), E.num_nodes, matrix::col_major> const>;
            return ReferenceCoordinates(Self::reference_nodes.data.data());
        }

        Components components;

        Neighbours neighbours;

        Indx count;

    };

    template<Element E, Domain D, template<Element, Domain, auto...> typename T, auto... A>
    requires(Point<E>)
    struct FiniteElementConnectivity<E, D, T, A...> : public ElementGeometry<E>
    {

    private:

        using Self = FiniteElementConnectivity;

        template<Element Eb>
        using FiniteElementPointer = SharedPointer<T<Eb, D, A...>>;

        template<Element Eb>
        using Neighbourggg = FiniteElementNeighbour<Eb, D, T, A...>;

    public:

        using Coordinates = SharedPointer<Vector<Real, D.dim>>;

        template<Element Eb>
        struct Neighbour
        {

            Neighbour()
            :
            ptr()
            {}

            Neighbour(
                    SharedPointer<T<Eb, D, A...>> const &
                    ptr_arg
            )
            :
            ptr(ptr_arg)
            {}

            Bool
            operator==(
                    Neighbour const &
                    other
            )
            const = default;

            Bool
            operator!=(
                    Neighbour const &
                    other
            )
            const = default;

            SharedPointer<T<Eb, D, A...>> ptr;

        };

        //using Neighbours = ElementNeighbourArray<E, D.dim, Neighbour>;

        using Neighbours = ElementNeighbourArray<E, D.dim, Neighbourggg>;

        FiniteElementConnectivity()
        :
        coordinates(),
        neighbours(),
        count(0)
        {}

        template<auto I, auto J>
        auto &
        getNeighbour(
                auto
                index_arg
        )
        {
            return neighbours.template get<I>().template get<J>().get(index_arg).get();
        }

        template<auto I, auto J>
        auto const &
        getNeighbour(
                auto
                index_arg
        )
        const
        {
            return neighbours.template get<I>().template get<J>().get(index_arg).get();
        }

        auto
        getCurrentCoordinates()
        {
            //return coordinates;
            return coordinates.get();
        }

        auto
        getCurrentCoordinates()
        const
        {
            //return coordinates;
            return coordinates.get();
        }

        static
        auto
        getReferenceCoordinates()
        {
            return Vector<Real, 1>{0.0};
        }

        Coordinates coordinates;

        Neighbours neighbours;

        Indx count;

    };

    /*
     * BASIS
     */

    template<Element E, Basis B>
    struct FiniteElementBasis;

    template<Element E, Basis B>
    requires(B.basis == BasisName::Monomial)
    struct FiniteElementBasis<E, B>
    {

        auto const static constexpr dim_basis = numerics::binomial(E.dim + B.ord, E.dim);

//        template<typename T>
//        struct Implementation : public T
        template<Domain D, template<Element, Domain, auto...> typename T, auto... A>
        struct Implementation : public FiniteElementGeometry<E, D, T, A...>
        {

        private:

            using Self = Implementation;

            static constexpr
            auto
            setExponents()
            {
                auto exponents_values = Array<Indx, dim_basis, E.dimPoint()>();
                if constexpr (E.dim == 0) {
                    exponents_values.get(0, 0) = 0;
                }
                else if constexpr (E.dim == 1) {
                    auto row = Indx(0);
                    for (Indx i = 0; i < B.ord + 1; ++i) {
                        exponents_values.get(row, 0) = i;
                        row += 1;
                    }
                }
                else if constexpr (E.dim == 2) {
                    auto row = Indx(0);
                    for (Indx i = 0; i < B.ord + 1; ++i) {
                        for (Indx j = 0; j < i + 1; ++j) {
                            exponents_values.get(row, 0) = i - j;
                            exponents_values.get(row, 1) = j;
                            row += 1;
                        }
                    }
                }
                else if constexpr (E.dim == 3) {
                    auto row = Indx(0);
                    for (Indx i = 0; i < B.ord + 1; ++i) {
                        for (Indx j = 0; j < i + 1; ++j) {
                            for (Indx k = 0; k < i + 1; ++k) {
                                if (j + k < i + 1) {
                                    exponents_values.get(row, 0) = i - (j + k);
                                    exponents_values.get(row, 1) = k;
                                    exponents_values.get(row, 2) = j;
                                    row += 1;
                                }
                            }
                        }
                    }
                }
                return exponents_values;
            }

            auto const static constexpr exponents = setExponents();

        public:

            auto
            evaluate(
                    auto const &
                    point_arg
            )
            const
            {
                auto basis_vector_values = Vector<Real, dim_basis>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (Indx i = 0; i < dim_basis; ++i) {
                    auto value = Real(1);
                    for (Indx j = 0; j < E.dim; ++j) {
                        auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                        //value *= numerics::pow(2.0 * dist / diameters(j), exponents.get(i, j));
                        value *= std::pow(2.0 * dist / diameters(j), exponents.get(i, j));
                    }
                    basis_vector_values(i) = value;
                }
                return basis_vector_values;
            }

            auto
            evaluate(
                    auto const &
                    point_arg,
                    auto
                    derivative_direction_arg
            )
            const
            {
                auto basis_vector_values = Vector<Real, dim_basis>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (Indx i = 0; i < dim_basis; ++i) {
                    auto value = Real(1);
                    for (Indx j = 0; j < E.dim; ++j) {
                        if (j != derivative_direction_arg) {
                            auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                            value *= std::pow(2.0 * (dist) / diameters(j), exponents.get(i, j));
                        }
                        else {
                            if (exponents.get(i, j) > 0) {
                                auto c = 2.0 * exponents.get(i, j) / diameters(j);
                                auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                                //value *= c * std::pow(2.0 * (dist) / diameters(j), exponents.get(i, j) - 1);
                                value *= c * numerics::pow(2.0 * (dist) / diameters(j), exponents.get(i, j) - 1);
                            }
                            else {
                                value *= 0.0;
                            }
                        }
                    }
                    basis_vector_values(i) = value;
                }
                return basis_vector_values;
            }

        };

    };

    template<Element E, Basis B>
    static constexpr
    auto
    dimBasis()
    {
        return FiniteElementBasis<E, B>::dim_basis;
    }

    /*
     * BUILD
     */

    template<Element E, Domain D, template<Element, Domain, auto...> typename T, auto... A>
    struct FiniteElementGeometry : public FiniteElementConnectivity<E, D, T, A...>
    {

    private:

        using Self = element::FiniteElementGeometry<E, D, T, A...>;

    public:

        FiniteElementGeometry()
        :
        FiniteElementConnectivity<E, D, T, A...>()
        {}

        template<Basis B>
        auto
        getBasisEvaluation(
                auto const &
                point_arg
        )
        const
        {
//            using FiniteElementBasis = typename EB<E, B>::template Implementation<FiniteElementGeometry>;
            using FiniteElementBasis = typename FiniteElementBasis<E, B>::template Implementation<D, T, A...>;
            return static_cast<FiniteElementBasis const *>(this)->evaluate(point_arg);
        }

        template<Basis B>
        auto
        getBasisDerivative(
                auto const &
                point_arg,
                auto
                derivative_direction_arg
        )
        const
        {
//            using FiniteElementBasis = typename EB<E, B>::template Implementation<FiniteElementGeometry>;
            using FiniteElementBasis = typename FiniteElementBasis<E, B>::template Implementation<D, T, A...>;
            return static_cast<FiniteElementBasis const *>(this)->evaluate(point_arg, derivative_direction_arg);
        }

        auto
        getShapeMappingDifferential(
                auto const &
                point_arg
        )
        const
        {
            auto const nds = this->getCurrentCoordinates();
            auto ru = Matrix<Real, 3, E.dimPoint()>().setZero();
            auto du = Real(0);
            for (Indx i = 0; i < D.dim; ++i) {
                for (Indx j = 0; j < E.dimPoint(); ++j) {
                    ru(i, j) = Self::getShapeMappingDerivative(nds.row(i), point_arg, j);
                }
            }
            if constexpr (E.dimPoint() == 1) {
                du = numerics::abs(ru.col(0).norm());
            }
            else if constexpr (E.dimPoint() == 2) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else if constexpr (E.dimPoint() == 3) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            if constexpr (D.frame == EuclideanFrame::AxiSymmetric) {
                Real r0 = getShapeMappingEvaluation(nds.row(0), point_arg);
                if (r0 < 1.e-10) {
                    r0 = 1.e-10;
                }
                du *= 2.0 * numerics::pi * r0;
            }
            return du;
        }

        auto
        getDistanceInCurrentConfiguration(
                auto const &
                first_point_arg,
                auto const &
                second_point_arg,
                auto
                direction_arg = Intg(-1)
        )
        const
        {
            assert(-1 <= direction_arg <= E.dim);
            auto distance = Real();
            if constexpr (Cell<E, D>) {
                distance = getCartesianDistance(first_point_arg, second_point_arg, direction_arg);
            }
            else {
                distance = getCurvedDistance(first_point_arg, second_point_arg, direction_arg);
            }
            return distance;
        }

        auto
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = Self::getReferenceCoordinates();
            auto current_diameters = Vector<Real, E.dim>().setZero();
            for (Indx i = 0; i < E.num_nodes; ++i) {
                for (Indx j = i + 1; j < E.num_nodes; ++j) {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (Indx k = 0; k < E.dimPoint(); ++k) {
                        auto new_value = numerics::abs(getDistanceInCurrentConfiguration(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }

        auto
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            return geometry::getBarycenter(current_nodes_coordinates);
        }

        static
        auto
        getReferenceDistance(
                auto const &
                first_point_arg,
                auto const &
                second_point_arg,
                auto
                direction_arg = Intg(-1)
        )
        {
            assert((-1 <= direction_arg <= static_cast<Intg>(E.dim)));
            if (direction_arg == -1) {
                return (second_point_arg - first_point_arg).norm();
            } else {
                return (second_point_arg - first_point_arg)(direction_arg);
            }
        }

//        static
//        Real
//        getReferenceDistance(
//                auto const &
//                first_point_arg,
//                auto const &
//                second_point_arg
//        )
//        {
//            Real distance;
//            distance = (second_point_arg - first_point_arg).norm();
//            return distance;
//        }

        static
        auto
        getReferenceCentroid()
        {
            auto nds = Self::getReferenceCoordinates();
            return geometry::getBarycenter(nds);
        }

        static
        auto
        getReferenceDiameters()
        {
            auto dts = Vector<Real, E.dim>().setZero();
            auto nds = Self::getReferenceCoordinates();
            for (Indx i = 0; i < E.num_nodes; ++i) {
                for (Indx j = i + 1; j < E.num_nodes; ++j) {
                    for (Indx k = 0; k < E.dim; ++k) {
                        auto & a = dts(k);
                        auto b = numerics::abs(nds(k, i) - nds(k, j));
                        if (b > a) {
//                            reference_diameters(k) = b;
                            a = b;
                        }
                    }
                }
            }
            return dts;
        }

        /*
         * QUADRATURE
         */

        template<Quadrature Q>
        static
        auto
        getReferenceQuadratureWeight(
                auto
                index_arg
        )
        {
            return ElementQuadrature<E, Q>::reference_weights.get(index_arg);
        }

        template<Quadrature Q>
        static
        auto
        getReferenceQuadraturePoint(
                auto
                index_arg
        )
        {
            auto offset = E.dimPoint() * index_arg;
            //Real const constexpr * const d = ShapeQuadrature<E.shape_description, Q>::reference_points.data.data();
            auto const constexpr * const d = ElementQuadrature<E, Q>::reference_points.data.data();
            return matrix::MatMap<Vector<Real, E.dimPoint()> const>(d + offset);
        }

        template<Quadrature Q, Indx I, Indx J>
        static
        auto
        getComponentReferenceQuadratureWeight(
                Indx
                index_arg
        )
        {
            using ComponentGeometry = FiniteElementGeometry<component<E, I, J>(), D, T, A...>;
            return ComponentGeometry::template getReferenceQuadratureWeight<Q>(index_arg);
        }

        template<Quadrature Q, Indx I, Indx J>
        static
        auto
        getComponentReferenceQuadraturePoint(
                auto
                component_index_arg,
                auto
                index_arg
        )
        {
            auto p = Vector<Real, E.dimPoint()>();
            using ComponentGeometry = FiniteElementGeometry<component<E, I, J>(), D, T, A...>;
            auto const & cpt_node_tags = ElementGeometry<E>::node_connectivity.template get<I>().template get<J>();
            auto const & elt_reference_nodes = ElementGeometry<E>::reference_nodes;
            for (Indx i = 0; i < E.dimPoint(); ++i) {
                auto cpt_coordinates = Vector<Real, component<E, I, J>().num_nodes>();
                for (Indx j = 0; j < component<E, I, J>().num_nodes; ++j) {
                    auto const node_tag = cpt_node_tags.get(component_index_arg).get(j);
                    cpt_coordinates(j) = elt_reference_nodes.get(node_tag, i);
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<Q>(index_arg);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }

        template<Quadrature Q>
        auto
        getCurrentQuadratureWeight(
                auto
                index_arg
        )
        const
        {
            auto w = getReferenceQuadratureWeight<Q>(index_arg);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<Q>(index_arg));
        }

        template<Quadrature Q>
        auto
        getCurrentQuadraturePoint(
                auto
                index_arg
        )
        const
        {
            auto p = Vector<Real, D.dim>();
            auto const nds = this->getCurrentCoordinates();
            for (Indx j = 0; j < D.dim; ++j) {
                p(j) = Self::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<Q>(index_arg));
            }
            return p;
        }

        template<Quadrature Q, Indx I, Indx J>
        auto
        getComponentCurrentQuadratureWeight(
                auto
                component_index_arg,
                Indx
                index_arg
        )
        {
            auto const & cmp =  this->components.template get<I>().template get<J>().get(component_index_arg).get();
            return cmp.template getCurrentQuadratureWeight<Q>(index_arg);
        }

        template<Quadrature Q, Indx I, Indx J>
        auto
        getComponentCurrentQuadraturePoint(
                auto
                component_index_arg,
                auto
                index_arg
        )
        const
        {
            auto p = Vector<Real, D.dim>();
            auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<Q, I, J>(component_index_arg, index_arg);
            auto const nds = this->getCurrentCoordinates();
            for (Indx j = 0; j < D.dim; ++j) {
                p(j) = Self::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

    private:

        auto
        getCurvedDistance(
                auto const &
                first_point_arg,
                auto const &
                second_point_arg,
                auto
                direction_arg = Intg(-1)
        )
        const
        requires(!Point<E>)
        {
            auto const constexpr q_seg = Quadrature(QuadratureRule::Gauss, 4);
            auto distance = Real(0);
            auto dt = Real();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (Indx q = 0; q < dimQuadrature<seg_02, q_seg>(); ++q) {
                auto pq = ElementQuadrature<seg_02, q_seg>::reference_points.get(q, 0);
                auto wq = ElementQuadrature<seg_02, q_seg>::reference_weights.get(q);
                auto ru = Matrix<Real, D.dim, E.dim>();
                auto difference = second_point_arg - first_point_arg;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                for (Indx i = 0; i < D.dim; ++i) {
                    for (Indx j = 0; j < E.dim; ++j) {
                        if (direction_arg == -1 || i == static_cast<Indx>(direction_arg)) {
                            auto du = (1.0 / 2.0) * (second_point_arg(j) - first_point_arg(j));
                            auto dx = Self::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                            ru(i, j) = dx * du;
                        }
                    }
                }
                if constexpr (Segment<E>) {
                    auto Eff = ru.col(0).template dot(ru.col(0));
                    dt = std::sqrt(Eff);
                }
                else if constexpr (Surface<E>) {
                    auto Eff = ru.col(0).template dot(ru.col(0));
                    auto Fff = ru.col(0).template dot(ru.col(1));
                    auto Gff = ru.col(1).template dot(ru.col(1));
                    dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                }
                distance += wq * dt;
            }
            return distance;
        }

        auto
        getCartesianDistance(
                auto const &
                first_point_arg,
                auto const &
                second_point_arg,
                auto
                direction_arg = Intg(-1)
        )
        const
        {
            auto const & nds = this->getCurrentCoordinates();
            auto distance = Real();
            auto mp0 = Vector<Real, D.dim>();
            auto mp1 = Vector<Real, D.dim>();
            for (Indx i = 0; i < D.dim; ++i) {
                mp0(i) = Self::getShapeMappingEvaluation(nds.row(i), first_point_arg);
                mp1(i) = Self::getShapeMappingEvaluation(nds.row(i), second_point_arg);
            }
            direction_arg == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction_arg);
            return distance;
        }

    public:

        Vector<Real, D.dim>
        getNormalVector(
                auto const &
                point_arg
        )
        const
        requires(Face<E, D>)
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto derivative_vector = Matrix<Real, 3, E.dimPoint()>().setZero();
            for (Indx i = 0; i < D.dim; ++i) {
                for (Indx j = 0; j < E.dimPoint(); ++j) {
                    derivative_vector(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, j);
                }
            }
            auto const & block = derivative_vector.template block<D.dim, E.dimPoint()>(0, 0);
            auto normal_vector = geometry::getNormalVector(block);
            return normal_vector;
        }

        Vector<Real, D.dim>
        getTangentVector(
                auto const &
                point_arg,
                Indx
                direction_arg
        )
        const
        requires(Edge<E, D>)
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = Vector<Real, D.dim>();
            for (Indx i = 0; i < D.dim; ++i) {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, direction_arg);
            }
            return tangent_vector;
        }

    };

    template<Element E, Domain D, auto F>
    struct FiniteElementPolicy2;

    template<Element E, Domain D, auto F>
    struct FiniteElementUnknownF
    {

        auto const static constexpr field = Field(F.field().ord_field, D.dim);

        static constexpr
        auto
        dimUnknowns()
        {
            return numerics::max(-1, FiniteElementPolicy2<E, D, F>::dimUnknowns());
        }

        static constexpr
        auto
        numUnknowns()
        {
            return numerics::max(-1, dimUnknowns() * field.size());
        }

        struct DegreeOfFreedom
        {

            DegreeOfFreedom()
            :
            coefficients(),
            index(-1)
            {}

            Vector<Real> coefficients;
            
            Intg index;

        };

        FiniteElementUnknownF()
        :
        unknowns(),
        bindings(),
        loads()
        {}

        void
        setUnknowns()
        {
            if constexpr (dimUnknowns() > 0) {
                for (int i = 0; i < field.rows(); ++i) {
                    for (int j = 0; j < field.cols(); ++j) {
                        unknowns.get(i, j).coefficients.resize(dimUnknowns());
                        unknowns.get(i, j).coefficients.setZero();
                    }
                }
            }
        }

        void
        setUnknowns(
                auto &
                degree_of_freedom_index_arg
        )
        {
            if constexpr (dimUnknowns() > 0) {
                for (int i = 0; i < field.rows(); ++i) {
                    for (int j = 0; j < field.cols(); ++j) {
                        unknowns.get(i, j).coefficients.resize(dimUnknowns());
                        unknowns.get(i, j).coefficients.setZero();
                        unknowns.get(i, j).index = degree_of_freedom_index_arg;
                        degree_of_freedom_index_arg += dimUnknowns();
                    }
                }
            }
        }

        void
        setLoad(
                auto &
                degree_of_freedom_index_arg,
                auto const &
                load_pointer_arg,
                auto
                row_arg,
                auto
                col_arg
        )
        {
            if (load_pointer_arg.get().load_type == LoadType::Natural) {
                loads.get(row_arg, col_arg) = load_pointer_arg;
            }
            else {
                loads.get(row_arg, col_arg) = load_pointer_arg;
                bindings.get(row_arg, col_arg).coefficients.resize(dimUnknowns());
                bindings.get(row_arg, col_arg).coefficients.setZero();
                bindings.get(row_arg, col_arg).index = degree_of_freedom_index_arg;
                degree_of_freedom_index_arg += dimUnknowns();
            }
        }

        template<Indx S>
        auto
        getUnknownCoefficients(
                std::integral auto
                row_arg,
                std::integral auto
                col_arg
        )
        const
        {
            return unknowns.get(row_arg, col_arg).coefficients.template segment<S>(0);
        };

        template<Indx S>
        auto
        getBindingCoefficients(
                auto
                row_arg,
                auto
                col_arg
        )
        const
        {
            return bindings.get(row_arg, col_arg).coefficients.template segment<S>(0);
        };

        Array<DegreeOfFreedom, field.rows(), field.cols()> unknowns;

        Array<DegreeOfFreedom, field.rows(), field.cols()> bindings;

        Array<SharedPointer<LoadComponent<D>>, field.rows(), field.cols()> loads;

    };

    template<Element E, Domain D, auto F>
    requires(Cell<E, D>)
    struct FiniteElementUnknownMod<E, D, F> : public FiniteElementUnknownF<E, D, F>, public FiniteElementGeometry<E, D, FiniteElementUnknownMod, F>
    {

    private:

        using Base = FiniteElementUnknownF<E, D, F>;

        using Base2 = FiniteElementGeometry<E, D, FiniteElementUnknownMod, F>;

        using Module = typename FiniteElementPolicy2<E, D, F>::Module;

    public:

        static constexpr
        auto
        rowsOperator(
                MappingOperator
                App
        )
        {
            return Field::fromMapping(Base::field, App).size();
        }

        static constexpr
        auto
        rowOperator(
                MappingOperator
                App
        )
        {
            auto value = Indx(0);
            for (auto m : F.field().mappings.data) {
                if (m == App) {
                    return value;
                }
                value += Field::fromMapping(Base::field, m).size();
            }
            return value;
        }

        static constexpr
        auto
        rowsOperator()
        {
            auto value = Indx(0);
            for (int i = 0; i < F.field().mappings.size(); ++i) {
                value += rowsOperator(F.field().mappings.get(i));
            }
            return value;
        }

        static constexpr
        auto
        colsOperator()
        {
            return FiniteElementPolicy2<E, D, F>::colsOperator();
        }

        auto const static constexpr quadrature = Quadrature(QuadratureRule::Gauss, 1);

        using Operators = Array<Matrix<Real, rowsOperator(), colsOperator()>, dimQuadrature<E, quadrature>()>;

        using MaterialPoints = Array<SharedPointer<mgis::behaviour::BehaviourData>, dimQuadrature<E, quadrature>()>;

        FiniteElementUnknownMod()
        :
        Base(),
        Base2()
        {}

        void
        getGrads()
        const
        {
            using Implementation = typename FiniteElementPolicy2<E, D, F>::Implementation;
            auto res = static_cast<Implementation const *>(this)->getGradients(0);
            print("res");
            print(res);
            auto extf = static_cast<Implementation const *>(this)->getExternalForces(0);
            print("extf");
            print(extf);
        }

        void
        setUnknowns(
                auto &
                degree_of_freedom_index_arg
        )
        {
            reinterpret_cast<FiniteElementUnknownF<E, D, F> *>(this)->setUnknowns();
        }

        void
        setMaterial(
                auto const &
                behaviour_arg
        )
        {
            for (int i = 0; i < dimQuadrature<E, quadrature>(); ++i) {
                material_points.get(i) = SharedPointer<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(behaviour_arg));
            }
        }

        void
        initialize()
        {
            using Implementation = typename FiniteElementPolicy2<E, D, F>::Implementation;
            static_cast<Implementation *>(this)->setModule();
            auto set_operator = [&] <auto I = 0> (auto & self)
            mutable
            {
                static_cast<Implementation *>(this)->template setOperator<F.field().mappings.get(I), I>();
                if constexpr (I < F.field().mappings.size() - 1) {
                    self.template operator()<I + 1>(self);
                }
            };
            set_operator(set_operator);
        }

        Operators operators;

        Module module;

        MaterialPoints material_points;

    };

    template<Element E, Domain D, auto F>
    requires(Face<E, D>)
    struct FiniteElementUnknownMod<E, D, F>
    :
    public FiniteElementUnknownF<E, D, F>,
    public FiniteElementGeometry<E, D, FiniteElementUnknownMod, F>
    {

    private:

        using Base = FiniteElementUnknownF<E, D, F>;

        using Base2 = FiniteElementGeometry<E, D, FiniteElementUnknownMod, F>;

        using Module = typename FiniteElementPolicy2<E, D, F>::Module;

    public:

        FiniteElementUnknownMod()
        :
        Base(),
        Base2()
        {}

        void
        initialize()
        {
            // this->setUnknowns(degree_of_freedom_index_arg, load_pointer_arg);
        }

        void
        setMaterial(
                auto const &
                behaviour_arg
        )
        {

        }

        void
        setUnknowns(
                auto &
                degree_of_freedom_index_arg
        )
        {
            reinterpret_cast<FiniteElementUnknownF<E, D, F> *>(this)->setUnknowns(degree_of_freedom_index_arg);
        }

    };

    template<Element E, Domain D, auto F>
    requires(Edge<E, D>)
    struct FiniteElementUnknownMod<E, D, F> : public FiniteElementUnknownF<E, D, F>,  public FiniteElementGeometry<E, D, FiniteElementUnknownMod, F>
    {

    private:

        using Base = FiniteElementUnknownF<E, D, F>;

        using Base2 = FiniteElementGeometry<E, D, FiniteElementUnknownMod, F>;

        using Module = typename FiniteElementPolicy2<E, D, F>::Module;

    public:

        FiniteElementUnknownMod()
                :
                Base(),
                Base2()
        {}

        void
        initialize()
        {
            
        }

        void
        setMaterial(
                auto const &
                behaviour_arg
        )
        {

        }

        void
        setUnknowns(
                auto &
                degree_of_freedom_index_arg
        )
        {
            
        }

    };

    template<Element E, Domain D, auto F>
    requires(Node<E, D>)
    struct FiniteElementUnknownMod<E, D, F> : public FiniteElementUnknownF<E, D, F>,  public FiniteElementGeometry<E, D, FiniteElementUnknownMod, F>
    {

    private:

        using Base = FiniteElementUnknownF<E, D, F>;

        using Base2 = FiniteElementGeometry<E, D, FiniteElementUnknownMod, F>;

        using Module = typename FiniteElementPolicy2<E, D, F>::Module;

    public:

        FiniteElementUnknownMod()
                :
                Base(),
                Base2()
        {}

        void
        initialize()
        {
            
        }

        void
        setMaterial(
                auto const &
                behaviour_arg
        )
        {

        }

        void
        setUnknowns(
                auto &
                degree_of_freedom_index_arg
        )
        {
            
        }

    };

    template<Element E, Domain D, auto F>
    requires(Cell<E, D> && F.method == fem_hho)
    struct FiniteElementPolicy2<E, D, F>
    {

    private:

        using Base = FiniteElementUnknownMod<E, D, F>;

    public:

        static constexpr
        auto
        dimUnknowns()
        {
            return Intg(FiniteElementBasis<E, Basis(BasisName::Monomial, F.discretization().ord_cell)>::dim_basis);
        }

        static constexpr
        auto
        colsOperator()
        {
            auto in = Intg(0);
            auto set_num_components = [&] <auto K> () constexpr mutable {
                auto nj = FiniteElementUnknownMod<component<E, 0, K>(), D, F>::numUnknowns();
                auto num_c = numComponents<E, 0, K>();
                in += num_c * nj;
            };
            collection::apply<numComponents<E, 0>()>(set_num_components);
            return FiniteElementUnknownMod<E, D, F>::numUnknowns() + in;
        }

        struct Module
        {

            Matrix<Real, colsOperator(), colsOperator()> stabilization;

        };

        struct Implementation : public FiniteElementUnknownMod<E, D, F>
        {

            void
            setModule()
            {
                this->module.stabilization.setOnes();
                this->module.stabilization(0, 0) = 3;
            }

            template<MappingOperator Map, auto I>
            void
            setOperator();

            template<MappingOperator Map, auto I>
            void
            setOperator()
            requires(Map == MappingOperator::Gradient)
            {
                auto const constexpr ord_max = numerics::max(F.discretization().ordMapping(Map), F.discretization().ord_cell, F.discretization().ord_face);
                /*
                 * Defining constants
                 */
                auto const constexpr ord_qad_grd = D.ordIntegration(ord_max);
                auto const constexpr qad_grd = Quadrature(QuadratureRule::Gauss, ord_qad_grd);
                auto const constexpr dim_qad_grd = dimQuadrature<E, qad_grd>();
                auto const constexpr bas_grd = Basis(BasisName::Monomial, F.discretization().ordMapping(Map));
                auto const constexpr bas_cel = Basis(BasisName::Monomial, F.discretization().ord_cell);
                auto const constexpr bas_fce = Basis(BasisName::Monomial, F.discretization().ord_face);
                auto const constexpr dim_bas_grd = dimBasis<E, bas_grd>();
                auto const constexpr dim_bas_cel = dimBasis<E, bas_cel>();
                /*
                 * LHS
                 */
                auto get_lhs = [&] ()
                mutable
                {
                    auto lhs = Matrix<Real, dim_bas_grd, dim_bas_grd>().setZero();
                    for (int i = 0; i < dim_qad_grd; ++i) {
                        auto pt = this->template getReferenceQuadraturePoint<qad_grd>(i);
                        auto wt = this->template getCurrentQuadratureWeight<qad_grd>(i);
                        auto vr = this->template getBasisEvaluation<bas_grd>(pt);
                        lhs += wt * vr * vr.transpose();
                    }
                    lhs = lhs.llt().solve(decltype(lhs)::Identity());
                    return lhs;
                };
                /*
                 * RHS
                 */
                auto get_rhs = [&] (
                        auto
                        i_arg,
                        auto
                        j_arg
                )
                mutable
                {
                    /*
                     * Initializing the RHS part of the gradient operator
                     */
                    auto rhs = Matrix<Real, dim_bas_grd, colsOperator()>().setZero();
                    /*
                     * Setting the cell part of the gradient operator
                     */
                    auto set_rhs_cell = [&] ()
                    mutable
                    {
                        for (int i = 0; i < dim_qad_grd; ++i) {
                            auto pc = this->template getReferenceQuadraturePoint<qad_grd>(i);
                            auto wc = this->template getCurrentQuadratureWeight<qad_grd>(i);
                            auto vr = this->template getBasisEvaluation<bas_grd>(pc);
                            auto vc = this->template getBasisDerivative<bas_cel>(pc, j_arg);
                            auto prt = rhs.template block<dim_bas_grd, dim_bas_cel>(0, dim_bas_cel * i_arg);
                            prt += wc * vr * vc.transpose();
                        }
                    };
                    /*
                     * Defining the face offset
                     */
                    auto get_faces_offset2 = [] <auto K> ()
                    constexpr
                    {
                        auto offset = FiniteElementUnknownMod<E, D, F>::numUnknowns();
                        auto get_faces_offset22 = [&] <auto L> (auto & self)
                        constexpr mutable
                        {
                            if constexpr (L > 0) {
                                auto nj = FiniteElementUnknownMod<component<E, 0, L - 1>(), D, F>::numUnknowns();
                                auto num_c = numComponents<E, 0, L - 1>();
                                offset += num_c * nj;
                            }
                        };
                        get_faces_offset22.template operator()<K>(get_faces_offset22);
                        return offset;
                    };
                    /*
                     * Setting the jump part of the gradient operator, at a given face
                     */
                    auto set_rhs_face = [&] <auto K = 0> (
                            auto &
                            self
                    )
                    mutable
                    {
                        auto const constexpr dim_bas_fce = dimBasis<component<E, 0, K>(), bas_fce>();
                        for (int num_f = 0; num_f < numComponents<E, 0, K>(); ++num_f) {
                            auto const & face = this->template getComponentPointer<0, K>(num_f).ptr.get();
                            auto const & n_dir = this->template getComponentPointer<0, K>(num_f).orientation();
                            for (int i = 0; i < dim_qad_grd; ++i) {
                                auto pf = face.template getReferenceQuadraturePoint<qad_grd>(i);
                                auto pc = this->template getComponentReferenceQuadraturePoint<qad_grd, 0, K>(num_f, i);
                                auto wf = face.template getCurrentQuadratureWeight<qad_grd>(i);
                                auto n = face.getNormalVector(pf);
                                auto vr = this->template getBasisEvaluation<bas_grd>(pc);
                                auto vc = this->template getBasisEvaluation<bas_cel>(pc);
                                auto vf = face.template getBasisEvaluation<bas_fce>(pf);
                                auto offset = dim_bas_cel * i_arg;
                                auto prt_cell = rhs.template block<dim_bas_grd, dim_bas_cel>(0, offset);
                                prt_cell -= wf * vr * vc.transpose() * n(j_arg) * n_dir;
                                //offset = get_faces_offset.template operator()<K>();
                                offset = get_faces_offset2.template operator()<K>();
                                offset += dim_bas_fce * (num_f * field<D, F>().size() + i_arg);
                                auto prt_face = rhs.template block<dim_bas_grd, dim_bas_fce>(0, offset);
                                prt_face += wf * vr * vf.transpose() * n(j_arg) * n_dir;
                            }
                        }
                        if constexpr (K < numComponents<E, 0>() - 1) {
                            self.template operator()<K + 1>(i_arg, j_arg, self);
                        }
                    };
                    set_rhs_face(set_rhs_face);
                    return rhs;
                };
                /*
                 *
                 */
                auto const constexpr fld = Field(F.field().ord_field, D.dim);
                auto const constexpr grd = Field::fromMapping(Base::field, Map);
                auto row = Indx(0);
                auto lhs = get_lhs();
                for (int i = 0; i < grd.cols(); ++i) {
                    for (int j = 0; j < grd.rows(); ++j) {
                        auto rhs = get_rhs(i, j);
                        for (int k = 0; k < dimQuadrature<E, this->quadrature>(); ++k) {
                            auto pnt = this->template getReferenceQuadraturePoint<this->quadrature>(k);
                            auto vct = this->template getBasisEvaluation<bas_grd>(pnt);
                            this->operators.get(k).row(this->rowOperator(Map) + row) = vct.transpose() * lhs * rhs;
                        }
                        row += 1;
                    }
                }
            }

            template<MappingOperator Map, auto I>
            void
            setOperator()
            requires(Map == MappingOperator::Identity)
            {
                /*
                 *
                 */
                auto const constexpr bas_ide = Basis(BasisName::Monomial, F.discretization().ordMapping(Map));
                auto const constexpr dim_bas_ide = dimBasis<E, bas_ide>();
                auto const constexpr bas_cel = Basis(BasisName::Monomial, F.discretization().ord_cell);
                auto const constexpr dim_bas_cel = dimBasis<E, bas_cel>();
                /*
                 *
                 */
                auto get_rhs = [&] (
                        auto
                        i_arg,
                        auto
                        j_arg
                )
                        mutable
                {
                    /*
                     * Initializing the RHS part of the identity operator
                     */
                    auto rhs = Matrix<Real, dim_bas_ide, colsOperator()>().setZero();
                    /*
                     * Setting the cell part of the identity operator
                     */
                    auto set_rhs_cell = [&] ()
                            mutable
                    {
                        auto prt = rhs.template block<dim_bas_ide, dim_bas_cel>(0, dim_bas_cel * i_arg);
                        prt = Matrix<Real, dim_bas_ide, dim_bas_cel>::Identity();
                    };
                    set_rhs_cell();
                    return rhs;
                };
                /*
                 *
                 */
                auto const constexpr grd = Field::fromMapping(Base::field, Map);
                auto row = Indx(0);
                for (int i = 0; i < grd.rows(); ++i) {
                    for (int j = 0; j < grd.cols(); ++j) {
                        auto rhs = get_rhs(i, j);
                        for (int k = 0; k < dimQuadrature<E, this->quadrature>(); ++k) {
                            auto pnt = this->template getReferenceQuadraturePoint<this->quadrature>(k);
                            auto vct = this->template getBasisEvaluation<bas_ide>(pnt);
                            this->operators.get(k).row(this->rowOperator(Map) + row) = vct.transpose() * rhs;
                        }
                        row += 1;
                    }
                }
            }

            auto
            getUnknowns()
            const
            {
                auto offset = Indx(0);
                auto vec = Vector<Real, Base::colsOperator()>().setZero();
                auto get_cell_unknowns = [&] ()
                        mutable
                {
                    auto const constexpr dim_cel_unk = Base::dimUnknowns();
                    //auto const & cel = this->fetch(finite_element_arg);
                    for (int i = 0; i < Base::field.rows(); ++i) {
                        for (int j = 0; j < Base::field.cols(); ++j) {
                            vec.template segment<dim_cel_unk>(offset) = this->unknowns.get(i, j).get().coefficients;
                            offset += dim_cel_unk;
                        }
                    }
                };
                auto get_face_unknowns = [&] <auto K> (auto num_f)
                        mutable
                {
                    auto const constexpr dim_fce_unk = FiniteElementUnknownMod<component<E, 0, K>(), D, F>::dimUnknowns();
                    //auto const & fce = this->fetch(finite_element_arg.template getComponentPointer<0, K>(num_f).get());
                    auto const & fce = this->template getComponentPointer<0, K>(num_f).get();
                    for (int i = 0; i < Base::field.rows(); ++i) {
                        for (int j = 0; j < Base::field.cols(); ++j) {
                            vec.template segment<dim_fce_unk>(offset) = fce.unknowns.get(i, j).get().coefficients;
                            offset += dim_fce_unk;
                        }
                    }
                };
                auto get_faces_unknowns = [&] <auto K> ()
                        mutable
                {
                    for (int i = 0; i < numComponents<E, 0, K>(); ++i) {
                        get_face_unknowns.template operator()<K>(i);
                    }
                };
                get_cell_unknowns();
                collection::apply<numComponents<E, 0>()>(get_faces_unknowns);
                return vec;
            }

            auto
            getGradients(
                    auto
                    index_arg
            )
            const
            {
                print("ops :");
                print(this->operators.get(index_arg));
                return Vector<Real, Base::rowsOperator()>(this->operators.get(index_arg) * getUnknowns());
            }

            auto
            getExternalForces(
                    auto const &
                    time_arg
            )
            const
            {
                auto const constexpr bas_cel = Basis(BasisName::Monomial, F.discretization.ord_cell);
                auto const constexpr bas_fce = Basis(BasisName::Monomial, F.discretization.ord_face);
                auto const constexpr dim_bas_cel = dimBasis<E, bas_cel>();
                auto const constexpr dim_cel_unk = Base::dimUnknowns();
                auto offset = Indx(0);
                auto vec = Vector<Real, Base::colsOperator()>().setZero();
                auto get_cell_unknowns = [&] ()
                        mutable
                {
                    //auto const & cel = this->fetch(finite_element_arg);
                    for (int i = 0; i < this->field.rows(); ++i) {
                        for (int j = 0; j < this->field.cols(); ++j) {
                            for (int k = 0; k < dimQuadrature<E, this->quadrature>(); ++k) {
                                auto pnt = this->template getReferenceQuadraturePoint<this->quadrature>(k);
                                auto pntc = this->template getCurrentQuadraturePoint<this->quadrature>(k);
                                auto wgt = this->template getCurrentQuadratureWeight<this->quadrature>(k);
                                auto vct = this->template getBasisEvaluation<bas_cel>(pnt);
                                auto val = Real();
                                if (this->loads.get(i, j).exists()) {
                                    val = this->loads.get(i, j).get().getImposedValue(pntc, time_arg);
                                }
                                else {
                                    val = 0;
                                }
                                vec.template segment<dim_bas_cel>(offset) += val * wgt * vct;
                            }
                            offset += dim_cel_unk;
                        }
                    }
                };
                auto get_face_unknowns = [&] <auto K> (auto num_f)
                        mutable
                {
                    auto const constexpr dim_fce_unk = FiniteElementUnknownMod<component<E, 0, K>(), D, F>::dimUnknowns();
                    //auto const & fce = this->fetch(finite_element_arg.template getComponentPointer<0, K>(num_f).get());
                    auto const & fce = this->template getComponentPointer<0, K>(num_f).get();
                    for (int i = 0; i < Base::field.rows(); ++i) {
                        for (int j = 0; j < Base::field.cols(); ++j) {
                            vec.template segment<dim_fce_unk>(offset) = fce.unknowns.get(i, j).get().coefficients;
                            offset += dim_fce_unk;
                        }
                    }
                };
                auto get_faces_unknowns = [&] <auto K> ()
                        mutable
                {
                    for (int i = 0; i < numComponents<E, 0, K>(); ++i) {
                        get_face_unknowns.template operator()<K>(i);
                    }
                };
                get_cell_unknowns();
                collection::apply<numComponents<E, 0>()>(get_faces_unknowns);
                return vec;
            }

        };

    };

    template<Element E, Domain D, auto F>
    requires(Face<E, D> && F.method == fem_hho)
    struct FiniteElementPolicy2<E, D, F>
    {

        static constexpr
        auto
        dimUnknowns()
        {
            return Intg(FiniteElementBasis<E, Basis(BasisName::Monomial, F.discretization().ord_face)>::dim_basis);
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementUnknownMod<E, D, F>
        {

        };

    };

    template<Element E, Domain D, auto F>
    requires(Edge<E, D> && F.method == fem_hho)
    struct FiniteElementPolicy2<E, D, F>
    {

        static constexpr
        auto
        dimUnknowns()
        {
            return Intg(0);
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementUnknownMod<E, D, F>
        {

        };

    };

    template<Element E, Domain D, auto F>
    requires(Node<E, D> && F.method == fem_hho)
    struct FiniteElementPolicy2<E, D, F>
    {

        static constexpr
        auto
        dimUnknowns()
        {
            return Intg(0);
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementUnknownMod<E, D, F>
        {

        };

    };

    /*
     * Unknowns
     */

    template<Element E, Domain D, auto C>
    using MixedElement = typename decltype(C)::template MixedElementType<FiniteElementUnknownMod, E, D>;

    template<Element E, Domain D, auto C>
    struct FiniteElement : public MixedElement<E, D, C>, public FiniteElementGeometry<E, D, FiniteElement, C>
    {

        FiniteElement()
        :
        MixedElement<E, D, C>(),
        FiniteElementGeometry<E, D, FiniteElement, C>()
        {}

        void
        initialize()
        {

        }

        void
        getGrads()
        {
            this->template get<0>().getGrads();
        }

        Bool
        operator==(
                FiniteElement const &
                other
        )
        const
        {
            return this->tag == other.tag;
        }

        Bool
        operator!=(
                FiniteElement const &
                other
        )
        const
        {
            return !(other == * this);
        }

        auto
        hash()
        const
        {
            StrgStream hash;
            if constexpr(Point<E>) {
                hash << tag + 1;
            }
            else {
                auto const & nodes = this->components.template get<E.dim - 1>().template get<0>();
                for (Indx i = 0; i < E.num_nodes; ++i) {
                    hash << nodes.get(i).ptr.get().hash();
                }
            }
            return hash.str();
        }

        Long tag;

    };

    template<auto I>
    struct ATest
    {

        SharedPointer<ATest<I - 1>> ptr;

        Indx j;

    };

    template<auto I>
    requires(I == 0)
    struct ATest<I>
    {

        SharedPointer<Indx> ptr;

        Indx j;

    };

//    template<Element E, Domain D, auto F>
//    struct FiniteElement : public FiniteElementUnknownMod<E, D, F>, public FiniteElementGeometry<E, D, FiniteElement, F>
//    {
//
//        FiniteElement()
//        :
//        FiniteElementUnknownMod<E, D, F>(),
//        tag()
//        {}
//
//        Indx tag;
//
//        auto
//        getHash()
//        const
//        {
//            StrgStream ss;
//            ss << tag + 1;
//            return ss.str();
//        }
//
//        void
//        initialize()
//        {
//            this->template get<0>().initialize();
//        }
//
//        void
//        getGrads()
//        {
//            this->template get<0>().getGrads();
//        }
//
//        Bool
//        operator==(
//                FiniteElement const &
//                other
//        )
//        const
//        {
//            return this->tag == other.tag;
//        }
//
//        Bool
//        operator!=(
//                FiniteElement const &
//                other
//        )
//        const
//        {
//            return !(other == * this);
//        }
//
//    };

}

#endif //LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
