//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
#define LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX

#include "lolita/lolita_pointers.hxx"
#include "lolita/lolita_matrix.hxx"

//#include "lolita/lolita_core_element_quadrature.hxx"
#include "lolita/lolita_core_element_geometry.hxx"
//#include "lolita/lolita_core_element_basis.hxx"
//#include "lolita/lolita_core_element_degree_of_freedom.hxx"

namespace lolita::core::element
{

    namespace sand_box
    {

        template<Indx D>
        struct GenericFiniteElement
        {

            virtual
            ~GenericFiniteElement() = default;

            virtual
            UniquePointer<GenericFiniteElement>
            copy()
            const = 0;

//        virtual
//        Matrix<Real>
//        getCurrentCoordinates()
//        const = 0;

        };

        template<Indx D, typename Derived>
        struct FiniteElement2 : public GenericFiniteElement<D>
        {

        protected:

            // We make clear Shape class needs to be inherited
            FiniteElement2() = default;

            FiniteElement2(FiniteElement2 const &) = default;

            FiniteElement2(FiniteElement2 &&) = default;

        public:

            UniquePointer<GenericFiniteElement<D>>
            copy()
            const override
            {
                return UniquePointer<GenericFiniteElement<D>>::template make<Derived>(* static_cast<Derived const *>(this));
            }

//        Matrix<Real>
//        getCurrentCoordinates()
//        const override
//        {
//            Matrix<Real> a = static_cast<Derived const *>(this)->getCurrentCoordinates();
//            return a;
//        }

        };

        template<Indx D>
        struct ShapeImpl : public FiniteElement2<D, ShapeImpl<D>>
        {

            ShapeImpl() = default;

//        auto
//        getCurrentCoordinates()
//        const
//        {
//            return Matrix<Real, 2, 2>().setZero();
//        }

        };

        template<Indx D>
        struct ShapeImpl2 : public FiniteElement2<D, ShapeImpl2<D>>
        {

            ShapeImpl2() = default;

//        auto
//        getCurrentCoordinates()
//        const
//        {
//            return Matrix<Real, 4, 4>().setZero();
//        }

        };

    }

    /*
     * CONN
     */

    template<Element E, Domain D, auto C>
    struct FiniteElementConnectivity;

    template<Element E, Domain D, auto C>
    requires(!Point<E>)
    struct FiniteElementConnectivity<E, D, C> : public ElementGeometry<E>
    {

    private:

        using Self = FiniteElementConnectivity;

        template<Element Eb>
        using FiniteElementPointer = SharedPointer<FiniteElement<Eb, D, C>>;

    public:

        template<Element Eb>
        struct Component : public FiniteElementPointer<Eb>
        {

        private:

            using Base = Self::FiniteElementPointer<Eb>;

        public:

            Component()
            :
            Base(),
            orientation()
            {}

            Component(
                    Base const &
                    base_arg,
                    Intg
                    orientation_arg
            )
            :
            Base(base_arg),
            orientation(orientation_arg)
            {}

            Component(
                    typename Base::Type const &
                    element_arg,
                    Intg
                    orientation_arg
            )
            :
            Base(element_arg),
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

            Intg orientation;

        };

        using Components = ElementComponents<E, Component>;

        using Neighbours = ElementNeighbourArray<E, D.dim, FiniteElementPointer>;

        FiniteElementConnectivity()
        :
        tag(),
        components(),
        neighbours()
        {}

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
                current_nodes_coordinates.col(i) = nodes.get(i).get().getCurrentCoordinates();
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

        auto
        getHash()
        const
        {
            StrgStream ss;
            auto const & nodes = components.template get<E.dim - 1>().template get<0>();
            for (Indx i = 0; i < E.num_nodes; ++i) {
                ss << nodes.get(i).get().getHash();
            }
            return ss.str();
        }

        Indx tag;

        Components components;

        Neighbours neighbours;

    };

    template<Element E, Domain D, auto C>
    requires(Point<E>)
    struct FiniteElementConnectivity<E, D, C> : public ElementGeometry<E>
    {

    private:

        using Self = FiniteElementConnectivity;

        template<Element N>
        using FiniteElementPointer = SharedPointer<FiniteElement<N, D, C>>;

    public:

        //using Coordinates = Vector<Real, D.dim>;
        using Coordinates = SharedPointer<Vector<Real, D.dim>>;

        using Neighbours = ElementNeighbourArray<E, D.dim, FiniteElementPointer>;

        FiniteElementConnectivity()
        :
        tag(),
        coordinates(),
        neighbours()
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

        auto
        getHash()
        const
        {
            StrgStream ss;
            ss << tag + 1;
            return ss.str();
        }

        Indx tag;

        Coordinates coordinates;

        Neighbours neighbours;

    };

    /*
     * BASIS
     */

    template<Element E, Basis B>
    struct EB;

    template<Element E, Basis B>
    requires(B.basis == BasisName::Monomial)
    struct EB<E, B>
    {

        auto const static constexpr dim_basis = numerics::binomial(E.dim + B.ord, E.dim);

        template<auto D, auto M>
        struct Implementation : public FiniteElement<E, D, M>
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
        return EB<E, B>::dim_basis;
    }

    /*
     * BUILD
     */

    template<Element E, Domain D, auto M>
    struct FiniteElementGeometry : public FiniteElementConnectivity<E, D, M>
    {

    private:

        using Self = element::FiniteElementGeometry<E, D, M>;

    public:

        FiniteElementGeometry()
        :
        FiniteElementConnectivity<E, D, M>()
        {}

        template<Basis B>
        auto
        getBasisEvaluation(
                auto const &
                point_arg
        )
        const
        {
            using FiniteElementBasis = typename element::EB<E, B>::template Implementation<D, M>;
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
            //using FiniteElementBasis = typename finite_element::FiniteElementBasis<E, B>::template Implementation<M>;
            using FiniteElementBasis = typename element::EB<E, B>::template Implementation<D, M>;
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
            using ComponentGeometry = FiniteElementGeometry<component<E, I, J>(), D, M>;
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
            using ComponentGeometry = FiniteElementGeometry<component<E, I, J>(), D, M>;
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
//            auto derivative_vector = Matrix<Real, 3, Self::dim_nodes>::Zero();
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
    struct FiniteElementUnknown2 : public FiniteElementGeometry<E, D, F>
    {

        auto const static constexpr field = Field(F.ord_field, D.dim);

        using Module = typename FiniteElementPolicy2<E, D, F>::Module;

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

            using Coefficients = Vector<Real, dimUnknowns()>;

            DegreeOfFreedom()
            requires(StaticMatrixType<Coefficients>)
                    :
                    coefficients(Coefficients::Zero()),
                    index(0)
            {}

            DegreeOfFreedom()
            requires(DynamicMatrixType<Coefficients>)
                    :
                    coefficients(Coefficients::Zero(0)),
                    index(0)
            {}

            explicit
            DegreeOfFreedom(
                    RealType auto &
                    index_arg
            )
            requires(StaticMatrixType<Coefficients>)
                    :
                    coefficients(Coefficients::Zero()),
                    index(index_arg)
            {
                index_arg += dimUnknowns();
            }

            DegreeOfFreedom(
                    RealType auto
                    size_arg,
                    RealType auto &
                    index_arg
            )
            requires(DynamicMatrixType<Coefficients>)
                    :
                    coefficients(Coefficients::Zero(size_arg)),
                    index(index_arg)
            {
                index_arg += size_arg;
            }

            Coefficients coefficients;

            Indx index;

        };

        FiniteElementUnknown2()
        :
        FiniteElementGeometry<E, D, F>(),
        unknowns(setUnknowns()),
        bindings(),
        loads()
        {}

    private:

        static
        auto
        setUnknowns()
        {
            auto dof = Array<UniquePointer<DegreeOfFreedom>, field.rows(), field.cols()>();
            for (int i = 0; i < field.rows(); ++i) {
                for (int j = 0; j < field.cols(); ++j) {
                    dof.get(i, j) = UniquePointer<DegreeOfFreedom>(DegreeOfFreedom());
                }
            }
            return dof;
        }

    public:

        Array<UniquePointer<DegreeOfFreedom>, field.rows(), field.cols()> unknowns;

        Array<UniquePointer<DegreeOfFreedom>, field.rows(), field.cols()> bindings;

        Array<SharedPointer<LoadComponent<D>>, field.rows(), field.cols()> loads;

    };

    template<Element E, Domain D, auto F>
    requires(Cell<E, D>)
    struct FiniteElementUnknownMod : public FiniteElementUnknown2<E, D, F>
    {

    };

    /*
     * Unknowns
     */

    template<Element E, Domain D, auto F, auto M, auto C>
    struct FiniteElementPolicy;

    template<Element E, Domain D, auto F, auto M, auto C>
    struct FiniteElementBase;

    template<Element E, Domain D, auto F, auto M, auto C>
    struct FiniteElementUnknown
    {

        auto const static constexpr mixed_element_index = aggregate::index(C, M);

        auto const static constexpr finite_element_index = aggregate::index(M, F);

        auto const static constexpr field = Field(F.ord_field, D.dim);

        static constexpr
        auto
        dimUnknowns()
        {
            return numerics::max(-1, FiniteElementPolicy<E, D, F, M, C>::dimUnknowns());
        }

        static constexpr
        auto
        numUnknowns()
        {
            return numerics::max(-1, dimUnknowns() * field.size());
        }

        static
        auto
        fetch(
                auto &
                finite_element_arg
        )
        {
            return finite_element_arg.template get<mixed_element_index>().template get<finite_element_index>();
        }

        struct DegreeOfFreedom
        {

            using Coefficients = Vector<Real, dimUnknowns()>;

            DegreeOfFreedom()
            requires(StaticMatrixType<Coefficients>)
            :
            coefficients(Coefficients::Zero()),
            index(0)
            {}

            DegreeOfFreedom()
            requires(DynamicMatrixType<Coefficients>)
            :
            coefficients(Coefficients::Zero(0)),
            index(0)
            {}

            explicit
            DegreeOfFreedom(
                    RealType auto &
                    index_arg
            )
            requires(StaticMatrixType<Coefficients>)
            :
            coefficients(Coefficients::Zero()),
            index(index_arg)
            {
                index_arg += dimUnknowns();
            }

            DegreeOfFreedom(
                    RealType auto
                    size_arg,
                    RealType auto &
                    index_arg
            )
            requires(DynamicMatrixType<Coefficients>)
            :
            coefficients(Coefficients::Zero(size_arg)),
            index(index_arg)
            {
                index_arg += size_arg;
            }

            Coefficients coefficients;

            Indx index;

        };

        FiniteElementUnknown()
        :
        unknowns(setUnknowns()),
        bindings(),
        loads()
        {}

    private:

        static
        auto
        setUnknowns()
        {
            auto dof = Array<UniquePointer<DegreeOfFreedom>, field.rows(), field.cols()>();
            for (int i = 0; i < field.rows(); ++i) {
                for (int j = 0; j < field.cols(); ++j) {
                    dof.get(i, j) = UniquePointer<DegreeOfFreedom>(DegreeOfFreedom());
                }
            }
            return dof;
        }

    public:

        Array<UniquePointer<DegreeOfFreedom>, field.rows(), field.cols()> unknowns;

        Array<UniquePointer<DegreeOfFreedom>, field.rows(), field.cols()> bindings;

        Array<SharedPointer<LoadComponent<D>>, field.rows(), field.cols()> loads;

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Cell<E, D>)
    struct FiniteElementBase<E, D, F, M, C> : public FiniteElementUnknown<E, D, F, M, C>
    {

    private:

        using Base = FiniteElementUnknown<E, D, F, M, C>;

        using Module = typename FiniteElementPolicy<E, D, F, M, C>::Module;

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
            for (auto m : F.mappings.data) {
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
            for (int i = 0; i < F.mappings.size(); ++i) {
                value += rowsOperator(F.mappings.get(i));
            }
            return value;
        }

        static constexpr
        auto
        colsOperator()
        {
            return FiniteElementPolicy<E, D, F, M, C>::colsOperator();
        }

        auto const static constexpr quadrature = Quadrature(C.quadrature_rule, D.ordIntegration(C.ordIntegration()));

        using Operators = Array<Matrix<Real, rowsOperator(), colsOperator()>, dimQuadrature<E, quadrature>()>;

//        FiniteElementBase()
//        :
//        Base()
//        {}

        Operators operators;

        Module module;

        void
        getGrads(
                auto const &
                finite_element_arg
        )
        {
            using Implementation = typename FiniteElementPolicy<E, D, F, M, C>::Implementation;
            auto res = static_cast<Implementation *>(this)->getGradients(finite_element_arg, 0);
            print("res");
            print(res);
            auto extf = static_cast<Implementation *>(this)->getExternalForces(finite_element_arg, 0);
            print("extf");
            print(extf);
        }

        void
        initialize(
                auto &
                finite_element_arg
        )
        {
            using Implementation = typename FiniteElementPolicy<E, D, F, M, C>::Implementation;
            static_cast<Implementation *>(this)->setModule(finite_element_arg);
            auto set_num_components = [&] <auto I> ()
            mutable
            {
                static_cast<Implementation *>(this)->template setOperator<F.mappings.get(I), I>(finite_element_arg);
            };
            collection::apply<F.mappings.size()>(set_num_components);
        }

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Face<E, D>)
    struct FiniteElementBase<E, D, F, M, C> : public FiniteElementUnknown<E, D, F, M, C>
    {

        void
        initialize(
                auto &
                finite_element_arg
        )
        {

        }

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Edge<E, D>)
    struct FiniteElementBase<E, D, F, M, C> : public FiniteElementUnknown<E, D, F, M, C>
    {

        void
        initialize(
                auto &
                finite_element_arg
        )
        {

        }

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Node<E, D>)
    struct FiniteElementBase<E, D, F, M, C> : public FiniteElementUnknown<E, D, F, M, C>
    {

        void
        initialize(
                auto &
                finite_element_arg
        )
        {

        }

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Cell<E, D> && F.method == fem_hho && M.model == Model::Solid)
    struct FiniteElementPolicy<E, D, F, M, C>
    {

        static constexpr
        auto
        dimUnknowns()
        {
            return Intg(EB<E, Basis(BasisName::Monomial, F.discretization.ord_cell)>::dim_basis);
        }

        static constexpr
        auto
        colsOperator()
        {
            auto in = Intg(0);
            auto set_num_components = [&] <auto K> () constexpr mutable {
                auto nj = FiniteElementBase<component<E, 0, K>(), D, F, M, C>::numUnknowns();
                auto num_c = numComponents<E, 0, K>();
                in += num_c * nj;
            };
            collection::apply<numComponents<E, 0>()>(set_num_components);
            return FiniteElementBase<E, D, F, M, C>::numUnknowns() + in;
        }

        struct Module
        {

            Matrix<Real, colsOperator(), colsOperator()> stabilization;

        };

        struct Implementation : public FiniteElementBase<E, D, F, M, C>
        {

        private:

            using Base = FiniteElementBase<E, D, F, M, C>;

        public:

            void
            setModule(
                    auto &
                    finite_element_arg
            )
            {
                this->module.stabilization.setOnes();
                this->module.stabilization(0, 0) = finite_element_arg.tag;
            }

            template<MappingOperator Map, auto I>
            void
            setOperator(
                    auto &
                    finite_element_arg
            );

            template<MappingOperator Map, auto I>
            void
            setOperator(
                    auto &
                    finite_element_arg
            )
            requires(Map == MappingOperator::Gradient)
            {
                auto const constexpr ord_max = numerics::max(
                        F.discretization.ordMapping(Map),
                        F.discretization.ord_cell,
                        F.discretization.ord_face
                );
                /*
                 * Defining constants
                 */
                auto const constexpr ord_qad_grd = D.ordIntegration(ord_max);
                auto const constexpr qad_grd = Quadrature(C.quadrature_rule, ord_qad_grd);
                auto const constexpr dim_qad_grd = dimQuadrature<E, qad_grd>();
                auto const constexpr bas_grd = Basis(BasisName::Monomial, F.discretization.ordMapping(Map));
                auto const constexpr bas_cel = Basis(BasisName::Monomial, F.discretization.ord_cell);
                auto const constexpr bas_fce = Basis(BasisName::Monomial, F.discretization.ord_face);
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
                        auto pt = finite_element_arg.template getReferenceQuadraturePoint<qad_grd>(i);
                        auto wt = finite_element_arg.template getCurrentQuadratureWeight<qad_grd>(i);
                        auto vr = finite_element_arg.template getBasisEvaluation<bas_grd>(pt);
                        lhs += wt * vr * vr.transpose();
                    }
                    //auto mymat1 = Matrix<Real, dim_bas_grd, dim_bas_grd>().setIdentity();
                    //auto res1 = mymat1.llt().solve(Matrix<Real, dim_bas_grd, dim_bas_grd>::Identity());
                    //print(mymat1);
                    //print(res1);
                    //print("lhs before inv");
                    //print(lhs);
                    //Matrix<Real, dim_bas_grd, dim_bas_grd> rhs = Matrix<Real, dim_bas_grd, dim_bas_grd>::Identity();
                    //Matrix<Real, dim_bas_grd, dim_bas_grd> res = lhs.llt().solve(rhs);
                    //print("lhs after inv");
                    //print(res);
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
                    auto set_rhs_cell = [&] (
                            auto
                            i_arg,
                            auto
                            j_arg
                    )
                    mutable
                    {
                        for (int i = 0; i < dim_qad_grd; ++i) {
                            auto pc = finite_element_arg.template getReferenceQuadraturePoint<qad_grd>(i);
                            auto wc = finite_element_arg.template getCurrentQuadratureWeight<qad_grd>(i);
                            auto vr = finite_element_arg.template getBasisEvaluation<bas_grd>(pc);
                            auto vc = finite_element_arg.template getBasisDerivative<bas_cel>(pc, j_arg);
                            auto prt = rhs.template block<dim_bas_grd, dim_bas_cel>(0, dim_bas_cel * i_arg);
                            prt += wc * vr * vc.transpose();
                        }
                    };
                    /*
                     * Defining the face offset
                     */
                    auto get_faces_offset = [] <auto N> ()
                    constexpr
                    {
                        auto in = FiniteElementBase<E, D, F, M, C>::numUnknowns();
                        auto set_num_components = [&] <auto K> () constexpr mutable {
                            auto nj = FiniteElementBase<component<E, 0, K - 1>(), D, F, M, C>::numUnknowns();
                            auto num_c = numComponents<E, 0, K - 1>();
                            in += num_c * nj;
                        };
                        if constexpr (N > 0) collection::apply<N>(set_num_components);
                        return in;
                    };
                    /*
                     * Setting the jump part of the gradient operator, at a given face
                     */
                    auto set_rhs_face = [&] <auto K> (
                            auto
                            i_arg,
                            auto
                            j_arg,
                            auto
                            num_f
                    )
                    mutable
                    {
                        auto const constexpr dim_bas_fce = dimBasis<component<E, 0, K>(), bas_fce>();
                        auto const & face = finite_element_arg.template getComponentPointer<0, K>(num_f).get();
                        auto const & n_dir = finite_element_arg.template getComponentPointer<0, K>(num_f).orientation;
                        for (int i = 0; i < dim_qad_grd; ++i) {
                            auto pf = face.template getReferenceQuadraturePoint<qad_grd>(i);
                            auto pc = finite_element_arg.template getComponentReferenceQuadraturePoint<qad_grd, 0, K>(num_f, i);
                            auto wf = face.template getCurrentQuadratureWeight<qad_grd>(i);
                            auto n = face.getNormalVector(pf);
                            auto vr = finite_element_arg.template getBasisEvaluation<bas_grd>(pc);
                            auto vc = finite_element_arg.template getBasisEvaluation<bas_cel>(pc);
                            auto vf = face.template getBasisEvaluation<bas_fce>(pf);
                            auto offset = dim_bas_cel * i_arg;
                            auto prt_cell = rhs.template block<dim_bas_grd, dim_bas_cel>(0, offset);
                            prt_cell -= wf * vr * vc.transpose() * n(j_arg) * n_dir;
                            offset = get_faces_offset.template operator()<K>();
                            offset += dim_bas_fce * (num_f * field<D, F>().size() + i_arg);
                            auto prt_face = rhs.template block<dim_bas_grd, dim_bas_fce>(0, offset);
                            prt_face += wf * vr * vf.transpose() * n(j_arg) * n_dir;
                        }
                    };
                    /*
                     * Setting the jump parts of the gradient operator
                     */
                    auto set_rhs_faces = [&] <auto K> (
                            auto
                            i_arg,
                            auto
                            j_arg
                    )
                    mutable
                    {
                        for (int i = 0; i < numComponents<E, 0, K>(); ++i) {
                            set_rhs_face.template operator()<K>(i_arg, j_arg, i);
                        }
                    };
                    set_rhs_cell(i_arg, j_arg);
                    collection::apply<numComponents<E, 0>()>(set_rhs_faces, i_arg, j_arg);
                    return rhs;
                };
                /*
                 *
                 */
                auto const constexpr fld = Field(F.ord_field, D.dim);
                auto const constexpr grd = Field::fromMapping(Base::field, Map);
                auto row = Indx(0);
                auto lhs = get_lhs();
                for (int i = 0; i < grd.rows(); ++i) {
                    for (int j = 0; j < grd.cols(); ++j) {
                        auto rhs = get_rhs(i, j);
                        for (int k = 0; k < dimQuadrature<E, this->quadrature>(); ++k) {
                            auto pnt = finite_element_arg.template getReferenceQuadraturePoint<this->quadrature>(k);
                            auto vct = finite_element_arg.template getBasisEvaluation<bas_grd>(pnt);
                            this->operators.get(k).row(this->rowOperator(Map) + row) = vct.transpose() * lhs * rhs;
                        }
                        row += 1;
                    }
                }
            }

            template<MappingOperator Map, auto I>
            void
            setOperator(
                    auto &
                    finite_element_arg
            )
            requires(Map == MappingOperator::Identity)
            {
                /*
                 *
                 */
                auto const constexpr bas_ide = Basis(BasisName::Monomial, F.discretization.ordMapping(Map));
                auto const constexpr dim_bas_ide = dimBasis<E, bas_ide>();
                auto const constexpr bas_cel = Basis(BasisName::Monomial, F.discretization.ord_cell);
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
                    auto set_rhs_cell = [&] (
                            auto
                            i_arg,
                            auto
                            j_arg
                    )
                    mutable
                    {
                        auto prt = rhs.template block<dim_bas_ide, dim_bas_cel>(0, dim_bas_cel * i_arg);
                        prt = Matrix<Real, dim_bas_ide, dim_bas_cel>::Identity();
                    };
                    set_rhs_cell(i_arg, j_arg);
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
                            auto pnt = finite_element_arg.template getReferenceQuadraturePoint<this->quadrature>(k);
                            auto vct = finite_element_arg.template getBasisEvaluation<bas_ide>(pnt);
                            this->operators.get(k).row(this->rowOperator(Map) + row) = vct.transpose() * rhs;
                        }
                        row += 1;
                    }
                }
            }

            auto
            getUnknowns(
                    auto const &
                    finite_element_arg
            )
            const
            {
                auto offset = Indx(0);
                auto vec = Vector<Real, Base::colsOperator()>().setZero();
                auto get_cell_unknowns = [&] ()
                mutable
                {
                    auto const constexpr dim_cel_unk = Base::dimUnknowns();
                    auto const & cel = this->fetch(finite_element_arg);
                    for (int i = 0; i < Base::field.rows(); ++i) {
                        for (int j = 0; j < Base::field.cols(); ++j) {
                            vec.template segment<dim_cel_unk>(offset) = cel.unknowns.get(i, j).get().coefficients;
                            offset += dim_cel_unk;
                        }
                    }
                };
                auto get_face_unknowns = [&] <auto K> (auto num_f)
                mutable
                {
                    auto const constexpr dim_fce_unk = FiniteElementBase<component<E, 0, K>(), D, F, M, C>::dimUnknowns();
                    auto const & fce = this->fetch(finite_element_arg.template getComponentPointer<0, K>(num_f).get());
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
                    auto const &
                    finite_element_arg,
                    auto
                    index_arg
            )
            const
            {
                print("ops :");
                print(this->operators.get(index_arg));
                return Vector<Real, Base::rowsOperator()>(this->operators.get(index_arg) * getUnknowns(finite_element_arg));
            }

            auto
            getExternalForces(
                    auto const &
                    finite_element_arg,
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
                    auto const & cel = this->fetch(finite_element_arg);
                    for (int i = 0; i < this->field.rows(); ++i) {
                        for (int j = 0; j < this->field.cols(); ++j) {
                            for (int k = 0; k < dimQuadrature<E, this->quadrature>(); ++k) {
                                auto pnt = finite_element_arg.template getReferenceQuadraturePoint<this->quadrature>(k);
                                auto pntc = finite_element_arg.template getCurrentQuadraturePoint<this->quadrature>(k);
                                auto wgt = finite_element_arg.template getCurrentQuadratureWeight<this->quadrature>(k);
                                auto vct = finite_element_arg.template getBasisEvaluation<bas_cel>(pnt);
                                auto val = Real();
                                if (cel.loads.get(i, j).exists()) {
                                    val = cel.loads.get(i, j).get().getImposedValue(pntc, time_arg);
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
                    auto const constexpr dim_fce_unk = FiniteElementBase<component<E, 0, K>(), D, F, M, C>::dimUnknowns();
                    auto const & fce = this->fetch(finite_element_arg.template getComponentPointer<0, K>(num_f).get());
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

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Face<E, D> && F.method == fem_hho && M.model == Model::Solid)
    struct FiniteElementPolicy<E, D, F, M, C>
    {

        static constexpr
        auto
        dimUnknowns()
        {
            return Intg(EB<E, Basis(BasisName::Monomial, F.discretization.ord_face)>::dim_basis);
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementBase<E, D, F, M, C>
        {

        };

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Edge<E, D> && F.method == fem_hho && M.model == Model::Solid)
    struct FiniteElementPolicy<E, D, F, M, C>
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

        struct Implementation : public FiniteElementBase<E, D, F, M, C>
        {

        };

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Node<E, D> && F.method == fem_hho && M.model == Model::Solid)
    struct FiniteElementPolicy<E, D, F, M, C>
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

        struct Implementation : public FiniteElementBase<E, D, F, M, C>
        {

        };

    };


    template<Element E, Domain D, auto M, auto C>
    struct MixedElementPolicy2 : public decltype(M)::template Elements2<FiniteElementBase, E, D, M, C>
    {

        SharedPointer<Behaviour> behaviour;

    };

    template<Element E, Domain D, auto C>
    struct CoupledElementPolicy2 : public decltype(C)::template Elements2<MixedElementPolicy2, E, D, C>
    {

    };

    template<Element E, Domain D, auto C>
    struct FiniteElement : public FiniteElementGeometry<E, D, C>, public CoupledElementPolicy2<E, D, C>
    {

        FiniteElement()
        :
        FiniteElementGeometry<E, D, C>(),
        CoupledElementPolicy2<E, D, C>()
        {}

        void
        initialize()
        {
            this->template get<0>().template get<0>().initialize(* this);
        }

        void
        getGrads()
        {
            this->template get<0>().template get<0>().getGrads(* this);
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

    };

}

#endif //LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
