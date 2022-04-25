//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
#define LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX

#include "lolita/lolita_pointers.hxx"
#include "lolita/lolita_matrix.hxx"

#include "lolita/lolita_core_element_quadrature.hxx"
#include "lolita/lolita_core_element_geometry.hxx"
#include "lolita/lolita_core_element_basis.hxx"
#include "lolita/lolita_core_element_degree_of_freedom.hxx"

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
            orientation(1)
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

        using DomainPointer = SharedPointer<core::MeshInteriorDomain<C>>;

        FiniteElementConnectivity()
        :
        tag(),
        components(Components()),
        neighbours(Neighbours()),
        domain(DomainPointer())
        {}

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
            using ReferenceCoordinates = Matrix<Matrix<Real, E.dimPoint(), E.num_nodes, matrix::col_major> const>;
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

        DomainPointer domain;

    };

    template<Element E, Domain D, auto M>
    requires(Point<E>)
    struct FiniteElementConnectivity<E, D, M> : public ElementGeometry<E>
    {

    private:

        using Self = FiniteElementConnectivity;

        template<Element N>
        using FiniteElementPointer = SharedPointer<FiniteElement<N, D, M>>;

    public:

        using Coordinates = Vector<Real, D.dim>;

        using Neighbours = ElementNeighbourArray<E, D.dim, FiniteElementPointer>;

        using DomainPointer = SharedPointer<core::MeshInteriorDomain<M>>;

        FiniteElementConnectivity()
                :
                tag(),
                coordinates(Coordinates::Zero()),
                neighbours(Neighbours()),
                domain(DomainPointer())
        {}

        auto
        getCurrentCoordinates()
        {
            return coordinates;
        }

        auto
        getCurrentCoordinates()
        const
        {
            return coordinates;
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

        DomainPointer domain;

    };

    /*
     * BUILD
     */

    template<Element E, Domain D, auto M>
    struct FiniteElementGeometry : public FiniteElementConnectivity<E, D, M>
    {

        using ElementConnectivity = FiniteElementConnectivity<E, D, M>;

    private:

        using Self = element::FiniteElementGeometry<E, D, M>;

    public:

        FiniteElementGeometry()
        :
        ElementConnectivity()
        {}

        template<BasisDescription B>
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

        template<BasisDescription B>
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
                        auto new_value = std::abs(getDistanceInCurrentConfiguration(pt0, pt1, k));
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
            return Matrix<Vector<Real, E.dimPoint()> const>(d + offset);
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
            for (Indx q = 0; q < ElementQuadrature<seg_02, q_seg>::dim_quadrature; ++q) {
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
            static_assert(E.dim == D.dim - 1);
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
            static_assert(E.dim == D.dim - 2);
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = Vector<Real, D.dim>();
            for (Indx i = 0; i < D.dim; ++i) {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, direction_arg);
            }
            return tangent_vector;
        }

    };

    /*
     * Unknowns
     */

    template<Element E, Domain D, auto F, auto M, auto C>
    struct FiniteElementPolicy;

    template<Element E, Domain D, auto F, auto M, auto C>
    struct FiniteElementBase;

    template<Element E, Domain D, auto F, auto M, auto C>
    struct FiniteElementUnknown //public : FiniteElementModule<E, F, M, C, D>::Module
    {

        auto const static constexpr mixed_element_index = aggregate::index(C, M);

        auto const static constexpr finite_element_index = aggregate::index(M, F);

        static constexpr
        auto
        dimField()
        {
            return FField(F.ord_field, D.dim).size();
        }

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
            return numerics::max(-1, dimUnknowns() * dimField());
        }

        template<Indx I, Indx J>
        static constexpr
        auto
        dimComponentUnknowns()
        requires(!Point<E>)
        {
            return numerics::max(-1, ElementComponent<E, I, J, FiniteElementPolicy, D, F, M, C>::dimUnknowns());
        }

        template<Indx I, Indx J>
        static constexpr
        auto
        numComponentUnknowns()
        requires(!Point<E>)
        {
            return numerics::max(-1, dimComponentUnknowns<I, J>() * dimField());
        }

        FiniteElementUnknown()
        :
        unknowns(),
        bindings()
        {}

        FiniteElementUnknown(
                auto &
                unknown_index_arg
        )
        requires(numUnknowns() > -1)
        :
        unknowns(numUnknowns()),
        bindings(0)
        {
            unknown_index_arg += numUnknowns();
        }

        FiniteElementUnknown(
                auto &
                unknown_index_arg,
                auto &
                binding_index_arg,
                auto &&...
                direction_args
        )
        requires(numUnknowns() > -1)
        :
        unknowns(numUnknowns()),
        bindings(numUnknowns() * sizeof...(direction_args))
        {
            unknown_index_arg += numUnknowns();
            binding_index_arg += numUnknowns() * sizeof...(direction_args);
        }

        Vector<Real, numUnknowns()> unknowns;
        Vector<Real> bindings;

        struct Implementation : public FiniteElement<E, D, C>
        {

//            using Self = Implementation;
//
//            using Module = typename Self::template Type<mixed_element_index>::template Type<finite_element_index>;

//            using Implll = typename FiniteElementBase<E, D, F, M, C>::ImplementationBase;
//
//            void
//            initialize()
//            {
//                static_cast<Implll *>(this)->initialize();
//            }

            auto &
            fetch()
            {
                return this->template get<mixed_element_index>().template get<finite_element_index>();
            }

            auto const &
            fetch()
            const
            {
                return this->template get<mixed_element_index>().template get<finite_element_index>();
            }

            template<Indx I, Indx J>
            auto &
            fetch(
                    auto
                    index_arg
            )
            requires(!Point<E>)
            {
                return this->components.template get<I>().template get<J>().template get(index_arg).get().template get<mixed_element_index>().template get<finite_element_index>();
            }

            template<Indx I, Indx J>
            auto const &
            fetch(
                    auto
                    index_arg
            )
            const
            requires(!Point<E>)
            {
                return this->components.template get<I>().template get<J>().template get(index_arg).get().template get<mixed_element_index>().template get<finite_element_index>();
            }

        };

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Cell<E, D>)
    struct FiniteElementBase<E, D, F, M, C> : public FiniteElementUnknown<E, D, F, M, C> // va constituer l'elemeent
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
        requires(Cell<E, D>)
        {
            return FField::fromMapping(FField(F.ord_field, D.dim), App).size();
        }

        static constexpr
        auto
        rangeOperator(
                MappingOperator
                App
        )
        requires(Cell<E, D>)
        {
            auto value = Indx(0);
            for (auto m : F.mappings.data) {
                if (m == App) {
                    return Pair<Indx>(value, value + rowsOperator(App));
                }
                value += FField::fromMapping(FField(F.ord_field, D.dim), m).size();
            }
            return Pair<Indx>(Indx(0), Indx(0));
        }

        static constexpr
        auto
        rowsOperator()
        requires(Cell<E, D>)
        {
            auto value = Indx(0);
            for (int i = 0; i < F.mappings.size(); ++i) {
                value += FField::fromMapping(FField(F.ord_field, D.dim), F.mappings.get(i)).size();
            }
            return value;
        }

        static constexpr
        auto
        colsOperator()
        requires(Cell<E, D>)
        {
            return FiniteElementPolicy<E, D, F, M, C>::colsOperator();
        }

    private:

        template<MappingOperator App>
        using OperatorMatrix = Matrix<Real, rowsOperator(App), colsOperator()>;

    public:

        auto const static constexpr quadrature = Quadrature(C.quadrature_rule, C.ordIntegration(D));

        using Operators = Array<collection::ArrayWrapper<F.mappings, OperatorMatrix>, dimQuadrature<E, quadrature>()>;

        FiniteElementBase()
        {
            // print("I am here right now :", this->tag);
        }

        Operators operators;

        Module module;

        void
        initialize(
                auto &
                finite_element_arg
        )
        {
            module.initialize(finite_element_arg);
        }

        struct Implementation : public FiniteElementUnknown<E, D, F, M, C>::Implementation // herite de l'element
        {

            using ImplementationCast = typename FiniteElementPolicy<E, D, F, M, C>::Implementation;

            void
            initialize()
            {
                static_cast<ImplementationCast *>(this)->initialize();
            }

            void
            initialize(
                    auto &
                    finite_element_arg
            )
            {
                static_cast<ImplementationCast *>(this)->initialize(finite_element_arg);
            }

            static
            void
            init()
            {
                ImplementationCast::init();
            }

        };

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Face<E, D>)
    struct FiniteElementBase<E, D, F, M, C> : public FiniteElementUnknown<E, D, F, M, C>
    {

        struct Implementation : public FiniteElementUnknown<E, D, F, M, C>::Implementation // herite de l'element
        {

            using ImplementationCast = typename FiniteElementPolicy<E, D, F, M, C>::Implementation;

            void
            initialize()
            {
                static_cast<ImplementationCast *>(this)->initialize();
            }

            static
            void
            init()
            {
                ImplementationCast::init();
            }

        };

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Edge<E, D>)
    struct FiniteElementBase<E, D, F, M, C> : public FiniteElementUnknown<E, D, F, M, C>
    {

        struct Implementation : public FiniteElementUnknown<E, D, F, M, C>::Implementation // herite de l'element
        {

            using ImplementationCast = typename FiniteElementPolicy<E, D, F, M, C>::Implementation;

            void
            initialize()
            {
                static_cast<ImplementationCast *>(this)->initialize();
            }

            static
            void
            init()
            {
                ImplementationCast::init();
            }

        };

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Node<E, D>)
    struct FiniteElementBase<E, D, F, M, C> : public FiniteElementUnknown<E, D, F, M, C>
    {

        struct Implementation : public FiniteElementUnknown<E, D, F, M, C>::Implementation // herite de l'element
        {

            using ImplementationCast = typename FiniteElementPolicy<E, D, F, M, C>::Implementation;

            void
            initialize()
            {
                static_cast<ImplementationCast *>(this)->initialize();
            }

            static
            void
            init()
            {
                ImplementationCast::init();
            }

        };

    };

    template<Element E, Domain D, auto F, auto M, auto C>
    requires(Cell<E, D> && F.method == fem_hho && M.model == Model::Solid)
    struct FiniteElementPolicy<E, D, F, M, C>
    {

        using Mother = FiniteElementBase<E, D, F, M, C>;

        static constexpr
        auto
        dimUnknowns()
        {
            return Intg(EB<E, BasisDescription(Basis::Monomial, F.discretization.ord_cell)>::dim_basis);
        }

        static constexpr
        auto
        colsOperator()
        {

            auto in = Intg(0);
            auto set_num_components = [& in] <auto K> () constexpr mutable {
                auto nj = FiniteElementBase<E, D, F, M, C>::template numComponentUnknowns<0, K>();
                auto num_c = numComponents<E, 0, K>();
                in += num_c * nj;
            };
            collection::apply<numComponents<E, 0>()>(set_num_components);
            return FiniteElementBase<E, D, F, M, C>::numUnknowns() + in;
        }

        struct Module
        {

        private:

            using Base = FiniteElementBase<E, D, F, M, C>;

        public:

            Indx value;

            Matrix<Real, Base::colsOperator(), Base::colsOperator()> stabilization;

            void
            setStabilization(
                    auto &
                    finite_element_arg
            )
            {
                stabilization.setOnes();
                stabilization(0, 0) = finite_element_arg.tag;
            }

            void
            initialize(
                    auto &
                    finite_element_arg
            )
            {
                setStabilization(finite_element_arg);
            }

        };

        struct Implementation : public FiniteElementBase<E, D, F, M, C>::Implementation
        {

            void
            setStabilization(
                    auto &
                    finite_element_arg
            )
            {
                auto & mod = finite_element_arg.template get<0>().template get<0>();
                mod.module.stabilization.setOnes();
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
                print("---- SETTING Gradient");
            }

            template<MappingOperator Map, auto I>
            void
            setOperator(
                    auto &
                    finite_element_arg
            )
            requires(Map == MappingOperator::Identity)
            {
                print("---- SETTING Identity");
            }

            void
            setOperators(
                    auto &
                    finite_element_arg
            )
            {
                auto set_num_components = [&] <auto I> () constexpr mutable {
                    setOperator<F.mappings.get(I), I>();
                };
                collection::apply<F.mappings.size()>(set_num_components);
            }

            Module &
            getModule()
            {
                return this->fetch().module;;
            }

            void
            initialize(
                    auto &
                    finite_element_arg
            )
            {
                setStabilization(finite_element_arg);
            }

            void
            initialize()
            {
                print("hey tag :", this->tag);
                print(this->getCurrentCoordinates());
                this->fetch().module.value = 31;
                print("modulevalue :", this->fetch().module.value);
//                this->fetch().module.stabilization = Matrix<Real, colsOperator(), colsOperator()>::Identity();
                this->fetch().module.stabilization.setOnes();
//                setOperators();
                getModule().stabilization.setOnes();
                print("stabilization print start");
                print(this->fetch().module.stabilization);
                print("stabilization print stop");
            }

            static
            void
            init()
            {
                print("hey static");
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
            return Intg(EB<E, BasisDescription(Basis::Monomial, F.discretization.ord_face)>::dim_basis);
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementBase<E, D, F, M, C>::Implementation
        {

            void
            initialize()
            {
                print("hey tag :", this->tag());
            }

            static
            void
            init()
            {
                print("hey static");
            }

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

        struct Module// : public FiniteElementUnknown<E, F, M, C, D>
        {

        };

        struct Implementation : public FiniteElementBase<E, D, F, M, C>::Implementation
        {

            void
            initialize()
            {
                print("hey tag :", this->tag());
            }

            static
            void
            init()
            {
                print("hey static");
            }

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

        struct Module// : public FiniteElementUnknown<E, F, M, C, D>
        {

        };

        struct Implementation : public FiniteElementBase<E, D, F, M, C>::Implementation
        {

            void
            initialize()
            {
                print("hey tag :", this->tag());
            }

            static
            void
            init()
            {
                print("hey static");
            }

        };

    };


    template<Element E, Domain D, auto M, auto C>
    struct MixedElementPolicy2 : public decltype(M)::template Elements2<FiniteElementBase, E, D, M, C>
    {

        using Base = typename decltype(M)::template Elements2<FiniteElementBase, E, D, M, C>;

        template<Indx I>
        using ImplementationCast = typename Base::template Type<I>::Implementation;

        template<Indx I>
        struct Iml : public Base::template Type<I>::Implementation
        {



        };

        void
        initialize()
        {
            this->template get<0>().initialize();
//            static_cast<Implementation<0> *>(this)->initialize();
        }

        MixedElementPolicy2()
        {}

        struct Implementation
        {

            void
            initialize()
            {
//                this->template get<0>().initialize();
                auto set_num_components = [&] <auto K> (auto & x) constexpr mutable {
                    x.fun();
                };
//            static_cast<Implementation<0> *>(this)->initialize();
            }



        };

    };

    template<Element E, Domain D, auto C>
    struct CoupledElementPolicy2 : public decltype(C)::template Elements2<MixedElementPolicy2, E, D, C>
    {

        using Base = typename decltype(C)::template Elements2<MixedElementPolicy2, E, D, C>;

        template<Indx I>
        using Implementation = typename Base::template Type<I>::ImplementationBase;

        void
        initialize()
        {
            this->template get<0>().initialize();
            auto set_num_components = [&] <auto K> () constexpr mutable {
                static_cast<Base *>(this)->implementation();
            };
//            static_cast<Implementation<0> *>(this)->initialize();
        }



    };

    template<Element E, Domain D, auto C>
    struct ComponentTest
    {

        ComponentTest(
                auto const &
                ref
        )
        :
//        _ref(UniquePointer<FiniteElementGeometry<E, D, C>>(ref)),
        _ref2(reinterpret_cast<FiniteElementGeometry<E, D, C> const * const>(ref)),
        _tag(ref->getCurrentCentroid()(0))
        {
//            sayCurrentCentroid();
        }

        void
        init()
        {
            _tag = _ref2->getCurrentCentroid()(0);
        }

        void
        sayCurrentCentroid()
        const
        {
            print(_ref2->getCurrentCentroid());
        }

        Indx _tag;

//        UniquePointer<FiniteElementGeometry<E, D, C>> _ref;

        FiniteElementGeometry<E, D, C> const * const _ref2;

    };

    template<Element E, Domain D, auto C>
    struct FiniteElement :
//            public FiniteElementConnectivity<E, D, C>,
            public FiniteElementGeometry<E, D, C>,
//                    public CoupledElementPolicy<E, D, C>::Type
            public CoupledElementPolicy2<E, D, C>
//            public ComponentTest<E, D, C>
//            public SubElement<E, D, C, C>
    {

        FiniteElement()
        :
        FiniteElementGeometry<E, D, C>(),
        CoupledElementPolicy2<E, D, C>()
//        ComponentTest<E, D, C>(static_cast<FiniteElementGeometry<E, D, C>>(* this))
//        ComponentTest<E, D, C>(this)
        {}

        void
        build()
        {

        }

        auto
        toGeom()
        {
            return UniquePointer<FiniteElementGeometry<E, D, C>>::make(static_cast<FiniteElementGeometry<E, D, C>>(* this));
//            return SharedPointer<FiniteElementGeometry<E, D, C>>::make(* this);
        }

        void
        initialize()
        {
            static_cast<typename CoupledElementPolicy2<E, D, C>::template Type<0>::template Type<0>::Implementation *>(this)->initialize();
        }

        void
        initialize2()
        {
            this->template get<0>().template get<0>().initialize(* this);
        }

//        FiniteElementGeometry<E, D, C> const &
//        toGeometry()
//        {
//            return static_cast<FiniteElementGeometry<E, D, C> *>(* this);
//        }

//        static
//        void
//        init()
//        {
//            CoupledElementPolicy2<E, D, C>::template Type<0>::template Type<0>::Implementation::init();
//        }

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

    struct CC;

    struct AA
    {

        AA() : mem_(0) {}

        AA(Indx arg) : mem_(arg) {}

        Indx mem_;

    };

    struct BB
    {

        BB() {}

        BB(auto & arg) : ref_(arg) {}

        void
        set(auto const & a)
        {
            ref_ = & a;
        }

        Indx memb_;

        CC * ref_;

    };

    struct CC : public AA, public Collection<BB, BB>
    {

        using Col = Collection<BB, BB>;

        CC() : AA() {}

    };






    /*
     * Operators
     */

//    template<Element E, auto F, auto M, auto C, auto D>
//    struct FiniteElementOperators;
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    struct FiniteElementOperatorsPolicy
//    {
//
//        auto const static constexpr index = aggregate::index(M, F);
//
////        using DegreesOfFreedom = DegOfFre<FField(F.ord_field, D.dim)>;
////        DegreesOfFreedom degrees_of_freedom;
//
//        struct ImplementationBase : public FiniteElement<E, C, D>
//        {
//
//            using Self = ImplementationBase;
//
//            using FiniteElementOperators = typename Self::template Type<index>;
//
//            auto
//            fetch()
//            {
//                return this->template get<index>();
//            }
//
//            auto
//            fetch()
//            const
//            {
//                return this->template get<index>();
//            }
//
//            template<Indx I, Indx J>
//            auto
//            fetch(
//                    auto
//                    index_arg
//            )
//            {
//                return this->components.template get<I>().template get<J>().template get(index_arg).get().template get<index>();
//            }
//
//            template<Indx I, Indx J>
//            auto
//            fetch(
//                    auto
//                    index_arg
//            )
//            const
//            {
//                return this->components.template get<I>().template get<J>().template get(index_arg).get().template get<index>();
//            }
//
//        };
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    struct CellModuleBase : public FiniteElementModuleBase<E, F, M, C, D>
//    {
//
//        using Unknowns = FiniteElementModule<E, F, M, C, D>;
//
//        using Self = CellModuleBase;
//
//        using Base = FiniteElementModuleBase<E, F, M, C, D>;
//
//        static constexpr
//        auto
//        rowsOperator(
//                MappingOperator
//                App
//        )
//        {
//            return FField::fromMapping(FField(F.ord_field, D.dim), App).size();
//        }
//
//        static constexpr
//        auto
//        rangeOperator(
//                MappingOperator
//                App
//        )
//        {
//            auto value = Indx(0);
//            for (auto m : F.mappings.data) {
//                if (m == App) {
//                    return Pair<Indx>(value, value + rowsOperator(App));
//                }
//                value += FField::fromMapping(FField(F.ord_field, D.dim), m).size();
//            }
//            return Pair<Indx>(Indx(0), Indx(0));
//        }
//
//        static constexpr
//        auto
//        rowsOperator()
//        {
//            auto value = Indx(0);
//            for (int i = 0; i < F.mappings.size(); ++i) {
//                value += FField::fromMapping(FField(F.ord_field, D.dim), F.mappings.get(i)).size();
//            }
//            return value;
//        }
//
//        using Operator = Matrix<Real, Self::rowsOperator(), Unknowns::colsOperator()>;
//
//        using Operators = Array<Operator, C.ordIntegration(D)>;
//
//        using IntegrationTags = Array<Indx, C.ordIntegration(D)>;
//
//        CellModuleBase()
//        :
//        Base(),
//        operators()
//        {}
//
//        template<MappingOperator App>
//        auto
//        getOperator(
//                Indx
//                num_op
//        )
//        const
//        {
//            auto const constexpr rng = rangeOperator(App);
//            auto const constexpr sze = rng.j - rng.i;
//            if constexpr(Unknowns::numUnknowns() != -1) {
//                return operators.get(num_op).template block<sze, Unknowns::numUnknowns()>(rng.i, 0);
//            }
//            else {
//                return operators.get(num_op).block(rng.i, 0, sze, operators.get(num_op).cols());
//            }
//        }
//
//        template<MappingOperator App>
//        auto
//        getOperator(
//                Indx
//                num_op
//        )
//        {
//            auto const constexpr rng = rangeOperator(App);
//            auto const constexpr sze = rng.j - rng.i;
//            if constexpr(Unknowns::dimUnknowns() != -1) {
//                return operators.get(num_op).template block<sze, Unknowns::dimUnknowns()>(rng.i, 0);
//            }
//            else {
//                return operators.get(num_op).block(rng.i, 0, sze, operators.get(num_op).cols());
//            }
//        }
//
//        Operators operators;
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    concept SolidHHOCell = D.dim - E.dim == 0 && F.method == fem_hho && M.model == Model::Solid;
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(D.dim - E.dim == 0 && F.method == fem_hho && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D> : public CellModuleBase<E, F, M, C, D>
//    {
//
//        using Base = FiniteElementOperatorsPolicy<E, F, M, C, D>;
//
//        using Self = FiniteElementOperators<E, F, M, C, D>;
//
//        using UUU = FiniteElementModule<E, F, M, C, D>;
//
//        using Stabilization = Matrix<Real, UUU::dimUnknowns(), UUU::dimUnknowns()>;
//
//        Stabilization stabilization;
//
////        FiniteElementOperators() = default;
////
////        FiniteElementOperators(
////                Indx
////                u_index_arg
////        )
////        {}
//
////        template<MappingOperator App>
////        struct OperatorFactory;
////
////        template<MappingOperator App>
////        requires(App == MappingOperator::Gradient)
////        struct OperatorFactory<App>
////        {
////
////            auto const static constexpr rows = Self::dimOperator(App);
////            auto const static constexpr cols = UUU::numUnknowns();
////
////            using Operator = Matrix<Real, rows, cols>;
////
////            static
////            auto
////            makeOp()
////            {
////                return Matrix<Real, rows, cols>().setOnes();
////            };
////
////        };
////
////        template<MappingOperator App>
////        using Operatorr = typename OperatorFactory<App>::Operator;
//
////        using Operatorf = typename collection::ArrayCollectionWrapper<F.mappings>::template Wrapper<Operatorr>;
//
////        using Grad = OperatorFactory<MappingOperator::Gradient>;
//
////        Matrix<Real, Grad::rows, Grad::cols> rss2;
//
////        DegreesOfFreedom degrees_of_freedom;
//
////        Stabilization stabilization;
//
//        struct Implementation : public FiniteElementOperatorsPolicy<E, F, M, C, D>::ImplementationBase
//        {
//
//            template<MappingOperator App>
//            struct OperatorFactory;
//
//            template<MappingOperator App>
//            requires(App == MappingOperator::Gradient)
//            struct OperatorFactory<App>
//            {
//
//                auto const static constexpr rows = Self::dimOperator(App);
//                auto const static constexpr cols = Self::numUnknowns();
//
//                static
//                auto
//                makeOp()
//                {
//                    return Matrix<Real, rows, cols>().setOnes();
//                };
//
//            };
//
//        };
//
//        static
//        void
//        heyCell()
//        {
//            print("hey cell");
//        }
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(D.dim - E.dim == 1 && F.method == fem_hho && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D>
//    {
//
//        using Pol = FiniteElementOperatorsPolicy<E, F, M, C, D>;
//
//        static constexpr
//        auto
//        numUnknowns()
//        {
//            return EB<E, BasisDescription(Basis::Monomial, F.discretization.ord_face)>::dim_basis;
//        }
//
//        using DegreesOfFreedom = DegOfFre<FField(F.ord_field, D.dim)>;
//
//        static
//        void
//        heyFace()
//        {
//            print("hey face");
//        }
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(D.dim - E.dim == 2 && F.method == fem_hho && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D>
//    {
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(D.dim - E.dim == 3 && F.method == fem_hho && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D>
//    {
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(E != pnt_00 && D.dim - E.dim == 0 && F.method == fem_lag && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D>
//    {
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(E != pnt_00 && D.dim - E.dim == 1 && F.method == fem_lag && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D>
//    {
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(E != pnt_00 && D.dim - E.dim == 2 && F.method == fem_lag && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D>
//    {
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(E != pnt_00 && D.dim - E.dim == 3 && F.method == fem_lag && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D>
//    {
//
//    };
//
//    template<Element E, auto F, auto M, auto C, auto D>
//    requires(E == pnt_00 && F.method == fem_lag && M.model == Model::Solid)
//    struct FiniteElementOperators<E, F, M, C, D>
//    {
//
//    };


//    template<auto F, auto M>
//    constexpr inline
//    auto
//    fieldDescription()
//    {
//        return FField{F.ord_field, M.dim_euclidean};
//    }
//
//    template<auto F, auto M>
//    constexpr inline
//    auto
//    fieldDescription(
//            MappingOperator
//            mapping_arg
//    )
//    {
//        return FField::fromMapping(fieldDescription<F, M>(), mapping_arg);
//    }

//    constexpr inline
//    auto
//    fieldSize(
//            FField
//            field_description_arg
//    )
//    {
//        return field_description_arg.size();
//    }

//    template<Element E, auto M>
//    constexpr inline
//    auto
//    modelDescription()
//    {
//        return ModelDescription(E, M.dim_euclidean, M.model);
//    }

//    template<auto F>
//    constexpr inline
//    auto
//    finiteElementMethod()
//    {
//        return F.fem;
//    }





















//    namespace detail
//    {
//
//        template<Element E, Indx, FiniteElementMethod D, auto F, auto M>
//        struct FiniteElementE;
//
//        template<Element E, FiniteElementMethod D, auto F, auto M>
//        struct FiniteElementEPolicy
//        {
//
//        private:
//
//            using Self = FiniteElementE<E, E.dim - M.dim_euclidean, D, F, M>;
//
//            template<Element Ec>
//            using ComponentT = FiniteElementE<Ec, Ec.dim - M.dim_euclidean, D, F, M>;
//
//        public:
//
//            using Components = typename element::ElementGeometry<E>::template Components<ComponentT>;
//
//            template<Indx I, Indx J>
//            using Component = typename Components::template Type<I>::template Type<J>::Type;
//
//            struct ImplementationBase : public FiniteElement<E, M>
//            {
//
//                using Self = ImplementationBase;
//
//                using FiniteElement = typename Self::template Type<M.index(F)>;
//
//                auto
//                fetch()
//                {
//                    return this->template get<M.index(F)>();
//                }
//
//                auto
//                fetch()
//                const
//                {
//                    return this->template get<M.index(F)>();
//                }
//
//                template<Indx I, Indx J>
//                auto
//                fetch(
//                        auto
//                        index_arg
//                )
//                {
//                    auto & p = this->components.template get<I>().template get<J>().template get(index_arg).get();
//                    return p.template get<M.index(F)>();
//                }
//
//                template<Indx I, Indx J>
//                auto
//                fetch(
//                        auto
//                        index_arg
//                )
//                const
//                {
//                    auto const & p = this->components.template get<I>().template get<J>().template get(index_arg).get();
//                    return p.template get<M.index(F)>();
//                }
//
//            };
//
//        };
//
//        /*
//         * CELL
//         */
//
//        template<Element E, auto F, auto M>
//        struct FiniteElementE<E, 0, fem_hho, F, M>
//        {
//
//        private:
//
//            using Self = FiniteElementE;
//
//            auto const static constexpr basis_description = BasisDescription(Basis::Monomial, F.ord_cell);
//
//        public:
//
//            using DegreesOfFreedom = DegreesOfFreedom<E, FField{F.ord_field, M.dim_euclidean}, BasisDescription(Basis::Monomial, F.ord_cell)>;
//
//        private:
//
//            using FiniteElementPolicy = detail::FiniteElementEPolicy<E, fem_hho, F, M>;
//
//            template<Indx I>
//            static constexpr
//            auto
//            numUnknowns()
//            {
//                return FiniteElementPolicy::template Component<0, I>::DegreesOfFreedom::num_unknowns;
//            }
//
//            template<Indx I = 0>
//            static constexpr
//            auto
//            numUnknowns(
//                    auto &
//                    value
//            )
//            {
//                value += numUnknowns<I>() * numComponents<E, 0, I>();
//                if constexpr (I < numComponents<E, 0>() - 1) {
//                    numUnknowns<I + 1>(value);
//                }
//            }
//
//            static constexpr
//            auto
//            numUnknowns()
//            {
//                auto value = DegreesOfFreedom::num_unknowns;
//                numUnknowns(value);
//                return value;
//            }
//
//            auto const static constexpr num_unknowns = numUnknowns();
//
//        public:
//
//            using Stabilization = Matrix<Real, num_unknowns, num_unknowns>;
//
////            template<Mapping O>
////            using MappingMatrix = Matrix<Real, mappingMatrixSize<O, F, M>(), num_unknowns>;
//
//            template<MappingOperator O>
//            using MappingMatrix = Matrix<Real, fieldSize(fieldDescription<F, M>(O)), num_unknowns>;
//
//            using Mappings = typename ArrayCollectionWrapper<FiniteElementPolicy::mappings>::template Wrapper<MappingMatrix>;
//
//            FiniteElement()
//            :
//            degrees_of_freedom(),
//            stabilization(Stabilization::Zero())
//            {}
//
//            FiniteElement(
//                    DegreesOfFreedom const &
//                    degrees_of_freedom_arg,
//                    Stabilization const &
//                    stabilization_arg
//            )
//            :
//            degrees_of_freedom(degrees_of_freedom_arg),
//            stabilization(stabilization_arg)
//            {}
//
//            FiniteElement(
//                    DegreesOfFreedom &&
//                    degrees_of_freedom_arg,
//                    Stabilization &&
//                    stabilization_arg
//            )
//            :
//            degrees_of_freedom(degrees_of_freedom_arg),
//            stabilization(stabilization_arg)
//            {}
//
//            DegreesOfFreedom degrees_of_freedom;
//
//            Stabilization stabilization;
//
//            struct Implementation : public FiniteElementPolicy::ImplementationBase
//            {
//
//                using Self = Implementation;
//
//                template<Mapping O>
//                auto
//                getMappingMatrix(
////                        auto const &
////                        point_arg
//                )
//                const
//                {
//                    return MappingMatrix<O>().setZero();
//                }
//
//                auto
//                setStabilization()
//                const
//                {
//                    return Stabilization::Identity();
//                }
//
//                template<Indx I>
//                auto
//                getUnknowns(
//                        auto &
//                        v,
//                        auto &
//                        in
//                )
//                const
//                {
//                    for (int i = 0; i < numComponents<E, 0>(); ++i) {
//                        auto const & n = this->template fetch<0, I>(i).degrees_of_freedom.unknowns.values;
//                        v.template segment<numUnknowns<I>()>(in) = n;
//                        in += numUnknowns<I>();
//                    }
//                    if constexpr (I < numComponents<E, 0>() - 1) {
//                        getUnknowns<I + 1>(v, in);
//                    }
//                }
//
//                auto
//                getUnknowns()
//                const
//                {
//                    auto v = Matrix<Real, num_unknowns>();
//                    auto in = Indx(0);
//                    auto const & n = this->template fetch().degrees_of_freedom.unknowns.values;
//                    v.template segment<DegreesOfFreedom::num_unknowns>(in) = n;
//                    in += DegreesOfFreedom::num_unknowns;
//                    getUnknowns(v, in);
//                    return v;
//                }
//
//                auto
//                make()
//                const
//                {
//                    auto unknown_index = Indx(0);
//                    auto binding_index = Indx(0);
////                    auto u = DegreesOfFreedom(unknown_index, binding_index, {{1, 2}});
//                    auto u = DegreesOfFreedom(unknown_index);
//                    auto s = setStabilization();
//                    return FiniteElement(u, s);
//                }
//
//                void
//                hey()
//                const
//                {
//                    print("hey a :", num_unknowns);
//                }
//
//            };
//
//        };
//
//        /*
//         * FACE
//         */
//
//        template<ElementDescription E, auto F, auto M>
//        struct FiniteElement<E, face_solid, FiniteElementMethod::HybridHighOrder, F, M>
//        {
//
//        private:
//
//            auto const static constexpr basis_description = BasisDescription(Basis::Monomial, F.face_order);
//
//        public:
//
//            using DegreesOfFreedom = finite_element::DegreesOfFreedom<E, fieldDescription<F, M>(), basis_description>;
//
//            DegreesOfFreedom degrees_of_freedom;
//
//            struct Implementation : public Element<E, M>
//            {
//
//                using Self = Implementation;
//
//                static
//                FiniteElement
//                make()
//                {
//                    return FiniteElement();
//                }
//
//            };
//
//        };
//
//        /*
//         * EDGE
//         */
//
//        template<ElementDescription E, auto F, auto M>
//        struct FiniteElement<E, edge_solid, FiniteElementMethod::HybridHighOrder, F, M>
//        {
//
//            using DegreesOfFreedom = Void;
//
//            struct Implementation : public Element<E, M>
//            {
//
//                using Self = Implementation;
//
//                static
//                FiniteElement
//                make()
//                {
//                    return FiniteElement();
//                }
//
//            };
//
//        };
//
//        /*
//         * NODE
//         */
//
//        template<ElementDescription E, auto F, auto M>
//        struct FiniteElement<E, node_solid, FiniteElementMethod::HybridHighOrder, F, M>
//        {
//
//            using DegreesOfFreedom = Void;
//
//            struct Implementation : public Element<E, M>
//            {
//
//                using Self = Implementation;
//
//                static
//                FiniteElement
//                make()
//                {
//                    return FiniteElement();
//                }
//
//            };
//
//        };
//
//    }
//
//    template<ElementDescription E, auto F, auto M>
//    using FiniteElement = detail::FiniteElement<E, modelDescription<E, M>(), finiteElementMethod<F>(), F, M>;
//
//}
//
//    namespace detail {
//
//        template<Element E, auto M>
//        struct MixedElementPolicy {
//
//        private:
//
//            template<auto F>
//            using FiniteElement = FiniteElement<E, F, M>;
//
//        public:
//
//            using Type = typename decltype(M)::template Elements<FiniteElement>;
//
//            struct MixedElementBase : public Type {
//
//            private:
//
//                using Base = Type;
//
//                using Self = MixedElementBase;
//
//                template<Indx I>
//                using Implementation = typename Self::template Type<I>::Implementation;
//
//            public:
//
//                MixedElementBase()
//                        :
//                        Base(make()) {}
//
//            private:
//
//                template<Indx I = 0>
//                auto
//                make(
//                        auto &
//                        b
//                )
//                const {
//                    auto e = static_cast<Implementation<I> const *>(this)->make();
//                    b.template get<I>() = e;
//                    if constexpr (I < Self::size() - 1) {
//                        make<I + 1>(b);
//                    }
//                }
//
//                auto
//                make()
//                const {
//                    Base b;
//                    make(b);
//                    return b;
//                }
//
//            };
//
//        };
//
//        template<Element E, Model Md, auto M>
//        struct MixedElement2;
//
//        template<Element E, auto M>
//        struct MixedElement2<E, cell_solid, M> : public MixedElementPolicy<E, M>::MixedElementBase
//        {
//
//        private:
//
//            using Base = typename MixedElementPolicy<E, M>::MixedElementBase;
//
//            using Self = MixedElement2;
//
//            template<Indx I>
//            using Implementation = typename Self::template Type<I>::Implementation;
//
//            template<Indx I>
//            using Mappings = typename Self::template Type<I>::Mappings;
//
//            using MappingMatrices = typename ArrayCollectionWrapper<range<Self::size()>(0)>::template Wrapper<Mappings>;
//
//        public:
//
//            using Base::Base;
//
//            MappingMatrices mapss;
//
//            template<Indx I, Mapping O>
//            auto
//            getMappingMatrix()
//            const
//            {
//                return static_cast<Implementation<I> const *>(this)->template getMappingMatrix<O>();
//            }
//
////            template<Indx I>
////            auto
////            getUnknowns()
////            const
////            {
////                return static_cast<Implementation<I> const *>(this)->getUnknowns();
////            }
//
////            template<Indx I, Mapping O>
////            auto
////            getMappedField()
////            const
////            {
////                print("getMappingMatrix<I, O>() :", getMappingMatrix<I, O>());
////                print("getUnknowns<I>() :", getUnknowns<I>());
////                return getMappingMatrix<I, O>() * getUnknowns<I>();
////            }
//
//            template<Indx I>
//            void
//            hey()
//            const
//            {
//                static_cast<Implementation<I> const *>(this)->hey();
//            }
//
//        };
//
//        template<ElementDescription E, auto M>
//        struct MixedElement2<E, face_solid, M> : public MixedElementPolicy<E, M>::MixedElementBase
//        {
//
//        private:
//
//            using Base = typename MixedElementPolicy<E, M>::MixedElementBase;
//
//            using Self = MixedElement2;
//
//            template<Indx I>
//            using Implementation = typename Self::template Type<I>::Implementation;
//
//        public:
//
//            using Base::Base;
//
//        };
//
//        template<ElementDescription E, auto M>
//        struct MixedElement2<E, edge_solid, M> : public MixedElementPolicy<E, M>::MixedElementBase
//        {
//
//        private:
//
//            using Base = typename MixedElementPolicy<E, M>::MixedElementBase;
//
//            using Self = MixedElement2;
//
//            template<Indx I>
//            using Implementation = typename Self::template Type<I>::Implementation;
//
//        public:
//
//            using Base::Base;
//
//        };
//
//        template<ElementDescription E, auto M>
//        struct MixedElement2<E, node_solid, M> : public MixedElementPolicy<E, M>::MixedElementBase
//        {
//
//        private:
//
//            using Base = typename MixedElementPolicy<E, M>::MixedElementBase;
//
//            using Self = MixedElement2;
//
//            template<Indx I>
//            using Implementation = typename Self::template Type<I>::Implementation;
//
//        public:
//
//            using Base::Base;
//
//        };
//
//    }
//
//    template<ElementDescription E, auto M>
//    struct MixedElement : public detail::MixedElement2<E, finite_element::modelDescription<E, M>(), M>
//    {
//
//    private:
//
//        using Base = detail::MixedElement2<E, finite_element::modelDescription<E, M>(), M>;
//
//    public:
//
//        using Base::Base;
//
//    };








}

#endif //LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
