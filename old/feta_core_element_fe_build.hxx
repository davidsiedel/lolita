//
// Created by dsiedel on 01/04/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_FE_BUILD_HXX
#define FETA_FETA_CORE_ELEMENT_FE_BUILD_HXX

#include "new/feta_core_element_element_connectivity.hxx"
#include "new/feta_core_element_quadrature_gauss_segment.hxx"
#include "new/feta_core_geometry_geometry.hxx"
#include "new/feta_core_element_basis.hxx"
#include "new/feta_core_element_fe_data.hxx"

namespace feta::core::element
{

    namespace detail
    {

        template<ElementDescription E, auto M>
        struct MixedElementPolicy
        {

        private:

            template<auto F>
            using FiniteElement = finite_element::FiniteElement<E, F, M>;

        public:

            using Type = typename TypeTraits<decltype(M)>::RawType::template Elements<FiniteElement>;

            struct MixedElementBase : public TypeTraits<decltype(M)>::RawType::template Elements<FiniteElement>
            {

            private:

                using Base = typename TypeTraits<decltype(M)>::RawType::template Elements<FiniteElement>;

                using Self = MixedElementBase;

                template<Indx I>
                using Implementation = typename Self::template Type<I>::Implementation;

            public:

                MixedElementBase()
                        :
                        Base(make())
                {}

            private:

                template<Indx I = 0>
                auto
                make(
                        auto &
                        b
                )
                const
                {
                    auto e = static_cast<Implementation<I> const *>(this)->make();
                    b.template get<I>() = e;
                    if constexpr (I < Self::size() - 1) {
                        make<I + 1>(b);
                    }
                }

                auto
                make()
                const
                {
                    Base b;
                    make(b);
                    return b;
                }

            };

        };

        template<ElementDescription E, ModelDescription Md, auto M>
        struct MixedElement2;

        template<ElementDescription E, auto M>
        struct MixedElement2<E, cell_solid, M> : public MixedElementPolicy<E, M>::MixedElementBase
        {

        private:

            using Base = typename MixedElementPolicy<E, M>::MixedElementBase;

            using Self = MixedElement2;

            template<Indx I>
            using Implementation = typename Self::template Type<I>::Implementation;

            template<Indx I>
            using Mappings = typename Self::template Type<I>::Mappings;

            using MappingMatrices = typename ArrayCollectionWrapper<range<Self::size()>(0)>::template Wrapper<Mappings>;

        public:

            using Base::Base;

            MappingMatrices mapss;

            template<Indx I, Mapping O>
            auto
            getMappingMatrix()
            const
            {
                return static_cast<Implementation<I> const *>(this)->template getMappingMatrix<O>();
            }

//            template<Indx I>
//            auto
//            getUnknowns()
//            const
//            {
//                return static_cast<Implementation<I> const *>(this)->getUnknowns();
//            }

//            template<Indx I, Mapping O>
//            auto
//            getMappedField()
//            const
//            {
//                print("getMappingMatrix<I, O>() :", getMappingMatrix<I, O>());
//                print("getUnknowns<I>() :", getUnknowns<I>());
//                return getMappingMatrix<I, O>() * getUnknowns<I>();
//            }

            template<Indx I>
            void
            hey()
            const
            {
                static_cast<Implementation<I> const *>(this)->hey();
            }

        };

        template<ElementDescription E, auto M>
        struct MixedElement2<E, face_solid, M> : public MixedElementPolicy<E, M>::MixedElementBase
        {

        private:

            using Base = typename MixedElementPolicy<E, M>::MixedElementBase;

            using Self = MixedElement2;

            template<Indx I>
            using Implementation = typename Self::template Type<I>::Implementation;

        public:

            using Base::Base;

        };

        template<ElementDescription E, auto M>
        struct MixedElement2<E, edge_solid, M> : public MixedElementPolicy<E, M>::MixedElementBase
        {

        private:

            using Base = typename MixedElementPolicy<E, M>::MixedElementBase;

            using Self = MixedElement2;

            template<Indx I>
            using Implementation = typename Self::template Type<I>::Implementation;

        public:

            using Base::Base;

        };

        template<ElementDescription E, auto M>
        struct MixedElement2<E, node_solid, M> : public MixedElementPolicy<E, M>::MixedElementBase
        {

        private:

            using Base = typename MixedElementPolicy<E, M>::MixedElementBase;

            using Self = MixedElement2;

            template<Indx I>
            using Implementation = typename Self::template Type<I>::Implementation;

        public:

            using Base::Base;

        };

    }

    template<ElementDescription E, auto M>
    struct MixedElement : public detail::MixedElement2<E, finite_element::modelDescription<E, M>(), M>
    {

    private:

        using Base = detail::MixedElement2<E, finite_element::modelDescription<E, M>(), M>;

    public:

        using Base::Base;

    };

    template<ElementDescription E, auto M>
    struct ElementBuild : public ElementConnectivity<E, M>, public MixedElement<E, M>
    {

        using ElementConnectivity = element::ElementConnectivity<E, M>;

    private:

        using Self = ElementBuild;

//    protected:
//
//        template<auto F>
//        using ElementUnknowns = element::ElementUnknowns<E, F, M>;
//
//        template<auto F>
//        using ElementImplementation = element::ElementImplementation<E, F, M>;

    public:

//        using Unknowns = typename std::remove_cvref_t<decltype(M)>::template Elements<ElementUnknowns>;
//
//        using Implementation = typename std::remove_cvref_t<decltype(M)>::template Elements<ElementImplementation>;
//
//        Unknowns unknowns;

        ElementBuild()
        :
        ElementConnectivity()
        {}

//        template<auto F>
//        static constexpr
//        auto
//        getNumUnknowns()
//        {
//            if constexpr (TypeTraits<decltype(F)>::template is_raw_in<Indx, Intg>) {
//                return Implementation::template Type<F>::getNumUnknowns();
//            } else if constexpr (M.has(F)) {
//                return ElementImplementation<F>::getNumUnknowns();
//            } else {
//
//            }
//        }

        template<BasisDescription B>
        static constexpr
        auto
        getDimBasis()
        {
            using FiniteElementBasis = typename finite_element::FiniteElementBasis<E, B>;
            return FiniteElementBasis::dim_basis;
//            return finite_element::FiniteElementBasis<E, B>::dim_basis;
        }

        template<BasisDescription B>
        auto
        getBasisEvaluation(
                auto const &
                point_arg
        )
        const
        {
            using FiniteElementBasis = typename finite_element::FiniteElementBasis<E, B>::template Implementation<M>;
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
            using FiniteElementBasis = typename finite_element::FiniteElementBasis<E, B>::template Implementation<M>;
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
            auto ru = Matrix<Real, 3, Self::dim_shape>().setZero();
            auto du = Real(0);
            for (Indx i = 0; i < Self::dim_euclidean; ++i) {
                for (Indx j = 0; j < Self::dim_shape; ++j) {
                    ru(i, j) = Self::getShapeMappingDerivative(nds.row(i), point_arg, j);
                }
            }
            if constexpr (Self::dim_shape == 1) {
                du = std::abs(ru.col(0).norm());
            }
            else if constexpr (Self::dim_shape == 2) {
                du = std::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else if constexpr (Self::dim_shape == 3) {
                du = std::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
//            if constexpr (Ft == FrameType::AxiSymmetric) {
//                Real r0 = getShapeMappingEvaluation(nds.row(0), point);
//                if (r0 < 1.e-10) {
//                    r0 = 1.e-10;
//                }
//                du *= 2.0 * PI_ * r0;
//            }
            return du;
        }

        auto
        getDistanceInCurrentConfiguration(
                auto const &
                first_point_arg,
                auto const &
                second_point_arg,
                Intg
                direction_arg = -1
        )
        const
        {
            assert(-1 <= direction_arg <= Self::ord_shape);
            auto distance = Real();
            if constexpr (Self::dim_euclidean == Self::ord_shape) {
                distance = getCartesianDistance(first_point_arg, second_point_arg, direction_arg);
            } else {
                distance = getCurvedDistance(first_point_arg, second_point_arg, direction_arg);
            }
            return distance;
        }

        auto
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = Self::getReferenceCoordinates();
            auto current_diameters = Matrix<Real, Self::ord_shape>().setZero();
            for (Indx i = 0; i < Self::num_nodes; ++i) {
                for (Indx j = i + 1; j < Self::num_nodes; ++j) {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (Indx k = 0; k < Self::dim_shape; ++k) {
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

        // --- REFERENCE CONFIG

//        static
//        ReferenceNodesCoordinatesViewT
//        getReferenceCoordinates()
//        {
//            return ElementConnectivityT::getReferenceCoordinates();
//        }

        static
        Real
        getReferenceDistance(
                auto const &
                first_point_arg,
                auto const &
                second_point_arg,
                Intg
                direction_arg = -1
        )
        {
            assert((direction_arg <= static_cast<Intg>(Self::ord_shape)));
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
        Matrix<Real, Self::ord_shape>
        getReferenceCentroid()
        {
            auto reference_nodes_coordinates = Self::getReferenceCoordinates();
            return core::geometry::getBarycenter(reference_nodes_coordinates);
        }

        static
        Matrix<Real, Self::ord_shape>
        getReferenceDiameters()
        {
            auto reference_nodes_coordinates = Self::getReferenceCoordinates();
            auto reference_diameters = Matrix<Real, Self::ord_shape>().setZero();
            for (Indx i = 0; i < Self::num_nodes; ++i) {
                for (Indx j = i + 1; j < Self::num_nodes; ++j) {
                    for (Indx k = 0; k < Self::ord_shape; ++k) {
                        auto & a = reference_diameters(k);
                        auto b = std::abs(reference_nodes_coordinates(k, i) - reference_nodes_coordinates(k, j));
                        if (b > a) {
//                            reference_diameters(k) = b;
                            a = b;
                        }
                    }
                }
            }
            return reference_diameters;
        }

        /*
         * QUADRATURE
         */

        template<QuadratureDescription Q>
        static constexpr
        Indx
        getDimQuadrature()
        {
            return ShapeQuadrature<E.shape_description, Q>::dim_quadrature;
        }

        template<QuadratureDescription Q>
        static
        Real
        getReferenceQuadratureWeight(
                Indx
                index_arg
        )
        {
            assert(index_arg <= getDimQuadrature<Q>());
            return ShapeQuadrature<E.shape_description, Q>::reference_weights(index_arg);
        }

        template<QuadratureDescription Q>
        static
        auto
        getReferenceQuadraturePoint(
                Indx
                index_arg
        )
        {
            assert(index_arg <= getDimQuadrature<Q>());
            auto offset = Self::dim_shape * index_arg;
//            Real const constexpr * const d = ShapeQuadrature<E.shape_description, Q>::reference_points.data.data();
            auto const constexpr * const d = ShapeQuadrature<E.shape_description, Q>::reference_points.data.data();
            return Matrix<Matrix<Real, Self::dim_shape> const>(d + offset);
        }

        template<QuadratureDescription Q, Indx I, Indx J>
        static
        Matrix<Real, Self::dim_shape>
        getComponentReferenceQuadraturePoint(
                Indx
                component_index_arg,
                Indx
                index_arg
        )
        {
//            using AAA = ComponentReferenceQuadraturePolicy<Q, I, J>;
//            using BBB = typename AAA::ComponentReferenceQuadrature;
//            assert(index_arg < BBB::size);
//            return AAA::getComponentReferencePoint(component_index_arg, index_arg);
////            Matrix<Real, Util::dim_shape> component_reference_quadrature_point;
////            for (Indx j = 0; j < Util::dim_nodes; ++j) {
////                auto const refpt = AAA::getElementComponentReferencePoint(component_index_arg, i);
////                component_reference_quadrature_point(j) = refpt(j);
////            }
////            return component_reference_quadrature_point;
            using Component = typename Self::Components::template Type<I>::template Type<J>::Type::Type;
            using ComponentGeometry = ElementBuild<Component::element_description, M>;
            auto component_reference_point1 = Matrix<Real, Self::dim_shape>();
            auto const & components_node_tags = Self::node_connectivity.template get<I>().template get<J>();
            auto const & element_reference_nodes = Self::reference_nodes;
            for (Indx i = 0; i < Self::dim_shape; ++i) {
                auto component_coordinates = Matrix<Real, Component::num_nodes>();
                for (Indx j = 0; j < Component::num_nodes; ++j) {
//                    auto const node_tag = components_node_tags.get(component_index_arg).get(j);
                    auto const node_tag = components_node_tags.get(component_index_arg, j);
                    component_coordinates(j) = element_reference_nodes.get(node_tag, i);
                }
                auto component_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<Q>(index_arg);
                component_reference_point1(i) = ComponentGeometry::getShapeMappingEvaluation(component_coordinates, component_reference_point);
            }
            return component_reference_point1;
        }

        template<QuadratureDescription Q>
        Real
        getCurrentQuadratureWeight(
                Indx
                index_arg
        )
        const
        {
            assert(index_arg < getDimQuadrature<Q>());
            auto reference_weight = getReferenceQuadratureWeight<Q>(index_arg);
            return reference_weight * getShapeMappingDifferential(getReferenceQuadraturePoint<Q>(index_arg));
        }

        template<QuadratureDescription Q>
        Matrix<Real, Self::dim_euclidean>
        getCurrentQuadraturePoint(
                Indx
                index_arg
        )
        const
        {
            assert(index_arg <= getDimQuadrature<Q>());
            auto p = Matrix<Real, Self::dim_euclidean>();
            auto const nds = this->getCurrentCoordinates();
            for (Indx j = 0; j < Self::dim_euclidean; ++j) {
                p(j) = Self::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<Q>(index_arg));
            }
            return p;
        }

        template<QuadratureDescription Q, Indx I, Indx J>
        Matrix<Real, Self::dim_euclidean>
        getComponentCurrentQuadraturePoint(
                Indx
                component_index_arg,
                Indx
                index_arg
        )
        const
        {
            auto const rp11 = getComponentReferenceQuadraturePoint<Q, I, J>(component_index_arg, index_arg);
            auto p = Matrix<Real, Self::dim_euclidean>();
            auto const nds = this->getCurrentCoordinates();
            for (Indx j = 0; j < Self::dim_euclidean; ++j) {
                p(j) = Self::getShapeMappingEvaluation(nds.row(j), rp11);
            }
            return p;
//            using AAA = ComponentReferenceQuadraturePolicy<Q, I, J>;
//            using BBB = typename AAA::ComponentReferenceQuadrature;
//            assert(index_arg < BBB::size);
////            auto const rp = AAA::getElementComponentReferencePoint(component_index_arg, i);
//            auto const rp1 = getComponentReferenceQuadraturePoint<Q>(component_index_arg, index_arg);
////            using AAA = typename ComponentReferenceQuadraturePointsPolicy<
////                    Qt,
////                    Q,
////                    I,
////                    J
////            >::ComponentQuadratureRule;
////            raise(i >= AAA::getQuadratureSize(), i, "is too large, max :", AAA::getQuadratureSize());
////            ReferencePointT component_reference_quadrature_point = getComponentReferenceQuadraturePoint<Qt, Q, I, J>(
////                    component_index_arg,
////                    i
////            );
//            auto p = Matrix<Real, Self::dim_euclidean>();
//            auto const nds = this->getCurrentCoordinates();
//            for (Indx j = 0; j < Self::dim_euclidean; ++j) {
//                p(j) = getShapeMappingEvaluation(nds.row(j), rp1);
//            }
//            return p;
//            using Component = typename Self::Components::template Type<I>::template Type<J>::Type::Type;
//            using ComponentGeometry = ElementBuild<Component::element_description, M>;
//            auto component_reference_point1 = Matrix<Real, Self::dim_shape>();
//            auto const & components_node_tags = Self::node_connectivity.template get<I>().template get<J>();
//            auto const & element_reference_nodes = Self::reference_nodes;
//            for (Indx i = 0; i < Self::dim_shape; ++i) {
//                auto component_coordinates = Matrix<Real, Component::num_nodes>();
//                for (Indx j = 0; j < Component::num_nodes; ++j) {
////                    auto const node_tag = components_node_tags.get(component_index_arg).get(j);
//                    auto const node_tag = components_node_tags.get(component_index_arg, j);
//                    component_coordinates(j) = element_reference_nodes.get(node_tag, i);
//                }
//                auto component_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<Q>(index_arg);
//                component_reference_point1(i) = ComponentGeometry::getShapeMappingEvaluation(component_coordinates, component_reference_point);
//            }
//            return component_reference_point1;
        }

    private:

        Real
        getCurvedDistance(
                auto const &
                first_point_arg,
                auto const &
                second_point_arg,
                Intg
                direction_arg = -1
        )
        const
        {
            static_assert(Self::element_description != pnt_0);
            auto const constexpr q_seg = QuadratureDescription(Quadrature::Gauss, 4);
//            using SegmentQuadratureRule = ShapeQuadrature<seg, q_seg>;
//            using DerivativeVectors = Matrix<Real, Self::dim_euclidean, Self::ord_shape>;
//                Indx const constexpr dim_shape = ElementPolicyT::getDimShape();
//                Indx const constexpr quadrature_order = SegmentQuadratureRule::getQuadratureOrder();
//                Indx const constexpr quadrature_size = SegmentQuadratureRule::getQuadratureSize();
//            Real distance = 0.0;
//            Real dt;
            auto distance = Real(0);
            auto dt = Real();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (Indx q = 0; q < ShapeQuadrature<shape::seg, q_seg>::dim_quadrature; ++q) {
                auto pq = ShapeQuadrature<shape::seg, q_seg>::reference_points(q, 0);
                auto wq = ShapeQuadrature<shape::seg, q_seg>::reference_weights(q);
//                auto ru = Matrix<Real, Self::dim_euclidean, Self::dim_shape>::Zero();
                auto ru = Matrix<Real, Self::dim_euclidean, Self::ord_shape>();
                auto difference = second_point_arg - first_point_arg;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                for (Indx i = 0; i < Self::dim_euclidean; ++i) {
                    for (Indx j = 0; j < Self::ord_shape; ++j) {
                        if (direction_arg == -1 || i == static_cast<Indx>(direction_arg)) {
                            auto du = (1.0 / 2.0) * (second_point_arg(j) - first_point_arg(j));
                            auto dx = Self::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                            ru(i, j) = dx * du;
                        }
                    }
                }
                if constexpr (Self::ord_shape == 1) {
                    auto Eff = ru.col(0).template dot(ru.col(0));
                    dt = std::sqrt(Eff);
                } else if constexpr (Self::ord_shape == 2) {
                    auto Eff = ru.col(0).template dot(ru.col(0));
                    auto Fff = ru.col(0).template dot(ru.col(1));
                    auto Gff = ru.col(1).template dot(ru.col(1));
                    dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                }
                distance += wq * dt;
            }
            return distance;
        }

        Real
        getCartesianDistance(
                auto const &
                first_point_arg,
                auto const &
                second_point_arg,
                Intg
                direction_arg = -1
        )
        const
        {
            auto const & nds = this->getCurrentCoordinates();
            auto distance = Real();
            auto mp0 = Matrix<Real, Self::dim_euclidean>();
            auto mp1 = Matrix<Real, Self::dim_euclidean>();
            for (Indx i = 0; i < Self::dim_euclidean; ++i) {
                mp0(i) = Self::getShapeMappingEvaluation(nds.row(i), first_point_arg);
                mp1(i) = Self::getShapeMappingEvaluation(nds.row(i), second_point_arg);
            }
            if (direction_arg == -1) {
                distance = (mp1 - mp0).norm();
            } else {
                distance = (mp1 - mp0)(direction_arg);
            }
            return distance;
        }

    public:

        Matrix<Real, Self::dim_euclidean>
        getNormalVector(
                auto const &
                point_arg
        )
        const
        {
            static_assert(Self::ord_shape == Self::dim_euclidean - 1);
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
//            auto derivative_vector = Matrix<Real, 3, Self::dim_nodes>::Zero();
            auto derivative_vector = Matrix<Real, 3, Self::dim_shape>().setZero();
            for (Indx i = 0; i < Self::dim_euclidean; ++i) {
                for (Indx j = 0; j < Self::dim_shape; ++j) {
                    derivative_vector(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, j);
                }
            }
            auto const & block = derivative_vector.template block<Self::dim_euclidean, Self::dim_shape>(0, 0);
            auto normal_vector = core::geometry::getNormalVector(block);
            return normal_vector;
        }

        Matrix<Real, Self::dim_euclidean>
        getTangentVector(
                auto const &
                point_arg,
                Indx
                direction_arg
        )
        const
        {
            static_assert(Self::ord_shape == Self::dim_euclidean - 2);
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = Matrix<Real, Self::dim_euclidean>();
            for (Indx i = 0; i < Self::dim_euclidean; ++i) {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, direction_arg);
            }
            return tangent_vector;
        }

    };

}

#endif //FETA_FETA_CORE_ELEMENT_FE_BUILD_HXX
