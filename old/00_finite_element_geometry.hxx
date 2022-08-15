//
// Created by dsiedel on 16/03/2022.
//

#ifndef FETA_00_FINITE_ELEMENT_GEOMETRY_HXX
#define FETA_00_FINITE_ELEMENT_GEOMETRY_HXX

#include "feta/core/00_finite_element_description_hho.hxx"
#include "feta/core/element_final.hxx"
#include "feta/core/frame_type.hxx"
//#include "feta/core/discrete_field.hxx"
#include "feta/core/finite_element_policy.hxx"
#include "feta/core/00_finite_element_connectivity.hxx"

namespace feta::core::internal_good
{

    template<
            ElementType ElementTypeArg,
            SupportType SupportTypeArg,
            auto MixedElementDescriptionArg
    >
    struct FiniteElementGeometry
    :
    public
    FiniteElementConnectivity<
            ElementTypeArg,
            SupportTypeArg,
            MixedElementDescriptionArg
    >
    {

        using Base = FiniteElementConnectivity<
                ElementTypeArg,
                SupportTypeArg,
                MixedElementDescriptionArg
        >;

        void
        sayHello()
        {
            print("Hi !!!");
            print(this->getReferenceCoordinates());
        }

//        using ElementConnectivityT = internal_good::FiniteElementConnectivity<
//                ElementTypeArg,
//                SupportTypeArg,
//                MixedElementDescriptionArg
//        >;

        using NodalValuesT = NodalValues<
                ElementTypeArg
        >;

        using ReferencePointT = ReferencePoint<
                ElementTypeArg
        >;

        using ElementPolicyT = ElementPE<
                ElementTypeArg
        >;

        using CurrentNodesCoordinatesT = CurrentNodesCoordinates<
                ElementTypeArg,
                MixedElementDescriptionArg.getDimEuclidean()
        >;

        using ReferenceNodesCoordinatesViewT = ReferenceNodesCoordinatesView<
                ElementTypeArg
        >;

        using NormalVectorT = StaticColVector<
                MixedElementDescriptionArg.getDimEuclidean()
        >;

        using TangentVectorT = StaticColVector<
                MixedElementDescriptionArg.getDimEuclidean()
        >;

        using CurrentPointT = CurrentPoint<
                MixedElementDescriptionArg.getDimEuclidean()
        >;

        using ReferencePointViewT = ReferencePointView<
                ElementTypeArg
        >;

    private:

        template<
                QuadratureType Qt,
                Indx Q,
                Indx I,
                Indx J
        >
        struct ComponentReferenceQuadraturePointsPolicy
        {

        private:

            template<
                    ElementType ElementTypeA
            >
            using FaceQuadratureRule = QuadratureRule<
                    ElementTypeA,
                    Qt,
                    Q
            >;

            using FaceQuadratureRuleComponents = typename ElementPolicyT::template ComponentsCollection<
                    FaceQuadratureRule
            >::template Type<
                    I
            >::template Type<
                    J
            >;

            static constexpr
            ElementType
            getElementComponentType()
            {
                return FaceQuadratureRuleComponents::getElementType();
            }

            using ElementComponentPolicyT = ElementPE<
                    getElementComponentType()
            >;

            using ComponentNodesCoordinatesArrayT = StaticArray<
                    Real,
                    ElementComponentPolicyT::getNumNodes(),
                    ElementPolicyT::getDimNodes()
            >;

            static constexpr
            ComponentNodesCoordinatesArrayT
            getElementComponentNodesCoordinatesArray(
                    Indx
                    component_index_arg
            )
            {
                ComponentNodesCoordinatesArrayT nodes_coordinates_array;
                for (Indx i = 0; i < ElementComponentPolicyT::getNumNodes(); ++i) {
//                    auto const & table = ElementPolicyT::getNodesConnectivity().template operator()<I>().template operator()<J>();
                    auto const & co = ElementPolicyT::getNodesConnectivity().template asMatrix<I>().
                            template asMatrix<J>();
                    auto const ghj = co.asMatrix(component_index_arg).asMatrix(i);
//                    Indx a = ElementPolicyT::getNodesConnectivity().template operator()<I>().template operator()<J>()(
//                            component_index_arg
//                    )(i);
                    for (Indx j = 0; j < ElementPolicyT::getDimNodes(); ++j) {
                        auto const reference_coordinates_array = ElementPolicyT::getReferenceNodesCoordinatesArray();
                        nodes_coordinates_array(i, j) = reference_coordinates_array.asMatrix(ghj, j);
//                        nodes_coordinates_array(i, j) = ElementPolicyT::getReferenceNodesCoordinatesArray()(ghj, j);
                    }
                }
                return nodes_coordinates_array;
            }

        public:

            using ComponentQuadratureRule = QuadratureRule<
                    getElementComponentType(),
                    Qt,
                    Q
            >;

            static constexpr
            ReferencePointT
            getElementComponentReferencePoint(
                    Indx
                    component_index_arg,
                    Indx
                    q
            )
            {
                ReferencePointT oneo;
                ComponentNodesCoordinatesArrayT component_nodes = getElementComponentNodesCoordinatesArray(
                        component_index_arg
                );
                for (Indx i = 0; i < ElementPolicyT::getDimNodes(); ++i) {
                    NodalValuesArray<getElementComponentType()> nodal_coordinates;
                    for (Indx j = 0; j < ElementComponentPolicyT::getNumNodes(); ++j) {
                        nodal_coordinates(j) = component_nodes(j, i);
                    }
                    ReferencePointArray<getElementComponentType()> reference_point_array;
                    for (Indx k = 0; k < reference_point_array.getSize(); ++k) {
                        reference_point_array(k) = ComponentQuadratureRule::getReferencePoints()(q, k);
                    }
                    oneo(i) = ElementComponentPolicyT::getMappingEvaluation(nodal_coordinates, reference_point_array);
                }
                return oneo;
            }

        };

    public:

//        ElementGeometry()
//                :
//                element_connectivity(ElementConnectivityPointerT(ElementConnectivityT()))
//        {}
//
//        explicit
//        ElementGeometry(
//                ElementConnectivityPointerT const &
//        a_element_geometry
//        )
//        :
//        element_connectivity(a_element_geometry)
//                {}
//
//        explicit
//        ElementGeometry(
//                ElementConnectivityPointerT &&
//        a_element_geometry
//        )
//        :
//        element_connectivity(a_element_geometry)
//                {}

        static
        Real
        getMappingEvaluation(
                NodalValuesT const &
                nodal_values_arg,
                ReferencePointT const &
                point_arg
        )
        {
            return ElementPolicyT::getMappingEvaluation(nodal_values_arg, point_arg);
        }

        static
        Real
        getMappingDerivative(
                NodalValuesT const &
                nodal_values_arg,
                ReferencePointT const &
                point_arg,
                Indx
                derivative_direction_arg
        )
        {
            return ElementPolicyT::getMappingDerivative(nodal_values_arg, point_arg, derivative_direction_arg);
        }

//        static constexpr
//        GeometryType
//        getGeometryType()
//        {
//            Indx const constexpr dim_shape = ElementPolicyT::getDimShape();
//            Indx const constexpr dim_euclidean = D;
//            GeometryType geometry_type;
//            if constexpr (dim_shape == D) {
//                geometry_type = GeometryType::Cell;
//            } else if constexpr (dim_shape == dim_euclidean - 1) {
//                geometry_type = GeometryType::Face;
//            } else if constexpr (dim_shape == dim_euclidean - 2) {
//                geometry_type = GeometryType::Edge;
//            } else {
//                geometry_type = GeometryType::Node;
//            }
//            return geometry_type;
//        }

//        GeometryType const static constexpr geometry_type = GeometryType::Cell;

        // --- CURRENT CONFIG

//        CurrentNodesCoordinatesT
//        getCurrentCoordinates()
//        const
//        {
//            return element_connectivity().getCurrentCoordinates();
//        }

        Real
        getMappingDifferential(
                ReferencePointT const &
                point_arg
        )
        const
        {
            Indx const constexpr ajj = ElementPolicyT::getDimShape();
            using HJKO = StaticMatrix<
                    3,
                    ajj
            >;
            CurrentNodesCoordinatesT const nds = this->getCurrentCoordinates();
            HJKO ru = HJKO::Zero();
            Real du;
            for (Indx i = 0; i < MixedElementDescriptionArg.getDimEuclidean(); ++i) {
                for (Indx j = 0; j < ajj; ++j) {
                    ru(i, j) = getMappingDerivative(nds.row(i), point_arg, j);
                }
            }
            if constexpr (ajj == 1) {
                du = std::abs(ru.col(0).norm());
            }
            else if constexpr (ajj == 2) {
                du = std::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else if constexpr (ajj == 3) {
                du = std::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
//            if constexpr (Ft == FrameType::AxiSymmetric) {
//                Real r0 = getMappingEvaluation(nds.row(0), point);
//                if (r0 < 1.e-10) {
//                    r0 = 1.e-10;
//                }
//                du *= 2.0 * PI_ * r0;
//            }
            return du;
        }

        Real
        getCurrentDistance(
                ReferencePointT const &
                first_point_arg,
                ReferencePointT const &
                second_point_arg,
                Indx
                direction_arg
        )
        const
        {
            raise((direction_arg >= ElementPolicyT::getDimShape()), "dir must be lesser than ", ElementPolicyT::getDimShape());
//            CurrentNodesCoordinates<ElementTypeArg, D> const nds = getCurrentCoordinates();
            Real distance;
            if constexpr (MixedElementDescriptionArg.getDimEuclidean() > ElementPolicyT::getDimShape()) {
                distance = getLocalFrameDistance(first_point_arg, second_point_arg, direction_arg);
            } else {
                distance = getGlobalFrameDistance(first_point_arg, second_point_arg, direction_arg);
            }
            return distance;
        }

        Real
        getCurrentDistance(
                ReferencePointT const &
                first_point_arg,
                ReferencePointT const &
                second_point_arg
        )
        const
        {
            Real distance;
            if constexpr (MixedElementDescriptionArg.getDimEuclidean() > ElementPolicyT::getDimShape()) {
                distance = getLocalFrameDistance(first_point_arg, second_point_arg, 3);
            } else {
                distance = getGlobalFrameDistance(first_point_arg, second_point_arg, 3);
            }
            return distance;
        }

        ReferencePointT
        getCurrentDiameters()
        const
        {
            ReferenceNodesCoordinatesViewT reference_nodes_coordinates = Base::getReferenceCoordinates();
            ReferencePointT current_diameters = ReferencePointT::Zero();
            for (Indx i = 0; i < ElementPolicyT::getNumNodes(); ++i) {
                for (Indx j = 0; j < ElementPolicyT::getNumNodes(); ++j) {
                    if (j > i) {
                        for (Indx k = 0; k < ElementPolicyT::getNumNodes(); ++k) {
                            Real b = std::abs(getCurrentDistance(
                                    reference_nodes_coordinates.col(i),
                                    reference_nodes_coordinates.col(j),
                                    k
                            ));
                            if (b > std::abs(current_diameters(k))) {
                                current_diameters(k) = b;
                            }
                        }
                    }
                }
            }
            return current_diameters;
        }

        CurrentPointT
        getCurrentCentroid()
        const
        {
            CurrentNodesCoordinatesT const current_nodes_coordinates = this->getCurrentCoordinates();
            return geometry::getBarycenter<MixedElementDescriptionArg.getDimEuclidean(), ElementPolicyT::getNumNodes()>(current_nodes_coordinates);
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
                ReferencePointT const &
                first_point_arg,
                ReferencePointT const &
                second_point_arg,
                Indx
                direction_arg
        )
        {
            raise((direction_arg >= ElementPolicyT::getDimShape()), "dir must be lesser than ", ElementPolicyT::getDimShape());
            Real distance;
            distance = (second_point_arg - first_point_arg)(direction_arg);
            return distance;
        }

        static
        Real
        getReferenceDistance(
                ReferencePointT const &
                first_point_arg,
                ReferencePointT const &
                second_point_arg
        )
        {
            Real distance;
            distance = (second_point_arg - first_point_arg).norm();
            return distance;
        }

        static
        ReferencePointT
        getReferenceCentroid()
        {
            ReferenceNodesCoordinatesViewT reference_nodes_coordinates = Base::getReferenceCoordinates();
            return geometry::getBarycenter<MixedElementDescriptionArg.getDimEuclidean(), ElementPolicyT::getDimShape()>(reference_nodes_coordinates);
        }

        static
        ReferencePointT
        getReferenceDiameters()
        {
            ReferenceNodesCoordinatesViewT reference_nodes_coordinates = Base::getReferenceCoordinates();
            ReferencePointT reference_diameters = ReferencePointT::Zero();
            for (Indx i = 0; i < ElementPolicyT::getNumNodes(); ++i) {
                for (Indx j = 0; j < ElementPolicyT::getNumNodes(); ++j) {
                    if (j > i) {
                        for (Indx k = 0; k < ElementPolicyT::getNumNodes(); ++k) {
                            Real b = std::abs(
                                    reference_nodes_coordinates(k, i) - reference_nodes_coordinates.col(k, j)
                            );
                            if (b > std::abs(reference_diameters(k))) {
                                reference_diameters(k) = b;
                            }
                        }
                    }
                }
            }
            return reference_diameters;
        }

        // --- QUADRATURE

        template<
                QuadratureType Qt,
                Indx Q
        >
        static constexpr
        Indx
        getQuadratureSize()
        {
            using QuadratureRuleT = QuadratureRule<
                    ElementTypeArg,
                    Qt,
                    Q
            >;
            return QuadratureRuleT::getQuadratureSize();
        }

//        /**
//         * @brief
//         * @tparam Qt
//         * @param Q
//         * @return
//         */
//        template<
//                QuadratureType Qt
//        >
//        static
//        Indx
//        getQuadratureSize(
//                Indx
//                Q
//        )
//        {
//            Indx res;
//            switch (Q) {
//                case 0:
//                    res = QuadratureRule<ElementTypeArg, Qt, 0>::getQuadratureSize();
//                    break;
//                case 1:
//                    res = QuadratureRule<ElementTypeArg, Qt, 1>::getQuadratureSize();
//                    break;
//                case 2:
//                    res = QuadratureRule<ElementTypeArg, Qt, 2>::getQuadratureSize();
//                    break;
//                case 3:
//                    res = QuadratureRule<ElementTypeArg, Qt, 3>::getQuadratureSize();
//                    break;
//                case 4:
//                    res = QuadratureRule<ElementTypeArg, Qt, 4>::getQuadratureSize();
//                    break;
//                case 5:
//                    res = QuadratureRule<ElementTypeArg, Qt, 5>::getQuadratureSize();
//                    break;
//                case 6:
//                    res = QuadratureRule<ElementTypeArg, Qt, 6>::getQuadratureSize();
//                    break;
//                case 7:
//                    res = QuadratureRule<ElementTypeArg, Qt, 7>::getQuadratureSize();
//                    break;
//                case 8:
//                    res = QuadratureRule<ElementTypeArg, Qt, 8>::getQuadratureSize();
//                    break;
//                default:
//                    raise(false, "unsupported quadrature order :", Q);
//            }
//            return res;
//        }

        template<
                QuadratureType Qt,
                Indx Q
        >
        static
        Real
        getReferenceQuadratureWeight(
                Indx
                i
        )
        {
            using QuadratureRuleT = QuadratureRule<
                    ElementTypeArg,
                    Qt,
                    Q
            >;
            raise(i >= getQuadratureSize<Qt, Q>(), i, "is too large, max :", getQuadratureSize<Qt, Q>());
            return QuadratureRuleT::reference_weights(i);
        }

        template<
                QuadratureType Qt,
                Indx Q
        >
        static
        ReferencePointViewT
        getReferenceQuadraturePoint(
                Indx
                i
        )
        {
            using QuadratureRuleT = QuadratureRule<
                    ElementTypeArg,
                    Qt,
                    Q
            >;
            raise(i >= getQuadratureSize<Qt, Q>(), i, "is too large, max :", getQuadratureSize<Qt, Q>());
            Indx offset = ElementPolicyT::getDimNodes() * i;
            Real const constexpr * const d = QuadratureRuleT::reference_points.data.data();
            return ReferencePointViewT(d + offset);
        }

        template<
                QuadratureType Qt,
                Indx Q,
                Indx I,
                Indx J
        >
        static
        ReferencePointT
        getComponentReferenceQuadraturePoint(
                Indx
                component_index_arg,
                Indx
                i
        )
        {
            using AAA = ComponentReferenceQuadraturePointsPolicy<
                    Qt,
                    Q,
                    I,
                    J
            >;
            using BBB = typename AAA::ComponentQuadratureRule;
            raise(i >= BBB::getQuadratureSize(), i, "is too large, max :", BBB::getQuadratureSize());
            ReferencePointT component_reference_quadrature_point;
            for (Indx j = 0; j < ElementPolicyT::getDimNodes(); ++j) {
                component_reference_quadrature_point(j) = AAA::getElementComponentReferencePoint(
                        component_index_arg,
                        i
                )(j);
            }
            return component_reference_quadrature_point;
        }

        template<
                QuadratureType Qt,
                Indx Q
        >
        Real
        getCurrentQuadratureWeight(
                Indx
                i
        )
        const
        {
            raise(i >= getQuadratureSize<Qt, Q>(), i, "is too large, max :", getQuadratureSize<Qt, Q>());
//            CurrentNodesCoordinates<ElementTypeArg, D> const nds = getCurrentCoordinates();
            Real w = getReferenceQuadratureWeight<Qt, Q>(i);
            return w * getMappingDifferential(getReferenceQuadraturePoint<Qt, Q>(i));
        }

        template<
                QuadratureType Qt,
                Indx Q
        >
        CurrentPointT
        getCurrentQuadraturePoint(
                Indx
                i
        )
        const
        {
            raise(i >= getQuadratureSize<Qt, Q>(), i, "is too large, max :", getQuadratureSize<Qt, Q>());
            CurrentPointT p;
            CurrentNodesCoordinatesT const nds = this->getCurrentCoordinates();
            for (Indx j = 0; j < MixedElementDescriptionArg.getDimEuclidean(); ++j) {
                p(j) = getMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<Qt, Q>(i));
            }
            return p;
        }

        template<
                QuadratureType Qt,
                Indx Q,
                Indx I,
                Indx J
        >
        ReferencePointT
        getComponentCurrentQuadraturePoint(
                Indx
                component_index_arg,
                Indx
                i
        )
        const
        {
            using AAA = typename ComponentReferenceQuadraturePointsPolicy<
                    Qt,
                    Q,
                    I,
                    J
            >::ComponentQuadratureRule;
            raise(i >= AAA::getQuadratureSize(), i, "is too large, max :", AAA::getQuadratureSize());
            ReferencePointT component_reference_quadrature_point = getComponentReferenceQuadraturePoint<Qt, Q, I, J>(
                    component_index_arg,
                    i
            );
            CurrentPointT p;
            CurrentNodesCoordinatesT const nds = this->getCurrentCoordinates();
            for (Indx j = 0; j < MixedElementDescriptionArg.getDimEuclidean(); ++j) {
                p(j) = getMappingEvaluation(nds.row(j), component_reference_quadrature_point);
            }
            return p;
        }

    private:

        Real
        getLocalFrameDistance(
                ReferencePointT const &
                first_point_arg,
                ReferencePointT const &
                second_point_arg,
                Indx
                direction_arg = 3
        )
        const
        {
            using SegmentQuadratureRule = QuadratureRule<
                    ElementType::LinearSegment,
                    QuadratureType::Gauss,
                    4
            >;
            using DerivativeVectors = StaticMatrix<
                    MixedElementDescriptionArg.getDimEuclidean(),
                    ElementPolicyT::getDimShape()
            >;
            Indx const constexpr dim_shape = ElementPolicyT::getDimShape();
            Indx const constexpr quadrature_order = SegmentQuadratureRule::getQuadratureOrder();
            Indx const constexpr quadrature_size = SegmentQuadratureRule::getQuadratureSize();
            CurrentNodesCoordinatesT const current_nodes_coordinates = this->getCurrentCoordinates();
            Real distance = 0.0;
            Real dt;
            for (Indx q = 0; q < quadrature_size; ++q) {
                Real pq = SegmentQuadratureRule::getReferencePoints()(q, 0);
                Real wq = SegmentQuadratureRule::getReferenceWeights()(q);
                DerivativeVectors ru = DerivativeVectors::Zero();
                Real difference = second_point_arg - first_point_arg;
                ReferencePointT uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                for (Indx i = 0; i < MixedElementDescriptionArg.getDimEuclidean(); ++i) {
                    for (Indx j = 0; j < dim_shape; ++j) {
                        if (direction_arg == 3 || i == direction_arg) {
                            Real du = (1.0 / 2.0) * (second_point_arg(j) - first_point_arg(j));
                            Real dx = getMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                            ru(i, j) = dx * du;
                        }
                    }
                }
                if constexpr (dim_shape == 1) {
                    Real E = ru.col(0).template dot(ru.col(0));
                    dt = std::sqrt(E);
                } else if constexpr (dim_shape == 2) {
                    Real E = ru.col(0).template dot(ru.col(0));
                    Real F = ru.col(0).template dot(ru.col(1));
                    Real G = ru.col(1).template dot(ru.col(1));
                    dt = std::sqrt(E + 2.0 * F + G);
                }
                distance += wq * dt;
            }
            return distance;
        }

        Real
        getGlobalFrameDistance(
                ReferencePointT const &
                first_point_arg,
                ReferencePointT const &
                second_point_arg,
                Indx
                direction_arg = 3
        )
        const
        {
            CurrentNodesCoordinatesT const nds = this->getCurrentCoordinates();
            Real distance;
            CurrentPointT mp0 = CurrentPointT::Zero();
            CurrentPointT mp1 = CurrentPointT::Zero();
            for (int i = 0; i < MixedElementDescriptionArg.getDimEuclidean(); ++i) {
                mp0(i) = getMappingEvaluation(nds.row(i), first_point_arg);
                mp1(i) = getMappingEvaluation(nds.row(i), second_point_arg);
            }
            if (direction_arg == 3) {
                distance = (mp1 - mp0).norm();
            } else {
                distance = (mp1 - mp0)(direction_arg);
            }
            return distance;
        }

        NormalVectorT
        getNormalVector(
                ReferencePointT const &
                point
        )
        const
        {

//            static_assert(getGeometryType() == GeometryType::Face);
            static_assert(ElementPolicyT::getDimShape() == MixedElementDescriptionArg.getDimEuclidean() - 1);

            Indx const constexpr dim_nodes = ElementPolicyT::getDimNodes();
            CurrentNodesCoordinatesT const current_nodes_coordinates = this->getCurrentCoordinates();
            Indx const constexpr dn = dim_nodes;
            StaticMatrix<3, dn> derivative_vector = StaticMatrix<3, dn>::Zero();
            for (Indx i = 0; i < MixedElementDescriptionArg.getDimEuclidean(); ++i) {
                for (Indx j = 0; j < dim_nodes; ++j) {
                    derivative_vector(i, j) = this->getMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            NormalVectorT normal_vector = geometry::getNormalVector<MixedElementDescriptionArg.getDimEuclidean()>(
                    derivative_vector.template block<MixedElementDescriptionArg.getDimEuclidean(), dn>(0, 0)
            );
            return normal_vector;
        }

        TangentVectorT
        getTangentVector(
                ReferencePointT const &
                point,
                Indx
                dir
        )
        const
        {

//            static_assert(getGeometryType() == GeometryType::Edge);
            static_assert(ElementPolicyT::getDimShape() == MixedElementDescriptionArg.getDimEuclidean() - 2);

            CurrentNodesCoordinatesT const current_nodes_coordinates = this->getCurrentCoordinates();
            TangentVectorT tangent_vector;
            for (Indx i = 0; i < MixedElementDescriptionArg.getDimEuclidean(); ++i) {
                tangent_vector(i) = this->getMappingDerivative(current_nodes_coordinates.row(i), point, dir);
            }
            return tangent_vector;
        }

    };

}

#endif //FETA_00_FINITE_ELEMENT_GEOMETRY_HXX
