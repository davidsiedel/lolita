//
// Created by dsiedel on 09/05/22.
//

#ifndef LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
#define LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_element.hxx"

namespace lolita::core::element
{

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    struct TensorPolicy;

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 0)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr cardinality_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 0),
                lolita::numerics::pow(_dim_euclidean, 0)
        };

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 1)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr cardinality_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 1),
                lolita::numerics::pow(_dim_euclidean, 0)
        };

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 2)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr cardinality_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 1),
                lolita::numerics::pow(_dim_euclidean, 1)
        };

    };

    template<Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto>
    struct FiniteElementPolicy;

    template<Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto>
    struct FiniteElementModule;

    template<template<Element, lolita::geometry::Domain, auto...> typename, Element, lolita::geometry::Domain, auto...>
    struct FiniteElementConnectivity;

    template<template<Element, lolita::geometry::Domain, auto...> typename, Element, lolita::geometry::Domain, auto...>
    struct FiniteElementGeometry;

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, Element _element, lolita::geometry::Domain _domain, auto... _arg>
    requires(!PointConcept<_element>)
    struct FiniteElementConnectivity<_T, _element, _domain, _arg...> : public ElementGeometry<_element>
    {

    private:

        template<Element __element, lolita::geometry::Domain __domain, auto... __arg>
        using _FiniteElementPointer = std::shared_ptr<_T<__element, __domain, __arg...>>;

    public:

        using Neighbours = typename ElementGeometry<_element>::template Neighbours<_FiniteElementPointer, _domain, _arg...>;

        using Components = typename ElementGeometry<_element>::template Components<_FiniteElementPointer, _domain, _arg...>;

        template<auto _i, auto _j>
        using NeighbourFamily =  std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>>;

        FiniteElementConnectivity()
        :
                neighbours_(),
                components_(),
                tag_()
        {}

        lolita::boolean
        operator==(
                FiniteElementConnectivity const &
                other
        )
        const
        {
            return this->tag_ == other.tag_;
        };

        lolita::boolean
        operator!=(
                FiniteElementConnectivity const &
                other
        )
        const
        {
            return !(* this == other);
        }

        /*
         *
         */

        template<auto _i, auto _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        template<auto _i, auto _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        template<auto _i, auto _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Components>> &
        getComponents()
        {
            return std::get<_j>(std::get<_i>(components_));
        }

        template<auto _i, auto _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Components>> const &
        getComponents()
        const
        {
            return std::get<_j>(std::get<_i>(components_));
        }

        template<lolita::index _i, lolita::index _j>
        static
        lolita::index
        getComponentNodeConnection(
                lolita::index i,
                lolita::index j
        )
        {
            return std::get<_j>(std::get<_i>(ElementGeometry<_element>::node_connectivity_))[i][j];
        }

//        template<auto _i, auto _j>
//        lolita::index
//        getComponentIndex(
//                lolita::index
//                index
//        )
//        const
//        {
//            auto const constexpr _neighbour = element::neighbour<0, _i, _j>();
//            auto const constexpr _pos_neighbour = element::neighbourPosition<_neighbour, _element>();
//            auto count = lolita::index(0);
//            for (auto const & item : getNeighbours2<_pos_neighbour[0], _pos_neighbour[1], _pos_neighbour[2]>()) {
////            for (auto const & item : std::get<_pos_neighbour[2]>(std::get<_pos_neighbour[1]>(std::get<_pos_neighbour[0]>(neighbours)))) {
//                if (* item == * this) {
//                    return count;
//                }
//                count ++;
//            }
//            return count;
//        }

        lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_>
        getCurrentCoordinates()
        const
        {
            auto current_nodes_coordinates = lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_>();
            auto count = lolita::index(0);
            for (auto const & node : getComponents<_element.dim_ - 1, 0>()) {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }

        static
        lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_, matrix::col_major> const>
        getReferenceCoordinates()
        {
            using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_, matrix::col_major> const>;
            return _ReferenceCoordinates(FiniteElementConnectivity::reference_nodes_.begin()->begin());
        }

        std::basic_string<lolita::character>
        hash()
        const
        {
            std::basic_stringstream<lolita::character> hash;
            auto const & nodes = this->template getComponents<_element.dim_ - 1, 0>();
//            for (lolita::index i = 0; i < _element.num_nodes; ++i) {
//                hash << nodes[i]->hash();
//            }
            for (auto const & node : nodes) {
                hash << node->hash();
            }
            return hash.str();
        }

        Neighbours neighbours_;

        Components components_;

        lolita::natural tag_;

    };

    template<template<Element, lolita::geometry::Domain, auto...> typename _T, Element _element, lolita::geometry::Domain _domain, auto... _arg>
    requires(PointConcept<_element>)
    struct FiniteElementConnectivity<_T, _element, _domain, _arg...> : public ElementGeometry<_element>
    {

    private:

        template<Element __element, lolita::geometry::Domain __domain, auto... __arg>
        using _FiniteElementPointer = std::shared_ptr<_T<__element, __domain, __arg...>>;

    public:

        using Neighbours = typename ElementGeometry<_element>::template Neighbours<_FiniteElementPointer, _domain, _arg...>;

        FiniteElementConnectivity()
        :
                neighbours_(),
                coordinates_(),
                tag_()
        {}

        lolita::boolean
        operator==(
                FiniteElementConnectivity const &
                other
        )
        const
        {
            return this->tag_ == other.tag_;
        };

        lolita::boolean
        operator!=(
                FiniteElementConnectivity const &
                other
        )
        const
        {
            return !(* this == other);
        }

        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        lolita::geometry::Point
        getCurrentCoordinates()
        {
            return * coordinates_;
        }

        lolita::geometry::Point
        getCurrentCoordinates()
        const
        {
            return * coordinates_;
        }

        static
        lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_> const>
        getReferenceCoordinates()
        {
            using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_> const>;
            return _ReferenceCoordinates(FiniteElementConnectivity::reference_nodes.begin()->begin());
        }

        std::basic_string<lolita::character>
        hash()
        const
        {
            return std::to_string(this->tag_ + 1);
        }

        Neighbours neighbours_;

        std::shared_ptr<lolita::geometry::Point> coordinates_;

        lolita::natural tag_;

    };

//    template<Element, Basis>
    template<Element _element, lolita::finite_element::Basis _basis, lolita::index _ord>
    struct FiniteElementBasis;

//    template<Element _element, Basis _basis>
    template<Element _element, lolita::finite_element::Basis _basis, lolita::index _ord>
    requires(_basis == lolita::finite_element::Basis::Monomial)
    struct FiniteElementBasis<_element, _basis, _ord>
    {

        lolita::finite_element::Basis const static constexpr basis_ = _basis;

        lolita::index const static constexpr ord_ = _ord;

        lolita::index const static constexpr dim_ = numerics::binomial(_element.dim_ + _ord, _element.dim_);

//        template<typename T>
//        struct Implementation : public T
        template<template<Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _arg>
        struct Implementation : public FiniteElementGeometry<_T, _element, _domain, _arg...>
        {

        private:

            static constexpr
            std::array<std::array<lolita::index, 3>, dim_>
            setExponents()
            {
                auto exponents_values = std::array<std::array<lolita::index, 3>, dim_>();
                auto row = lolita::index(0);
                if constexpr (_element.dim_ == 0) {
                    exponents_values[row][0] = 0;
                    exponents_values[row][1] = 0;
                    exponents_values[row][2] = 0;
                }
                else if constexpr (_element.dim_ == 1) {
                    for (lolita::index i = 0; i < _ord + 1; ++i) {
                        exponents_values[row][0] = i;
                        exponents_values[row][1] = 0;
                        exponents_values[row][2] = 0;
                        row += 1;
                    }
                }
                else if constexpr (_element.dim_ == 2) {
                    for (lolita::index i = 0; i < _ord + 1; ++i) {
                        for (lolita::index j = 0; j < i + 1; ++j) {
                            exponents_values[row][0] = i - j;
                            exponents_values[row][1] = j;
                            exponents_values[row][2] = 0;
                            row += 1;
                        }
                    }
                }
                else if constexpr (_element.dim_ == 3) {
                    for (lolita::index i = 0; i < _ord + 1; ++i) {
                        for (lolita::index j = 0; j < i + 1; ++j) {
                            for (lolita::index k = 0; k < i + 1; ++k) {
                                if (j + k < i + 1) {
                                    exponents_values[row][0] = i - (j + k);
                                    exponents_values[row][1] = k;
                                    exponents_values[row][2] = j;
                                    row += 1;
                                }
                            }
                        }
                    }
                }
                return exponents_values;
            }

            std::array<std::array<lolita::index, 3>, dim_> const static constexpr exponents = setExponents();

        public:

            lolita::matrix::Vector<lolita::real, dim_>
            evaluate(
                    lolita::geometry::Point const &
                    point_arg
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::index i = 0; i < dim_; ++i) {
                    auto value = lolita::real(1);
                    for (lolita::index j = 0; j < _element.dim_; ++j) {
                        auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                        //value *= numerics::pow(2.0 * dist / diameters(j), exponents.get(i, j));
                        value *= std::pow(2.0 * dist / diameters(j), exponents[i][j]);
                    }
                    basis_vector_values(i) = value;
                }
                return basis_vector_values;
            }

            lolita::matrix::Vector<lolita::real, dim_>
            evaluate(
                    lolita::geometry::Point const &
                    point_arg,
                    lolita::index
                    derivative_direction_arg
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::index i = 0; i < dim_; ++i) {
                    auto value = lolita::real(1);
                    for (lolita::index j = 0; j < _element.dim_; ++j) {
                        if (j != derivative_direction_arg) {
                            auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                            value *= std::pow(2.0 * (dist) / diameters(j), exponents[i][j]);
                        }
                        else {
                            if (exponents.get(i, j) > 0) {
                                auto c = 2.0 * exponents.get(i, j) / diameters(j);
                                auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                                value *= c * std::pow(2.0 * (dist) / diameters(j), exponents[i][j] - 1);
//                                value *= c * numerics::pow(2.0 * (dist) / diameters(j), exponents.get(i, j) - 1);
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

    /*
     * BUILD
     */

    template<template<Element, lolita::geometry::Domain, auto...> typename T, Element E, lolita::geometry::Domain D, auto... A>
    struct FiniteElementGeometry : public FiniteElementConnectivity<T, E, D, A...>
    {

    private:

        using Self = element::FiniteElementGeometry<T, E, D, A...>;

    public:

        FiniteElementGeometry()
        :
        FiniteElementConnectivity<T, E, D, A...>()
        {}

//        template<Basis B>
        template<lolita::finite_element::Basis _basis, lolita::index _ord>
        lolita::matrix::Vector<lolita::real, FiniteElementBasis<E, _basis, _ord>::dim_>
        getBasisEvaluation(
                lolita::geometry::Point const &
                point_arg
        )
        const
        {
//            using FiniteElementBasis = typename EB<E, B>::template Implementation<FiniteElementGeometry>;
            using FiniteElementBasis = typename FiniteElementBasis<E, _basis, _ord>::template Implementation<T, D, A...>;
            return static_cast<FiniteElementBasis const *>(this)->evaluate(point_arg);
        }

//        template<Basis B>
        template<lolita::finite_element::Basis _basis, lolita::index _ord>
        lolita::matrix::Vector<lolita::real, FiniteElementBasis<E, _basis, _ord>::dim_>
        getBasisDerivative(
                lolita::geometry::Point const &
                point_arg,
                lolita::index
                derivative_direction_arg
        )
        const
        {
//            using FiniteElementBasis = typename EB<E, B>::template Implementation<FiniteElementGeometry>;
            using FiniteElementBasis = typename FiniteElementBasis<E, _basis, _ord>::template Implementation<T, D, A...>;
            return static_cast<FiniteElementBasis const *>(this)->evaluate(point_arg, derivative_direction_arg);
        }

        lolita::real
        getShapeMappingDifferential(
                lolita::geometry::Point const &
                point_arg
        )
        const
        {
            auto const nds = this->getCurrentCoordinates();
//            auto ru = lolita::matrix::Matrix<lolita::real, 3, E.dimPoint()>().setZero();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>().setZero();
            auto du = lolita::real(0);
            for (lolita::index i = 0; i < D.dim_; ++i) {
//                for (lolita::index j = 0; j < E.dimPoint(); ++j) {
                for (lolita::index j = 0; j < E.dim_; ++j) {
                    ru(i, j) = Self::getShapeMappingDerivative(nds.row(i), point_arg, j);
                }
            }
//            if constexpr (E.dimPoint() == 1) {
//                du = numerics::abs(ru.col(0).norm());
//            }
//            else if constexpr (E.dimPoint() == 2) {
//                du = numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
//            }
//            else if constexpr (E.dimPoint() == 3) {
//                du = numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
//            }
//            if constexpr (E.dim == 1) {
//                du = numerics::abs(ru.col(0).norm());
//            }
            if constexpr (E.dim_ == 3) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            else if constexpr (E.dim_ == 2) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else {
                du = numerics::abs(ru.col(0).norm());
            }
            if constexpr (D.frame_ == lolita::geometry::Frame::AxiSymmetric) {
                lolita::real r0 = getShapeMappingEvaluation(nds.row(0), point_arg);
                if (r0 < 1.e-10) {
                    r0 = 1.e-10;
                }
                du *= 2.0 * numerics::pi * r0;
            }
            return du;
        }

        lolita::real
        getDistanceInCurrentConfiguration(
                lolita::geometry::Point const &
                first_point_arg,
                lolita::geometry::Point const &
                second_point_arg,
                lolita::integer
                direction_arg = -1
        )
        const
        {
            assert(-1 <= direction_arg <= E.dim_);
            auto distance = lolita::real();
            if constexpr (CellConcept<E, D>) {
                distance = getCartesianDistance(first_point_arg, second_point_arg, direction_arg);
            }
            else {
                distance = getCurvedDistance(first_point_arg, second_point_arg, direction_arg);
            }
            return distance;
        }

        lolita::geometry::Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = Self::getReferenceCoordinates();
//            auto current_diameters = lolita::matrix::Vector<lolita::real, E.dim>().setZero();
            auto current_diameters = lolita::geometry::Point().setZero();
            for (lolita::index i = 0; i < E.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < E.num_nodes_; ++j) {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
//                    for (lolita::index k = 0; k < E.dimPoint(); ++k) {
                    for (lolita::index k = 0; k < 3; ++k) {
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

        lolita::geometry::Point
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            return lolita::geometry::getBarycenter(current_nodes_coordinates);
        }

        static
        lolita::real
        getReferenceDistance(
                lolita::geometry::Point const &
                first_point_arg,
                lolita::geometry::Point const &
                second_point_arg,
                lolita::integer
                direction_arg = -1
        )
        {
            assert((-1 <= direction_arg <= static_cast<lolita::integer>(E.dim_)));
            if (direction_arg == -1) {
                return (second_point_arg - first_point_arg).norm();
            } else {
                return (second_point_arg - first_point_arg)(direction_arg);
            }
        }

        static
        lolita::geometry::Point
        getReferenceCentroid()
        {
            auto nds = Self::getReferenceCoordinates();
            return lolita::geometry::getBarycenter(nds);
        }

        static
        lolita::geometry::Point
        getReferenceDiameters()
        {
            auto dts = lolita::geometry::Point().setZero();//Vector<lolita::real, E.dim>().setZero();
            auto nds = Self::getReferenceCoordinates();
            for (lolita::index i = 0; i < E.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < E.num_nodes_; ++j) {
                    for (lolita::index k = 0; k < 3; ++k) {
//                        auto & a = dts(k);
//                        auto b = numerics::abs(nds(k, i) - nds(k, j));
//                        if (b > a) {
////                            reference_diameters(k) = b;
//                            a = b;
//                        }
                        auto new_value = numerics::abs(nds(k, i) - nds(k, j));
                        auto & current_value = dts(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return dts;
        }

        /*
         * QUADRATURE
         */

        template<lolita::finite_element::Quadrature Q, lolita::index _ord>
        static
        lolita::real
        getReferenceQuadratureWeight(
                lolita::index
                index_arg
        )
        {
            return ElementQuadrature<E, Q, _ord>::reference_weights_[index_arg];
        }

        template<lolita::finite_element::Quadrature Q, lolita::index _ord>
        static
        lolita::matrix::Span<lolita::geometry::Point const>
        getReferenceQuadraturePoint(
                lolita::index
                index_arg
        )
        {
//            auto offset = E.dimPoint() * index_arg;
//            //lolita::real const constexpr * const d = ShapeQuadrature<E.shape_description, Q>::reference_points.data.data();
//            auto const constexpr * const d = ElementQuadrature<E, Q>::reference_points.data.data();
//            return matrix::MatMap<Vector<lolita::real, E.dimPoint()> const>(d + offset);
            return lolita::matrix::Span<lolita::geometry::Point const>(ElementQuadrature<E, Q, _ord>::reference_points_[index_arg].begin());
        }

        template<lolita::finite_element::Quadrature Q, lolita::index _ord, lolita::index I, lolita::index J>
//        template<Quadrature Q, lolita::index I, lolita::index J>
        static
        lolita::real
        getComponentReferenceQuadratureWeight(
                lolita::index
                index_arg
        )
        {
            using ComponentGeometry = FiniteElementGeometry<T, element::component<E, D, I, J>(), D, A...>;
            return ComponentGeometry::template getReferenceQuadratureWeight<Q, _ord>(index_arg);
        }

        template<lolita::finite_element::Quadrature Q, lolita::index _ord, lolita::index I, lolita::index J>
//        template<Quadrature Q, lolita::index I, lolita::index J>
        static
        lolita::geometry::Point
        getComponentReferenceQuadraturePoint(
                lolita::index
                component_index_arg,
                lolita::index
                index_arg
        )
        {
//            auto p = Vector<lolita::real, E.dimPoint()>();
            auto p = lolita::geometry::Point();
            using ComponentGeometry = FiniteElementGeometry<T, element::component<E, D, I, J>(), D, A...>;
//            auto const & cpt_node_tags = ElementGeometry<E>::node_connectivity.template get<I>().template get<J>();
            auto const & elt_reference_nodes = ElementGeometry<E>::reference_nodes;
            for (lolita::index i = 0; i < 3; ++i) {
                auto cpt_coordinates = lolita::matrix::Vector<lolita::real, element::component<E, D, I, J>().num_nodes_>();
                for (lolita::index j = 0; j < element::component<E, D, I, J>().num_nodes_; ++j) {
                    auto const node_tag = Self::template getComponentNodeConnection<I, J>(component_index_arg, j);//.get(component_index_arg).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<Q, _ord>(index_arg);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }

//        template<Quadrature Q>
        template<lolita::finite_element::Quadrature Q, lolita::index _ord>
        lolita::real
        getCurrentQuadratureWeight(
                lolita::index
                index_arg
        )
        const
        {
            auto w = getReferenceQuadratureWeight<Q, _ord>(index_arg);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<Q, _ord>(index_arg));
        }

        template<lolita::finite_element::Quadrature Q, lolita::index _ord>
        lolita::geometry::Point
        getCurrentQuadraturePoint(
                lolita::index
                index_arg
        )
        const
        {
            auto p = lolita::geometry::Point();
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = Self::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<Q, _ord>(index_arg));
            }
            return p;
        }

        template<lolita::finite_element::Quadrature Q, lolita::index _ord, lolita::index I, lolita::index J>
        lolita::real
        getComponentCurrentQuadratureWeight(
                lolita::index
                component_index_arg,
                lolita::index
                index_arg
        )
        {
            auto const & cmp =  this->template getComponents<I, J>()[component_index_arg];//.template get<I>().template get<J>().get(component_index_arg).get();
            return cmp->template getCurrentQuadratureWeight<Q, _ord>(index_arg);
        }

        template<lolita::finite_element::Quadrature Q, lolita::index _ord, lolita::index I, lolita::index J>
        lolita::geometry::Point
        getComponentCurrentQuadraturePoint(
                lolita::index
                component_index_arg,
                lolita::index
                index_arg
        )
        const
        {
            auto p = lolita::geometry::Point();
            auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<Q, _ord, I, J>(component_index_arg, index_arg);
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = Self::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

    private:

        lolita::real
        getCurvedDistance(
                lolita::geometry::Point const &
                first_point_arg,
                lolita::geometry::Point const &
                second_point_arg,
                lolita::integer
                direction_arg = -1
        )
        const
        requires(!PointConcept<E>)
        {
//            auto const constexpr q_seg = Quadrature(QuadratureRule::Gauss, 4);
//            auto distance = lolita::real(0);
//            auto dt = lolita::real();
//            auto const current_nodes_coordinates = this->getCurrentCoordinates();
//            for (lolita::index q = 0; q < dimQuadrature<seg_02, lolita::finite_element::Quadrature::Gauss, 4>(); ++q) {
//                auto pq = ElementQuadrature<seg_02, lolita::finite_element::Quadrature::Gauss, 4>::reference_points[q][0];
//                auto wq = ElementQuadrature<seg_02, lolita::finite_element::Quadrature::Gauss, 4>::reference_weights[q];
//                auto ru = Matrix<lolita::real, D.dim_, E.dim>();
//                auto difference = second_point_arg - first_point_arg;
//                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
//                for (lolita::index i = 0; i < D.dim_; ++i) {
//                    for (lolita::index j = 0; j < E.dim; ++j) {
//                        if (direction_arg == -1 || i == static_cast<lolita::index>(direction_arg)) {
//                            auto du = (1.0 / 2.0) * (second_point_arg(j) - first_point_arg(j));
//                            auto dx = Self::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
//                            ru(i, j) = dx * du;
//                        }
//                    }
//                }
//                if constexpr (SegmentConcept<E>) {
//                    auto Eff = ru.col(0).template dot(ru.col(0));
//                    dt = std::sqrt(Eff);
//                }
//                else if constexpr (SurfaceConcept<E>) {
//                    auto Eff = ru.col(0).template dot(ru.col(0));
//                    auto Fff = ru.col(0).template dot(ru.col(1));
//                    auto Gff = ru.col(1).template dot(ru.col(1));
//                    dt = std::sqrt(Eff + 2.0 * Fff + Gff);
//                }
//                distance += wq * dt;
//            }
//            return distance;
            using SegmentQuadrature = lolita::core::element::ElementQuadrature<seg_02, lolita::finite_element::Quadrature::Gauss, 4>;
            auto distance = lolita::real(0);
            auto dt = lolita::real();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (lolita::index q = 0; q < SegmentQuadrature::dim_; ++q) {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>().setZero();
                auto difference = second_point_arg - first_point_arg;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                for (lolita::index i = 0; i < D.dim_; ++i) {
                    for (lolita::index j = 0; j < E.dim_; ++j) {
                        if (direction_arg == -1 || i == static_cast<lolita::index>(direction_arg)) {
                            auto du = (1.0 / 2.0) * (second_point_arg(j) - first_point_arg(j));
                            auto dx = Self::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                            ru(i, j) = dx * du;
                        }
                    }
                }
                if constexpr (SegmentConcept<E>) {
                    auto Eff = ru.col(0).template dot(ru.col(0));
                    dt = std::sqrt(Eff);
                }
                else if constexpr (SurfaceConcept<E>) {
                    auto Eff = ru.col(0).template dot(ru.col(0));
                    auto Fff = ru.col(0).template dot(ru.col(1));
                    auto Gff = ru.col(1).template dot(ru.col(1));
                    dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                }
                distance += wq * dt;
            }
            return distance;
        }

        lolita::real
        getCartesianDistance(
                lolita::geometry::Point const &
                first_point_arg,
                lolita::geometry::Point const &
                second_point_arg,
                lolita::integer
                direction_arg = -1
        )
        const
        {
            auto const & nds = this->getCurrentCoordinates();
            auto distance = lolita::real();
//            auto mp0 = Vector<lolita::real, D.dim>();
//            auto mp1 = Vector<lolita::real, D.dim>();
            auto mp0 = lolita::geometry::Point();
            auto mp1 = lolita::geometry::Point();
//            for (lolita::index i = 0; i < D.dim; ++i) {
            for (lolita::index i = 0; i < 3; ++i) {
                mp0(i) = Self::getShapeMappingEvaluation(nds.row(i), first_point_arg);
                mp1(i) = Self::getShapeMappingEvaluation(nds.row(i), second_point_arg);
            }
            direction_arg == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction_arg);
            return distance;
        }

    public:

//        Vector<lolita::real, D.dim>
        lolita::geometry::Point
        getNormalVector(
                lolita::geometry::Point const &
                point_arg
        )
        const
        requires(FaceConcept<E, D>)
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto derivative_vector = lolita::matrix::Matrix<lolita::real, 3, E.dimPoint()>().setZero();
            for (lolita::index i = 0; i < 3; ++i) {
                for (lolita::index j = 0; j < E.dim_; ++j) {
                    derivative_vector(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, j);
                }
            }
//            auto const & block = derivative_vector.template block<D.dim_, E.dimPoint()>(0, 0);
//            auto normal_vector = lolita::geometry::getNormalVector(block);
            auto normal_vector = lolita::geometry::getNormalVector2(derivative_vector);
            return normal_vector;
        }

//        Vector<lolita::real, D.dim>
        lolita::geometry::Point
        getTangentVector(
                lolita::geometry::Point const &
                point_arg,
                lolita::index
                direction_arg
        )
        const
        requires(EdgeConcept<E, D>)
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = lolita::geometry::Point();
            for (lolita::index i = 0; i < 3; ++i) {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point_arg, direction_arg);
            }
            return tangent_vector;
        }

    };







    struct DegreeOfFreedom
    {

        DegreeOfFreedom()
        :
        coefficients_(),
        index_(-1)
        {}

        lolita::boolean
        isLocal()
        const
        {
            return index_ == -1 ? true : false;
        }

        lolita::boolean
        isActive()
        const
        {
            return coefficients_.size() > 0 ? true : false;
        }

        lolita::matrix::Vector<lolita::real> coefficients_;

        lolita::integer index_;

    };

    template<Element _element, lolita::geometry::Domain _domain, lolita::field::Tensor _tensor>
    struct FiniteElementUnknownF
    {

        lolita::matrix::Shape const static constexpr shp_field_ = _tensor.shape(_domain);

        FiniteElementUnknownF()
        :
        unknowns_(),
        bindings_(),
        loads_()
        {}

        void
        setUnknowns(
                lolita::index
                dim_unknowns
        )
        {
            for (int i = 0; i < shp_field_.rows_; ++i) {
                for (int j = 0; j < shp_field_.cols_; ++j) {
                    unknowns_[i][j].coefficients_.resize(dim_unknowns);
                    unknowns_[i][j].coefficients_.setZero();
                }
            }
        }

        void
        setUnknowns(
                lolita::index
                dim_unknowns,
                lolita::index &
                unknown_index
        )
        {
            for (int i = 0; i < shp_field_.rows_; ++i) {
                for (int j = 0; j < shp_field_.cols_; ++j) {
                    unknowns_[i][j].coefficients_.resize(dim_unknowns);
                    unknowns_[i][j].coefficients_.setZero();
                    unknowns_[i][j].index_ = unknown_index;
                    unknown_index += dim_unknowns;
                }
            }
        }

        void
        setUnknown(
                lolita::index
                dim_unknowns,
                lolita::index
                i,
                lolita::index
                j
        )
        {
            unknowns_[i][j].coefficients_.resize(dim_unknowns);
            unknowns_[i][j].coefficients_.setZero();
        }

        void
        setUnknown(
                lolita::index
                dim_unknowns,
                lolita::index &
                unknown_index,
                lolita::index
                i,
                lolita::index
                j
        )
        {
            unknowns_[i][j].coefficients_.resize(dim_unknowns);
            unknowns_[i][j].coefficients_.setZero();
            unknowns_[i][j].index_ = unknown_index;
            unknown_index += dim_unknowns;
        }

        void
        setBinding(
                lolita::index
                dim_unknowns,
                lolita::index &
                binding_index,
                lolita::index
                i,
                lolita::index
                j
        )
        {
            bindings_[i][j].coefficients_.resize(dim_unknowns);
            bindings_[i][j].coefficients_.setZero();
            bindings_[i][j].index_ = binding_index;
            binding_index += dim_unknowns;
        }

        void
        setLoad(
                lolita::index
                dim_unknowns,
                lolita::index &
                binding_index,
                std::shared_ptr<lolita::finite_element::LoadComponent> const &
                load_pointer,
                lolita::index
                i,
                lolita::index
                j
        )
        {
            loads_[i][j] = load_pointer;
            if (load_pointer->loading_ == lolita::finite_element::Loading::Constraint) {
                bindings_[i][j].coefficients_.resize(dim_unknowns);
                bindings_[i][j].coefficients_.setZero();
                bindings_[i][j].index_ = binding_index;
                binding_index += dim_unknowns;
            }
        }

        template<lolita::index _size>
        auto
        getUnknownCoefficients(
                lolita::index
                i,
                lolita::index
                j
        )
        const
        {
            return unknowns_[i][j].coefficients_.template segment<_size>(0);
        };

        template<lolita::index _size>
        auto
        getBindingCoefficients(
                lolita::index
                i,
                lolita::index
                j
        )
        const
        {
            return bindings_[i][j].coefficients_.template segment<_size>(0);
        };

        lolita::boolean
        isActive(
                lolita::index
                i,
                lolita::index
                j
        )
        const
        {
            return unknowns_[i][j].coefficients_.size() > 0 ? true : false;
        }

        lolita::boolean
        isBound(
                lolita::index
                i,
                lolita::index
                j
        )
        const
        {
            return bindings_[i][j].coefficients_.size() > 0 ? true : false;
        }

        std::array<std::array<DegreeOfFreedom, shp_field_.rows_>, shp_field_.cols_> unknowns_;

        std::array<std::array<DegreeOfFreedom, shp_field_.rows_>, shp_field_.cols_> bindings_;

        std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, shp_field_.rows_>, shp_field_.cols_> loads_;

    };

    namespace utility
    {

        struct MappingPosition
        {

            lolita::index i_;

            lolita::index j_;

        };

        template<lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        static constexpr
        lolita::index
        getDimField()
        {
            return _finite_element.unknown_.tensor_.shape(_domain).size_;
        }

        template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        static constexpr
        lolita::index
        getNumUnknowns()
        {
            return FiniteElementPolicy<_element, _domain, _finite_element>::dimUnknowns() * getDimField<_domain, _finite_element>();
        }

        template<lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element, lolita::field::Mapping _mapping>
        static constexpr
        lolita::index
        getDimMapping()
        {
            return lolita::field::MappingPolicy<_finite_element.unknown_.tensor_, _domain, _mapping>::cardinality_.size_;
        }

        template<lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element, lolita::field::Mapping _mapping>
        static constexpr
        lolita::index
        getRowMapping()
        {
            auto value = lolita::index(0);
            for (auto m : _finite_element.unknown_.mappings_) {
                if (m == _mapping) {
                    return value;
                }
                value += getDimMapping<_domain, _finite_element, _mapping>();
            }
            return value;
        }

        template<lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        static constexpr
        lolita::index
        getRowsOperator()
        {
            auto dim_mapping = lolita::index(0);
            auto set_dim_mapping = [&] <lolita::index _i = 0> (auto & self)
                    constexpr mutable
            {
                dim_mapping += getDimMapping<_domain, _finite_element, _finite_element.unknown_.mappings_[_i]>();
                if constexpr (_i < _finite_element.unknown_.mappings_.size() - 1) {
                    self.template operator()<_i + 1>(self);
                }
            };
            set_dim_mapping(set_dim_mapping);
            return dim_mapping;
        }

    }

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(CellConcept<_element, _domain> && _finite_element.finite_element_method_.method_ == lolita::finite_element::Method::HHO)
    struct FiniteElementPolicy<_element, _domain, _finite_element>
    {

    private:

        using Base = FiniteElementModule<_element, _domain, _finite_element>;

    public:

        static constexpr
        lolita::index
        dimUnknowns()
        {
            return FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.finite_element_method_.ord_cell_>::dim_;
        }

        static constexpr
        lolita::index
        getColsOperator()
        {
            using _CellBasis = FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.finite_element_method_.ord_cell_>;
            auto dim_field_mapping = _CellBasis::dim_ * utility::getDimField<_domain, _finite_element>();
            auto set_dim_field_mapping = [&] <lolita::index K = 0> (auto & self) constexpr mutable {
                auto const constexpr _component = element::component<_element, _domain, 0, K>();
                using _FaceBasis = FiniteElementBasis<_component, lolita::finite_element::Basis::Monomial, _finite_element.finite_element_method_.ord_face_>;
                dim_field_mapping += _FaceBasis::dim_ * utility::getDimField<_domain, _finite_element>();
                if constexpr (K < element::numComponents<_element, _domain, 0>() - 1) {
                    self.template operator ()<K + 1>(self);
                }
            };
            set_dim_field_mapping(set_dim_field_mapping);
            return dim_field_mapping;
        }

        struct Module
        {

//            using CellBasis = FiniteElementBasis<E, lolita::finite_element::Basis::Monomial, F.finite_element_method_.ord_cell_>;

//            static constexpr
//            lolita::index
//            dimCellUnknowns()
//            {
//                return FiniteElementBasis<E, lolita::finite_element::Basis::Monomial, F.finite_element_method_.ord_cell_>::dim_;
//            }

            static constexpr
            lolita::index
            numCellUnknowns()
            {
                using _CellBasis = FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.finite_element_method_.ord_cell_>;
                return _CellBasis::dim_ * _finite_element.unknown_.tensor_.shape(_domain).size_;
            }

            static constexpr
            lolita::index
            dimOperator()
            {
                using _CellBasis = FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.finite_element_method_.ord_cell_>;
                auto dim_field_mapping = _CellBasis::dim_ * utility::getDimField<_domain, _finite_element>();
                auto set_dim_field_mapping = [&] <lolita::index K = 0> (
                        auto & self
                )
                constexpr mutable
                {
                    auto const constexpr _component = element::component<_element, _domain, 0, K>();
                    using _FaceBasis = FiniteElementBasis<_component, lolita::finite_element::Basis::Monomial, _finite_element.finite_element_method_.ord_face_>;
                    dim_field_mapping += _FaceBasis::dim_ * utility::getDimField<_domain, _finite_element>();
                    if constexpr (K < element::numComponents<_element, _domain, 0>() - 1) {
                        self.template operator ()<K + 1>(self);
                    }
                };
                set_dim_field_mapping(set_dim_field_mapping);
                return dim_field_mapping;
            }

            lolita::matrix::Matrix<lolita::real, dimOperator(), dimOperator()> stabilization_;

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

            void
            setModule()
            {
                this->module_.stabilization_.setOnes();
                this->module_.stabilization_(0, 0) = 3;
            }

            template<lolita::field::Mapping Map, auto I>
            void
            setOperator();

            template<lolita::field::Mapping Map, auto I>
            void
            setOperator()
            requires(Map == lolita::field::Mapping::Gradient)
            {
                /*
                 * Defining constants
                 */
                auto const constexpr ord_cell = _finite_element.finite_element_method_.ord_cell_;
                auto const constexpr ord_face = _finite_element.finite_element_method_.ord_face_;
                auto const constexpr ord_grad = _finite_element.finite_element_method_.ordMapping(Map);
                auto const constexpr ord_max = lolita::numerics::max(ord_cell, ord_face, ord_grad);
                auto const constexpr ord_qad = _domain.ordIntegration(ord_max);
//                auto const constexpr dim_qad = dimQuadrature<E, lolita::finite_element::Quadrature::Gauss, ord_qad>();
//                auto const constexpr dim_bas_grad = dimBasis<E, BasisName::Monomial, ord_grad>();
//                auto const constexpr dim_bas_cell = dimBasis<E, BasisName::Monomial, ord_cell>();
//                auto const constexpr dim_bas_face = dimBasis<E, BasisName::Monomial, ord_face>();
                using _Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, ord_qad>;
                using _Quadrature2 = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _finite_element.ord_quadrature_>;
                using _CellBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_cell>;
                using _FaceBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_face>;
                using _GradBasis = lolita::core::element::FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, ord_grad>;
                /*
                 *
                 */
                using _MapPol = lolita::field::MappingPolicy<_finite_element.unknown_.tensor_, _domain, Map>;
                /*
                 * LHS
                 */
                auto get_lhs = [&] ()
                        mutable
                {
                    auto lhs = lolita::matrix::Matrix<lolita::real, _GradBasis::dim_, _GradBasis::dim_>().setZero();
                    for (int i = 0; i < _Quadrature::dim_quadrature; ++i) {
                        auto pt = this->template getReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                        auto wt = this->template getCurrentQuadratureWeight<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                        auto vr = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pt);
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
                    auto rhs = lolita::matrix::Matrix<lolita::real, _GradBasis::dim_, getColsOperator()>().setZero();
                    /*
                     * Setting the cell part of the gradient operator
                     */
                    auto set_rhs_cell = [&] ()
                            mutable
                    {
                        for (int i = 0; i < _Quadrature::dim_quadrature; ++i) {
                            auto pc = this->template getReferenceQuadraturePoint<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                            auto wc = this->template getCurrentQuadratureWeight<_Quadrature::quadrature_, _Quadrature::ord_>(i);
                            auto vr = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pc);
                            auto vc = this->template getBasisDerivative<_CellBasis::basis_, _CellBasis::ord_>(pc, j_arg);
                            auto prt = rhs.template block<_GradBasis::dim_, _CellBasis::dim_>(0, _CellBasis::dim_ * i_arg);
                            prt += wc * vr * vc.transpose();
                        }
                    };
                    /*
                     * Defining the face offset
                     */
                    auto get_faces_offset2 = [] <auto K> ()
                            constexpr
                    {
                        auto offset = _CellBasis::dim_ * utility::getDimField<_domain, _finite_element>();
                        auto get_faces_offset22 = [&] <auto L = 0> (auto & self)
                                constexpr mutable
                        {
                            if constexpr (L > 0) {
                                auto nj = FiniteElementPolicy<lolita::core::element::component<_element, _domain, 0, L - 1>(), _domain, _finite_element>::numUnknowns();
                                auto num_c = numComponents<_element, _domain, 0, L - 1>();
                                offset += num_c * nj;
                            }
                            if constexpr (L < K) {
                                self.template operator()<K + 1>(self);
                            }
                        };
                        get_faces_offset22(get_faces_offset22);
                        return offset;
                    };
                    /*
                     * Setting the jump part of the gradient operator, at a given face
                     */
//                    auto set_rhs_face = [&] <auto K = 0> (
//                            auto &
//                            self
//                    )
//                    mutable
//                    {
//                        auto const constexpr dim_bas_fce = dimBasis<component<E, 0, K>(), bas_fce>();
//                        for (int num_f = 0; num_f < numComponents<E, 0, K>(); ++num_f) {
//                            auto const & face = this->template getComponentPointer<0, K>(num_f).ptr.get();
//                            auto const & n_dir = this->template getComponentPointer<0, K>(num_f).orientation();
//                            for (int i = 0; i < dim_qad_grd; ++i) {
//                                auto pf = face.template getReferenceQuadraturePoint<qad_grd>(i);
//                                auto pc = this->template getComponentReferenceQuadraturePoint<qad_grd, 0, K>(num_f, i);
//                                auto wf = face.template getCurrentQuadratureWeight<qad_grd>(i);
//                                auto n = face.getNormalVector(pf);
//                                auto vr = this->template getBasisEvaluation<bas_grd>(pc);
//                                auto vc = this->template getBasisEvaluation<bas_cel>(pc);
//                                auto vf = face.template getBasisEvaluation<bas_fce>(pf);
//                                auto offset = dim_bas_cel * i_arg;
//                                auto prt_cell = rhs.template block<dim_bas_grd, dim_bas_cel>(0, offset);
//                                prt_cell -= wf * vr * vc.transpose() * n(j_arg) * n_dir;
//                                //offset = get_faces_offset.template operator()<K>();
//                                offset = get_faces_offset2.template operator()<K>();
//                                offset += dim_bas_fce * (num_f * field<D, F>().size() + i_arg);
//                                auto prt_face = rhs.template block<dim_bas_grd, dim_bas_fce>(0, offset);
//                                prt_face += wf * vr * vf.transpose() * n(j_arg) * n_dir;
//                            }
//                        }
//                        if constexpr (K < numComponents<E, 0>() - 1) {
//                            self.template operator()<K + 1>(i_arg, j_arg, self);
//                        }
//                    };
//                    set_rhs_face(set_rhs_face);
                    return rhs;
                };
                /*
                 *
                 */
//                auto const constexpr fld = Field(F.field().ord_field, D.dim);
//                auto const constexpr grd = Field::fromMapping(Base::field, Map);
                auto row = lolita::index(0);
                auto lhs = get_lhs();
                for (lolita::field::MappingValues const & item : _MapPol::cardinality_) {
                    auto rhs = get_rhs(item.col_, item.row_);
                    for (int k = 0; k < _Quadrature::dim; ++k) {
                        auto pnt = this->template getReferenceQuadraturePoint<_Quadrature2::quadrature_, _Quadrature2::ord_>(k);
                        auto vct = this->template getBasisEvaluation<_GradBasis::basis_, _GradBasis::ord_>(pnt);
                        this->operators_[k].row(this->rowMapping(Map) + row) = vct.transpose() * lhs * rhs;
                    }
                    row += 1;
                }
//                for (int i = 0; i < _MapPol::cardinality_.cols_; ++i) {
//                    for (int j = 0; j < _MapPol::cardinality_.rows_; ++j) {
//                        auto rhs = get_rhs(i, j);
//                        for (int k = 0; k < dimQuadrature<E, this->quadrature>(); ++k) {
//                            auto pnt = this->template getReferenceQuadraturePoint<this->quadrature>(k);
//                            auto vct = this->template getBasisEvaluation<bas_grd>(pnt);
//                            this->operators.get(k).row(this->rowOperator(Map) + row) = vct.transpose() * lhs * rhs;
//                        }
//                        row += 1;
//                    }
//                }
            }

        };

    };

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(FaceConcept<_element, _domain> && _finite_element.finite_element_method_.method_ == lolita::finite_element::Method::HHO)
    struct FiniteElementPolicy<_element, _domain, _finite_element>
    {

        static constexpr
        lolita::index
        dimUnknowns()
        {
            return FiniteElementBasis<_element, lolita::finite_element::Basis::Monomial, _finite_element.finite_element_method_.ord_face_>::dim_;
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

        };

    };

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(EdgeConcept<_element, _domain> && _finite_element.finite_element_method_.method_ == lolita::finite_element::Method::HHO)
    struct FiniteElementPolicy<_element, _domain, _finite_element>
    {

        static constexpr
        lolita::index
        dimUnknowns()
        {
            return 0;
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

        };

    };

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(NodeConcept<_element, _domain> && _finite_element.finite_element_method_.method_ == lolita::finite_element::Method::HHO)
    struct FiniteElementPolicy<_element, _domain, _finite_element>
    {

        static constexpr
        lolita::index
        dimUnknowns()
        {
            return 0;
        }

        struct Module
        {

        };

        struct Implementation : public FiniteElementModule<_element, _domain, _finite_element>
        {

        };

    };

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(CellConcept<_element, _domain>)
    struct FiniteElementModule<_element, _domain, _finite_element> : public FiniteElementUnknownF<_element, _domain, _finite_element.unknown_.tensor_>, public FiniteElementGeometry<FiniteElementModule, _element, _domain, _finite_element>
    {

    private:

        using _Module = typename FiniteElementPolicy<_element, _domain, _finite_element>::Module;

        using _Policy = FiniteElementPolicy<_element, _domain, _finite_element>;

    public:

        using Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _finite_element.ord_quadrature_>;

        using Operator = lolita::matrix::Matrix<lolita::real, utility::getRowsOperator<_domain, _finite_element>(), _Policy::getColsOperator()>;

//        using Operators = std::array<Operator, Quadrature::dim_>;
//
//        using MaterialPoints = std::array<std::shared_ptr<mgis::behaviour::BehaviourData>, Quadrature::dim_>;

        FiniteElementModule()
        :
        FiniteElementUnknownF<_element, _domain, _finite_element.unknown_.tensor_>(),
        FiniteElementGeometry<FiniteElementModule, _element, _domain, _finite_element>()
        {}

        void
        setUnknowns(
                lolita::index &
                degree_of_freedom_index_arg
        )
        {
            reinterpret_cast<FiniteElementUnknownF<_element, _domain, _finite_element> *>(this)->setUnknowns(FiniteElementPolicy<_element, _domain, _finite_element>::numUnknowns());
        }

        void
        setMaterial(
                mgis::behaviour::Behaviour const &
                behaviour_arg
        )
        {
            for (int i = 0; i < Quadrature::dim_quadrature; ++i) {
                material_points_[i] = std::make_shared<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(behaviour_arg));
            }
        }

        void
        initialize()
        {
            using Implementation = typename FiniteElementPolicy<_element, _domain, _finite_element>::Implementation;
            static_cast<Implementation *>(this)->setModule();
            auto set_operator = [&] <lolita::index I = 0> (auto & self)
            mutable
            {
                static_cast<Implementation *>(this)->template setOperator<_finite_element.unknown_.mappings_[I], I>();
                if constexpr (I < _finite_element.unknown_.mappings_.size() - 1) {
                    self.template operator()<I + 1>(self);
                }
            };
            set_operator(set_operator);
        }

        _Module module_;

        std::array<Operator, Quadrature::dim_> operators_;

        std::array<std::shared_ptr<mgis::behaviour::BehaviourData>, Quadrature::dim_> material_points_;

    };

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(FaceConcept<_element, _domain>)
    struct FiniteElementModule<_element, _domain, _finite_element> : public FiniteElementUnknownF<_element, _domain, _finite_element.unknown_.tensor_>, public FiniteElementGeometry<FiniteElementModule, _element, _domain, _finite_element>
    {

    private:

        using _Module = typename FiniteElementPolicy<_element, _domain, _finite_element>::Module;

    public:

        FiniteElementModule()
        :
        FiniteElementUnknownF<_element, _domain, _finite_element.unknown_.tensor_>(),
        FiniteElementGeometry<FiniteElementModule, _element, _domain, _finite_element>()
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
            reinterpret_cast<FiniteElementUnknownF<_element, _domain, _finite_element> *>(this)->setUnknowns(degree_of_freedom_index_arg);
        }

    };

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(EdgeConcept<_element, _domain>)
    struct FiniteElementModule<_element, _domain, _finite_element> : public FiniteElementUnknownF<_element, _domain, _finite_element.unknown_.tensor_>, public FiniteElementGeometry<FiniteElementModule, _element, _domain, _finite_element>
    {

    private:

        using _Module = typename FiniteElementPolicy<_element, _domain, _finite_element>::Module;

    public:

        FiniteElementModule()
        :
        FiniteElementUnknownF<_element, _domain, _finite_element.unknown_.tensor_>(),
        FiniteElementGeometry<FiniteElementModule, _element, _domain, _finite_element>()
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

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    requires(NodeConcept<_element, _domain>)
    struct FiniteElementModule<_element, _domain, _finite_element> : public FiniteElementUnknownF<_element, _domain, _finite_element.unknown_.tensor_>, public FiniteElementGeometry<FiniteElementModule, _element, _domain, _finite_element>
    {

    private:

        using _Module = typename FiniteElementPolicy<_element, _domain, _finite_element>::Module;

    public:

        FiniteElementModule()
        :
        FiniteElementUnknownF<_element, _domain, _finite_element.unknown_.tensor_>(),
        FiniteElementGeometry<FiniteElementModule, _element, _domain, _finite_element>()
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

    /*
     *
     */

    template<Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    struct FiniteElementFinal : public FiniteElementGeometry<FiniteElementFinal, _element, _domain, _finite_element...>
    {

    private:

        using _ElementsPtr = std::tuple<std::shared_ptr<FiniteElementModule<_element, _domain, _finite_element>>...>;

        using _Elements = std::tuple<FiniteElementModule<_element, _domain, _finite_element>...>;

    public:

        FiniteElementFinal()
        :
        FiniteElementGeometry<FiniteElementFinal, _element, _domain, _finite_element...>(),
        elements_()
        {}

        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementsPtr> const &
        getElement()
        const
        {
            return std::get<_i>(elements_);
        }

        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementsPtr> &
        getElement()
        {
            return std::get<_i>(elements_);
        }

        template<lolita::index _i>
        void
        makeElement()
        {
            getElement<_i>() = std::make_shared<std::tuple_element_t<_i, _Elements>>(std::tuple_element_t<_i, _Elements>());
        }

        void
        init()
        {
            auto initialize_element = [&] <lolita::index _i = 0> (auto & self)
            constexpr mutable
            {
                using _Element = typename std::tuple_element_t<_i, _ElementsPtr>::element_type;
                getElement<_i>() = std::make_shared<_Element>(_Element());
                getElement<_i>()->initialize();
            };
            initialize_element(initialize_element);
        }

        _ElementsPtr elements_;

    };
//
//    /*
//     * Unknowns
//     */
//
//    template<Element E, lolita::geometry::Domain D, auto C>
//    using MixedElement = typename decltype(C)::template MixedElementType<FiniteElementUnknownMod, E, D>;
//
//    template<Element E, lolita::geometry::Domain D, auto C>
//    struct FiniteElement : public MixedElement<E, D, C>, public FiniteElementGeometry<E, D, FiniteElement, C>
//    {
//
//        FiniteElement()
//        :
//        MixedElement<E, D, C>(),
//        FiniteElementGeometry<E, D, FiniteElement, C>()
//        {}
//
//        void
//        initialize()
//        {
//
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
//        auto
//        hash()
//        const
//        {
//            StrgStream hash;
//            if constexpr(Point<E>) {
//                hash << tag + 1;
//            }
//            else {
//                auto const & nodes = this->components.template get<E.dim - 1>().template get<0>();
//                for (lolita::index i = 0; i < E.num_nodes; ++i) {
//                    hash << nodes.get(i).ptr.get().hash();
//                }
//            }
//            return hash.str();
//        }
//
//        Long tag;
//
//    };

//    template<template<Element, lolita::geometry::Domain, auto...> typename T, Element E, lolita::geometry::Domain D, auto... A>
//    requires(PointConcept<E>)
//    struct FiniteElementConnectivity<T, E, D, A...> : public ElementGeometry<E>
//    {
//
//    private:
//
//        using Self = FiniteElementConnectivity;
//
//        template<Element Eb>
//        using FiniteElementPointer = SharedPointer<T<Eb, D, A...>>;
//
//        template<Element Eb>
//        using Neighbourggg = FiniteElementNeighbour<Eb, D, T, A...>;
//
//    public:
//
//        using Coordinates = SharedPointer<Vector<lolita::real, D.dim>>;
//
//        template<Element Eb>
//        struct Neighbour
//        {
//
//            Neighbour()
//            :
//            ptr()
//            {}
//
//            Neighbour(
//                    SharedPointer<T<Eb, D, A...>> const &
//                    ptr_arg
//            )
//            :
//            ptr(ptr_arg)
//            {}
//
//            Bool
//            operator==(
//                    Neighbour const &
//                    other
//            )
//            const = default;
//
//            Bool
//            operator!=(
//                    Neighbour const &
//                    other
//            )
//            const = default;
//
//            SharedPointer<T<Eb, D, A...>> ptr;
//
//        };
//
//        //using Neighbours = ElementNeighbourArray<E, D.dim, Neighbour>;
//
//        using Neighbours = ElementNeighbourArray<E, D.dim, Neighbourggg>;
//
//        FiniteElementConnectivity()
//        :
//        coordinates(),
//        neighbours(),
//        count(0)
//        {}
//
//        template<auto I, auto J>
//        auto &
//        getNeighbour(
//                auto
//                index_arg
//        )
//        {
//            return neighbours.template get<I>().template get<J>().get(index_arg).get();
//        }
//
//        template<auto I, auto J>
//        auto const &
//        getNeighbour(
//                auto
//                index_arg
//        )
//        const
//        {
//            return neighbours.template get<I>().template get<J>().get(index_arg).get();
//        }
//
//        auto
//        getCurrentCoordinates()
//        {
//            //return coordinates;
//            return coordinates.get();
//        }
//
//        auto
//        getCurrentCoordinates()
//        const
//        {
//            //return coordinates;
//            return coordinates.get();
//        }
//
//        static
//        auto
//        getReferenceCoordinates()
//        {
//            return Vector<lolita::real, 1>{0.0};
//        }
//
//        Coordinates coordinates;
//
//        Neighbours neighbours;
//
//        lolita::index count;
//
//    };

    /*
     * BASIS
     */

}

#endif //LOLITA_LOLITA_CORE_FINITE_ELEMENT_HXX
