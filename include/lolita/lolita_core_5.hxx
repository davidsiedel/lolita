//
// Created by dsiedel on 04/06/22.
//

#ifndef LOLITA_LOLITA_CORE_5_HXX
#define LOLITA_LOLITA_CORE_5_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_4.hxx"

namespace lolita::core2::finite_element
{

    /**
     * @brief
     */
    auto const static null_load_ptr = std::make_shared<lolita::finite_element::ScalarLoad>(lolita::finite_element::ScalarLoad());

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementTraits;

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, auto t_finite_element>
    struct FiniteElement;

    /**
     * @brief
     * @tparam t_T
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, auto t_finite_element>
    struct FiniteElementConnectivity
    {

    private:

        /**
         * @brief
         */
        using t_ElementTraits = lolita::core2::geometry::ElementTraits<t_element, t_domain>;

        /**
         * @brief
         */
        template<lolita::core2::geometry::Element t__element, lolita::domain::Domain t__domain, auto... t__finite_element>
        using t_ElementPointer = std::shared_ptr<FiniteElement<t__element, t__domain, t__finite_element...>>;

    public:

        /**
         * @brief
         */
        using Components = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer, t_finite_element>;

        /**
         * @brief
         */
        using Neighbours = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer, t_finite_element>;

        /**
         * @brief
         */
        Neighbours neighbours_;

        /**
         * @brief
         */
        Components components_;

        /**
         * @brief
         */
        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

        /**
         * @brief
         */
        lolita::natural tag_;

        /**
         * @brief
         */
        std::shared_ptr<lolita::domain::Point> coordinates_;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        lolita::boolean
        operator==(
                FiniteElementConnectivity const & other
        )
        const = default;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        lolita::boolean
        operator!=(
                FiniteElementConnectivity const & other
        )
        const = default;

        /**
         * @brief
         * @return
         */
        std::basic_string<lolita::character>
        hash()
        const
        requires(t_element.isPoint())
        {
            return std::to_string(this->tag_);
        }

        /**
         * @brief
         * @return
         */
        std::basic_string<lolita::character>
        hash()
        const
        {
            std::basic_stringstream<lolita::character> hash;
            auto const & nodes = getComponents<t_element.dim_ - 1, 0>();
            for (auto const & node : nodes) {
                hash << node->hash();
            }
            return hash.str();
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @param i
         * @param j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        static constexpr
        lolita::index
        getComponentNodeConnection(
                lolita::index i,
                lolita::index j
        )
        requires(!t_element.isPoint())
        {
            return std::get<t_j>(std::get<t_i>(lolita::core2::geometry::ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> &
        getComponents()
        requires(!t_element.isPoint())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> const &
        getComponents()
        const
        requires(!t_element.isPoint())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @param i
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        lolita::index
        getComponentIndex(
                lolita::index i
        )
        const
        requires(!t_element.isPoint())
        {
            auto constexpr t_component = t_ElementTraits::template getComponent<t_i, t_j>();
            using t_NeighbourTraits = lolita::core2::geometry::ElementTraits<t_component, t_domain>;
            auto constexpr t_coordinates = t_NeighbourTraits::template getNeighbourCoordinates<t_element>();
            auto const & items = getComponents<t_i, t_j>()[i]->template getNeighbours<t_coordinates.dim_, t_coordinates.tag_>();
            auto is_equal = [&] (std::shared_ptr<FiniteElementConnectivity> const & ptr_element) {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @param i
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        lolita::integer
        getComponentOrientation(
                lolita::index i
        )
        const
        requires(!t_element.isPoint())
        {
            return getComponentIndex<t_i, t_j>(i) == 0 ? 1 : -1;
        }

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementGeometry : virtual FiniteElementConnectivity<t_element, t_domain, t_finite_element>
    {

    private:

        /**
         * @brief
         */
        using t_ElementTraits = lolita::core2::geometry::ElementTraits<t_element, t_domain>;

        /**
         * @brief
         */
        using t_FiniteElementTraits = lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         */
        using t_Quadrature = typename t_FiniteElementTraits::Quadrature;

    public:

        /**
         * @brief Get the current coordinates of the point/node in the mesh
         * @return the vector with x, y and z coordinates
         */
        lolita::domain::Point &
        getCurrentCoordinates()
        requires(t_element.isPoint())
        {
            return * this->coordinates_;
        }

        /**
         * @brief Get the current coordinates of the point/node in the mesh
         * @return the vector with x, y and z coordinates
         */
        lolita::domain::Point const &
        getCurrentCoordinates()
        const
        requires(t_element.isPoint())
        {
            return * this->coordinates_;
        }

        /**
         * @brief Get the current coordinates of the element in the mesh
         * @return the matrix composed of the column vectors with x, y and z coordinates of each node
         */
        lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>
        getCurrentCoordinates()
        const
//        requires(!t_element.isPoint())
        {
            auto current_nodes_coordinates = lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>();
            auto count = lolita::index(0);
            for (auto const & node : this->template getComponents<t_element.dim_ - 1, 0>()) {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }

        /**
         * @brief
         * @return
         */
        static
        lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.numNodes(), matrix::col_major> const>
        getReferenceCoordinates()
        {
            using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.numNodes(), matrix::col_major> const>;
            return _ReferenceCoordinates(lolita::core2::geometry::ElementTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
        }

//        /**
//         * @brief
//         * @tparam _basis
//         * @tparam _ord
//         * @param point
//         * @return
//         */
//        template<lolita::core2::finite_element::basis::Basis _basis, lolita::index _ord>
//        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, _basis, _ord>::dim_>
//        getBasisEvaluation(
//                lolita::domain::Point const & point
//        )
//        const
//        {
//            using t_FiniteElementBasis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, _basis, _ord>;
//            using t_FiniteElementBasisImplementation = typename t_FiniteElementBasis::template Implementation<t_domain, t_finite_element>;
//            return static_cast<t_FiniteElementBasisImplementation const *>(this)->evaluate(point);
//        }
//
//        /**
//         * @brief
//         * @tparam _basis
//         * @tparam _ord
//         * @param point
//         * @param derivative_direction
//         * @return
//         */
//        template<lolita::core2::finite_element::basis::Basis _basis, lolita::index _ord>
//        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, _basis, _ord>::dim_>
//        getBasisDerivative(
//                lolita::domain::Point const & point,
//                lolita::index derivative_direction
//        )
//        const
//        {
//            using t_FiniteElementBasis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, _basis, _ord>;
//            using t_FiniteElementBasisImplementation = typename t_FiniteElementBasis::template Implementation<t_domain, t_finite_element>;
//            return static_cast<t_FiniteElementBasisImplementation const *>(this)->evaluate(point, derivative_direction);
//        }

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, t_element.numNodes()> const & nodal_field_values,
                lolita::domain::Point const & reference_point
        )
        {
            return t_ElementTraits::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @param derivative_direction
         * @return
         */
        static
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, t_element.numNodes()> const & nodal_field_values,
                lolita::domain::Point const & reference_point,
                lolita::index derivative_direction
        )
        {
            return t_ElementTraits::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }

        /**
         * @brief
         * @param point
         * @return
         */
        lolita::real
        getShapeMappingDifferential(
                lolita::domain::Point const & point
        )
        const
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>();
            auto du = lolita::real(0);
            ru.setZero();
            for (lolita::index i = 0; i < t_domain.dim_; ++i) {
                for (lolita::index j = 0; j < t_element.dim_; ++j) {
                    ru(i, j) = FiniteElementGeometry::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.dim_ == 3) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            else if constexpr (t_element.dim_ == 2) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else {
                du = numerics::abs(ru.col(0).norm());
            }
            if constexpr (t_domain.frame_ == lolita::domain::Frame::AxiSymmetric()) {
                lolita::real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10) {
                    r0 = 1.e-10;
                }
                du *= 2.0 * numerics::pi * r0;
            }
            return du;
        }

        /**
         * @brief
         * @param first_point
         * @param second_point
         * @param direction
         * @return
         */
        lolita::real
        getRiemannianDistance(
                lolita::domain::Point const & first_point,
                lolita::domain::Point const & second_point,
                lolita::integer direction = -1
        )
        const
        {
            if constexpr (t_ElementTraits::hasDim(0)) {
                auto const & current_coordinates = this->getCurrentCoordinates();
                auto distance = lolita::real();
                auto mp0 = lolita::domain::Point();
                auto mp1 = lolita::domain::Point();
                for (lolita::index i = 0; i < t_element.dim_; ++i) {
                    mp0(i) = FiniteElementGeometry::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                    mp1(i) = FiniteElementGeometry::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
                }
                direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
                return distance;
            }
            else {
                auto constexpr _quadrature = lolita::finite_element::Quadrature::Gauss();
                using SegmentQuadrature = lolita::core2::geometry::ElementQuadratureTraits<lolita::core2::geometry::Element::LinearSegment(), _quadrature, 4>;
                auto distance = lolita::real(0);
                auto dt = lolita::real();
                auto const current_nodes_coordinates = this->getCurrentCoordinates();
                for (lolita::index q = 0; q < SegmentQuadrature::getDim(); ++q) {
                    auto pq = SegmentQuadrature::reference_points_[q][0];
                    auto wq = SegmentQuadrature::reference_weights_[q];
                    auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>().setZero();
                    auto difference = second_point - first_point;
                    auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                    for (lolita::index i = 0; i < t_domain.dim_; ++i) {
                        for (lolita::index j = 0; j < t_element.dim_; ++j) {
                            if (direction == -1 || i == static_cast<lolita::index>(direction)) {
                                auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                                auto dx = FiniteElementGeometry::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                                ru(i, j) = dx * du;
                            }
                        }
                    }
                    if constexpr (t_element.isFacet()) {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        auto Fff = ru.col(0).template dot(ru.col(1));
                        auto Gff = ru.col(1).template dot(ru.col(1));
                        dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                    }
                    else if constexpr (t_element.isCurve()) {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        dt = std::sqrt(Eff);
                    }
                    else {
                        dt = 0;
                    }
                    distance += wq * dt;
                }
                return distance;
            }
        }

        /*
         * GEOMETRY
         */

        /**
         * @brief
         * @return
         */
        lolita::domain::Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElementGeometry::getReferenceCoordinates();
            auto current_diameters = lolita::domain::Point().setZero();
            for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (lolita::index k = 0; k < 3; ++k) {
                        auto new_value = numerics::abs(getRiemannianDistance(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }

        /**
         * @brief
         * @return
         */
        lolita::domain::Point
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            return lolita::domain::getBarycenter(current_nodes_coordinates);
        }

        /**
         * @brief
         * @return
         */
        static
        lolita::domain::Point
        getReferenceCentroid()
        {
            auto nds = FiniteElementGeometry::getReferenceCoordinates();
            return lolita::domain::getBarycenter(nds);
        }

        /**
         * @brief
         * @return
         */
        static
        lolita::domain::Point
        getReferenceDiameters()
        {
            auto dts = lolita::domain::Point();
            auto nds = FiniteElementGeometry::getReferenceCoordinates();
            dts.setZero();
            for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
                    for (lolita::index k = 0; k < 3; ++k) {
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

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
        static
        lolita::real
        getReferenceQuadratureWeight(
                lolita::index index
        )
        {
            return lolita::core2::geometry::ElementQuadratureTraits<t_element, _quadrature, _ord>::reference_weights_[index];
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
        static
        lolita::matrix::Span<lolita::domain::Point const>
        getReferenceQuadraturePoint(
                lolita::index index
        )
        {
            return lolita::matrix::Span<lolita::domain::Point const>(
                    lolita::core2::geometry::ElementQuadratureTraits<t_element, _quadrature, _ord>::reference_points_[index].begin()
            );
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
        lolita::real
        getCurrentQuadratureWeight(
                lolita::index index
        )
        const
        {
            auto w = getReferenceQuadratureWeight<_quadrature, _ord>(index);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<_quadrature, _ord>(index));
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
        lolita::domain::Point
        getCurrentQuadraturePoint(
                lolita::index index
        )
        const
        {
            auto p = lolita::domain::Point();
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<_quadrature, _ord>(index));
            }
            return p;
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @tparam _i
         * @tparam _j
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
        static
        lolita::real
        getComponentReferenceQuadratureWeight(
                lolita::index index
        )
        requires(!t_element.isPoint())
        {
            auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
            using ComponentGeometry = FiniteElementGeometry<_component, t_domain, t_finite_element>;
            return ComponentGeometry::template getReferenceQuadratureWeight<_quadrature, _ord>(index);
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @tparam _i
         * @tparam _j
         * @param component_index
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
        static
        lolita::domain::Point
        getComponentReferenceQuadraturePoint(
                lolita::index component_index,
                lolita::index index
        )
        requires(!t_element.isPoint())
        {
//            auto const constexpr _component = lolita::core::component<_element, _domain, _i, _j>();
            auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
            auto p = lolita::domain::Point();
            using ComponentGeometry = FiniteElementGeometry<_component, t_domain, t_finite_element>;
            auto const & elt_reference_nodes = lolita::core2::geometry::ElementTraits<t_element, t_domain>::reference_nodes_;
            for (lolita::index i = 0; i < 3; ++i) {
                auto cpt_coordinates = lolita::matrix::Vector<lolita::real, _component.num_nodes_>();
                for (lolita::index j = 0; j < _component.num_nodes_; ++j) {
                    auto const node_tag = getComponentNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<_quadrature, _ord>(index);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @tparam _i
         * @tparam _j
         * @param component_index
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
        lolita::real
        getComponentCurrentQuadratureWeight(
                lolita::index component_index,
                lolita::index index
        )
        const
        requires(!t_element.isPoint())
        {
            auto const & cmp =  this->template getComponents<_i, _j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<_quadrature, _ord>(index);
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @tparam _i
         * @tparam _j
         * @param component_index
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
        lolita::domain::Point
        getComponentCurrentQuadraturePoint(
                lolita::index component_index,
                lolita::index index
        )
        const
        requires(!t_element.isPoint())
        {
            auto p = lolita::domain::Point();
            auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<_quadrature, _ord, _i, _j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

        /**
         * @brief
         * @param point
         * @return
         */
        lolita::domain::Point
        getNormalVector(
                lolita::domain::Point const & point
        )
        const
        requires(t_ElementTraits::isFace())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, t_element.dim_>();
            ru.setZero();
            for (lolita::index i = 0; i < 3; ++i) {
                for (lolita::index j = 0; j < t_element.dim_; ++j) {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.isPoint()) {
                return lolita::domain::Point{0, 0, 0};
            }
            else if constexpr (t_element.isCurve()) {
                return lolita::domain::Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
            }
            else {
                return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            }
        }

        /**
         * @brief
         * @param point
         * @param direction
         * @return
         */
        lolita::domain::Point
        getTangentVector(
                lolita::domain::Point const & point,
                lolita::index direction
        )
        const
        requires(t_element.isCurve())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = lolita::domain::Point();
            for (lolita::index i = 0; i < 3; ++i) {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
            }
            return tangent_vector;
        }

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementBehaviour : virtual FiniteElementConnectivity<t_element, t_domain, t_finite_element>
    {

        /**
         * @brief
         * @tparam t_element_group
         * @param mesh
         */
        template<auto t_element_group>
        void
        setBehaviour(
                lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
        )
        {
            for (auto const & behaviour : mesh.behaviours_) {
                auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & domain) {
                    return * domain == behaviour.domain_tag_;
                };
                auto const has_domain = std::find_if(this->domains_.begin(), this->domains_.end(), is_equal) != this->domains_.end();
                auto const has_unknown = behaviour.unknown_tag_ == t_finite_element.unknown_.tensor_;
                if (has_domain && has_unknown) {
                    std::cout << "setting behaviour : " << behaviour.domain_tag_ << std::endl;
                    behaviour_ = behaviour.behaviour_data_;
                }
            }
        }

        /**
         * @brief
         */
        std::shared_ptr<lolita::behaviour::MgisBehaviourData> behaviour_;

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementFieldLoad : virtual FiniteElementConnectivity<t_element, t_domain, t_finite_element>
    {

    private:

        /**
         * @brief
         */
        using t_Field = lolita::core2::field::TensorPolicy<t_finite_element.unknown_.tensor_, t_domain.dim_>;

        /**
         * @brief loads initializer
         * @return loads member initialized with the zero natural load
         */
        static
        std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_>
        def_loads()
        {
            auto loads = std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_>();
            for (int i = 0; i < t_Field::shape_.rows_; ++i) {
                for (int j = 0; j < t_Field::shape_.cols_; ++j) {
                    loads[j][i] = lolita::core2::finite_element::null_load_ptr;
                }
            }
            return loads;
        }

    public:

        /**
         * @brief
         * @param row
         * @param col
         * @return
         */
        std::shared_ptr<lolita::finite_element::ScalarLoad> const &
        getLoadComponent(
                lolita::integer row,
                lolita::integer col
        )
        const
        {
            return field_load_[col][row];
        }

        /**
         * @brief
         * @param row
         * @param col
         * @return
         */
        std::shared_ptr<lolita::finite_element::ScalarLoad> &
        getLoadComponent(
                lolita::integer row,
                lolita::integer col
        )
        {
            return field_load_[col][row];
        }

        /**
         * @brief
         * @tparam t_element_group
         * @param mesh
         */
        template<auto t_element_group>
        void
        setLoads(
                lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
        )
        {
            for (auto const & load : mesh.loads_) {
                auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & domain) {
                    return * domain == load.domain_tag_;
                };
                auto const has_domain = std::find_if(this->domains_.begin(), this->domains_.end(), is_equal) != this->domains_.end();
                auto const has_unknown = load.unknown_tag_ == t_finite_element.unknown_.tensor_;
                auto const has_dimension = load.element_dim_ == t_element.dim_;
                if (has_domain && has_unknown && has_dimension) {
                    std::cout << "setting load : " << load.domain_tag_ << std::endl;
                    auto const i = load.components_.row_;
                    auto const j = load.components_.col_;
                    getLoadComponent(i, j) = load.scalar_load_;
                }
            }
        }

        /**
         * @brief
         */
        std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_> field_load_ = def_loads();

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementBase;

    namespace basis
    {

        /**
         * @brief
         */
        struct Basis : public lolita::utility::Enumeration<Basis>
        {

            /**
             * @brief
             * @param tag
             */
            constexpr
            Basis(
                    std::basic_string_view<lolita::character> && tag
            )
            :
            lolita::utility::Enumeration<Basis>(std::forward<std::basic_string_view<lolita::character>>(tag))
            {}

            /**
             * @brief
             * @return
             */
            static constexpr
            Basis
            Monomial()
            {
                return Basis("Monomial");
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::boolean
            isMonomial()
            const
            {
                return * this == Monomial();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            Basis
            Lagrange()
            {
                return Basis("Lagrange");
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::boolean
            isLagrange()
            const
            {
                return * this == Lagrange();
            }

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        struct FiniteElementBasis;

        /**
         * @brief Implementation object for the definition of a polynomial basis function
         * @tparam t_element the element on which the monomial is defined
         * @tparam t_basis the basis enum
         * @tparam t_ord the polynomial order of the basis
         */
        template<lolita::core2::geometry::Element t_element, lolita::core2::finite_element::basis::Basis t_basis, lolita::index t_ord>
        struct FiniteElementBasisTraits;

        /**
         * @brief Implementation object for the definition of a monomial basis function
         * @tparam t_element the element on which the monomial is defined
         * @tparam t_basis the monomial basis enum
         * @tparam t_ord the polynomial order of the basis
         */
        template<lolita::core2::geometry::Element t_element, lolita::core2::finite_element::basis::Basis t_basis, lolita::index t_ord>
        requires(t_basis.isMonomial())
        struct FiniteElementBasisTraits<t_element, t_basis, t_ord>
        {

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::integer
            dim()
            {
                return numerics::binomial(t_element.dim_ + t_ord, t_element.dim_);
            }

            /**
             * @brief The basis dimension or cardinality
             */
            lolita::index const static constexpr dim_ = dim();

            /**
             * @brief Implementation of the evaluation and derivative function for the given polynomial basis
             * @tparam t_domain
             * @tparam t_finite_element
             */
            template<lolita::domain::Domain t_domain, auto t_finite_element>
            struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
            {

            private:

                /**
                 * @brief
                 * @return
                 */
                static constexpr
                std::array<std::array<lolita::index, 3>, dim_>
                make_exponents()
                {
                    auto exponents = std::array<std::array<lolita::index, 3>, dim_>();
                    auto row = lolita::index(0);
                    if constexpr (t_element.dim_ == 0) {
                        exponents[row][0] = 0;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                    }
                    else if constexpr (t_element.dim_ == 1) {
                        for (lolita::index i = 0; i < t_ord + 1; ++i) {
                            exponents[row][0] = i;
                            exponents[row][1] = 0;
                            exponents[row][2] = 0;
                            row += 1;
                        }
                    }
                    else if constexpr (t_element.dim_ == 2) {
                        for (lolita::index i = 0; i < t_ord + 1; ++i) {
                            for (lolita::index j = 0; j < i + 1; ++j) {
                                exponents[row][0] = i - j;
                                exponents[row][1] = j;
                                exponents[row][2] = 0;
                                row += 1;
                            }
                        }
                    }
                    else if constexpr (t_element.dim_ == 3) {
                        for (lolita::index i = 0; i < t_ord + 1; ++i) {
                            for (lolita::index j = 0; j < i + 1; ++j) {
                                for (lolita::index k = 0; k < i + 1; ++k) {
                                    if (j + k < i + 1) {
                                        exponents[row][0] = i - (j + k);
                                        exponents[row][1] = k;
                                        exponents[row][2] = j;
                                        row += 1;
                                    }
                                }
                            }
                        }
                    }
                    return exponents;
                }

                /**
                 * @brief
                 */
                std::array<std::array<lolita::index, 3>, dim_> const static constexpr exponents_ = make_exponents();

            public:

                /**
                 * @brief evaluate the polynomial basis function at some point
                 * @param point the point where to evaluate the polynomial basis function
                 * @return a scalar-valued vector of size the dimension of the polynomial basis
                 */
                lolita::matrix::Vector<lolita::real, dim_>
                getBasisEvaluation(
                        lolita::domain::Point const & point
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::index i = 0; i < dim_; ++i) {
                        auto value = lolita::real(1);
                        for (lolita::index j = 0; j < t_element.dim_; ++j) {
                            auto dist = this->getRiemannianDistance(centroid, point, j);
                            //value *= numerics::pow(2.0 * dist / diameters(j), exponents.get(i, j));
                            value *= std::pow(2.0 * dist / diameters(j), exponents_[i][j]);
                        }
                        basis_vector_values(i) = value;
                    }
                    return basis_vector_values;
                }

                /**
                 * @brief evaluate the polynomial basis derivative at some point
                 * @param point the point where to evaluate the polynomial basis function
                 * @param derivative_direction the derivative direction, that must be less than or equal to the element dimension
                 * @return a scalar-valued vector of size the dimension of the polynomial basis
                 */
                lolita::matrix::Vector<lolita::real, dim_>
                getBasisDerivative(
                        lolita::domain::Point const & point,
                        lolita::index derivative_direction
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::index i = 0; i < dim_; ++i) {
                        auto value = lolita::real(1);
                        for (lolita::index j = 0; j < t_element.dim_; ++j) {
                            if (j != derivative_direction) {
                                auto dist = this->getRiemannianDistance(centroid, point, j);
                                value *= std::pow(2.0 * (dist) / diameters(j), exponents_[i][j]);
                            }
                            else {
                                if (exponents_[i][j] > 0) {
                                    auto c = 2.0 * exponents_[i][j] / diameters(j);
                                    auto dist = this->getRiemannianDistance(centroid, point, j);
                                    value *= c * std::pow(2.0 * (dist) / diameters(j), exponents_[i][j] - 1);
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

    }

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementBasis
    {

        /**
         * @brief
         * @tparam t_basis
         * @tparam t_ord
         * @param point
         * @return
         */
        template<lolita::core2::finite_element::basis::Basis t_basis, lolita::index t_ord>
        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, t_basis, t_ord>::dim_>
        getBasisEvaluation(
                lolita::domain::Point const & point
        )
        const
        {
            using t_Implementation = typename basis::FiniteElementBasisTraits<t_element, t_basis, t_ord>::template Implementation<t_domain, t_finite_element>;
            return static_cast<t_Implementation const *>(this)->getBasisEvaluation(point);
        }

        /**
         * @brief
         * @tparam t_basis
         * @tparam t_ord
         * @param point
         * @param derivative_direction
         * @return
         */
        template<lolita::core2::finite_element::basis::Basis t_basis, lolita::index t_ord>
        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, t_basis, t_ord>::dim_>
        getBasisDerivative(
                lolita::domain::Point const & point,
                lolita::index derivative_direction
        )
        const
        {
            using t_Implementation = typename basis::FiniteElementBasisTraits<t_element, t_basis, t_ord>::template Implementation<t_domain, t_finite_element>;
            return static_cast<t_Implementation const *>(this)->getBasisDerivative(point, derivative_direction);
        }

    };

    namespace unknown
    {

        /**
         * @brief
         */
        struct Unknown : public lolita::utility::Enumeration<Unknown>
        {

            /**
             * @brief
             * @param tag
             */
            constexpr
            Unknown(
                    std::basic_string_view<lolita::character> && tag
            )
            :
            lolita::utility::Enumeration<Unknown>(std::forward<std::basic_string_view<lolita::character>>(tag))
            {}

            /**
             * @brief
             * @return
             */
            static constexpr
            Unknown
            Structural()
            {
                return Unknown("Structural");
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::boolean
            isStructural()
            const
            {
                return * this == Structural();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            Unknown
            Subsidiary()
            {
                return Unknown("Subsidiary");
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::boolean
            isSubsidiary()
            const
            {
                return * this == Subsidiary();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            std::array<Unknown, 2>
            UnknownTypes()
            {
                return std::array<Unknown, 2>{
                        Unknown::Structural(),
                        Unknown::Subsidiary(),
                };
            }

        };

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_dim
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_dim>
        struct ScalarUnknown;

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_dim
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_dim>
        requires(t_unknown.isStructural())
        struct ScalarUnknown<t_unknown, t_dim>
        {

            /**
             * @brief The binding coefficient vector type
             */
            using CoordinateVector = lolita::matrix::Vector<lolita::integer, t_dim>;

            /**
             * @brief The binding coefficient vector type
             */
            using CoefficientVector = lolita::matrix::Vector<lolita::real, t_dim>;

            /**
             * @brief
             * @return
             */
            lolita::boolean
            isBound()
            const
            {
                return !(bindings_coefficients_ == nullptr);
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector &
            getUnknownsCoefficients()
            {
                return * unknowns_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector const &
            getUnknownsCoefficients()
            const
            {
                return * unknowns_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector &
            getBindingsCoefficients()
            {
                return * bindings_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector const &
            getBindingsCoefficients()
            const
            {
                return * bindings_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoordinateVector &
            getUnknownsCoordinates()
            {
                return * unknowns_coordinates_;
            }

            /**
             * @brief
             * @return
             */
            CoordinateVector const &
            getUnknownsCoordinates()
            const
            {
                return * unknowns_coordinates_;
            }

            /**
             * @brief
             * @return
             */
            CoordinateVector &
            getBindingsCoordinates()
            {
                return * bindings_coordinates_;
            }

            /**
             * @brief
             * @return
             */
            CoordinateVector const &
            getBindingsCoordinates()
            const
            {
                return * bindings_coordinates_;
            }

            /**
             * @brief
             */
            std::unique_ptr<CoefficientVector> unknowns_coefficients_;

            /**
             * @brief
             */
            std::unique_ptr<CoefficientVector> bindings_coefficients_;

            /**
             * @brief
             */
            std::unique_ptr<CoordinateVector> unknowns_coordinates_;

            /**
             * @brief
             */
            std::unique_ptr<CoordinateVector> bindings_coordinates_;

        };

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_dim
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_dim>
        requires(t_unknown.isSubsidiary())
        struct ScalarUnknown<t_unknown, t_dim>
        {

            /**
             * @brief The binding coefficient vector type
             */
            using CoefficientVector = lolita::matrix::Vector<lolita::real, t_dim>;

            /**
             * @brief
             * @return
             */
            lolita::boolean
            isBound()
            const
            {
                return !(bindings_coefficients_ == nullptr);
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector &
            getUnknownsCoefficients()
            {
                return * unknowns_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector const &
            getUnknownsCoefficients()
            const
            {
                return * unknowns_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector &
            getBindingsCoefficients()
            {
                return * bindings_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector const &
            getBindingsCoefficients()
            const
            {
                return * bindings_coefficients_;
            }

            /**
             * @brief
             */
            std::unique_ptr<CoefficientVector> unknowns_coefficients_;

            /**
             * @brief
             */
            std::unique_ptr<CoefficientVector> bindings_coefficients_;

        };

        /**
         * @brief
         * @tparam t_domain
         * @tparam t_field
         * @tparam t_dim
         */
        template<
                lolita::domain::Domain t_domain,
                lolita::field::Field t_field,
                lolita::integer t_dim,
                lolita::core2::finite_element::unknown::Unknown t_unknown
        >
        struct FieldUnknown
        {

            /**
             * @brief The unknown type object
             */
            lolita::core2::finite_element::unknown::Unknown const static constexpr unknown_ = t_unknown;

            /**
             * @brief The field object
             */
            lolita::field::Field const static constexpr field_ = t_field;

            /**
             * @brief The dimension or cardinality of the unknown
             */
            lolita::integer const static constexpr dim_ = t_dim;

        private:

            /**
             * @brief The field type
             */
            using t_Field = lolita::core2::field::TensorPolicy<t_field, t_domain.dim_>;

        public:

            /**
             * @brief
             * @param row
             * @param col
             * @return
             */
            ScalarUnknown<t_unknown, t_dim> &
            getUnknownComponent(
                    lolita::integer row,
                    lolita::integer col
            )
            {
                return field_unknown_[col][row];
            }

            /**
             * @brief
             * @param row
             * @param col
             * @return
             */
            ScalarUnknown<t_unknown, t_dim> const &
            getUnknownComponent(
                    lolita::integer row,
                    lolita::integer col
            )
            const
            {
                return field_unknown_[col][row];
            }

            /**
             * @brief
             * @param i
             * @param j
             */
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j
            )
            requires(t_dim > 0 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero());
            }

            /**
             * @brief
             * @param i
             * @param j
             * @param dim_unknown
             */
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim_unknown
            )
            requires(t_dim == -1 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim_unknown));
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim > 0 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero());
                auto & k = mesh.systems_[_finite_element_index].num_unknowns_;
                getUnknownComponent(i, j).unknowns_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(k, k + t_dim));
                k += t_dim;
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param dim_unknown
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim_unknown,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim == -1 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim_unknown));
                auto & k = mesh.systems_[_finite_element_index].num_unknowns_;
                getUnknownComponent(i, j).unknowns_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(dim_unknown, k, k + dim_unknown));
                k += dim_unknown;
            }

            /**
             * @brief
             * @param i
             * @param j
             */
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j
            )
            requires(t_dim > 0 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero());
            }

            /**
             * @brief
             * @param i
             * @param j
             * @param dim_unknown
             */
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim_unknown
            )
            requires(t_dim == -1 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim_unknown));
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim > 0 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero());
                auto & k = mesh.systems_[_finite_element_index].num_bindings_;
                getUnknownComponent(i, j).bindings_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(k, k + t_dim));
                k += t_dim;
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param dim_unknown
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim_unknown,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim == -1 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim_unknown));
                auto & k = mesh.systems_[_finite_element_index].num_bindings_;
                getUnknownComponent(i, j).bindings_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(dim_unknown, k, k + dim_unknown));
                k += dim_unknown;
            }

            /**
             * @brief
             */
            std::array<std::array<ScalarUnknown<t_unknown, t_dim>, t_Field::shape_.rows_>, t_Field::shape_.cols_> field_unknown_;

        };

        namespace detail
        {

            template<typename t_T>
            struct is_unknown : public std::false_type {};

            template<
                    lolita::domain::Domain t_domain,
                    lolita::field::Field t_field,
                    lolita::integer t_dim,
                    lolita::core2::finite_element::unknown::Unknown t_unknown
            >
            struct is_unknown<lolita::core2::finite_element::unknown::FieldUnknown<t_domain, t_field, t_dim, t_unknown>> : public std::true_type {};

        }

        /**
         * @brief
         * @tparam t_T
         */
        template<typename t_T>
        concept FieldUnknownConcept = lolita::core2::finite_element::unknown::detail::is_unknown<t_T>::value;

        /**
         * @brief
         * @tparam t_FieldUnknown
         */
        template<lolita::core2::finite_element::unknown::FieldUnknownConcept... t_FieldUnknown>
        struct FieldUnknownCollection
        {

            /**
             * @brief
             */
            using Unknowns = std::tuple<t_FieldUnknown...>;

            /**
             * @brief
             */
            std::array<lolita::integer, sizeof...(t_FieldUnknown)> const static constexpr dim_unknowns_ = {t_FieldUnknown::dim_...};

            /**
             * @brief
             */
            lolita::integer const static constexpr num_unknowns_ = sizeof...(t_FieldUnknown);

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::boolean
            hasStructuralUnknowns()
            {
                return ((t_FieldUnknown::unknown_type_.isStructural()) || ...);
            }

            /**
             * @brief
             * @tparam t_domain
             * @tparam _unknown_type
             * @return
             */
            template<lolita::domain::Domain t_domain, lolita::core2::finite_element::unknown::Unknown t_unknown>
            static constexpr
            lolita::integer
            getNumUnknowns()
            {
                auto num_unknowns = lolita::integer(0);
                auto set_num_unknowns  = [&] <lolita::integer _i = 0> (auto & self) constexpr mutable {
                    if constexpr (std::tuple_size_v<Unknowns> > 0) {
                        using _Unknown = std::tuple_element_t<_i, Unknowns>;
                        if constexpr (_Unknown::unknown_ == t_unknown) {
                            num_unknowns += _Unknown::dim_ * lolita::core2::field::TensorPolicy<_Unknown::field_, t_domain.dim_>::shape_.size_;
                        }
                        if constexpr (_i < std::tuple_size_v<Unknowns> - 1) {
                            self.template operator()<_i + 1>(self);
                        }
                    }
                };
                set_num_unknowns(set_num_unknowns);
                return num_unknowns;
            }

            /**
             * @brief
             * @tparam t_i
             * @return
             */
            template<lolita::integer t_i>
            std::tuple_element_t<t_i, Unknowns> const &
            getUnknown()
            const
            {
                return std::get<t_i>(unknowns_);
            }

            /**
             * @brief
             * @tparam t_i
             * @return
             */
            template<lolita::integer t_i>
            std::tuple_element_t<t_i, Unknowns> &
            getUnknown()
            {
                return std::get<t_i>(unknowns_);
            }

            /**
             * @brief
             */
            Unknowns unknowns_;

        };

        /**
         * @brief Default initialization is zero unknowns
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        struct FiniteElementFieldUnknownsTraits
        {

            /**
             * @brief
             */
            using FieldUnknowns = lolita::core2::finite_element::unknown::FieldUnknownCollection<>;

            /**
             * @brief
             */
            struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
            {

                /**
                 * @brief
                 * @tparam t_element_group
                 * @param mesh
                 */
                template<auto t_element_group>
                void
                setUnknowns(
                        lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
                )
                {
                    std::cout << "----> setting unknown empty" << std::endl;
                }

            };

//            /**
//             * @brief
//             * @tparam t_element_group
//             * @param mesh
//             */
//            template<auto t_element_group>
//            void
//            setUnknowns(
//                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
//            )
//            {
//                std::cout << "----> setting unknown empty" << std::endl;
//            }
//
//            /**
//             * @brief
//             */
//            FieldUnknowns field_unknowns_;

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        requires(t_finite_element.discretization_.isHHO() && t_element.isSub(t_domain, 0))
        struct FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>
        {

        private:

            /**
             * @brief
             */
            using t_Basis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<
                    t_element,
                    lolita::core2::finite_element::basis::Basis::Monomial(),
                    t_finite_element.discretization_.ord_cell_
            >;

            /**
             * @brief
             */
            using t_Field = typename lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>::Field;

        public:

            /**
             * @brief
             */
            using FieldUnknowns = lolita::core2::finite_element::unknown::FieldUnknownCollection<
                    lolita::core2::finite_element::unknown::FieldUnknown<
                            t_domain,
                            t_finite_element.unknown_.tensor_,
                            t_Basis::dim_,
                            lolita::core2::finite_element::unknown::Unknown::Subsidiary()
                    >
            >;

            /**
             * @brief
             */
            struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
            {

                /**
                 * @brief
                 * @tparam t_element_group
                 * @param mesh
                 */
                template<auto t_element_group>
                void
                setUnknowns(
                        lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
                )
                {
                    for (int i = 0; i < t_Field::shape_.rows_; ++i) {
                        for (int j = 0; j < t_Field::shape_.cols_; ++j) {
                            this->field_unknowns_.template getUnknown<0>().setUnknownCoefficients(i, j);
                            if (this->getLoadComponent(i, j)->load_.isConstraint()) {
                                this->field_unknowns_.template getUnknown<0>().setBindingCoefficients(i, j);
                            }
                        }
                    }
                }

            };

//            /**
//             * @brief
//             * @tparam t_element_group
//             * @param mesh
//             */
//            template<auto t_element_group>
//            void
//            setUnknowns(
//                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
//            )
//            {
//                for (int i = 0; i < t_Field::shape_.rows_; ++i) {
//                    for (int j = 0; j < t_Field::shape_.cols_; ++j) {
//                        field_unknowns_.template getUnknown<0>().setUnknownCoefficients(i, j);
//                        if (this->getLoad(i, j)->load_ == lolita::finite_element::Load::Constraint()) {
//                            field_unknowns_.template getUnknown<0>().setBindingCoefficients(i, j);
//                        }
//                    }
//                }
//            }
//
//            /**
//             * @brief
//             */
//            FieldUnknowns field_unknowns_;

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        requires(t_finite_element.discretization_.isHHO() && t_element.isSub(t_domain, 1))
        struct FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>
        {

        private:

            /**
             * @brief
             */
            using _Basis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<
                    t_element,
                    lolita::core2::finite_element::basis::Basis::Monomial(),
                    t_finite_element.discretization_.ord_face_
            >;

            /**
             * @brief
             */
            using t_Field = typename lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>::Field;

        public:

            /**
             * @brief
             */
            using FieldUnknowns = lolita::core2::finite_element::unknown::FieldUnknownCollection<
                    lolita::core2::finite_element::unknown::FieldUnknown<
                            t_domain,
                            t_finite_element.unknown_.tensor_,
                            _Basis::dim_,
                            lolita::core2::finite_element::unknown::Unknown::Structural()
                    >
            >;

            /**
             * @brief
             */
            struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
            {

                /**
                 * @brief
                 * @tparam t_element_group
                 * @param mesh
                 */
                template<auto t_element_group>
                void
                setUnknowns(
                        lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
                )
                {
                    for (int i = 0; i < t_Field::shape_.rows_; ++i) {
                        for (int j = 0; j < t_Field::shape_.cols_; ++j) {
                            this->field_unknowns_.template getUnknown<0>().template setUnknownCoefficients<t_finite_element>(i, j, mesh);
                            if (this->getLoadComponent(i, j)->load_.isConstraint()) {
                                std::cout << "setting const " << std::endl;
                                this->field_unknowns_.template getUnknown<0>().template setBindingCoefficients<t_finite_element>(i, j, mesh);
                            }
                        }
                    }
                }

            };

//            /**
//             * @brief
//             * @tparam t_element_group
//             * @param mesh
//             */
//            template<auto t_element_group>
//            void
//            setUnknowns(
//                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
//            )
//            {
//                for (int i = 0; i < t_Field::shape_.rows_; ++i) {
//                    for (int j = 0; j < t_Field::shape_.cols_; ++j) {
//                        field_unknowns_.template getUnknown<0>().template setUnknownCoefficients<t_finite_element>(i, j, mesh);
//                        if (this->getLoad(i, j)->load_.isConstraint()) {
//                            field_unknowns_.template getUnknown<0>().template setBindingCoefficients<t_finite_element>(i, j, mesh);
//                        }
//                    }
//                }
//            }
//
//            /**
//             * @brief
//             */
//            FieldUnknowns field_unknowns_;

        };

    }

    /**
     * @brief Default initialization is zero unknowns
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementFieldUnknowns
    {

        /**
         * @brief
         */
        using t_FieldUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>::FieldUnknowns;

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown>
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            using t_FiniteElementFieldUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>::FieldUnknowns;
            return t_FiniteElementFieldUnknowns::template getNumUnknowns<t_domain, t_unknown>();
        }

        /**
         * @brief
         * @tparam t_element_group
         * @param mesh
         */
        template<auto t_element_group>
        void
        setUnknowns(
                lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
        )
        {
            using t_Implementation = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>::Implementation;
            static_cast<t_Implementation *>(this)->setUnknowns(mesh);
        }

        /**
         * @brief
         */
        t_FieldUnknowns field_unknowns_;

    };

    namespace cell
    {



    }

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementTraits
    {

    private:

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getOrdQuadrature()
        {
            auto ord_quadrature = t_finite_element.ord_quadrature_;
            return t_domain.frame_ == lolita::domain::Frame::AxiSymmetric() ? 2 * ord_quadrature + 1 : 2 * ord_quadrature;
        }


    public:

        /**
         * @brief
         */
        lolita::finite_element::FiniteElementConcept auto const static constexpr finite_element_ = t_finite_element;

        /**
         * @brief
         */
        using ElementDescription = lolita::core2::geometry::ElementTraits<t_element, t_domain>;

        /**
         * @brief
         */
        using Quadrature = lolita::core2::geometry::ElementQuadratureTraits<t_element, t_finite_element.quadrature_, getOrdQuadrature()>;

        /**
         * @brief
         */
        using Field = lolita::core2::field::TensorPolicy<t_finite_element.unknown_.tensor_, t_domain.dim_>;

        /**
         * @brief
         * @tparam t_method
         * @return
         */
        template<lolita::finite_element::FiniteElementMethod t_method>
        static constexpr
        lolita::boolean
        hasMethod()
        {
            return finite_element_.discretization_ == t_method;
        }

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown>
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            auto num_element_unknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getNumUnknowns<t_unknown>();
            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                auto const constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
                auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getNumUnknowns<t_unknown>();
                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
                num_element_unknowns += num_component_unknowns * num_components;
                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            get_num_components_unknowns(get_num_components_unknowns);
            return num_element_unknowns;
        }

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_i
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_i>
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            if constexpr (t_i == 0) {
                return FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getNumUnknowns<t_unknown>();
            }
            else {
                auto num_unknowns = 0;
                auto get_num_components_unknowns = [&] <lolita::integer t_j = 0> (auto & self) constexpr mutable {
                    auto const constexpr t_component = ElementDescription::template getComponent<t_i - 1, t_j>();
                    auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getNumUnknowns<t_unknown>();
                    auto num_components = ElementDescription::template getNumComponents<t_i - 1, t_j>();
                    num_unknowns += num_component_unknowns * num_components;
                    if constexpr (t_j < ElementDescription::template getNumComponents<t_i - 1>() - 1) {
                        self.template operator()<t_j + 1>(self);
                    }
                };
                get_num_components_unknowns(get_num_components_unknowns);
                return num_unknowns;
            }
        }

        /**
         * @brief
         * @tparam t_mapping
         * @return
         */
        template<lolita::field::Mapping t_mapping>
        static constexpr
        lolita::integer
        getMappingSize()
        {
            return lolita::core2::field::MappingPolicy<t_finite_element.unknown_.tensor_, t_domain, t_mapping>::shape_.size_;
        }

        /**
         * @brief
         * @tparam t_mapping
         * @return
         */
        template<lolita::field::Mapping t_mapping>
        static constexpr
        lolita::matrix::VectorBlock
        getMappingBlock()
        {

            auto mapping_row = lolita::integer(0);
            for (auto mapping : t_finite_element.unknown_.mappings_) {
                if (mapping == t_mapping) {
                    return lolita::matrix::VectorBlock(mapping_row, mapping_row + getMappingSize<t_mapping>());
                }
                mapping_row += getMappingSize<t_mapping>();
            }
            return lolita::matrix::VectorBlock();
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getOperatorNumRows()
        {
            auto mapping_size = lolita::integer(0);
            auto set_dim_mapping = [&] <lolita::integer _i = 0> (auto & self) constexpr mutable {
                mapping_size += getMappingSize<t_finite_element.unknown_.mappings_[_i]>();
                if constexpr (_i < t_finite_element.unknown_.mappings_.size() - 1) {
                    self.template operator()<_i + 1>(self);
                }
            };
            set_dim_mapping(set_dim_mapping);
            return mapping_size;
        }

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_element_group
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, auto t_element_group>
    struct FiniteElement : FiniteElementConnectivity<t_element, t_domain, t_element_group>
    {

        /**
         * @brief
         */
        lolita::boolean const static constexpr is_finite_element_ = false;

    private:

        /**
         * @brief
         */
        struct ElementGroupTraits
        {

            /**
             * @brief
             */
            using ElementPointers = typename std::remove_cvref_t<decltype(t_element_group)>::template ElementPointers<FiniteElement, t_element, t_domain>;

            /**
             * @brief
             */
            using Elements = typename std::remove_cvref_t<decltype(t_element_group)>::template Elements<FiniteElement, t_element, t_domain>;

        };

        /**
         * @brief
         */
        using t_ElementPointers = typename ElementGroupTraits::ElementPointers;

        /**
         * @brief
         */
        using t_Elements = typename ElementGroupTraits::Elements;

    public:

        /**
         * @brief
         * @tparam t_i
         * @return
         */
        template<lolita::index t_i>
        std::tuple_element_t<t_i, t_ElementPointers> const &
        getElement()
        const
        {
            return std::get<t_i>(elements_);
        }

        /**
         * @brief
         * @tparam t_i
         * @return
         */
        template<lolita::index t_i>
        std::tuple_element_t<t_i, t_ElementPointers> &
        getElement()
        {
            return std::get<t_i>(elements_);
        }

        /**
         * @brief
         */
        void
        make()
        {
            auto make_elements = [&] <lolita::index t_k = 0> (auto & make_elements_imp) mutable {
                using t_Element = std::tuple_element_t<t_k, t_Elements>;
                this->template getElement<t_k>() = std::make_shared<t_Element>(t_Element());
                if constexpr (!t_Element::is_finite_element_) {
                    this->template getElement<t_k>()->make();
                }
                if constexpr (t_k < std::tuple_size_v<t_Elements> - 1) {
                    make_elements_imp.template operator()<t_k + 1u>(make_elements_imp);
                }
            };
            make_elements(make_elements);
        }

        template<lolita::core2::geometry::Element t__element = t_element, auto t_arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<t_element> const & initialization_data,
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
        );

        template<lolita::core2::geometry::Element ___element = t_element, auto t_arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<t_element> const & initialization_data,
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
        )
        requires(!___element.isPoint())
        {
            /*
             *
             */
            auto get_element_hash = [&] <lolita::core2::geometry::Element __element> (
                    std::array<lolita::index, __element.num_nodes_> node_tags
            )
            {
                auto element_hash = std::basic_stringstream<lolita::character>();
                if constexpr (node_tags.size() > 0) {
                    std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
                }
                for (int i = 0; i < node_tags.size(); ++i) {
                    element_hash << node_tags[i];
                }
                return element_hash.str();
            };
            /*
             * make element
             */
            auto make_element = [&] <lolita::core2::geometry::Element __element = t_element, lolita::index _i = 0, lolita::index _j = 0> (
                    std::shared_ptr<lolita::core2::finite_element::FiniteElement<__element, t_domain, t_element_group>> & ptr_element,
                    lolita::core2::mesh::ElementInitializationData<__element> const & element_initialization_data,
                    auto & make_element_imp
            )
                    mutable
            {
                if constexpr (!__element.isPoint()) {
                    using __ElementDescription = lolita::core2::geometry::ElementTraits<__element, t_domain>;
                    using __ComponentDescription = lolita::core2::geometry::ElementTraits<__ElementDescription::template getComponent<_i, _j>(), t_domain>;
                    using __MeshDescription = lolita::core2::geometry::DomainTraits<t_domain>;
                    auto const constexpr _is_initialized = _i == 0 && _j == 0;
                    auto const constexpr _component = __ElementDescription::template getComponent<_i, _j>();
                    auto const constexpr _component_coordinates = __MeshDescription::template getElementCoordinates<_component>();
                    auto const constexpr _neighbour_coordinates = __ComponentDescription::template getNeighbourCoordinates<__element>();
                    auto const constexpr _element_coordinates = __MeshDescription::template getElementCoordinates<__element>();
                    auto const constexpr _node_coordinates = __ElementDescription::template getComponentCoordinates<lolita::core2::geometry::Element::Node()>();
                    using __Component = lolita::core2::finite_element::FiniteElement<_component, t_domain, t_element_group>;
                    using __Self = lolita::core2::finite_element::FiniteElement<__element, t_domain, t_element_group>;
                    auto & components = mesh_data.elements_.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
                    auto & element_component_array = ptr_element->template getComponents<_i, _j>();
                    for (auto i = 0; i < element_component_array.size(); ++i) {
                        auto component_hash = std::basic_string<lolita::character>();
                        if constexpr(!_component.isPoint()) {
                            auto component_initialization_data = lolita::core2::mesh::ElementInitializationData<_component>();
                            for (int j = 0; j < _component.num_nodes_; ++j) {
                                auto const k = __Self::template getComponentNodeConnection<_i, _j>(i, j);
                                component_initialization_data.node_tags_[j] = element_initialization_data.node_tags_[k];
                            }
                            component_initialization_data.tag_ = components.size();
                            component_hash = get_element_hash.template operator ()<_component>(component_initialization_data.node_tags_);
                            if (!components.contains(component_hash)) {
                                auto ptr_component = std::make_shared<__Component>(__Component());
                                make_element_imp.template operator ()<_component, 0, 0>(ptr_component, component_initialization_data, make_element_imp);
                            }
                        }
                        else {
                            component_hash = std::to_string(element_initialization_data.node_tags_[__Self::template getComponentNodeConnection<_i, _j>(i, 0)]);
                        }
                        element_component_array[i] = components[component_hash];
                        components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
                    }
                    if constexpr (_j < __ElementDescription::template getNumComponents<_i>() - 1) {
                        make_element_imp.template operator()<__element, _i, _j + 1u>(ptr_element, element_initialization_data, make_element_imp);
                    }
                    else if constexpr (_i < __ElementDescription::getNumComponents() - 1) {
                        make_element_imp.template operator()<__element, _i + 1u, 0u>(ptr_element, element_initialization_data, make_element_imp);
                    }
                    if constexpr (_is_initialized) {
                        auto & elements = mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>();
                        auto const & nodes = ptr_element->template getComponents<_node_coordinates.dim_, _node_coordinates.tag_>();
                        auto domains = std::unordered_set<std::shared_ptr<std::basic_string<lolita::character>>>();
                        for (auto const & domain : nodes[0]->domains_) {
                            auto has_domain = true;
                            for (int j = 1; j < __element.num_nodes_; ++j) {
                                has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                            }
                            if (has_domain) {
                                domains.insert(domain);
                            }
                        }
                        ptr_element->tag_ = element_initialization_data.tag_;
                        ptr_element->domains_.assign(domains.begin(), domains.end());
                        ptr_element->make();
                        auto element_hash = get_element_hash.template operator ()<__element>(element_initialization_data.node_tags_);
                        elements.insert({element_hash, ptr_element});
                    }
                }
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FiniteElement>(FiniteElement());
            make_element(ptr_element, initialization_data, make_element);
        }

        template<lolita::core2::geometry::Element ___element = t_element, auto t_arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<t_element> const & initialization_data,
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
        )
        requires(___element.isPoint())
        {
            /*
             *
             */
            auto make_element = [&] (
                    std::shared_ptr<lolita::core2::finite_element::FiniteElement<t_element, t_domain, t_element_group>> & ptr_element
            )
                    mutable
            {
                auto const constexpr _element_coordinates = lolita::core2::geometry::DomainTraits<t_domain>::template getElementCoordinates<t_element>();
                ptr_element->tag_ = initialization_data.tag_;
                ptr_element->domains_ = initialization_data.domains_;
                ptr_element->coordinates_ = std::make_shared<lolita::domain::Point>(initialization_data.coordinates_);
//                ptr_element->coordinates_ = std::make_unique<lolita::domain::Point>(initialization_data.coordinates_);
                ptr_element->make();
                auto elem_hash = std::to_string(initialization_data.tag_);
                mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>().insert({elem_hash, ptr_element});
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FiniteElement>(FiniteElement());
            make_element(ptr_element);
        }

        template<lolita::core2::geometry::Element __element = t_element, auto t_arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
        );

        template<lolita::core2::geometry::Element __element = t_element, auto t_arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
        )
        requires(!__element.isPoint())
        {
            /*
             *
             */
            auto initialize_element = [&] <lolita::index _k = 0u> (
                    auto & initialize_element_imp
            )
                    mutable
            {
                /*
                 *
                 */
                auto initialize_components = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_components_imp
                )
                        mutable
                {
                    for (int i = 0; i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<__i, __j>(); ++i) {
                        auto & rhs = this->template getComponents<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getComponents<__i, __j>()[i];
                        lhs = rhs;
                    }
                    if constexpr (__j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<__i>() - 1) {
                        initialize_components_imp.template operator()<__i, __j + 1u>(initialize_components_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::getNumComponents() - 1) {
                        initialize_components_imp.template operator()<__i + 1u, 0u>(initialize_components_imp);
                    }
                };
                /*
                 *
                 */
                auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_neighbours_imp
                )
                        mutable
                {
                    for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                        auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (__j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::getNumNeighbours() - 1) {
                        initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                    }
                };
                /*
                 *
                 */
                initialize_components(initialize_components);
                initialize_neighbours(initialize_neighbours);
                this->template getElement<_k>()->tag_ = this->tag_;
                this->template getElement<_k>()->domains_ = this->domains_;
                this->template getElement<_k>()->initialize(mesh_data);
                if constexpr (_k < std::tuple_size_v<t_Elements> - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        template<lolita::core2::geometry::Element __element = t_element, auto _arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<t_domain, _arg> & mesh_data
        )
        requires(__element.isPoint())
        {
            /*
              *
              */
            auto initialize_element = [&] <lolita::index _k = 0u> (
                    auto & initialize_element_imp
            )
                    constexpr mutable
            {
                /*
                 *
                 */
                auto initialize_coordinates = [&] ()
                {
                    this->template getElement<_k>()->components_ = this->components_;
                    this->template getElement<_k>()->coordinates_ = this->coordinates_;
                };
                /*
                 *
                 */
                auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_neighbours_imp
                )
                        mutable
                {
                    for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                        auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (__j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::getNumNeighbours() - 1) {
                        initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                    }
                };
                /*
                 *
                 */
                initialize_coordinates();
                initialize_neighbours(initialize_neighbours);
                this->template getElement<_k>()->tag_ = this->tag_;
                this->template getElement<_k>()->domains_ = this->domains_;
                this->template getElement<_k>()->initialize(mesh_data);
                if constexpr (_k < std::tuple_size_v<t_Elements> - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        t_ElementPointers elements_;

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElement<t_element, t_domain, t_finite_element>
    :
    FiniteElementGeometry<t_element, t_domain, t_finite_element>,
    FiniteElementBehaviour<t_element, t_domain, t_finite_element>,
    FiniteElementFieldLoad<t_element, t_domain, t_finite_element>,
    FiniteElementBasis<t_element, t_domain, t_finite_element>,
    FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>
    {

        lolita::boolean const static constexpr is_finite_element_ = true;

        template<auto t_element_group>
        void
        initialize(
                lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
        )
        {
            std::cout << "making element : " << this->hash() << std::endl;
            this->setBehaviour(mesh);
            this->setLoads(mesh);
            this->setUnknowns(mesh);
            this->getCurrentCoordinates();
            this->action();
        }

        void
        action()
        const
        {
            auto pt = lolita::domain::Point();
            this->template getBasisEvaluation<lolita::core2::finite_element::basis::Basis::Monomial(), 1>(pt);
        }

    };

}

#endif //LOLITA_LOLITA_CORE_5_HXX



//namespace temp
//{
//
//    /**
//     * @brief
//     * @tparam t_element
//     * @tparam t_domain
//     * @tparam t_finite_element
//     */
//    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//    struct FiniteElementGeometry2
//    {
//
//    private:
//
//        /**
//         * @brief
//         */
//        using t_ElementTraits = lolita::core2::geometry::ElementGeometryTraits<t_element, t_domain>;
//
//        /**
//         * @brief
//         */
//        using t_FiniteElementTraits = lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>;
//
//        /**
//         * @brief
//         */
//        using t_Quadrature = typename t_FiniteElementTraits::Quadrature;
//
//    public:
//
//        /**
//         * @brief Get the current coordinates of the point/node in the mesh
//         * @return the vector with x, y and z coordinates
//         */
//        lolita::domain::Point &
//        getCurrentCoordinates()
//        requires(t_element.isPoint())
//        {
//            return connectivity_->coordinates_;
//        }
//
//        /**
//         * @brief Get the current coordinates of the point/node in the mesh
//         * @return the vector with x, y and z coordinates
//         */
//        lolita::domain::Point const &
//        getCurrentCoordinates()
//        const
//        requires(t_element.isPoint())
//        {
//            return connectivity_->coordinates_;
//        }
//
//        /**
//         * @brief Get the current coordinates of the element in the mesh
//         * @return the matrix composed of the column vectors with x, y and z coordinates of each node
//         */
//        lolita::matrix::Matrix<lolita::real, 3, t_element.numNodes()>
//        getCurrentCoordinates()
//        const
//        requires(!t_element.isPoint())
//        {
//            auto current_nodes_coordinates = lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>();
//            auto count = lolita::index(0);
//            for (auto const & node : connectivity_->template getComponents<t_element.dim_ - 1, 0>()) {
//                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
//                count ++;
//            }
//            return current_nodes_coordinates;
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        static
//        lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.numNodes(), matrix::col_major> const>
//        getReferenceCoordinates()
//        {
//            using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.numNodes(), matrix::col_major> const>;
//            return _ReferenceCoordinates(lolita::core2::geometry::ElementGeometryTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
//        }
//
//        /**
//         * @brief
//         * @param nodal_field_values
//         * @param reference_point
//         * @return
//         */
//        static
//        lolita::real
//        getShapeMappingEvaluation(
//                lolita::matrix::Vector<lolita::real, t_element.numNodes()> const & nodal_field_values,
//                lolita::domain::Point const & reference_point
//        )
//        {
//            return lolita::core2::geometry::ElementGeometryTraits<t_element, t_domain>::getShapeMappingEvaluation(nodal_field_values, reference_point);
//        }
//
//        /**
//         * @brief
//         * @param nodal_field_values
//         * @param reference_point
//         * @param derivative_direction
//         * @return
//         */
//        static
//        lolita::real
//        getShapeMappingDerivative(
//                lolita::matrix::Vector<lolita::real, t_element.numNodes()> const & nodal_field_values,
//                lolita::domain::Point const & reference_point,
//                lolita::index derivative_direction
//        )
//        {
//            return lolita::core2::geometry::ElementGeometryTraits<t_element, t_domain>::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
//        }
//
//        /*
//         *
//         */
//
//        /**
//         * @brief
//         * @param point
//         * @return
//         */
//        lolita::real
//        getShapeMappingDifferential(
//                lolita::domain::Point const & point
//        )
//        const
//        {
//            auto const current_coordinates = this->getCurrentCoordinates();
//            auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>();
//            auto du = lolita::real(0);
//            ru.setZero();
//            for (lolita::index i = 0; i < t_domain.dim_; ++i) {
//                for (lolita::index j = 0; j < t_element.dim_; ++j) {
//                    ru(i, j) = FiniteElementGeometry2::getShapeMappingDerivative(current_coordinates.row(i), point, j);
//                }
//            }
//            if constexpr (t_element.dim_ == 3) {
//                du = numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
//            }
//            else if constexpr (t_element.dim_ == 2) {
//                du = numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
//            }
//            else {
//                du = numerics::abs(ru.col(0).norm());
//            }
//            if constexpr (t_domain.frame_ == lolita::domain::Frame::AxiSymmetric()) {
//                lolita::real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
//                if (r0 < 1.e-10) {
//                    r0 = 1.e-10;
//                }
//                du *= 2.0 * numerics::pi * r0;
//            }
//            return du;
//        }
//
//        /**
//         * @brief
//         * @param first_point
//         * @param second_point
//         * @param direction
//         * @return
//         */
//        lolita::real
//        getRiemannianDistance(
//                lolita::domain::Point const & first_point,
//                lolita::domain::Point const & second_point,
//                lolita::integer direction = -1
//        )
//        const
//        {
//            if constexpr (t_ElementTraits::isCell()) {
//                auto const & current_coordinates = this->getCurrentCoordinates();
//                auto distance = lolita::real();
//                auto mp0 = lolita::domain::Point();
//                auto mp1 = lolita::domain::Point();
//                for (lolita::index i = 0; i < t_element.dim_; ++i) {
//                    mp0(i) = FiniteElementGeometry2::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
//                    mp1(i) = FiniteElementGeometry2::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
//                }
//                direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
//                return distance;
//            }
//            else {
//                auto const constexpr _quadrature = lolita::finite_element::Quadrature::Gauss();
//                using SegmentQuadrature = lolita::core2::geometry::ElementQuadratureTraits<lolita::core2::geometry::Element::LinearSegment(), _quadrature, 4>;
//                auto distance = lolita::real(0);
//                auto dt = lolita::real();
//                auto const current_nodes_coordinates = this->getCurrentCoordinates();
//                for (lolita::index q = 0; q < SegmentQuadrature::getDim(); ++q) {
//                    auto pq = SegmentQuadrature::reference_points_[q][0];
//                    auto wq = SegmentQuadrature::reference_weights_[q];
//                    auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>().setZero();
//                    auto difference = second_point - first_point;
//                    auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
//                    for (lolita::index i = 0; i < t_domain.dim_; ++i) {
//                        for (lolita::index j = 0; j < t_element.dim_; ++j) {
//                            if (direction == -1 || i == static_cast<lolita::index>(direction)) {
//                                auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
//                                auto dx = FiniteElementGeometry2::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
//                                ru(i, j) = dx * du;
//                            }
//                        }
//                    }
//                    if constexpr (t_element.isFacet()) {
//                        auto Eff = ru.col(0).template dot(ru.col(0));
//                        auto Fff = ru.col(0).template dot(ru.col(1));
//                        auto Gff = ru.col(1).template dot(ru.col(1));
//                        dt = std::sqrt(Eff + 2.0 * Fff + Gff);
//                    }
//                    else if constexpr (t_element.isCurve()) {
//                        auto Eff = ru.col(0).template dot(ru.col(0));
//                        dt = std::sqrt(Eff);
//                    }
//                    else {
//                        dt = 0;
//                    }
//                    distance += wq * dt;
//                }
//                return distance;
//            }
//        }
//
//        /*
//         * GEOMETRY
//         */
//
//        /**
//         * @brief
//         * @return
//         */
//        lolita::domain::Point
//        getCurrentDiameters()
//        const
//        {
////            auto reference_coordinates = FiniteElementGeometry2::getReferenceCoordinates();
//            auto const current_nodes_coordinates = this->getCurrentCoordinates();
//            auto current_diameters = lolita::domain::Point().setZero();
//            for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
//                for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
//                    auto const & pt0 = current_nodes_coordinates.col(i);
//                    auto const & pt1 = current_nodes_coordinates.col(j);
//                    for (lolita::index k = 0; k < 3; ++k) {
//                        auto new_value = numerics::abs(getRiemannianDistance(pt0, pt1, k));
//                        auto & current_value = current_diameters(k);
//                        if (new_value > current_value) {
//                            current_value = new_value;
//                        }
//                    }
//                }
//            }
//            return current_diameters;
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        lolita::domain::Point
//        getCurrentCentroid()
//        const
//        {
//            auto const current_nodes_coordinates = this->getCurrentCoordinates();
//            return lolita::domain::getBarycenter(current_nodes_coordinates);
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        static
//        lolita::domain::Point
//        getReferenceCentroid()
//        {
//            auto nds = FiniteElementGeometry2::getReferenceCoordinates();
//            return lolita::domain::getBarycenter(nds);
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        static
//        lolita::domain::Point
//        getReferenceDiameters()
//        {
////            auto dts = lolita::domain::Point();
////            auto nds = FiniteElementGeometry2::getReferenceCoordinates();
////            dts.setZero();
////            for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
////                for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
////                    for (lolita::index k = 0; k < 3; ++k) {
////                        auto new_value = numerics::abs(nds(k, i) - nds(k, j));
////                        auto & current_value = dts(k);
////                        if (new_value > current_value) {
////                            current_value = new_value;
////                        }
////                    }
////                }
////            }
////            return dts;
//            auto const reference_nodes_coordinates = FiniteElementGeometry2::getReferenceCoordinates();
////            auto const reference_nodes_coordinates = this->getCurrentCoordinates();
//            auto current_diameters = lolita::domain::Point().setZero();
//            for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
//                for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
//                    auto const & pt0 = reference_nodes_coordinates.col(i);
//                    auto const & pt1 = reference_nodes_coordinates.col(j);
//                    for (lolita::index k = 0; k < 3; ++k) {
//                        auto new_value = numerics::abs((pt0 - pt1)(k));
//                        auto & current_value = current_diameters(k);
//                        if (new_value > current_value) {
//                            current_value = new_value;
//                        }
//                    }
//                }
//            }
//            return current_diameters;
//        }
//
//        /*
//         * QUADRATURE
//         */
//
//        /**
//         * @brief
//         * @tparam _quadrature
//         * @tparam _ord
//         * @param index
//         * @return
//         */
//        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
//        static
//        lolita::real
//        getReferenceQuadratureWeight(
//                lolita::index index
//        )
//        {
//            return lolita::core2::geometry::ElementQuadratureTraits<t_element, _quadrature, _ord>::reference_weights_[index];
//        }
//
//        /**
//         * @brief
//         * @tparam _quadrature
//         * @tparam _ord
//         * @param index
//         * @return
//         */
//        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
//        static
//        lolita::matrix::Span<lolita::domain::Point const>
//        getReferenceQuadraturePoint(
//                lolita::index index
//        )
//        {
//            return lolita::matrix::Span<lolita::domain::Point const>(
//                    lolita::core2::geometry::ElementQuadratureTraits<t_element, _quadrature, _ord>::reference_points_[index].begin()
//            );
//        }
//
//        /**
//         * @brief
//         * @tparam _quadrature
//         * @tparam _ord
//         * @param index
//         * @return
//         */
//        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
//        lolita::real
//        getCurrentQuadratureWeight(
//                lolita::index index
//        )
//        const
//        {
//            auto w = getReferenceQuadratureWeight<_quadrature, _ord>(index);
//            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<_quadrature, _ord>(index));
//        }
//
//        /**
//         * @brief
//         * @tparam _quadrature
//         * @tparam _ord
//         * @param index
//         * @return
//         */
//        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
//        lolita::domain::Point
//        getCurrentQuadraturePoint(
//                lolita::index index
//        )
//        const
//        {
//            auto p = lolita::domain::Point();
//            auto const nds = this->getCurrentCoordinates();
//            for (lolita::index j = 0; j < 3; ++j) {
//                p(j) = FiniteElementGeometry2::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<_quadrature, _ord>(index));
//            }
//            return p;
//        }
//
//        /**
//         * @brief
//         * @tparam _quadrature
//         * @tparam _ord
//         * @tparam _i
//         * @tparam _j
//         * @param index
//         * @return
//         */
//        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
//        static
//        lolita::real
//        getComponentReferenceQuadratureWeight(
//                lolita::index index
//        )
//        requires(!t_element.isPoint())
//        {
//            auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
//            using ComponentGeometry = FiniteElementGeometry2<_component, t_domain, t_finite_element>;
//            return ComponentGeometry::template getReferenceQuadratureWeight<_quadrature, _ord>(index);
//        }
//
//        /**
//         * @brief
//         * @tparam _quadrature
//         * @tparam _ord
//         * @tparam _i
//         * @tparam _j
//         * @param component_index
//         * @param index
//         * @return
//         */
//        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
//        static
//        lolita::domain::Point
//        getComponentReferenceQuadraturePoint(
//                lolita::index component_index,
//                lolita::index index
//        )
//        requires(!t_element.isPoint())
//        {
////            auto const constexpr _component = lolita::core::component<_element, _domain, _i, _j>();
//            auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
//            auto p = lolita::domain::Point();
//            using ComponentGeometry = FiniteElementGeometry2<_component, t_domain, t_finite_element>;
//            auto const & elt_reference_nodes = lolita::core2::geometry::ElementGeometryTraits<t_element, t_domain>::reference_nodes_;
//            for (lolita::index i = 0; i < 3; ++i) {
//                auto cpt_coordinates = lolita::matrix::Vector<lolita::real, _component.num_nodes_>();
//                for (lolita::index j = 0; j < _component.num_nodes_; ++j) {
//                    auto const node_tag = getComponentNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
//                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
//                }
//                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<_quadrature, _ord>(index);
//                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
//            }
//            return p;
//        }
//
//        /**
//         * @brief
//         * @tparam _quadrature
//         * @tparam _ord
//         * @tparam _i
//         * @tparam _j
//         * @param component_index
//         * @param index
//         * @return
//         */
//        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
//        lolita::real
//        getComponentCurrentQuadratureWeight(
//                lolita::index component_index,
//                lolita::index index
//        )
//        const
//        requires(!t_element.isPoint())
//        {
//            auto const & cmp =  connectivity_->template getComponents<_i, _j>()[component_index];
//            return cmp->template getCurrentQuadratureWeight<_quadrature, _ord>(index);
//        }
//
//        /**
//         * @brief
//         * @tparam _quadrature
//         * @tparam _ord
//         * @tparam _i
//         * @tparam _j
//         * @param component_index
//         * @param index
//         * @return
//         */
//        template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
//        lolita::domain::Point
//        getComponentCurrentQuadraturePoint(
//                lolita::index component_index,
//                lolita::index index
//        )
//        const
//        requires(!t_element.isPoint())
//        {
//            auto p = lolita::domain::Point();
//            auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<_quadrature, _ord, _i, _j>(component_index, index);
//            auto const nds = this->getCurrentCoordinates();
//            for (lolita::index j = 0; j < 3; ++j) {
//                p(j) = FiniteElementGeometry2::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
//            }
//            return p;
//        }
//
//        /**
//         * @brief
//         * @param point
//         * @return
//         */
//        lolita::domain::Point
//        getNormalVector(
//                lolita::domain::Point const & point
//        )
//        const
//        requires(t_ElementTraits::isFace())
//        {
//            auto const current_nodes_coordinates = this->getCurrentCoordinates();
//            auto ru = lolita::matrix::Matrix<lolita::real, 3, t_element.dim_>();
//            ru.setZero();
//            for (lolita::index i = 0; i < 3; ++i) {
//                for (lolita::index j = 0; j < t_element.dim_; ++j) {
//                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
//                }
//            }
//            if constexpr (t_element.isPoint()) {
//                return lolita::domain::Point{0, 0, 0};
//            }
//            else if constexpr (t_element.isCurve()) {
//                return lolita::domain::Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
//            }
//            else {
//                return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
//            }
//        }
//
//        /**
//         * @brief
//         * @param point
//         * @param direction
//         * @return
//         */
//        lolita::domain::Point
//        getTangentVector(
//                lolita::domain::Point const & point,
//                lolita::index direction
//        )
//        const
//        requires(t_element.isCurve())
//        {
//            auto const current_nodes_coordinates = this->getCurrentCoordinates();
//            auto tangent_vector = lolita::domain::Point();
//            for (lolita::index i = 0; i < 3; ++i) {
//                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
//            }
//            return tangent_vector;
//        }
//
//        /**
//         * @brief
//         */
//        std::shared_ptr<FiniteElementConnectivity<FiniteElement, t_element, t_domain, t_finite_element>> connectivity_;
//
//    };
//
//    /**
//     * @brief
//     * @tparam t_element
//     * @tparam t_domain
//     * @tparam t_finite_element
//     */
//    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//    struct FiniteElementLoad2
//    {
//
//    private:
//
//        /**
//         * @brief
//         */
//        using t_Field = lolita::core2::field::TensorPolicy<t_finite_element.unknown_.tensor_, t_domain.dim_>;
//
//        /**
//         * @brief loads initializer
//         * @return loads member initialized with the zero natural load
//         */
//        static
//        std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, t_Field::shape_.rows_>, t_Field::shape_.cols_>
//        def_loads()
//        {
//            auto loads = std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, t_Field::shape_.rows_>, t_Field::shape_.cols_>();
//            for (int i = 0; i < t_Field::shape_.rows_; ++i) {
//                for (int j = 0; j < t_Field::shape_.cols_; ++j) {
//                    loads[j][i] = lolita::core2::finite_element::null_load_ptr;
//                }
//            }
//            return loads;
//        }
//
//    public:
//
//        /**
//         * @brief
//         * @param row
//         * @param col
//         * @return
//         */
//        std::shared_ptr<lolita::finite_element::LoadComponent> const &
//        getLoad(
//                lolita::integer row,
//                lolita::integer col
//        )
//        const
//        {
//            return loads_[col][row];
//        }
//
//        /**
//         * @brief
//         * @param row
//         * @param col
//         * @return
//         */
//        std::shared_ptr<lolita::finite_element::LoadComponent> &
//        getLoad(
//                lolita::integer row,
//                lolita::integer col
//        )
//        {
//            return loads_[col][row];
//        }
//
//        /**
//         * @brief
//         * @tparam t_element_group
//         * @param mesh
//         */
//        template<auto t_element_group>
//        void
//        setLoads(
//                lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
//        )
//        {
//            for (auto const & load : mesh.loads_) {
//                auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & domain) {
//                    return * domain == load.domain_tag_;
//                };
//                auto const has_domain = std::find_if(connectivity_->domains_.begin(), connectivity_->domains_.end(), is_equal) != connectivity_->domains_.end();
//                auto const has_unknown = load.unknown_tag_ == t_finite_element.unknown_.tensor_;
//                auto const has_dimension = load.element_dim_ == t_element.dim_;
//                if (has_domain && has_unknown && has_dimension) {
//                    std::cout << "setting load : " << load.domain_tag_ << std::endl;
//                    auto const i = load.components_.row_;
//                    auto const j = load.components_.col_;
//                    getLoad(i, j) = load.load_;
//                }
//            }
//        }
//
//        /**
//         * @brief
//         */
//        std::shared_ptr<FiniteElementConnectivity<FiniteElement, t_element, t_domain, t_finite_element>> connectivity_;
//
//        /**
//         * @brief
//         */
//        std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, t_Field::shape_.rows_>, t_Field::shape_.cols_> loads_ = def_loads();
//
//    };
//
//    /**
//     * @brief
//     * @tparam t_element
//     * @tparam t_domain
//     * @tparam t_finite_element
//     */
//    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
//    struct FiniteElementBehaviour2
//    {
//
//        /**
//         * @brief
//         * @tparam t_element_group
//         * @param mesh
//         */
//        template<auto t_element_group>
//        void
//        setBehaviour(
//                lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
//        )
//        {
//            for (auto const & behaviour : mesh.behaviours_) {
//                auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & domain) {
//                    return * domain == behaviour.domain_tag_;
//                };
//                auto const has_domain = std::find_if(connectivity_->domains_.begin(), connectivity_->domains_.end(), is_equal) != connectivity_->domains_.end();
//                auto const has_unknown = behaviour.unknown_tag_ == t_finite_element.unknown_.tensor_;
//                if (has_domain && has_unknown) {
//                    std::cout << "setting behaviour : " << behaviour.domain_tag_ << std::endl;
//                    behaviour_ = behaviour.behaviour_data_;
//                }
//            }
//        }
//
//        /**
//         * @brief
//         */
//        std::shared_ptr<FiniteElementConnectivity<FiniteElement, t_element, t_domain, t_finite_element>> connectivity_;
//
//        /**
//         * @brief
//         */
//        std::shared_ptr<lolita::behaviour::MgisBehaviourData> behaviour_;
//
//    };
//
//}

/**
 * @brief
 * @tparam t_T
 * @tparam t_element
 * @tparam t_domain
 * @tparam t_finite_element
 */
//template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, auto t_finite_element>
//struct FiniteElementBaseM
//{
//
//private:
//
//    /**
//     * @brief
//     */
//    template<lolita::core2::geometry::Element t__element, lolita::domain::Domain t__domain, auto... t__finite_element>
//    using t_FiniteElementPointer = std::shared_ptr<FiniteElement<t__element, t__domain, t__finite_element...>>;
//
//    /**
//     * @brief
//     */
//    using t_ElementTraits = lolita::core2::geometry::ElementTraits<t_element, t_domain>;
//
//    /**
//     * @brief
//     */
//    using t_FiniteElementTraits = lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>;
//
//    /**
//     * @brief
//     */
//    using t_FiniteElement = FiniteElement<t_element, t_domain, t_finite_element>;
//
//    /**
//     * @brief
//     */
//    using t_Quadrature = typename t_FiniteElementTraits::Quadrature;
//
//    /**
//     * @brief
//     */
//    using t_Field = lolita::core2::field::TensorPolicy<t_finite_element.unknown_.tensor_, t_domain.dim_>;
//
//    /**
//     * @brief loads initializer
//     * @return loads member initialized with the zero natural load
//     */
//    static
//    std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_>
//    def_loads()
//    {
//        auto loads = std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_>();
//        for (int i = 0; i < t_Field::shape_.rows_; ++i) {
//            for (int j = 0; j < t_Field::shape_.cols_; ++j) {
//                loads[j][i] = lolita::core2::finite_element::null_load_ptr;
//            }
//        }
//        return loads;
//    }
//
//public:
//
//    /**
//     * @brief
//     */
//    using Components = typename t_ElementTraits::template InnerConnectivity<t_FiniteElementPointer, t_finite_element>;
//
//    /**
//     * @brief
//     */
//    using Neighbours = typename t_ElementTraits::template OuterConnectivity<t_FiniteElementPointer, t_finite_element>;
//
//    /**
//     * @brief Basic equality comparison operator
//     * @param other
//     * @return
//     */
//    lolita::boolean
//    operator==(
//            FiniteElementBaseM const & other
//    )
//    const = default;
//
//    /**
//     * @brief Basic equality comparison operator
//     * @param other
//     * @return
//     */
//    lolita::boolean
//    operator!=(
//            FiniteElementBaseM const & other
//    )
//    const = default;
//
//    /**
//     * @brief
//     * @return
//     */
//    std::basic_string<lolita::character>
//    hash()
//    const
//    requires(t_element.isPoint())
//    {
//        return std::to_string(this->tag_);
//    }
//
//    /**
//     * @brief
//     * @return
//     */
//    std::basic_string<lolita::character>
//    hash()
//    const
//    {
//        std::basic_stringstream<lolita::character> hash;
//        auto const & nodes = getComponents<t_element.dim_ - 1, 0>();
//        for (auto const & node : nodes) {
//            hash << node->hash();
//        }
//        return hash.str();
//    }
//
//    /**
//     * @brief
//     * @tparam t_i
//     * @tparam t_j
//     * @return
//     */
//    template<lolita::index t_i, lolita::index t_j>
//    std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> &
//    getNeighbours()
//    {
//        return std::get<t_j>(std::get<t_i>(neighbours_));
//    }
//
//    /**
//     * @brief
//     * @tparam t_i
//     * @tparam t_j
//     * @return
//     */
//    template<lolita::index t_i, lolita::index t_j>
//    std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> const &
//    getNeighbours()
//    const
//    {
//        return std::get<t_j>(std::get<t_i>(neighbours_));
//    }
//
//    /**
//     * @brief
//     * @tparam t_i
//     * @tparam t_j
//     * @return
//     */
//    template<lolita::index t_i, lolita::index t_j>
//    std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> &
//    getComponents()
//    requires(!t_element.isPoint())
//    {
//        return std::get<t_j>(std::get<t_i>(components_));
//    }
//
//    /**
//     * @brief
//     * @tparam t_i
//     * @tparam t_j
//     * @return
//     */
//    template<lolita::index t_i, lolita::index t_j>
//    std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> const &
//    getComponents()
//    const
//    requires(!t_element.isPoint())
//    {
//        return std::get<t_j>(std::get<t_i>(components_));
//    }
//
//    /**
//     * @brief
//     * @tparam t_i
//     * @tparam t_j
//     * @param i
//     * @param j
//     * @return
//     */
//    template<lolita::index t_i, lolita::index t_j>
//    static constexpr
//    lolita::index
//    getComponentNodeConnection(
//            lolita::index i,
//            lolita::index j
//    )
//    requires(!t_element.isPoint())
//    {
//        return std::get<t_j>(std::get<t_i>(t_ElementTraits::node_connectivity_))[i][j];
//    }
//
//    /**
//     * @brief
//     * @tparam t_i
//     * @tparam t_j
//     * @param i
//     * @return
//     */
//    template<lolita::index t_i, lolita::index t_j>
//    lolita::index
//    getComponentIndex(
//            lolita::index i
//    )
//    const
//    requires(!t_element.isPoint())
//    {
//        auto constexpr t_component = t_ElementTraits::template getComponent<t_i, t_j>();
//        using t_NeighbourTraits = lolita::core2::geometry::ElementTraits<t_component, t_domain>;
//        auto constexpr t_coordinates = t_NeighbourTraits::template getNeighbourCoordinates<t_element>();
//        auto const & items = getComponents<t_i, t_j>()[i]->template getNeighbours<t_coordinates.dim_, t_coordinates.tag_>();
//        auto is_equal = [&] (std::shared_ptr<t_FiniteElement> const & ptr_element) {
//            return * ptr_element == * this;
//        };
//        return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
//    }
//
//    /**
//     * @brief
//     * @tparam t_i
//     * @tparam t_j
//     * @param i
//     * @return
//     */
//    template<lolita::index t_i, lolita::index t_j>
//    lolita::integer
//    getComponentOrientation(
//            lolita::index i
//    )
//    const
//    requires(!t_element.isPoint())
//    {
//        return getComponentIndex<t_i, t_j>(i) == 0 ? 1 : -1;
//    }
//
//    /**
//     * @brief Get the current coordinates of the point/node in the mesh
//     * @return the vector with x, y and z coordinates
//     */
//    lolita::domain::Point &
//    getCurrentCoordinates()
//    requires(t_element.isPoint())
//    {
//        return * this->coordinates_;
//    }
//
//    /**
//     * @brief Get the current coordinates of the point/node in the mesh
//     * @return the vector with x, y and z coordinates
//     */
//    lolita::domain::Point const &
//    getCurrentCoordinates()
//    const
//    requires(t_element.isPoint())
//    {
//        return * this->coordinates_;
//    }
//
//    /**
//     * @brief Get the current coordinates of the element in the mesh
//     * @return the matrix composed of the column vectors with x, y and z coordinates of each node
//     */
//    lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>
//    getCurrentCoordinates()
//    const
//    {
//        auto current_nodes_coordinates = lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>();
//        auto count = lolita::index(0);
//        for (auto const & node : this->template getComponents<t_element.dim_ - 1, 0>()) {
//            current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
//            count ++;
//        }
//        return current_nodes_coordinates;
//    }
//
//    /**
//     * @brief
//     * @return
//     */
//    static
//    lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.numNodes(), matrix::col_major> const>
//    getReferenceCoordinates()
//    {
//        using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.numNodes(), matrix::col_major> const>;
//        return _ReferenceCoordinates(t_ElementTraits::reference_nodes_.begin()->begin());
//    }
//
////        /**
////         * @brief
////         * @tparam _basis
////         * @tparam _ord
////         * @param point
////         * @return
////         */
////        template<lolita::core2::finite_element::basis::Basis _basis, lolita::index _ord>
////        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, _basis, _ord>::dim_>
////        getBasisEvaluation(
////                lolita::domain::Point const & point
////        )
////        const
////        {
////            using t_FiniteElementBasis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, _basis, _ord>;
////            using t_FiniteElementBasisImplementation = typename t_FiniteElementBasis::template Implementation<t_domain, t_finite_element>;
////            return static_cast<t_FiniteElementBasisImplementation const *>(this)->evaluate(point);
////        }
////
////        /**
////         * @brief
////         * @tparam _basis
////         * @tparam _ord
////         * @param point
////         * @param derivative_direction
////         * @return
////         */
////        template<lolita::core2::finite_element::basis::Basis _basis, lolita::index _ord>
////        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, _basis, _ord>::dim_>
////        getBasisDerivative(
////                lolita::domain::Point const & point,
////                lolita::index derivative_direction
////        )
////        const
////        {
////            using t_FiniteElementBasis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, _basis, _ord>;
////            using t_FiniteElementBasisImplementation = typename t_FiniteElementBasis::template Implementation<t_domain, t_finite_element>;
////            return static_cast<t_FiniteElementBasisImplementation const *>(this)->evaluate(point, derivative_direction);
////        }
//
//    /**
//     * @brief
//     * @param nodal_field_values
//     * @param reference_point
//     * @return
//     */
//    static
//    lolita::real
//    getShapeMappingEvaluation(
//            lolita::matrix::Vector<lolita::real, t_element.numNodes()> const & nodal_field_values,
//            lolita::domain::Point const & reference_point
//    )
//    {
//        return t_ElementTraits::getShapeMappingEvaluation(nodal_field_values, reference_point);
//    }
//
//    /**
//     * @brief
//     * @param nodal_field_values
//     * @param reference_point
//     * @param derivative_direction
//     * @return
//     */
//    static
//    lolita::real
//    getShapeMappingDerivative(
//            lolita::matrix::Vector<lolita::real, t_element.numNodes()> const & nodal_field_values,
//            lolita::domain::Point const & reference_point,
//            lolita::index derivative_direction
//    )
//    {
//        return t_ElementTraits::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
//    }
//
//    /*
//     *
//     */
//
//    /**
//     * @brief
//     * @param point
//     * @return
//     */
//    lolita::real
//    getShapeMappingDifferential(
//            lolita::domain::Point const & point
//    )
//    const
//    {
//        auto const current_coordinates = this->getCurrentCoordinates();
//        auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>();
//        auto du = lolita::real(0);
//        ru.setZero();
//        for (lolita::index i = 0; i < t_domain.dim_; ++i) {
//            for (lolita::index j = 0; j < t_element.dim_; ++j) {
//                ru(i, j) = FiniteElementBaseM::getShapeMappingDerivative(current_coordinates.row(i), point, j);
//            }
//        }
//        if constexpr (t_element.dim_ == 3) {
//            du = numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
//        }
//        else if constexpr (t_element.dim_ == 2) {
//            du = numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
//        }
//        else {
//            du = numerics::abs(ru.col(0).norm());
//        }
//        if constexpr (t_domain.frame_ == lolita::domain::Frame::AxiSymmetric()) {
//            lolita::real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
//            if (r0 < 1.e-10) {
//                r0 = 1.e-10;
//            }
//            du *= 2.0 * numerics::pi * r0;
//        }
//        return du;
//    }
//
//    /**
//     * @brief
//     * @param first_point
//     * @param second_point
//     * @param direction
//     * @return
//     */
//    lolita::real
//    getRiemannianDistance(
//            lolita::domain::Point const & first_point,
//            lolita::domain::Point const & second_point,
//            lolita::integer direction = -1
//    )
//    const
//    {
//        if constexpr (t_ElementTraits::isCell()) {
//            auto const & current_coordinates = this->getCurrentCoordinates();
//            auto distance = lolita::real();
//            auto mp0 = lolita::domain::Point();
//            auto mp1 = lolita::domain::Point();
//            for (lolita::index i = 0; i < t_element.dim_; ++i) {
//                mp0(i) = FiniteElementBaseM::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
//                mp1(i) = FiniteElementBaseM::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
//            }
//            direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
//            return distance;
//        }
//        else {
//            auto const constexpr _quadrature = lolita::finite_element::Quadrature::Gauss();
//            using SegmentQuadrature = lolita::core2::geometry::ElementQuadratureTraits<lolita::core2::geometry::Element::LinearSegment(), _quadrature, 4>;
//            auto distance = lolita::real(0);
//            auto dt = lolita::real();
//            auto const current_nodes_coordinates = this->getCurrentCoordinates();
//            for (lolita::index q = 0; q < SegmentQuadrature::getDim(); ++q) {
//                auto pq = SegmentQuadrature::reference_points_[q][0];
//                auto wq = SegmentQuadrature::reference_weights_[q];
//                auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>().setZero();
//                auto difference = second_point - first_point;
//                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
//                for (lolita::index i = 0; i < t_domain.dim_; ++i) {
//                    for (lolita::index j = 0; j < t_element.dim_; ++j) {
//                        if (direction == -1 || i == static_cast<lolita::index>(direction)) {
//                            auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
//                            auto dx = FiniteElementBaseM::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
//                            ru(i, j) = dx * du;
//                        }
//                    }
//                }
//                if constexpr (t_element.isFacet()) {
//                    auto Eff = ru.col(0).template dot(ru.col(0));
//                    auto Fff = ru.col(0).template dot(ru.col(1));
//                    auto Gff = ru.col(1).template dot(ru.col(1));
//                    dt = std::sqrt(Eff + 2.0 * Fff + Gff);
//                }
//                else if constexpr (t_element.isCurve()) {
//                    auto Eff = ru.col(0).template dot(ru.col(0));
//                    dt = std::sqrt(Eff);
//                }
//                else {
//                    dt = 0;
//                }
//                distance += wq * dt;
//            }
//            return distance;
//        }
//    }
//
//    /*
//     * GEOMETRY
//     */
//
//    /**
//     * @brief
//     * @return
//     */
//    lolita::domain::Point
//    getCurrentDiameters()
//    const
//    {
//        auto reference_coordinates = FiniteElementBaseM::getReferenceCoordinates();
//        auto current_diameters = lolita::domain::Point().setZero();
//        for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
//            for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
//                auto const & pt0 = reference_coordinates.col(i);
//                auto const & pt1 = reference_coordinates.col(j);
//                for (lolita::index k = 0; k < 3; ++k) {
//                    auto new_value = numerics::abs(getRiemannianDistance(pt0, pt1, k));
//                    auto & current_value = current_diameters(k);
//                    if (new_value > current_value) {
//                        current_value = new_value;
//                    }
//                }
//            }
//        }
//        return current_diameters;
//    }
//
//    /**
//     * @brief
//     * @return
//     */
//    lolita::domain::Point
//    getCurrentCentroid()
//    const
//    {
//        auto const current_nodes_coordinates = this->getCurrentCoordinates();
//        return lolita::domain::getBarycenter(current_nodes_coordinates);
//    }
//
//    /**
//     * @brief
//     * @return
//     */
//    static
//    lolita::domain::Point
//    getReferenceCentroid()
//    {
//        auto nds = FiniteElementBaseM::getReferenceCoordinates();
//        return lolita::domain::getBarycenter(nds);
//    }
//
//    /**
//     * @brief
//     * @return
//     */
//    static
//    lolita::domain::Point
//    getReferenceDiameters()
//    {
//        auto dts = lolita::domain::Point();
//        auto nds = FiniteElementBaseM::getReferenceCoordinates();
//        dts.setZero();
//        for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
//            for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
//                for (lolita::index k = 0; k < 3; ++k) {
//                    auto new_value = numerics::abs(nds(k, i) - nds(k, j));
//                    auto & current_value = dts(k);
//                    if (new_value > current_value) {
//                        current_value = new_value;
//                    }
//                }
//            }
//        }
//        return dts;
//    }
//
//    /*
//     * QUADRATURE
//     */
//
//    /**
//     * @brief
//     * @tparam _quadrature
//     * @tparam _ord
//     * @param index
//     * @return
//     */
//    template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
//    static
//    lolita::real
//    getReferenceQuadratureWeight(
//            lolita::index index
//    )
//    {
//        return lolita::core2::geometry::ElementQuadratureTraits<t_element, _quadrature, _ord>::reference_weights_[index];
//    }
//
//    /**
//     * @brief
//     * @tparam _quadrature
//     * @tparam _ord
//     * @param index
//     * @return
//     */
//    template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
//    static
//    lolita::matrix::Span<lolita::domain::Point const>
//    getReferenceQuadraturePoint(
//            lolita::index index
//    )
//    {
//        return lolita::matrix::Span<lolita::domain::Point const>(
//                lolita::core2::geometry::ElementQuadratureTraits<t_element, _quadrature, _ord>::reference_points_[index].begin()
//        );
//    }
//
//    /**
//     * @brief
//     * @tparam _quadrature
//     * @tparam _ord
//     * @param index
//     * @return
//     */
//    template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
//    lolita::real
//    getCurrentQuadratureWeight(
//            lolita::index index
//    )
//    const
//    {
//        auto w = getReferenceQuadratureWeight<_quadrature, _ord>(index);
//        return w * getShapeMappingDifferential(getReferenceQuadraturePoint<_quadrature, _ord>(index));
//    }
//
//    /**
//     * @brief
//     * @tparam _quadrature
//     * @tparam _ord
//     * @param index
//     * @return
//     */
//    template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_>
//    lolita::domain::Point
//    getCurrentQuadraturePoint(
//            lolita::index index
//    )
//    const
//    {
//        auto p = lolita::domain::Point();
//        auto const nds = this->getCurrentCoordinates();
//        for (lolita::index j = 0; j < 3; ++j) {
//            p(j) = FiniteElementBaseM::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<_quadrature, _ord>(index));
//        }
//        return p;
//    }
//
//    /**
//     * @brief
//     * @tparam _quadrature
//     * @tparam _ord
//     * @tparam _i
//     * @tparam _j
//     * @param index
//     * @return
//     */
//    template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
//    static
//    lolita::real
//    getComponentReferenceQuadratureWeight(
//            lolita::index index
//    )
//    requires(!t_element.isPoint())
//    {
//        auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
//        using ComponentGeometry = FiniteElementBaseM<_component, t_domain, t_finite_element>;
//        return ComponentGeometry::template getReferenceQuadratureWeight<_quadrature, _ord>(index);
//    }
//
//    /**
//     * @brief
//     * @tparam _quadrature
//     * @tparam _ord
//     * @tparam _i
//     * @tparam _j
//     * @param component_index
//     * @param index
//     * @return
//     */
//    template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
//    static
//    lolita::domain::Point
//    getComponentReferenceQuadraturePoint(
//            lolita::index component_index,
//            lolita::index index
//    )
//    requires(!t_element.isPoint())
//    {
////            auto const constexpr _component = lolita::core::component<_element, _domain, _i, _j>();
//        auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
//        auto p = lolita::domain::Point();
//        using ComponentGeometry = FiniteElementBaseM<_component, t_domain, t_finite_element>;
//        auto const & elt_reference_nodes = lolita::core2::geometry::ElementTraits<t_element, t_domain>::reference_nodes_;
//        for (lolita::index i = 0; i < 3; ++i) {
//            auto cpt_coordinates = lolita::matrix::Vector<lolita::real, _component.num_nodes_>();
//            for (lolita::index j = 0; j < _component.num_nodes_; ++j) {
//                auto const node_tag = getComponentNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
//                cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
//            }
//            auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<_quadrature, _ord>(index);
//            p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
//        }
//        return p;
//    }
//
//    /**
//     * @brief
//     * @tparam _quadrature
//     * @tparam _ord
//     * @tparam _i
//     * @tparam _j
//     * @param component_index
//     * @param index
//     * @return
//     */
//    template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
//    lolita::real
//    getComponentCurrentQuadratureWeight(
//            lolita::index component_index,
//            lolita::index index
//    )
//    const
//    requires(!t_element.isPoint())
//    {
//        auto const & cmp =  this->template getComponents<_i, _j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
//        return cmp->template getCurrentQuadratureWeight<_quadrature, _ord>(index);
//    }
//
//    /**
//     * @brief
//     * @tparam _quadrature
//     * @tparam _ord
//     * @tparam _i
//     * @tparam _j
//     * @param component_index
//     * @param index
//     * @return
//     */
//    template<lolita::finite_element::Quadrature _quadrature = t_Quadrature::quadrature_, lolita::index _ord = t_Quadrature::ord_, lolita::index _i, lolita::index _j>
//    lolita::domain::Point
//    getComponentCurrentQuadraturePoint(
//            lolita::index component_index,
//            lolita::index index
//    )
//    const
//    requires(!t_element.isPoint())
//    {
//        auto p = lolita::domain::Point();
//        auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<_quadrature, _ord, _i, _j>(component_index, index);
//        auto const nds = this->getCurrentCoordinates();
//        for (lolita::index j = 0; j < 3; ++j) {
//            p(j) = FiniteElementBaseM::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
//        }
//        return p;
//    }
//
//    /**
//     * @brief
//     * @param point
//     * @return
//     */
//    lolita::domain::Point
//    getNormalVector(
//            lolita::domain::Point const & point
//    )
//    const
//    requires(t_ElementTraits::isFace())
//    {
//        auto const current_nodes_coordinates = this->getCurrentCoordinates();
//        auto ru = lolita::matrix::Matrix<lolita::real, 3, t_element.dim_>();
//        ru.setZero();
//        for (lolita::index i = 0; i < 3; ++i) {
//            for (lolita::index j = 0; j < t_element.dim_; ++j) {
//                ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
//            }
//        }
//        if constexpr (t_element.isPoint()) {
//            return lolita::domain::Point{0, 0, 0};
//        }
//        else if constexpr (t_element.isCurve()) {
//            return lolita::domain::Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
//        }
//        else {
//            return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
//        }
//    }
//
//    /**
//     * @brief
//     * @param point
//     * @param direction
//     * @return
//     */
//    lolita::domain::Point
//    getTangentVector(
//            lolita::domain::Point const & point,
//            lolita::index direction
//    )
//    const
//    requires(t_element.isCurve())
//    {
//        auto const current_nodes_coordinates = this->getCurrentCoordinates();
//        auto tangent_vector = lolita::domain::Point();
//        for (lolita::index i = 0; i < 3; ++i) {
//            tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
//        }
//        return tangent_vector;
//    }
//
//    /**
//     * @brief
//     * @param row
//     * @param col
//     * @return
//     */
//    std::shared_ptr<lolita::finite_element::ScalarLoad> const &
//    getLoad(
//            lolita::integer row,
//            lolita::integer col
//    )
//    const
//    {
//        return loads_[col][row];
//    }
//
//    /**
//     * @brief
//     * @param row
//     * @param col
//     * @return
//     */
//    std::shared_ptr<lolita::finite_element::ScalarLoad> &
//    getLoad(
//            lolita::integer row,
//            lolita::integer col
//    )
//    {
//        return loads_[col][row];
//    }
//
//    /**
//     * @brief
//     * @tparam t_element_group
//     * @param mesh
//     */
//    template<auto t_element_group>
//    void
//    setLoads(
//            lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
//    )
//    {
//        for (auto const & load : mesh.loads_) {
//            auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & domain) {
//                return * domain == load.domain_tag_;
//            };
//            auto const has_domain = std::find_if(domains_.begin(), domains_.end(), is_equal) != domains_.end();
//            auto const has_unknown = load.unknown_tag_ == t_finite_element.unknown_.tensor_;
//            auto const has_dimension = load.element_dim_ == t_element.dim_;
//            if (has_domain && has_unknown && has_dimension) {
//                std::cout << "setting load : " << load.domain_tag_ << std::endl;
//                auto const i = load.components_.row_;
//                auto const j = load.components_.col_;
//                getLoad(i, j) = load.load_;
//            }
//        }
//    }
//
//    /**
//     * @brief
//     * @tparam t_element_group
//     * @param mesh
//     */
//    template<auto t_element_group>
//    void
//    setBehaviour(
//            lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
//    )
//    {
//        for (auto const & behaviour : mesh.behaviours_) {
//            auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & domain) {
//                return * domain == behaviour.domain_tag_;
//            };
//            auto const has_domain = std::find_if(domains_.begin(), domains_.end(), is_equal) != domains_.end();
//            auto const has_unknown = behaviour.unknown_tag_ == t_finite_element.unknown_.tensor_;
//            if (has_domain && has_unknown) {
//                std::cout << "setting behaviour : " << behaviour.domain_tag_ << std::endl;
//                behaviour_ = behaviour.behaviour_data_;
//            }
//        }
//    }
//
//    /**
//     * @brief
//     */
//    Neighbours neighbours_;
//
//    /**
//     * @brief
//     */
//    Components components_;
//
//    /**
//     * @brief
//     */
//    std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;
//
//    /**
//     * @brief
//     */
//    lolita::natural tag_;
//
//    /**
//     * @brief
//     */
//    std::shared_ptr<lolita::domain::Point> coordinates_;
//
//    /**
//     * @brief
//     */
//    std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_> loads_ = def_loads();
//
//    /**
//     * @brief
//     */
//    std::shared_ptr<lolita::behaviour::MgisBehaviourData> behaviour_;
//
//};
