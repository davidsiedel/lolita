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
            Unknowns()
            {
                return std::array<Unknown, 2>{
                        Unknown::Subsidiary(),
                        Unknown::Structural(),
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
             * @brief The field type
             */
            using Field = lolita::core2::field::TensorPolicy<t_field, t_domain.dim_>;

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
                //getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::LinSpaced(0, 0 + t_dim - 1));
            }

            /**
             * @brief
             * @param i
             * @param j
             * @param dim
             */
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim
            )
            requires(t_dim == -1 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim));
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
                //getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::LinSpaced(0, 0 + t_dim - 1));
                auto & k = mesh.systems_[_finite_element_index].num_unknowns_;
                getUnknownComponent(i, j).unknowns_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(k, k + t_dim - 1));
                k += t_dim;
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param dim
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim == -1 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim));
                auto & k = mesh.systems_[_finite_element_index].num_unknowns_;
                getUnknownComponent(i, j).unknowns_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(dim, k, k + dim - 1));
                k += dim;
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
             * @param dim
             */
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim
            )
            requires(t_dim == -1 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim));
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
                getUnknownComponent(i, j).bindings_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(k, k + t_dim - 1));
                k += t_dim;
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param dim
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim == -1 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim));
                auto & k = mesh.systems_[_finite_element_index].num_bindings_;
                getUnknownComponent(i, j).bindings_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(dim, k, k + dim - 1));
                k += dim;
            }

            /**
             * @brief
             */
            std::array<std::array<ScalarUnknown<t_unknown, t_dim>, Field::shape_.rows_>, Field::shape_.cols_> field_unknown_;

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
                return ((t_FieldUnknown::unknown_.isStructural()) || ...);
            }

            /**
             * @brief
             * @tparam t_unknown
             * @return
             */
            template<lolita::core2::finite_element::unknown::Unknown... t_unknown>
            static constexpr
            lolita::integer
            getDimUnknowns()
            {
                auto constexpr unknowns = std::array<lolita::core2::finite_element::unknown::Unknown, sizeof...(t_unknown)>{t_unknown...};
                auto dim_unknowns = lolita::integer(0);
                auto set_dim_unknowns  = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                    if constexpr (std::tuple_size_v<Unknowns> > 0) {
                        using t_Unknown = std::tuple_element_t<t_j, Unknowns>;
                        if constexpr (t_Unknown::unknown_ == unknowns[t_i]) {
                            dim_unknowns += t_Unknown::dim_;
                        }
                        if constexpr (t_j < std::tuple_size_v<Unknowns> - 1) {
                            self.template operator()<t_i, t_j + 1>(self);
                        }
                        else if constexpr (t_i < sizeof...(t_unknown) - 1) {
                            self.template operator()<t_i + 1, 0>(self);
                        }
                    }
                };
                set_dim_unknowns(set_dim_unknowns);
                return dim_unknowns;
            }

            /**
             * @brief
             * @tparam t_domain
             * @tparam t_unknown
             * @return
             */
            template<lolita::domain::Domain t_domain, lolita::core2::finite_element::unknown::Unknown... t_unknown>
            static constexpr
            lolita::integer
            getNumUnknowns()
            {
//                auto num_unknowns = lolita::integer(0);
//                auto set_num_unknowns  = [&] <lolita::integer t_i = 0> (auto & self) constexpr mutable {
//                    if constexpr (std::tuple_size_v<Unknowns> > 0) {
//                        using t_Unknown = std::tuple_element_t<t_i, Unknowns>;
//                        if constexpr (t_Unknown::unknown_ == t_unknown) {
//                            num_unknowns += t_Unknown::dim_ * t_Unknown::Field::shape_.size_;
//                        }
//                        if constexpr (t_i < std::tuple_size_v<Unknowns> - 1) {
//                            self.template operator()<t_i + 1>(self);
//                        }
//                    }
//                };
//                set_num_unknowns(set_num_unknowns);
//                return num_unknowns;
                auto constexpr unknowns = std::array<lolita::core2::finite_element::unknown::Unknown, sizeof...(t_unknown)>{t_unknown...};
                auto num_unknowns = lolita::integer(0);
                auto set_num_unknowns  = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                    if constexpr (std::tuple_size_v<Unknowns> > 0) {
                        using t_Unknown = std::tuple_element_t<t_j, Unknowns>;
                        if constexpr (t_Unknown::unknown_ == unknowns[t_i]) {
                            num_unknowns += t_Unknown::dim_ * t_Unknown::Field::shape_.size_;
                        }
                        if constexpr (t_j < std::tuple_size_v<Unknowns> - 1) {
                            self.template operator()<t_i, t_j + 1>(self);
                        }
                        else if constexpr (t_i < sizeof...(t_unknown) - 1) {
                            self.template operator()<t_i + 1, 0>(self);
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
            using t_Field = typename lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>::Field;

        public:

            /**
             * @brief
             */
            using Basis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<
                    t_element,
                    lolita::core2::finite_element::basis::Basis::Monomial(),
                    t_finite_element.discretization_.ord_cell_
            >;

            /**
             * @brief
             */
            using FieldUnknowns = lolita::core2::finite_element::unknown::FieldUnknownCollection<
                    lolita::core2::finite_element::unknown::FieldUnknown<
                            t_domain,
                            t_finite_element.unknown_.tensor_,
                            Basis::dim_,
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
            using t_Field = typename lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>::Field;

        public:

            /**
             * @brief
             */
            using Basis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<
                    t_element,
                    lolita::core2::finite_element::basis::Basis::Monomial(),
                    t_finite_element.discretization_.ord_face_
            >;

            /**
             * @brief
             */
            using FieldUnknowns = lolita::core2::finite_element::unknown::FieldUnknownCollection<
                    lolita::core2::finite_element::unknown::FieldUnknown<
                            t_domain,
                            t_finite_element.unknown_.tensor_,
                            Basis::dim_,
                            lolita::core2::finite_element::unknown::Unknown::Structural()
                    >
//                    ,lolita::core2::finite_element::unknown::FieldUnknown<
//                            t_domain,
//                            t_finite_element.unknown_.tensor_,
//                            Basis::dim_,
//                            lolita::core2::finite_element::unknown::Unknown::Structural()
//                    >
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
//                    for (int i = 0; i < t_Field::shape_.rows_; ++i) {
//                        for (int j = 0; j < t_Field::shape_.cols_; ++j) {
//                            this->field_unknowns_.template getUnknown<1>().template setUnknownCoefficients<t_finite_element>(i, j, mesh);
//                            if (this->getLoadComponent(i, j)->load_.isConstraint()) {
//                                std::cout << "setting const " << std::endl;
//                                this->field_unknowns_.template getUnknown<1>().template setBindingCoefficients<t_finite_element>(i, j, mesh);
//                            }
//                        }
//                    }
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
    struct FiniteElementFieldUnknowns
    {

        /**
         * @brief
         */
        using t_FiniteElementFieldUnknownsTraits = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         */
        using t_FieldUnknowns = typename t_FiniteElementFieldUnknownsTraits::FieldUnknowns;

        /**
         * @brief
         */
        using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown... t_unknown>
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            return t_FieldUnknowns::template getNumUnknowns<t_domain, t_unknown...>();
        }

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown... t_unknown>
        static constexpr
        lolita::integer
        getDimUnknowns()
        {
            return t_FieldUnknowns::template getDimUnknowns<t_unknown...>();
        }

        /**
         * @brief
         */
        struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
        {

            /**
             * @brief
             * @tparam t_unknown
             * @return
             */
            template<lolita::core2::finite_element::unknown::Unknown t_unknown>
            lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns<t_unknown>()>
            getUnknowns()
            const
            {
                using t_UnknownVector = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns<t_unknown>()>;
                auto unknowns = lolita::matrix::Zero<t_UnknownVector>();
                auto count = 0;
                auto set_element_unknowns = [&] <lolita::integer t_k = 0> (auto & self) constexpr mutable {
                    using t_ElementUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>::FieldUnknowns;
                    if constexpr (std::tuple_size_v<typename t_ElementUnknowns::Unknowns> > 0) {
                        using t_ElementUnknown = std::tuple_element_t<t_k, typename t_ElementUnknowns::Unknowns>;
                        if constexpr (t_ElementUnknown::unknown_ == t_unknown) {
                            for (int i = 0; i < t_ElementUnknown::Field::shape_.rows_; ++i) {
                                for (int j = 0; j < t_ElementUnknown::Field::shape_.cols_; ++j) {
                                    auto const & rhs = this->field_unknowns_.template getUnknown<t_k>().getUnknownComponent(i, j).getUnknownsCoefficients();
                                    unknowns.template segment<t_ElementUnknown::dim_>(count) = rhs;
                                    count += t_ElementUnknown::dim_;
                                }
                            }
                        }
                    }
                    if constexpr (t_k < lolita::integer(std::tuple_size_v<typename t_ElementUnknowns::Unknowns>) - 1) {
                        self.template operator()<t_k + 1>(self);
                    }
                };
                set_element_unknowns(set_element_unknowns);
                auto set_neighbours_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0, lolita::integer t_k = 0> (auto & self) constexpr mutable {
                    //std::cout << t_element << " " << "t_i : " << t_i << " " << "t_j : " << t_j << " " << "t_k : " << t_k << std::endl;
                    auto constexpr t_face = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<t_i, t_j>();
                    using t_FaceUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<t_face, t_domain, t_finite_element>::FieldUnknowns;
                    if constexpr (std::tuple_size_v<typename t_FaceUnknowns::Unknowns> > 0) {
                        using t_FaceUnknown = std::tuple_element_t<t_k, typename t_FaceUnknowns::Unknowns>;
                        if constexpr (t_FaceUnknown::unknown_ == t_unknown) {
                            auto const & faces = this->template getComponents<t_i, t_j>();
                            for (auto const & face : faces) {
                                for (int i = 0; i < t_FaceUnknown::Field::shape_.rows_; ++i) {
                                    for (int j = 0; j < t_FaceUnknown::Field::shape_.cols_; ++j) {
                                        auto const & rhs = face->field_unknowns_.template getUnknown<t_k>().getUnknownComponent(i, j).getUnknownsCoefficients();
                                        unknowns.template segment<t_FaceUnknown::dim_>(count) = rhs;
                                        count += t_FaceUnknown::dim_;
                                    }
                                }
                            }
                        }
                    }
                    if constexpr (t_k < lolita::integer(std::tuple_size_v<typename t_FaceUnknowns::Unknowns>) - 1) {
                        self.template operator()<t_i, t_j, t_k + 1>(self);
                    }
                    else if constexpr (t_j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<t_i>() - 1) {
                        self.template operator()<t_i, t_j + 1, 0>(self);
                    }
                    else if constexpr (t_i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<>() - 1) {
                        self.template operator()<t_i + 1, 0, 0>(self);
                    }
                };
                if constexpr (!t_element.isPoint()) {
                    set_neighbours_unknowns(set_neighbours_unknowns);
                }
                //std::cout << t_element << " " << t_unknown << " unknowns :" << std::endl;
                //std::cout << unknowns << std::endl;
                return unknowns;
            }

            /**
             * @brief
             * @return
             */
            lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns()>
            getUnknowns()
            const
            {
                auto unknowns = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns()>();
                auto count = lolita::integer(0);
                auto set_unknown_block = [&] <lolita::integer t_i = 0> (auto & self) mutable {
                    auto constexpr t_unknown = lolita::core2::finite_element::unknown::Unknown::Unknowns()[t_i];
                    auto constexpr num_unknowns = t_FiniteElementTraits::template getNumUnknowns<t_unknown>();
                    unknowns.template segment<num_unknowns>(count) = getUnknowns<t_unknown>();
                    count += num_unknowns;
                    if constexpr (t_i < lolita::core2::finite_element::unknown::Unknown::Unknowns().size() - 1) {
                        self.template operator()<t_i + 1>(self);
                    }
                };
                set_unknown_block(set_unknown_block);
                return unknowns;
            }

        };

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
            static_cast<typename t_FiniteElementFieldUnknownsTraits::Implementation *>(this)->setUnknowns(mesh);
        }

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown>
        lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns<t_unknown>()>
        getUnknowns()
        const
        {
            return static_cast<Implementation const *>(this)->template getUnknowns<t_unknown>();
        }

        /**
         * @brief
         * @return
         */
        lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns()>
        getUnknowns()
        const
        {
            return static_cast<Implementation const *>(this)->getUnknowns();
        }

        /**
         * @brief
         */
        t_FieldUnknowns field_unknowns_;

    };

    namespace face
    {

        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        struct FiniteElementFaceTraits;

        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        requires(t_finite_element.discretization_.isHHO())
        struct FiniteElementFaceTraits<t_element, t_domain, t_finite_element>
        {

            struct Face
            {};

            struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
            {

                template<auto t_element_group>
                void
                makeDirichlet(
                        lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
                )
                {

                }

            };

        };

    }

    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementFace
    {

    private:

        /**
         * @brief
         */
        using t_FiniteElementTraits = lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>;

    public:

        /**
         * @brief
         */
        lolita::integer static constexpr dim_structural_unknowns = t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Structural()>();

        /**
         * @brief
         */
        lolita::integer static constexpr dim_subsidiary_unknowns = t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Subsidiary()>();

        /**
         * @brief
         */
        lolita::integer static constexpr num_structural_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Structural()>();

        /**
         * @brief
         */
        lolita::integer static constexpr num_subsidiary_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Subsidiary()>();

        /**
         * @brief
         */
        struct IntegrationPoint
        {

            /**
             * @brief
             */
            lolita::domain::Point coordinates_;

            /**
             * @brief
             */
            lolita::real weight_;

            /**
             * @brief
             */
            std::unique_ptr<mgis::behaviour::BehaviourData> material_point_;

            /**
             * @brief
             */
            lolita::matrix::Vector<lolita::real, dim_structural_unknowns> structural_load_operator_;

            /**
             * @brief
             */
            lolita::matrix::Vector<lolita::real, dim_subsidiary_unknowns> subsidiary_load_operator_;

        };

        /**
         * @brief
         */
        std::array<IntegrationPoint, t_FiniteElementTraits::Quadrature::dim_> integration_points_;

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementCell;

    namespace cell
    {

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        struct FiniteElementCellTraits;

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        requires(t_finite_element.discretization_.isHHO())
        struct FiniteElementCellTraits<t_element, t_domain, t_finite_element>
        {

        private:

            /**
             * @brief
             */
            using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element>;

            /**
             * @brief
             */
            using t_FiniteElementCell = FiniteElementCell<t_element, t_domain, t_finite_element>;

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::integer
            getNumUnknowns()
            {
                return t_FiniteElementTraits::getNumUnknowns();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::integer
            getNumCellUnknowns()
            {
                return t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::integer
            getNumFaceUnknowns()
            {
                return t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Structural()>();
            }

        public:

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::integer
            getLoadVectorSize()
            {
                return t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Subsidiary()>();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::integer
            getUnknownVectorSize()
            {
                auto num_cell_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
                auto num_face_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Structural()>();
                return t_FiniteElementTraits::getNumUnknowns();
            }

            struct Cell
            {

                using Stabilization = lolita::matrix::Matrix<lolita::real, getUnknownVectorSize(), getUnknownVectorSize()>;

                using KTTinv = lolita::matrix::Matrix<lolita::real, getNumCellUnknowns(), getNumCellUnknowns()>;

                using KFT = lolita::matrix::Matrix<lolita::real, getNumFaceUnknowns(), getNumCellUnknowns()>;

                Stabilization stabilization_;

                KTTinv ktt_inv_;

                KFT kft_;

            };

            struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
            {

//                lolita::matrix::Vector<lolita::real, getUnknownVectorSize()>
//                getUnknowns()
//                const
//                {
//                    using t_UnknownVector = lolita::matrix::Vector<lolita::real, getUnknownVectorSize()>;
//                    auto unknowns = lolita::matrix::Zero<t_UnknownVector>();
//                    auto count = 0;
//                    auto set_cell_block = [&] () {
//                        using t_CellUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>::FieldUnknowns ;
//                        using t_CellUnknown = std::tuple_element_t<0, typename t_CellUnknowns::Unknowns>;
//                        for (int i = 0; i < t_FiniteElementTraits::Field::shape_.rows_; ++i) {
//                            for (int j = 0; j < t_FiniteElementTraits::Field::shape_.cols_; ++j) {
//                                auto const & rhs = this->unknowns_.template getUnknown<0>().getDegreeOfFreedom(i, j).getUnknownsCoefficients();
//                                unknowns.template segment<t_CellUnknown::dim_>(count) = rhs;
//                                count += t_CellUnknown::dim_;
//                            }
//                        }
//                    };
//                    auto set_faces_block = [&] <lolita::integer t_i = 0> (auto & self) mutable {
//                        auto const constexpr _face = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<0, t_i>();
//                        using t_FaceUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<_face, t_domain, t_finite_element>::FieldUnknowns;
//                        using t_FaceUnknown = std::tuple_element_t<0, typename t_FaceUnknowns::Unknowns>;
//                        auto const & faces = this->template getComponents<0, t_i>();
//                        for (auto const & face : faces) {
//                            for (int i = 0; i < t_FiniteElementTraits::Field::shape_.rows_; ++i) {
//                                for (int j = 0; j < t_FiniteElementTraits::Field::shape_.cols_; ++j) {
//                                    auto const & rhs = face->unknowns_.template getUnknown<0>().getDegreeOfFreedom(i, j).getUnknownsCoefficients();
//                                    unknowns.template segment<t_FaceUnknown::dim_>(count) = rhs;
//                                    count += t_FaceUnknown::dim_;
//                                }
//                            }
//                        }
//                        if constexpr (t_i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<0>() - 1) {
//                            self.template operator()<t_i + 1>(self);
//                        }
//                    };
//                    set_cell_block();
//                    set_faces_block(set_faces_block);
//                    return unknowns;
//                }

                void
                setGeneralizedGradients()
                {
                    using t_Strain = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows()>;
//                    auto unknowns = this->template getUnknowns<unknown::Unknown::Structural()>();
                    auto structural_unknowns = this->template getUnknowns<unknown::Unknown::Structural()>();
                    auto subsidiary_unknowns = this->template getUnknowns<unknown::Unknown::Subsidiary()>();
                    for (int i = 0; i < t_FiniteElementTraits::Quadrature::dim_; ++i) {
                        /*
                         * update the generalized gradient vector for the current integration point
                         */
//                        auto generalized_gradient_values = t_Strain(this->integration_points_[i].generalized_strain_operator_ * unknowns);
                        auto structural_strain_values = t_Strain(this->integration_points_[i].structural_strain_operator_ * structural_unknowns);
                        auto subsidiary_strain_values = t_Strain(this->integration_points_[i].subsidiary_strain_operator_ * subsidiary_unknowns);
                        auto strain = t_Strain(structural_strain_values + subsidiary_strain_values);
//                        auto generalized_gradient_values = _Gradients::Zero();
//                        generalized_gradient_values.setZero();
                        auto set_generalized_gradient = [&] <lolita::integer t_i = 0> (auto & self) mutable {
                            auto constexpr t_mapping = t_finite_element.unknown_.mappings_[t_i];
                            using t_MappingPolicy = lolita::core2::field::MappingPolicy<t_finite_element.unknown_.tensor_, t_domain, t_mapping>;
                            auto constexpr t_mapping_vector_block = t_FiniteElementTraits::template getMappingBlock<t_mapping>();
                            auto constexpr t_mapping_block_size = t_mapping_vector_block.j_ - t_mapping_vector_block.i_;
                            auto strain_block = strain.template segment<t_mapping_block_size>(t_mapping_vector_block.i_);
                            t_MappingPolicy::non_linear(strain_block);
                            if constexpr (t_i < t_finite_element.unknown_.mappings_.size() - 1) {
                                self.template operator()<t_i + 1>(self);
                            }
                        };
                        set_generalized_gradient(set_generalized_gradient);
                        auto generalized_gradients = lolita::matrix::Span<t_Strain>(this->integration_points_[i].material_point_->s1.gradients.data());
                        generalized_gradients = strain;
                        /*
                         * update the generalized flux, tangent operator, and internal state variables for the current integration point
                         */
                        auto const & behaviour = this->behaviour_->behaviour_;
                        auto behaviour_view = mgis::behaviour::make_view(* this->integration_points_[i].material_point_);
                        auto result = mgis::behaviour::integrate(behaviour_view, behaviour);
                    }
                }

                lolita::matrix::Vector<lolita::real, getUnknownVectorSize()>
                getInternalForces()
                const
                {
                    using t_InternalForcesVector = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getNumUnknowns()>;
                    using t_Stress = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows()>;
                    auto internal_forces = lolita::matrix::Zero<t_InternalForcesVector>();
                    for (int i = 0; i < t_FiniteElementTraits::Quadrature::dim_; ++i) {
                        auto generalized_flux = lolita::matrix::Span<t_Stress>(this->integration_points_[i].material_point_->s1.thermodynamic_forces.data());
                        internal_forces += this->integration_points_[i].weight_ * this->integration_points_[i].operator_.transpose() * generalized_flux;
                    }
                    internal_forces += this->cell_.stabilization_ * this->getUnknowns();
                    return internal_forces;
                }

                lolita::matrix::Vector<lolita::real, getUnknownVectorSize()>
                getExternalForces()
                const
                {
                    using t_ExternalForcesVector = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getNumUnknowns()>;
                    auto external_forces_vector = lolita::matrix::Zero<t_ExternalForcesVector>();
                    auto count = 0;
                    auto TIME___ = 0.0;
                    auto set_cell_block = [&] () {
                        using t_CellFieldUnknownsTraits = unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>;
                        using t_CellUnknowns = typename t_CellFieldUnknownsTraits::FieldUnknowns;
                        using t_CellUnknown = std::tuple_element_t<0, typename t_CellUnknowns::Unknowns>;
                        for (int i = 0; i < t_FiniteElementTraits::Field::shape_.rows_; ++i) {
                            for (int j = 0; j < t_FiniteElementTraits::Field::shape_.cols_; ++j) {
                                for (int k = 0; k < t_FiniteElementTraits::Quadrature::dim_; ++k) {
                                    auto external_cell_forces_vector = external_forces_vector.template segment<t_CellUnknown::dim_>(count);
                                    auto cell_load_value = this->getLoad(i, j)->getImposedValue(this->integration_points_[k].coordinates_, TIME___);
                                    auto const & generalized_load_operator = this->integration_points_[k].generalized_load_operator_;
                                    external_cell_forces_vector += this->integration_points_[k].weight_ * cell_load_value * generalized_load_operator;
                                }
                                count += t_CellUnknown::dim_;
                            }
                        }
                    };
                    auto set_faces_block = [&] <lolita::integer _i = 0> (auto & self) mutable {
                        auto const constexpr t_face = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<0, _i>();
                        using t_FaceFieldUnknownsTraits = unknown::FiniteElementFieldUnknownsTraits<t_face, t_domain, t_finite_element>;
                        using t_FaceUnknowns = typename t_FaceFieldUnknownsTraits::FieldUnknowns;
                        using t_FaceUnknown = std::tuple_element_t<0, typename t_FaceUnknowns::Unknowns>;
                        using t_FaceElementTraits = lolita::core2::finite_element::FiniteElementTraits<t_face, t_domain, t_finite_element>;
                        auto const & faces = this->template getComponents<0, _i>();
                        for (auto const & face : faces) {
                            for (int i = 0; i < t_FaceElementTraits::Field::shape_.rows_; ++i) {
                                for (int j = 0; j < t_FaceElementTraits::Field::shape_.cols_; ++j) {
                                    for (int k = 0; k < t_FaceElementTraits::Quadrature::dim_; ++k) {
                                        auto external_face_forces_vector = external_forces_vector.template segment<t_FaceUnknown::dim_>(count);
                                        auto face_load_value = face->getLoad(i, j)->getImposedValue(face->integration_points_[k].coordinates_, TIME___);
                                        auto const & generalized_load_operator = face->integration_points_[k].generalized_load_operator_;
                                        external_face_forces_vector += face->integration_points_[k].weight_ * face_load_value * generalized_load_operator;
                                    }
                                    count += t_FaceUnknown::dim_;
                                }
                            }
                        }
                        if constexpr (_i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<0>() - 1) {
                            self.template operator()<_i + 1>(self);
                        }
                    };
                    set_cell_block();
                    set_faces_block(set_faces_block);
                    return external_forces_vector;
                }
//
//                lolita::matrix::Matrix<lolita::real, getUnknownVectorSize(), getUnknownVectorSize()>
//                getTangentOperator()
//                const
//                {
//                    using _TangentOperatorMatrix = lolita::matrix::Matrix<lolita::real, t_FiniteElementTraits::getNumUnknowns(), t_FiniteElementTraits::getNumUnknowns()>;
//                    using _TangentOperatorTensor = lolita::matrix::Matrix<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows(), t_FiniteElementTraits::getGeneralizedStrainNumRows()>;
//                    auto tangent_operator_matrix = lolita::matrix::Zero<_TangentOperatorMatrix>();
//                    for (int i = 0; i < t_FiniteElementTraits::Quadrature::dim_; ++i) {
//                        auto count = 0;
//                        auto tangent_operator_tensor = lolita::matrix::Zero<_TangentOperatorTensor>();
//                        auto set_tangent_operator_tensor = [&] <lolita::integer _i = 0> (auto & self) mutable {
//                            auto const constexpr _mapping = t_finite_element.unknown_.mappings_[_i];
//                            auto const constexpr _mapping_vector_block = t_FiniteElementTraits::template getMappingBlock<_mapping>();
//                            auto const constexpr _mapping_block_size = _mapping_vector_block.j_ - _mapping_vector_block.i_;
//                            using _MappingPolicy = lolita::core2::field::MappingPolicy<t_finite_element.unknown_.tensor_, t_domain, _mapping>;
//                            using _MappingMatrix = lolita::matrix::Matrix<lolita::real, _mapping_block_size, _mapping_block_size>;
//                            auto tangent_operator_block = tangent_operator_tensor.template block<_mapping_block_size, _mapping_block_size>(
//                                    _mapping_vector_block.i_,
//                                    _mapping_vector_block.i_
//                            );
//                            auto tangent_operator_values = lolita::matrix::Span<_MappingMatrix const>(this->integration_points_[i].material_point_->K.data() + count);
//                            tangent_operator_block += tangent_operator_values;
//                            count += _mapping_block_size * _mapping_block_size;
//                            if constexpr (_i < t_finite_element.unknown_.mappings_.size() - 1) {
//                                self.template operator()<_i + 1>(self);
//                            }
//                        };
//                        set_tangent_operator_tensor(set_tangent_operator_tensor);
//                        tangent_operator_matrix += this->integration_points_[i].weight_ * this->integration_points_[i].generalized_strain_operator_.transpose() * tangent_operator_tensor * this->integration_points_[i].generalized_strain_operator_;
//                    }
//                    tangent_operator_matrix += this->cell_.stabilization_;
//                    return tangent_operator_matrix;
//                }
//
//                template<auto _arg>
//                void
//                assemble(
//                        lolita::core2::mesh::Mesh<t_domain, _arg> & mesh
//                )
//                {
//                    auto const constexpr _cs = getNumCellUnknowns();
//                    auto const constexpr _fs = getNumFaceUnknowns();
//                    //
//                    using FullMatrix = lolita::matrix::Matrix<lolita::real, getNumUnknowns(), getNumUnknowns()>;
//                    using FullVector = lolita::matrix::Vector<lolita::real, getNumUnknowns()>;
//                    //
//                    using CondMatrix = lolita::matrix::Matrix<lolita::real, getNumFaceUnknowns(), getNumFaceUnknowns()>;
//                    using CondVector = lolita::matrix::Vector<lolita::real, getNumFaceUnknowns()>;
//                    //
//                    using CellCellMatrixBlock = lolita::matrix::Matrix<lolita::real, getNumCellUnknowns(), getNumCellUnknowns()>;
//                    using CellFaceMatrixBlock = lolita::matrix::Matrix<lolita::real, getNumCellUnknowns(), getNumFaceUnknowns()>;
//                    using FaceCellMatrixBlock = lolita::matrix::Matrix<lolita::real, getNumFaceUnknowns(), getNumCellUnknowns()>;
//                    using FaceFaceMatrixBlock = lolita::matrix::Matrix<lolita::real, getNumFaceUnknowns(), getNumFaceUnknowns()>;
//                    //
//                    using CellVectorBlock = lolita::matrix::Vector<lolita::real, getNumCellUnknowns()>;
//                    using FaceVectorBlock = lolita::matrix::Vector<lolita::real, getNumFaceUnknowns()>;
//                    //
//                    auto residual_vector = getInternalForces() - getExternalForces();
//                    auto tangent_matrix = getTangentOperator();
//                    auto Ktt = tangent_matrix.template block<_cs, _cs>(0, 0);
//                    auto Ktf = tangent_matrix.template block<_cs, _fs>(0, _cs);
//                    auto Kft = tangent_matrix.template block<_fs, _cs>(_cs, 0);
//                    auto Kff = tangent_matrix.template block<_fs, _fs>(_cs, _cs);
//                    auto Rt = residual_vector.template segment<_cs>(0);
//                    auto Rf = residual_vector.template segment<_fs>(_cs);
//                    auto Ktt_inv = Ktt.llt().solve(decltype(Ktt)::Identity());
//                    auto Kc = Kff - Kft * Ktt_inv * Ktf;
//                    auto Rc = Rf - Kft * Ktt_inv * Rt;
//                    this->kft_ = Kft;
//                    this->ktt_inv_ = Ktt_inv;
//                    auto constexpr _finite_element_index = _arg.template getFiniteElementIndex<t_finite_element>();
//                    auto & system = mesh.systems_[_finite_element_index];
//                    auto row_offset = 0;
//                    auto col_offset = 0;
//                    auto set_rows = [&] <lolita::integer _i = 0> (auto & t_set_rows) mutable {
//                        auto const & faces_row = this->template getComponents<0, _i>();
//                        auto const constexpr _face_r = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<0, _i>();
//                        using _RowFaceUnknowns = typename lolita::core2::finite_element::unknown::FiniteElementFieldUnknownsTraits<_face_r, t_domain, t_finite_element>::FieldUnknowns;
//                        using _RowFaceUnknown = std::tuple_element_t<0, typename _RowFaceUnknowns::Unknowns>;
////                        using _RowFiniteElementFace = lolita::core2::finite_element::face::FiniteElementFace<_face_r, t_domain, t_finite_element>;
//                        using _RowFiniteElementFace = lolita::core2::finite_element::FiniteElementTraits<_face_r, t_domain, t_finite_element>;
//                        auto set_cols = [&] <lolita::integer _j = 0> (auto & t_set_cols) mutable {
//                            auto const & faces_col = this->template getComponents<0, _j>();
//                            auto const constexpr _face_c = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<0, _j>();
//                            using _ColFaceUnknowns = typename lolita::core2::finite_element::unknown::FiniteElementFieldUnknownsTraits<_face_c, t_domain, t_finite_element>::FieldUnknowns;
//                            using _ColFaceUnknown = std::tuple_element_t<0, typename _ColFaceUnknowns::Unknowns>;
////                            using _ColFiniteElementFace = lolita::core2::finite_element::face::FiniteElementFace<_face_c, t_domain, t_finite_element>;
//                            using _ColFiniteElementFace = lolita::core2::finite_element::FiniteElementTraits<_face_c, t_domain, t_finite_element>;
//                            for (auto const & face_row : faces_row) {
//                                for (auto const & face_col : faces_col) {
//                                    for (int i = 0; i < t_FiniteElementTraits::Field::shape_.rows_; ++i) {
//                                        for (int j = 0; j < t_FiniteElementTraits::Field::shape_.cols_; ++j) {
//                                            auto const & i_rows = face_row->unknwons_.template getUnknown<0>().getUnknownComponent(i, j)->unknowns_coordinates_;
//                                            auto const & i_cols = face_col->unknwons_.template getUnknown<0>().getUnknownComponent(i, j)->unknowns_coordinates_;
//                                            for (auto const & idx_row : i_rows) {
//                                                for (auto const & idx_col : i_cols) {
//                                                    auto triplet = Eigen::Triplet<lolita::real>(row_offset, col_offset, Kc(row_offset, col_offset));
//                                                    system.lhs_values.push_back(Eigen::Triplet<lolita::real>(idx_row, idx_col, Kc(row_offset, col_offset)));
//                                                    system.lhs_values.push_back(lolita::matrix::SparseMatrixInput(idx_row, idx_col, Kc(row_offset, col_offset)));
//                                                    Kc(row_offset, col_offset);
//                                                    col_offset += 1;
//                                                }
//                                                system.rhs_values.push_back(lolita::matrix::SparseVectorInput(idx_row, Rc(row_offset)));
//                                                Rc(row_offset);
//                                                row_offset += 1;
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                            if constexpr (_j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<0>() - 1) {
//                                t_set_cols.template operator()<_j + 1>(t_set_cols);
//                            }
//                        };
//                        if constexpr (_i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<0>() - 1) {
//                            t_set_rows.template operator()<_i + 1>(t_set_rows);
//                        }
//                    };
//                }
//
//                void
//                setElement()
//                {
//                    for (int i = 0; i < t_FiniteElementTraits::Quadrature::dim_; ++i) {
//                        auto const & bhv = this->behaviour_->behaviour_;
//                        this->integration_points_[i].coordinates_ = this->template getCurrentQuadraturePoint<t_FiniteElementTraits::Quadrature::quadrature_, t_FiniteElementTraits::Quadrature::ord_>(i);
//                        std::cout << "my coordinates are : " << i << std::endl;
//                        std::cout << this->integration_points_[i].coordinates_ << std::endl;
//                        this->integration_points_[i].material_point_ = std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(bhv));
//                        this->integration_points_[i].material_point_->K[0] = 4;
//                        auto behaviour_view = mgis::behaviour::make_view(* this->integration_points_[i].material_point_);
//                        auto result = mgis::behaviour::integrate(behaviour_view, bhv);
//                        getUnknowns();
//                        auto res = this->integration_points_[i].operator_ * getUnknowns();
//                        setGeneralizedGradients();
//                        getInternalForces();
//                        getExternalForces();
//                        getTangentOperator();
//                    }
//                }

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
    struct FiniteElementCell
    {

    private:

        /**
         * @brief
         */
        using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         */
        using t_FiniteElementCellTraits = cell::FiniteElementCellTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         */
        using t_Quadrature = typename FiniteElementTraits<t_element, t_domain, t_finite_element>::Quadrature;

    public:

        /**
         * @brief
         */
        lolita::integer static constexpr dim_structural_unknowns = t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Structural()>();

        /**
         * @brief
         */
        lolita::integer static constexpr dim_subsidiary_unknowns = t_FiniteElementTraits::template getDimUnknowns<unknown::Unknown::Subsidiary()>();

        /**
         * @brief
         */
        lolita::integer static constexpr num_structural_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Structural()>();

        /**
         * @brief
         */
        lolita::integer static constexpr num_subsidiary_unknowns = t_FiniteElementTraits::template getNumUnknowns<unknown::Unknown::Subsidiary()>();

        /**
         * @brief
         */
        struct IntegrationPoint
        {

            /**
             * @brief
             */
            lolita::domain::Point coordinates_;

            /**
             * @brief
             */
            lolita::real weight_;

            /**
             * @brief
             */
            std::unique_ptr<mgis::behaviour::BehaviourData> material_point_;

            /**
             * @brief
             */
            lolita::matrix::Matrix<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows(), num_structural_unknowns> structural_strain_operator_;

            /**
             * @brief
             */
            lolita::matrix::Matrix<lolita::real, t_FiniteElementTraits::getGeneralizedStrainNumRows(), num_subsidiary_unknowns> subsidiary_strain_operator_;

            /**
             * @brief
             */
            lolita::matrix::Vector<lolita::real, dim_structural_unknowns> structural_load_operator_;

            /**
             * @brief
             */
            lolita::matrix::Vector<lolita::real, dim_subsidiary_unknowns> subsidiary_load_operator_;

        };

        /**
         * @brief
         */
        struct Implementation : FiniteElement<t_element, t_domain, t_finite_element>
        {

            void
            setIntegrationPoint()
            {
                auto const & behaviour = this->behaviour_->behaviour_;
                for (auto & integration_point : this->integration_points_) {
                    integration_point.material_point_ = std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(behaviour));
                    integration_point.material_point_->K[0] = 4;
                }
            }

        };

        void
        setIntegrationPoint()
        {
            static_cast<Implementation *>(this)->setIntegrationPoint();
        }

        void
        setGeneralizedGradients()
        {
            setIntegrationPoint();
            using t_Implementation = typename t_FiniteElementCellTraits::Implementation;
            static_cast<t_Implementation *>(this)->setGeneralizedGradients();
        }

        /**
         * @brief
         */
        std::array<IntegrationPoint, t_Quadrature::dim_> integration_points_;

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementModule
    {};

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    requires(t_element.isSub(t_domain, 0))
    struct FiniteElementModule<t_element, t_domain, t_finite_element> : FiniteElementCell<t_element, t_domain, t_finite_element>
    {};

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    requires(t_element.isSub(t_domain, 1))
    struct FiniteElementModule<t_element, t_domain, t_finite_element> : FiniteElementFace<t_element, t_domain, t_finite_element>
    {};

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
         * @return
         */
        static constexpr
        lolita::integer
        getDimUnknowns()
        {
            using t_ElementFieldUnknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>;
            auto num_structural_unknowns = t_ElementFieldUnknowns::template getDimUnknowns<unknown::Unknown::Structural()>();
            auto num_subsidiary_unknowns = t_ElementFieldUnknowns::template getDimUnknowns<unknown::Unknown::Subsidiary()>();
            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                auto constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
                using t_NeighbourFieldUnknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>;
                auto num_structural_component_unknowns = t_NeighbourFieldUnknowns::template getDimUnknowns<unknown::Unknown::Structural()>();
                auto num_subsidiary_component_unknowns = t_NeighbourFieldUnknowns::template getDimUnknowns<unknown::Unknown::Subsidiary()>();
                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
                num_structural_unknowns += num_structural_component_unknowns * num_components;
                num_subsidiary_unknowns += num_subsidiary_component_unknowns * num_components;
                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            if constexpr (!t_element.isPoint()) {
                get_num_components_unknowns(get_num_components_unknowns);
            }
            return num_structural_unknowns + num_subsidiary_unknowns;
        }

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown>
        static constexpr
        lolita::integer
        getDimUnknowns()
        {
            auto num_element_unknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                auto const constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
                auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
                num_element_unknowns += num_component_unknowns * num_components;
                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            if constexpr (!t_element.isPoint()) {
                get_num_components_unknowns(get_num_components_unknowns);
            }
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
        getDimUnknowns()
        {
            if constexpr (t_i == 0) {
                return FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
            }
            else {
                auto num_unknowns = 0;
                auto get_num_components_unknowns = [&] <lolita::integer t_j = 0> (auto & self) constexpr mutable {
                    auto const constexpr t_component = ElementDescription::template getComponent<t_i - 1, t_j>();
                    auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
                    auto num_components = ElementDescription::template getNumComponents<t_i - 1, t_j>();
                    num_unknowns += num_component_unknowns * num_components;
                    if constexpr (t_j < ElementDescription::template getNumComponents<t_i - 1>() - 1) {
                        self.template operator()<t_j + 1>(self);
                    }
                };
                if constexpr (!t_element.isPoint()) {
                    get_num_components_unknowns(get_num_components_unknowns);
                }
                return num_unknowns;
            }
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            using t_ElementFieldUnknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>;
            auto num_structural_unknowns = t_ElementFieldUnknowns::template getNumUnknowns<unknown::Unknown::Structural()>();
            auto num_subsidiary_unknowns = t_ElementFieldUnknowns::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                auto constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
                using t_NeighbourFieldUnknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>;
                auto num_structural_component_unknowns = t_NeighbourFieldUnknowns::template getNumUnknowns<unknown::Unknown::Structural()>();
                auto num_subsidiary_component_unknowns = t_NeighbourFieldUnknowns::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
                num_structural_unknowns += num_structural_component_unknowns * num_components;
                num_subsidiary_unknowns += num_subsidiary_component_unknowns * num_components;
                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            if constexpr (!t_element.isPoint()) {
                get_num_components_unknowns(get_num_components_unknowns);
            }
            return num_structural_unknowns + num_subsidiary_unknowns;
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
            if constexpr (!t_element.isPoint()) {
                get_num_components_unknowns(get_num_components_unknowns);
            }
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
                if constexpr (!t_element.isPoint()) {
                    get_num_components_unknowns(get_num_components_unknowns);
                }
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
        getGeneralizedStrainNumRows()
        {
            auto mapping_size = lolita::integer(0);
            auto set_dim_mapping = [&] <lolita::integer t_i = 0> (auto & self) constexpr mutable {
                mapping_size += getMappingSize<t_finite_element.unknown_.mappings_[t_i]>();
                if constexpr (t_i < t_finite_element.unknown_.mappings_.size() - 1) {
                    self.template operator()<t_i + 1>(self);
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
    FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>,
    FiniteElementModule<t_element, t_domain, t_finite_element>
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
            this->template getUnknowns<unknown::Unknown::Subsidiary()>();
            this->template getUnknowns<unknown::Unknown::Structural()>();
            if constexpr (t_element.isSub(t_domain, 0)) {
                this->setGeneralizedGradients();
            }
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