//
// Created by dsiedel on 14/06/22.
//

#ifndef LOLITA_LOLITA_CORE_5_001_BASE_HXX
#define LOLITA_LOLITA_CORE_5_001_BASE_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_4.hxx"
#include "lolita/lolita_core_5_000_connectivity.hxx"

namespace lolita::core::finite_element
{

    /**
     * @brief
     * 
     */
    auto const static null_load_ptr = std::make_shared<lolita::finite_element::ScalarLoad>(lolita::finite_element::ScalarLoad());

    /**
     * @brief
     * 
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementGeometry : virtual FiniteElementConnectivity<t_element, t_domain, t_finite_element>
    {

    private:

        /**
         * @brief
         */
        using t_ElementTraits = lolita::core::geometry::ElementTraits<t_element, t_domain>;

        /**
         * @brief
         */
        using t_FiniteElementTraits = lolita::core::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>;

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
            return _ReferenceCoordinates(lolita::core::geometry::ElementTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
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
                using SegmentQuadrature = lolita::core::geometry::ElementQuadratureTraits<lolita::core::geometry::Element::LinearSegment(), _quadrature, 4>;
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
            return lolita::core::geometry::ElementQuadratureTraits<t_element, _quadrature, _ord>::reference_weights_[index];
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
                    lolita::core::geometry::ElementQuadratureTraits<t_element, _quadrature, _ord>::reference_points_[index].begin()
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
            auto const & elt_reference_nodes = lolita::core::geometry::ElementTraits<t_element, t_domain>::reference_nodes_;
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
     * 
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
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
                lolita::core::mesh::Mesh<t_domain, t_element_group> & mesh
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
     * 
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementFieldLoad : virtual FiniteElementConnectivity<t_element, t_domain, t_finite_element>
    {

    private:

        /**
         * @brief
         */
        using t_Field = lolita::core::field::TensorPolicy<t_finite_element.unknown_.tensor_, t_domain.dim_>;

        /**
         * @brief loads initializer
         * @return loads member initialized with the zero natural load
         */
        static
        std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_>
        setLoads()
        {
            auto loads = std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_>();
            for (int i = 0; i < t_Field::shape_.rows_; ++i) {
                for (int j = 0; j < t_Field::shape_.cols_; ++j) {
                    loads[j][i] = lolita::core::finite_element::null_load_ptr;
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
                lolita::core::mesh::Mesh<t_domain, t_element_group> & mesh
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
        std::array<std::array<std::shared_ptr<lolita::finite_element::ScalarLoad>, t_Field::shape_.rows_>, t_Field::shape_.cols_> field_load_ = setLoads();

    };

    /**
     * @brief 
     * 
     * @tparam t_element 
     * @tparam t_domain 
     * @tparam t_finite_element 
     */
    template<lolita::core::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementUnknowns : virtual FiniteElementConnectivity<t_element, t_domain, t_finite_element>
    {

        /**
         * @brief
         * 
         * @return
         */
        lolita::boolean
        isBound()
        const
        {
            return binding_index_ == -1;
        }

        /**
         * @brief 
         * 
         */
        lolita::natural unknown_index_ = -1;
        
        /**
         * @brief 
         * 
         */
        lolita::natural binding_index_ = -1;

        /**
         * @brief 
         * 
         */
        std::shared_ptr<lolita::core::system::FiniteElementLinearSystem> system_;

    };

    /**
     * @brief
     * 
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core::geometry::Element t_element, lolita::domain::Domain t_domain, auto t_finite_element>
    struct FiniteElementBase
    :
    FiniteElementGeometry<t_element, t_domain, t_finite_element>,
    FiniteElementBehaviour<t_element, t_domain, t_finite_element>,
    FiniteElementFieldLoad<t_element, t_domain, t_finite_element>,
    FiniteElementUnknowns<t_element, t_domain, t_finite_element>
    {};

}

#endif //LOLITA_LOLITA_CORE_5_001_BASE_HXX
