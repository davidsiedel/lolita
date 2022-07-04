#ifndef C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC
#define C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC

#include <execution>
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_core_n_1.hxx"

namespace lolita2::geometry
{
    
    template<Domain _domain>
    struct Mesh;

    template<Element t_element, Domain t_domain>
    struct FiniteElement : FiniteElementConnectivity<t_element, t_domain>
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;

    public:
    
        lolita2::Point &
        getCurrentCoordinates()
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        lolita2::Point const &
        getCurrentCoordinates()
        const
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>
        getCurrentCoordinates()
        const
//        requires(!t_element.isNode())
        {
            auto current_nodes_coordinates = lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>();
            auto count = lolita::index(0);
            for (auto const & node : this->template getComponents<t_element.dim_ - 1, 0>()) {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }
        
        static
        lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_, lolita::matrix::col_major> const>
        getReferenceCoordinates()
        {
            using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_, lolita::matrix::col_major> const>;
            return _ReferenceCoordinates(ElementTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
        }
        
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const & nodal_field_values,
                lolita2::Point const & reference_point
        )
        {
            return t_ElementTraits::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }
        
        static
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const & nodal_field_values,
                lolita2::Point const & reference_point,
                lolita::index derivative_direction
        )
        {
            return t_ElementTraits::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }
        
        lolita::real
        getShapeMappingDifferential(
                lolita2::Point const & point
        )
        const
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>();
            auto du = lolita::real(0);
            ru.setZero();
            for (lolita::index i = 0; i < t_domain.dim_; ++i) {
                for (lolita::index j = 0; j < t_element.dim_; ++j) {
                    ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.dim_ == 3) {
                du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            else if constexpr (t_element.dim_ == 2) {
                du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else {
                du = lolita::numerics::abs(ru.col(0).norm());
            }
            if constexpr (t_domain.frame_ == lolita2::Domain::Frame::AxiSymmetric) {
                lolita::real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10) {
                    r0 = 1.e-10;
                }
                du *= 2.0 * lolita::numerics::pi * r0;
            }
            return du;
        }
        
        lolita::real
        getRiemannianDistance(
                lolita2::Point const & first_point,
                lolita2::Point const & second_point,
                lolita::integer direction = -1
        )
        const
        {
            if constexpr (t_ElementTraits::hasDim(0)) {
                auto const & current_coordinates = this->getCurrentCoordinates();
                auto distance = lolita::real();
                auto mp0 = lolita2::Point();
                auto mp1 = lolita2::Point();
                for (lolita::index i = 0; i < t_element.dim_; ++i) {
                    mp0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                    mp1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
                }
                direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
                return distance;
            }
            else {
                // auto constexpr _quadrature = lolita2::Quadrature::gauss(4);
                using SegmentQuadrature = ElementQuadratureTraits<Element::segment(1), lolita2::Quadrature::gauss(4)>;
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
                                auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
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
        
        lolita2::Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElement::getReferenceCoordinates();
            auto current_diameters = lolita2::Point().setZero();
            for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (lolita::index k = 0; k < 3; ++k) {
                        auto new_value = lolita::numerics::abs(getRiemannianDistance(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }
        
        lolita2::Point
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            return lolita::domain::getBarycenter(current_nodes_coordinates);
        }
        
        static
        lolita2::Point
        getReferenceCentroid()
        {
            auto nds = FiniteElement::getReferenceCoordinates();
            return lolita::domain::getBarycenter(nds);
        }
        
        static
        lolita2::Point
        getReferenceDiameters()
        {
            auto dts = lolita2::Point();
            auto nds = FiniteElement::getReferenceCoordinates();
            dts.setZero();
            for (lolita::index i = 0; i < t_element.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < t_element.num_nodes_; ++j) {
                    for (lolita::index k = 0; k < 3; ++k) {
                        auto new_value = lolita::numerics::abs(nds(k, i) - nds(k, j));
                        auto & current_value = dts(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return dts;
        }
        
        template<lolita2::Quadrature _quadrature>
        static
        lolita::real
        getReferenceQuadratureWeight(
                lolita::index index
        )
        {
            return ElementQuadratureTraits<t_element, _quadrature>::reference_weights_[index];
        }
        
        template<lolita2::Quadrature _quadrature>
        static
        lolita::matrix::Span<lolita2::Point const>
        getReferenceQuadraturePoint(
                lolita::index index
        )
        {
            return lolita::matrix::Span<lolita2::Point const>(
                    ElementQuadratureTraits<t_element, _quadrature>::reference_points_[index].begin()
            );
        }
        
        template<lolita2::Quadrature _quadrature>
        lolita::real
        getCurrentQuadratureWeight(
                lolita::index index
        )
        const
        {
            auto w = getReferenceQuadratureWeight<_quadrature>(index);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<_quadrature>(index));
        }
        
        template<lolita2::Quadrature _quadrature>
        lolita2::Point
        getCurrentQuadraturePoint(
                lolita::index index
        )
        const
        {
            auto p = lolita2::Point();
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = FiniteElement::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<_quadrature>(index));
            }
            return p;
        }
        
        template<lolita2::Quadrature _quadrature, lolita::index _i, lolita::index _j>
        static
        lolita::real
        getComponentReferenceQuadratureWeight(
                lolita::index index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
            using ComponentGeometry = FiniteElement<_component, t_domain>;
            return ComponentGeometry::template getReferenceQuadratureWeight<_quadrature>(index);
        }
        
        template<lolita2::Quadrature _quadrature, lolita::index _i, lolita::index _j>
        static
        lolita2::Point
        getComponentReferenceQuadraturePoint(
                lolita::index component_index,
                lolita::index index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
            auto p = lolita2::Point();
            using ComponentGeometry = FiniteElement<_component, t_domain>;
            auto const & elt_reference_nodes = ElementTraits<t_element, t_domain>::reference_nodes_;
            for (lolita::index i = 0; i < 3; ++i) {
                auto cpt_coordinates = lolita::matrix::Vector<lolita::real, _component.num_nodes_>();
                for (lolita::index j = 0; j < _component.num_nodes_; ++j) {
                    auto const node_tag = getComponentNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<_quadrature>(index);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }
        
        template<lolita2::Quadrature _quadrature, lolita::index _i, lolita::index _j>
        lolita::real
        getComponentCurrentQuadratureWeight(
                lolita::index component_index,
                lolita::index index
        )
        const
        requires(!t_element.isNode())
        {
            auto const & cmp =  this->template getComponents<_i, _j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<_quadrature>(index);
        }
        
        template<lolita2::Quadrature _quadrature, lolita::index _i, lolita::index _j>
        lolita2::Point
        getComponentCurrentQuadraturePoint(
                lolita::index component_index,
                lolita::index index
        )
        const
        requires(!t_element.isNode())
        {
            auto p = lolita2::Point();
            auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<_quadrature, _i, _j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = FiniteElement::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }
        
        lolita2::Point
        getNormalVector(
                lolita2::Point const & point
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
            if constexpr (t_element.isNode()) {
                return lolita2::Point{0, 0, 0};
            }
            else if constexpr (t_element.isCurve()) {
                return lolita2::Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
            }
            else {
                return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            }
        }
        
        lolita2::Point
        getTangentVector(
                lolita2::Point const & point,
                lolita::index direction
        )
        const
        requires(t_element.isCurve())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = lolita2::Point();
            for (lolita::index i = 0; i < 3; ++i) {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
            }
            return tangent_vector;
        }

    };

}


#endif /* C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC */
