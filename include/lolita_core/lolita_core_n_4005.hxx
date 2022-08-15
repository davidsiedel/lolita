#ifndef FC0FBEA8_1C07_47E6_97F1_3D709687C47A
#define FC0FBEA8_1C07_47E6_97F1_3D709687C47A

#include "lolita_core/lolita.hxx"
#include "lolita_core/lolita_core_n_0000.hxx"
#include "lolita_core/lolita_core_n_1000.hxx"
#include "lolita_core/lolita_core_n_2000.hxx"
#include "lolita_core/lolita_core_n_3000.hxx"
#include "lolita_core/lolita_core_n_4001.hxx"
#include "lolita_core/lolita_core_n_4002.hxx"
#include "lolita_core/lolita_core_n_4003.hxx"
#include "lolita_core/lolita_core_n_4004.hxx"

namespace lolita
{

    template<Element t_element, Domain t_domain>
    struct FiniteElementHolder
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;
        
        template<Element t__element, Domain t__domain>
        using t_ElementPointer = std::shared_ptr<FiniteElementHolder<t__element, t__domain>>;

    public:
    
        using t_InnerNeighbors = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer>;
        
        using t_OuterNeighbors = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer>;

        template<Basis t_basis>
        using t_Basis = typename FiniteElementBasisTraits<t_basis>::template Implementation<t_element, t_domain>;

        template<auto t_discretization>
        using t_Disc = typename HybridDiscontinuousGalerkinTraits<t_discretization>::template Implementation<t_element, t_domain>;
        
        Boolean
        operator==(
            FiniteElementHolder const & other
        )
        const = default;
        
        Boolean
        operator!=(
            FiniteElementHolder const & other
        )
        const = default;

        Boolean
        isIn(
            std::basic_string_view<Character> domain
        )
        const
        {
            auto has_domain = [&] (
                std::shared_ptr<MeshDomain> const & mesh_domain
            )
            {
                return mesh_domain->tag_ == domain;
            };
            return std::find_if(domains_.begin(), domains_.end(), has_domain) != domains_.end();
        }
        
        std::basic_string<Character>
        getHash()
        const
        requires(t_element.isNode())
        {
            return std::to_string(this->tag_);
        }
        
        std::basic_string<Character>
        getHash()
        const
        {
            std::basic_stringstream<Character> hash;
            auto const & nodes = getInnerNeighbors<t_element.dim_ - 1, 0>();
            for (auto const & node : nodes)
            {
                hash << node->getHash();
            }
            return hash.str();
        }
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer i,
            Integer j
        )
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_OuterNeighbors>> &
        getOuterNeighbors()
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_OuterNeighbors>> const &
        getOuterNeighbors()
        const
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_InnerNeighbors>> &
        getInnerNeighbors()
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_InnerNeighbors>> const &
        getInnerNeighbors()
        const
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborIndex(
            Integer i
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_inner_neighbor = t_ElementTraits::template getInnerNeighbor<t_i, t_j>();
            auto constexpr t_coordinates = ElementTraits<t_inner_neighbor, t_domain>::template getOuterNeighborCoordinates<t_element>();
            auto const & items = getInnerNeighbors<t_i, t_j>()[i]->template getOuterNeighbors<t_coordinates.dim_, t_coordinates.tag_>();
            auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element)
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborIndex(
            std::shared_ptr<FiniteElementHolder<ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>();
            auto const & inner_neighbors = getInnerNeighbors<t_i, t_j>();
            auto is_equal = [&] (std::shared_ptr<FiniteElementHolder<t_inner_neighbor, t_domain>> const & neighbor)
            {
                return * neighbor == * ptr_neighbor;
            };
            auto neighbor_index = std::distance(inner_neighbors.begin(), std::find_if(inner_neighbors.begin(), inner_neighbors.end(), is_equal));
            return getInnerNeighborIndex<t_i, t_j>(neighbor_index);
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborOrientation(
            Integer i
        )
        const
        requires(!t_element.isNode())
        {
            return getInnerNeighborIndex<t_i, t_j>(i) == 0 ? 1 : -1;
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborOrientation(
            std::shared_ptr<FiniteElementHolder<ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        )
        const
        requires(!t_element.isNode())
        {
            return getInnerNeighborIndex<t_i, t_j>(ptr_neighbor) == 0 ? 1 : -1;
        }
    
        Point &
        getCurrentCoordinates()
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        Point const &
        getCurrentCoordinates()
        const
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        lolita::algebra::Matrix<Real, 3, t_element.getNumNodes()>
        getCurrentCoordinates()
        const
        requires(!t_element.isNode())
        {
            auto current_nodes_coordinates = lolita::algebra::Matrix<Real, 3, t_element.getNumNodes()>();
            auto count = Integer(0);
            for (auto const & node : this->template getInnerNeighbors<t_element.dim_ - 1, 0>())
            {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }
        
        static
        lolita::algebra::Span<lolita::algebra::Matrix<Real, 3, t_element.getNumNodes(), lolita::algebra::colMajor()> const>
        getReferenceCoordinates()
        {
            using t_ReferenceCoordinates = lolita::algebra::Span<lolita::algebra::Matrix<Real, 3, t_element.getNumNodes(), lolita::algebra::colMajor()> const>;
            return t_ReferenceCoordinates(ElementTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
        }
        
        static
        Real
        getShapeMappingEvaluation(
            lolita::algebra::Vector<Real, t_element.getNumNodes()> const & nodal_field_values,
            Point const & reference_point
        )
        {
            return t_ElementTraits::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }
        
        static
        Real
        getShapeMappingDerivative(
            lolita::algebra::Vector<Real, t_element.getNumNodes()> const & nodal_field_values,
            Point const & reference_point,
            Integer derivative_direction
        )
        {
            return t_ElementTraits::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::algebra::Matrix<Real, 3, 3>();
            auto du = Real(0);
            ru.setZero();
            for (auto i = 0; i < t_domain.dim_; ++i)
            {
                for (auto j = 0; j < t_element.dim_; ++j)
                {
                    ru(i, j) = FiniteElementHolder::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.dim_ == 3)
            {
                du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            else if constexpr (t_element.dim_ == 2)
            {
                du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else
            {
                du = lolita::numerics::abs(ru.col(0).norm());
            }
            if constexpr (t_domain.frame_ == Domain::Frame::AxiSymmetric)
            {
                Real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10)
                {
                    r0 = 1.e-10;
                }
                du *= 2.0 * lolita::numerics::pi() * r0;
            }
            return du;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction = -1
        )
        const
        {
            if constexpr (t_ElementTraits::hasDim(0))
            {
                auto const & current_coordinates = this->getCurrentCoordinates();
                auto distance = Real();
                auto mp0 = Point();
                auto mp1 = Point();
                for (auto i = 0; i < t_element.dim_; ++i)
                {
                    mp0(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                    mp1(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
                }
                direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
                return distance;
            }
            else
            {
                using SegmentQuadrature = ElementQuadratureRuleTraits<Element::segment(1), Quadrature::gauss(4)>;
                auto distance = Real(0);
                auto dt = Real();
                auto const current_nodes_coordinates = this->getCurrentCoordinates();
                for (auto q = 0; q < SegmentQuadrature::dim_; ++q)
                {
                    auto pq = SegmentQuadrature::reference_points_[q][0];
                    auto wq = SegmentQuadrature::reference_weights_[q];
                    auto ru = lolita::algebra::Matrix<Real, 3, 3>().setZero();
                    auto difference = second_point - first_point;
                    auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                    for (auto i = 0; i < t_domain.dim_; ++i)
                    {
                        for (auto j = 0; j < t_element.dim_; ++j)
                        {
                            if (direction == -1 || i == static_cast<Integer>(direction))
                            {
                                auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                                auto dx = FiniteElementHolder::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                                ru(i, j) = dx * du;
                            }
                        }
                    }
                    if constexpr (t_element.isFacet())
                    {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        auto Fff = ru.col(0).template dot(ru.col(1));
                        auto Gff = ru.col(1).template dot(ru.col(1));
                        dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                    }
                    else if constexpr (t_element.isCurve())
                    {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        dt = std::sqrt(Eff);
                    }
                    else
                    {
                        dt = 0;
                    }
                    distance += wq * dt;
                }
                return distance;
            }
        }
        
        Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElementHolder::getReferenceCoordinates();
            auto current_diameters = Point().setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i)
            {
                for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
                {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (auto k = 0; k < 3; ++k)
                    {
                        auto new_value = lolita::numerics::abs(getRiemannianDistance(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }
        
        Point
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i) {
                barycenter += current_nodes_coordinates.col(i);
            }
            barycenter /= Real(t_element.getNumNodes());
            return barycenter;
        }
        
        static
        Point
        getReferenceCentroid()
        {
            auto reference_nodes_coordinates = FiniteElementHolder::getReferenceCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i) {
                barycenter += reference_nodes_coordinates.col(i);
            }
            barycenter /= Real(t_element.getNumNodes());
            return barycenter;
        }
        
        static
        Point
        getReferenceDiameters()
        {
            auto dts = Point();
            auto nds = FiniteElementHolder::getReferenceCoordinates();
            dts.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i)
            {
                for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
                {
                    for (auto k = 0; k < 3; ++k)
                    {
                        auto new_value = lolita::numerics::abs(nds(k, i) - nds(k, j));
                        auto & current_value = dts(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return dts;
        }
        
        Point
        getNormalVector(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1))
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::algebra::Matrix<Real, 3, t_element.dim_>();
            ru.setZero();
            for (auto i = 0; i < 3; ++i)
            {
                for (auto j = 0; j < t_element.dim_; ++j)
                {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.isNode()) {
                return Point{0, 0, 0};
            }
            else if constexpr (t_element.isCurve())
            {
                return Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
            }
            else
            {
                return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            }
        }
        
        Point
        getTangentVector(
            Point const & point,
            Integer direction
        )
        const
        requires(t_element.isCurve())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = Point();
            for (auto i = 0; i < 3; ++i)
            {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
            }
            return tangent_vector;
        }
        
        template<Quadrature t_quadrature>
        static
        Real
        getReferenceQuadratureWeight(
            Integer index
        )
        {
            return ElementQuadratureRuleTraits<t_element, t_quadrature>::reference_weights_[index];
        }
        
        template<Quadrature t_quadrature>
        static
        lolita::algebra::Span<Point const>
        getReferenceQuadraturePoint(
            Integer index
        )
        {
            return lolita::algebra::Span<Point const>(ElementQuadratureRuleTraits<t_element, t_quadrature>::reference_points_[index].begin());
        }
        
        template<Quadrature t_quadrature>
        Real
        getCurrentQuadratureWeight(
            Integer index
        )
        const
        {
            auto w = getReferenceQuadratureWeight<t_quadrature>(index);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<t_quadrature>(index));
        }
        
        template<Quadrature t_quadrature>
        Point
        getCurrentQuadraturePoint(
            Integer index
        )
        const
        {
            auto p = Point();
            auto const nds = this->getCurrentCoordinates();
            for (auto j = 0; j < 3; ++j)
            {
                p(j) = FiniteElementHolder::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<t_quadrature>(index));
            }
            return p;
        }
        
        template<Quadrature t_quadrature, Integer _i, Integer _j>
        static
        Real
        getInnerNeighborReferenceQuadratureWeight(
            Integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits::template getInnerNeighbor<_i, _j>();
            using ComponentGeometry = FiniteElementHolder<_component, t_domain>;
            return ComponentGeometry::template getReferenceQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, Integer _i, Integer _j>
        static
        Point
        getInnerNeighborReferenceQuadraturePoint(
            Integer component_index,
            Integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits ::template getInnerNeighbor<_i, _j>();
            auto p = Point();
            using ComponentGeometry = FiniteElementHolder<_component, t_domain>;
            auto const & elt_reference_nodes = ElementTraits<t_element, t_domain>::reference_nodes_;
            for (auto i = 0; i < 3; ++i)
            {
                auto cpt_coordinates = lolita::algebra::Vector<Real, _component.getNumNodes()>();
                for (auto j = 0; j < _component.getNumNodes(); ++j)
                {
                    auto const node_tag = getInnerNeighborNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<t_quadrature>(index);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }
        
        template<Quadrature t_quadrature, Integer t_i, Integer t_j>
        Real
        getInnerNeighborCurrentQuadratureWeight(
            Integer component_index,
            Integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto const & cmp =  this->template getInnerNeighbors<t_i, t_j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, Integer t_i, Integer t_j>
        Point
        getInnerNeighborCurrentQuadraturePoint(
            Integer component_index,
            Integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto p = Point();
            auto const cpt_ref_pnt = getInnerNeighborReferenceQuadraturePoint<t_quadrature, t_i, t_j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            for (auto j = 0; j < 3; ++j)
            {
                p(j) = FiniteElementHolder::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // NEW
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        Boolean
        hasFiniteElement(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::shared_ptr<FiniteElement<t_element, t_domain>> const & finite_element
            )
            {
                return finite_element->getLabel() == label;
            };
            return std::find_if(finite_elements_.begin(), finite_elements_.end(), has_label) != finite_elements_.end();
        }

        Boolean
        hasBehavior(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::shared_ptr<ElementIntegrationPoints<t_domain>> const & element_integration_points
            )
            {
                return element_integration_points->getLabel() == label;
            };
            return std::find_if(integration_points_.begin(), integration_points_.end(), has_label) != integration_points_.end();
        }

        Integer
        getFiniteElementIndex(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::shared_ptr<FiniteElement<t_element, t_domain>> const & finite_element
            )
            {
                return finite_element->getLabel() == label;
            };
            if (std::find_if(finite_elements_.begin(), finite_elements_.end(), has_label) != finite_elements_.end())
            {
                return std::distance(finite_elements_.begin(), std::find_if(finite_elements_.begin(), finite_elements_.end(), has_label));
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }

        Integer
        getBehaviorIndex(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::shared_ptr<ElementIntegrationPoints<t_domain>> const & element_integration_points
            )
            {
                return element_integration_points->getLabel() == label;
            };
            if (std::find_if(integration_points_.begin(), integration_points_.end(), has_label) != integration_points_.end())
            {
                return std::distance(integration_points_.begin(), std::find_if(integration_points_.begin(), integration_points_.end(), has_label));
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }
        
        std::shared_ptr<FiniteElement<t_element, t_domain>> const &
        getFiniteElement(
            std::basic_string_view<Character> label
        )
        const
        {
            return finite_elements_[getFiniteElementIndex(label)];
        }
        
        std::shared_ptr<FiniteElement<t_element, t_domain>> &
        getFiniteElement(
            std::basic_string_view<Character> label
        )
        {
            return finite_elements_[getFiniteElementIndex(label)];
        }
        
        std::shared_ptr<ElementIntegrationPoints<t_domain>> const &
        getIntegrationPoints(
            std::basic_string_view<Character> label
        )
        const
        {
            return integration_points_[getBehaviorIndex(label)];
        }
        
        std::shared_ptr<ElementIntegrationPoints<t_domain>> &
        getIntegrationPoints(
            std::basic_string_view<Character> label
        )
        {
            return integration_points_[getBehaviorIndex(label)];
        }

        void
        addElement(
            std::shared_ptr<std::basic_string<Character>> const & label
        )
        {
            if (!hasFiniteElement(* label))
            {
                finite_elements_.push_back(std::make_shared<FiniteElement<t_element, t_domain>>(label));
            }
        }

        template<Quadrature t_quadrature>
        void
        addBehavior(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            if (!hasBehavior(behavior->behaviour))
            {
                auto element_integration_points = std::make_shared<ElementIntegrationPoints<t_domain>>(behavior);
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
                {
                    auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto integration_point = IntegrationPoint<t_domain>(point, weight, behavior);
                    element_integration_points->integration_points_.push_back(std::move(integration_point));
                }
                std::cout << "dt : " << element_integration_points->integration_points_[0].behavior_data_->dt << std::endl;
                integration_points_.push_back(element_integration_points);
            }
        }

        void
        addBehavior(
            std::basic_string_view<Character> finite_element_label,
            std::basic_string_view<Character> behavior_label
        )
        {
            getFiniteElement(finite_element_label)->element_integration_points_ = getIntegrationPoints(behavior_label);
        }
        
        template<Basis t_basis>
        RealVector<t_Basis<t_basis>::getSize()>
        getBasisEvaluation(
            Point const & point
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisEvaluation(point);
        }
        
        template<Basis t_basis>
        RealVector<t_Basis<t_basis>::getSize()>
        getBasisDerivative(
            Point const & point,
            Integer derivative_direction
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisDerivative(point, derivative_direction);
        }

        template<Field t_field, Mapping t_mapping, auto t_discretization>
        RealMatrix<MappingTraits<t_mapping>::template getSize<t_domain, t_field>(), t_Disc<t_discretization>::template getNumElementUnknowns<t_field>()>
        getMapping(
            Point const & point
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getMapping<t_mapping, t_field>(point);
        }

        template<Field t_field, Quadrature t_quadrature, auto t_discretization>
        void
        makeQuadrature(
            std::basic_string_view<Character> label
        )
        {
            for (auto i = 0; i < getFiniteElement(label)->element_integration_points_->getSize(); i++)
            {
                auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
                auto mapp = this->template getMapping<t_field, Mapping::gradient(), t_discretization>(point);
                // element_integration_points_->integration_points_[i].quadrature_operators_[0] = mapp;
            }
            // for (auto const & integration_point : element_integration_points_->integration_points_)
            // {
            //     integration_point.quadrature_operators_[0] = 
            // }
            // std::cout << "mapp : " << std::endl;
            // std::cout << mapp << std::endl;
        }
        
        t_OuterNeighbors outer_neighbors_;
        
        t_InnerNeighbors inner_neighbors_;
        
        Natural tag_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;
        
        std::shared_ptr<Point> coordinates_;

        std::vector<std::shared_ptr<FiniteElement<t_element, t_domain>>> finite_elements_;

        std::vector<std::shared_ptr<ElementIntegrationPoints<t_domain>>> integration_points_;

        std::vector<std::shared_ptr<ElementIntegrationPoints<t_domain>>> integration_points2_;

    };

} // namespace lolita


#endif /* FC0FBEA8_1C07_47E6_97F1_3D709687C47A */
