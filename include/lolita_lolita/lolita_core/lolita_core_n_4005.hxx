#ifndef FC0FBEA8_1C07_47E6_97F1_3D709687C47A
#define FC0FBEA8_1C07_47E6_97F1_3D709687C47A

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_3000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4001.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4002.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4003.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4004.hxx"

namespace lolita
{

    template<Domain t_domain>
    struct IntegrationPoint2
    {

    private:

        template<FiniteElementMethodConcept auto t_finite_element_method>
        static constexpr
        Integer
        getGeneralizedStrainSize()
        {
            return FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
        }

        template<BehaviorConcept auto t_behavior>
        static constexpr
        Integer
        getGeneralizedStrainSize()
        {
            return BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
        }

    public:

        IntegrationPoint2(
            algebra::View<Point const> ref_pt,
            Point const & coordinates,
            Real const & weight,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        :
        reference_coordinates_(ref_pt),
        coordinates_(coordinates),
        weight_(weight),
        behavior_(behavior),
        behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behavior)))
        {
            behavior_data_->K[0] = 4;
        }
    
        inline
        Boolean
        operator==(
            IntegrationPoint2 const & other
        )
        const = default;
        
        inline
        Boolean
        operator!=(
            IntegrationPoint2 const & other
        )
        const = default;

        void
        integrate()
        {
            auto behaviour_data_view = mgis::behaviour::make_view(* behavior_data_);
            auto res = mgis::behaviour::integrate(behaviour_data_view, * behavior_);
            // auto strain_view = lolita::algebra::View<Vector<Real> const>(behavior_data_->s1.gradients.data(), behavior_data_->s1.gradients.size());
            // auto stress_view = lolita::algebra::View<Vector<Real> const>(behavior_data_->s1.thermodynamic_forces.data(), behavior_data_->s1.thermodynamic_forces.size());
            // auto K = lolita::algebra::View<Vector<Real> const>(behavior_data_->K.data(), behavior_data_->K.size());
            // std::cout << "strain : " << strain_view << std::endl;
            // std::cout << "stress : " << stress_view << std::endl;
            // std::cout << "K : " << K << std::endl;
            // std::cout << "res : " << res << std::endl;
        }

        void
        setMaterialProperty(
            std::basic_string_view<Character> material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto value = std::forward<std::function<Real(Point const &)>>(function)(coordinates_);
            mgis::behaviour::setMaterialProperty(behavior_data_->s0, std::string(material_property_label), value);
            mgis::behaviour::setMaterialProperty(behavior_data_->s1, std::string(material_property_label), value);
        }

        void
        setExternalVariable(
            std::basic_string_view<Character> material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto value = std::forward<std::function<Real(Point const &)>>(function)(coordinates_);
            mgis::behaviour::setExternalStateVariable(behavior_data_->s0, std::string(material_property_label), value);
            mgis::behaviour::setExternalStateVariable(behavior_data_->s1, std::string(material_property_label), value);
        }
        
        template<BehaviorConcept auto t_behavior>
        lolita::algebra::Span<Vector<Real, BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
            return lolita::algebra::Span<Vector<Real, size> const>(behavior_data_->s1.gradients.data());
        }
        
        template<auto t_finite_element_method>
        lolita::algebra::Span<Vector<Real, FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
            return lolita::algebra::Span<Vector<Real, size> const>(behavior_data_->s1.gradients.data() + offset);
        }
        
        template<auto t_finite_element_method>
        lolita::algebra::Span<Vector<Real, FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()>>
        getGeneralizedStrain()
        {
            auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
            return lolita::algebra::Span<Vector<Real, size>>(behavior_data_->s1.gradients.data() + offset);
        }

        template<FiniteElementMethodConcept auto t_finite_element_method>
        Matrix<Real, getGeneralizedStrainSize<t_finite_element_method>(), getGeneralizedStrainSize<t_finite_element_method>()>
        getJacobian()
        const
        {
            auto constexpr strain_operator_num_rows = getGeneralizedStrainSize<t_finite_element_method>();
            using IntegrationPointJacobian = Matrix<Real, strain_operator_num_rows, strain_operator_num_rows>;
            auto integration_point_jacobian = IntegrationPointJacobian();
            integration_point_jacobian.setZero();
            auto set_mapping_block = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr mapping = t_finite_element_method.template getMapping<t_i>();
                auto constexpr mapping_size = FiniteElementMethodTraits<t_finite_element_method>::template getMappingSize<t_domain, mapping>();
                auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getMappingOffset<t_domain, mapping>();
                auto jacobian_block_view = algebra::View<Matrix<Real, mapping_size, mapping_size> const>(behavior_data_->K.data() + offset * offset);
                auto integration_point_jacobian_block_view = integration_point_jacobian.template block<mapping_size, mapping_size>(offset, offset);
                integration_point_jacobian_block_view = jacobian_block_view;
                if constexpr (t_i < t_finite_element_method.getGeneralizedStrain().getNumMappings() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_mapping_block(set_mapping_block);
            return integration_point_jacobian;
        }

        Point &
        getCurrentCoordinates()
        {
            return coordinates_;
        }
        
        Point const &
        getCurrentCoordinates()
        const
        {
            return coordinates_;
        }
        
        algebra::View<Point const>
        getReferenceCoordinates()
        const
        {
            return reference_coordinates_;
        }

        algebra::View<Point const> reference_coordinates_;
        
        Real reference_weight_;
    
        Point coordinates_;
        
        Real weight_;
    
        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

        std::vector<Matrix<Real>> operators_;

        std::map<std::basic_string<Character>, Matrix<Real>> ops_;

    };

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

        template<Basis t_basis>
        static constexpr
        Integer
        getBasisSize()
        {
            return t_Basis<t_basis>::getSize();
        }

        template<Field t_field>
        static constexpr
        Integer
        getFieldSize()
        {
            return FieldTraits<t_field>::template getSize<t_domain>();
        }

        template<Field t_field, Basis t_basis>
        static constexpr
        Integer
        getFieldSize()
        {
            return t_Basis<t_basis>::getSize() * FieldTraits<t_field>::template getSize<t_domain>();
        }
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer component_index,
            Integer node_index
        )
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[component_index][node_index];
        }
        
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
        {
            auto hash = std::basic_stringstream<Character>();
            for (auto const & node : getInnerNeighbors<t_element.getDim() - 1, 0>())
            {
                hash << node->getHash();
            }
            return hash.str();
        }
        
        std::basic_string<Character>
        getHash()
        const
        requires(t_element.isNode())
        {
            return std::to_string(this->tag_);
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
            // auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element)
            auto is_equal = [&] (std::shared_ptr<FiniteElementHolder<t_element, t_domain>> const & ptr_element)
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }
        
        // template<Integer t_i, Integer t_j>
        // Integer
        // getInnerNeighborIndex(
        //     std::shared_ptr<FiniteElementHolder<ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        // )
        // const
        // requires(!t_element.isNode())
        // {
        //     auto const & inner_neighbors = getInnerNeighbors<t_i, t_j>();
        //     // auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>();
        //     // auto is_equal = [&] (std::shared_ptr<FiniteElementHolder<t_inner_neighbor, t_domain>> const & neighbor)
        //     // {
        //     //     return * neighbor == * ptr_neighbor;
        //     // };
        //     // auto neighbor_index = std::distance(inner_neighbors.begin(), std::find_if(inner_neighbors.begin(), inner_neighbors.end(), is_equal));
        //     auto neighbor_index = std::distance(inner_neighbors.begin(), std::find(inner_neighbors.begin(), inner_neighbors.end(), ptr_neighbor));
        //     return getInnerNeighborIndex<t_i, t_j>(neighbor_index);
        // }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborOrientation(
            Integer i
        )
        const
        requires(!t_element.isNode())
        {
            // return getInnerNeighborIndex<t_i, t_j>(i) == 0 ? 1 : -1;
            auto constexpr t_inner_neighbor = t_ElementTraits::template getInnerNeighbor<t_i, t_j>();
            auto constexpr ggg = t_inner_neighbor.getDim() - 1;
            // auto constexpr t_inner_neighbor_num_nodes = ElementTraits<t_inner_neighbor, t_domain>::template getNumInnerNeighbors<ggg, 0>();
            auto ori = 1;
            for (auto node_tag = 0; node_tag < ElementTraits<t_inner_neighbor, t_domain>::template getNumInnerNeighbors<ggg, 0>(); node_tag++)
            {
                // auto mmm = FiniteElementHolder<t_inner_neighbor, t_domain>::template getInnerNeighborNodeConnection<ggg, 0>(node_tag, 0);
                auto lll = getInnerNeighbors<t_i, t_j>()[i]->getCurrentCoordinates(node_tag);
                auto kkk = getCurrentCoordinates(getInnerNeighborNodeConnection<t_i, t_j>(i, node_tag));
                if (!lll.isApprox(kkk))
                {
                    ori = -1;
                }
            }
            return ori;
        }
        
        // template<Integer t_i, Integer t_j>
        // Integer
        // getInnerNeighborOrientation(
        //     std::shared_ptr<FiniteElementHolder<ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        // )
        // const
        // requires(!t_element.isNode())
        // {
        //     return getInnerNeighborIndex<t_i, t_j>(ptr_neighbor) == 0 ? 1 : -1;
        // }
    
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
    
        Point &
        getCurrentCoordinates(
            Integer node_tag
        )
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        Point const &
        getCurrentCoordinates(
            Integer node_tag
        )
        const
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
    
        Point &
        getCurrentCoordinates(
            Integer node_tag
        )
        {
            return this->template getInnerNeighbors<t_element.getDim() - 1, 0>()[node_tag]->getCurrentCoordinates();
        }
        
        Point const &
        getCurrentCoordinates(
            Integer node_tag
        )
        const
        {
            return this->template getInnerNeighbors<t_element.getDim() - 1, 0>()[node_tag]->getCurrentCoordinates();
        }
    
        static
        algebra::View<Point const>
        getReferenceCoordinates(
            Integer node_tag
        )
        {
            return algebra::View<Point const>(ElementTraits<t_element, t_domain>::reference_nodes_[node_tag].data());
        }
        
        Matrix<Real, 3, t_element.getNumNodes()>
        getCurrentCoordinates()
        const
        requires(!t_element.isNode())
        {
            auto current_nodes_coordinates = Matrix<Real, 3, t_element.getNumNodes()>();
            auto count = 0;
            for (auto const & node : this->template getInnerNeighbors<t_element.getDim() - 1, 0>())
            {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }
        
        static
        lolita::algebra::Span<Matrix<Real, 3, t_element.getNumNodes(), lolita::algebra::colMajor()> const>
        getReferenceCoordinates()
        {
            using t_ReferenceCoordinates = lolita::algebra::Span<Matrix<Real, 3, t_element.getNumNodes(), lolita::algebra::colMajor()> const>;
            return t_ReferenceCoordinates(ElementTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
        }
        
        static
        Real
        getShapeMappingEvaluation(
            Vector<Real, t_element.getNumNodes()> const & nodal_field_values,
            Point const & reference_point
        )
        {
            return t_ElementTraits::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }
        
        static
        Real
        getShapeMappingDerivative(
            Vector<Real, t_element.getNumNodes()> const & nodal_field_values,
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
            auto ru = Matrix<Real, 3, 3>();
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
            if constexpr (t_element.isSub(t_domain, 0))
            {
                auto const & current_coordinates = this->getCurrentCoordinates();
                auto distance = Real();
                auto mp0 = Point();
                auto mp1 = Point();
                mp0.setZero();
                mp1.setZero();
                for (auto i = 0; i < t_domain.getDim(); ++i)
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
                for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
                {
                    auto pq = SegmentQuadrature::reference_points_[q][0];
                    auto wq = SegmentQuadrature::reference_weights_[q];
                    auto ru = Matrix<Real, 3, 3>();
                    auto difference = second_point - first_point;
                    auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                    ru.setZero();
                    for (auto i = 0; i < t_domain.dim_; ++i)
                    {
                        for (auto j = 0; j < t_element.dim_; ++j)
                        {
                            if (direction == -1 || i == direction)
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
        
        Real
        getLocalFrameDistance(
            Point const & first_point,
            Point const & second_point,
            Integer kkk
        )
        const
        requires(t_element.isSub(t_domain, 0))
        {
            return (second_point - first_point)(kkk);
        }
        
        Real
        getLocalFrameDistance(
            Point const & first_point,
            Point const & second_point,
            Integer kkk
        )
        const
        requires(t_element.isSub(t_domain, 1))
        {
            auto rot = getRotationMatrix(getReferenceCentroid());
            return (rot * (second_point - first_point))(kkk);
        }
        
        // Real
        // getLocalFrameDistance(
        //     Point const & first_point,
        //     Point const & second_point,
        //     Integer kkk
        // )
        // const
        // {
        //     auto mp_0 = Point();
        //     auto mp_1 = Point();
        //     auto const & current_coordinates = this->getCurrentCoordinates();
        //     mp_0.setZero();
        //     mp_1.setZero();
        //     for (auto i = 0; i < t_domain.getDim(); ++i)
        //     {
        //         mp_0(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
        //         mp_1(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
        //     }
        //     if constexpr (t_element.isSub(t_domain, 0))
        //     {
        //         return (mp_1 - mp_0)(kkk);
        //     }
        //     else if constexpr (t_element.isSub(t_domain, 1))
        //     {
        //         auto rot = getRotationMatrix(getReferenceCentroid());
        //         return (rot * (mp_1 - mp_0))(kkk);
        //     }
        // }
        
        Point
        getLocalFrameDiameters()
        const
        requires(t_element.isSub(t_domain, 0))
        {
            return getCurrentDiameters();
        }
        
        Point
        getLocalFrameDiameters()
        const
        requires(t_element.isSub(t_domain, 1))
        {
            auto rotation_matrix = getRotationMatrix(getReferenceCentroid());
            auto coordinates = getCurrentCoordinates();
            auto projected_coordinates = rotation_matrix * coordinates;
            auto current_diameters = Point();
            current_diameters.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i)
            {
                for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
                {
                    // auto pt0 = projected_coordinates.col(i);
                    // auto pt1 = projected_coordinates.col(j);
                    for (auto k = 0; k < t_element.getDim(); ++k)
                    {
                        // auto new_value = (pt1 - pt0)(k);
                        auto new_value = projected_coordinates(k, j) - projected_coordinates(k, i);
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
        
        // Point
        // getLocalFrameDiameters()
        // const
        // {
        //     if constexpr (t_element.isSub(t_domain, 0))
        //     {
        //         return getCurrentDiameters();
        //     }
        //     else if constexpr (t_element.isSub(t_domain, 1))
        //     {
        //         auto proj_v = getRotationMatrix(getReferenceCentroid()) * getCurrentCoordinates();
        //         auto current_diameters = Point();
        //         current_diameters.setZero();
        //         for (auto i = 0; i < t_element.getNumNodes(); ++i)
        //         {
        //             for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
        //             {
        //                 auto pt0 = proj_v.col(i);
        //                 auto pt1 = proj_v.col(j);
        //                 for (auto k = 0; k < t_element.getDim(); ++k)
        //                 {
        //                     auto new_value = (pt1 - pt0)(k);
        //                     auto & current_value = current_diameters(k);
        //                     if (new_value > current_value)
        //                     {
        //                         current_value = new_value;
        //                     }
        //                 }
        //             }
        //         }
        //         return current_diameters;
        //     }
        // }
        
        Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElementHolder::getReferenceCoordinates();
            auto current_diameters = Point();
            current_diameters.setZero();
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
            auto ru = Matrix<Real, 3, t_element.dim_>();
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

        Matrix<Real, 3, 3>
        getRotationMatrix(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(3))
        {
            // getting the element centroid and the first node coordinates
            auto centroid = getCurrentCentroid();
            auto node_coordinates = getCurrentCoordinates(0);
            // setting basis vectors
            auto basis_vector_0 = getNormalVector(point);
            auto basis_vector_1 = (node_coordinates - centroid) / (node_coordinates - centroid).norm();
            auto basis_vector_2 = basis_vector_0.cross(basis_vector_1);
            // setting rotation matrix
            auto rotation_matrix = Matrix<Real, 3, 3>();
            rotation_matrix.setZero();
            rotation_matrix(0, 0) = + basis_vector_0(0);
            rotation_matrix(0, 1) = + basis_vector_0(1);
            rotation_matrix(0, 2) = + basis_vector_0(2);
            rotation_matrix(1, 0) = + basis_vector_1(0);
            rotation_matrix(1, 1) = + basis_vector_1(1);
            rotation_matrix(1, 2) = + basis_vector_1(2);
            rotation_matrix(2, 0) = + basis_vector_2(0);
            rotation_matrix(2, 1) = + basis_vector_2(1);
            rotation_matrix(2, 2) = + basis_vector_2(2);
            return rotation_matrix;
        }

        Matrix<Real, 3, 3>
        getRotationMatrix(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(2))
        {
            // setting the only basis vector needed
            auto basis_vector_0 = getNormalVector(point);
            // setting rotation matrix
            auto rotation_matrix = Matrix<Real, 3, 3>();
            rotation_matrix.setZero();
            rotation_matrix(0, 0) = - basis_vector_0(1);
            rotation_matrix(0, 1) = + basis_vector_0(0);
            rotation_matrix(1, 0) = + basis_vector_0(0);
            rotation_matrix(1, 1) = + basis_vector_0(1);
            return rotation_matrix;
        }

        Matrix<Real, 3, 3>
        getRotationMatrix(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(1))
        {
            // setting rotation matrix
            auto rotation_matrix = Matrix<Real, 3, 3>();
            rotation_matrix.setIdentity();
            return rotation_matrix;
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
            auto const constexpr t_component = t_ElementTraits::template getInnerNeighbor<_i, _j>();
            using ComponentGeometry = FiniteElementHolder<t_component, t_domain>;
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
            auto constexpr t_component = t_ElementTraits ::template getInnerNeighbor<_i, _j>();
            auto p = Point();
            using ComponentGeometry = FiniteElementHolder<t_component, t_domain>;
            auto const & elt_reference_nodes = ElementTraits<t_element, t_domain>::reference_nodes_;
            for (auto i = 0; i < 3; ++i)
            {
                auto cpt_coordinates = Vector<Real, t_component.getNumNodes()>();
                for (auto j = 0; j < t_component.getNumNodes(); ++j)
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
        // ICI
        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        
        template<Basis t_basis>
        Vector<Real, t_Basis<t_basis>::getSize()>
        getBasisEvaluation(
            auto const &... args
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisEvaluation(args...);
        }
        
        template<Basis t_basis>
        Vector<Real, t_Basis<t_basis>::getSize()>
        getBasisDerivative(
            auto const &... args
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisDerivative(args...);
        }

        // template<Basis t_row_basis, Basis t_col_basis, Quadrature t_quadrature>
        // Matrix<Real, t_Basis<t_row_basis>::getSize(), t_Basis<t_col_basis>::getSize()>
        // getMassMatrix()
        // const
        // {
        //     auto mass_matrix = Matrix<Real, t_Basis<t_row_basis>::getSize(), t_Basis<t_col_basis>::getSize()>();
        //     mass_matrix.setZero();
        //     for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
        //     {
        //         auto point = getCurrentQuadraturePoint<t_quadrature>(i);
        //         auto weight = getCurrentQuadratureWeight<t_quadrature>(i);
        //         auto row_vector = getBasisEvaluation<t_row_basis>(point);
        //         auto col_vector = getBasisEvaluation<t_col_basis>(point);
        //         mass_matrix += weight * row_vector * col_vector.transpose();
        //     }
        //     return mass_matrix;
        // }

        // template<Basis t_row_basis, Basis t_col_basis>
        // Matrix<Real, t_Basis<t_row_basis>::getSize(), t_Basis<t_col_basis>::getSize()>
        // getMassMatrix(
        //     std::basic_string<Character> const & label
        // )
        // const
        // {
        //     auto mass_matrix = Matrix<Real, t_Basis<t_row_basis>::getSize(), t_Basis<t_col_basis>::getSize()>();
        //     mass_matrix.setZero();
        //     for (auto const & ip : quadrature_.at(label).ips_)
        //     {
        //         auto point = ip.coordinates_;
        //         auto weight = ip.weight_;
        //         auto row_vector = getBasisEvaluation<t_row_basis>(point);
        //         auto col_vector = getBasisEvaluation<t_col_basis>(point);
        //         mass_matrix += weight * row_vector * col_vector.transpose();
        //     }
        //     return mass_matrix;
        // }

        template<Field t_field, Mapping t_mapping, auto t_discretization>
        Matrix<Real, MappingTraits<t_mapping>::template getSize<t_domain, t_field>(), t_Disc<t_discretization>::template getNumElementUnknowns<t_field>()>
        getMapping(
            Point const & point
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getMapping<t_mapping, t_field>(point);
        }

        template<Field t_field, auto t_discretization>
        Matrix<Real, t_Disc<t_discretization>::template getNumElementUnknowns<t_field>(), t_Disc<t_discretization>::template getNumElementUnknowns<t_field>()>
        getStabilization()
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getStabilization<t_field>();
        }

        template<Field t_field, auto t_discretization>
        auto
        getGradientRhs(
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getGradientRhs<t_field>(row, col);
        }

        template<Field t_field, auto t_discretization>
        auto
        getPotential(
            Point const & point
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getPotential<t_field>(point);
        }

        template<Field t_field, auto t_discretization>
        auto
        getSymmetricGradientRhs(
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getSymmetricGradientRhs<t_field>(row, col);
        }

        template<Field t_field, auto t_discretization>
        auto
        getSymmetricGradientRhsDEBUG(
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getSymmetricGradientRhsDEBUG<t_field>(row, col);
        }

        template<Field t_field, auto t_discretization>
        Vector<Real, t_Disc<t_discretization>::template getNumElementUnknowns<t_field>()>
        getUnknowns(
            std::basic_string_view<Character> label
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getUnknowns<t_field>(label);
        }

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // LOAD
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        // struct Prescription
        // {

        //     Prescription()
        //     {}

        //     Prescription(
        //         std::shared_ptr<Loading> const & function,
        //         Integer row,
        //         Integer col
        //     )
        //     :
        //     function_(function),
        //     row_(row),
        //     col_(col)
        //     {}
            
        //     inline
        //     Boolean
        //     operator==(
        //         Prescription const & other
        //     )
        //     const = default;
            
        //     inline
        //     Boolean
        //     operator!=(
        //         Prescription const & other
        //     )
        //     const = default;
            
        //     inline
        //     Real
        //     getImposedValue(
        //         Point const & point,
        //         Real const & time
        //     )
        //     const
        //     {
        //         return function_->operator()(point, time);
        //     }
            
        //     Integer
        //     getRow()
        //     const
        //     {
        //         return row_;
        //     }
            
        //     Integer
        //     getCol()
        //     const
        //     {
        //         return col_;
        //     }

        //     Index row_;

        //     Index col_;

        //     std::shared_ptr<Loading> function_;

        // };
        
        // void
        // setLoad(
        //     std::basic_string<Character> const & label,
        //     std::shared_ptr<Loading> const & function,
        //     Integer row,
        //     Integer col
        // )
        // {
        //     loads_[label] = Prescription(function, row, col);
        // }
        
        // void
        // setConstraint(
        //     std::basic_string<Character> const & label,
        //     std::shared_ptr<Loading> const & function,
        //     Integer row,
        //     Integer col
        // )
        // {
        //     constraints_[label] = Prescription(function, row, col);
        // }
        
        // std::map<std::basic_string<Character>, Prescription> loads_;
        
        // std::map<std::basic_string<Character>, Prescription> constraints_;

        struct Load2
        {

            Load2()
            {}

            Load2(
                std::shared_ptr<Function> const & function
            )
            :
            function_(function)
            {}
            
            Boolean
            operator==(
                Load2 const & other
            )
            const = default;
            
            Boolean
            operator!=(
                Load2 const & other
            )
            const = default;

            Function const &
            getFunction()
            const
            {
                return * function_;
            }

            Function &
            getFunction()
            {
                return * function_;
            }

            std::shared_ptr<Function> function_;

        };
        
        void
        setLoad(
            std::basic_string<Character> const & label,
            std::shared_ptr<Function> const & function
        )
        {
            loads_[label] = Load2(function);
        }
        
        void
        setConstraint(
            std::basic_string<Character> const & label,
            std::shared_ptr<Function> const & function
        )
        {
            constraints_[label] = Load2(function);
        }
        
        std::map<std::basic_string<Character>, Load2> loads_;
        
        std::map<std::basic_string<Character>, Load2> constraints_;

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // DOF
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        struct DegreeOfFreedom
        {

            template<Field t_field, Basis t_basis>
            static constexpr
            Integer
            getSize()
            {
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
                return t_field_size * t_basis_size;
            }

            template<Field t_field, Basis t_basis>
            static
            DegreeOfFreedom
            make(
                std::shared_ptr<Vector<Real>> const & dof
            )
            {
                auto element_degree_of_freedom = DegreeOfFreedom(dof->size(), dof);
                dof->resize(dof->size() + getSize<t_field, t_basis>());
                return element_degree_of_freedom;
            }

            DegreeOfFreedom()
            {}

            DegreeOfFreedom(
                Natural tag,
                std::shared_ptr<Vector<Real>> const & dof
            )
            :
            tag_(tag),
            dof_(dof)
            {}
            
            Boolean
            operator==(
                DegreeOfFreedom const & other
            )
            const = default;
            
            Boolean
            operator!=(
                DegreeOfFreedom const & other
            )
            const = default;
            
            Natural
            getTag()
            const
            {
                return tag_;
            }

            template<Field t_field, Basis t_basis>
            lolita::algebra::View<Vector<Real, getSize<t_field, t_basis>()>>
            getCoefficients()
            {
                auto const & data = dof_->data() + tag_;
                return lolita::algebra::View<Vector<Real, getSize<t_field, t_basis>()>>(data);
            }

            template<Field t_field, Basis t_basis>
            lolita::algebra::View<Vector<Real, getSize<t_field, t_basis>()> const>
            getCoefficients()
            const
            {
                auto const & data = dof_->data() + tag_;
                return lolita::algebra::View<Vector<Real, getSize<t_field, t_basis>()> const>(data);
            }

            template<Field t_field, Basis t_basis>
            lolita::algebra::View<Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>
            getCoefficients(
                Integer row,
                Integer col
            )
            {
                auto constexpr field_num_cols = FieldTraits<t_field>::template getCols<t_domain>();
                auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
                auto const & data = dof_->data() + tag_ + (field_num_cols * row  + col) * t_basis_size;
                return lolita::algebra::View<Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
            }

            template<Field t_field, Basis t_basis>
            lolita::algebra::View<Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()> const>
            getCoefficients(
                Integer row,
                Integer col
            )
            const
            {
                auto constexpr field_num_cols = FieldTraits<t_field>::template getCols<t_domain>();
                auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
                auto const & data = dof_->data() + tag_ + (field_num_cols * row  + col) * t_basis_size;
                return lolita::algebra::View<Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()> const>(data);
            }

            template<Field t_field, Basis t_basis>
            Vector<Real, getSize<t_field, t_basis>()>
            getCoefficientsCopy()
            const
            {
                return Vector<Real, getSize<t_field, t_basis>()>(getCoefficients<t_field, t_basis>());
            }

            template<Field t_field, Basis t_basis>
            Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>
            getCoefficientsCopy(
                Integer row,
                Integer col
            )
            const
            {
                return Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>(getCoefficients<t_field, t_basis>(row, col));
            }

            template<Field t_field, Basis t_basis>
            void
            addCoefficients(
                lolita::algebra::View<Vector<Real> const> const & vector
            )
            {
                getCoefficients<t_field, t_basis>() += vector.template segment<getSize<t_field, t_basis>()>(tag_);
            }

            template<Field t_field, Basis t_basis>
            void
            setCoefficients(
                lolita::algebra::View<Vector<Real> const> const & vector
            )
            {
                getCoefficients<t_field, t_basis>() = vector.template segment<getSize<t_field, t_basis>()>(tag_);
            }

            Natural tag_;

            std::shared_ptr<Vector<Real>> dof_;

        };

        template<Field t_field, Basis t_basis>
        void
        setDegreeOfFreedom(
            std::basic_string<Character> const & label,
            std::shared_ptr<Vector<Real>> const & dof
        )
        {
            degrees_of_freedom_[label] = DegreeOfFreedom::template make<t_field, t_basis>(dof);
        }

        template<Field t_field, Basis t_basis>
        void
        addDegreeOfFreedomCoefficients(
            std::basic_string<Character> const & label,
            lolita::algebra::View<Vector<Real> const> const & vector
        )
        {
            degrees_of_freedom_.at(label).template addCoefficients<t_field, t_basis>(vector);
        }

        template<Field t_field, Basis t_basis>
        void
        setDegreeOfFreedomCoefficients(
            std::basic_string<Character> const & label,
            lolita::algebra::View<Vector<Real> const> const & vector
        )
        {
            degrees_of_freedom_.at(label).template setCoefficients<t_field, t_basis>(vector);
        }

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        updateUnknown(
            std::basic_string_view<Character> unknown_label,
            std::unique_ptr<System> const & system
        )
        {
            static_cast<t_Disc<t_discretization> *>(this)->template updateUnknown<t_finite_element_method>(unknown_label, system);
        }

        template<Field t_field, Basis t_basis>
        void
        updateUnknown(
            std::basic_string_view<Character> unknown_label,
            std::unique_ptr<System> const & system
        )
        {
            degrees_of_freedom_.at(std::string(unknown_label)).template addCoefficients<t_field, t_basis>(system->getUnknownCorrection(unknown_label));
        }

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        updateBinding(
            std::basic_string_view<Character> binding_label,
            std::unique_ptr<System> const & system
        )
        {
            static_cast<t_Disc<t_discretization> *>(this)->template updateBinding<t_finite_element_method>(binding_label, system);
            // constexpr bool has_method = requires(t_Disc<t_discretization> & t) {
            //     t.hello();
            // };
            // if constexpr (has_method)
            // {
            //     std::cout << "HAS !\n";
            // }
            // else
            // {
            //     std::cout << "HAS NOT !\n";
            // }
        }

        template<Field t_field, Basis t_basis>
        void
        updateBinding(
            std::basic_string_view<Character> binding_label,
            std::unique_ptr<System> const & system
        )
        {
            degrees_of_freedom_.at(std::string(binding_label)).template addCoefficients<t_field, t_basis>(system->getBindingCorrection(binding_label));
        }

        // template<Field t_field, Basis t_basis>
        // Real
        // getUnknownValue(
        //     std::basic_string_view<Character> binding_label,
        //     Point const & point,
        //     Integer row,
        //     Integer col
        // )
        // const
        // {
        //     auto coefficients = degrees_of_freedom_.at(std::string(binding_label)).template getCoefficients<t_field, t_basis>(row, col);
        //     auto basis_vector = this->getBasisEvaluation<t_basis>(point);
        //     return coefficients.dot(basis_vector);
        // }

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        Real
        getUnknownValue(
            std::basic_string_view<Character> binding_label,
            Point const & point,
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getUnknownValue<t_finite_element_method>(binding_label, point, row, col);
        }

        // template<Field t_field, Basis t_basis>
        // Real
        // getBindingValue(
        //     std::basic_string_view<Character> binding_label,
        //     Point const & point,
        //     Integer row,
        //     Integer col
        // )
        // const
        // {
        //     auto coefficients = degrees_of_freedom_.at(std::string(binding_label)).template getCoefficients<t_field, t_basis>(row, col);
        //     auto basis_vector = this->getBasisEvaluation<t_basis>(point);
        //     return coefficients.dot(basis_vector);
        // }

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        Real
        getBindingValue(
            std::basic_string_view<Character> binding_label,
            Point const & point,
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getBindingValue<t_finite_element_method>(binding_label, point, row, col);
        }
        
        std::map<std::basic_string<Character>, DegreeOfFreedom> degrees_of_freedom_;

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // QUAD
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        struct QuadratureElement
        {

            struct IntegrationPoint
            {

            private:

                template<FiniteElementMethodConcept auto t_finite_element_method>
                static constexpr
                Integer
                getGeneralizedStrainSize()
                {
                    return FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
                }

                template<BehaviorConcept auto t_behavior>
                static constexpr
                Integer
                getGeneralizedStrainSize()
                {
                    return BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
                }

            public:

                IntegrationPoint(
                    algebra::View<Point const> ref_pt,
                    Point const & coordinates,
                    Real const & weight,
                    std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
                )
                :
                reference_coordinates_(ref_pt),
                coordinates_(coordinates),
                weight_(weight),
                behavior_(behavior),
                behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behavior))),
                behavior_data_view_(std::make_unique<mgis::behaviour::BehaviourDataView>(mgis::behaviour::make_view(* behavior_data_)))
                {
                    behavior_data_->K[0] = 4;
                }
            
                inline
                Boolean
                operator==(
                    IntegrationPoint const & other
                )
                const = default;
                
                inline
                Boolean
                operator!=(
                    IntegrationPoint const & other
                )
                const = default;

                void
                integrate()
                {
                    auto behaviour_data_view = mgis::behaviour::make_view(* behavior_data_);
                    auto res = mgis::behaviour::integrate(behaviour_data_view, * behavior_);
                    // auto strain_view = lolita::algebra::View<Vector<Real> const>(behavior_data_->s1.gradients.data(), behavior_data_->s1.gradients.size());
                    // auto stress_view = lolita::algebra::View<Vector<Real> const>(behavior_data_->s1.thermodynamic_forces.data(), behavior_data_->s1.thermodynamic_forces.size());
                    // auto K = lolita::algebra::View<Vector<Real> const>(behavior_data_->K.data(), behavior_data_->K.size());
                    // std::cout << "strain : " << strain_view << std::endl;
                    // std::cout << "stress : " << stress_view << std::endl;
                    // std::cout << "K : " << K << std::endl;
                    // std::cout << "res : " << res << std::endl;
                }

                void
                setMaterialProperty(
                    std::basic_string_view<Character> material_property_label,
                    std::function<Real(Point const &)> && function
                )
                {
                    auto value = std::forward<std::function<Real(Point const &)>>(function)(coordinates_);
                    mgis::behaviour::setMaterialProperty(behavior_data_->s0, std::string(material_property_label), value);
                    mgis::behaviour::setMaterialProperty(behavior_data_->s1, std::string(material_property_label), value);
                }

                void
                setExternalVariable(
                    std::basic_string_view<Character> material_property_label,
                    std::function<Real(Point const &)> && function
                )
                {
                    auto value = std::forward<std::function<Real(Point const &)>>(function)(coordinates_);
                    mgis::behaviour::setExternalStateVariable(behavior_data_->s0, std::string(material_property_label), value);
                    mgis::behaviour::setExternalStateVariable(behavior_data_->s1, std::string(material_property_label), value);
                }
                
                template<BehaviorConcept auto t_behavior>
                lolita::algebra::Span<Vector<Real, BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>()> const>
                getGeneralizedStrain()
                const
                {
                    auto constexpr size = BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
                    return lolita::algebra::Span<Vector<Real, size> const>(behavior_data_->s1.gradients.data());
                }
                
                template<auto t_finite_element_method>
                lolita::algebra::Span<Vector<Real, FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()> const>
                getGeneralizedStrain()
                const
                {
                    auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
                    auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
                    return lolita::algebra::Span<Vector<Real, size> const>(behavior_data_->s1.gradients.data() + offset);
                }
                
                template<auto t_finite_element_method>
                lolita::algebra::Span<Vector<Real, FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()>>
                getGeneralizedStrain()
                {
                    auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
                    auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
                    return lolita::algebra::Span<Vector<Real, size>>(behavior_data_->s1.gradients.data() + offset);
                }

                template<FiniteElementMethodConcept auto t_finite_element_method>
                Matrix<Real, getGeneralizedStrainSize<t_finite_element_method>(), getGeneralizedStrainSize<t_finite_element_method>()>
                getJacobian()
                const
                {
                    auto constexpr strain_operator_num_rows = getGeneralizedStrainSize<t_finite_element_method>();
                    using Jac = Matrix<Real, strain_operator_num_rows, strain_operator_num_rows>;
                    auto jac = Jac();
                    jac.setZero();
                    auto set_mapping_block = [&] <Integer t_i = 0> (
                        auto & self
                    )
                    constexpr mutable
                    {
                        auto constexpr mapping = t_finite_element_method.template getMapping<t_i>();
                        auto constexpr mapping_size = FiniteElementMethodTraits<t_finite_element_method>::template getMappingSize<t_domain, mapping>();
                        auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getMappingOffset<t_domain, mapping>();
                        auto jacobian_block_view = algebra::View<Matrix<Real, mapping_size, mapping_size> const>(behavior_data_->K.data() + offset * offset);
                        auto tangent_block_view = jac.template block<mapping_size, mapping_size>(offset, offset);
                        tangent_block_view = jacobian_block_view;
                        if constexpr (t_i < t_finite_element_method.getGeneralizedStrain().getNumMappings() - 1)
                        {
                            self.template operator ()<t_i + 1>(self);
                        }
                    };
                    set_mapping_block(set_mapping_block);
                    return jac;
                }
    
                Point &
                getCurrentCoordinates()
                {
                    return coordinates_;
                }
                
                Point const &
                getCurrentCoordinates()
                const
                {
                    return coordinates_;
                }
                
                algebra::View<Point const>
                getReferenceCoordinates()
                const
                {
                    return reference_coordinates_;
                }

                algebra::View<Point const> reference_coordinates_;
            
                Point coordinates_;
                
                Real weight_;
            
                std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

                std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

                std::unique_ptr<mgis::behaviour::BehaviourDataView> behavior_data_view_;

                std::map<std::basic_string<Character>, Matrix<Real>> ops_;

            };

            static inline
            QuadratureElement
            make(
                std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
            )
            {
                return QuadratureElement(behavior);
            }

            QuadratureElement()
            :
            ips_(),
            behavior_()
            {}

            explicit
            QuadratureElement(
                std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
            )
            :
            ips_(),
            behavior_(behavior)
            {}
            
            inline
            Boolean
            operator==(
                QuadratureElement const & other
            )
            const = default;
            
            inline
            Boolean
            operator!=(
                QuadratureElement const & other
            )
            const = default;

            std::vector<IntegrationPoint> ips_;
            
            std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        };

        template<Quadrature t_quadrature>
        void
        setBehavior(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            auto behavior_label = behavior->behaviour;
            quadrature_[behavior_label] = QuadratureElement::make(behavior);
            for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
            {
                auto point = getCurrentQuadraturePoint<t_quadrature>(i);
                auto r_point = getReferenceQuadraturePoint<t_quadrature>(i);
                auto weight = getCurrentQuadratureWeight<t_quadrature>(i);
                auto ip = typename QuadratureElement::IntegrationPoint(r_point, point, weight, quadrature_.at(behavior_label).behavior_);
                quadrature_.at(behavior_label).ips_.push_back(std::move(ip));
            }
            
        }

        // template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        // void
        // setStrainOperators(
        //     std::basic_string_view<Character> label,
        //     std::basic_string_view<Character> label2
        // )
        // {
        //     auto constexpr strain_operator_num_rows = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
        //     auto constexpr strain_operator_num_cols = t_Disc<t_discretization>::template getNumElementUnknowns<t_finite_element_method.getField()>();
        //     auto quadrature_point_count = 0;
        //     for (auto & ip : quadrature_.at(std::string(label)).ips_)
        //     {
        //         auto strain_operator = Matrix<Real, strain_operator_num_rows, strain_operator_num_cols>();
        //         strain_operator.setZero();
        //         auto set_mapping_block = [&] <Integer t_i = 0> (
        //             auto & self
        //         )
        //         constexpr mutable
        //         {
        //             auto constexpr mapping = t_finite_element_method.template getMapping<t_i>();
        //             auto constexpr mapping_size = FiniteElementMethodTraits<t_finite_element_method>::template getMappingSize<t_domain, mapping>();
        //             auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getMappingOffset<t_domain, mapping>();
        //             auto point = ip.reference_coordinates_;
        //             auto mapping_operator = this->template getMapping<t_finite_element_method.getField(), mapping, t_discretization>(point);
        //             auto mapping_block = strain_operator.template block<mapping_size, strain_operator_num_cols>(offset, 0);
        //             mapping_block = mapping_operator;
        //             if constexpr (t_i < t_finite_element_method.getGeneralizedStrain().getNumMappings() - 1)
        //             {
        //                 self.template operator ()<t_i + 1>(self);
        //             }
        //         };
        //         set_mapping_block(set_mapping_block);
        //         ip.ops_[std::string(label2)] = strain_operator;
        //         quadrature_point_count ++;
        //     }
        // }

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        setStrainOperators(
            std::basic_string_view<Character> label,
            std::basic_string_view<Character> label2
        )
        {
            auto constexpr strain_operator_num_rows = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            auto constexpr strain_operator_num_cols = t_Disc<t_discretization>::template getNumElementUnknowns<t_finite_element_method.getField()>();
            auto quadrature_point_count = 0;
            for (auto & ip : quadrature_.at(std::string(label)).ips_)
            {
                auto strain_operator = Matrix<Real, strain_operator_num_rows, strain_operator_num_cols>();
                strain_operator.setZero();
                auto set_mapping_block = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr mapping = t_finite_element_method.template getMapping<t_i>();
                    auto constexpr mapping_size = FiniteElementMethodTraits<t_finite_element_method>::template getMappingSize<t_domain, mapping>();
                    auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getMappingOffset<t_domain, mapping>();
                    auto const & point = ip.coordinates_;
                    auto mapping_operator = this->template getMapping<t_finite_element_method.getField(), mapping, t_discretization>(point);
                    auto mapping_block = strain_operator.template block<mapping_size, strain_operator_num_cols>(offset, 0);
                    mapping_block = mapping_operator;
                    if constexpr (t_i < t_finite_element_method.getGeneralizedStrain().getNumMappings() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_mapping_block(set_mapping_block);
                ip.ops_[std::string(label2)] = strain_operator;
                quadrature_point_count ++;
            }
        }

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        setStrainValues(
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> degree_of_freedom_label
        )
        {
            auto constexpr strain_operator_num_rows = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            auto constexpr strain_operator_num_cols = t_Disc<t_discretization>::template getNumElementUnknowns<t_finite_element_method.getField()>();
            using StrainOperatorView = lolita::algebra::View<Matrix<Real, strain_operator_num_rows, strain_operator_num_cols> const>;
            using StrainView = lolita::algebra::View<Vector<Real, strain_operator_num_rows>>;
            using StrainValue = Vector<Real, strain_operator_num_rows>;
            auto unknown = this->template getUnknowns<t_finite_element_method.getField(), t_discretization>(degree_of_freedom_label);
            for (auto & ip : quadrature_.at(std::string(behavior_label)).ips_)
            {
                auto strain_operator = StrainOperatorView(ip.ops_.at(std::string(degree_of_freedom_label)).data());
                auto strain_value = StrainValue(strain_operator * unknown);
                auto strain_view = StrainView(ip.behavior_data_->s1.gradients.data());
                strain_view = strain_value;
                strain_view.setZero();
            }
        }

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        assembleUnknownBlock(
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> degree_of_freedom_label,
            std::unique_ptr<System> const & system
        )
        {
            static_cast<t_Disc<t_discretization> *>(this)->template assembleUnknownBlock<t_finite_element_method>(behavior_label, degree_of_freedom_label, system);
        }

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        assembleBindingBlock(
            // std::basic_string_view<Character> behavior_label,
            // std::basic_string_view<Character> degree_of_freedom_label,
            // std::unique_ptr<System> const & system
            std::basic_string_view<Character> binding_label,
            std::basic_string_view<Character> unknown_label,
            std::basic_string_view<Character> constraint_label,
            std::unique_ptr<System> const & system
        )
        const
        {
            static_cast<t_Disc<t_discretization> const *>(this)->template assembleBindingBlock<t_finite_element_method>(binding_label, unknown_label, constraint_label, system);
        }

        void
        integrate(
            std::basic_string_view<Character> behavior_label
        )
        {
            for (auto & ip : quadrature_.at(std::string(behavior_label)).ips_)
            {
                ip.integrate();
            }
        }

        void
        setMaterialProperty(
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            for (auto & ip : quadrature_.at(std::string(behavior_label)).ips_)
            {
                ip.setMaterialProperty(material_property_label, std::forward<std::function<Real(Point const &)>>(function));
            }
        }

        void
        setExternalVariable(
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            for (auto & ip : quadrature_.at(std::string(behavior_label)).ips_)
            {
                ip.setExternalVariable(material_property_label, std::forward<std::function<Real(Point const &)>>(function));
            }
        }

        void
        setParameter(
            std::basic_string<Character> const & parameter_label,
            std::function<Real(Point const &)> && function
        )
        {
            parameters_[parameter_label] = std::forward<std::function<Real(Point const &)>>(function)(* coordinates_);
        }

        Integer
        getNumIntegrationPoints(
            std::basic_string_view<Character> behavior_label
        )
        const
        {
            return quadrature_.contains(std::string(behavior_label)) ? quadrature_.at(std::string(behavior_label)).ips_.size() : 0;
        }

        Integer
        getNumIntegrationPoints()
        const
        {
            auto num_integration_points = 0;
            for (auto const & quadrature : quadrature_)
            {
                num_integration_points += quadrature.second.ips_.size();
            }
            return num_integration_points;
        }

        std::map<std::basic_string<Character>, QuadratureElement> quadrature_;

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // ELEM OPS
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        template<FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        setElementOperator(
            std::basic_string_view<Character> label
        )
        {
            operators_[std::string(label)] = this->template getStabilization<t_finite_element_method.getField(), t_discretization>();
        }

        std::map<std::basic_string<Character>, Matrix<Real>> operators_;

        std::map<std::basic_string<Character>, Real> parameters_;

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // NEW
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        Natural const &
        getTag()
        const
        {
            return tag_;
        }
        
        t_OuterNeighbors outer_neighbors_;
        
        t_InnerNeighbors inner_neighbors_;
        
        Natural tag_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;
        
        std::shared_ptr<Point> coordinates_;

    };

} // namespace lolita


#endif /* FC0FBEA8_1C07_47E6_97F1_3D709687C47A */
