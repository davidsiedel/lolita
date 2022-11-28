#ifndef B44BF9EE_1549_420E_9AF4_5820C5464430
#define B44BF9EE_1549_420E_9AF4_5820C5464430

#include "2/core/_include.hxx"
#include "2/core/region.hxx"
#include "2/core/dof.hxx"
#include "2/core/implementation/hdg.hxx"
#include "2/core/implementation/basis.hxx"

// #include "2/core/frm.hxx"
#include "2/core/element_integration_point.hxx"
#include "2/core/element_potential.hxx"
#include "2/core/element_lagrangian_implementation.hxx"
#include "2/core/element_lagrangian_interface.hxx"

#include "2/core/region_potential.hxx"
#include "2/core/region_lagrangian_implementation.hxx"
#include "2/core/region_lagrangian_interface.hxx"

namespace lolita::core
{
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    struct FiniteElement
    {

    private:
    
        using ElementTraits_ = ShapeTraits<t_element>;

        template<BasisConcept auto t_basis>
        using BasisImplementation_ = BasisImplementation<t_element, t_domain, t_basis>;

        template<auto t_discretization>
        using t_Disc = typename DiscretizationTraits<t_discretization>::template Implementation<t_element, t_domain>;

        template<auto t_field>
        using t_Disc2 = typename DiscretizationTraits<t_field.getDiscretization()>::template Implementation<t_element, t_domain, t_field>;

        template<FieldConcept auto field_>
        using DiscretizationImplementation_ = DiscretizationImplementation<t_element, t_domain, field_>;
        
        template<ShapeConcept auto t__element, MeshConcept auto t__domain>
        using NeighborPointer_ = std::shared_ptr<FiniteElement<t__element, t__domain>>;
        
        template<DomainConcept auto t_dim, MeshConcept auto t__domain>
        using DomainPtr_ = std::shared_ptr<FiniteDomain<t_dim, t__domain>>;

        using Domain_ = typename ElementTraits_::template DomainConnectivity<DomainPtr_, t_domain>;
    
        using InnerNeighbors_ = typename ElementTraits_::template InnerConnectivity<NeighborPointer_, t_domain>;
        
        using OuterNeighbors_ = typename ElementTraits_::template OuterConnectivity<NeighborPointer_, t_domain>;

        using t_Fld = DegreesOfFreedomInterface<t_element, t_domain>;

        template<FieldConcept auto t_lag>
        using t_FldImpl = DegreesOfFreedomImplementation<t_element, t_domain, t_lag>;

        using t_Lag = AbstractElementLagrangian<t_element, t_domain>;

        template<LagrangianConcept auto t_lag>
        using t_LagImpl = ElementLagrangian<t_element, t_domain, t_lag>;
        
        Natural tag_;
        
        Point coordinates_;

        Domain_ domain_;
        
        OuterNeighbors_ outer_neighbors_;
        
        InnerNeighbors_ inner_neighbors_;

    public:

        explicit
        FiniteElement(
            Natural const & tag
        )
        :
        tag_(tag),
        coordinates_(),
        domain_(),
        outer_neighbors_(),
        inner_neighbors_()
        {}
        
        Boolean
        operator==(
            FiniteElement const & other
        )
        const
        = default;
        
        Boolean
        operator!=(
            FiniteElement const & other
        )
        const
        = default;
        
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
        requires(t_element == Node())
        {
            return std::to_string(this->tag_);
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, OuterNeighbors_>> &
        getOuterNeighbors()
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, OuterNeighbors_>> const &
        getOuterNeighbors()
        const
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, InnerNeighbors_>> &
        getInnerNeighbors()
        requires(t_element != Node())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, InnerNeighbors_>> const &
        getInnerNeighbors()
        const
        requires(t_element != Node())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }

        Point const &
        getCoordinates()
        const
        {
            return coordinates_;
        }
        
        template<BasisConcept auto t_basis>
        DenseVector<Real, BasisTraits<t_basis>::template getSize<t_element>()>
        getBasisEvaluation(
            auto const &... args
        )
        const
        {
            return static_cast<BasisImplementation_<t_basis> const *>(this)->getBasisEvaluation(args...);
        }
        
        template<BasisConcept auto t_basis>
        DenseVector<Real, BasisTraits<t_basis>::template getSize<t_element>()>
        getBasisDerivative(
            auto const &... args
        )
        const
        {
            return static_cast<BasisImplementation_<t_basis> const *>(this)->getBasisDerivative(args...);
        }

        Natural const &
        getTag()
        const
        {
            return tag_;
        }

        void
        setCoordinates(
            Point const & point
        )
        {
            coordinates_ = point;
        }

        void
        setDomain(
            Domain_ const & domain
        )
        {
            domain_ = domain;
        }

        Boolean
        hasDomain()
        const
        {
            return domain_ != nullptr;
        }

        Domain_ const &
        getDomain()
        const
        {
            return domain_;
        }

        /**
         * This part is dedicated to the Dof implementation
         * *****************************************************************************************************************************************************
         */

        template<FieldConcept auto t_lag>
        void
        setDiscreteField()
        {
            if (fields_ == nullptr)
            {
                fields_ = std::make_unique<std::vector<std::unique_ptr<t_Fld>>>();
            }
            for (auto & kk : * fields_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    kk = std::make_unique<t_FldImpl<t_lag>>(* this);
                    return;
                }
            }
            fields_->push_back(std::make_unique<t_FldImpl<t_lag>>(* this));
        }

        template<FieldConcept auto t_lag>
        Boolean
        hasDiscreteField()
        const
        {
            if (fields_ == nullptr)
            {
                return false;
            }
            for (auto const & kk : * fields_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return true;
                }
            }
            return false;
        }

        template<FieldConcept auto t_lag>
        t_Fld const &
        getDiscreteField()
        const
        {
            if (fields_ == nullptr)
            {
                throw std::runtime_error("NO");
            }
            for (auto const & kk : * fields_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return * kk;
                }
            }
            throw std::runtime_error("NO");
        }

        template<FieldConcept auto t_lag>
        t_Fld &
        getDiscreteField()
        {
            if (fields_ == nullptr)
            {
                throw std::runtime_error("NO");
            }
            for (auto & kk : * fields_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return * kk;
                }
            }
            throw std::runtime_error("NO");
        }

        std::unique_ptr<std::vector<std::unique_ptr<t_Fld>>> fields_;

        /**
         * Formulation
         * *****************************************************************************************************************************************************
         */

        template<LagrangianConcept auto t_lag>
        void
        setLagrangian()
        {
            if (lags_ == nullptr)
            {
                lags_ = std::make_unique<std::vector<std::unique_ptr<t_Lag>>>();
            }
            for (auto & kk : * lags_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    kk = std::make_unique<t_LagImpl<t_lag>>(* this);
                    return;
                }
            }
            lags_->push_back(std::make_unique<t_LagImpl<t_lag>>(* this));
        }

        template<LagrangianConcept auto t_lag>
        Boolean
        hasLagrangian()
        const
        {
            if (lags_ == nullptr)
            {
                return false;
            }
            for (auto const & kk : * lags_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return true;
                }
            }
            return false;
        }

        template<LagrangianConcept auto t_lag>
        t_Lag const &
        getLagrangian()
        const
        {
            if (lags_ == nullptr)
            {
                throw std::runtime_error("NO");
            }
            for (auto const & kk : * lags_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return * kk;
                }
            }
            throw std::runtime_error("NO");
        }

        template<LagrangianConcept auto t_lag>
        t_Lag &
        getLagrangian()
        {
            if (lags_ == nullptr)
            {
                throw std::runtime_error("NO");
            }
            for (auto & kk : * lags_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return * kk;
                }
            }
            throw std::runtime_error("NO");
        }

        std::unique_ptr<std::vector<std::unique_ptr<t_Lag>>> lags_;

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer component_index,
            Integer node_index
        )
        requires(t_element != Node())
        {
            return std::get<t_j>(std::get<t_i>(ShapeTraits<t_element>::node_connectivity_))[component_index][node_index];
        }

        template<Integer t_i>
        Boolean
        isIn(
            std::basic_string<Character> && domain
        )
        const
        requires(t_i == t_element.getDim())
        {
            if (this->hasDomain())
            {
                return this->domain_->getLabel() == std::forward<std::basic_string<Character>>(domain);
            }
            return false;
        }

        template<Integer t_i>
        Boolean
        isIn(
            std::basic_string<Character> && domain
        )
        const
        requires(t_i > t_element.getDim())
        {
            auto is_in_domain = Boolean(false);
            auto find_domain = [&] <Integer t_j = 0> (
                auto & self
            )
            constexpr mutable
            {
                for (auto const & neighbor : this->template getOuterNeighbors<t_i - t_element.getDim(), t_j>())
                {
                    if (neighbor->template isIn<t_i>(std::forward<std::basic_string<Character>>(domain)))
                    {
                        is_in_domain = true;
                    }
                }
                if constexpr (t_j < ElementTraits_::template getNumOuterNeighbors<t_domain, t_i - t_element.getDim()>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }                
            };
            find_domain(find_domain);
            return is_in_domain;
        }

        template<Integer t_i>
        Boolean
        isIn(
            std::basic_string<Character> && domain
        )
        const
        requires(t_i < t_element.getDim())
        {
            return false;
        }

        Boolean
        isIn(
            std::basic_string<Character> && domain
        )
        const
        {
            auto is_in_domain = Boolean(false);
            auto find_domain = [&] <Integer t_i = t_element.getDim()> (
                auto & self
            )
            constexpr mutable
            {
                if (this->template isIn<t_i>(std::forward<std::basic_string<Character>>(domain)))
                {
                    is_in_domain = true;
                }
                if constexpr (t_i < t_domain.getDim())
                {
                    self.template operator()<t_i + 1>(self);
                }
            };
            find_domain(find_domain);
            return is_in_domain;
        }
        
        // template<Integer t_i, Integer t_j>
        // Integer
        // getInnerNeighborIndex(
        //     Integer i
        // )
        // const
        // requires(t_element != Node())
        // {
        //     auto constexpr t_inner_neighbor = ElementTraits_::template getInnerNeighbor<t_i, t_j>();
        //     auto constexpr t_coordinates = ShapeTraits<t_inner_neighbor, t_domain>::template getOuterNeighborCoordinates<t_element>();
        //     auto const & items = getInnerNeighbors<t_i, t_j>()[i]->template getOuterNeighbors<t_coordinates.dim_, t_coordinates.tag_>();
        //     // auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element)
        //     auto is_equal = [&] (std::shared_ptr<FiniteElement<t_element, t_domain>> const & ptr_element)
        //     {
        //         return * ptr_element == * this;
        //     };
        //     return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        // }
        
        // template<Integer t_i, Integer t_j>
        // Integer
        // getInnerNeighborIndex(
        //     std::shared_ptr<FiniteElement<ShapeTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        // )
        // const
        // requires(t_element != Node())
        // {
        //     auto const & inner_neighbors = getInnerNeighbors<t_i, t_j>();
        //     // auto constexpr t_inner_neighbor = ShapeTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>();
        //     // auto is_equal = [&] (std::shared_ptr<FiniteElement<t_inner_neighbor, t_domain>> const & neighbor)
        //     // {
        //     //     return * neighbor == * ptr_neighbor;
        //     // };
        //     // auto neighbor_index = std::distance(inner_neighbors.begin(), std::find_if(inner_neighbors.begin(), inner_neighbors.end(), is_equal));
        //     auto neighbor_index = std::distance(inner_neighbors.begin(), std::find(inner_neighbors.begin(), inner_neighbors.end(), ptr_neighbor));
        //     return getInnerNeighborIndex<t_i, t_j>(neighbor_index);
        // }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborSign(
            Integer i
        )
        const
        requires(t_element != Node())
        {
            // return getInnerNeighborIndex<t_i, t_j>(i) == 0 ? 1 : -1;
            auto constexpr t_inner_neighbor = ElementTraits_::template getInnerNeighbor<t_i, t_j>();
            auto constexpr ggg = t_inner_neighbor.getDim() - 1;
            auto constexpr t_num_inner_neighbor_nodes = ShapeTraits<t_inner_neighbor>::template getNumInnerNeighbors<t_inner_neighbor.getDim() - 1, 0>();
            // auto constexpr t_inner_neighbor_num_nodes = ShapeTraits<t_inner_neighbor, t_domain>::template getNumInnerNeighbors<ggg, 0>();
            auto ori = 1;
            for (auto node_tag = 0; node_tag < t_num_inner_neighbor_nodes; node_tag++)
            {
                // auto mmm = FiniteElement<t_inner_neighbor, t_domain>::template getInnerNeighborNodeConnection<ggg, 0>(node_tag, 0);
                // auto lll = getInnerNeighbors<t_i, t_j>()[i]->getCurrentCoordinates(node_tag);
                // auto kkk = getCurrentCoordinates(getInnerNeighborNodeConnection<t_i, t_j>(i, node_tag));
                // if (!lll.isApprox(kkk))
                // {
                //     ori = -1;
                // }
                // //
                auto lll = getInnerNeighbors<t_i, t_j>()[i]->template getInnerNeighbors<t_inner_neighbor.getDim() - 1, 0>()[node_tag]->getTag();
                auto kkk = getInnerNeighbors<t_element.getDim() - 1, 0>()[getInnerNeighborNodeConnection<t_i, t_j>(i, node_tag)]->getTag();
                if (lll != kkk)
                {
                    ori = -1;
                    break;
                }
            }
            return ori;
            // auto inner_neighbor = getInnerNeighbors<t_i, t_j>()[i];
            // // auto inner_neighbor_reference_centroid = inner_neighbor->getReferenceCentroid();
            // // auto inner_neighbor_current_centroid = inner_neighbor->getCurrentCentroid();
            // // auto current_centroid = this->getCurrentCentroid();
            // auto inner_neighbor_rotation_matrix = inner_neighbor->getRotationMatrix(inner_neighbor->getReferenceCentroid());
            // return (inner_neighbor_rotation_matrix * (inner_neighbor->getCurrentCentroid() - this->getCurrentCentroid()))(t_element.getDim() - 1) > 0 ? 1 : -1;
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborOrientation(
            Integer i
        )
        const
        requires(t_element != Node())
        {
            // // return getInnerNeighborIndex<t_i, t_j>(i) == 0 ? 1 : -1;
            // auto constexpr t_inner_neighbor = ElementTraits_::template getInnerNeighbor<t_i, t_j>();
            // auto constexpr ggg = t_inner_neighbor.getDim() - 1;
            // // auto constexpr t_inner_neighbor_num_nodes = ShapeTraits<t_inner_neighbor, t_domain>::template getNumInnerNeighbors<ggg, 0>();
            // auto ori = 1;
            // for (auto node_tag = 0; node_tag < ShapeTraits<t_inner_neighbor, t_domain>::template getNumInnerNeighbors<ggg, 0>(); node_tag++)
            // {
            //     // auto mmm = FiniteElement<t_inner_neighbor, t_domain>::template getInnerNeighborNodeConnection<ggg, 0>(node_tag, 0);
            //     auto lll = getInnerNeighbors<t_i, t_j>()[i]->getCurrentCoordinates(node_tag);
            //     auto kkk = getCurrentCoordinates(getInnerNeighborNodeConnection<t_i, t_j>(i, node_tag));
            //     if (!lll.isApprox(kkk))
            //     {
            //         ori = -1;
            //     }
            // }
            auto inner_neighbor = getInnerNeighbors<t_i, t_j>()[i];
            // auto inner_neighbor_reference_centroid = inner_neighbor->getReferenceCentroid();
            // auto inner_neighbor_current_centroid = inner_neighbor->getCurrentCentroid();
            // auto current_centroid = this->getCurrentCentroid();
            auto inner_neighbor_rotation_matrix = inner_neighbor->getRotationMatrix(inner_neighbor->getReferenceCentroid());
            return (inner_neighbor_rotation_matrix * (inner_neighbor->getCurrentCentroid() - this->getCurrentCentroid()))(t_element.getDim() - 1) > 0 ? 1 : -1;
        }
        
        // template<Integer t_i, Integer t_j>
        // Integer
        // getInnerNeighborOrientation(
        //     std::shared_ptr<FiniteElement<ShapeTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        // )
        // const
        // requires(t_element != Node())
        // {
        //     return getInnerNeighborIndex<t_i, t_j>(ptr_neighbor) == 0 ? 1 : -1;
        // }
    
        // Point &
        // getCurrentCoordinates()
        // requires(t_element == Node())
        // {
        //     return this->coordinates_;
        // }
        
        Point const &
        getCurrentCoordinates()
        const
        requires(t_element == Node())
        {
            return this->coordinates_;
        }
    
        // Point &
        // getCurrentCoordinates(
        //     Integer node_tag
        // )
        // requires(t_element == Node())
        // {
        //     return this->coordinates_;
        // }
        
        Point const &
        getCurrentCoordinates(
            Integer node_tag
        )
        const
        requires(t_element == Node())
        {
            return this->coordinates_;
        }
    
        // Point &
        // getCurrentCoordinates(
        //     Integer node_tag
        // )
        // {
        //     return this->template getInnerNeighbors<t_element.getDim() - 1, 0>()[node_tag]->getCurrentCoordinates();
        // }
        
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
            return algebra::View<Point const>(ShapeTraits<t_element>::reference_nodes_[node_tag].data());
        }
        
        DenseMatrix<Real, 3, t_element.getNumNodes()>
        getCurrentCoordinates()
        const
        requires(t_element != Node())
        {
            auto current_nodes_coordinates = DenseMatrix<Real, 3, t_element.getNumNodes()>();
            auto count = 0;
            for (auto const & node : this->template getInnerNeighbors<t_element.getDim() - 1, 0>())
            {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }
        
        static
        lolita::algebra::Span<DenseMatrix<Real, 3, t_element.getNumNodes(), lolita::algebra::colMajor()> const>
        getReferenceCoordinates()
        {
            using t_ReferenceCoordinates = lolita::algebra::Span<DenseMatrix<Real, 3, t_element.getNumNodes(), lolita::algebra::colMajor()> const>;
            return t_ReferenceCoordinates(ShapeTraits<t_element>::reference_nodes_.begin()->begin());
        }
        
        static
        Real
        getShapeMappingEvaluation(
            DenseVector<Real, t_element.getNumNodes()> const & nodal_field_values,
            Point const & reference_point
        )
        {
            return ElementTraits_::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }
        
        static
        Real
        getShapeMappingDerivative(
            DenseVector<Real, t_element.getNumNodes()> const & nodal_field_values,
            Point const & reference_point,
            Integer derivative_direction
        )
        {
            return ElementTraits_::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        requires(t_element.hasDim(3))
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, 3>();
            auto du = Real(0);
            ru.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            du = lolita::numerics::abs((ru.col(0).cross(ru.col(1))).dot(ru.col(2)));
            return du;
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        requires(t_element.hasDim(2))
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, 3>();
            auto du = Real(0);
            ru.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            du = lolita::numerics::abs((ru.col(0).cross(ru.col(1))).norm());
            if constexpr (MeshTraits<t_domain>::hasAxiSymmetricCoordinates())
            {
                auto r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10)
                {
                    r0 = 1.e-10;
                }
                du *= 2.0 * lolita::numerics::pi() * r0;
            }
            return du;
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        requires(t_element.hasDim(1))
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, 3>();
            auto du = Real(0);
            ru.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            du = lolita::numerics::abs(ru.col(0).norm());
            if constexpr (MeshTraits<t_domain>::hasAxiSymmetricCoordinates())
            {
                auto r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10)
                {
                    r0 = 1.e-10;
                }
                du *= 2.0 * lolita::numerics::pi() * r0;
            }
            return du;
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        requires(t_element == Node())
        {
            return 1.0;
        }
        
        // Real
        // getShapeMappingDifferential(
        //     Point const & point
        // )
        // const
        // {
        //     auto const current_coordinates = this->getCurrentCoordinates();
        //     auto ru = DenseMatrix<Real, 3, 3>();
        //     auto du = Real(0);
        //     ru.setZero();
        //     for (auto i = 0; i < t_domain.dim_; ++i)
        //     {
        //         for (auto j = 0; j < t_element.dim_; ++j)
        //         {
        //             ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
        //         }
        //     }
        //     if constexpr (t_element.dim_ == 3)
        //     {
        //         du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
        //     }
        //     else if constexpr (t_element.dim_ == 2)
        //     {
        //         du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
        //     }
        //     else
        //     {
        //         du = lolita::numerics::abs(ru.col(0).norm());
        //     }
        //     if constexpr (t_domain.frame_ == Domain::Frame::AxiSymmetric)
        //     {
        //         Real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
        //         if (r0 < 1.e-10)
        //         {
        //             r0 = 1.e-10;
        //         }
        //         du *= 2.0 * lolita::numerics::pi() * r0;
        //     }
        //     return du;
        // }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point
        )
        const
        requires(t_element.isSub(t_domain, 0))
        {
            auto const & current_coordinates = this->getCurrentCoordinates();
            auto mp0 = Point();
            auto mp1 = Point();
            mp0.setZero();
            mp1.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                mp0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                mp1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            }
            return (mp1 - mp0).norm();
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction
        )
        const
        requires(t_element.isSub(t_domain, 0))
        {
            auto const & current_coordinates = this->getCurrentCoordinates();
            auto mp0 = Point();
            auto mp1 = Point();
            mp0.setZero();
            mp1.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                mp0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                mp1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            }
            return (mp1 - mp0)(direction);
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.hasDim(2))
        {
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature::gauss(4)>::template Rule<Element::segment(1)>;
            using SegmentQuadrature = typename QuadratureTraits<GaussQuadrature(4)>::template Rule<LinearSegment{}>;
            auto distance = Real(0);
            auto dt = Real();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
            {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = DenseMatrix<Real, 3, 3>();
                auto difference = second_point - first_point;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                ru.setZero();
                for (auto i = 0; i < t_domain.dim_; ++i)
                {
                    for (auto j = 0; j < t_element.dim_; ++j)
                    {
                        auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                        auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                        ru(i, j) = dx * du;
                    }
                }
                auto Eff = ru.col(0).dot(ru.col(0));
                auto Fff = ru.col(0).dot(ru.col(1));
                auto Gff = ru.col(1).dot(ru.col(1));
                dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                distance += wq * dt;
            }
            return distance;
        }
        
        // DenseMatrix<Real, 2, 2>
        // getMetricTensor(
        //     Point const & first_point,
        //     Point const & second_point,
        //     Integer direction
        // )
        // const
        // requires(!t_element.isSub(t_domain, 0))
        // {
        //     using SegmentQuadrature = QuadratureTraits<Element::segment(1), Quadrature::gauss(4)>;
        //     auto const current_nodes_coordinates = this->getCurrentCoordinates();
        //     for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
        //     {
        //         auto pq = SegmentQuadrature::reference_points_[q][0];
        //         auto wq = SegmentQuadrature::reference_weights_[q];
        //         auto ru = DenseMatrix<Real, 3, 3>();
        //         auto difference = second_point - first_point;
        //         auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
        //         ru.setZero();
        //         // for (auto i = 0; i < t_domain.dim_; ++i)
        //         // {
        //         for (auto j = 0; j < t_element.dim_; ++j)
        //         {
        //             // if (i == direction)
        //             // {
        //             auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
        //             auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(direction), uq, j);
        //             ru(direction, j) = dx * du;
        //             // }
        //         }
        //         // }
        //         auto Eff = ru.col(0).dot(ru.col(0));
        //         auto Fff = ru.col(0).dot(ru.col(1));
        //         auto Gff = ru.col(1).dot(ru.col(1));
        //         dt = std::sqrt(Eff + 2.0 * Fff + Gff);
        //         distance += wq * dt;
        //     }
        //     return sign * distance;
        // }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.hasDim(2))
        {
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature::gauss(4)>::template Rule<Element::segment(1)>;
            using SegmentQuadrature = typename QuadratureTraits<GaussQuadrature(4)>::template Rule<LinearSegment{}>;
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature(4, Quadrature::Gauss)>::template Rule<Element::segment(1)>;
            auto distance = Real(0);
            auto dt = Real();
            // -> TEST
            auto sign = (second_point - first_point)(direction) > 0 ? 1 : -1;
            // <_ TEST
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
            {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = DenseMatrix<Real, 3, 3>();
                auto difference = second_point - first_point;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                ru.setZero();
                // for (auto i = 0; i < t_domain.dim_; ++i)
                // {
                for (auto j = 0; j < t_element.dim_; ++j)
                {
                    // if (i == direction)
                    // {
                    auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                    auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(direction), uq, j);
                    ru(direction, j) = dx * du;
                    // }
                }
                // }
                auto Eff = ru.col(0).dot(ru.col(0));
                auto Fff = ru.col(0).dot(ru.col(1));
                auto Gff = ru.col(1).dot(ru.col(1));
                dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                distance += wq * dt;
            }
            return sign * distance;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.hasDim(1))
        {
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature::gauss(4)>::template Rule<Element::segment(1)>;
            using SegmentQuadrature = typename QuadratureTraits<GaussQuadrature(4)>::template Rule<LinearSegment{}>;
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature(4, Quadrature::Gauss)>::template Rule<Element::segment(1)>;
            auto distance = Real(0);
            auto dt = Real();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
            {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = DenseMatrix<Real, 3, 3>();
                auto difference = second_point - first_point;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                ru.setZero();
                for (auto i = 0; i < t_domain.dim_; ++i)
                {
                    for (auto j = 0; j < t_element.dim_; ++j)
                    {
                        auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                        auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                        ru(i, j) = dx * du;
                    }
                }
                auto Eff = ru.col(0).dot(ru.col(0));
                dt = std::sqrt(Eff);
                distance += wq * dt;
            }
            return distance;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.hasDim(1))
        {
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature::gauss(4)>::template Rule<Element::segment(1)>;
            using SegmentQuadrature = typename QuadratureTraits<GaussQuadrature(4)>::template Rule<LinearSegment{}>;
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature(4, Quadrature::Gauss)>::template Rule<Element::segment(1)>;
            auto distance = Real(0);
            auto dt = Real();
            // -> TEST
            auto sign = (second_point - first_point)(direction) > 0 ? 1 : -1;
            // <_ TEST
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
            {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = DenseMatrix<Real, 3, 3>();
                auto difference = second_point - first_point;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                ru.setZero();
                // for (auto i = 0; i < t_domain.dim_; ++i)
                // {
                for (auto j = 0; j < t_element.dim_; ++j)
                {
                    // if (i == direction)
                    // {
                    auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                    auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(direction), uq, j);
                    ru(direction, j) = dx * du;
                    // }
                }
                // }
                auto Eff = ru.col(0).dot(ru.col(0));
                dt = std::sqrt(Eff);
                distance += wq * dt;
            }
            return sign * distance;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element == Node())
        {
            return 0.0;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element == Node())
        {
            return 0.0;
        }
        
        // Real
        // getRiemannianDistance(
        //     Point const & first_point,
        //     Point const & second_point,
        //     Integer direction = -1
        // )
        // const
        // {
        //     if constexpr (t_element.isSub(t_domain, 0))
        //     {
        //         auto const & current_coordinates = this->getCurrentCoordinates();
        //         auto distance = Real();
        //         auto mp0 = Point();
        //         auto mp1 = Point();
        //         mp0.setZero();
        //         mp1.setZero();
        //         for (auto i = 0; i < t_domain.getDim(); ++i)
        //         {
        //             mp0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
        //             mp1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
        //         }
        //         direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
        //         return distance;
        //     }
        //     else
        //     {
        //         using SegmentQuadrature = QuadratureTraits<Element::segment(1), Quadrature::gauss(4)>;
        //         auto distance = Real(0);
        //         auto dt = Real();
        //         auto const current_nodes_coordinates = this->getCurrentCoordinates();
        //         for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
        //         {
        //             auto pq = SegmentQuadrature::reference_points_[q][0];
        //             auto wq = SegmentQuadrature::reference_weights_[q];
        //             auto ru = DenseMatrix<Real, 3, 3>();
        //             auto difference = second_point - first_point;
        //             auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
        //             ru.setZero();
        //             for (auto i = 0; i < t_domain.dim_; ++i)
        //             {
        //                 for (auto j = 0; j < t_element.dim_; ++j)
        //                 {
        //                     if (direction == -1 || i == direction)
        //                     {
        //                         auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
        //                         auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
        //                         ru(i, j) = dx * du;
        //                     }
        //                 }
        //             }
        //             if constexpr (t_element.hasDim(2))
        //             {
        //                 auto Eff = ru.col(0).template dot(ru.col(0));
        //                 auto Fff = ru.col(0).template dot(ru.col(1));
        //                 auto Gff = ru.col(1).template dot(ru.col(1));
        //                 dt = std::sqrt(Eff + 2.0 * Fff + Gff);
        //             }
        //             else if constexpr (t_element.hasDim(1))
        //             {
        //                 auto Eff = ru.col(0).template dot(ru.col(0));
        //                 dt = std::sqrt(Eff);
        //             }
        //             else
        //             {
        //                 dt = 0;
        //             }
        //             distance += wq * dt;
        //         }
        //         return distance;
        //     }
        // }
        
        Real
        getLocalFrameDistance(
            Point const & first_point,
            Point const & second_point,
            Integer kkk
        )
        const
        requires(t_element.isSub(t_domain, 0))
        {
            auto first_point_mapping = Point();
            auto second_point_mapping = Point();
            auto const & current_coordinates = this->getCurrentCoordinates();
            first_point_mapping.setZero();
            second_point_mapping.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                first_point_mapping(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                second_point_mapping(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            }
            return (second_point_mapping - first_point_mapping)(kkk);
            // return (second_point - first_point)(kkk);
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
            auto rotation_matrix = getRotationMatrix(getReferenceCentroid());
            auto first_point_mapping = Point();
            auto second_point_mapping = Point();
            auto const & current_coordinates = this->getCurrentCoordinates();
            first_point_mapping.setZero();
            second_point_mapping.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                first_point_mapping(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                second_point_mapping(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            }
            return (rotation_matrix * (second_point_mapping - first_point_mapping))(kkk);
            // return (rot * (second_point - first_point))(kkk);
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
        //         mp_0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
        //         mp_1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
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
                    for (auto k = 0; k < t_element.getDim(); ++k)
                    {
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
            // auto reference_coordinates = FiniteElement::getReferenceCoordinates();
            // auto current_diameters = Point();
            // current_diameters.setZero();
            // for (auto i = 0; i < t_element.getNumNodes(); ++i)
            // {
            //     for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
            //     {
            //         auto const & pt0 = reference_coordinates.col(i);
            //         auto const & pt1 = reference_coordinates.col(j);
            //         for (auto k = 0; k < t_element.getDim(); ++k)
            //         {
            //             auto new_value = lolita::numerics::abs(getRiemannianDistance(pt0, pt1, k));
            //             auto & current_value = current_diameters(k);
            //             if (new_value > current_value)
            //             {
            //                 current_value = new_value;
            //             }
            //         }
            //     }
            // }
            // return current_diameters;
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
            auto reference_coordinates = FiniteElement::getReferenceCoordinates();
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
        
        static
        Point
        getReferenceCentroid()
        {
            auto reference_nodes_coordinates = FiniteElement::getReferenceCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i)
            {
                barycenter += reference_nodes_coordinates.col(i);
            }
            barycenter /= Real(t_element.getNumNodes());
            return barycenter;
        }
        
        Point
        getCurrentCentroid()
        const
        {
            // auto const current_nodes_coordinates = this->getCurrentCoordinates();
            // auto barycenter = Point();
            // barycenter.setZero();
            // for (auto i = 0; i < t_element.getNumNodes(); ++i) {
            //     barycenter += current_nodes_coordinates.col(i);
            // }
            // barycenter /= Real(t_element.getNumNodes());
            // return barycenter;
            auto const reference_centroid = FiniteElement::getReferenceCentroid();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_domain.getDim(); i++)
            {
                barycenter(i) = FiniteElement::getShapeMappingEvaluation(current_nodes_coordinates.row(i), reference_centroid);
            }
            return barycenter;
        }
        
        static
        Point
        getReferenceDiameters()
        {
            auto dts = Point();
            auto nds = FiniteElement::getReferenceCoordinates();
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
        
        // Point
        // getNormalVector(
        //     Point const & point
        // )
        // const
        // requires(t_element.isSub(t_domain, 1))
        // {
        //     auto const current_nodes_coordinates = this->getCurrentCoordinates();
        //     auto ru = DenseMatrix<Real, 3, t_element.getDim()>();
        //     ru.setZero();
        //     for (auto i = 0; i < t_domain.getDim(); ++i)
        //     {
        //         for (auto j = 0; j < t_element.getDim(); ++j)
        //         {
        //             ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
        //         }
        //     }
        //     if constexpr (t_element == Node()) {
        //         return Point{0, 0, 0};
        //     }
        //     else if constexpr (t_element.hasDim(1))
        //     {
        //         return Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
        //     }
        //     else
        //     {
        //         return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
        //     }
        // }
        
        Point
        getNormalVector(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(2))
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, t_element.getDim()>();
            auto normal_vector = Point();
            ru.setZero();
            normal_vector.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            normal_vector(0) = + ru(1)/ru.norm();
            normal_vector(1) = - ru(0)/ru.norm();
            return normal_vector;
        }
        
        Point
        getNormalVector(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(3))
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, t_element.getDim()>();
            auto normal_vector = Point();
            ru.setZero();
            normal_vector.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            // normal_vector = (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            normal_vector = ru.col(0).cross(ru.col(1));
            normal_vector.normalize();
            return normal_vector;
        }

        DenseMatrix<Real, 3, 3>
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
            auto rotation_matrix = DenseMatrix<Real, 3, 3>();
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

        DenseMatrix<Real, 3, 3>
        getRotationMatrix(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(2))
        {
            // setting the only basis vector needed
            auto basis_vector_0 = getNormalVector(point);
            // setting rotation matrix
            auto rotation_matrix = DenseMatrix<Real, 3, 3>();
            rotation_matrix.setZero();
            rotation_matrix(0, 0) = - basis_vector_0(1);
            rotation_matrix(0, 1) = + basis_vector_0(0);
            rotation_matrix(1, 0) = + basis_vector_0(0);
            rotation_matrix(1, 1) = + basis_vector_0(1);
            return rotation_matrix;
        }

        DenseMatrix<Real, 3, 3>
        getRotationMatrix(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(1))
        {
            // setting rotation matrix
            auto rotation_matrix = DenseMatrix<Real, 3, 3>();
            rotation_matrix.setZero();
            rotation_matrix(0, 0) = 1.0;
            return rotation_matrix;
        }
        
        Point
        getTangentVector(
            Point const & point,
            Integer direction
        )
        const
        requires(t_element.hasDim(1))
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = Point();
            tangent_vector.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
            }
            return tangent_vector;
        }
        
        template<QuadratureConcept auto t_quadrature>
        static
        Real
        getReferenceQuadratureWeight(
            Integer index
        )
        {
            return QuadratureTraits<t_quadrature>::template Rule<t_element>::reference_weights_[index];
        }
        
        template<QuadratureConcept auto t_quadrature>
        static
        lolita::algebra::Span<Point const>
        getReferenceQuadraturePoint(
            Integer index
        )
        {
            return lolita::algebra::Span<Point const>(QuadratureTraits<t_quadrature>::template Rule<t_element>::reference_points_[index].begin());
        }
        
        template<QuadratureConcept auto t_quadrature>
        Real
        getCurrentQuadratureWeight(
            Integer index
        )
        const
        {
            auto w = getReferenceQuadratureWeight<t_quadrature>(index);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<t_quadrature>(index));
        }
        
        template<QuadratureConcept auto t_quadrature>
        Point
        getCurrentQuadraturePoint(
            Integer index
        )
        const
        {
            auto p = Point();
            auto const nds = this->getCurrentCoordinates();
            p.setZero();
            for (auto j = 0; j < t_domain.getDim(); ++j)
            {
                p(j) = FiniteElement::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<t_quadrature>(index));
            }
            return p;
        }
        
        template<QuadratureConcept auto t_quadrature, Integer _i, Integer _j>
        static
        Real
        getInnerNeighborReferenceQuadratureWeight(
            Integer index
        )
        requires(t_element != Node())
        {
            auto const constexpr t_component = ElementTraits_::template getInnerNeighbor<_i, _j>();
            return FiniteElement<t_component, t_domain>::template getReferenceQuadratureWeight<t_quadrature>(index);
        }
        
        template<QuadratureConcept auto t_quadrature, Integer _i, Integer _j>
        static
        Point
        getInnerNeighborReferenceQuadraturePoint(
            Integer component_index,
            Integer index
        )
        requires(t_element != Node())
        {
            auto constexpr t_component = ElementTraits_ ::template getInnerNeighbor<_i, _j>();
            auto p = Point();
            p.setZero();
            auto const & elt_reference_nodes = ShapeTraits<t_element>::reference_nodes_;
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                auto cpt_coordinates = DenseVector<Real, t_component.getNumNodes()>();
                for (auto j = 0; j < t_component.getNumNodes(); ++j)
                {
                    auto const node_tag = getInnerNeighborNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = FiniteElement<t_component, t_domain>::template getReferenceQuadraturePoint<t_quadrature>(index);
                p(i) = FiniteElement<t_component, t_domain>::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }
        
        template<QuadratureConcept auto t_quadrature, Integer t_i, Integer t_j>
        Real
        getInnerNeighborCurrentQuadratureWeight(
            Integer component_index,
            Integer index
        )
        const
        requires(t_element != Node())
        {
            auto const & cmp =  this->template getInnerNeighbors<t_i, t_j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<t_quadrature>(index);
        }
        
        template<QuadratureConcept auto t_quadrature, Integer t_i, Integer t_j>
        Point
        getInnerNeighborCurrentQuadraturePoint(
            Integer component_index,
            Integer index
        )
        const
        requires(t_element != Node())
        {
            auto p = Point();
            auto const cpt_ref_pnt = getInnerNeighborReferenceQuadraturePoint<t_quadrature, t_i, t_j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            p.setZero();
            for (auto j = 0; j < t_domain.getDim(); ++j)
            {
                p(j) = FiniteElement::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

    };

} // namespace lolita::core


#endif /* B44BF9EE_1549_420E_9AF4_5820C5464430 */
