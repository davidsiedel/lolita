#ifndef F0FF180A_BC96_4D84_87EF_E2EAB877F599
#define F0FF180A_BC96_4D84_87EF_E2EAB877F599

#include "2/core/_include.hxx"

namespace lolita::core
{
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    struct FiniteElement;
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain, FieldConcept auto t_field>
    struct DiscretizationImplementation;

    struct ExternalLoad
    {

        ExternalLoad(
            Integer row,
            Integer col,
            std::function<Real(Point const &, Real const &)> const & function
        )
        :
        row_(row),
        col_(col),
        function_(function)
        {}

        ExternalLoad(
            Integer row,
            Integer col,
            std::function<Real(Point const &, Real const &)> && function
        )
        :
        row_(row),
        col_(col),
        function_(std::move(function))
        {}

        inline
        Integer
        getRow()
        const
        {
            return row_;
        }

        inline
        Integer
        getCol()
        const
        {
            return col_;
        }

        inline
        Real
        getValue(
            Point const & point,
            Real const & time
        )
        const
        {
            return function_(point, time);
        }

    private:

        Integer row_;

        Integer col_;

        std::function<Real(Point const &, Real const &)> function_;

    };

    template<Integer num_coefficients_>
    struct UknImplementation;
    
    struct UknInterface
    {

    protected:

        UknInterface()
        {}

    public:

        template<Integer num_coefficients_>
        DenseVector<Real, num_coefficients_> const &
        getCoefficients()
        const
        {
            return static_cast<UknImplementation<num_coefficients_> const *>(this)->getCoefficients();
        }

    };

    template<Integer num_coefficients_>
    struct UknImplementation : UknInterface
    {

        static constexpr
        Integer
        getNumCoefficients()
        {
            return num_coefficients_;
        }

        UknImplementation()
        :
        current_coefficients_(DenseVector<Real, num_coefficients_>::Zero()),
        previous_coefficients_(DenseVector<Real, num_coefficients_>::Zero())
        {}

        DenseVector<Real, num_coefficients_> const &
        getCoefficients()
        const
        {
            return current_coefficients_;
        }
        
    private:

        DenseVector<Real, num_coefficients_> current_coefficients_;

        DenseVector<Real, num_coefficients_> previous_coefficients_;
        
    };

    struct DegreesOfFreedomCoordinates
    {
        
        constexpr
        DegreesOfFreedomCoordinates(
            Natural const & rhs_offset,
            Natural const & lhs_offset
        )
        :
        rhs_offset_(rhs_offset),
        lhs_offset_(lhs_offset)
        {}

        constexpr
        DegreesOfFreedomCoordinates(
            Natural && rhs_offset,
            Natural && lhs_offset
        )
        :
        rhs_offset_(std::move(rhs_offset)),
        lhs_offset_(std::move(lhs_offset))
        {}
    
        constexpr
        Natural const &
        getRhsOffset()
        const
        {
            return rhs_offset_;
        }
    
        constexpr
        Natural const &
        getLhsOffset()
        const
        {
            return lhs_offset_;
        }

        constexpr
        void
        setOffsets(
            Natural && rhs_offset,
            Natural && lhs_offset
        )
        {
            rhs_offset_ = std::forward<Natural>(rhs_offset);
            lhs_offset_ = std::forward<Natural>(lhs_offset);
        }

    private:

        Natural rhs_offset_;

        Natural lhs_offset_;

    };

    template<LagrangeShapeConcept auto t_element, MeshConcept auto t_domain, FieldConcept auto t_field>
    struct DegreesOfFreedomImplementation;

    template<LagrangeShapeConcept auto t_element, MeshConcept auto t_domain>
    struct DegreesOfFreedomInterface
    {

        template<FieldConcept auto t_field>
        static constexpr
        Integer
        getScalarSize()
        {
            return DegreesOfFreedomImplementation<t_element, t_domain, t_field>::getScalarSize();
        }

        template<FieldConcept auto t_field>
        static constexpr
        Integer
        getTensorSize()
        {
            return DegreesOfFreedomImplementation<t_element, t_domain, t_field>::getTensorSize();
        }

    private:

        template<FieldConcept auto t_field>
        using Implementation_ = DegreesOfFreedomImplementation<t_element, t_domain, t_field>;

        using FiniteElement_ = FiniteElement<t_element, t_domain>;
        
        template<FieldConcept auto t_field>
        using DiscretizationImplementation_ = DiscretizationImplementation<t_element, t_domain, t_field>;

    public:

        DegreesOfFreedomInterface(
            FieldConcept auto const & field,
            FiniteElement_ const & finite_element
        )
        :
        label_(field.getLabel()),
        finite_element_(finite_element),
        coordinates_()
        {}
        
        Boolean
        operator==(
            DegreesOfFreedomInterface const & other
        )
        const
        = default;
        
        Boolean
        operator!=(
            DegreesOfFreedomInterface const & other
        )
        const
        = default;

        Label const &
        getLabel()
        const
        {
            return label_;
        }
        
        Boolean
        isLinked()
        const
        {
            return coordinates_ != nullptr;
        }

        // virtual
        // Integer
        // getScalarSize2()
        // const
        // =0;

        // virtual
        // Integer
        // getTensorSize2()
        // const
        // =0;

        void
        link(
            Natural && rhs_offset,
            Natural && lhs_offset
        )
        {
            coordinates_ = std::make_unique<DegreesOfFreedomCoordinates>(std::forward<Natural>(rhs_offset), std::forward<Natural>(lhs_offset));
        }

        template<FieldConcept auto t_field>
        auto
        getCoefficients()
        const
        {
            return static_cast<Implementation_<t_field> const *>(this)->getCoefficients();
        }

        template<FieldConcept auto t_field>
        auto
        getCoefficients(
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<Implementation_<t_field> const *>(this)->getCoefficients(row, col);
        }

        template<FieldConcept auto t_field>
        void
        upgradeCoefficients(
            DenseVectorConcept<Real> auto && input
        )
        {
            static_cast<Implementation_<t_field> *>(this)->upgradeCoefficients(std::forward<decltype(input)>(input));
        }

        template<FieldConcept auto t_field>
        void
        reserveCoefficients()
        {
            static_cast<Implementation_<t_field> *>(this)->reserveCoefficients();
        }

        template<FieldConcept auto t_field>
        void
        recoverCoefficients()
        {
            static_cast<Implementation_<t_field> *>(this)->recoverCoefficients();
        }

        template<FieldConcept auto t_field, auto t_s>
        void
        addToLinearSystem(
            LinearSystem<t_s> & linear_system
        )
        {
            coordinates_ = std::make_unique<DegreesOfFreedomCoordinates>(linear_system.getRhsSize(), linear_system.getLhsSize());
            // this->link(linear_system.getRhsSize(), linear_system.getLhsSize());
            linear_system.addRhsSize(DiscretizationTraits<t_field>::template getTensorBasisSize<t_element, t_domain>());
            auto ci = this->template getBandWidth<t_field>();
            std::cout << "field size : " << FieldTraits<t_field>::template getSize<t_domain>() << std::endl;
            std::cout << this->getFiniteElement().getHash() << " : " << ci << std::endl;
            linear_system.addLhsSize(this->template getBandWidth<t_field>());
        }

        template<FieldConcept auto t_field>
        Real
        getDiscreteFieldDegreeOfFreedomValue(
            PointConcept auto const & point,
            Integer row,
            Integer col
        )
        const
        {
            auto const coefficients = this->template getUnknownCoefficients<t_field>();
            auto const basis_vector = this->template getFieldDualVector<t_field>(point, row, col);
            return coefficients.dot(basis_vector);
        }

        template<FieldConcept auto t_field>
        DenseVector<Real, DiscretizationTraits<t_field>::template getTensorSpaceSize<t_element, t_domain>()>
        getUnknownCoefficients()
        const
        {
            auto constexpr cell_range = DiscretizationTraits<t_field>::template getTensorBasisSize<t_element, t_domain>();
            auto offset = 0;
            auto unknown = DenseVector<Real, DiscretizationTraits<t_field>::template getTensorSpaceSize<t_element, t_domain>()>();
            auto cell_block = unknown.template segment<cell_range>(offset);
            cell_block = this->template getCoefficients<t_field>();
            offset += cell_range;
            auto set_faces_unknowns = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_inner_neighbor = ShapeTraits<t_element>::template getInnerNeighbor<t_i, t_j>();
                auto constexpr range = DiscretizationTraits<t_field>::template getTensorBasisSize<t_inner_neighbor, t_domain>();
                if constexpr (range > 0)
                {
                    for (auto const & face : this->getFiniteElement().template getInnerNeighbors<t_i, t_j>())
                    {
                        auto face_block = unknown.template segment<range>(offset);
                        face_block = face->template getDiscreteField<t_field>().template getCoefficients<t_field>();
                        offset += range;
                    }
                }
                if constexpr (t_j < ShapeTraits<t_element>::template getNumInnerNeighbors<t_i>() - 1)
                {
                    self.template operator ()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ShapeTraits<t_element>::template getNumInnerNeighbors<>() - 1)
                {
                    self.template operator ()<t_i + 1, 0>(self);
                }
            };
            set_faces_unknowns(set_faces_unknowns);
            return unknown;
        }
        
        template<FieldConcept auto t_field, Integer t_a, Integer t_b>
        static constexpr
        Integer
        getNeighborOffset(
            Integer kkk
        )
        {
            auto num_unknowns = DiscretizationTraits<t_field>::template getTensorBasisSize<t_element, t_domain>();
            auto set_num_unknowns = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_num_unknowns
            )
            constexpr mutable
            {
                auto constexpr inner_neighbor = ShapeTraits<t_element>::template getInnerNeighbor<t_i, t_j>();
                auto constexpr num_inner_neighbors = ShapeTraits<t_element>::template getNumInnerNeighbors<t_i, t_j>();
                for (auto iii = 0; iii < num_inner_neighbors; iii++)
                {
                    if (kkk == iii && t_a == t_i && t_b == t_j)
                    {
                        return num_unknowns;
                    }
                    else
                    {
                        num_unknowns += DiscretizationTraits<t_field>::template getTensorBasisSize<inner_neighbor, t_domain>();
                    }
                }
                if constexpr (t_j < ShapeTraits<t_element>::template getNumInnerNeighbors<t_i>() - 1)
                {
                    t_set_num_unknowns.template operator ()<t_i, t_j + 1>(t_set_num_unknowns);
                }
                else if constexpr (t_i < ShapeTraits<t_element>::template getNumInnerNeighbors<>() - 1)
                {
                    t_set_num_unknowns.template operator ()<t_i + 1, 0>(t_set_num_unknowns);
                }
            };
            set_num_unknowns(set_num_unknowns);
            return num_unknowns;
        }
        
        template<FieldConcept auto t_field, Integer t_a, Integer t_b>
        Integer
        getNeighborOffset(
            std::shared_ptr<FiniteElement<ShapeTraits<t_element>::template getInnerNeighbor<t_a, t_b>(), t_domain>> const & inner_neighbor
        )
        {
            auto num_unknowns = DiscretizationTraits<t_field>::template getTensorBasisSize<t_element, t_domain>();
            auto set_num_unknowns = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_num_unknowns
            )
            constexpr mutable
            {
                auto constexpr inner_neighbor = ShapeTraits<t_element>::template getInnerNeighbor<t_i, t_j>();
                auto constexpr num_inner_neighbors = ShapeTraits<t_element>::template getNumInnerNeighbors<t_i, t_j>();
                for (auto const & in : this->getFiniteElement().template getInnerNeighbors<t_i, t_j>())
                {
                    if (utility::areEqual(inner_neighbor, in))
                    {
                        return num_unknowns;
                    }
                    else
                    {
                        num_unknowns += DiscretizationTraits<t_field>::template getTensorBasisSize<inner_neighbor, t_domain>();
                    }
                }
                if constexpr (t_j < ShapeTraits<t_element>::template getNumInnerNeighbors<t_i>() - 1)
                {
                    t_set_num_unknowns.template operator ()<t_i, t_j + 1>(t_set_num_unknowns);
                }
                else if constexpr (t_i < ShapeTraits<t_element>::template getNumInnerNeighbors<>() - 1)
                {
                    t_set_num_unknowns.template operator ()<t_i + 1, 0>(t_set_num_unknowns);
                }
            };
            set_num_unknowns(set_num_unknowns);
            return num_unknowns;
        }
        
        template<FieldConcept auto t_field>
        Integer
        getBandWidth()
        const
        {
            auto band_width = DiscretizationTraits<t_field>::template getTensorBasisSize<t_element, t_domain>();
            auto constexpr cell_offset = t_field.getDimDomain() - t_element.getDim();
            auto set_faces_unknowns2 = [&] <Integer t_i = 0, Integer t_j = 0, Integer t_k = 0> (
                auto & self2
            )
            constexpr mutable
            {
                auto constexpr t_cell = ShapeTraits<t_element>::template getOuterNeighbor<t_domain, cell_offset, t_i>();
                auto constexpr t_cell_neighbor = ShapeTraits<t_cell>::template getInnerNeighbor<t_j, t_k>();
                for (auto const & cell : this->getFiniteElement().template getOuterNeighbors<cell_offset, t_i>())
                {
                    for (auto const & cell_neighbor : cell->template getInnerNeighbors<t_j, t_k>())
                    {
                        if (!utility::areEqual(* cell_neighbor, this->getFiniteElement()))
                        {
                            // std::cout << t_cell_neighbor << std::endl;
                            // band_width += DiscretizationTraits<t_field>::template getTensorSpaceSize<t_cell_neighbor, t_domain>();
                            band_width += DiscretizationTraits<t_field>::template getTensorBasisSize<t_cell_neighbor, t_domain>();
                        }
                        else
                        {
                            // std::cout << "HI : " << this->getFiniteElement().getHash() << std::endl;
                        }
                    }
                }
                if constexpr (t_k < ShapeTraits<t_cell>::template getNumInnerNeighbors<t_j>() - 1)
                {
                    self2.template operator ()<t_i, t_j, t_k + 1>(self2);
                }
                else if constexpr (t_j < ShapeTraits<t_cell>::template getNumInnerNeighbors<>() - 1)
                {
                    self2.template operator ()<t_i, t_j + 1, 0>(self2);
                }
                else if constexpr (t_i < ShapeTraits<t_element>::template getNumOuterNeighbors<t_domain, cell_offset>() - 1)
                {
                    self2.template operator ()<t_i + 1, 0, 0>(self2);
                }
            };
            if constexpr (t_element.getDim() < t_field.getDimDomain())
            {
                set_faces_unknowns2(set_faces_unknowns2);
            }
            return band_width;
        }

        template<FieldConcept auto t_field>
        DenseVector<Real, DiscretizationTraits<t_field>::template getTensorSpaceSize<t_element, t_domain>()>
        getFieldDualVector(
            PointConcept auto const & point
        )
        const
        {
            return static_cast<DiscretizationImplementation_<t_field> const *>(this)->getFieldDualVector(point);
        };

        template<FieldConcept auto t_field>
        DenseVector<Real, DiscretizationTraits<t_field>::template getTensorSpaceSize<t_element, t_domain>()>
        getFieldDualVector(
            PointConcept auto const & point,
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<DiscretizationImplementation_<t_field> const *>(this)->getFieldDualVector(point, row, col);
        };

        template<FieldConcept auto t_field>
        DenseVector<Real, DiscretizationTraits<t_field>::template getTensorSpaceSize<t_element, t_domain>()>
        getFieldPrimalVector(
            PointConcept auto const & point
        )
        const
        {
            return static_cast<DiscretizationImplementation_<t_field> const *>(this)->getFieldPrimalVector(point);
        };

        template<FieldConcept auto t_field>
        DenseVector<Real, DiscretizationTraits<t_field>::template getTensorSpaceSize<t_element, t_domain>()>
        getFieldPrimalVector(
            PointConcept auto const & point,
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<DiscretizationImplementation_<t_field> const *>(this)->getFieldPrimalVector(point, row, col);
        };

        template<FieldConcept auto t_field, LinearOperatorConcept auto t_mapping>
        auto
        letLinearOperator(
            Point const & point
        )
        const
        {
            return static_cast<DiscretizationImplementation_<t_field> const *>(this)->template letLinearOperator<t_mapping>(point);
        }

        template<LinearOperatorConcept auto t_mapping>
        auto
        letLinearOperator(
            Point const & point
        )
        const
        {
            return static_cast<DiscretizationImplementation_<t_mapping.getField()> const *>(this)->template letLinearOperator<t_mapping>(point);
        }

        FiniteElement_ const &
        getFiniteElement()
        const
        {
            return finite_element_;
        }

    protected:

        Label const & label_;

        FiniteElement_ const & finite_element_;

        std::unique_ptr<DegreesOfFreedomCoordinates> coordinates_;

    };

    template<LagrangeShapeConcept auto t_element, MeshConcept auto t_domain, FieldConcept auto t_field>
    struct DegreesOfFreedomImplementation : DegreesOfFreedomInterface<t_element, t_domain>
    {

    private:

        static constexpr
        Integer
        getNumRows()
        {
            return FieldTraits<t_field>::template getNumRows<t_domain>();
        }

        static constexpr
        Integer
        getNumCols()
        {
            return FieldTraits<t_field>::template getNumCols<t_domain>();
        }

        static constexpr
        Integer
        getScalarSize()
        {
            return DiscretizationTraits<t_field>::template getScalarBasisSize<t_element, t_domain>();
        }

        static constexpr
        Integer
        getTensorSize()
        {
            return DiscretizationTraits<t_field>::template getTensorBasisSize<t_element, t_domain>();
        }

        using Base_ = DegreesOfFreedomInterface<t_element, t_domain>;

        using FiniteElement_ = FiniteElement<t_element, t_domain>;
        
        using DiscretizationImplementation_ = DiscretizationImplementation<t_element, t_domain, t_field>;

    public:

        explicit
        DegreesOfFreedomImplementation(
            FiniteElement_ const & finite_element
        )
        :
        Base_(t_field, finite_element),
        s0_(),
        s1_()
        {
            algebra::View<DenseVector<Real, getTensorSize()>>(s1_[0][0].data()) = DenseVector<Real, getTensorSize()>::Zero();
            algebra::View<DenseVector<Real, getTensorSize()>>(s0_[0][0].data()) = DenseVector<Real, getTensorSize()>::Zero();
            std::cout << algebra::View<DenseVector<Real, getTensorSize()> const>(s1_[0][0].data()) << std::endl;
            std::cout << "---" << std::endl;
        }

        algebra::View<DenseVector<Real, getScalarSize()> const>
        getCoefficients(
            Integer row,
            Integer col
        )
        const
        {
            return algebra::View<DenseVector<Real, getScalarSize()> const>(s1_[row][col].data());
        }

        algebra::View<DenseVector<Real, getTensorSize()> const>
        getCoefficients()
        const
        {
            return algebra::View<DenseVector<Real, getTensorSize()> const>(s1_[0][0].data());
        }

        void
        upgradeCoefficients(
            DenseVectorConcept<Real> auto && input
        )
        {
            algebra::View<DenseVector<Real, getTensorSize()>>(s1_[0][0].data()) += std::forward<decltype(input)>(input);
        }

        void
        reserveCoefficients()
        {
            s0_ = s1_;
        }

        void
        recoverCoefficients()
        {
            s1_ = s0_;
        }

    private:
        
        std::array<std::array<DenseVector<Real, getScalarSize()>, getNumRows()>, getNumCols()> s1_;
        
        std::array<std::array<DenseVector<Real, getScalarSize()>, getNumRows()>, getNumCols()> s0_;

        std::unique_ptr<StaticCondensationInterface> condensation_;

    };

    template<MeshConcept auto t_domain, FieldConcept auto t_field>
    struct MeshFieldImplementation;
    
    template<MeshConcept auto t_domain>
    struct MeshDiscreteField
    {

        explicit
        MeshDiscreteField(
            FieldConcept auto const & field
        )
        :
        label_(field.getLabel())
        {}
        
        Boolean
        operator==(
            MeshDiscreteField const & other
        )
        const
        = default;
        
        Boolean
        operator!=(
            MeshDiscreteField const & other
        )
        const
        = default;

        Label const &
        getLabel()
        const
        {
            return label_;
        }

        void
        addLoad(
            Integer row,
            Integer col,
            std::function<Real(Point const &, Real const &)> && function
        )
        {
            if (loads_ == nullptr)
            {
                loads_ = std::make_unique<std::vector<ExternalLoad>>();
            }
            for (auto & load : * loads_)
            {
                if (load.getRow() == row && load.getCol() == col)
                {
                    load = ExternalLoad(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function));
                    return;
                }
            }
            loads_->push_back(ExternalLoad(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function)));
        }

        Boolean
        hasLoads()
        const
        {
            return loads_ != nullptr;
        }

        std::vector<ExternalLoad> const &
        getLoads()
        const
        {
            return * loads_;
        }

    private:

        Label const & label_;

        std::unique_ptr<std::vector<ExternalLoad>> loads_;

    };
    
} // namespace lolita::core

#endif /* F0FF180A_BC96_4D84_87EF_E2EAB877F599 */
