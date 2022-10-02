#ifndef A28EDC0F_E4ED_44BB_A20D_555BFCABEC71
#define A28EDC0F_E4ED_44BB_A20D_555BFCABEC71

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/100_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/200_quadrature.hxx"

namespace lolita
{

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

    struct DofCoordinates
    {
        
        DofCoordinates(
            Natural const & rhs_offset,
            Natural const & lhs_offset
        )
        :
        rhs_offset_(rhs_offset),
        lhs_offset_(lhs_offset)
        {}

        DofCoordinates(
            Natural && rhs_offset,
            Natural && lhs_offset
        )
        :
        rhs_offset_(std::move(rhs_offset)),
        lhs_offset_(std::move(lhs_offset))
        {}
    
        inline
        Natural const &
        getRhsOffset()
        const
        {
            return rhs_offset_;
        }
    
        inline
        Natural const &
        getLhsOffset()
        const
        {
            return lhs_offset_;
        }

        inline
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

    template<Domain t_domain>
    struct DegreeOfFreedom
    {

        // template<Integer t_size>
        // static constexpr
        // Integer
        // getSize()
        // {
        //     return t_size;
        // }

        // template<FieldConcept auto t_field, Integer t_size>
        // static constexpr
        // Integer
        // getSize()
        // {
        //     return FieldTraits<t_field>::template getSize<t_domain>() * t_size;
        // }

        // template<Element t_element, Basis t_basis>
        // static constexpr
        // Integer
        // getSize()
        // {
        //     return BasisTraits<t_basis>::template getSize<t_element>();
        // }

        // template<Element t_element, FieldConcept auto t_field, Basis t_basis>
        // static constexpr
        // Integer
        // getSize()
        // {
        //     return FieldTraits<t_field>::template getSize<t_domain>() * BasisTraits<t_basis>::template getSize<t_element>();
        // }

        // template<FieldConcept auto t_field, Integer t_size, Strategy t_s>
        // static
        // DegreeOfFreedom
        // make(
        //     std::unique_ptr<LinearSystem<t_s>> const & linear_system
        // )
        // {
        //     auto element_degree_of_freedom = DegreeOfFreedom(getSize<t_field, t_size>(), linear_system->getSize());
        //     linear_system->getSize() += getSize<t_field, t_size>();
        //     return element_degree_of_freedom;
        // }

        // template<FieldConcept auto t_field, Integer t_size>
        // static
        // DegreeOfFreedom
        // make()
        // {
        //     auto element_degree_of_freedom = DegreeOfFreedom(t_size);
        //     return element_degree_of_freedom;
        // }

        // template<Element t_element, FieldConcept auto t_field, Basis t_basis, Strategy t_s>
        // static
        // DegreeOfFreedom
        // make(
        //     std::unique_ptr<LinearSystem<t_s>> const & linear_system
        // )
        // {
        //     auto element_degree_of_freedom = DegreeOfFreedom(getSize<t_element, t_field, t_basis>(), linear_system->getSize());
        //     linear_system->getSize() += getSize<t_element, t_field, t_basis>();
        //     return element_degree_of_freedom;
        // }

        // template<Element t_element, FieldConcept auto t_field, Basis t_basis>
        // static
        // DegreeOfFreedom
        // make()
        // {
        //     auto element_degree_of_freedom = DegreeOfFreedom(getSize<t_element, t_field, t_basis>());
        //     return element_degree_of_freedom;
        // }

    // private:

        explicit
        DegreeOfFreedom(
            Integer size
        )
        :
        coordinates_(),
        s0(DenseVector<Real>::Zero(size)),
        s1(DenseVector<Real>::Zero(size))
        {}
        
        // DegreeOfFreedom(
        //     Integer size,
        //     Integer offset
        // )
        // :
        // offset_(offset),
        // s0(DenseVector<Real>::Zero(size)),
        // s1(DenseVector<Real>::Zero(size))
        // {}

    // public:
        
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
        
        Boolean
        isLinked()
        const
        {
            return coordinates_ != nullptr;
        }

        void
        link(
            Natural && rhs_offset,
            Natural && lhs_offset
        )
        {
            coordinates_ = std::make_unique<DofCoordinates>(std::forward<Natural>(rhs_offset), std::forward<Natural>(lhs_offset));
        }
        
        // Natural &
        // getOffset()
        // {
        //     return rhs_offset_;
        // }
        
        // Natural const &
        // getOffset()
        // const
        // {
        //     return rhs_offset_;
        // }
        
        // Natural &
        // getLhsOffset()
        // {
        //     return lhs_offset_;
        // }
        
        // Natural const &
        // getLhsOffset()
        // const
        // {
        //     return lhs_offset_;
        // }

        Integer
        getSize()
        const
        {
            return s1.size();
        }
        
        void
        addCoefficients(
            DenseVectorConcept<Real> auto && input
        )
        {
            s1 += std::forward<decltype(input)>(input);
        }

        // template<FieldConcept auto t_field, Integer t_size>
        // void
        // addCoefficients(
        //     DenseVectorConcept<Real> auto && input
        // )
        // {
        //     auto const & offset = coordinates_->getRhsOffset();
        //     s1 += std::forward<decltype(input)>(input).template segment<DegreeOfFreedom::template getSize<t_field, t_size>()>(offset);
        // }

        // template<Element t_element, FieldConcept auto t_field, Basis t_basis>
        // void
        // addCoefficients(
        //     DenseVectorConcept<Real> auto && input
        // )
        // {
        //     auto const & offset = coordinates_->getRhsOffset();
        //     s1 += std::forward<decltype(input)>(input).template segment<DegreeOfFreedom::template getSize<t_element, t_field, t_basis>()>(offset);
        // }

        // template<Element t_element, FieldConcept auto t_field, Basis t_basis>
        // void
        // addCoefficients(
        //     DenseVectorConcept<Real> auto && input,
        //     Integer offset
        // )
        // {
        //     s1 += std::forward<decltype(input)>(input).template segment<DegreeOfFreedom::template getSize<t_element, t_field, t_basis>()>(offset);
        // }

        // template<FieldConcept auto t_field, Integer t_size>
        // algebra::View<DenseVector<Real, DegreeOfFreedom::template getSize<t_size>()> const>
        // getCoefficients(
        //     Integer row,
        //     Integer col
        // )
        // const
        // {
        //     auto const & data = s1.data() + FieldTraits<t_field>::template getCols<t_domain>() * row + col;
        //     return algebra::View<DenseVector<Real, DegreeOfFreedom::template getSize<t_size>()> const>(data);
        // }

        // template<FieldConcept auto t_field, Integer t_size>
        // algebra::View<DenseVector<Real, DegreeOfFreedom::template getSize<t_field, t_size>()> const>
        // getCoefficients()
        // const
        // {
        //     auto const & data = s1.data();
        //     return algebra::View<DenseVector<Real, DegreeOfFreedom::template getSize<t_field, t_size>()> const>(data);
        // }

        template<Integer t_size>
        algebra::View<DenseVector<Real, t_size> const>
        getCoefficients(
            Integer offset
        )
        const
        {
            auto const & data = s1.data() + offset;
            return algebra::View<DenseVector<Real, t_size> const>(data);
        }

        template<FieldConcept auto t_field, Integer t_size>
        algebra::View<DenseVector<Real, t_size> const>
        getCoefficients(
            Integer row,
            Integer col
        )
        const
        {
            auto const & data = s1.data() + FieldTraits<t_field>::template getCols<t_domain>() * row + col;
            return algebra::View<DenseVector<Real, t_size> const>(data);
        }

        template<Integer t_size>
        algebra::View<DenseVector<Real, t_size> const>
        getCoefficients()
        const
        {
            auto const & data = s1.data();
            return algebra::View<DenseVector<Real, t_size> const>(data);
        }

        template<Integer t_size1, Integer t_size>
        auto
        packUp(
            DenseSystem<Real, t_size> const & system
        )
        {
            condensation_ = StaticCondensation::template make<Real, t_size, t_size1>(system);
            auto j = condensation_->template packUp<Real, t_size, t_size1>(system);
            return j;
        }

        template<Integer t_size1, Integer t_size>
        auto
        unPack(
            DenseVectorConcept<Real, t_size> auto const & vector
        )
        {
            auto j = condensation_->template unPack<Real, t_size, t_size1>(vector);
            condensation_ = nullptr;
            return j;
        }

        void
        reserve()
        {
            s0 = s1;
        }

        void
        recover()
        {
            s1 = s0;
        }

        Boolean
        isCondensed()
        const
        {
            return condensation_ != nullptr;
        }

    private:

        std::unique_ptr<StaticCondensation> condensation_;

        std::unique_ptr<DofCoordinates> coordinates_;

        DenseVector<Real> s0;

        DenseVector<Real> s1;

    };

    template<typename t_Key, typename t_Value>
    struct ElementaryOperator
    {

        explicit
        ElementaryOperator(
            t_Key const & tag
        )
        :
        tag_(tag),
        operator_()
        {}
        
        ElementaryOperator(
            t_Key const & tag,
            auto const & value
        )
        :
        tag_(tag),
        operator_(value)
        {}
        
        ElementaryOperator(
            t_Key const & tag,
            auto && value
        )
        :
        tag_(tag),
        operator_(std::move(value))
        {}

        t_Key const &
        getTag()
        const
        {
            return tag_;
        }

        t_Value const &
        getOperator()
        const
        {
            return operator_;
        }

        t_Value &
        getOperator()
        {
            return operator_;
        }

    private:

        t_Key const & tag_;

        t_Value operator_;

    };
    
    template<Domain t_domain>
    struct DiscreteFieldBase
    {

        explicit
        DiscreteFieldBase(
            FieldConcept auto const & field
        )
        :
        label_(field.getLabel())
        {}
        
        Boolean
        operator==(
            DiscreteFieldBase const & other
        )
        const = default;
        
        Boolean
        operator!=(
            DiscreteFieldBase const & other
        )
        const = default;

        Label const &
        getLabel()
        const
        {
            return label_;
        }

        void
        addScalar(
            Label const & tag,
            auto &&... scalar
        )
        {
            if (scalar_items_ == nullptr)
            {
                scalar_items_ = std::make_unique<std::vector<ElementaryOperator<Label, Real>>>();
            }
            for (auto & m : * scalar_items_)
            {
                if (m.getTag() == tag)
                {
                    return;
                }
            }
            scalar_items_->push_back(ElementaryOperator<Label, Real>(tag, std::forward<decltype(scalar)>(scalar)...));
        }

        void
        addVector(
            Label const & tag,
            DenseVectorConcept<Real> auto &&... vector
        )
        {
            if (vector_items_ == nullptr)
            {
                vector_items_ = std::make_unique<std::vector<ElementaryOperator<Label, DenseVector<Real>>>>();
            }
            for (auto & m : * vector_items_)
            {
                if (m.getTag() == tag)
                {
                    return;
                }
            }
            vector_items_->push_back(ElementaryOperator<Label, DenseVector<Real>>(tag, std::forward<decltype(vector)>(vector)...));
        }

        void
        addMatrix(
            Label const & tag,
            DenseMatrixConcept<Real> auto &&... matrix
        )
        {
            if (matrix_items_ == nullptr)
            {
                matrix_items_ = std::make_unique<std::vector<ElementaryOperator<Label, DenseMatrix<Real>>>>();
            }
            for (auto & m : * matrix_items_)
            {
                if (m.getTag() == tag)
                {
                    return;
                }
            }
            matrix_items_->push_back(ElementaryOperator<Label, DenseMatrix<Real>>(tag, std::forward<decltype(matrix)>(matrix)...));
        }

        Real const &
        getScalar(
            std::basic_string_view<Character> && tag
        )
        const
        {
            if (scalar_items_ == nullptr)
            {
                throw std::runtime_error("No scalar field data");
            }
            else
            {
                for (auto const & m : * scalar_items_)
                {
                    if (m.getTag() == std::forward<std::basic_string_view<Character>>(tag))
                    {
                        return m.getOperator();
                    }
                }
                throw std::runtime_error("No scalar with this tag");
            }
        }

        Real &
        getScalar(
            std::basic_string_view<Character> && tag
        )
        {
            if (scalar_items_ == nullptr)
            {
                throw std::runtime_error("No scalar field data");
            }
            else
            {
                for (auto & m : * scalar_items_)
                {
                    if (m.getTag() == std::forward<std::basic_string_view<Character>>(tag))
                    {
                        return m.getOperator();
                    }
                }
                throw std::runtime_error("No scalar with this tag");
            }
        }

        DenseVector<Real> const &
        getVector(
            std::basic_string_view<Character> && tag
        )
        const
        {
            if (vector_items_ == nullptr)
            {
                throw std::runtime_error("No vector field data");
            }
            else
            {
                for (auto const & m : * vector_items_)
                {
                    if (m.getTag() == std::forward<std::basic_string_view<Character>>(tag))
                    {
                        return m.getOperator();
                    }
                }
                throw std::runtime_error("No vector with this tag");
            }
        }

        DenseVector<Real> &
        getVector(
            std::basic_string_view<Character> && tag
        )
        {
            if (vector_items_ == nullptr)
            {
                throw std::runtime_error("No vector field data");
            }
            else
            {
                for (auto & m : * vector_items_)
                {
                    if (m.getTag() == std::forward<std::basic_string_view<Character>>(tag))
                    {
                        return m.getOperator();
                    }
                }
                throw std::runtime_error("No vector with this tag");
            }
        }

        DenseMatrix<Real> const &
        getMatrix(
            std::basic_string_view<Character> && tag
        )
        const
        {
            if (matrix_items_ == nullptr)
            {
                throw std::runtime_error("No matrix field data");
            }
            else
            {
                for (auto const & m : * matrix_items_)
                {
                    if (m.getTag() == std::forward<std::basic_string_view<Character>>(tag))
                    {
                        return m.getOperator();
                    }
                }
                throw std::runtime_error("No matrix with this tag");
            }
        }

        DenseMatrix<Real> &
        getMatrix(
            std::basic_string_view<Character> && tag
        )
        {
            if (matrix_items_ == nullptr)
            {
                throw std::runtime_error("No matrix field data");
            }
            else
            {
                for (auto & m : * matrix_items_)
                {
                    if (m.getTag() == std::forward<std::basic_string_view<Character>>(tag))
                    {
                        return m.getOperator();
                    }
                }
                throw std::runtime_error("No matrix with this tag");
            }
        }

    protected:

        Label const & label_;

        std::unique_ptr<std::vector<ElementaryOperator<Label, Real>>> scalar_items_;

        std::unique_ptr<std::vector<ElementaryOperator<Label, DenseVector<Real>>>> vector_items_;

        std::unique_ptr<std::vector<ElementaryOperator<Label, DenseMatrix<Real>>>> matrix_items_;

    };

    template<Domain t_domain>
    struct MeshDiscreteField : DiscreteFieldBase<t_domain>
    {

        explicit
        MeshDiscreteField(
            FieldConcept auto const & field
        )
        :
        DiscreteFieldBase<t_domain>(field)
        {}
        
        Boolean
        operator==(
            MeshDiscreteField const & other
        )
        const = default;
        
        Boolean
        operator!=(
            MeshDiscreteField const & other
        )
        const = default;

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

        // template<FieldConcept auto t_field, Integer t_size, Strategy t_s>
        // void
        // addDegreeOfFreedom(
        //     std::unique_ptr<LinearSystem<t_s>> const & linear_system
        // )
        // {
        //     if (dof_ == nullptr)
        //     {
        //         dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(DegreeOfFreedom<t_domain>::template make<t_field, t_size>(linear_system));
        //     }
        // }

        // template<FieldConcept auto t_field, Integer t_size>
        // void
        // addDegreeOfFreedom()
        // {
        //     if (dof_ == nullptr)
        //     {
        //         // dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(DegreeOfFreedom<t_domain>::template make<t_field, t_size>());
        //         dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(t_size);
        //     }
        // }

        // template<FieldConcept auto t_field, Integer t_size>
        void
        addDegreeOfFreedom(
            Integer size
        )
        {
            if (dof_ == nullptr)
            {
                // dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(DegreeOfFreedom<t_domain>::template make<t_field, t_size>());
                dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(size);
            }
        }

        DegreeOfFreedom<t_domain> const &
        getDegreeOfFreedom()
        const
        {
            if (dof_ == nullptr)
            {
                throw std::runtime_error("No such field data");
            }
            return * dof_;
        }

        DegreeOfFreedom<t_domain> &
        getDegreeOfFreedom()
        {
            if (dof_ == nullptr)
            {
                throw std::runtime_error("No such field data");
            }
            return * dof_;
        }

    private:

        std::unique_ptr<DegreeOfFreedom<t_domain>> dof_;

        std::unique_ptr<std::vector<ExternalLoad>> loads_;

    };

    template<Domain t_domain>
    struct ElementDiscreteField : DiscreteFieldBase<t_domain>
    {

        explicit
        ElementDiscreteField(
            FieldConcept auto const & field
        )
        :
        DiscreteFieldBase<t_domain>(field)
        {}
        
        Boolean
        operator==(
            ElementDiscreteField const & other
        )
        const = default;
        
        Boolean
        operator!=(
            ElementDiscreteField const & other
        )
        const = default;

        Boolean
        hasDegreeOfFreedom()
        const
        {
            return dof_ != nullptr;
        }

        // template<FieldConcept auto t_field, Integer t_size>
        // void
        // addDegreeOfFreedom()
        // {
        //     if (dof_ == nullptr)
        //     {
        //         dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(DegreeOfFreedom<t_domain>::template make<t_field, t_size>());
        //     }
        // }

        // template<FieldConcept auto t_field, Integer t_size>
        void
        addDegreeOfFreedom(
            Integer size
        )
        {
            if (dof_ == nullptr)
            {
                dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(size);
            }
        }

        // template<Element t_element, FieldConcept auto t_field, Basis t_basis, Strategy t_s>
        // void
        // addDegreeOfFreedom(
        //     std::unique_ptr<LinearSystem<t_s>> const & linear_system
        // )
        // {
        //     if (dof_ == nullptr)
        //     {
        //         dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(DegreeOfFreedom<t_domain>::template make<t_element, t_field, t_basis>(linear_system));
        //     }
        // }

        // template<Element t_element, FieldConcept auto t_field, Basis t_basis>
        // void
        // addDegreeOfFreedom()
        // {
        //     if (dof_ == nullptr)
        //     {
        //         dof_ = std::make_unique<DegreeOfFreedom<t_domain>>(DegreeOfFreedom<t_domain>::template make<t_element, t_field, t_basis>());
        //     }
        // }

        DegreeOfFreedom<t_domain> const &
        getDegreeOfFreedom()
        const
        {
            if (dof_ == nullptr)
            {
                throw std::runtime_error("No such field data");
            }
            return * dof_;
        }

        DegreeOfFreedom<t_domain> &
        getDegreeOfFreedom()
        {
            if (dof_ == nullptr)
            {
                throw std::runtime_error("No such field data");
            }
            return * dof_;
        }

    private:

        std::unique_ptr<DegreeOfFreedom<t_domain>> dof_;

    };
    
} // namespace lolita


#endif /* A28EDC0F_E4ED_44BB_A20D_555BFCABEC71 */
