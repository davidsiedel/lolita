#ifndef A28EDC0F_E4ED_44BB_A20D_555BFCABEC71
#define A28EDC0F_E4ED_44BB_A20D_555BFCABEC71

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/100_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/200_quadrature.hxx"

namespace lolita
{
        
    template<Basis... t_basis>
    struct BasisTraits;

    template<auto t_discretization>
    struct DiscretizationTraits;

    
    template<Element t_element, Domain t_domain>
    struct Coefficients
    {

        template<Basis t_basis>
        static constexpr
        Integer
        getSize()
        {
            return BasisTraits<t_basis>::template getSize<t_element>();
        }

        template<Basis t_basis>
        static
        Coefficients
        make()
        {
            return Coefficients(getSize<t_basis>());
        }

    private:

        explicit
        Coefficients(
            Integer size
        )
        :
        s1(Vector<Real>::Zero(size)),
        s0(Vector<Real>::Zero(size))
        {}

    public:
    
        Boolean
        operator==(
            Coefficients const & other
        )
        const = default;
        
        Boolean
        operator!=(
            Coefficients const & other
        )
        const = default;

        Vector<Real> const &
        get()
        const
        {
            return s1;
        }

        Vector<Real> &
        get()
        {
            return s1;
        }

        template<Basis t_basis>
        algebra::View<Vector<Real, getSize<t_basis>()> const>
        get()
        const
        {
            return algebra::View<Vector<Real, getSize<t_basis>()> const>(s1.data());
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

    private:

        Vector<Real> s0;

        Vector<Real> s1;

    };

    template<Element t_element, Domain t_domain>
    struct ElementDegreeOfFreedom
    {

        template<Field t_field, Basis t_basis>
        static constexpr
        Integer
        getSize()
        {
            return FieldTraits<t_field>::template getSize<t_domain>() * BasisTraits<t_basis>::template getSize<t_element>();
        }

        template<Field t_field, Basis t_basis, Strategy t_s>
        static
        ElementDegreeOfFreedom
        make(
            std::unique_ptr<LinearSystem<t_s>> const & linear_system
        )
        {
            auto element_degree_of_freedom = ElementDegreeOfFreedom();
            element_degree_of_freedom.getOffset() = linear_system->getSize();
            linear_system->getSize() += getSize<t_field, t_basis>();
            for (auto i = 0; i < FieldTraits<t_field>::template getSize<t_domain>(); i++)
            {
                element_degree_of_freedom.coefficients_.push_back(Coefficients<t_element, t_domain>::template make<t_basis>());
            }
            return element_degree_of_freedom;
        }

        template<Field t_field, Basis t_basis>
        static
        ElementDegreeOfFreedom
        make()
        {
            auto element_degree_of_freedom = ElementDegreeOfFreedom();
            for (auto i = 0; i < FieldTraits<t_field>::template getSize<t_domain>(); i++)
            {
                element_degree_of_freedom.coefficients_.push_back(Coefficients<t_element, t_domain>::template make<t_basis>());
            }
            return element_degree_of_freedom;
        }

    private:

        explicit
        ElementDegreeOfFreedom()
        :
        offset_(0),
        coefficients_()
        {}

    public:
        
        Boolean
        operator==(
            ElementDegreeOfFreedom const & other
        )
        const = default;
        
        Boolean
        operator!=(
            ElementDegreeOfFreedom const & other
        )
        const = default;
        
        Natural &
        getOffset()
        {
            return offset_;
        }
        
        Natural const &
        getOffset()
        const
        {
            return offset_;
        }

        template<Field t_field, Basis t_basis>
        void
        add(
            VectorConcept<Real> auto && input
        )
        {
            auto constexpr offset = BasisTraits<t_basis>::template getSize<t_element>();
            for (auto row = 0; row < FieldTraits<t_field>::template getRows<t_domain>(); row++)
            {
                for (auto col = 0; col < FieldTraits<t_field>::template getCols<t_domain>(); col++)
                {
                    auto index = FieldTraits<t_field>::template getCols<t_domain>() * row + col;
                    coefficients_[index].add(std::forward<decltype(input)>(input).template segment<offset>(offset_ + index * offset));
                }
            }
        }

        template<Field t_field, Basis t_basis>
        algebra::View<Vector<Real, Coefficients<t_element, t_domain>::template getSize<t_basis>()> const>
        get(
            Integer row,
            Integer col
        )
        const
        {
            return coefficients_[FieldTraits<t_field>::template getCols<t_domain>() * row + col].template get<t_basis>();
        }

        template<Field t_field, Basis t_basis>
        Vector<Real, getSize<t_field, t_basis>()>
        get()
        const
        {
            auto constexpr offset = BasisTraits<t_basis>::template getSize<t_element>();
            auto output = Vector<Real, getSize<t_field, t_basis>()>();
            for (auto row = 0; row < FieldTraits<t_field>::template getRows<t_domain>(); row++)
            {
                for (auto col = 0; col < FieldTraits<t_field>::template getCols<t_domain>(); col++)
                {
                    auto index = FieldTraits<t_field>::template getCols<t_domain>() * row + col;
                    output.template segment<offset>(index * offset) = coefficients_[index].template get<t_basis>();
                }
            }
            return output;
        }

        void
        reserve()
        {
            for (auto & c : coefficients_)
            {
                c.reserve();
            }
        }

        void
        recover()
        {
            for (auto & c : coefficients_)
            {
                c.recover();
            }
        }

    private:

        Natural offset_;

        std::vector<Coefficients<t_element, t_domain>> coefficients_;

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

    template<Element t_element, Domain t_domain>
    struct ElementDiscreteField
    {

        using HHH = ElementDegreeOfFreedom<t_element, t_domain>;

        explicit
        ElementDiscreteField(
            Field const & field
        )
        :
        label_(field.getLabel())
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

        utility::Label const &
        getLabel()
        const
        {
            return label_;
        }

        std::basic_string_view<Character>
        getLabelView()
        const
        {
            return label_.view();
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

        template<Field t_field, Basis t_basis, Strategy t_s>
        void
        addDegreeOfFreedom(
            std::unique_ptr<LinearSystem<t_s>> const & linear_system
        )
        {
            dof_ = std::make_unique<HHH>(HHH::template make<t_field, t_basis>(linear_system));
        }

        template<Field t_field, Basis t_basis>
        void
        addDegreeOfFreedom()
        {
            dof_ = std::make_unique<HHH>(HHH::template make<t_field, t_basis>());
        }

        HHH const &
        getDegreeOfFreedom()
        const
        {
            if (dof_ == nullptr)
            {
                throw std::runtime_error("No such field data");
            }
            return * dof_;
        }

        HHH &
        getDegreeOfFreedom()
        {
            if (dof_ == nullptr)
            {
                throw std::runtime_error("No such field data");
            }
            return * dof_;
        }

        void
        addScalar(
            utility::Label const & tag,
            auto &&... scalar
        )
        {
            if (scalar_items_ == nullptr)
            {
                scalar_items_ = std::make_unique<std::vector<ElementaryOperator<utility::Label, Real>>>();
            }
            for (auto & m : * scalar_items_)
            {
                if (m.getTag() == tag)
                {
                    return;
                }
            }
            scalar_items_->push_back(ElementaryOperator<utility::Label, Real>(tag, std::forward<decltype(scalar)>(scalar)...));
        }

        void
        addVector(
            utility::Label const & tag,
            VectorConcept<Real> auto &&... vector
        )
        {
            if (vector_items_ == nullptr)
            {
                vector_items_ = std::make_unique<std::vector<ElementaryOperator<utility::Label, Vector<Real>>>>();
            }
            for (auto & m : * vector_items_)
            {
                if (m.getTag() == tag)
                {
                    return;
                }
            }
            vector_items_->push_back(ElementaryOperator<utility::Label, Vector<Real>>(tag, std::forward<decltype(vector)>(vector)...));
        }

        void
        addMatrix(
            utility::Label const & tag,
            MatrixConcept<Real> auto &&... matrix
        )
        {
            if (matrix_items_ == nullptr)
            {
                matrix_items_ = std::make_unique<std::vector<ElementaryOperator<utility::Label, Matrix<Real>>>>();
            }
            for (auto & m : * matrix_items_)
            {
                if (m.getTag() == tag)
                {
                    return;
                }
            }
            matrix_items_->push_back(ElementaryOperator<utility::Label, Matrix<Real>>(tag, std::forward<decltype(matrix)>(matrix)...));
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

        Vector<Real> const &
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

        Vector<Real> &
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

        Matrix<Real> const &
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

        Matrix<Real> &
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

    private:

        utility::Label const & label_;

        std::unique_ptr<HHH> dof_;

        std::unique_ptr<std::vector<ExternalLoad>> loads_;

        std::unique_ptr<std::vector<ElementaryOperator<utility::Label, Real>>> scalar_items_;

        std::unique_ptr<std::vector<ElementaryOperator<utility::Label, Vector<Real>>>> vector_items_;

        std::unique_ptr<std::vector<ElementaryOperator<utility::Label, Matrix<Real>>>> matrix_items_;

    };

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------

    struct Output
    {

        enum Type
        {
            Success,
            Failure,
            Warning
        };

        static constexpr
        Output
        success()
        {
            return Output(Output::Success);
        }

        static constexpr
        Output
        failure()
        {
            return Output(Output::Failure);
        }

        static constexpr
        Output
        warning()
        {
            return Output(Output::Warning);
        }

        constexpr explicit
        Output(
            Type type
        )
        :
        type_(type)
        {}

        constexpr
        Boolean
        isSuccess()
        const
        {
            return type_ == Output::Success;
        }

        constexpr
        Boolean
        isFailure()
        const
        {
            return type_ == Output::Failure;
        }

        constexpr
        Boolean
        isWarning()
        const
        {
            return type_ == Output::Warning;
        }

        Type type_;

    };

    struct OutputHandler
    {

        OutputHandler()
        :
        output_(Output::success())
        {}

        inline
        void
        setSuccess()
        {
            auto lock = std::scoped_lock<std::mutex>(mutex_);
            output_ = Output::success();
        }

        inline
        void
        setFailure()
        {
            auto lock = std::scoped_lock<std::mutex>(mutex_);
            output_ = Output::failure();
        }

        inline
        void
        setWarning()
        {
            auto lock = std::scoped_lock<std::mutex>(mutex_);
            output_ = Output::warning();
        }

        inline
        Output
        getOutput()
        {
            return output_;
        }

    private:

        Output output_;

        std::mutex mutex_;

    };
    
    template<typename t_Tag, typename t_Value>
    struct ElementaryDataHandler
    {

    private:

        struct Item
        {

            Item(
                t_Tag tag,
                t_Value const & value
            )
            :
            tag_(tag),
            value_(value)
            {}

            Item(
                t_Tag tag,
                t_Value && value
            )
            :
            tag_(tag),
            value_(std::move(value))
            {}
        
            Boolean
            operator==(
                Item const & other
            )
            const = default;
            
            Boolean
            operator!=(
                Item const & other
            )
            const = default;

            t_Tag
            getTag()
            const
            {
                return tag_;
            }

            t_Value const &
            getValue()
            const
            {
                return value_;
            }

            t_Value &
            getValue()
            {
                return value_;
            }

            t_Tag tag_;

            t_Value value_;

        };

    public:

        ElementaryDataHandler()
        :
        items_()
        {}

        void
        add(
            t_Tag tag,
            t_Value && value
        )
        {
            auto find_item = [&] (Item const & item)
            {
                return item.getTag() == tag;
            };
            if(std::find_if(items_.begin(), items_.end(), find_item) == items_.end())
            {
                items_.push_back(Item(tag, std::forward<t_Value>(value)));
            }
        }

        t_Value const &
        get(
            t_Tag tag
        )
        const
        {
            auto find_item = [&] (Item const & item)
            {
                return item.getTag() == tag;
            };
            return std::find_if(items_.begin(), items_.end(), find_item)->getValue();
        }

        t_Value &
        get(
            t_Tag tag
        )
        {
            auto find_item = [&] (Item const & item)
            {
                return item.getTag() == tag;
            };
            return std::find_if(items_.begin(), items_.end(), find_item)->getValue();
        }

    private:

        std::vector<Item> items_;

    };
    
} // namespace lolita


#endif /* A28EDC0F_E4ED_44BB_A20D_555BFCABEC71 */
