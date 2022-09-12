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
    
    template<Element t_element, Domain t_domain>
    struct ElementDegreeOfFreedom
    {

    private:

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

    public:

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
            Integer tag,
            std::unique_ptr<LinearSystem<t_s>> const & linear_system
        )
        {
            auto element_degree_of_freedom = ElementDegreeOfFreedom(tag);
            element_degree_of_freedom.getOffset() = linear_system->getSize();
            linear_system->getSize() += getSize<t_field, t_basis>();
            for (auto i = 0; i < FieldTraits<t_field>::template getSize<t_domain>(); i++)
            {
                element_degree_of_freedom.coefficients_.push_back(Coefficients::template make<t_basis>());
            }
            return element_degree_of_freedom;
        }

        template<Field t_field, Basis t_basis>
        static
        ElementDegreeOfFreedom
        make(
            Integer tag
        )
        {
            auto element_degree_of_freedom = ElementDegreeOfFreedom(tag);
            for (auto i = 0; i < FieldTraits<t_field>::template getSize<t_domain>(); i++)
            {
                element_degree_of_freedom.coefficients_.push_back(Coefficients::template make<t_basis>());
            }
            return element_degree_of_freedom;
        }

    private:

        ElementDegreeOfFreedom(
            Integer tag
        )
        :
        tag_(tag),
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

        Integer
        getTag()
        const
        {
            return tag_;
        }
        
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
        algebra::View<Vector<Real, Coefficients::template getSize<t_basis>()> const>
        get(
            Integer row,
            Integer col
        )
        const
        {
            auto index = FieldTraits<t_field>::template getCols<t_domain>() * row + col;
            return coefficients_[index].template get<t_basis>();
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

        void
        addElementMatrix(
            MatrixConcept<Real> auto && input
        )
        {
            matrix_auxiliaries_.push_back(std::forward<decltype(input)>(input));
        }

        void
        addElementVector(
            VectorConcept<Real> auto && input
        )
        {
            vector_auxiliaries_.push_back(std::forward<decltype(input)>(input));
        }

        void
        addRealParameter(
            Real && input
        )
        {
            parameter_auxiliaries_.push_back(std::forward<Real>(input));
        }

        Matrix<Real> const &
        getElementMatrix(
            Integer tag
        )
        const
        {
            return matrix_auxiliaries_[tag];
        }

        Matrix<Real> &
        getElementMatrix(
            Integer tag
        )
        {
            return matrix_auxiliaries_[tag];
        }

        Vector<Real> const &
        getElementVector(
            Integer tag
        )
        const
        {
            return vector_auxiliaries_[tag];
        }

        Vector<Real> &
        getElementVector(
            Integer tag
        )
        {
            return vector_auxiliaries_[tag];
        }

        Real const &
        getRealParameter(
            Integer tag
        )
        const
        {
            return parameter_auxiliaries_[tag];
        }

        Real &
        getRealParameter(
            Integer tag
        )
        {
            return parameter_auxiliaries_[tag];
        }

    private:

        Integer tag_;

        Natural offset_;

        std::vector<Coefficients> coefficients_;

        std::vector<Matrix<Real>> matrix_auxiliaries_;

        std::vector<Vector<Real>> vector_auxiliaries_;

        std::vector<Real> parameter_auxiliaries_;

        std::vector<Boolean> binary_auxiliaries_;

    };
    
} // namespace lolita


#endif /* A28EDC0F_E4ED_44BB_A20D_555BFCABEC71 */
