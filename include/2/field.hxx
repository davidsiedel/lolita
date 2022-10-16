#ifndef F3E97EE6_11E2_4C64_8FD9_08023BC931D6
#define F3E97EE6_11E2_4C64_8FD9_08023BC931D6

#include "config.hxx"

#include "2/label.hxx"
#include "2/basis.hxx"
#include "2/discretization.hxx"

namespace lolita
{

    template<typename Implementation_>
    struct FieldBase;

    namespace detail
    {
        
        struct FieldType
        {

        private:

            constexpr
            FieldType()
            {}

            constexpr
            Boolean
            operator==(
                FieldType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                FieldType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::FieldBase;

        };
        
    } // namespace detail

    template<typename T>
    concept FieldConcept = std::derived_from<std::decay_t<T>, detail::FieldType>;

    template<typename Implementation_>
    struct FieldBase : detail::FieldType
    {

    protected:

        constexpr
        FieldBase(
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        detail::FieldType(),
        label_(),
        dim_domain_(dim_domain),
        dim_tensor_(dim_tensor)
        {}
        
        constexpr
        FieldBase(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        detail::FieldType(),
        label_(std::forward<std::basic_string_view<Character>>(label)),
        dim_domain_(dim_domain),
        dim_tensor_(dim_tensor)
        {}

    public:

        constexpr
        Boolean
        operator==(
            FieldBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            FieldBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            FieldBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            FieldBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

        constexpr
        Integer
        getDimDomain()
        const
        {
            return dim_domain_;
        }

        constexpr
        Integer
        getDimTensor()
        const
        {
            return dim_tensor_;
        }

        constexpr
        Label const &
        getLabel()
        const
        {
            return label_;
        }

        constexpr
        Boolean
        isTensor(
            Integer dim
        )
        const
        {
            return dim_tensor_ == dim;
        }

        Label label_;

        Integer dim_domain_;

        Integer dim_tensor_;

    };

    struct ImposedField : FieldBase<ImposedField>
    {

    private:

        using Base_ = FieldBase<ImposedField>;

    public:
        
        constexpr
        ImposedField(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor)
        {}

        constexpr
        Boolean
        operator==(
            ImposedField const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            ImposedField const & other
        )
        const = default;

    };

    namespace detail
    {

        template<typename... T>
        struct IsImposedFieldTraits : std::false_type
        {};

        template<>
        struct IsImposedFieldTraits<ImposedField> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept ImposedFieldConcept = detail::IsImposedFieldTraits<T>::value;

    template<DiscretizationConcept Discretization_>
    struct UnknownField : FieldBase<UnknownField<Discretization_>>
    {

        using Discretization = Discretization_;

    private:

        using Base_ = FieldBase<UnknownField<Discretization_>>;

    public:

        constexpr
        UnknownField(
            Integer dim_domain,
            Integer dim_tensor,
            Discretization const & discretization
        )
        :
        Base_(dim_domain, dim_tensor),
        discretization_(discretization)
        {}
        
        constexpr
        UnknownField(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor,
            Discretization const & discretization
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        discretization_(discretization)
        {}

        constexpr
        Boolean
        operator==(
            UnknownField const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            UnknownField const & other
        )
        const = default;

        constexpr
        Discretization const &
        getDiscretization()
        const
        {
            return discretization_;
        }

        Discretization discretization_;

    };

    namespace detail
    {

        template<typename... T>
        struct IsUnknownFieldTraits : std::false_type
        {};

        template<typename T>
        struct IsUnknownFieldTraits<UnknownField<T>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept UnknownFieldConcept = detail::IsUnknownFieldTraits<T>::value;

} // namespace lolita


#endif /* F3E97EE6_11E2_4C64_8FD9_08023BC931D6 */
