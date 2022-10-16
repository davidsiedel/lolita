#ifndef C058C37E_1D45_4C6F_8ECD_DD6B61E78080
#define C058C37E_1D45_4C6F_8ECD_DD6B61E78080

#include "config.hxx"
#include "2/field.hxx"

namespace lolita
{

    template<typename Implementation_, typename UnknownField_>
    struct LinearOperatorBase;

    namespace detail
    {
        
        struct LinearOperatorType
        {

        private:

            constexpr
            LinearOperatorType()
            {}

            constexpr
            Boolean
            operator==(
                LinearOperatorType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                LinearOperatorType const & other
            )
            const
            = default;

            template<typename Implementation_, typename UnknownField_>
            friend class lolita::LinearOperatorBase;

        };
        
    } // namespace detail

    template<typename T>
    concept LinearOperatorConcept = std::derived_from<std::decay_t<T>, detail::LinearOperatorType>;

    template<typename Implementation_, typename UnknownField_>
    struct LinearOperatorBase : detail::LinearOperatorType, FieldBase<LinearOperatorBase<Implementation_, UnknownField_>>
    {

    private:

        using Base_ = FieldBase<LinearOperatorBase<Implementation_, UnknownField_>>;

    protected:

        constexpr 
        LinearOperatorBase(
            Integer dim_domain,
            Integer dim_tensor,
            UnknownField_ const & unknown_field
        )
        :
        detail::LinearOperatorType(),
        Base_(dim_domain, dim_tensor),
        row_(-1),
        col_(-1),
        unknown_field_(unknown_field)
        {}

        constexpr 
        LinearOperatorBase(
            Integer dim_domain,
            Integer dim_tensor,
            Integer row,
            Integer col,
            UnknownField_ const & unknown_field
        )
        :
        detail::LinearOperatorType(),
        Base_(dim_domain, dim_tensor),
        row_(row),
        col_(col),
        unknown_field_(unknown_field)
        {}

        constexpr 
        LinearOperatorBase(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor,
            UnknownField_ const & unknown_field
        )
        :
        detail::LinearOperatorType(),
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        row_(-1),
        col_(-1),
        unknown_field_(unknown_field)
        {}

        constexpr 
        LinearOperatorBase(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor,
            Integer row,
            Integer col,
            UnknownField_ const & unknown_field
        )
        :
        detail::LinearOperatorType(),
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        row_(row),
        col_(col),
        unknown_field_(unknown_field)
        {}

    public:

        constexpr
        Boolean
        operator==(
            LinearOperatorBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            LinearOperatorBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            LinearOperatorBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            LinearOperatorBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

        constexpr
        Integer
        getRow()
        const
        {
            return row_;
        }

        constexpr
        Integer
        getCol()
        const
        {
            return col_;
        }

        constexpr
        UnknownField_ const &
        getField()
        const
        {
            return unknown_field_;
        }

        Integer row_;

        Integer col_;

        UnknownField_ unknown_field_;

    };

    template<typename CellField_>
    struct GradientOperator : LinearOperatorBase<GradientOperator<CellField_>, CellField_>
    {

    private:

        using Base_ = LinearOperatorBase<GradientOperator<CellField_>, CellField_>;

    public:

        explicit constexpr
        GradientOperator(
            CellField_ const & cell_field
        )
        :
        Base_("Gradient", cell_field.getDimDomain(), cell_field.getDimTensor() + 1, cell_field)
        {}

        constexpr 
        GradientOperator(
            std::basic_string_view<Character> && label,
            CellField_ const & cell_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), cell_field.getDimDomain(), cell_field.getDimTensor() + 1, cell_field)
        {}

        constexpr
        Boolean
        operator==(
            GradientOperator const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            GradientOperator const & other
        )
        const = default;
        
    };

    namespace detail
    {

        template<typename... T>
        struct IsGradientOperatorTraits : std::false_type
        {};

        template<typename... T>
        struct IsGradientOperatorTraits<GradientOperator<T...>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept GradientOperatorConcept = detail::IsGradientOperatorTraits<std::decay_t<T>>::value;

    template<typename CellFieldT>
    struct TraceOperator : LinearOperatorBase<TraceOperator<CellFieldT>, CellFieldT>
    {

    private:

        using Base_ = LinearOperatorBase<TraceOperator<CellFieldT>, CellFieldT>;

    public:

        explicit constexpr
        TraceOperator(
            CellFieldT const & cell_field
        )
        :
        Base_(cell_field.getDomainDim(), cell_field.getTensorDim(), cell_field)
        {}
        
        constexpr
        TraceOperator(
            std::basic_string_view<Character> && label,
            CellFieldT const & cell_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), cell_field.getDomainDim(), cell_field.getTensorDim(), cell_field)
        {}

        constexpr
        Boolean
        operator==(
            TraceOperator const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            TraceOperator const & other
        )
        const = default;
        
    };

    namespace detail
    {

        template<typename... T>
        struct IsTraceOperatorTraits : std::false_type
        {};

        template<typename... T>
        struct IsTraceOperatorTraits<TraceOperator<T...>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept TraceOperatorConcept = detail::IsTraceOperatorTraits<std::decay_t<T>>::value;

    template<typename CellField_>
    struct SmallStrainOperator : LinearOperatorBase<SmallStrainOperator<CellField_>, CellField_>
    {

    private:

        using Base_ = LinearOperatorBase<SmallStrainOperator<CellField_>, CellField_>;

    public:

        explicit constexpr
        SmallStrainOperator(
            CellField_ const & cell_field
        )
        :
        Base_("SmallStrain", cell_field.getDimDomain(), cell_field.getDimTensor() + 1, cell_field)
        {}

        constexpr 
        SmallStrainOperator(
            std::basic_string_view<Character> && label,
            CellField_ const & cell_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), cell_field.getDimDomain(), cell_field.getDimTensor() + 1, cell_field)
        {}

        constexpr
        Boolean
        operator==(
            SmallStrainOperator const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            SmallStrainOperator const & other
        )
        const = default;
        
    };

    namespace detail
    {

        template<typename... T>
        struct IsSmallStrainOperatorTraits : std::false_type
        {};

        template<typename... T>
        struct IsSmallStrainOperatorTraits<SmallStrainOperator<T...>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept SmallStrainOperatorConcept = detail::IsSmallStrainOperatorTraits<std::decay_t<T>>::value;

    template<typename CellField_>
    struct LargeStrainOperator : LinearOperatorBase<LargeStrainOperator<CellField_>, CellField_>
    {

    private:

        using Base_ = LinearOperatorBase<LargeStrainOperator<CellField_>, CellField_>;

    public:

        explicit constexpr
        LargeStrainOperator(
            CellField_ const & cell_field
        )
        :
        Base_("LargeStrain", cell_field.getDimDomain(), cell_field.getDimTensor() + 1, cell_field)
        {}

        constexpr 
        LargeStrainOperator(
            std::basic_string_view<Character> && label,
            CellField_ const & cell_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), cell_field.getDimDomain(), cell_field.getDimTensor() + 1, cell_field)
        {}

        constexpr
        Boolean
        operator==(
            LargeStrainOperator const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            LargeStrainOperator const & other
        )
        const = default;
        
    };

    namespace detail
    {

        template<typename... T>
        struct IsLargeStrainOperatorTraits : std::false_type
        {};

        template<typename... T>
        struct IsLargeStrainOperatorTraits<LargeStrainOperator<T...>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept LargeStrainOperatorConcept = detail::IsLargeStrainOperatorTraits<std::decay_t<T>>::value;

    // template<typename CellFieldT, typename FaceFieldT>
    // struct HybridDiscontinuousGalerkinGradientOperator : LinearOperatorBase<HybridDiscontinuousGalerkinGradientOperator<CellFieldT, FaceFieldT>, CellFieldT, FaceFieldT>
    // {

    // private:

    //     using Base_ = LinearOperatorBase<HybridDiscontinuousGalerkinGradientOperator<CellFieldT, FaceFieldT>, CellFieldT, FaceFieldT>;

    // public:

    //     constexpr
    //     HybridDiscontinuousGalerkinGradientOperator(
    //         CellFieldT const & cell_field,
    //         CellFieldT const & face_field
    //     )
    //     :
    //     Base_(cell_field.getDomainDim(), cell_field.getTensorDim() + 1, cell_field, face_field)
    //     {}

    //     constexpr
    //     HybridDiscontinuousGalerkinGradientOperator(
    //         std::basic_string_view<Character> && label,
    //         CellFieldT const & cell_field,
    //         CellFieldT const & face_field
    //     )
    //     :
    //     Base_(std::forward<std::basic_string_view<Character>>(label), cell_field.getDomainDim(), cell_field.getTensorDim() + 1, cell_field, face_field)
    //     {}

    //     constexpr
    //     Boolean
    //     operator==(
    //         HybridDiscontinuousGalerkinGradientOperator const & other
    //     )
    //     const = default;

    //     constexpr
    //     Boolean
    //     operator!=(
    //         HybridDiscontinuousGalerkinGradientOperator const & other
    //     )
    //     const = default;
        
    // };

    // namespace detail
    // {

    //     template<typename... T>
    //     struct IsHybridDiscontinuousGalerkinGradientOperatorTraits : std::false_type
    //     {};

    //     template<typename... T>
    //     struct IsHybridDiscontinuousGalerkinGradientOperatorTraits<HybridDiscontinuousGalerkinGradientOperator<T...>> : std::true_type
    //     {};
        
    // } // namespace detail

    // template<typename T>
    // concept HybridDiscontinuousGalerkinGradientOperatorConcept = detail::IsHybridDiscontinuousGalerkinGradientOperatorTraits<T>::value;

    // template<typename CellFieldT, typename FaceFieldT>
    // struct HybridDiscontinuousGalerkinStabilizationOperator : LinearOperatorBase<HybridDiscontinuousGalerkinStabilizationOperator<CellFieldT, FaceFieldT>, CellFieldT, FaceFieldT>
    // {

    // private:

    //     using Base_ = LinearOperatorBase<HybridDiscontinuousGalerkinStabilizationOperator<CellFieldT, FaceFieldT>, CellFieldT, FaceFieldT>;

    // public:

    //     constexpr
    //     HybridDiscontinuousGalerkinStabilizationOperator(
    //         CellFieldT const & cell_field,
    //         CellFieldT const & face_field
    //     )
    //     :
    //     Base_(face_field.getDomainDim(), cell_field.getTensorDim(), cell_field, face_field)
    //     {}

    //     constexpr
    //     HybridDiscontinuousGalerkinStabilizationOperator(
    //         std::basic_string_view<Character> && label,
    //         CellFieldT const & cell_field,
    //         CellFieldT const & face_field
    //     )
    //     :
    //     Base_(std::forward<std::basic_string_view<Character>>(label), face_field.getDomainDim(), cell_field.getTensorDim(), cell_field, face_field)
    //     {}

    //     constexpr
    //     Boolean
    //     operator==(
    //         HybridDiscontinuousGalerkinStabilizationOperator const & other
    //     )
    //     const = default;

    //     constexpr
    //     Boolean
    //     operator!=(
    //         HybridDiscontinuousGalerkinStabilizationOperator const & other
    //     )
    //     const = default;
        
    // };

    // namespace detail
    // {

    //     template<typename... T>
    //     struct IsHybridDiscontinuousGalerkinStabilizationOperatorTraits : std::false_type
    //     {};

    //     template<typename... T>
    //     struct IsHybridDiscontinuousGalerkinStabilizationOperatorTraits<HybridDiscontinuousGalerkinStabilizationOperator<T...>> : std::true_type
    //     {};
        
    // } // namespace detail

    // template<typename T>
    // concept HybridDiscontinuousGalerkinStabilizationOperatorConcept = detail::IsHybridDiscontinuousGalerkinStabilizationOperatorTraits<T>::value;

} // namespace lolita

#endif /* C058C37E_1D45_4C6F_8ECD_DD6B61E78080 */
