#ifndef A1D91AB1_91FC_4E21_AE3D_AE1379B054B0
#define A1D91AB1_91FC_4E21_AE3D_AE1379B054B0

#include "config.hxx"
#include "2/basis.hxx"
#include "2/field.hxx"

namespace lolita
{

    template<typename Implementation_>
    struct DiscretizationBase;

    namespace detail
    {
        
        struct DiscretizationType
        {

        private:

            constexpr
            DiscretizationType()
            {}

            constexpr
            Boolean
            operator==(
                DiscretizationType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                DiscretizationType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::DiscretizationBase;

        };
        
    } // namespace detail

    template<typename T>
    concept DiscretizationConcept = std::derived_from<std::decay_t<T>, detail::DiscretizationType>;

    template<typename Implementation_>
    struct DiscretizationBase : detail::DiscretizationType
    {

    protected:
        
        constexpr
        DiscretizationBase()
        :
        detail::DiscretizationType()
        {}

    public:

        constexpr
        Boolean
        operator==(
            DiscretizationBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            DiscretizationBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            DiscretizationBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            DiscretizationBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

    };

    template<BasisConcept CellBasis_, BasisConcept FaceBasis_>
    struct HybridDiscontinuousGalerkinDiscretization : DiscretizationBase<HybridDiscontinuousGalerkinDiscretization<CellBasis_, FaceBasis_>>
    {

    private:

        using Base_ = DiscretizationBase<HybridDiscontinuousGalerkinDiscretization<CellBasis_, FaceBasis_>>;

    public:
    
        constexpr
        HybridDiscontinuousGalerkinDiscretization(
            CellBasis_ const & cell_basis,
            FaceBasis_ const & face_basis
        )
        :
        Base_(),
        cell_basis_(cell_basis),
        face_basis_(face_basis)
        {}

        constexpr
        CellBasis_ const &
        getCellBasis()
        const
        {
            return cell_basis_;
        }

        constexpr
        FaceBasis_ const &
        getFaceBasis()
        const
        {
            return face_basis_;
        }

        CellBasis_ cell_basis_;

        FaceBasis_ face_basis_;

    };

    namespace detail
    {

        template<typename... T>
        struct IsHybridDiscontinuousGalerkinDiscretizationTraits : std::false_type
        {};

        template<typename... T>
        struct IsHybridDiscontinuousGalerkinDiscretizationTraits<HybridDiscontinuousGalerkinDiscretization<T...>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept HybridDiscontinuousGalerkinDiscretizationConcept = detail::IsHybridDiscontinuousGalerkinDiscretizationTraits<std::decay_t<T>>::value;

    struct LagrangeDiscretization : DiscretizationBase<LagrangeDiscretization>
    {

    private:

        using Base_ = DiscretizationBase<LagrangeDiscretization>;

    public:
    
        constexpr
        LagrangeDiscretization(
            Integer order
        )
        :
        Base_(),
        basis_(order)
        {}

        constexpr
        LagrangeBasis const &
        getBasis()
        const
        {
            return basis_;
        }

        LagrangeBasis basis_;

    };

    namespace detail
    {

        template<typename... T>
        struct IsLagrangeDiscretizationTraits : std::false_type
        {};

        template<>
        struct IsLagrangeDiscretizationTraits<LagrangeDiscretization> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept LagrangeDiscretizationConcept = detail::IsLagrangeDiscretizationTraits<std::decay_t<T>>::value;

    template<typename Implementation_>
    struct ElementBase;

    namespace detail
    {
        
        struct ElementType
        {

        private:

            constexpr
            ElementType()
            {}

            constexpr
            Boolean
            operator==(
                ElementType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                ElementType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::ElementBase;

        };
        
    } // namespace detail

    template<typename T>
    concept ElementConcept = std::derived_from<std::decay_t<T>, detail::ElementType>;

    template<typename Implementation_>
    struct ElementBase : detail::ElementType, FieldBase<ElementBase<Implementation_>>
    {

    protected:

        using FieldBase_ = FieldBase<ElementBase<Implementation_>>;
        
        constexpr
        ElementBase(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        detail::ElementType(),
        FieldBase_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor)
        {}
        
        constexpr
        ElementBase(
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        detail::ElementType(),
        FieldBase_(dim_domain, dim_tensor)
        {}

    public:

        constexpr
        Boolean
        operator==(
            ElementBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            ElementBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            ElementBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            ElementBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

    };

    template<typename CellBasis_, typename FaceBasis_>
    struct HdgElement : ElementBase<HdgElement<CellBasis_, FaceBasis_>>
    {

    private:

        using ElementBase_ = ElementBase<HdgElement<CellBasis_, FaceBasis_>>;

    public:
    
        constexpr
        HdgElement(
            Integer dim_domain,
            Integer dim_tensor,
            CellBasis_ const & cell_basis,
            FaceBasis_ const & face_basis
        )
        :
        ElementBase_(dim_domain, dim_tensor),
        cell_basis_(cell_basis),
        face_basis_(face_basis)
        {}
    
        constexpr
        HdgElement(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor,
            CellBasis_ const & cell_basis,
            FaceBasis_ const & face_basis
        )
        :
        ElementBase_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        cell_basis_(cell_basis),
        face_basis_(face_basis)
        {}

        constexpr
        CellBasis_ const &
        getCellBasis()
        const
        {
            return cell_basis_;
        }

        constexpr
        FaceBasis_ const &
        getFaceBasis()
        const
        {
            return face_basis_;
        }

        CellBasis_ cell_basis_;

        FaceBasis_ face_basis_;

    };

    namespace detail
    {

        template<typename... T>
        struct IsHdgElementTraits : std::false_type
        {};

        template<typename... T>
        struct IsHdgElementTraits<HdgElement<T...>> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept HdgElementConcept = detail::IsHdgElementTraits<std::decay_t<T>>::value;
    
} // namespace lolita

#endif /* A1D91AB1_91FC_4E21_AE3D_AE1379B054B0 */
