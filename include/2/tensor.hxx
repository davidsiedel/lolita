#ifndef C8A31CE9_BBB4_4A85_9C83_422B5B663D52
#define C8A31CE9_BBB4_4A85_9C83_422B5B663D52

#include "config.hxx"

namespace lolita
{

    template<typename Implementation_>
    struct TensorBase;

    namespace detail
    {
        
        struct TensorType
        {

        private:

            constexpr
            TensorType()
            {}

            constexpr
            Boolean
            operator==(
                TensorType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                TensorType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::TensorBase;

        };
        
    } // namespace detail

    template<typename T>
    concept TensorConcept = std::derived_from<std::decay_t<T>, detail::TensorType>;
    
    template<typename Implementation_>
    struct TensorBase : detail::TensorType
    {

    protected:

        constexpr
        TensorBase(
            Integer dim
        )
        :
        detail::TensorType(),
        dim_(dim)
        {}

    public:

        constexpr
        Boolean
        operator==(
            TensorBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            TensorBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            TensorBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            TensorBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }

        Integer dim_;

    };

    struct Scalar : TensorBase<Scalar>
    {

    private:

        using Base_ = TensorBase<Scalar>;

    public:

        constexpr
        Scalar()
        :
        Base_(0)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsScalarTraits : std::false_type
        {};

        template<>
        struct IsScalarTraits<Scalar> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept ScalarConcept = detail::IsScalarTraits<T>::value;

    struct VectorTensorField : TensorBase<VectorTensorField>
    {

    private:

        using Base_ = TensorBase<VectorTensorField>;

    public:

        constexpr
        VectorTensorField()
        :
        Base_(1)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsVectorTensorFieldTraits : std::false_type
        {};

        template<>
        struct IsVectorTensorFieldTraits<VectorTensorField> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept VectorTensorFieldConcept = detail::IsVectorTensorFieldTraits<T>::value;

    struct VectorTensorField : TensorBase<VectorTensorField>
    {

    private:

        using Base_ = TensorBase<VectorTensorField>;

    public:

        constexpr
        VectorTensorField()
        :
        Base_(1)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsVectorTensorFieldTraits : std::false_type
        {};

        template<>
        struct IsVectorTensorFieldTraits<VectorTensorField> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept VectorTensorFieldConcept = detail::IsVectorTensorFieldTraits<T>::value;

    struct SolidTensor : TensorBase<SolidTensor>
    {

    private:

        using Base_ = TensorBase<SolidTensor>;

    public:

        constexpr
        SolidTensor()
        :
        Base_(1)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsSolidTensorTraits : std::false_type
        {};

        template<>
        struct IsSolidTensorTraits<SolidTensor> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept SolidTensorConcept = detail::IsSolidTensorTraits<T>::value;

    struct TensorLibrary
    {

        template<template<TensorConcept auto, auto...> typename T, auto... t_arg>
        using Tensors = std::tuple<
            T<Scalar{}, t_arg...>,
            T<VectorTensorField{}, t_arg...>,
            T<VectorTensorField{}, t_arg...>,
            T<SolidTensor{}, t_arg...>
        >;

    };

} // namespace lolita

#endif /* C8A31CE9_BBB4_4A85_9C83_422B5B663D52 */
