/**
 * @file field.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-04-05
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef D40C81AD_38A1_4446_A043_2CDF6516250E
#define D40C81AD_38A1_4446_A043_2CDF6516250E

#include "config.hxx"
#include "tensor.hxx"
#include "geometry/frame.hxx"

namespace lolita::physics
{

    namespace internal
    {

        template<Integer... i_>
        using FieldStaticTensor = tensor::StaticTensor<Real, i_...>;

        template<Integer rank_, Integer dim_>
        struct FieldStaticTensorTraits
        {

            using Type = tuple_cat_t<FieldStaticTensor<dim_>, typename FieldStaticTensorTraits<rank_ - 1, dim_>::Type>;

        };

        template<Integer dim_>
        struct FieldStaticTensorTraits<1, dim_>
        {

            using Type = FieldStaticTensor<dim_>;

        };

        template<Integer dim_>
        struct FieldStaticTensorTraits<0, dim_>
        {

            using Type = FieldStaticTensor<>;

        };
        
    } // namespace internal
    
    
    template<Integer dim_, Integer rank_>
    struct Field
    {

        using DynamicTensor = tensor::DynamicTensor<Real, rank_>;

        using StaticTensor = internal::FieldStaticTensorTraits<rank_, dim_>::Type;

        static constexpr
        Integer
        getDimField()
        {
            return dim_;
        }

        static constexpr
        Integer
        getRank()
        {
            return rank_;
        }

    };
    
    namespace internal
    {

        template<typename T_>
        struct IsField : std::false_type
        {};

        template<Integer... i_>
        struct IsField<Field<i_...>> : std::true_type
        {};

    } // namespace internal
    
    template<typename T_>
    concept FieldConcept = internal::IsField<T_>::value;

    template<geometry::FrameConcept Frame_, Integer rank_>
    using CartesianField = Field<Frame_::getDimEuclidean(), rank_>;

    template<FieldConcept Field_>
    using Gradient = Field<Field_::getDim(), Field_::getRank() + 1>;

    template<FieldConcept Field_>
    using Divergence = Field<Field_::getDim(), Field_::getRank() - 1>;

    template<FieldConcept Field_>
    using Identity = Field_;

} // namespace lolita::physics


#endif /* D40C81AD_38A1_4446_A043_2CDF6516250E */
