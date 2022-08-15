//
// Created by dsiedel on 23/03/2022.
//

#ifndef FETA_00_FINITE_ELEMENT_DATA_HXX
#define FETA_00_FINITE_ELEMENT_DATA_HXX

#include "feta/core/00_finite_element_geometry.hxx"

namespace feta::core::internal_good
{

    template<
            ElementType ElementTypeArg,
            SupportType SupportTypeArg,
            FiniteElementType FiniteElementTypeArg,
            auto FiniteElementDescriptionArg,
            auto MixedElementDescriptionArg
    >
    struct FiniteElementGeometry3;

    template<
            ElementType ElementTypeArg,
            SupportType SupportTypeArg,
            FiniteElementType FiniteElementTypeArg,
            auto FiniteElementDescriptionArg,
            auto MixedElementDescriptionArg
    >
    struct FiniteElementGeometry3Base
    {

        FiniteElementType const static constexpr finite_element_type = FiniteElementTypeArg;

        static_assert(FiniteElementDescriptionArg.getFiniteElementType() == finite_element_type);

    };

    template<
            ElementType ElementTypeArg,
            auto FiniteElementDescriptionArg,
            auto MixedElementDescriptionArg
    >
    struct FiniteElementGeometry3<
            ElementTypeArg,
            SupportType::Cell,
            FiniteElementType::HybridHighOrder,
            FiniteElementDescriptionArg,
            MixedElementDescriptionArg
    >
    :
    public
    FiniteElementGeometry3Base<
            ElementTypeArg,
            SupportType::Cell,
            FiniteElementType::HybridHighOrder,
            FiniteElementDescriptionArg,
            MixedElementDescriptionArg
    >
    {

        using Base = FiniteElementGeometry3Base<
                ElementTypeArg,
                SupportType::Cell,
                FiniteElementType::HybridHighOrder,
                FiniteElementDescriptionArg,
                MixedElementDescriptionArg
        >;

        BasisType const static constexpr basis_type = BasisType::Monomial;

    };

    template<
            ElementType ElementTypeArg,
            auto FiniteElementDescriptionArg,
            auto MixedElementDescriptionArg
    >
    struct FiniteElementData
    {

        using MixedElementDescription = StructuralType<
                MixedElementDescriptionArg
        >;

        using FiniteElementDescription = typename MixedElementDescription::template FiniteElementDescription<
                FiniteElementDescriptionArg
        >;

    };

//    template<
//            ElementType ElementTypeArg,
//            auto FiniteElementDescriptionArg,
//            auto MixedElementDescriptionArg
//    >
//    struct Element<
//            ElementTypeArg,
//            SupportType::Cell,
//            FiniteElementDescriptionArg,
//            MixedElementDescriptionArg
//    >
//    :
//    public
//    FiniteElementGeometry<
//            ElementTypeArg,
//            SupportType::Cell,
//            FiniteElementDescriptionArg,
//            MixedElementDescriptionArg
//    >
//    {
//
//        using BaseT = FiniteElementGeometry<
//                ElementTypeArg,
//                SupportType::Cell,
//                FiniteElementDescriptionArg,
//                MixedElementDescriptionArg
//        >;
//
//    };

}

#endif //FETA_00_FINITE_ELEMENT_DATA_HXX
