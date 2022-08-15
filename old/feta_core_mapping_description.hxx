//
// Created by dsiedel on 06/04/2022.
//

#ifndef FETA_FETA_CORE_MAPPING_DESCRIPTION_HXX
#define FETA_FETA_CORE_MAPPING_DESCRIPTION_HXX

#include "new/_feta.hxx"
#include "new/_feta_collection.hxx"

namespace feta::core
{

    enum struct Mapping
    {

        Gradient,
        SymmetricGradient,
        Divergence,
        SmallStrain,
        LargeStrain,
        Identity,

    };

    auto const static constexpr map_grd = Mapping::Gradient;
    auto const static constexpr map_grd_sym = Mapping::SymmetricGradient;
    auto const static constexpr map_div = Mapping::Divergence;
    auto const static constexpr map_eye = Mapping::Identity;
    auto const static constexpr map_str_sym = Mapping::SmallStrain;
    auto const static constexpr map_str = Mapping::LargeStrain;

    template<template<Mapping> typename T = Self<Mapping>::TemplatedType>
    using Mappings = Collection<
            T<Mapping::Gradient>,
            T<Mapping::SymmetricGradient>,
            T<Mapping::Divergence>,
            T<Mapping::SmallStrain>,
            T<Mapping::LargeStrain>,
            T<Mapping::Identity>
    >;

    struct MappingDescription
    {

        constexpr
        MappingDescription()
        :
        mapping(),
        ord_mapping(-1)
        {}

        constexpr
        MappingDescription(
                Mapping
                mpa_arg
        )
        :
        mapping(mpa_arg),
        ord_mapping(-1)
        {}

        constexpr
        MappingDescription(
                Mapping
                mpa_arg,
                Intg
                order_arg
        )
        :
        mapping(mpa_arg),
        ord_mapping(order_arg)
        {}

        constexpr
        Bool
        operator==(
                MappingDescription const &
                other
        )
        const
        {
            auto eq_0 = mapping == other.mapping;
            auto eq_1 = ord_mapping == other.ord_mapping;
            return eq_0 && eq_1;
        }

        constexpr
        Bool
        operator!=(
                MappingDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

        Mapping mapping;

        Intg ord_mapping;

    };

}

#endif //FETA_FETA_CORE_MAPPING_DESCRIPTION_HXX
