//
// Created by dsiedel on 29/03/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_QUADRATURE_DESCRIPTION_HXX
#define FETA_FETA_CORE_ELEMENT_QUADRATURE_DESCRIPTION_HXX

#include "new/_feta.hxx"

namespace feta::core
{

    enum struct Quadrature
    {

        Gauss

    };

    struct QuadratureDescription
    {

        constexpr
        QuadratureDescription()
        :
        quadrature(Quadrature::Gauss),
        order(0)
        {}

        constexpr
        QuadratureDescription(
                Quadrature
                quadrature_arg,
                Indx
                order_arg
        )
        :
        quadrature(quadrature_arg),
        order(order_arg)
        {}

        constexpr
        Bool
        operator==(
                QuadratureDescription const &
                other
        )
        const
        {
            auto eq_0 = quadrature == other.quadrature;
            auto eq_1 = order == other.order;
            return eq_0 && eq_1;
        }

        constexpr
        Bool
        operator!=(
                QuadratureDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

        Quadrature quadrature;

        Indx order;

    };

}

#endif //FETA_FETA_CORE_ELEMENT_QUADRATURE_DESCRIPTION_HXX
