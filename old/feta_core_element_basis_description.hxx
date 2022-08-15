//
// Created by dsiedel on 02/04/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_BASIS_DESCRIPTION_HXX
#define FETA_FETA_CORE_ELEMENT_BASIS_DESCRIPTION_HXX

#include "new/_feta.hxx"

namespace feta::core::element
{

    enum struct Basis
    {

        Monomial,
        Lagrange

    };

    struct BasisDescription
    {

        constexpr
        BasisDescription()
        :
        basis(),
        order()
        {}

        constexpr
        BasisDescription(
                Basis
                basis_arg,
                Indx
                order_arg
        )
        :
        basis(basis_arg),
        order(order_arg)
        {}

        constexpr
        Bool
        operator==(
                BasisDescription const &
                other
        )
        const
        {
            auto eq_0 = basis == other.basis;
            auto eq_1 = order == other.order;
            return eq_0 && eq_1;
        }

        constexpr
        Bool
        operator!=(
                BasisDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

        Basis basis;

        Indx order;

    };

}

#endif //FETA_FETA_CORE_ELEMENT_BASIS_DESCRIPTION_HXX
