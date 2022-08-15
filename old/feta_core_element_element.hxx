//
// Created by dsiedel on 01/04/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_ELEMENT_HXX
#define FETA_FETA_CORE_ELEMENT_ELEMENT_HXX

#include "new/feta_core_element_fe_build.hxx"

namespace feta::core::element
{

    template<ElementDescription E, auto M>
    struct Element : public ElementBuild<E, M>
    {

        using Base = ElementBuild<E, M>;

        Element()
        :
        Base()
        {}

        Bool
        operator==(
                Element const &
                other
        )
        const
        {
            auto eq_0 = this->tag == other.tag;
            return eq_0;
        }

        Bool
        operator!=(
                Element const &
                other
        )
        const
        {
            return !(other == * this);
        }

    };

}

#endif //FETA_FETA_CORE_ELEMENT_ELEMENT_HXX
