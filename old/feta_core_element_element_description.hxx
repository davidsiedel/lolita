//
// Created by dsiedel on 30/03/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_ELEMENT_DESCRIPTION_HXX
#define FETA_FETA_CORE_ELEMENT_ELEMENT_DESCRIPTION_HXX

#include "new/_feta.hxx"
#include "new/feta_core_element_shape_description.hxx"

namespace feta::core::element
{

    struct ElementDescription
    {

        constexpr
        ElementDescription()
        :
                shape_description(ShapeDescription()),
                ord_element(),
                tag()
        {}

        constexpr
        ElementDescription(
                Shape
                shape_arg,
                Indx
                order_arg,
                Indx
                tag_arg
        )
        :
                shape_description(ShapeDescription(shape_arg)),
                ord_element(order_arg),
                tag(tag_arg)
        {}

        constexpr inline
        Bool
        operator==(
                ElementDescription const &
                other
        )
        const
        {
            auto eq_0 = shape_description == other.shape_description;
            auto eq_1 = ord_element == other.ord_element;
            auto eq_3 = tag == other.tag;
            return eq_0 && eq_1 && eq_3;
        }

        constexpr inline
        Bool
        operator!=(
                ElementDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

        constexpr inline
        auto
        shape()
        const
        {
            return shape_description.shape;
        }

        constexpr inline
        auto
        ordShape()
        const
        {
            return shape_description.ord_shape;
        }

        constexpr inline
        auto
        dimShape()
        const
        {
            return shape_description.dim_shape;
        }

        ShapeDescription shape_description;

        Indx ord_element;

        Indx tag;

    };

    auto const static constexpr pnt_0 = ElementDescription(Shape::Point, 0, 0);
    auto const static constexpr seg_2 = ElementDescription(Shape::Segment, 1, 0);
    auto const static constexpr tri_3 = ElementDescription(Shape::Triangle, 1, 0);
    auto const static constexpr qua_4 = ElementDescription(Shape::Quadrangle, 1, 0);
    auto const static constexpr tet_4 = ElementDescription(Shape::Tetrahedron, 1, 0);

}

#endif //FETA_FETA_CORE_ELEMENT_ELEMENT_DESCRIPTION_HXX
