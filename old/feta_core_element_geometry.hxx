//
// Created by dsiedel on 13/04/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_GEOMETRY_HXX
#define FETA_FETA_CORE_ELEMENT_GEOMETRY_HXX

#include "new/feta_feta_types.hxx"

namespace feta::core::element
{

    namespace shape
    {

        enum struct Shape
        {

            Point,
            Segment,
            Triangle,
            Quadrangle,
            Polygon,
            Tetrahedron,
            Pyramid,
            Brick,
            Prism,
            Polyhedron,

        };

        namespace detail
        {

            static constexpr inline
            auto
            setShape(
                    Shape
                    shape_arg
            )
            {
                if (isInn(shape_arg, Shape::Point)) {
                    return std::pair(0, 1);
                }
                else if (isInn(shape_arg, Shape::Segment)) {
                    return std::pair(1, 1);
                }
                else if (isInn(shape_arg, Shape::Triangle, Shape::Quadrangle, Shape::Polygon)) {
                    return std::pair(2, 2);
                }
                else {
                    return std::pair(3, 3);
                }
            }

            static constexpr inline
            auto
            ordShape(
                    Shape
                    shape_arg
            )
            {
                if (isInn(shape_arg, Shape::Point)) {
                    return 0;
                }
                else if (isInn(shape_arg, Shape::Segment)) {
                    return 1;
                }
                else if (isInn(shape_arg, Shape::Triangle, Shape::Quadrangle, Shape::Polygon)) {
                    return 2;
                }
                else {
                    return 3;
                }
            }

            static constexpr inline
            auto
            dimShape(
                    Shape
                    shape_arg
            )
            {
                if (isInn(shape_arg, Shape::Point)) {
                    return 1;
                }
                else if (isInn(shape_arg, Shape::Segment)) {
                    return 1;
                }
                else if (isInn(shape_arg, Shape::Triangle, Shape::Quadrangle, Shape::Polygon)) {
                    return 2;
                }
                else {
                    return 3;
                }
            }

        }

        struct ShapeDescription
        {

            constexpr
            ShapeDescription()
                    :
                    shape(),
                    ord_shape(),
                    dim_shape()
            {}

            constexpr explicit
            ShapeDescription(
                    Shape
                    shape_arg
            )
                    :
                    shape(shape_arg),
                    ord_shape(shape::detail::ordShape(shape_arg)),
                    dim_shape(shape::detail::dimShape(shape_arg))
            {}

            constexpr
            Bool
            operator==(
                    ShapeDescription const &
                    other
            )
            const
            {
                auto eq_0 = shape == other.shape;
                auto eq_1 = ord_shape == other.ord_shape;
                auto eq_2 = dim_shape == other.dim_shape;
                return eq_0 && eq_1 && eq_2;
            }

            constexpr
            Bool
            operator!=(
                    ShapeDescription const &
                    other
            )
            const
            {
                return !(other == * this);
            }

            Shape shape;

            Indx ord_shape;

            Indx dim_shape;

        };

        auto const static constexpr pnt = ShapeDescription(Shape::Point);
        auto const static constexpr seg = ShapeDescription(Shape::Segment);
        auto const static constexpr tri = ShapeDescription(Shape::Triangle);
        auto const static constexpr qua = ShapeDescription(Shape::Quadrangle);
        auto const static constexpr tet = ShapeDescription(Shape::Tetrahedron);

    }

    using Shape = shape::Shape;

    using ShapeDescription = shape::ShapeDescription;

}

#endif //FETA_FETA_CORE_ELEMENT_GEOMETRY_HXX
