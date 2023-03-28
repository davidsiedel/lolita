#ifndef E68CB6A8_1D3D_4E96_928D_FB405436F896
#define E68CB6A8_1D3D_4E96_928D_FB405436F896

#include "config.hxx"
#include "numerics.hxx"
#include "geometry/frame.hxx"

namespace lolita::geometry
{

    namespace shape
    {
        
        /**
         * @brief 
         * 
         */
        struct Point
        {
            
            static constexpr
            Integer
            getDimShape()
            {
                return 0;
            }
            
            static constexpr
            Integer
            getNumNodes()
            {
                return 1;
            }

        };

        namespace internal
        {
        
            template<typename T_>
            struct PointTraits : std::false_type
            {};
            
            template<>
            struct PointTraits<Point> : std::true_type
            {};
            
        } // namespace internal

        /**
         * @brief 
         * 
         */
        struct Curve
        {

            template<template<typename, typename...> typename T_, typename... U_>
            using InnerNeighborhood = std::tuple<
                std::tuple<
                    std::array<T_<Point, U_...>, 2>
                >
            >;
            
            static constexpr
            Integer
            getDimShape()
            {
                return 1;
            }
            
            static constexpr
            Integer
            getNumNodes()
            {
                return 2;
            }

        };

        namespace internal
        {

            template<typename T_>
            struct CurveTraits : std::false_type
            {};
            
            template<>
            struct CurveTraits<Curve> : std::true_type
            {};
            
        } // namespace internal

        /**
         * @brief 
         * 
         * @tparam num_points_ 
         */
        template<Integer num_points_>
        requires(num_points_ >= 3)
        struct Facet
        {

            template<template<typename, typename...> typename T_, typename... U_>
            using InnerNeighborhood = std::tuple<
                std::tuple<
                    std::array<T_<Curve, U_...>, num_points_>
                >,
                std::tuple<
                    std::array<T_<Point, U_...>, num_points_>
                >
            >;
            
            static constexpr
            Integer
            getDimShape()
            {
                return 2;
            }
            
            static constexpr
            Integer
            getNumNodes()
            {
                return num_points_;
            }

        };

        namespace internal
        {

            template<typename T_>
            struct FacetTraits : std::false_type
            {};
            
            template<Integer arg_>
            struct FacetTraits<Facet<arg_>> : std::true_type
            {};
            
        } // namespace internal

        namespace internal
        {

            struct Faces
            {
                
                constexpr
                Faces(
                    Integer num_nodes,
                    Integer num_faces
                )
                :
                num_nodes_(num_nodes),
                num_faces_(num_faces)
                {}

                Integer const num_nodes_;

                Integer const num_faces_;

            };

        } // namespace internal
        
        /**
         * @brief 
         * 
         * @tparam num_points_ 
         * @tparam num_curves_ 
         * @tparam faces_ 
         */
        template<Integer num_points_, Integer num_curves_, auto... faces_>
        struct Solid
        {

            template<template<typename, typename...> typename T_, typename... U_>
            using InnerNeighborhood = std::tuple<
                std::tuple<
                    std::array<T_<Facet<faces_.num_nodes_>, U_...>, faces_.num_faces_>...
                >,
                std::tuple<
                    std::array<T_<Curve, U_...>, num_curves_>
                >,
                std::tuple<
                    std::array<T_<Point, U_...>, num_points_>
                >
            >;
            
            static constexpr
            Integer
            getDimShape()
            {
                return 3;
            }
            
            static constexpr
            Integer
            getNumNodes()
            {
                return num_points_;
            }

        };

        namespace internal
        {
        
            template<typename T_>
            struct SolidTraits : std::false_type
            {};
            
            template<auto... args_>
            struct SolidTraits<Solid<args_...>> : std::true_type
            {};
            
        } // namespace internal

    } // namespace shape

    template<typename T_>
    concept PointConcept = shape::internal::PointTraits<T_>::value;

    template<typename T_>
    concept CurveConcept = shape::internal::CurveTraits<T_>::value;

    template<typename T_>
    concept FacetConcept = shape::internal::FacetTraits<T_>::value;

    template<typename T_>
    concept SolidConcept = shape::internal::SolidTraits<T_>::value;

    template<typename T_>
    concept ShapeConcept = PointConcept<T_> || CurveConcept<T_> || FacetConcept<T_> || SolidConcept<T_>;

    using Node = shape::Point;
    using Segment = shape::Curve;
    using Triangle = shape::Facet<3>;
    using Quadrangle = shape::Facet<4>;
    using Tetrahedron = shape::Solid<4, 6, shape::internal::Faces(3, 4)>;
    using Hexahedron = shape::Solid<8, 12, shape::internal::Faces(4, 8)>;

    template<template<typename...> typename T_, typename... U_>
    using Shapes = std::tuple<
        std::tuple<
            T_<Node, U_...>
        >,
        std::tuple<
            T_<Segment, U_...>
        >,
        std::tuple<
            T_<Triangle, U_...>,
            T_<Quadrangle, U_...>
        >,
        std::tuple<
            T_<Tetrahedron, U_...>
        >
    >;

    struct ShapeCoordinates
    {
        
        constexpr
        ShapeCoordinates(
            Integer i,
            Integer j
        )
        :
        i_(i),
        j_(j)
        {}

        Integer const i_;

        Integer const j_;

    };

    struct ShapeLibraryTraits
    {

        template<Integer i_, Integer j_>
        using Shape = std::tuple_element_t<j_, std::tuple_element_t<i_, Shapes<TypeView>>>::type;

        template<ShapeConcept Shape_>
        static constexpr
        Boolean
        hasShape()
        {
            using Shapes_ = std::tuple_element_t<Shape_::getDimShape(), Shapes<TypeView>>;
            auto tag = false;
            auto set_tag = [&] <Integer i_ = 0> (
                auto & set_tag_
            )
            constexpr mutable
            {
                if (std::same_as<Shape_, typename std::tuple_element_t<i_, Shapes_>::type>)
                {
                    tag = true;
                }
                if constexpr (i_ < std::tuple_size_v<Shapes_> - 1)
                {
                    set_tag_.template operator()<i_ + 1>(set_tag_);
                }
            };
            set_tag(set_tag);
            return tag;
        }

        template<ShapeConcept Shape_>
        static constexpr
        ShapeCoordinates
        getShapeCoordinates()
        requires(hasShape<Shape_>())
        {
            using Shapes_ = std::tuple_element_t<Shape_::getDimShape(), Shapes<TypeView>>;
            auto tag = -1;
            auto set_tag = [&] <Integer i_ = 0> (
                auto & set_tag_
            )
            constexpr mutable
            {
                if (std::same_as<Shape_, typename std::tuple_element_t<i_, Shapes_>::type>)
                {
                    tag = i_;
                }
                if constexpr (i_ < std::tuple_size_v<Shapes_> - 1)
                {
                    set_tag_.template operator()<i_ + 1>(set_tag_);
                }
            };
            set_tag(set_tag);
            return ShapeCoordinates(Shape_::getDimShape(), tag);
        }

        template<ShapeConcept Shape_>
        static constexpr
        Integer
        getShapeTag()
        requires(hasShape<Shape_>())
        {
            using Shapes_ = std::tuple_element_t<Shape_::getDimShape(), Shapes<TypeView>>;
            auto tag = -1;
            auto set_tag = [&] <Integer i_ = 0> (
                auto & set_tag_
            )
            constexpr mutable
            {
                if (std::same_as<Shape_, typename std::tuple_element_t<i_, Shapes_>::type>)
                {
                    tag = i_;
                }
                if constexpr (i_ < std::tuple_size_v<Shapes_> - 1)
                {
                    set_tag_.template operator()<i_ + 1>(set_tag_);
                }
            };
            set_tag(set_tag);
            return tag;
        }

    };

    template<ShapeConcept Shape_>
    requires(!PointConcept<Shape_>)
    struct ShapeInnerNeighborhoodTraits
    {

        template<template<typename, typename...> typename T_, typename... U_>
        using InnerNeighborhood = typename Shape_::template InnerNeighborhood<T_, U_...>;

        template<Integer i_, Integer j_>
        using InnerNeighbor = std::tuple_element_t<j_, std::tuple_element_t<i_, InnerNeighborhood<TypeView>>>::value_type::type;

        template<ShapeConcept Other_>
        static constexpr
        Boolean
        hasShape()
        {
            if (Other_::getDimShape() >= Shape_::getDimShape())
            {
                return false;
            }
            else
            {
                using Shapes_ = std::tuple_element_t<Shape_::getDimShape() - Other_::getDimShape() - 1, InnerNeighborhood<TypeView>>;
                auto tag = false;
                auto set_tag = [&] <Integer i_ = 0> (
                    auto & set_tag_
                )
                constexpr mutable
                {
                    if (std::same_as<Other_, typename std::tuple_element_t<i_, Shapes_>::value_type::type>)
                    {
                        tag = true;
                    }
                    if constexpr (i_ < std::tuple_size_v<Shapes_> - 1)
                    {
                        set_tag_.template operator()<i_ + 1>(set_tag_);
                    }
                };
                set_tag(set_tag);
                return tag;
            }
        }

        template<ShapeConcept Other_>
        static constexpr
        ShapeCoordinates
        getShapeCoordinates()
        requires(hasShape<Other_>())
        {
            using Shapes_ = std::tuple_element_t<Shape_::getDimShape() - Other_::getDimShape() - 1, InnerNeighborhood<TypeView>>;
            auto tag = -1;
            auto set_tag = [&] <Integer i_ = 0> (
                auto & set_tag_
            )
            constexpr mutable
            {
                if (std::same_as<Other_, typename std::tuple_element_t<i_, Shapes_>::value_type::type>)
                {
                    tag = i_;
                }
                if constexpr (i_ < std::tuple_size_v<Shapes_> - 1)
                {
                    set_tag_.template operator()<i_ + 1>(set_tag_);
                }
            };
            set_tag(set_tag);
            return ShapeCoordinates(Shape_::getDimShape() - Other_::getDimShape() - 1, tag);
        }
        
        template<Integer... i_>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(sizeof...(i_) == 2)
        {
            auto constexpr _ = std::array<Integer, 2>{i_...};
            return std::tuple_size_v<std::tuple_element_t<_[1], std::tuple_element_t<_[0], InnerNeighborhood<TypeView>>>>;
        }
        
        template<Integer... i_>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(sizeof...(i_) == 1)
        {
            auto constexpr _ = std::array<Integer, 1>{i_...};
            return std::tuple_size_v<std::tuple_element_t<_[0], InnerNeighborhood<TypeView>>>;
        }
        
        template<Integer... i_>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(sizeof...(i_) == 0)
        {
            return std::tuple_size_v<InnerNeighborhood<TypeView>>;
        }

    };

    template<ShapeConcept Shape_, FrameConcept Frame_>
    struct ShapeOuterNeighborhoodTraits
    {

        template<template<typename, FrameConcept, typename...> typename T_, typename... U_>
        using OuterNeighborhood = tuple_slice_t<Shapes<T_, Frame_, U_...>, Shape_::getDimShape(), Frame_::getDimEuclidean() + 1>;

        template<Integer i_, Integer j_>
        using OuterNeighbor = std::tuple_element_t<j_, std::tuple_element_t<i_, OuterNeighborhood<TypeView>>>::type;

        template<ShapeConcept Other_>
        static constexpr
        Boolean
        hasShape()
        {
            if (Other_::getDimShape() < Shape_::getDimShape())
            {
                return false;
            }
            else
            {
                using Shapes_ = std::tuple_element_t<Other_::getDimShape() - Shape_::getDimShape(), OuterNeighborhood<TypeView>>;
                auto tag = false;
                auto set_tag = [&] <Integer i_ = 0> (
                    auto & set_tag_
                )
                constexpr mutable
                {
                    if (std::same_as<Other_, typename std::tuple_element_t<i_, Shapes_>::type>)
                    {
                        tag = true;
                    }
                    if constexpr (i_ < std::tuple_size_v<Shapes_> - 1)
                    {
                        set_tag_.template operator()<i_ + 1>(set_tag_);
                    }
                };
                set_tag(set_tag);
                return tag;
            }
        }

        template<ShapeConcept Other_>
        static constexpr
        ShapeCoordinates
        getShapeCoordinates()
        requires(hasShape<Other_>())
        {
            using Shapes_ = std::tuple_element_t<Other_::getDimShape() - Shape_::getDimShape(), OuterNeighborhood<TypeView>>;
            auto tag = -1;
            auto set_tag = [&] <Integer i_ = 0> (
                auto & set_tag_
            )
            constexpr mutable
            {
                if (std::same_as<Other_, typename std::tuple_element_t<i_, Shapes_>::type>)
                {
                    tag = i_;
                }
                if constexpr (i_ < std::tuple_size_v<Shapes_> - 1)
                {
                    set_tag_.template operator()<i_ + 1>(set_tag_);
                }
            };
            set_tag(set_tag);
            return ShapeCoordinates(Other_::getDimShape() - Shape_::getDimShape(), tag);
        }
        
        template<Integer... i_>
        static constexpr
        Integer
        getNumOuterNeighbors()
        requires(sizeof...(i_) == 1)
        {
            auto constexpr _ = std::array<Integer, 1>{i_...};
            return std::tuple_size_v<std::tuple_element_t<_[0], OuterNeighborhood<TypeView>>>;
        }
        
        template<Integer... i_>
        static constexpr
        Integer
        getNumOuterNeighbors()
        requires(sizeof...(i_) == 0)
        {
            return std::tuple_size_v<OuterNeighborhood<TypeView>>;
        }

    };

} // namespace lolita::geometry


#endif /* E68CB6A8_1D3D_4E96_928D_FB405436F896 */
