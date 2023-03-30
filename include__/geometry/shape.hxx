/**
 * @file shape.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-03-29
 * 
 * @copyright Copyright (c) 2023
 * 
 */

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

    /**
     * @brief 
     * 
     * @tparam T_ 
     */
    template<typename T_>
    concept PointConcept = shape::internal::PointTraits<T_>::value;
    
    /**
     * @brief 
     * 
     * @tparam T_ 
     */
    template<typename T_>
    concept CurveConcept = shape::internal::CurveTraits<T_>::value;

    /**
     * @brief 
     * 
     * @tparam T_ 
     */
    template<typename T_>
    concept FacetConcept = shape::internal::FacetTraits<T_>::value;

    /**
     * @brief 
     * 
     * @tparam T_ 
     */
    template<typename T_>
    concept SolidConcept = shape::internal::SolidTraits<T_>::value;
    
    /**
     * @brief 
     * 
     * @tparam T_ 
     */
    template<typename T_>
    concept ShapeConcept = PointConcept<T_> || CurveConcept<T_> || FacetConcept<T_> || SolidConcept<T_>;

    using Node = shape::Point;
    using Segment = shape::Curve;
    using Triangle = shape::Facet<3>;
    using Quadrangle = shape::Facet<4>;
    using Tetrahedron = shape::Solid<4, 6, shape::internal::Faces(3, 4)>;
    using Hexahedron = shape::Solid<8, 12, shape::internal::Faces(4, 8)>;

    /**
     * @brief 
     * 
     * @tparam T_ 
     * @tparam U_ 
     */
    template<template<ShapeConcept, typename...> typename T_, typename... U_>
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

    /**
     * @brief 
     * 
     */
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
    
    template<FrameConcept Frame_>
    struct ShapeLibraryTraits
    {

        template<template<typename, typename...> typename T_, typename... U_>
        using MyShapes = tuple_slice_t<Shapes<T_, U_...>, 0, Frame_::getDimEuclidean() + 1>;
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         */
        template<Integer i_, Integer j_>
        using Shape = std::tuple_element_t<j_, std::tuple_element_t<i_, MyShapes<TypeView>>>::type;

        /**
         * @brief 
         * 
         * @tparam Shape_ 
         * @return constexpr Boolean 
         */
        template<ShapeConcept Shape_>
        static constexpr
        Boolean
        hasShape()
        {
            using Shapes_ = std::tuple_element_t<Shape_::getDimShape(), MyShapes<TypeView>>;
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

        /**
         * @brief 
         * 
         * @tparam Shape_ 
         */
        template<ShapeConcept Shape_>
        static constexpr
        ShapeCoordinates
        getShapeCoordinates()
        requires(hasShape<Shape_>())
        {
            using Shapes_ = std::tuple_element_t<Shape_::getDimShape(), MyShapes<TypeView>>;
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
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumElements()
        requires(sizeof...(t_i) == 0)
        {
            return std::tuple_size_v<MyShapes<TypeView>>;
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumElements()
        requires(sizeof...(t_i) == 1)
        {
            auto constexpr _ = std::array<Integer, 1>{t_i...};
            return std::tuple_size_v<std::tuple_element_t<_[0], MyShapes<TypeView>>>;
        }

    };

    /**
     * @brief 
     * 
     * @tparam Shape_ 
     */
    template<ShapeConcept Shape_>
    requires(!PointConcept<Shape_>)
    struct ShapeInnerNeighborhoodTraits
    {
        
        /**
         * @brief 
         * 
         * @tparam T_ 
         * @tparam U_ 
         */
        template<template<typename, typename...> typename T_, typename... U_>
        using InnerNeighborhood = typename Shape_::template InnerNeighborhood<T_, U_...>;

        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         */
        template<Integer i_, Integer j_>
        using InnerNeighbor = std::tuple_element_t<j_, std::tuple_element_t<i_, InnerNeighborhood<TypeView>>>::value_type::type;

        /**
         * @brief 
         * 
         * @tparam Other_ 
         * @return constexpr Boolean 
         */
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

        /**
         * @brief 
         * 
         * @tparam Other_ 
         */
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
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         */
        template<Integer... i_>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(sizeof...(i_) == 2)
        {
            auto constexpr _ = std::array<Integer, 2>{i_...};
            return std::tuple_size_v<std::tuple_element_t<_[1], std::tuple_element_t<_[0], InnerNeighborhood<TypeView>>>>;
        }
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         */
        template<Integer... i_>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(sizeof...(i_) == 1)
        {
            auto constexpr _ = std::array<Integer, 1>{i_...};
            return std::tuple_size_v<std::tuple_element_t<_[0], InnerNeighborhood<TypeView>>>;
        }
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         */
        template<Integer... i_>
        static constexpr
        Integer
        getNumInnerNeighbors()
        requires(sizeof...(i_) == 0)
        {
            return std::tuple_size_v<InnerNeighborhood<TypeView>>;
        }

    };

    /**
     * @brief 
     * 
     * @tparam Shape_ 
     * @tparam Frame_ 
     */
    template<ShapeConcept Shape_, FrameConcept Frame_>
    struct ShapeOuterNeighborhoodTraits
    {

        /**
         * @brief 
         * 
         * @tparam T_ 
         * @tparam U_ 
         */
        template<template<typename, typename...> typename T_, typename... U_>
        using OuterNeighborhood = tuple_slice_t<Shapes<T_, U_...>, Shape_::getDimShape(), Frame_::getDimEuclidean() + 1>;

        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         */
        template<Integer i_, Integer j_>
        using OuterNeighbor = std::tuple_element_t<j_, std::tuple_element_t<i_, OuterNeighborhood<TypeView>>>::type;

        /**
         * @brief 
         * 
         * @tparam Other_ 
         * @return constexpr Boolean 
         */
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

        /**
         * @brief 
         * 
         * @tparam Other_ 
         */
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
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         */
        template<Integer... i_>
        static constexpr
        Integer
        getNumOuterNeighbors()
        requires(sizeof...(i_) == 1)
        {
            auto constexpr _ = std::array<Integer, 1>{i_...};
            return std::tuple_size_v<std::tuple_element_t<_[0], OuterNeighborhood<TypeView>>>;
        }
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         */
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
