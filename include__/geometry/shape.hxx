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
            
            /**
             * @brief 
             * 
             */
            using Domain = geometry::Domain<0>;
            
            /**
             * @brief Get the Dim Shape object
             * 
             * @return constexpr Integer 
             */
            static constexpr
            Integer
            getDimShape()
            {
                return 0;
            }
            
            /**
             * @brief Get the Num Nodes object
             * 
             * @return constexpr Integer 
             */
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
            
            /**
             * @brief 
             * 
             */
            using Domain = geometry::Domain<1>;

            /**
             * @brief 
             * 
             * @tparam T_ 
             * @tparam U_ 
             */
            template<template<typename, typename...> typename T_, typename... U_>
            using InnerNeighborhood = std::tuple<
                std::tuple<
                    std::array<T_<Point, U_...>, 2>
                >
            >;
            
            /**
             * @brief Get the Dim Shape object
             * 
             * @return constexpr Integer 
             */
            static constexpr
            Integer
            getDimShape()
            {
                return 1;
            }
            
            /**
             * @brief Get the Num Nodes object
             * 
             * @return constexpr Integer 
             */
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
            
            /**
             * @brief 
             * 
             */
            using Domain = geometry::Domain<2>;

            /**
             * @brief 
             * 
             * @tparam T_ 
             * @tparam U_ 
             */
            template<template<typename, typename...> typename T_, typename... U_>
            using InnerNeighborhood = std::tuple<
                std::tuple<
                    std::array<T_<Curve, U_...>, num_points_>
                >,
                std::tuple<
                    std::array<T_<Point, U_...>, num_points_>
                >
            >;
            
            /**
             * @brief Get the Dim Shape object
             * 
             * @return constexpr Integer 
             */
            static constexpr
            Integer
            getDimShape()
            {
                return 2;
            }
            
            /**
             * @brief Get the Num Nodes object
             * 
             * @return constexpr Integer 
             */
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

            /**
             * @brief 
             * 
             */
            using Domain = geometry::Domain<3>;

            /**
             * @brief 
             * 
             * @tparam T_ 
             * @tparam U_ 
             */
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
            
            /**
             * @brief Get the Dim Shape object
             * 
             * @return constexpr Integer 
             */
            static constexpr
            Integer
            getDimShape()
            {
                return 3;
            }
            
            /**
             * @brief Get the Num Nodes object
             * 
             * @return constexpr Integer 
             */
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

    namespace internal
    {
        
        /**
         * @brief 
         * 
         * @tparam Frame_ 
         * @tparam Shape_ 
         */
        template<FrameConcept Frame_, ShapeConcept Shape_>
        struct ShapeTopologyTraits
        {

            static constexpr
            Boolean is_cell_ = Frame_::getDimEuclidean() == Shape_::getDimShape();

            static constexpr
            Boolean is_face_ = Frame_::getDimEuclidean() == Shape_::getDimShape() + 1;
            
        };
        
    } // namespace internal
    
    /**
     * @brief 
     * 
     * @tparam F_ 
     * @tparam T_ 
     */
    template<typename F_, typename T_>
    concept CellConcept = internal::ShapeTopologyTraits<F_, T_>::is_cell_;
    
    /**
     * @brief 
     * 
     * @tparam F_ 
     * @tparam T_ 
     */
    template<typename F_, typename T_>
    concept FaceConcept = internal::ShapeTopologyTraits<F_, T_>::is_face_;

    /**
     * @brief 
     * 
     */
    using Node = shape::Point;

    /**
     * @brief 
     * 
     */
    using Segment = shape::Curve;

    /**
     * @brief 
     * 
     */
    using Triangle = shape::Facet<3>;

    /**
     * @brief 
     * 
     */
    using Quadrangle = shape::Facet<4>;

    /**
     * @brief 
     * 
     */
    using Tetrahedron = shape::Solid<4, 6, shape::internal::Faces(3, 4)>;

    /**
     * @brief 
     * 
     */
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
    
    template<FrameConcept Frame_>
    struct ShapeCollectionTraits
    {

        template<template<typename, typename...> typename T_ = TypeView, typename... U_>
        using Collection = tuple_slice_t<Shapes<T_, U_...>, 0, Frame_::getDimEuclidean() + 1>;
        
        template<Integer i_, Integer j_>
        using Shape = std::tuple_element_t<j_, std::tuple_element_t<i_, Collection<>>>;
        
        template<ShapeConcept Shape_>
        struct ShapeTraits
        {

        private:
            
            static constexpr
            Integer
            getIndex(
                Integer i
            )
            {
                auto constexpr first_index = Shape_::getDimShape();
                if (i == 0)
                {
                    return first_index;
                }
                else
                {
                    using Shapes_ = std::tuple_element_t<first_index, Collection<>>;
                    auto tag = -1;
                    auto set_tag = [&] <Integer i_ = 0> (
                        auto & set_tag_
                    )
                    constexpr mutable
                    {
                        if (std::same_as<Shape_, typename std::tuple_element_t<i_, Shapes_>>)
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
            }

        public:
        
            static constexpr
            Boolean
            hasCoordinates()
            {
                return getIndex(1) != -1;
            }
            
            static constexpr
            std::array<Integer, 2> coordinates_ = {getIndex(0), getIndex(1)};

        };
    
        template<ShapeConcept Shape_>
        static constexpr
        Integer
        getShapeIndex(
            Integer i
        )
        requires(ShapeTraits<Shape_>::hasCoordinates())
        {
            return ShapeTraits<Shape_>::coordinates_[i];
        }

        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<Collection<>>;
        }

        template<Integer i_>
        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<std::tuple_element_t<i_, Collection<>>>;
        }

        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<Shape<i_, j_>>();
                if constexpr (j_ < getNumComponents<i_>() - 1)
                {
                    apply_.template operator()<i_, j_ + 1>(apply_);
                }
                else if constexpr (i_ < getNumComponents() - 1)
                {
                    apply_.template operator()<i_ + 1, 0>(apply_);
                }
            };
            apply(apply);
        }

        template<Integer i_>
        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer j_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<Shape<i_, j_>>();
                if constexpr (j_ < getNumComponents<i_>() - 1)
                {
                    apply_.template operator()<j_ + 1>(apply_);
                }
            };
            apply(apply);
        }

    };

    template<ShapeConcept Shape_>
    struct ShapeInnerNeighborhoodTraits
    {

        template<template<typename, typename...> typename T_ = TypeView, typename... U_>
        using InnerNeighborhood = typename Shape_::template InnerNeighborhood<T_, U_...>;
        
        template<Integer i_, Integer j_>
        using Shape = std::tuple_element_t<j_, std::tuple_element_t<i_, InnerNeighborhood<>>>::value_type;
        
        template<ShapeConcept Neighbor_>
        struct ShapeTraits
        {

        private:
            
            static constexpr
            Integer
            getIndex(
                Integer i
            )
            {
                auto constexpr first_index = Shape_::getDimShape() - Neighbor_::getDimShape() - 1;
                if (i == 0)
                {
                    return first_index;
                }
                else
                {
                    using Shapes_ = std::tuple_element_t<first_index, InnerNeighborhood<>>;
                    auto tag = -1;
                    auto set_tag = [&] <Integer i_ = 0> (
                        auto & set_tag_
                    )
                    constexpr mutable
                    {
                        if (std::same_as<Neighbor_, typename std::tuple_element_t<i_, Shapes_>::value_type>)
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
            }

        public:
        
            static constexpr
            Boolean
            hasCoordinates()
            {
                return getIndex(1) != -1;
            }
            
            static constexpr
            std::array<Integer, 2> coordinates_ = {getIndex(0), getIndex(1)};

        };
    
        template<ShapeConcept Neighbor_>
        static constexpr
        Integer
        getShapeIndex(
            Integer i
        )
        requires(ShapeTraits<Neighbor_>::hasCoordinates())
        {
            return ShapeTraits<Neighbor_>::coordinates_[i];
        }

        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<InnerNeighborhood<>>;
        }

        template<Integer i_>
        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<std::tuple_element_t<i_, InnerNeighborhood<>>>;
        }

        template<ShapeConcept Neighbor_>
        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<std::tuple_element_t<getShapeIndex<Neighbor_>(1), std::tuple_element_t<getShapeIndex<Neighbor_>(0), InnerNeighborhood<>>>>;
        }

        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<Shape<i_, j_>>();
                if constexpr (j_ < getNumComponents<i_>() - 1)
                {
                    apply_.template operator()<i_, j_ + 1>(apply_);
                }
                else if constexpr (i_ < getNumComponents() - 1)
                {
                    apply_.template operator()<i_ + 1, 0>(apply_);
                }
            };
            apply(apply);
        }

        template<Integer i_>
        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer j_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<Shape<i_, j_>>();
                if constexpr (j_ < getNumComponents<i_>() - 1)
                {
                    apply_.template operator()<j_ + 1>(apply_);
                }
            };
            apply(apply);
        }

    };

    template<ShapeConcept Shape_, FrameConcept Frame_>
    struct ShapeOuterNeighborhoodTraits
    {

        template<template<typename, typename...> typename T_ = TypeView, typename... U_>
        using OuterNeighborhood = tuple_slice_t<Shapes<T_, U_...>, Shape_::getDimShape() + 1, Frame_::getDimEuclidean() + 1>;

        template<Integer i_, Integer j_>
        using Shape = std::tuple_element_t<j_, std::tuple_element_t<i_, OuterNeighborhood<>>>;
        
        template<ShapeConcept Neighbor_>
        struct ShapeTraits
        {

        private:
        
            static constexpr
            Integer
            getIndex(
                Integer i
            )
            {
                auto constexpr first_index = Neighbor_::getDimShape() - Shape_::getDimShape() - 1;
                if (i == 0)
                {
                    return first_index;
                }
                else
                {
                    using Shapes_ = std::tuple_element_t<first_index, OuterNeighborhood<>>;
                    auto tag = -1;
                    auto set_tag = [&] <Integer i_ = 0> (
                        auto & set_tag_
                    )
                    constexpr mutable
                    {
                        if (std::same_as<Neighbor_, typename std::tuple_element_t<i_, Shapes_>>)
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
            }

        public:
        
            static constexpr
            Boolean
            hasCoordinates()
            {
                return getIndex(1) != -1;
            }
            
            static constexpr
            std::array<Integer, 2> coordinates_ = {getIndex(0), getIndex(1)};

        };
    
        template<ShapeConcept Neighbor_>
        static constexpr
        Integer
        getShapeIndex(
            Integer i
        )
        requires(ShapeTraits<Neighbor_>::hasCoordinates())
        {
            return ShapeTraits<Neighbor_>::coordinates_[i];
        }
        
        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<OuterNeighborhood<>>;
        }

        template<Integer i_>
        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<std::tuple_element_t<i_, OuterNeighborhood<>>>;
        }

        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<Shape<i_, j_>>();
                if constexpr (j_ < getNumComponents<i_>() - 1)
                {
                    apply_.template operator()<i_, j_ + 1>(apply_);
                }
                else if constexpr (i_ < getNumComponents() - 1)
                {
                    apply_.template operator()<i_ + 1, 0>(apply_);
                }
            };
            apply(apply);
        }

        template<Integer i_>
        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer j_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<Shape<i_, j_>>();
                if constexpr (j_ < getNumComponents<i_>() - 1)
                {
                    apply_.template operator()<j_ + 1>(apply_);
                }
            };
            apply(apply);
        }

    };

} // namespace lolita::geometry


#endif /* E68CB6A8_1D3D_4E96_928D_FB405436F896 */
