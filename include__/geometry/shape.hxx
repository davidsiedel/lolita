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

    /**
     * @brief 
     * 
     * @tparam Frame_ 
     * @tparam T_ 
     * @tparam U_ 
     */
    template<FrameConcept Frame_, template<typename, typename...> typename T_ = TypeView2, typename... U_>
    struct ShapeCollection
    {
        
        /**
         * @brief 
         * 
         */
        using Components = tuple_slice_t<Shapes<T_, U_...>, 0, Frame_::getDimEuclidean() + 1>;
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         */
        template<Integer i_, Integer j_>
        using Component = std::tuple_element_t<j_, std::tuple_element_t<i_, Components>>;

        /**
         * @brief 
         * 
         * @tparam Shape_ 
         */
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
                    using Shapes_ = std::tuple_element_t<first_index, typename ShapeCollection<Frame_>::Components>;
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
            
            /**
             * @brief 
             * 
             * @return constexpr Boolean 
             */
            static constexpr
            Boolean
            hasCoordinates()
            {
                return getIndex(1) != -1;
            }

            /**
             * @brief 
             * 
             */
            static constexpr
            std::array<Integer, 2> coordinates_ = {getIndex(0), getIndex(1)};

        };

        // template<ShapeConcept Shape_>
        // static constexpr
        // std::array<Integer, 2>
        // getShapeCoordinates()
        // {
        //     return ShapeTraits<Shape_>::coordinates_;
        // }

        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<Components>;
        }

        static constexpr
        Integer
        getNumComponents(
            Integer i
        )
        {
            auto size = -1;
            auto set_size = [&] <Integer i_ = 0> (
                auto & set_size_
            )
            constexpr mutable
            {
                if (i_ == i)
                {
                    size = std::tuple_size_v<std::tuple_element_t<i_, Components>>;
                }
                if constexpr (i_ < std::tuple_size_v<Components> - 1)
                {
                    set_size_.template operator()<i_ + 1>(set_size_);
                }
            };
            set_size(set_size);
            return size;
        }
    
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> const &
        getComponent()
        const
        {
            return std::get<j_>(std::get<i_>(components_));
        }
    
        template<ShapeConcept Shape_>
        std::tuple_element_t<ShapeTraits<Shape_>::coordinates_[1], std::tuple_element_t<ShapeTraits<Shape_>::coordinates_[0], Components>> const &
        getComponent()
        const
        requires(ShapeTraits<Shape_>::hasCoordinates())
        {
            return std::get<ShapeTraits<Shape_>::coordinates_[1]>(std::get<ShapeTraits<Shape_>::coordinates_[0]>(components_));
        }
        
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> &
        getComponent()
        {
            return std::get<j_>(std::get<i_>(components_));
        }
    
        template<ShapeConcept Shape_>
        std::tuple_element_t<ShapeTraits<Shape_>::coordinates_[1], std::tuple_element_t<ShapeTraits<Shape_>::coordinates_[0], Components>> &
        getComponent()
        requires(ShapeTraits<Shape_>::hasCoordinates())
        {
            return std::get<ShapeTraits<Shape_>::coordinates_[1]>(std::get<ShapeTraits<Shape_>::coordinates_[0]>(components_));
        }

    private:

        Components components_;

    };

    /**
     * @brief 
     * 
     * @tparam Shape_ 
     * @tparam T_ 
     * @tparam U_ 
     */
    template<ShapeConcept Shape_, template<typename, typename...> typename T_ = TypeView2, typename... U_>
    requires(!PointConcept<Shape_>)
    struct ShapeInnerNeighborhood
    {
        
        /**
         * @brief 
         * 
         */
        using Components = typename Shape_::template InnerNeighborhood<T_, U_...>;
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         */
        template<Integer i_, Integer j_>
        using Component = std::tuple_element_t<j_, std::tuple_element_t<i_, Components>>::value_type;

        /**
         * @brief 
         * 
         * @tparam Neighbor_ 
         */
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
                    using Shapes_ = std::tuple_element_t<first_index, typename ShapeInnerNeighborhood<Shape_>::Components>;
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
            
            /**
             * @brief 
             * 
             * @return constexpr Boolean 
             */
            static constexpr
            Boolean
            hasCoordinates()
            {
                return getIndex(1) != -1;
            }

            /**
             * @brief 
             * 
             */
            static constexpr
            std::array<Integer, 2> coordinates_ = {getIndex(0), getIndex(1)};

        };

        // template<ShapeConcept Neighbor_>
        // static constexpr
        // Boolean
        // hasNeighbor()
        // {
        //     return ShapeTraits<Neighbor_>::hasCoordinates();
        // }

        // template<ShapeConcept Neighbor_>
        // static constexpr
        // std::array<Integer, 2>
        // getNeighborCoordinates()
        // requires(hasNeighbor<Neighbor_>())
        // {
        //     return ShapeTraits<Neighbor_>::coordinates_;
        // }

        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<Components>;
        }

        static constexpr
        Integer
        getNumComponents(
            Integer i
        )
        {
            auto size = -1;
            auto set_size = [&] <Integer i_ = 0> (
                auto & set_size_
            )
            constexpr mutable
            {
                if (i_ == i)
                {
                    size = std::tuple_size_v<std::tuple_element_t<i_, Components>>;
                }
                if constexpr (i_ < std::tuple_size_v<Components> - 1)
                {
                    set_size_.template operator()<i_ + 1>(set_size_);
                }
            };
            set_size(set_size);
            return size;
        }

        static constexpr
        Integer
        getNumComponents(
            Integer i,
            Integer j
        )
        {
            auto size = -1;
            auto set_size = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & set_size_
            )
            constexpr mutable
            {
                if (i_ == i && j_ == j)
                {
                    size = std::tuple_size_v<std::tuple_element_t<j_, std::tuple_element_t<i_, Components>>>;
                }
                if constexpr (j_ < std::tuple_size_v<std::tuple_element_t<i_, Components>> - 1)
                {
                    set_size_.template operator()<i_, j_ + 1>(set_size_);
                }
                else if constexpr (i_ < std::tuple_size_v<Components> - 1)
                {
                    set_size_.template operator()<i_ + 1, 0>(set_size_);
                }
            };
            set_size(set_size);
            return size;
        }
    
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> const &
        getComponent()
        const
        {
            return std::get<j_>(std::get<i_>(components_));
        }
    
        template<ShapeConcept Neighbor_>
        std::tuple_element_t<ShapeTraits<Neighbor_>::coordinates_[1], std::tuple_element_t<ShapeTraits<Neighbor_>::coordinates_[0], Components>> const &
        getComponent()
        const
        requires(ShapeTraits<Neighbor_>::hasCoordinates())
        {
            return std::get<ShapeTraits<Neighbor_>::coordinates_[1]>(std::get<ShapeTraits<Neighbor_>::coordinates_[0]>(components_));
        }
        
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> &
        getComponent()
        {
            return std::get<j_>(std::get<i_>(components_));
        }
    
        template<ShapeConcept Neighbor_>
        std::tuple_element_t<ShapeTraits<Neighbor_>::coordinates_[1], std::tuple_element_t<ShapeTraits<Neighbor_>::coordinates_[0], Components>> &
        getComponent()
        requires(ShapeTraits<Neighbor_>::hasCoordinates())
        {
            return std::get<ShapeTraits<Neighbor_>::coordinates_[1]>(std::get<ShapeTraits<Neighbor_>::coordinates_[0]>(components_));
        }

    private:

        Components components_;

    };

    /**
     * @brief 
     * 
     * @tparam Shape_ 
     * @tparam Frame_ 
     * @tparam T_ 
     * @tparam U_ 
     */
    template<ShapeConcept Shape_, FrameConcept Frame_, template<typename, typename...> typename T_ = TypeView2, typename... U_>
    struct ShapeOuterNeighborhood
    {
        
        /**
         * @brief 
         * 
         */
        using Components = tuple_slice_t<Shapes<T_, U_...>, Shape_::getDimShape(), Frame_::getDimEuclidean() + 1>;;
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         */
        template<Integer i_, Integer j_>
        using Component = std::tuple_element_t<j_, std::tuple_element_t<i_, Components>>;

        /**
         * @brief 
         * 
         * @tparam Neighbor_ 
         */
        template<ShapeConcept Neighbor_>
        struct ShapeTraits
        {

        private:
            
            /**
             * @brief Get the Index object
             * 
             * @param i 
             * @return constexpr Integer 
             */
            static constexpr
            Integer
            getIndex(
                Integer i
            )
            {
                auto constexpr first_index = Neighbor_::getDimShape() - Shape_::getDimShape();
                if (i == 0)
                {
                    return first_index;
                }
                else
                {
                    using Shapes_ = std::tuple_element_t<first_index, typename ShapeOuterNeighborhood<Shape_, Frame_>::Components>;
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
            
            /**
             * @brief 
             * 
             * @return constexpr Boolean 
             */
            static constexpr
            Boolean
            hasCoordinates()
            {
                return getIndex(1) != -1;
            }

            /**
             * @brief 
             * 
             */
            static constexpr
            std::array<Integer, 2> coordinates_ = {getIndex(0), getIndex(1)};

        };

        // template<ShapeConcept Neighbor_>
        // static constexpr
        // std::array<Integer, 2>
        // getShapeCoordinates()
        // {
        //     return ShapeTraits<Neighbor_>::coordinates_;
        // }

        /**
         * @brief Get the Num Components object
         * 
         * @return constexpr Integer 
         */
        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<Components>;
        }

        /**
         * @brief Get the Num Components object
         * 
         * @param i 
         * @return constexpr Integer 
         */
        static constexpr
        Integer
        getNumComponents(
            Integer i
        )
        {
            auto size = -1;
            auto set_size = [&] <Integer i_ = 0> (
                auto & set_size_
            )
            constexpr mutable
            {
                if (i_ == i)
                {
                    size = std::tuple_size_v<std::tuple_element_t<i_, Components>>;
                }
                if constexpr (i_ < std::tuple_size_v<Components> - 1)
                {
                    set_size_.template operator()<i_ + 1>(set_size_);
                }
            };
            set_size(set_size);
            return size;
        }

        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         * @return std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> const& 
         */
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> const &
        getComponent()
        const
        {
            return std::get<j_>(std::get<i_>(components_));
        }
        
        /**
         * @brief 
         * 
         * @tparam Neighbor_ 
         * @return std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> const& 
         */
        template<ShapeConcept Neighbor_>
        std::tuple_element_t<ShapeTraits<Neighbor_>::coordinates_[1], std::tuple_element_t<ShapeTraits<Neighbor_>::coordinates_[0], Components>> const &
        getComponent()
        const
        requires(ShapeTraits<Neighbor_>::hasCoordinates())
        {
            return std::get<ShapeTraits<Neighbor_>::coordinates_[1]>(std::get<ShapeTraits<Neighbor_>::coordinates_[0]>(components_));
        }
        
        /**
         * @brief 
         * 
         * @tparam i_ 
         * @tparam j_ 
         * @return std::tuple_element_t<j_, std::tuple_element_t<i_, Components>>& 
         */
        template<Integer i_, Integer j_>
        std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> &
        getComponent()
        {
            return std::get<j_>(std::get<i_>(components_));
        }
        
        /**
         * @brief 
         * 
         * @tparam Neighbor_ 
         * @return std::tuple_element_t<j_, std::tuple_element_t<i_, Components>> const& 
         */
        template<ShapeConcept Neighbor_>
        std::tuple_element_t<ShapeTraits<Neighbor_>::coordinates_[1], std::tuple_element_t<ShapeTraits<Neighbor_>::coordinates_[0], Components>> &
        getComponent()
        requires(ShapeTraits<Neighbor_>::hasCoordinates())
        {
            return std::get<ShapeTraits<Neighbor_>::coordinates_[1]>(std::get<ShapeTraits<Neighbor_>::coordinates_[0]>(components_));
        }

    private:

        /**
         * @brief 
         * 
         */
        Components components_;

    };

} // namespace lolita::geometry


#endif /* E68CB6A8_1D3D_4E96_928D_FB405436F896 */
