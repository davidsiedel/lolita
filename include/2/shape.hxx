#ifndef AF1D9E47_C7A3_4252_B958_B6C481143805
#define AF1D9E47_C7A3_4252_B958_B6C481143805

#include "config.hxx"

namespace lolita
{

    template<typename Implementation_>
    struct ShapeBase;

    namespace detail
    {
        
        struct ShapeType
        {

        private:

            constexpr
            ShapeType()
            {}

            constexpr
            Boolean
            operator==(
                ShapeType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                ShapeType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::ShapeBase;

        };
        
    } // namespace detail

    template<typename T>
    concept ShapeConcept = std::derived_from<std::decay_t<T>, detail::ShapeType>;

    template<typename Implementation_>
    struct ShapeBase : detail::ShapeType
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return Implementation_::getLabel();
        }

    protected:

        constexpr
        ShapeBase(
            Integer dim,
            Integer num_points,
            Integer num_curves,
            Integer num_facets
        )
        :
        detail::ShapeType(),
        dim_(dim),
        num_points_(num_points),
        num_curves_(num_curves),
        num_facets_(num_facets)
        {}

    public:

        constexpr
        Boolean
        operator==(
            ShapeBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            ShapeBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            ShapeBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            ShapeBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }
        
        constexpr
        Boolean
        hasDim(
            Integer dim
        )
        const
        {
            return dim_ == dim;
        }
        
        constexpr
        Boolean
        isSub(
            auto const & frame,
            Integer level
        )
        const
        {
            return frame.getDim() - dim_ == level;
        }

        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }

        constexpr
        Integer
        getNumPoints()
        const
        {
            return num_points_;
        }

        constexpr
        Integer
        getNumCurves()
        const
        {
            return num_curves_;
        }

        constexpr
        Integer
        getNumFacets()
        const
        {
            return num_facets_;
        }

        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            ShapeBase const & shape
        )
        {
            os << "{ " << ShapeBase::getLabel() << " | ";
            os << shape.num_points_ << " points | ";
            os << shape.num_curves_ << " curves | ";
            os << shape.num_facets_ << " facets }";
            return os;
        }

        Integer dim_;

        Integer num_points_;

        Integer num_curves_;
        
        Integer num_facets_;

    };

    struct PointShape : ShapeBase<PointShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "Point";
        }

    private:

        using Base_ = ShapeBase<PointShape>;

    protected:

        constexpr
        PointShape()
        :
        Base_(0, 1, 0, 0)
        {}

    };

    template<typename T>
    concept PointShapeConcept = std::derived_from<std::decay_t<T>, PointShape>;
    
    struct SegmentShape : ShapeBase<SegmentShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "Segment";
        }

    private:

        using Base_ = ShapeBase<SegmentShape>;

    protected:

        constexpr
        SegmentShape()
        :
        Base_(1, 2, 0, 0)
        {}

    };

    template<typename T>
    concept SegmentShapeConcept = std::derived_from<std::decay_t<T>, SegmentShape>;
    
    struct TriangleShape : ShapeBase<TriangleShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "Triangle";
        }

    private:

        using Base_ = ShapeBase<TriangleShape>;

    protected:

        constexpr
        TriangleShape()
        :
        Base_(2, 3, 3, 0)
        {}

    };

    template<typename T>
    concept TriangleShapeConcept = std::derived_from<std::decay_t<T>, TriangleShape>;
    
    struct QuadrangleShape : ShapeBase<QuadrangleShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "Quadrangle";
        }

    private:

        using Base_ = ShapeBase<QuadrangleShape>;

    protected:

        constexpr
        QuadrangleShape()
        :
        Base_(2, 4, 4, 0)
        {}

    };

    template<typename T>
    concept QuadrangleShapeConcept = std::derived_from<std::decay_t<T>, QuadrangleShape>;
    
    struct TetrahedronShape : ShapeBase<TetrahedronShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "Tetrahedron";
        }

    private:

        using Base_ = ShapeBase<TetrahedronShape>;

    protected:

        constexpr
        TetrahedronShape()
        :
        Base_(3, 4, 6, 4)
        {}

    };
    
    template<typename T>
    concept TetrahedronShapeConcept = std::derived_from<std::decay_t<T>, TetrahedronShape>;
    
    template<typename Implementation_, typename Shape_>
    struct LagrangeShapeBase;

    namespace detail
    {
        
        struct LagrangeShapeType
        {

        private:

            constexpr
            LagrangeShapeType()
            {}

            constexpr
            Boolean
            operator==(
                LagrangeShapeType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                LagrangeShapeType const & other
            )
            const
            = default;

            template<typename Implementation_, typename Shape_>
            friend class lolita::LagrangeShapeBase;

        };
        
    } // namespace detail

    template<typename T>
    concept LagrangeShapeConcept = std::derived_from<std::decay_t<T>, detail::LagrangeShapeType>;
    
    template<typename Implementation_, typename Shape_>
    struct LagrangeShapeBase : detail::LagrangeShapeType, Shape_
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return Implementation_::getLabel();
        }

        using Shape = Shape_;

    protected:

        constexpr
        LagrangeShapeBase(
            Integer order,
            Integer num_nodes
        )
        :
        detail::LagrangeShapeType(),
        Shape(),
        order_(order),
        num_nodes_(num_nodes)
        {}

    public:

        constexpr
        Boolean
        operator==(
            LagrangeShapeBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            LagrangeShapeBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            LagrangeShapeBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        };

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            LagrangeShapeBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        };

        constexpr
        Integer
        getOrder()
        const
        {
            return order_;
        }

        constexpr
        Integer
        getNumNodes()
        const
        {
            return num_nodes_;
        }
        
        constexpr
        Boolean
        hasOrd(
            Integer ord
        )
        const
        {
            return order_ == ord;
        }

        Integer order_;

        Integer num_nodes_;

    };

    struct ShapeCoordinates
    {

        constexpr
        ShapeCoordinates()
        :
        dim_(-1),
        tag_(-1)
        {}
        
        constexpr
        ShapeCoordinates(
            Integer dim,
            Integer tag
        )
        :
        dim_(dim),
        tag_(tag)
        {}
        
        constexpr
        void
        setDim(
            Integer dim
        )
        {
            dim_ = dim;
        }
        
        constexpr
        void
        setTag(
            Integer tag
        )
        {
            tag_ = tag;
        }
        
        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }
        
        constexpr
        Integer
        getTag()
        const
        {
            return tag_;
        }
        
        Integer dim_;
        
        Integer tag_;

    };

    struct Node : LagrangeShapeBase<Node, PointShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "Node";
        }
        
        static constexpr
        Integer
        getMshTag()
        {
            return 15;
        }

    private:

        using Base_ = LagrangeShapeBase<Node, PointShape>;

    public:

        constexpr
        Node()
        :
        Base_(0, 1)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsNodeTraits : std::false_type
        {};

        template<>
        struct IsNodeTraits<Node> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept NodeConcept = detail::IsNodeTraits<std::decay_t<T>>::value;
    
    struct LinearSegment : LagrangeShapeBase<LinearSegment, SegmentShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "LinearSegment";
        }
        
        static constexpr
        Integer
        getMshTag()
        {
            return 1;
        }

    private:

        using Base_ = LagrangeShapeBase<LinearSegment, SegmentShape>;

    public:

        constexpr
        LinearSegment()
        :
        Base_(1, 2)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsLinearSegmentTraits : std::false_type
        {};

        template<>
        struct IsLinearSegmentTraits<LinearSegment> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept LinearSegmentConcept = detail::IsLinearSegmentTraits<std::decay_t<T>>::value;
    
    struct LinearTriangle : LagrangeShapeBase<LinearTriangle, TriangleShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "LinearTriangle";
        }
        
        static constexpr
        Integer
        getMshTag()
        {
            return 2;
        }

    private:

        using Base_ = LagrangeShapeBase<LinearTriangle, TriangleShape>;

    public:

        constexpr
        LinearTriangle()
        :
        Base_(1, 3)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsLinearTriangleTraits : std::false_type
        {};

        template<>
        struct IsLinearTriangleTraits<LinearTriangle> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept LinearTriangleConcept = detail::IsLinearTriangleTraits<std::decay_t<T>>::value;

    struct LinearQuadrangle : LagrangeShapeBase<LinearQuadrangle, QuadrangleShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "LinearQuadrangle";
        }
        
        static constexpr
        Integer
        getMshTag()
        {
            return 3;
        }

    private:

        using Base_ = LagrangeShapeBase<LinearQuadrangle, QuadrangleShape>;

    public:

        constexpr
        LinearQuadrangle()
        :
        Base_(1, 4)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsLinearQuadrangleTraits : std::false_type
        {};

        template<>
        struct IsLinearQuadrangleTraits<LinearQuadrangle> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept LinearQuadrangleConcept = detail::IsLinearQuadrangleTraits<std::decay_t<T>>::value;

    struct LinearTetrahedron : LagrangeShapeBase<LinearTetrahedron, TetrahedronShape>
    {
        
        static constexpr
        std::basic_string_view<Character>
        getLabel()
        {
            return "LinearTetrahedron";
        }
        
        static constexpr
        Integer
        getMshTag()
        {
            return 4;
        }

    private:

        using Base_ = LagrangeShapeBase<LinearTetrahedron, TetrahedronShape>;

    public:

        constexpr
        LinearTetrahedron()
        :
        Base_(1, 4)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsLinearTetrahedronTraits : std::false_type
        {};

        template<>
        struct IsLinearTetrahedronTraits<LinearTetrahedron> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept LinearTetrahedronConcept = detail::IsLinearTetrahedronTraits<std::decay_t<T>>::value;

    struct ShapeLibrary
    {

    private:

        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using Points_ = std::tuple<
            T<Node{}, t_args...>
        >;
        
        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using Curves_ = std::tuple<
            T<LinearSegment{}, t_args...>
        >;
        
        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using Facets_ = std::tuple<
            T<LinearTriangle{}, t_args...>,
            T<LinearQuadrangle{}, t_args...>
        >;
        
        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using Solids_ = std::tuple<
            T<LinearTetrahedron{}, t_args...>
        >;

    public:

        template<template<LagrangeShapeConcept auto, auto...> typename T, auto... t_args>
        using Elements = std::tuple<
            Points_<T, t_args...>,
            Curves_<T, t_args...>,
            Facets_<T, t_args...>,
            Solids_<T, t_args...>
        >;

    };

} // namespace lolita


#endif /* AF1D9E47_C7A3_4252_B958_B6C481143805 */
