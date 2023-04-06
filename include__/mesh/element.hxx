/**
 * @file element.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-03-29
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A
#define AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A

#include "geometry/frame.hxx"
#include "geometry/domain.hxx"
#include "geometry/point.hxx"
#include "geometry/shape.hxx"
#include "mesh/region.hxx"

namespace lolita::mesh
{
    
    template<geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_, typename... T_>
    struct Element;
    
    template<geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_, typename... T_>
    struct ElementConnectivity
    {

    private:

        template<typename... U_>
        using InnerNeighbor_ = std::shared_ptr<Element<U_...>>;

        template<typename... U_>
        using OuterNeighbors_ = std::vector<std::shared_ptr<Element<U_...>>>;

        using InnerNeighborhoodTraits_ = geometry::ShapeInnerNeighborhoodTraits<Shape_>;

        using OuterNeighborhoodTraits_ = geometry::ShapeOuterNeighborhoodTraits<Shape_, Frame_>;

    public:

        using InnerNeighborhood = typename InnerNeighborhoodTraits_::template InnerNeighborhood<InnerNeighbor_, Frame_, T_...>;

        using OuterNeighborhood = typename OuterNeighborhoodTraits_::template OuterNeighborhood<OuterNeighbors_, Frame_, T_...>;

        template<geometry::ShapeConcept Neighbor_>
        using InnerNeighbors = std::array<InnerNeighbor_<Neighbor_, Frame_, T_...>, InnerNeighborhoodTraits_::template getNumComponents<Neighbor_>()>;

        template<geometry::ShapeConcept Neighbor_>
        using OuterNeighbors = OuterNeighbors_<Neighbor_, Frame_, T_...>;

        template<geometry::ShapeConcept Neighbor_>
        InnerNeighbors<Neighbor_> const &
        getInnerNeighbors()
        const
        {
            auto constexpr i = InnerNeighborhoodTraits_::template getShapeIndex<Neighbor_>(0);
            auto constexpr j = InnerNeighborhoodTraits_::template getShapeIndex<Neighbor_>(1);
            return std::get<j>(std::get<i>(inner_neighborhood_));
        }

        template<geometry::ShapeConcept Neighbor_>
        InnerNeighbors<Neighbor_> &
        getInnerNeighbors()
        {
            auto constexpr i = InnerNeighborhoodTraits_::template getShapeIndex<Neighbor_>(0);
            auto constexpr j = InnerNeighborhoodTraits_::template getShapeIndex<Neighbor_>(1);
            return std::get<j>(std::get<i>(inner_neighborhood_));
        }

        template<geometry::ShapeConcept Neighbor_>
        OuterNeighbors<Neighbor_> const &
        getOuterNeighbors()
        const
        {
            auto constexpr i = OuterNeighborhoodTraits_::template getShapeIndex<Neighbor_>(0);
            auto constexpr j = OuterNeighborhoodTraits_::template getShapeIndex<Neighbor_>(1);
            return std::get<j>(std::get<i>(outer_neighborhood_));
        }

        template<geometry::ShapeConcept Neighbor_>
        OuterNeighbors<Neighbor_> &
        getOuterNeighbors()
        {
            auto constexpr i = OuterNeighborhoodTraits_::template getShapeIndex<Neighbor_>(0);
            auto constexpr j = OuterNeighborhoodTraits_::template getShapeIndex<Neighbor_>(1);
            return std::get<j>(std::get<i>(outer_neighborhood_));
        }

    private:

        OuterNeighborhood outer_neighborhood_;

        InnerNeighborhood inner_neighborhood_;

    };
    
    template<geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_, typename... T_>
    requires(geometry::PointConcept<Shape_>)
    struct ElementConnectivity<Shape_, Frame_, T_...>
    {

    private:

        template<typename... U_>
        using OuterNeighbors_ = std::vector<std::shared_ptr<Element<U_...>>>;

        using OuterNeighborhoodTraits_ = geometry::ShapeOuterNeighborhoodTraits<Shape_, Frame_>;

    public:

        using OuterNeighborhood = typename OuterNeighborhoodTraits_::template OuterNeighborhood<OuterNeighbors_, Frame_, T_...>;

        template<geometry::ShapeConcept Neighbor_>
        using OuterNeighbors = OuterNeighbors_<Neighbor_, Frame_, T_...>;

        template<geometry::ShapeConcept Neighbor_>
        OuterNeighbors<Neighbor_> const &
        getOuterNeighbors()
        const
        {
            auto constexpr i = OuterNeighborhoodTraits_::template getShapeIndex<Neighbor_>(0);
            auto constexpr j= OuterNeighborhoodTraits_::template getShapeIndex<Neighbor_>(1);
            return std::get<j>(std::get<i>(outer_neighborhood_));
        }

        template<geometry::ShapeConcept Neighbor_>
        OuterNeighbors<Neighbor_> &
        getOuterNeighbors()
        {
            auto constexpr i = OuterNeighborhoodTraits_::template getShapeIndex<Neighbor_>(0);
            auto constexpr j = OuterNeighborhoodTraits_::template getShapeIndex<Neighbor_>(1);
            return std::get<j>(std::get<i>(outer_neighborhood_));
        }

        geometry::Point<Frame_::getDimEuclidean()> const &
        getCoordinates()
        const
        {
            return coordinates_;
        }

        geometry::Point<Frame_::getDimEuclidean()> &
        getCoordinates()
        {
            return coordinates_;
        }

        void
        setCoordinates(
            geometry::Point<Frame_::getDimEuclidean()> const & point
        )
        {
            coordinates_ = point;
        }

    private:
        
        geometry::Point<Frame_::getDimEuclidean()> coordinates_;

        OuterNeighborhood outer_neighborhood_;

    };
    
    template<geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_, typename... T_>
    requires(geometry::CellConcept<Frame_, Shape_>)
    struct ElementConnectivity<Shape_, Frame_, T_...>
    {

    private:

        template<typename... U_>
        using InnerNeighbor_ = std::shared_ptr<Element<U_...>>;

        using InnerNeighborhoodTraits_ = geometry::ShapeInnerNeighborhoodTraits<Shape_>;

    public:

        using InnerNeighborhood = typename InnerNeighborhoodTraits_::template InnerNeighborhood<InnerNeighbor_, Frame_, T_...>;

        template<geometry::ShapeConcept Neighbor_>
        using InnerNeighbors = std::array<InnerNeighbor_<Neighbor_, Frame_, T_...>, InnerNeighborhoodTraits_::template getNumComponents<Neighbor_>()>;

        template<geometry::ShapeConcept Neighbor_>
        InnerNeighbors<Neighbor_> const &
        getInnerNeighbors()
        const
        {
            auto constexpr i = InnerNeighborhoodTraits_::template getShapeIndex<Neighbor_>(0);
            auto constexpr j = InnerNeighborhoodTraits_::template getShapeIndex<Neighbor_>(1);
            return std::get<j>(std::get<i>(inner_neighborhood_));
        }

        template<geometry::ShapeConcept Neighbor_>
        InnerNeighbors<Neighbor_> &
        getInnerNeighbors()
        {
            auto constexpr i = InnerNeighborhoodTraits_::template getShapeIndex<Neighbor_>(0);
            auto constexpr j = InnerNeighborhoodTraits_::template getShapeIndex<Neighbor_>(1);
            return std::get<j>(std::get<i>(inner_neighborhood_));
        }

    private:

        InnerNeighborhood inner_neighborhood_;

    };
    
    template<geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_, typename... T_>
    struct Element
    {

        using Connectivity_ = ElementConnectivity<Shape_, Frame_, T_...>;

        using Region_ = std::shared_ptr<Region<typename Shape_::Domain, Frame_>>;

        using Regions_ = std::vector<Region_>;

        explicit
        Element(
            Natural const & tag
        )
        :
        tag_(tag)
        {}

        Natural const &
        getTag()
        const
        {
            return tag_;
        }
        
        std::basic_string<Character>
        getHash()
        const
        {
            if constexpr (geometry::PointConcept<Shape_>)
            {
                return std::to_string(tag_);
            }
            else
            {
                auto hash = std::basic_stringstream<Character>();
                for (auto const & node : getConnectivity().template getInnerNeighbors<geometry::Node>())
                {
                    hash << node->getHash();
                }
                return hash.str();
            }
        }

        void
        addRegion(
            Region_ const & region
        )
        {
            regions_.push_back(region);
        }

        Connectivity_ const &
        getConnectivity()
        const
        {
            return connectivity_;
        }

        Connectivity_ &
        getConnectivity()
        {
            return connectivity_;
        }

    private:

        Connectivity_ connectivity_;
        
        Natural tag_;

        Regions_ regions_;

    };
    
} // namespace lolita::mesh


#endif /* AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A */
