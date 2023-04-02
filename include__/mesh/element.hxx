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
    struct Element
    {

    private:

        template<typename... U_>
        using InnerNeighbor_ = std::shared_ptr<Element<U_...>>;

        using InnerNeighborhood_ = geometry::ShapeInnerNeighborhood<Shape_, InnerNeighbor_, Frame_, T_...>;

        using Region_ = std::shared_ptr<Region<typename Shape_::Domain, Frame_>>;

        using Regions_ = std::vector<Region_>;

    public:

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
            auto hash = std::basic_stringstream<Character>();
            for (auto const & node : getInnerNeighborhood().template getComponent<geometry::Node>())
            {
                hash << node->getHash();
            }
            return hash.str();
        }

        void
        addRegion(
            Region_ const & region
        )
        {
            regions_.push_back(region);
        }
        
        InnerNeighborhood_ const &
        getInnerNeighborhood()
        const
        {
            return inner_neighborhood_;
        }
        
        InnerNeighborhood_ &
        getInnerNeighborhood()
        {
            return inner_neighborhood_;
        }

    private:

        InnerNeighborhood_ inner_neighborhood_;
        
        Natural tag_;

        Regions_ regions_;

    };

    template<geometry::FrameConcept Frame_, typename... T_>
    struct Element<geometry::Node, Frame_, T_...>
    {

    private:

        template<typename... U_>
        using OuterNeighbors_ = std::vector<std::shared_ptr<Element<U_...>>>;

        using Shape_ = geometry::Node;

        using OuterNeighborhood_ = geometry::ShapeOuterNeighborhood<Shape_, Frame_, OuterNeighbors_, Frame_, T_...>;

        using Region_ = std::shared_ptr<Region<typename Shape_::Domain, Frame_>>;

        using Regions_ = std::vector<Region_>;

    public:

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
            return std::to_string(this->tag_);
        }

        void
        addRegion(
            Region_ const & region
        )
        {
            regions_.push_back(region);
        }

        void
        setCoordinates(
            geometry::Point<Frame_::getDimEuclidean()> const & point
        )
        {
            coordinates_ = point;
        }
        
        OuterNeighborhood_ const &
        getOuterNeighborhood()
        const
        {
            return outer_neighborhood_;
        }
        
        OuterNeighborhood_ &
        getOuterNeighborhood()
        {
            return outer_neighborhood_;
        }

    private:

        OuterNeighborhood_ outer_neighborhood_;
        
        Natural tag_;
        
        geometry::Point<Frame_::getDimEuclidean()> coordinates_;

        Regions_ regions_;
        
    };
    
} // namespace lolita::mesh


#endif /* AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A */
