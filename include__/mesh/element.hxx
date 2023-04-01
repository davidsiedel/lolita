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
    
    template<geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_>
    struct Element
    {

    private:

        template<typename... U_>
        using ElementPointer_ = std::shared_ptr<Element<U_...>>;

        template<typename... U_>
        using ElementPointers_ = std::vector<std::shared_ptr<Element<U_...>>>;

        using InnerNeighborhood_ = geometry::ShapeInnerNeighborhood<Shape_, ElementPointer_, Frame_>;

        using OuterNeighborhood_ = geometry::ShapeOuterNeighborhood<Shape_, Frame_, ElementPointers_, Frame_>;

        using Region_ = Region<typename Shape_::Domain, Frame_>;

        using Regions_ = std::vector<std::shared_ptr<Region_>>;

    public:

        // explicit
        // Element(
        //     Natural const & tag
        // )
        // :
        // tag_(tag)
        // {}

        Element()
        {}

        void
        setTag(
            Natural const & tag
        )
        {
            tag_ = tag;
        }
        
        std::basic_string<Character>
        getHash()
        const
        {
            auto hash = std::basic_stringstream<Character>();
            // for (auto const & node : getInnerNeighborhood().template getComponent<geometry::Node>()<Shape_::getDimShape() - 1, 0>())
            for (auto const & node : getInnerNeighborhood().template getComponent<geometry::Node>())
            {
                hash << node->getHash();
            }
            return hash.str();
        }

        void
        addRegion(
            std::shared_ptr<Region_> const & region
        )
        {
            regions_.push_back(region);
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

        InnerNeighborhood_ inner_neighborhood_;

        OuterNeighborhood_ outer_neighborhood_;
        
        Natural tag_;

        // Region_ domain_;

        Regions_ regions_;

    };

    template<geometry::PointConcept Shape_, geometry::FrameConcept Frame_>
    struct Element<Shape_, Frame_>
    {

    private:

        template<typename... U_>
        using ElementPointers_ = std::vector<std::shared_ptr<Element<U_...>>>;

        using OuterNeighborhood_ = geometry::ShapeOuterNeighborhood<Shape_, Frame_, ElementPointers_, Frame_>;

        using Region_ = Region<typename Shape_::Domain, Frame_>;

        using Regions_ = std::vector<std::shared_ptr<Region_>>;

    public:

        // explicit
        // Element(
        //     Natural const & tag
        // )
        // :
        // tag_(tag)
        // {}

        Element()
        {}
        
        std::basic_string<Character>
        getHash()
        const
        {
            return std::to_string(this->tag_);
        }

        void
        addRegion(
            std::shared_ptr<Region_> const & region
        )
        {
            regions_.push_back(region);
        }

        void
        setTag(
            Natural const & tag
        )
        {
            tag_ = tag;
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

        OuterNeighborhood_ outer_neighborhood_;
        
        Natural tag_;
        
        geometry::Point<Frame_::getDimEuclidean()> coordinates_;

        // Region_ domain_;

        Regions_ regions_;
        
    };
    
} // namespace lolita::mesh


#endif /* AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A */
