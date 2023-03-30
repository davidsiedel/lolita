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

        using InnerNeighborhood_ = typename geometry::ShapeInnerNeighborhoodTraits<Shape_>::template InnerNeighborhood<ElementPointer_, Frame_>;

        using OuterNeighborhood_ = typename geometry::ShapeOuterNeighborhoodTraits<Shape_, Frame_>::template OuterNeighborhood<ElementPointers_, Frame_>;

        using Region_ = Region<geometry::Domain<Shape_::getDimShape()>, Frame_>;

    public:
        
        std::basic_string<Character>
        getHash()
        const
        {
            auto hash = std::basic_stringstream<Character>();
            for (auto const & node : getInnerNeighbors<Shape_::getDimShape() - 1, 0>())
            {
                hash << node->getHash();
            }
            return hash.str();
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, OuterNeighborhood_>> &
        getOuterNeighbors()
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, OuterNeighborhood_>> const &
        getOuterNeighbors()
        const
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, InnerNeighborhood_>> &
        getInnerNeighbors()
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, InnerNeighborhood_>> const &
        getInnerNeighbors()
        const
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }

        InnerNeighborhood_ inner_neighbors_;

        OuterNeighborhood_ outer_neighbors_;
        
        Natural tag_;

        Region_ domain_;

    };

    template<geometry::PointConcept Shape_, geometry::FrameConcept Frame_>
    struct Element<Shape_, Frame_>
    {

    private:

        template<typename... U_>
        using ElementPointers_ = std::vector<std::shared_ptr<Element<U_...>>>;

        using OuterNeighborhood_ = typename geometry::ShapeOuterNeighborhoodTraits<Shape_, Frame_>::template OuterNeighborhood<ElementPointers_, Frame_>;

        using Region_ = Region<geometry::Domain<Shape_::getDimShape()>, Frame_>;

    public:
        
        std::basic_string<Character>
        getHash()
        const
        {
            return std::to_string(this->tag_);
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, OuterNeighborhood_>> &
        getOuterNeighbors()
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, OuterNeighborhood_>> const &
        getOuterNeighbors()
        const
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }

        OuterNeighborhood_ outer_neighbors_;
        
        Natural tag_;
        
        geometry::Point<Frame_::getDimEuclidean()> coordinates_;

        Region_ domain_;
        
    };
    
} // namespace lolita::mesh


#endif /* AC1D9F5B_9AFE_4022_8BAA_92E35CEF446A */
